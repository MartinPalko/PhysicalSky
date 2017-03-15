using System.Collections.Generic;
using UnityEngine;
using System;

[CreateAssetMenu(fileName = "NewAtmosphereModel", menuName = "AtmosphereModel")]
public class AtmosphereModel : ScriptableObject
{
    [NonSerialized]
    bool needsRecompute = true;
    [NonSerialized]
    bool lookupTexturesDirty = true;

    enum PrecomputePass
    {
        Transmittance = 0,
        DirectIrradiance = 1,
        SingleScattering = 2,
        ScatteringDensity = 3,
        IndirectIrradiance = 4,
        MultipleScattering = 5
    }    

    // Externally visible properties
    public double SunRadius
    {
        get
        {
            return sunAngularRadius * Math.PI * 2.0;
        }
        set
        {
            sunAngularRadius = value / (Math.PI * 2.0);
        }
    }
    
    public double PlanetaryRadius
    {
        get
        {
            return bottomRadius;
        }
        set
        {
            double thickness = AtmosphereThickness;
            bottomRadius = value;
            AtmosphereThickness = thickness;
        }
    }
    
    public double AtmosphereThickness
    {
        get
        {
            return topRadius - bottomRadius;
        }
        set
        {
            topRadius = bottomRadius + Math.Max(value, 0);
        }
    }

    // Internal variables; serialized and fed to shaders
    [SerializeField] private double sunAngularRadius = 0.00872665;
    [SerializeField] private double constantSolarIrradiance = 1.5;
    [SerializeField] private double bottomRadius = 6360000.0;
    [SerializeField] private double topRadius = 6420000.0;
    [SerializeField] private double rayleigh = 1.24062e-6;
    [SerializeField] private double rayleighScaleHeight = 8000.0;
    [SerializeField] private double mieScaleHeight = 1200.0;
    [SerializeField] private double mieAngstromAlpha = 0.0;
    [SerializeField] private double mieAngstromBeta = 5.328e-3;
    [SerializeField] private double mieSingleScatteringAlbedo = 0.9;
    [SerializeField] private double miePhaseFunctionG = 0.8;
    [SerializeField] private double groundAlbedo = 0.1;
    [SerializeField] private double maxSunZenithAngle = 102.0 / 180.0 * Mathf.PI;

    // Constants
    const bool use_constant_solar_spectrum_ = false;

    const double kLambdaR = 680.0;
    const double kLambdaG = 550.0;
    const double kLambdaB = 440.0;

    // TODO: Try changing to half or double
    const RenderTextureFormat LUT_FORMAT = RenderTextureFormat.ARGBFloat;

    // Same as defined in PhysicalSkyCommon
    const int TRANSMITTANCE_TEXTURE_WIDTH = 256;
    const int TRANSMITTANCE_TEXTURE_HEIGHT = 64;
    const int SCATTERING_TEXTURE_R_SIZE = 32;
    const int SCATTERING_TEXTURE_MU_SIZE = 128;
    const int SCATTERING_TEXTURE_MU_S_SIZE = 32;
    const int SCATTERING_TEXTURE_NU_SIZE = 8;
    const int SCATTERING_TEXTURE_WIDTH = SCATTERING_TEXTURE_NU_SIZE * SCATTERING_TEXTURE_MU_S_SIZE;
    const int SCATTERING_TEXTURE_HEIGHT = SCATTERING_TEXTURE_MU_SIZE;
    const int SCATTERING_TEXTURE_DEPTH = SCATTERING_TEXTURE_R_SIZE;
    const int IRRADIANCE_TEXTURE_WIDTH = 64;
    const int IRRADIANCE_TEXTURE_HEIGHT = 16;

    // TODO: Make customizable?
    const int NUM_SCATTERING_ORDERS = 4;
    const double kLengthUnitInMeters = 1000.0;
    const int kLambdaMin = 360;
    const int kLambdaMax = 830;
    [NonSerialized]
    double[] kSolarIrradiance = {
    1.11776, 1.14259, 1.01249, 1.14716, 1.72765, 1.73054, 1.6887, 1.61253,
    1.91198, 2.03474, 2.02042, 2.02212, 1.93377, 1.95809, 1.91686, 1.8298,
    1.8685, 1.8931, 1.85149, 1.8504, 1.8341, 1.8345, 1.8147, 1.78158, 1.7533,
    1.6965, 1.68194, 1.64654, 1.6048, 1.52143, 1.55622, 1.5113, 1.474, 1.4482,
    1.41018, 1.36775, 1.34188, 1.31429, 1.28303, 1.26758, 1.2367, 1.2082,
    1.18737, 1.14683, 1.12362, 1.1058, 1.07124, 1.04992
    };

    // Computed values
    [NonSerialized]
    private double kSunSolidAngle;
    [NonSerialized]
    private List<double> wavelengths = new List<double>();
    [NonSerialized]
    private List<double> solar_irradiance = new List<double>();
    [NonSerialized]
    private List<double> rayleigh_scattering = new List<double>();
    [NonSerialized]
    private List<double> mie_scattering = new List<double>();
    [NonSerialized]
    private List<double> mie_extinction = new List<double>();
    [NonSerialized]
    private List<double> ground_albedo = new List<double>();

    // Computed textures
    [NonSerialized]
    private RenderTexture transmittanceLUT = null;
    public RenderTexture TransmittanceLUT { get { return transmittanceLUT; } }

    [NonSerialized]
    private RenderTexture scatteringLUT = null;
    public RenderTexture ScatteringLUT { get { return scatteringLUT; } }

    [NonSerialized]
    private RenderTexture irradianceLUT = null;
    public RenderTexture IrradianceLUT { get { return irradianceLUT; } }

    // Shgader and material used for precompute
    [SerializeField] [HideInInspector]
    public Shader PrecomputeShader = null;
    [NonSerialized]
    private Material PrecomputeMaterial = null;

    private double Interpolate(List<double> wavelengths, List<double> wavelength_function, double wavelength)
    {
        Debug.Assert(wavelengths.Count > 0);
        Debug.Assert(wavelength_function.Count == wavelengths.Count);

        if (wavelength < wavelengths[0])
        {
            return wavelength_function[0];
        }
        for (int i = 0; i < wavelengths.Count - 1; ++i)
        {
            if (wavelength < wavelengths[i + 1])
            {
                double u = (wavelength - wavelengths[i]) / (wavelengths[i + 1] - wavelengths[i]);
                return wavelength_function[i] * (1.0 - u) + wavelength_function[i + 1] * u;
            }
        }
        return wavelength_function[wavelength_function.Count - 1];
    }

    private Vector4 ScaleToWavelengths(List<double> v, double scale)
    {
        return new Vector4(
            (float)(Interpolate(wavelengths, v, kLambdaR) * scale),
            (float)(Interpolate(wavelengths, v, kLambdaG) * scale),
            (float)(Interpolate(wavelengths, v, kLambdaB) * scale),
            1.0f);
    }

    public void AllocateLookupTextures()
    {
        if (!transmittanceLUT)
            transmittanceLUT = new RenderTexture(TRANSMITTANCE_TEXTURE_WIDTH, TRANSMITTANCE_TEXTURE_HEIGHT, 0, LUT_FORMAT);
        if (!transmittanceLUT.IsCreated())
        {
            transmittanceLUT.useMipMap = false;
            transmittanceLUT.Create();
        }

        if (!scatteringLUT)
            scatteringLUT = new RenderTexture(SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT, 0, LUT_FORMAT);
        if (!scatteringLUT.IsCreated())
        {
            scatteringLUT.volumeDepth = SCATTERING_TEXTURE_DEPTH;
            scatteringLUT.useMipMap = false;
            scatteringLUT.dimension = UnityEngine.Rendering.TextureDimension.Tex3D;
            scatteringLUT.enableRandomWrite = true;
            scatteringLUT.Create();
        }

        if (!irradianceLUT)
            irradianceLUT = new RenderTexture(IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT, 0, LUT_FORMAT);
        if (!irradianceLUT.IsCreated())
        {
            irradianceLUT.useMipMap = false;
            irradianceLUT.Create();
        }
    }    

    private void Blit(RenderTexture dest, Material mat, int pass)
    {
        Graphics.Blit(null, dest, mat, pass);
    }

    private void BlitWithDummy(Material mat, int pass, int width, int height)
    {
        RenderTexture dummy = RenderTexture.GetTemporary(width, height, 0, RenderTextureFormat.R8);
        Blit(dummy, mat, pass);
        RenderTexture.ReleaseTemporary(dummy);
    }

    public void SetAtmosphereUniforms(Material mat)
    {
        mat.SetVector("_solar_irradiance", ScaleToWavelengths(solar_irradiance, 1.0));
        mat.SetFloat("_sun_angular_radius", (float)sunAngularRadius);
        mat.SetFloat("_bottom_radius", (float)(bottomRadius / kLengthUnitInMeters));
        mat.SetFloat("_top_radius", (float)(topRadius / kLengthUnitInMeters));
        mat.SetFloat("_rayleigh_scale_height", (float)(rayleighScaleHeight / kLengthUnitInMeters));
        mat.SetVector("_rayleigh_scattering", ScaleToWavelengths(rayleigh_scattering, kLengthUnitInMeters));
        mat.SetFloat("_mie_scale_height", (float)(mieScaleHeight / kLengthUnitInMeters));
        mat.SetVector("_mie_scattering", ScaleToWavelengths(mie_scattering, kLengthUnitInMeters));
        mat.SetVector("_mie_extinction", ScaleToWavelengths(mie_extinction, kLengthUnitInMeters));
        mat.SetFloat("_mie_phase_function_g", (float)miePhaseFunctionG);
        mat.SetVector("_ground_albedo", ScaleToWavelengths(ground_albedo, 1.0));
        mat.SetFloat("_mu_s_min", (float)Math.Cos(maxSunZenithAngle));

        mat.SetVector("sun_radiance", new Vector3((float)kSolarIrradiance[0], (float)kSolarIrradiance[1], (float)kSolarIrradiance[2]) / (float)kSunSolidAngle);
        mat.SetVector("sun_size", new Vector3(Mathf.Tan((float)sunAngularRadius), Mathf.Cos((float)sunAngularRadius), (float)sunAngularRadius));
    }

    public void ComputeLookupTextures()
    {
        Debug.Log("Computing LUT");
        float timerStartCompute = Time.realtimeSinceStartup;

        if (SystemInfo.graphicsShaderLevel < 50)
        {
            Debug.LogError("Computing atmospheric lookup textures requires shader model 5.0 or higher!");
            return;
        }

        wavelengths.Clear();
        solar_irradiance.Clear();
        rayleigh_scattering.Clear();
        mie_scattering.Clear();
        mie_extinction.Clear();
        ground_albedo.Clear();

        kSunSolidAngle = Math.PI * sunAngularRadius * sunAngularRadius;

        for (int l = kLambdaMin; l <= kLambdaMax; l += 10)
        {
            double lambda = l * 1e-3;  // micro-meters
            double mie = mieAngstromBeta / mieScaleHeight * Math.Pow(lambda, -mieAngstromAlpha);
            wavelengths.Add(l);

            if (use_constant_solar_spectrum_)
                solar_irradiance.Add(constantSolarIrradiance);
            else
                solar_irradiance.Add(kSolarIrradiance[(l - kLambdaMin) / 10]);

            rayleigh_scattering.Add(rayleigh * Math.Pow(lambda, -4));
            mie_scattering.Add(mie * mieSingleScatteringAlbedo);
            mie_extinction.Add(mie);
            ground_albedo.Add(groundAlbedo);
        }

        if (!PrecomputeMaterial)
            PrecomputeMaterial = new Material(PrecomputeShader);

        SetAtmosphereUniforms(PrecomputeMaterial);

        RenderTexture DeltaIrradianceTexture = new RenderTexture(IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT, 0, LUT_FORMAT, RenderTextureReadWrite.Linear);
        DeltaIrradianceTexture.useMipMap = false;
        DeltaIrradianceTexture.enableRandomWrite = true;
        DeltaIrradianceTexture.Create();

        RenderTexture DeltaRayleighScatteringTexture = new RenderTexture(SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT, 0, LUT_FORMAT, RenderTextureReadWrite.Linear);
        DeltaRayleighScatteringTexture.volumeDepth = SCATTERING_TEXTURE_DEPTH;
        DeltaRayleighScatteringTexture.useMipMap = false;
        DeltaRayleighScatteringTexture.dimension = UnityEngine.Rendering.TextureDimension.Tex3D;
        DeltaRayleighScatteringTexture.enableRandomWrite = true;
        DeltaRayleighScatteringTexture.Create();

        RenderTexture DeltaMieScatteringTexture = new RenderTexture(SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT, 0, LUT_FORMAT, RenderTextureReadWrite.Linear);
        DeltaMieScatteringTexture.volumeDepth = SCATTERING_TEXTURE_DEPTH;
        DeltaMieScatteringTexture.useMipMap = false;
        DeltaMieScatteringTexture.dimension = UnityEngine.Rendering.TextureDimension.Tex3D;
        DeltaMieScatteringTexture.enableRandomWrite = true;
        DeltaMieScatteringTexture.Create();

        RenderTexture DeltaScatteringDensityTexture = new RenderTexture(SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT, 0, LUT_FORMAT, RenderTextureReadWrite.Linear);
        DeltaScatteringDensityTexture.volumeDepth = SCATTERING_TEXTURE_DEPTH;
        DeltaScatteringDensityTexture.useMipMap = false;
        DeltaScatteringDensityTexture.dimension = UnityEngine.Rendering.TextureDimension.Tex3D;
        DeltaScatteringDensityTexture.enableRandomWrite = true;
        DeltaScatteringDensityTexture.Create();

        RenderTexture DeltaMultipleScatteringTexture = DeltaRayleighScatteringTexture;

        // Allocate final lookup textures
        AllocateLookupTextures();

        // Compute Transmittance LUT
        Blit(transmittanceLUT, PrecomputeMaterial, (int)PrecomputePass.Transmittance);
        PrecomputeMaterial.SetTexture("transmittance_texture", transmittanceLUT); // Set for subsequent shaders to read

        // Compute Direct Irradiance into DeltaIrradianceTexture and Initialize irradianceLUT with 0
        Graphics.ClearRandomWriteTargets();
        Graphics.SetRandomWriteTarget(1, DeltaIrradianceTexture);
        Blit(irradianceLUT, PrecomputeMaterial, (int)PrecomputePass.DirectIrradiance);

        // Compute Raylie and Mie Single Scattering and store them in DeltaRayleighScatteringTexture, DeltaMieScatteringTexture as well as scatteringLUT
        Graphics.ClearRandomWriteTargets();
        Graphics.SetRandomWriteTarget(1, DeltaRayleighScatteringTexture);
        Graphics.SetRandomWriteTarget(2, DeltaMieScatteringTexture);
        Graphics.SetRandomWriteTarget(3, scatteringLUT);
        PrecomputeMaterial.SetTexture("transmittance_texture", transmittanceLUT);
        BlitWithDummy(PrecomputeMaterial, (int)PrecomputePass.SingleScattering, SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT);

        //Compute to the nth order of scattering, in sequence
        for (int scatteringOrder = 2; scatteringOrder <= NUM_SCATTERING_ORDERS; ++scatteringOrder)
        {
            PrecomputeMaterial.SetTexture("transmittance_texture", transmittanceLUT);
            PrecomputeMaterial.SetTexture("single_rayleigh_scattering_texture", DeltaRayleighScatteringTexture);
            PrecomputeMaterial.SetTexture("single_mie_scattering_texture", DeltaMieScatteringTexture);
            PrecomputeMaterial.SetTexture("multiple_scattering_texture", DeltaMultipleScatteringTexture);
            PrecomputeMaterial.SetTexture("irradiance_texture", DeltaIrradianceTexture);
            PrecomputeMaterial.SetInt("scattering_order", scatteringOrder);

            // Compute the scattering density, and store it in DeltaScatteringDensityTexture.
            Graphics.ClearRandomWriteTargets();
            Graphics.SetRandomWriteTarget(1, DeltaScatteringDensityTexture);
            BlitWithDummy(PrecomputeMaterial, (int)PrecomputePass.ScatteringDensity, SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT);

            // Compute the indirect irradiance, store it in DeltaIrradianceTexture and accumulate it in irradianceLUT.
            Graphics.ClearRandomWriteTargets();
            Graphics.SetRandomWriteTarget(1, DeltaIrradianceTexture);
            Graphics.SetRandomWriteTarget(2, irradianceLUT);
            BlitWithDummy(PrecomputeMaterial, (int)PrecomputePass.IndirectIrradiance, IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT);

            // Compute the multiple scattering, store it in DeltaMultipleScatteringTexture, and accumulate it in scatteringLUT.
            Graphics.ClearRandomWriteTargets();
            Graphics.SetRandomWriteTarget(1, DeltaMultipleScatteringTexture);
            Graphics.SetRandomWriteTarget(2, scatteringLUT);
            PrecomputeMaterial.SetTexture("scattering_density_texture", DeltaScatteringDensityTexture);
            BlitWithDummy(PrecomputeMaterial, (int)PrecomputePass.MultipleScattering, SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT);
        }

        Graphics.ClearRandomWriteTargets();
        RenderTexture.active = null;
        
        // Release temporary textures.
        DeltaIrradianceTexture.Release();
        DeltaRayleighScatteringTexture.Release();
        DeltaMieScatteringTexture.Release();
        DeltaScatteringDensityTexture.Release();

        float timerEndCompute = Time.realtimeSinceStartup;
        Debug.Log("Computed atmospheric lookup textures in " + (timerEndCompute - timerStartCompute) * 1000.0f + "ms");
        needsRecompute = false;
    }

    void ReleaseLookupTextures()
    {
        if (transmittanceLUT && transmittanceLUT.IsCreated())
            transmittanceLUT.Release();

        if (scatteringLUT && scatteringLUT.IsCreated())
            scatteringLUT.Release();

        if (irradianceLUT && irradianceLUT.IsCreated())
            irradianceLUT.Release();
    }

    private void Awake()
    {
        needsRecompute = true;
    }

    private void OnDestroy()
    {
        ReleaseLookupTextures();
    }
}
