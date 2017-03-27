using System.Collections.Generic;
using UnityEngine;
using System;

namespace PhysicalSky
{
    [CreateAssetMenu(fileName = "NewAtmosphereModel", menuName = "AtmosphereModel")]
    public class AtmosphereModel : ScriptableObject, IAtmosphereModel
    {
        enum PrecomputePass
        {
            Transmittance = 0,
            DirectIrradiance = 1,
            SingleScattering = 2,
            ScatteringDensity = 3,
            IndirectIrradiance = 4,
            MultipleScattering = 5
        }
        
        private bool needsRecompute = true;
        public bool NeedsRecompute { get { return needsRecompute || TexturesInvalid(); } }

        [SerializeField]
        private float constantSolarIrradiance = 1.5f;
        public float ConstantSolarIrradiance
        {
            get { return constantSolarIrradiance; }
            set
            {
                if (constantSolarIrradiance != value)
                {
                    constantSolarIrradiance = value;
                    needsRecompute = true;
                }
            }
        }

        [SerializeField]
        private float sunAngularRadius = 0.00872665f;
        public float SunAngularRadius
        {
            get { return sunAngularRadius; }
            set
            {
                if (sunAngularRadius != value)
                {
                    sunAngularRadius = value;
                    needsRecompute = true;
                }
            }
        }

        [SerializeField]
        private float planetaryRadius = 6360000.0f;
        public float PlanetaryRadius
        {
            get
            {
                return planetaryRadius;
            }
            set
            {
                float thickness = AtmosphereThickness;
                planetaryRadius = value;
                AtmosphereThickness = thickness;
            }
        }

        [SerializeField]
        private float atmosphereThickness = 60000.0f;
        public float AtmosphereThickness
        {
            get { return atmosphereThickness; }
            set
            {
                float newValue = Math.Max(value, 0);
                if (newValue != atmosphereThickness)
                {
                    atmosphereThickness = newValue;
                    needsRecompute = true;
                }

            }
        }

        [SerializeField]
        private float rayleigh = 1.24062e-6f;
        public float Rayleigh
        {
            get { return rayleigh; }
            set
            {
                if (rayleigh != value)
                {
                    rayleigh = value;
                    needsRecompute = true;
                }
            }
        }

        [SerializeField]
        private float rayleighScaleHeight = 8000.0f;
        public float RayleighScaleHeight
        {
            get { return rayleighScaleHeight; }
            set
            {
                if (rayleighScaleHeight != value)
                {
                    rayleighScaleHeight = value;
                    needsRecompute = true;
                }
            }
        }

        [SerializeField]
        private float mieScaleHeight = 1200.0f;
        public float MieScaleHeight
        {
            get { return mieScaleHeight; }
            set
            {
                if (mieScaleHeight != value)
                {
                    mieScaleHeight = value;
                    needsRecompute = true;
                }
            }
        }

        [SerializeField]
        private float mieAngstromAlpha = 0.0f;
        public float MieAngstromAlpha
        {
            get { return mieAngstromAlpha; }
            set
            {
                if (mieAngstromAlpha != value)
                {
                    mieAngstromAlpha = value;
                    needsRecompute = true;
                }
            }
        }

        [SerializeField]
        private float mieAngstromBeta = 5.328e-3f;
        public float MieAngstromBeta
        {
            get { return mieAngstromBeta; }
            set
            {
                if (mieAngstromBeta != value)
                {
                    mieAngstromBeta = value;
                    needsRecompute = true;
                }
            }
        }

        [SerializeField]
        private float mieSingleScatteringAlbedo = 0.9f;
        public float MieSingleScatteringAlbedo
        {
            get { return mieSingleScatteringAlbedo; }
            set
            {
                if (mieSingleScatteringAlbedo != value)
                {
                    mieSingleScatteringAlbedo = value;
                    needsRecompute = true;
                }
            }
        }

        [SerializeField]
        private float miePhaseFunctionG = 0.8f;
        public float MiePhaseFunctionG
        {
            get { return miePhaseFunctionG; }
            set
            {
                if (miePhaseFunctionG != value)
                {
                    miePhaseFunctionG = value;
                    needsRecompute = true;
                }
            }
        }

        [SerializeField]
        private float groundAlbedo = 0.1f;
        public float GroundAlbedo
        {
            get { return groundAlbedo; }
            set
            {
                if (groundAlbedo != value)
                {
                    groundAlbedo = value;
                    needsRecompute = true;
                }
            }
        }

        [SerializeField]
        private float maxSunZenithAngle = 102.0f / 180.0f * Mathf.PI;
        public float MaxSunZenithAngle
        {
            get { return maxSunZenithAngle; }
            set
            {
                if (maxSunZenithAngle != value)
                {
                    maxSunZenithAngle = value;
                    needsRecompute = true;
                }
            }
        }

        // Constants
        const bool use_constant_solar_spectrum_ = true;

        const float kLambdaR = 680.0f;
        const float kLambdaG = 550.0f;
        const float kLambdaB = 440.0f;

        // TODO: Try changing to half or float
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
        const float kLengthUnitInMeters = 1000.0f;
        const int kLambdaMin = 360;
        const int kLambdaMax = 830;
        [NonSerialized]
        float[] kSolarIrradiance = {
    1.11776f, 1.14259f, 1.01249f, 1.14716f, 1.72765f, 1.73054f, 1.6887f, 1.61253f,
    1.91198f, 2.03474f, 2.02042f, 2.02212f, 1.93377f, 1.95809f, 1.91686f, 1.8298f,
    1.8685f, 1.8931f, 1.85149f, 1.8504f, 1.8341f, 1.8345f, 1.8147f, 1.78158f, 1.7533f,
    1.6965f, 1.68194f, 1.64654f, 1.6048f, 1.52143f, 1.55622f, 1.5113f, 1.474f, 1.4482f,
    1.41018f, 1.36775f, 1.34188f, 1.31429f, 1.28303f, 1.26758f, 1.2367f, 1.2082f,
    1.18737f, 1.14683f, 1.12362f, 1.1058f, 1.07124f, 1.04992f
    };

        // Computed values
        [NonSerialized]
        private float kSunSolidAngle;
        [NonSerialized]
        private List<float> wavelengths = new List<float>();
        [NonSerialized]
        private List<float> solar_irradiance = new List<float>();
        [NonSerialized]
        private List<float> rayleigh_scattering = new List<float>();
        [NonSerialized]
        private List<float> mie_scattering = new List<float>();
        [NonSerialized]
        private List<float> mie_extinction = new List<float>();
        [NonSerialized]
        private List<float> ground_albedo = new List<float>();

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
        [SerializeField]
        [HideInInspector]
        public Shader PrecomputeShader = null;
        [NonSerialized]
        private Material PrecomputeMaterial = null;

        private float Interpolate(List<float> wavelengths, List<float> wavelength_function, float wavelength)
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
                    float u = (wavelength - wavelengths[i]) / (wavelengths[i + 1] - wavelengths[i]);
                    return wavelength_function[i] * (1.0f - u) + wavelength_function[i + 1] * u;
                }
            }
            return wavelength_function[wavelength_function.Count - 1];
        }

        private Vector4 ScaleToWavelengths(List<float> v, float scale)
        {
            return new Vector4(
                (float)(Interpolate(wavelengths, v, kLambdaR) * scale),
                (float)(Interpolate(wavelengths, v, kLambdaG) * scale),
                (float)(Interpolate(wavelengths, v, kLambdaB) * scale),
                1.0f);
        }

        private void AllocateLookupTextures()
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

        private void ReleaseLookupTextures()
        {
            if (transmittanceLUT && transmittanceLUT.IsCreated())
                transmittanceLUT.Release();

            if (scatteringLUT && scatteringLUT.IsCreated())
                scatteringLUT.Release();

            if (irradianceLUT && irradianceLUT.IsCreated())
                irradianceLUT.Release();
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

        public void SetShaderUniforms(Material m)
        {
            m.SetVector("_solar_irradiance", ScaleToWavelengths(solar_irradiance, 1.0f));
            m.SetFloat("_sun_angular_radius", sunAngularRadius);
            m.SetFloat("_bottom_radius", planetaryRadius / kLengthUnitInMeters);
            m.SetFloat("_top_radius", (planetaryRadius + atmosphereThickness) / kLengthUnitInMeters);
            m.SetFloat("_rayleigh_scale_height", rayleighScaleHeight / kLengthUnitInMeters);
            m.SetVector("_rayleigh_scattering", ScaleToWavelengths(rayleigh_scattering, kLengthUnitInMeters));
            m.SetFloat("_mie_scale_height", mieScaleHeight / kLengthUnitInMeters);
            m.SetVector("_mie_scattering", ScaleToWavelengths(mie_scattering, kLengthUnitInMeters));
            m.SetVector("_mie_extinction", ScaleToWavelengths(mie_extinction, kLengthUnitInMeters));
            m.SetFloat("_mie_phase_function_g", miePhaseFunctionG);
            m.SetVector("_ground_albedo", ScaleToWavelengths(ground_albedo, 1.0f));
            m.SetFloat("_mu_s_min", Mathf.Cos(maxSunZenithAngle));
            m.SetVector("sun_radiance", new Vector3(kSolarIrradiance[0], kSolarIrradiance[1], kSolarIrradiance[2]) / kSunSolidAngle);
            m.SetVector("sun_size", new Vector3(Mathf.Tan(sunAngularRadius), Mathf.Cos(sunAngularRadius), sunAngularRadius));
        }

        public void Compute()
        {
            Debug.Log("Computing Atmospheric Lookup Textures");
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

            kSunSolidAngle = Mathf.PI * sunAngularRadius * sunAngularRadius;

            for (int l = kLambdaMin; l <= kLambdaMax; l += 10)
            {
                float lambda = l * 1e-3f;  // micro-meters
                float mie = mieAngstromBeta / mieScaleHeight * Mathf.Pow(lambda, -mieAngstromAlpha);
                wavelengths.Add(l);

                if (use_constant_solar_spectrum_)
                    solar_irradiance.Add(constantSolarIrradiance);
                else
                    solar_irradiance.Add(kSolarIrradiance[(l - kLambdaMin) / 10]);

                rayleigh_scattering.Add(rayleigh * Mathf.Pow(lambda, -4));
                mie_scattering.Add(mie * mieSingleScatteringAlbedo);
                mie_extinction.Add(mie);
                ground_albedo.Add(groundAlbedo);
            }

            if (!PrecomputeMaterial)
                PrecomputeMaterial = new Material(PrecomputeShader);

            SetShaderUniforms(PrecomputeMaterial);

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

        public void ReleaseResources()
        {
            Debug.Log("Released Atmosphere Resources");
            ReleaseLookupTextures();
            needsRecompute = true;
        }

        public bool TexturesInvalid()
        {
            if (!transmittanceLUT || !transmittanceLUT.IsCreated())
                return true;
            else if (!scatteringLUT || !scatteringLUT.IsCreated())
                return true;

            else if (!irradianceLUT || !irradianceLUT.IsCreated())
                return true;
            else
                return false;
        }

        private void Awake()
        {
            needsRecompute = true;
        }

        private void OnDestroy()
        {
            ReleaseResources();
        }
    }
}