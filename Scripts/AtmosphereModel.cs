using System.Collections.Generic;
using UnityEngine;
using System;

namespace PhysicalSky
{
    [CreateAssetMenu(fileName = "NewAtmosphereModel", menuName = "PhysicalSky/AtmosphereModel")]
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

        private bool m_needsRecompute = true;
        public bool NeedsRecompute { get { return m_needsRecompute || TexturesInvalid(); } }

        [SerializeField]
        private float m_constantSolarIrradiance = 1.5f;
        public float ConstantSolarIrradiance
        {
            get { return m_constantSolarIrradiance; }
            set
            {
                if (m_constantSolarIrradiance != value)
                {
                    m_constantSolarIrradiance = value;
                    m_needsRecompute = true;
                }
            }
        }

        [SerializeField]
        private float m_sunAngularRadius = 0.00872665f;
        public float SunAngularRadius
        {
            get { return m_sunAngularRadius; }
            set
            {
                if (m_sunAngularRadius != value)
                {
                    m_sunAngularRadius = value;
                    m_needsRecompute = true;
                }
            }
        }

        [SerializeField]
        private float m_planetaryRadius = 6360000.0f;
        public float PlanetaryRadius
        {
            get
            {
                return m_planetaryRadius;
            }
            set
            {
                float thickness = AtmosphereThickness;
                m_planetaryRadius = value;
                AtmosphereThickness = thickness;
            }
        }

        [SerializeField]
        private float m_atmosphereThickness = 60000.0f;
        public float AtmosphereThickness
        {
            get { return m_atmosphereThickness; }
            set
            {
                float newValue = Math.Max(value, 0);
                if (newValue != m_atmosphereThickness)
                {
                    m_atmosphereThickness = newValue;
                    m_needsRecompute = true;
                }

            }
        }

        [SerializeField]
        private float m_rayleigh = 1.24062e-6f;
        public float Rayleigh
        {
            get { return m_rayleigh; }
            set
            {
                if (m_rayleigh != value)
                {
                    m_rayleigh = value;
                    m_needsRecompute = true;
                }
            }
        }

        [SerializeField]
        private float m_rayleighScaleHeight = 8000.0f;
        public float RayleighScaleHeight
        {
            get { return m_rayleighScaleHeight; }
            set
            {
                if (m_rayleighScaleHeight != value)
                {
                    m_rayleighScaleHeight = value;
                    m_needsRecompute = true;
                }
            }
        }

        [SerializeField]
        private float m_mieScaleHeight = 1200.0f;
        public float MieScaleHeight
        {
            get { return m_mieScaleHeight; }
            set
            {
                if (m_mieScaleHeight != value)
                {
                    m_mieScaleHeight = value;
                    m_needsRecompute = true;
                }
            }
        }

        [SerializeField]
        private float m_mieAngstromAlpha = 0.0f;
        public float MieAngstromAlpha
        {
            get { return m_mieAngstromAlpha; }
            set
            {
                if (m_mieAngstromAlpha != value)
                {
                    m_mieAngstromAlpha = value;
                    m_needsRecompute = true;
                }
            }
        }

        [SerializeField]
        private float m_mieAngstromBeta = 5.328e-3f;
        public float MieAngstromBeta
        {
            get { return m_mieAngstromBeta; }
            set
            {
                if (m_mieAngstromBeta != value)
                {
                    m_mieAngstromBeta = value;
                    m_needsRecompute = true;
                }
            }
        }

        [SerializeField]
        private float m_mieSingleScatteringAlbedo = 0.9f;
        public float MieSingleScatteringAlbedo
        {
            get { return m_mieSingleScatteringAlbedo; }
            set
            {
                if (m_mieSingleScatteringAlbedo != value)
                {
                    m_mieSingleScatteringAlbedo = value;
                    m_needsRecompute = true;
                }
            }
        }

        [SerializeField]
        private float m_miePhaseFunctionG = 0.8f;
        public float MiePhaseFunctionG
        {
            get { return m_miePhaseFunctionG; }
            set
            {
                if (m_miePhaseFunctionG != value)
                {
                    m_miePhaseFunctionG = value;
                    m_needsRecompute = true;
                }
            }
        }

        [SerializeField]
        private float m_groundAlbedo = 0.1f;
        public float GroundAlbedo
        {
            get { return m_groundAlbedo; }
            set
            {
                if (m_groundAlbedo != value)
                {
                    m_groundAlbedo = value;
                    m_needsRecompute = true;
                }
            }
        }

        [SerializeField]
        private float m_maxSunZenithAngle = 102.0f / 180.0f * Mathf.PI;
        public float MaxSunZenithAngle
        {
            get { return m_maxSunZenithAngle; }
            set
            {
                if (m_maxSunZenithAngle != value)
                {
                    m_maxSunZenithAngle = value;
                    m_needsRecompute = true;
                }
            }
        }

        // Constants
        const bool USE_CONSTANT_SOLAR_SPECTRUM = true;

        const float LAMBDA_R = 680.0f;
        const float LAMBDA_G = 550.0f;
        const float LAMBDA_B = 440.0f;

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

        // TODO: Make num scattering orders customizable?
        const int NUM_SCATTERING_ORDERS = 4;
        const float LENGTH_UNIT_IN_METERS = 1000.0f;
        const int LAMBDA_MIN = 360;
        const int LAMBDA_MAX = 830;
        static readonly float[] SOLOAR_IRRADIANCE = {
    1.11776f, 1.14259f, 1.01249f, 1.14716f, 1.72765f, 1.73054f, 1.6887f, 1.61253f,
    1.91198f, 2.03474f, 2.02042f, 2.02212f, 1.93377f, 1.95809f, 1.91686f, 1.8298f,
    1.8685f, 1.8931f, 1.85149f, 1.8504f, 1.8341f, 1.8345f, 1.8147f, 1.78158f, 1.7533f,
    1.6965f, 1.68194f, 1.64654f, 1.6048f, 1.52143f, 1.55622f, 1.5113f, 1.474f, 1.4482f,
    1.41018f, 1.36775f, 1.34188f, 1.31429f, 1.28303f, 1.26758f, 1.2367f, 1.2082f,
    1.18737f, 1.14683f, 1.12362f, 1.1058f, 1.07124f, 1.04992f
    };

        // Computed values
        [NonSerialized]
        private float m_sunSolidAngle;
        [NonSerialized]
        private List<float> m_wavelengths = new List<float>();
        [NonSerialized]
        private List<float> m_solarIrradiance = new List<float>();
        [NonSerialized]
        private List<float> m_rayleighScattering = new List<float>();
        [NonSerialized]
        private List<float> m_mieScattering = new List<float>();
        [NonSerialized]
        private List<float> m_mieExtinction = new List<float>();

        // Computed textures
        [NonSerialized]
        private RenderTexture m_transmittanceLUT = null;
        public RenderTexture TransmittanceLUT { get { return m_transmittanceLUT; } }

        [NonSerialized]
        private RenderTexture m_scatteringLUT = null;
        public RenderTexture ScatteringLUT { get { return m_scatteringLUT; } }

        [NonSerialized]
        private RenderTexture m_irradianceLUT = null;
        public RenderTexture IrradianceLUT { get { return m_irradianceLUT; } }

        // Shgader and material used for precompute
        [SerializeField]
        private Shader m_PrecomputeShader = null;
        public Shader PrecomputeShader { get { return m_PrecomputeShader; } set { m_PrecomputeShader = value; } }

        [NonSerialized]
        private Material m_precomputeMaterial = null;

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
                (float)(Interpolate(m_wavelengths, v, LAMBDA_R) * scale),
                (float)(Interpolate(m_wavelengths, v, LAMBDA_G) * scale),
                (float)(Interpolate(m_wavelengths, v, LAMBDA_B) * scale),
                1.0f);
        }

        private void AllocateLookupTextures()
        {
            if (!m_transmittanceLUT)
                m_transmittanceLUT = new RenderTexture(TRANSMITTANCE_TEXTURE_WIDTH, TRANSMITTANCE_TEXTURE_HEIGHT, 0, LUT_FORMAT, RenderTextureReadWrite.Linear);
            if (!m_transmittanceLUT.IsCreated())
            {
                m_transmittanceLUT.useMipMap = false;
                m_transmittanceLUT.Create();
            }

            if (!m_scatteringLUT)
                m_scatteringLUT = new RenderTexture(SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT, 0, LUT_FORMAT, RenderTextureReadWrite.Linear);
            if (!m_scatteringLUT.IsCreated())
            {
                m_scatteringLUT.volumeDepth = SCATTERING_TEXTURE_DEPTH;
                m_scatteringLUT.useMipMap = false;
                m_scatteringLUT.dimension = UnityEngine.Rendering.TextureDimension.Tex3D;
                m_scatteringLUT.enableRandomWrite = true;
                m_scatteringLUT.Create();
            }

            if (!m_irradianceLUT)
                m_irradianceLUT = new RenderTexture(IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT, 0, LUT_FORMAT, RenderTextureReadWrite.Linear);
            if (!m_irradianceLUT.IsCreated())
            {
                m_irradianceLUT.useMipMap = false;
                m_irradianceLUT.enableRandomWrite = true;
                m_irradianceLUT.Create();
            }
        }

        private void ReleaseLookupTextures()
        {
            if (m_transmittanceLUT && m_transmittanceLUT.IsCreated())
                m_transmittanceLUT.Release();

            if (m_scatteringLUT && m_scatteringLUT.IsCreated())
                m_scatteringLUT.Release();

            if (m_irradianceLUT && m_irradianceLUT.IsCreated())
                m_irradianceLUT.Release();
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
            m.SetVector("_solar_irradiance", ScaleToWavelengths(m_solarIrradiance, 1.0f));
            m.SetFloat("_sun_angular_radius", m_sunAngularRadius);
            m.SetFloat("_bottom_radius", m_planetaryRadius / LENGTH_UNIT_IN_METERS);
            m.SetFloat("_top_radius", (m_planetaryRadius + m_atmosphereThickness) / LENGTH_UNIT_IN_METERS);
            m.SetFloat("_rayleigh_scale_height", m_rayleighScaleHeight / LENGTH_UNIT_IN_METERS);
            m.SetVector("_rayleigh_scattering", ScaleToWavelengths(m_rayleighScattering, LENGTH_UNIT_IN_METERS));
            m.SetFloat("_mie_scale_height", m_mieScaleHeight / LENGTH_UNIT_IN_METERS);
            m.SetVector("_mie_scattering", ScaleToWavelengths(m_mieScattering, LENGTH_UNIT_IN_METERS));
            m.SetVector("_mie_extinction", ScaleToWavelengths(m_mieExtinction, LENGTH_UNIT_IN_METERS));
            m.SetFloat("_mie_phase_function_g", m_miePhaseFunctionG);
            m.SetVector("_ground_albedo", new Vector3(m_groundAlbedo, m_groundAlbedo, m_groundAlbedo));
            m.SetFloat("_mu_s_min", Mathf.Cos(m_maxSunZenithAngle));
            m.SetVector("sun_radiance", new Vector3(SOLOAR_IRRADIANCE[0], SOLOAR_IRRADIANCE[1], SOLOAR_IRRADIANCE[2]) / m_sunSolidAngle);
            m.SetVector("sun_size", new Vector3(Mathf.Tan(m_sunAngularRadius), Mathf.Cos(m_sunAngularRadius), m_sunAngularRadius));
        }

        public void Compute()
        {
#if PHYSICAL_SKY_DEBUG
            Debug.Log("Computing Atmospheric Lookup Textures");
#endif
            float timerStartCompute = Time.realtimeSinceStartup;

            if (SystemInfo.graphicsShaderLevel < 50)
            {
                Debug.LogError("Computing atmospheric lookup textures requires shader model 5.0 or higher!");
                return;
            }

            m_wavelengths.Clear();
            m_solarIrradiance.Clear();
            m_rayleighScattering.Clear();
            m_mieScattering.Clear();
            m_mieExtinction.Clear();

            m_sunSolidAngle = Mathf.PI * m_sunAngularRadius * m_sunAngularRadius;

            for (int l = LAMBDA_MIN; l <= LAMBDA_MAX; l += 10)
            {
                float lambda = l * 1e-3f;  // micro-meters
                float mie = m_mieAngstromBeta / m_mieScaleHeight * Mathf.Pow(lambda, -m_mieAngstromAlpha);
                m_wavelengths.Add(l);

                if (USE_CONSTANT_SOLAR_SPECTRUM)
                    m_solarIrradiance.Add(m_constantSolarIrradiance);
                else
                    m_solarIrradiance.Add(SOLOAR_IRRADIANCE[(l - LAMBDA_MIN) / 10]);

                m_rayleighScattering.Add(m_rayleigh * Mathf.Pow(lambda, -4));
                m_mieScattering.Add(mie * m_mieSingleScatteringAlbedo);
                m_mieExtinction.Add(mie);
            }

            if (!m_precomputeMaterial)
                m_precomputeMaterial = new Material(PrecomputeShader);

            SetShaderUniforms(m_precomputeMaterial);

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
            Blit(m_transmittanceLUT, m_precomputeMaterial, (int)PrecomputePass.Transmittance);
            m_precomputeMaterial.SetTexture("transmittance_texture", m_transmittanceLUT); // Set for subsequent shaders to read

            // Compute Direct Irradiance into DeltaIrradianceTexture and Initialize irradianceLUT with 0
            Graphics.ClearRandomWriteTargets();
            Graphics.SetRandomWriteTarget(1, DeltaIrradianceTexture);
            Blit(m_irradianceLUT, m_precomputeMaterial, (int)PrecomputePass.DirectIrradiance);

            // Compute Raylie and Mie Single Scattering and store them in DeltaRayleighScatteringTexture, DeltaMieScatteringTexture as well as scatteringLUT
            Graphics.ClearRandomWriteTargets();
            Graphics.SetRandomWriteTarget(1, DeltaRayleighScatteringTexture);
            Graphics.SetRandomWriteTarget(2, DeltaMieScatteringTexture);
            Graphics.SetRandomWriteTarget(3, m_scatteringLUT);
            m_precomputeMaterial.SetTexture("transmittance_texture", m_transmittanceLUT);
            BlitWithDummy(m_precomputeMaterial, (int)PrecomputePass.SingleScattering, SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT);

            //Compute to the nth order of scattering, in sequence
            for (int scatteringOrder = 2; scatteringOrder <= NUM_SCATTERING_ORDERS; ++scatteringOrder)
            {
                m_precomputeMaterial.SetTexture("transmittance_texture", m_transmittanceLUT);
                m_precomputeMaterial.SetTexture("single_rayleigh_scattering_texture", DeltaRayleighScatteringTexture);
                m_precomputeMaterial.SetTexture("single_mie_scattering_texture", DeltaMieScatteringTexture);
                m_precomputeMaterial.SetTexture("multiple_scattering_texture", DeltaMultipleScatteringTexture);
                m_precomputeMaterial.SetTexture("irradiance_texture", DeltaIrradianceTexture);
                m_precomputeMaterial.SetInt("scattering_order", scatteringOrder);

                // Compute the scattering density, and store it in DeltaScatteringDensityTexture.
                Graphics.ClearRandomWriteTargets();
                Graphics.SetRandomWriteTarget(1, DeltaScatteringDensityTexture);
                BlitWithDummy(m_precomputeMaterial, (int)PrecomputePass.ScatteringDensity, SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT);

                // Compute the indirect irradiance, store it in DeltaIrradianceTexture and accumulate it in irradianceLUT.
                Graphics.ClearRandomWriteTargets();
                Graphics.SetRandomWriteTarget(1, DeltaIrradianceTexture);
                Graphics.SetRandomWriteTarget(2, m_irradianceLUT);
                BlitWithDummy(m_precomputeMaterial, (int)PrecomputePass.IndirectIrradiance, IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT);

                // Compute the multiple scattering, store it in DeltaMultipleScatteringTexture, and accumulate it in scatteringLUT.
                Graphics.ClearRandomWriteTargets();
                Graphics.SetRandomWriteTarget(1, DeltaMultipleScatteringTexture);
                Graphics.SetRandomWriteTarget(2, m_scatteringLUT);
                m_precomputeMaterial.SetTexture("scattering_density_texture", DeltaScatteringDensityTexture);
                BlitWithDummy(m_precomputeMaterial, (int)PrecomputePass.MultipleScattering, SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT);
            }

            Graphics.ClearRandomWriteTargets();
            RenderTexture.active = null;

            // Release temporary textures.
            DeltaIrradianceTexture.Release();
            DeltaRayleighScatteringTexture.Release();
            DeltaMieScatteringTexture.Release();
            DeltaScatteringDensityTexture.Release();

            float timerEndCompute = Time.realtimeSinceStartup;
#if PHYSICAL_SKY_DEBUG
            Debug.Log("Computed atmospheric lookup textures in " + (timerEndCompute - timerStartCompute) * 1000.0f + "ms");
#endif
            m_needsRecompute = false;
        }

        public void ReleaseResources()
        {
#if PHYSICAL_SKY_DEBUG
            Debug.Log("Released Atmosphere Resources");
#endif
            ReleaseLookupTextures();
            m_needsRecompute = true;
        }

        public bool TexturesInvalid()
        {
            if (!m_transmittanceLUT || !m_transmittanceLUT.IsCreated())
                return true;
            else if (!m_scatteringLUT || !m_scatteringLUT.IsCreated())
                return true;

            else if (!m_irradianceLUT || !m_irradianceLUT.IsCreated())
                return true;
            else
                return false;
        }

        private void Awake()
        {
            m_needsRecompute = true;
        }

        private void OnDestroy()
        {
            ReleaseResources();
        }

#if UNITY_EDITOR
        /// <summary>
        /// Invoked by Unity whenever our values change in the inspector. Only
        /// called in the Editor
        /// </summary>
        protected void OnValidate()
        {
            ConstantSolarIrradiance = Mathf.Max(0.001f, ConstantSolarIrradiance);
            SunAngularRadius = Mathf.Max(0.001f, SunAngularRadius);
            PlanetaryRadius = Mathf.Max(1.0f, PlanetaryRadius);
            AtmosphereThickness = Mathf.Max(1.0f, AtmosphereThickness);
            Rayleigh = Mathf.Max(1.0e-22f, Rayleigh);
            RayleighScaleHeight = Mathf.Max(0.0f, RayleighScaleHeight);
            MieScaleHeight = Mathf.Max(0.0f, MieScaleHeight);
            MieAngstromAlpha = Mathf.Max(0.0f, MieAngstromAlpha);
            MieAngstromBeta = Mathf.Max(0.0f, MieAngstromBeta);
            MieSingleScatteringAlbedo = Mathf.Max(0.0f, MieSingleScatteringAlbedo);
            MiePhaseFunctionG = Mathf.Max(0.0f, MiePhaseFunctionG);
            GroundAlbedo = Mathf.Max(0.0f, GroundAlbedo);
            MaxSunZenithAngle = Mathf.Max(0.0f, MaxSunZenithAngle);
        }
#endif
    }
}