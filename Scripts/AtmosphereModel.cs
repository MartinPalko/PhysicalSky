using System.Collections.Generic;
using UnityEngine;
using System;

namespace PhysicalSky
{
    [CreateAssetMenu(fileName = "NewAtmosphereModel", menuName = "PhysicalSky/AtmosphereModel")]
    public class AtmosphereModel : ScriptableObject, IAtmosphereModel
    {
        // Matches that in PhysicalSkyCommon.cginc
        public struct DensityProfileLayer
        {
            float width;
            float exp_term;
            float exp_scale;
            float linear_term;
            float constant_term;

            public DensityProfileLayer(float width, float exp_term, float exp_scale, float linear_term, float constant_term)
            {
                this.width = width;
                this.exp_term = exp_term;
                this.exp_scale = exp_scale;
                this.linear_term = linear_term;
                this.constant_term = constant_term;
            }

            public static bool operator ==(DensityProfileLayer x, DensityProfileLayer y)
            {
                return x.width == y.width &&
                    x.exp_term == y.exp_term &&
                    x.exp_scale == y.exp_scale &&
                    x.linear_term == y.linear_term &&
                    x.constant_term == y.constant_term;
            }
            public static bool operator !=(DensityProfileLayer x, DensityProfileLayer y) { return !x.Equals(y); }
            public override bool Equals(object obj) { return obj is DensityProfileLayer && this == (DensityProfileLayer)obj; }
            public override int GetHashCode() { return width.GetHashCode() ^ exp_term.GetHashCode() ^ exp_scale.GetHashCode() ^ linear_term.GetHashCode() ^ constant_term.GetHashCode(); }

            public void SetInMaterial(string name, Material m)
            {
                m.SetFloatArray(name, new float[5] {
                    width / LENGTH_UNIT_IN_METERS,
                    exp_term,
                    exp_scale * LENGTH_UNIT_IN_METERS,
                    linear_term * LENGTH_UNIT_IN_METERS,
                    constant_term });
            }
        };

        public struct DensityProfile
        {
            DensityProfileLayer layer0;
            DensityProfileLayer layer1;

            public DensityProfile(DensityProfileLayer layer0, DensityProfileLayer layer1)
            {
                this.layer0 = layer0;
                this.layer1 = layer1;
            }

            public DensityProfile(DensityProfileLayer layer)
            {
                this.layer0 = layer;
                this.layer1 = layer;
            }

            public static bool operator ==(DensityProfile x, DensityProfile y)
            {
                return x.layer0 == y.layer0 && x.layer1 == y.layer1;
            }
            public static bool operator !=(DensityProfile x, DensityProfile y) { return !x.Equals(y); }
            public override bool Equals(object obj) { return obj is DensityProfile && this == (DensityProfile)obj; }
            public override int GetHashCode() { return layer0.GetHashCode() ^ layer1.GetHashCode(); }

            public void SetInMaterial(string name, Material m)
            {
                layer0.SetInMaterial(name + "0", m);
                layer1.SetInMaterial(name + "1", m);
            }
        }

        enum PrecomputePass
        {
            Transmittance = 0,
            DirectIrradiance = 1,
            SingleScattering = 2,
            ScatteringDensity = 3,
            ComputeIndirectIrradiance = 4,
            AccumulateIndirectIrradiance = 5,
            ComputeMultipleScattering = 6,
            AccumulateMultipleScattering = 7
        }

        private bool m_needsRecompute = true;
        public bool NeedsRecompute { get { return m_needsRecompute || TexturesInvalid(); } }

        [SerializeField]
        int m_scatteringOrders = 1;
        public int ScatteringOrders
        {
            get { return m_scatteringOrders; }
            set
            {
                if (m_scatteringOrders != value)
                {
                    m_scatteringOrders = value;
                    m_needsRecompute = true;
                }
            }
        }

        [SerializeField]
        private bool m_useOzone = true;
        public bool UseOzone
        {
            get { return m_useOzone; }
            set
            {
                if (m_useOzone != value)
                {
                    m_useOzone = value;
                    m_needsRecompute = true;
                }
            }
        }

        [SerializeField]
        private bool m_useConstantSolarSpectrum = false;
        public bool UseConstantSolarSpectrum
        {
            get { return m_useConstantSolarSpectrum; }
            set
            {
                if (m_useConstantSolarSpectrum != value)
                {
                    m_useConstantSolarSpectrum = value;
                    m_needsRecompute = true;
                }
            }
        }

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
        const float LAMBDA_R = 680.0f;
        const float LAMBDA_G = 550.0f;
        const float LAMBDA_B = 440.0f;

        const RenderTextureFormat LUT_FORMAT = RenderTextureFormat.ARGBFloat; // TODO: Try changing to half or float

        const TextureWrapMode LUT_WRAP_MODE = TextureWrapMode.Clamp;
        const FilterMode LUT_FILTER_MODE = FilterMode.Trilinear;

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

        const float LENGTH_UNIT_IN_METERS = 1000.0f;

        // (see http://rredc.nrel.gov/solar/spectra/am1.5/ASTMG173/ASTMG173.html),
        // summed and averaged in each bin (e.g. the value for 360nm is the average
        // of the ASTM G-173 values for all wavelengths between 360 and 370nm).
        // Values in W.m^-2.
        const int LAMBDA_MIN = 360;
        const int LAMBDA_MAX = 830;
        static readonly float[] SOLAR_IRRADIANCE = {
    1.11776f, 1.14259f, 1.01249f, 1.14716f, 1.72765f, 1.73054f, 1.6887f, 1.61253f,
    1.91198f, 2.03474f, 2.02042f, 2.02212f, 1.93377f, 1.95809f, 1.91686f, 1.8298f,
    1.8685f, 1.8931f, 1.85149f, 1.8504f, 1.8341f, 1.8345f, 1.8147f, 1.78158f, 1.7533f,
    1.6965f, 1.68194f, 1.64654f, 1.6048f, 1.52143f, 1.55622f, 1.5113f, 1.474f, 1.4482f,
    1.41018f, 1.36775f, 1.34188f, 1.31429f, 1.28303f, 1.26758f, 1.2367f, 1.2082f,
    1.18737f, 1.14683f, 1.12362f, 1.1058f, 1.07124f, 1.04992f
        };

        // Values from http://www.iup.uni-bremen.de/gruppen/molspec/databases/
        // referencespectra/o3spectra2011/index.html for 233K, summed and averaged in
        // each bin (e.g. the value for 360nm is the average of the original values
        // for all wavelengths between 360 and 370nm). Values in m^2.
        static readonly float[] OZONE_CROSS_SECTION = {
    1.18e-27f, 2.182e-28f, 2.818e-28f, 6.636e-28f, 1.527e-27f, 2.763e-27f, 5.52e-27f,
    8.451e-27f, 1.582e-26f, 2.316e-26f, 3.669e-26f, 4.924e-26f, 7.752e-26f, 9.016e-26f,
    1.48e-25f, 1.602e-25f, 2.139e-25f, 2.755e-25f, 3.091e-25f, 3.5e-25f, 4.266e-25f,
    4.672e-25f, 4.398e-25f, 4.701e-25f, 5.019e-25f, 4.305e-25f, 3.74e-25f, 3.215e-25f,
    2.662e-25f, 2.238e-25f, 1.852e-25f, 1.473e-25f, 1.209e-25f, 9.423e-26f, 7.455e-26f,
    6.566e-26f, 5.105e-26f, 4.15e-26f, 4.228e-26f, 3.237e-26f, 2.451e-26f, 2.801e-26f,
    2.534e-26f, 1.624e-26f, 1.465e-26f, 2.078e-26f, 1.383e-26f, 7.105e-27f
        };

        // From https://en.wikipedia.org/wiki/Dobson_unit, in molecules.m^-2.
        const double kDobsonUnit = 2.687e20f;
        // Maximum number density of ozone molecules, in m^-3 (computed so at to get
        // 300 Dobson units of ozone - for this we divide 300 DU by the integral of
        // the ozone density profile defined below, which is equal to 15km).
        const float kMaxOzoneNumberDensity = (float)(300.0 * kDobsonUnit / 15000.0);


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
        [NonSerialized]
        private List<float> m_absorptionExtinction = new List<float>();
        [NonSerialized]
        private DensityProfile m_rayleighDensity;
        [NonSerialized]
        private DensityProfile m_mieDensity;
        [NonSerialized]
        private DensityProfile m_absorptionDensity;

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
            Vector4 ret = new Vector4(
                (float)(Interpolate(m_wavelengths, v, LAMBDA_R) * scale),
                (float)(Interpolate(m_wavelengths, v, LAMBDA_G) * scale),
                (float)(Interpolate(m_wavelengths, v, LAMBDA_B) * scale),
                1.0f);
            return ret;

            //return new Vector4(
            //    (float)(Interpolate(m_wavelengths, v, LAMBDA_R) * scale),
            //    (float)(Interpolate(m_wavelengths, v, LAMBDA_G) * scale),
            //    (float)(Interpolate(m_wavelengths, v, LAMBDA_B) * scale),
            //    1.0f);
        }

        private void AllocateLookupTextures()
        {
            if (!m_transmittanceLUT)
                m_transmittanceLUT = new RenderTexture(TRANSMITTANCE_TEXTURE_WIDTH, TRANSMITTANCE_TEXTURE_HEIGHT, 0, LUT_FORMAT, RenderTextureReadWrite.Linear);
            if (!m_transmittanceLUT.IsCreated())
            {
                m_transmittanceLUT.useMipMap = false;
                m_transmittanceLUT.wrapMode = LUT_WRAP_MODE;
                m_transmittanceLUT.filterMode = LUT_FILTER_MODE;
                m_transmittanceLUT.Create();
            }

            if (!m_scatteringLUT)
                m_scatteringLUT = new RenderTexture(SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT, 0, LUT_FORMAT, RenderTextureReadWrite.Linear);
            if (!m_scatteringLUT.IsCreated())
            {
                m_scatteringLUT.volumeDepth = SCATTERING_TEXTURE_DEPTH;
                m_scatteringLUT.useMipMap = false;
                m_scatteringLUT.dimension = UnityEngine.Rendering.TextureDimension.Tex3D;
                m_scatteringLUT.wrapMode = LUT_WRAP_MODE;
                m_scatteringLUT.filterMode = LUT_FILTER_MODE;
                m_scatteringLUT.Create();
            }

            if (!m_irradianceLUT)
                m_irradianceLUT = new RenderTexture(IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT, 0, LUT_FORMAT, RenderTextureReadWrite.Linear);
            if (!m_irradianceLUT.IsCreated())
            {
                m_irradianceLUT.useMipMap = false;
                m_irradianceLUT.wrapMode = LUT_WRAP_MODE;
                m_irradianceLUT.filterMode = LUT_FILTER_MODE;
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

        public void SetShaderUniforms(Material m)
        {
            m.SetVector("_solar_irradiance", ScaleToWavelengths(m_solarIrradiance, 1.0f));
            m.SetFloat("_sun_angular_radius", m_sunAngularRadius);
            m.SetFloat("_bottom_radius", m_planetaryRadius / LENGTH_UNIT_IN_METERS);
            m.SetFloat("_top_radius", (m_planetaryRadius + m_atmosphereThickness) / LENGTH_UNIT_IN_METERS);
            m_rayleighDensity.SetInMaterial("_rayleigh_density", m);
            m.SetVector("_rayleigh_scattering", ScaleToWavelengths(m_rayleighScattering, LENGTH_UNIT_IN_METERS));
            m_mieDensity.SetInMaterial("_mie_density", m);
            m.SetVector("_mie_scattering", ScaleToWavelengths(m_mieScattering, LENGTH_UNIT_IN_METERS));
            m.SetVector("_mie_extinction", ScaleToWavelengths(m_mieExtinction, LENGTH_UNIT_IN_METERS));
            m.SetFloat("_mie_phase_function_g", m_miePhaseFunctionG);
            m_absorptionDensity.SetInMaterial("_absorption_density", m);
            m.SetVector("_absorption_extinction", ScaleToWavelengths(m_absorptionExtinction, LENGTH_UNIT_IN_METERS));
            m.SetVector("_ground_albedo", new Vector3(m_groundAlbedo, m_groundAlbedo, m_groundAlbedo));
            m.SetFloat("_mu_s_min", Mathf.Cos(m_maxSunZenithAngle));
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
            m_absorptionExtinction.Clear();

            m_sunSolidAngle = Mathf.PI * m_sunAngularRadius * m_sunAngularRadius;

            for (int l = LAMBDA_MIN; l <= LAMBDA_MAX; l += 10)
            {
                float lambda = l * 1e-3f;  // micro-meters
                float mie = m_mieAngstromBeta / m_mieScaleHeight * Mathf.Pow(lambda, -m_mieAngstromAlpha);
                m_wavelengths.Add(l);

                if (m_useConstantSolarSpectrum)
                    m_solarIrradiance.Add(m_constantSolarIrradiance);
                else
                    m_solarIrradiance.Add(SOLAR_IRRADIANCE[(l - LAMBDA_MIN) / 10]);

                m_rayleighScattering.Add(m_rayleigh * Mathf.Pow(lambda, -4));
                m_mieScattering.Add(mie * m_mieSingleScatteringAlbedo);
                m_mieExtinction.Add(mie);
                m_absorptionExtinction.Add(m_useOzone ? kMaxOzoneNumberDensity * OZONE_CROSS_SECTION[(l - LAMBDA_MIN) / 10] : 0.0f);
            }

            m_rayleighDensity = new DensityProfile(new DensityProfileLayer(0.0f, 1.0f, -1.0f / m_rayleighScaleHeight, 0.0f, 0.0f));
            m_mieDensity = new DensityProfile(new DensityProfileLayer(0.0f, 1.0f, -1.0f / m_mieScaleHeight, 0.0f, 0.0f));

            // Density profile increasing linearly from 0 to 1 between 10 and 25km, and
            // decreasing linearly from 1 to 0 between 25 and 40km. This is an approximate
            // profile from http://www.kln.ac.lk/science/Chemistry/Teaching_Resources/
            // Documents/Introduction%20to%20atmospheric%20chemistry.pdf (page 10).
            // (ozone density)
            m_absorptionDensity = new DensityProfile(
                new DensityProfileLayer(25000.0f, 0.0f, 0.0f, 1.0f / 15000.0f, -2.0f / 3.0f),
                new DensityProfileLayer(0.0f, 0.0f, 0.0f, -1.0f / 15000.0f, 8.0f / 3.0f));

            if (!m_precomputeMaterial)
                m_precomputeMaterial = new Material(PrecomputeShader);

            SetShaderUniforms(m_precomputeMaterial);

            RenderTexture DeltaIrradianceTexture = new RenderTexture(IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT, 0, LUT_FORMAT, RenderTextureReadWrite.Linear);
            DeltaIrradianceTexture.useMipMap = false;
            DeltaIrradianceTexture.wrapMode = LUT_WRAP_MODE;
            DeltaIrradianceTexture.filterMode = LUT_FILTER_MODE;
            DeltaIrradianceTexture.Create();

            RenderTexture DeltaRayleighScatteringTexture = new RenderTexture(SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT, 0, LUT_FORMAT, RenderTextureReadWrite.Linear);
            DeltaRayleighScatteringTexture.volumeDepth = SCATTERING_TEXTURE_DEPTH;
            DeltaRayleighScatteringTexture.useMipMap = false;
            DeltaRayleighScatteringTexture.dimension = UnityEngine.Rendering.TextureDimension.Tex3D;
            DeltaRayleighScatteringTexture.wrapMode = LUT_WRAP_MODE;
            DeltaRayleighScatteringTexture.filterMode = LUT_FILTER_MODE;
            DeltaRayleighScatteringTexture.Create();

            RenderTexture DeltaMieScatteringTexture = new RenderTexture(SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT, 0, LUT_FORMAT, RenderTextureReadWrite.Linear);
            DeltaMieScatteringTexture.volumeDepth = SCATTERING_TEXTURE_DEPTH;
            DeltaMieScatteringTexture.useMipMap = false;
            DeltaMieScatteringTexture.dimension = UnityEngine.Rendering.TextureDimension.Tex3D;
            DeltaMieScatteringTexture.wrapMode = LUT_WRAP_MODE;
            DeltaMieScatteringTexture.filterMode = LUT_FILTER_MODE;
            DeltaMieScatteringTexture.Create();

            RenderTexture DeltaScatteringDensityTexture = new RenderTexture(SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT, 0, LUT_FORMAT, RenderTextureReadWrite.Linear);
            DeltaScatteringDensityTexture.volumeDepth = SCATTERING_TEXTURE_DEPTH;
            DeltaScatteringDensityTexture.useMipMap = false;
            DeltaScatteringDensityTexture.dimension = UnityEngine.Rendering.TextureDimension.Tex3D;
            DeltaScatteringDensityTexture.wrapMode = LUT_WRAP_MODE;
            DeltaScatteringDensityTexture.filterMode = LUT_FILTER_MODE;
            DeltaScatteringDensityTexture.Create();

            // DeltaMultipleScatteringTexture is only needed to compute scattering
            // order 3 or more, while DeltaRayleighScatteringTexture and
            // DeltaMieScatteringTexture are only needed to compute double scattering.
            // Therefore, to save memory, we can store DeltaRayleighScatteringTexture
            // and DeltaMultipleScatteringTexture in the same GPU texture.
            RenderTexture DeltaMultipleScatteringTexture = DeltaRayleighScatteringTexture;

            // Allocate final lookup textures
            AllocateLookupTextures();

            // Compute Transmittance LUT.
            Utilities.GraphicsHelpers.Blit(m_transmittanceLUT, m_precomputeMaterial, (int)PrecomputePass.Transmittance);
            m_precomputeMaterial.SetTexture("transmittance_texture", m_transmittanceLUT); // Set for subsequent shaders to read

            // Compute Direct Irradiance into DeltaIrradianceTexture and Initialize m_irradianceLUT with 0 (we don't want the direct irradiance in m_irradianceLUT, but only the irradiance from the sky)
            Utilities.GraphicsHelpers.Blit3D(new RenderTexture[2] { DeltaIrradianceTexture, m_irradianceLUT }, m_precomputeMaterial, (int)PrecomputePass.DirectIrradiance);

            // Compute Raylie and Mie Single Scattering, and store them in DeltaRayleighScatteringTexture, DeltaMieScatteringTexture as well as scatteringLUT
            Utilities.GraphicsHelpers.Blit3D(new RenderTexture[3] { DeltaRayleighScatteringTexture, DeltaMieScatteringTexture, m_scatteringLUT}, m_precomputeMaterial, (int)PrecomputePass.SingleScattering);

            // Compute to the nth order of scattering, in sequence
            for (int scatteringOrder = 2; scatteringOrder <= m_scatteringOrders; ++scatteringOrder)
            {
                m_precomputeMaterial.SetTexture("transmittance_texture", m_transmittanceLUT);
                m_precomputeMaterial.SetTexture("single_rayleigh_scattering_texture", DeltaRayleighScatteringTexture);
                m_precomputeMaterial.SetTexture("single_mie_scattering_texture", DeltaMieScatteringTexture);
                m_precomputeMaterial.SetTexture("multiple_scattering_texture", DeltaMultipleScatteringTexture);
                m_precomputeMaterial.SetTexture("irradiance_texture", DeltaIrradianceTexture);
                m_precomputeMaterial.SetTexture("scattering_density_texture", DeltaScatteringDensityTexture);
                m_precomputeMaterial.SetInt("scattering_order", scatteringOrder);

                // Compute the scattering density, and store it in DeltaScatteringDensityTexture.
                Utilities.GraphicsHelpers.Blit3D(DeltaScatteringDensityTexture, m_precomputeMaterial, (int)PrecomputePass.ScatteringDensity);

                // Compute the indirect irradiance, store it in DeltaIrradianceTexture and accumulate it in irradianceLUT.
                m_precomputeMaterial.SetInt("scattering_order", scatteringOrder - 1);
                Utilities.GraphicsHelpers.Blit3D(DeltaIrradianceTexture, m_precomputeMaterial, (int)PrecomputePass.ComputeIndirectIrradiance);
                Utilities.GraphicsHelpers.Blit3D(m_irradianceLUT, m_precomputeMaterial, (int)PrecomputePass.AccumulateIndirectIrradiance);

                // Compute the multiple scattering, store it in DeltaMultipleScatteringTexture, and accumulate it in scatteringLUT.
                Utilities.GraphicsHelpers.Blit3D(DeltaMultipleScatteringTexture, m_precomputeMaterial, (int)PrecomputePass.ComputeMultipleScattering);
                Utilities.GraphicsHelpers.Blit3D(m_scatteringLUT, m_precomputeMaterial, (int)PrecomputePass.AccumulateMultipleScattering);
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
            ScatteringOrders = Mathf.Clamp(ScatteringOrders, 0, 10);
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