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

        private enum PrecomputePass
        {
            ComputeTransmittance = 0,
            ComputeDirectIrradiance = 1,
            ComputeSingleScattering = 2,
            AccumulateSingleScattering = 3,
            ComputeScatteringDensity = 4,
            ComputeIndirectIrradiance = 5,
            AccumulateIndirectIrradiance = 6,
            ComputeMultipleScattering = 7,
            AccumulateMultipleScattering = 8
        }

        public enum LuminanceType
        {
            // Render the spectral radiance at kLambdaR, kLambdaG, kLambdaB.
            none = 0,
            // Render the sRGB luminance, using an approximate (on the fly) conversion
            // from 3 spectral radiance values only (see section 14.3 in "A Qualitative
            // and Quantitative Evaluation of 8 Clear Sky Models"
            // https://arxiv.org/pdf/1612.04336.pdf).
            approximate = 1,
            // Render the sRGB luminance, precomputed from 15 spectral radiance values
            // (see section 4.4 in "Real-time Spectral Scattering in Large-scale 
            // Natural Participating Media"
            // http://www.oskee.wz.cz/stranka/uploads/SCCG10ElekKmoch.pdf).
            precomputed = 2
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
        LuminanceType m_luminance = LuminanceType.precomputed;
        public LuminanceType Luminance
        {
            get { return m_luminance; }
            set
            {
                if (m_luminance != value)
                {
                    m_luminance = value;
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

        // The conversion factor between watts and lumens.
        const float MAX_LUMINOUS_EFFICACY = 683.0f;

        const float LAMBDA_R = 680.0f;
        const float LAMBDA_G = 550.0f;
        const float LAMBDA_B = 440.0f;

        // Values from "CIE (1931) 2-deg color matching functions", see
        // "http://web.archive.org/web/20081228084047/
        //    http://www.cvrl.org/database/data/cmfs/ciexyz31.txt".
        static readonly double[] CIE_2_DEG_COLOR_MATCHING_FUNCTIONS = {
    360, 0.000129900000, 0.000003917000, 0.000606100000,
    365, 0.000232100000, 0.000006965000, 0.001086000000,
    370, 0.000414900000, 0.000012390000, 0.001946000000,
    375, 0.000741600000, 0.000022020000, 0.003486000000,
    380, 0.001368000000, 0.000039000000, 0.006450001000,
    385, 0.002236000000, 0.000064000000, 0.010549990000,
    390, 0.004243000000, 0.000120000000, 0.020050010000,
    395, 0.007650000000, 0.000217000000, 0.036210000000,
    400, 0.014310000000, 0.000396000000, 0.067850010000,
    405, 0.023190000000, 0.000640000000, 0.110200000000,
    410, 0.043510000000, 0.001210000000, 0.207400000000,
    415, 0.077630000000, 0.002180000000, 0.371300000000,
    420, 0.134380000000, 0.004000000000, 0.645600000000,
    425, 0.214770000000, 0.007300000000, 1.039050100000,
    430, 0.283900000000, 0.011600000000, 1.385600000000,
    435, 0.328500000000, 0.016840000000, 1.622960000000,
    440, 0.348280000000, 0.023000000000, 1.747060000000,
    445, 0.348060000000, 0.029800000000, 1.782600000000,
    450, 0.336200000000, 0.038000000000, 1.772110000000,
    455, 0.318700000000, 0.048000000000, 1.744100000000,
    460, 0.290800000000, 0.060000000000, 1.669200000000,
    465, 0.251100000000, 0.073900000000, 1.528100000000,
    470, 0.195360000000, 0.090980000000, 1.287640000000,
    475, 0.142100000000, 0.112600000000, 1.041900000000,
    480, 0.095640000000, 0.139020000000, 0.812950100000,
    485, 0.057950010000, 0.169300000000, 0.616200000000,
    490, 0.032010000000, 0.208020000000, 0.465180000000,
    495, 0.014700000000, 0.258600000000, 0.353300000000,
    500, 0.004900000000, 0.323000000000, 0.272000000000,
    505, 0.002400000000, 0.407300000000, 0.212300000000,
    510, 0.009300000000, 0.503000000000, 0.158200000000,
    515, 0.029100000000, 0.608200000000, 0.111700000000,
    520, 0.063270000000, 0.710000000000, 0.078249990000,
    525, 0.109600000000, 0.793200000000, 0.057250010000,
    530, 0.165500000000, 0.862000000000, 0.042160000000,
    535, 0.225749900000, 0.914850100000, 0.029840000000,
    540, 0.290400000000, 0.954000000000, 0.020300000000,
    545, 0.359700000000, 0.980300000000, 0.013400000000,
    550, 0.433449900000, 0.994950100000, 0.008749999000,
    555, 0.512050100000, 1.000000000000, 0.005749999000,
    560, 0.594500000000, 0.995000000000, 0.003900000000,
    565, 0.678400000000, 0.978600000000, 0.002749999000,
    570, 0.762100000000, 0.952000000000, 0.002100000000,
    575, 0.842500000000, 0.915400000000, 0.001800000000,
    580, 0.916300000000, 0.870000000000, 0.001650001000,
    585, 0.978600000000, 0.816300000000, 0.001400000000,
    590, 1.026300000000, 0.757000000000, 0.001100000000,
    595, 1.056700000000, 0.694900000000, 0.001000000000,
    600, 1.062200000000, 0.631000000000, 0.000800000000,
    605, 1.045600000000, 0.566800000000, 0.000600000000,
    610, 1.002600000000, 0.503000000000, 0.000340000000,
    615, 0.938400000000, 0.441200000000, 0.000240000000,
    620, 0.854449900000, 0.381000000000, 0.000190000000,
    625, 0.751400000000, 0.321000000000, 0.000100000000,
    630, 0.642400000000, 0.265000000000, 0.000049999990,
    635, 0.541900000000, 0.217000000000, 0.000030000000,
    640, 0.447900000000, 0.175000000000, 0.000020000000,
    645, 0.360800000000, 0.138200000000, 0.000010000000,
    650, 0.283500000000, 0.107000000000, 0.000000000000,
    655, 0.218700000000, 0.081600000000, 0.000000000000,
    660, 0.164900000000, 0.061000000000, 0.000000000000,
    665, 0.121200000000, 0.044580000000, 0.000000000000,
    670, 0.087400000000, 0.032000000000, 0.000000000000,
    675, 0.063600000000, 0.023200000000, 0.000000000000,
    680, 0.046770000000, 0.017000000000, 0.000000000000,
    685, 0.032900000000, 0.011920000000, 0.000000000000,
    690, 0.022700000000, 0.008210000000, 0.000000000000,
    695, 0.015840000000, 0.005723000000, 0.000000000000,
    700, 0.011359160000, 0.004102000000, 0.000000000000,
    705, 0.008110916000, 0.002929000000, 0.000000000000,
    710, 0.005790346000, 0.002091000000, 0.000000000000,
    715, 0.004109457000, 0.001484000000, 0.000000000000,
    720, 0.002899327000, 0.001047000000, 0.000000000000,
    725, 0.002049190000, 0.000740000000, 0.000000000000,
    730, 0.001439971000, 0.000520000000, 0.000000000000,
    735, 0.000999949300, 0.000361100000, 0.000000000000,
    740, 0.000690078600, 0.000249200000, 0.000000000000,
    745, 0.000476021300, 0.000171900000, 0.000000000000,
    750, 0.000332301100, 0.000120000000, 0.000000000000,
    755, 0.000234826100, 0.000084800000, 0.000000000000,
    760, 0.000166150500, 0.000060000000, 0.000000000000,
    765, 0.000117413000, 0.000042400000, 0.000000000000,
    770, 0.000083075270, 0.000030000000, 0.000000000000,
    775, 0.000058706520, 0.000021200000, 0.000000000000,
    780, 0.000041509940, 0.000014990000, 0.000000000000,
    785, 0.000029353260, 0.000010600000, 0.000000000000,
    790, 0.000020673830, 0.000007465700, 0.000000000000,
    795, 0.000014559770, 0.000005257800, 0.000000000000,
    800, 0.000010253980, 0.000003702900, 0.000000000000,
    805, 0.000007221456, 0.000002607800, 0.000000000000,
    810, 0.000005085868, 0.000001836600, 0.000000000000,
    815, 0.000003581652, 0.000001293400, 0.000000000000,
    820, 0.000002522525, 0.000000910930, 0.000000000000,
    825, 0.000001776509, 0.000000641530, 0.000000000000,
    830, 0.000001251141, 0.000000451810, 0.000000000000,
        };

        // The conversion matrix from XYZ to linear sRGB color spaces.
        // Values from https://en.wikipedia.org/wiki/SRGB.
        static readonly float[] XYZ_TO_SRGB = {
    +3.2406f, -1.5372f, -0.4986f,
    -0.9689f, +1.8758f, +0.0415f,
    +0.0557f, -0.2040f, +1.0570f
        };

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
        Vector3 m_sky_k;
        [NonSerialized]
        Vector3 m_sun_k;
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

        private float CieColorMatchingFunctionTableValue(float wavelength, int column)
        {
            if (wavelength <= LAMBDA_MIN || wavelength >= LAMBDA_MAX)
            {
                return 0.0f;
            }
            float u = (wavelength - LAMBDA_MIN) / 5.0f;
            int row = (int)Mathf.Floor(u);
            Debug.Assert(row >= 0 && row + 1 < 95);
            Debug.Assert(CIE_2_DEG_COLOR_MATCHING_FUNCTIONS[4 * row] <= wavelength &&
                   CIE_2_DEG_COLOR_MATCHING_FUNCTIONS[4 * (row + 1)] >= wavelength);
            u -= row;
            double result = CIE_2_DEG_COLOR_MATCHING_FUNCTIONS[4 * row + column] * (1.0 - u) + CIE_2_DEG_COLOR_MATCHING_FUNCTIONS[4 * (row + 1) + column] * u;
            return (float)result;
        }

        // The returned constants are in lumen.nm / watt.
        private Vector3 ComputeSpectralRadianceToLuminanceFactors(
            List<float> wavelengths,
            List<float> solar_irradiance,
            float lambda_power)
        {
            Vector3 k = Vector3.zero;
            float solar_r = Interpolate(wavelengths, solar_irradiance, LAMBDA_R);
            float solar_g = Interpolate(wavelengths, solar_irradiance, LAMBDA_G);
            float solar_b = Interpolate(wavelengths, solar_irradiance, LAMBDA_B);
            int dlambda = 1;
            for (int lambda = LAMBDA_MIN; lambda < LAMBDA_MAX; lambda += dlambda)
            {
                float x_bar = CieColorMatchingFunctionTableValue(lambda, 1);
                float y_bar = CieColorMatchingFunctionTableValue(lambda, 2);
                float z_bar = CieColorMatchingFunctionTableValue(lambda, 3);
                float r_bar = XYZ_TO_SRGB[0] * x_bar + XYZ_TO_SRGB[1] * y_bar + XYZ_TO_SRGB[2] * z_bar;
                float g_bar = XYZ_TO_SRGB[3] * x_bar + XYZ_TO_SRGB[4] * y_bar + XYZ_TO_SRGB[5] * z_bar;
                float b_bar = XYZ_TO_SRGB[6] * x_bar + XYZ_TO_SRGB[7] * y_bar + XYZ_TO_SRGB[8] * z_bar;
                float irradiance = Interpolate(wavelengths, solar_irradiance, lambda);
                k.x += r_bar * irradiance / solar_r * 
                    Mathf.Pow(lambda / LAMBDA_R, lambda_power);
                k.y += g_bar * irradiance / solar_g *
                    Mathf.Pow(lambda / LAMBDA_G, lambda_power);
                k.z += b_bar * irradiance / solar_b *
                    Mathf.Pow(lambda / LAMBDA_B, lambda_power);
            }
            k.x *= MAX_LUMINOUS_EFFICACY * dlambda;
            k.y *= MAX_LUMINOUS_EFFICACY * dlambda;
            k.z *= MAX_LUMINOUS_EFFICACY * dlambda;
            return k;
        }

        private Vector4 ScaleToWavelengths(List<float> v, Vector3 lambdas, float scale)
        {
            return new Vector4(
                (Interpolate(m_wavelengths, v, lambdas.x) * scale),
                (Interpolate(m_wavelengths, v, lambdas.y) * scale),
                (Interpolate(m_wavelengths, v, lambdas.z) * scale),
                1.0f);
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
            SetShaderUniforms(m, new Vector3(LAMBDA_R, LAMBDA_G, LAMBDA_B));
        }

        private void SetShaderUniforms(Material m, Vector3 lambdas)
        {
            m.SetInt("_use_luminance", Luminance == LuminanceType.none ? 0 : 1);
            m.SetVector("_sky_spectral_radiance_to_luminance", Luminance != LuminanceType.none ? m_sky_k : Vector3.one);
            m.SetVector("_sun_spectral_radiance_to_luminance", Luminance != LuminanceType.none ? m_sun_k : Vector3.one);

            m.SetVector("_solar_irradiance", ScaleToWavelengths(m_solarIrradiance, lambdas, 1.0f));
            m.SetFloat("_sun_angular_radius", m_sunAngularRadius);
            m.SetFloat("_bottom_radius", m_planetaryRadius / LENGTH_UNIT_IN_METERS);
            m.SetFloat("_top_radius", (m_planetaryRadius + m_atmosphereThickness) / LENGTH_UNIT_IN_METERS);
            m_rayleighDensity.SetInMaterial("_rayleigh_density", m);
            m.SetVector("_rayleigh_scattering", ScaleToWavelengths(m_rayleighScattering, lambdas, LENGTH_UNIT_IN_METERS));
            m_mieDensity.SetInMaterial("_mie_density", m);
            m.SetVector("_mie_scattering", ScaleToWavelengths(m_mieScattering, lambdas, LENGTH_UNIT_IN_METERS));
            m.SetVector("_mie_extinction", ScaleToWavelengths(m_mieExtinction, lambdas, LENGTH_UNIT_IN_METERS));
            m.SetFloat("_mie_phase_function_g", m_miePhaseFunctionG);
            m_absorptionDensity.SetInMaterial("_absorption_density", m);
            m.SetVector("_absorption_extinction", ScaleToWavelengths(m_absorptionExtinction, lambdas, LENGTH_UNIT_IN_METERS));
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
                string currentLevelString = (SystemInfo.graphicsShaderLevel / 10).ToString() + "." + (SystemInfo.graphicsShaderLevel % 10).ToString();
                Debug.LogError("Computing atmospheric lookup textures requires shader model 5.0 or higher! Current is " + currentLevelString);
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

            int numPrecomputedWavelengths = m_luminance == LuminanceType.precomputed ? 15 : 3;
            bool precomputeIlluminance = numPrecomputedWavelengths > 3;

            // Compute the values for sky_spectral_radiance_to_luminance. In theory
            // this should be 1 in precomputed illuminance mode (because the precomputed
            // textures already contain illuminance values). In practice, however, storing
            // true illuminance values in half precision textures yields artefacts
            // (because the values are too large), so we store illuminance values divided
            // by MAX_LUMINOUS_EFFICACY instead. This is why, in precomputed illuminance
            // mode, we set sky_spectral_radiance_to_luminance to MAX_LUMINOUS_EFFICACY.
            if (precomputeIlluminance)
                m_sky_k = new Vector3(MAX_LUMINOUS_EFFICACY, MAX_LUMINOUS_EFFICACY, MAX_LUMINOUS_EFFICACY);
            else
                m_sky_k = ComputeSpectralRadianceToLuminanceFactors(m_wavelengths, m_solarIrradiance, -3 /* lambda_power */);

            m_sun_k = ComputeSpectralRadianceToLuminanceFactors(m_wavelengths, m_solarIrradiance, 0 /* lambda_power */);

            if (numPrecomputedWavelengths <= 3)
            {
                Vector3 lambdas = new Vector3(LAMBDA_R, LAMBDA_G, LAMBDA_B);
                float[] luminanceFromRadiance = new float[9]{
                    1.0f, 0.0f, 0.0f,
                    0.0f, 1.0f, 0.0f,
                    0.0f, 0.0f, 1.0f};

                Compute(lambdas, luminanceFromRadiance, false);
            }
            else
            {
                int num_iterations = (numPrecomputedWavelengths + 2) / 3;
                float dlambda = (LAMBDA_MAX - LAMBDA_MIN) / (3 * num_iterations);
                for (int i = 0; i < num_iterations; ++i)
                {
                    Vector3 lambdas = new Vector3(
                        LAMBDA_MIN + (3 * i + 0.5f) * dlambda,
                        LAMBDA_MIN + (3 * i + 1.5f) * dlambda,
                        LAMBDA_MIN + (3 * i + 2.5f) * dlambda);

                    Func<float, int, float> coeff = (lambda, component) =>
                    {
                        // Note that we don't include MAX_LUMINOUS_EFFICACY here, to avoid
                        // artefacts due to too large values when using half precision on GPU.
                        // We add this term back in kAtmosphereShader, via
                        // SKY_SPECTRAL_RADIANCE_TO_LUMINANCE (see also the comments in the
                        // Model constructor).
                        double x = CieColorMatchingFunctionTableValue(lambda, 1);
                        double y = CieColorMatchingFunctionTableValue(lambda, 2);
                        double z = CieColorMatchingFunctionTableValue(lambda, 3);
                        return (float)((
                            XYZ_TO_SRGB[component * 3] * x +
                            XYZ_TO_SRGB[component * 3 + 1] * y +
                            XYZ_TO_SRGB[component * 3 + 2] * z) * dlambda);
                    };

                    float[] luminanceFromRadiance = new float[9]{
                        coeff(lambdas[0], 0), coeff(lambdas[1], 0), coeff(lambdas[2], 0),
                        coeff(lambdas[0], 1), coeff(lambdas[1], 1), coeff(lambdas[2], 1),
                        coeff(lambdas[0], 2), coeff(lambdas[1], 2), coeff(lambdas[2], 2)
                    };
                    Compute(lambdas, luminanceFromRadiance, i == 0 ? false : true);
                }

                // After the above iterations, the transmittance texture contains the
                // transmittance for the 3 wavelengths used at the last iteration. But we
                // want the transmittance at kLambdaR, kLambdaG, kLambdaB instead, so we
                // must recompute it here for these 3 wavelengths:
                SetShaderUniforms(m_precomputeMaterial, new Vector3(LAMBDA_R, LAMBDA_G, LAMBDA_B));
                Utilities.GraphicsHelpers.Blit(m_transmittanceLUT, m_precomputeMaterial, (int)PrecomputePass.ComputeTransmittance);
            }

            float timerEndCompute = Time.realtimeSinceStartup;
#if PHYSICAL_SKY_DEBUG
            Debug.Log("Computed atmospheric lookup textures in " + (timerEndCompute - timerStartCompute) * 1000.0f + "ms");
#endif
            m_needsRecompute = false;
        }

        private void Compute(Vector3 lambdas, float[] luminanceFromRadiance, bool accumulate)
        {
            Debug.Assert(luminanceFromRadiance.Length == 9);

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

            SetShaderUniforms(m_precomputeMaterial, lambdas);
            m_precomputeMaterial.SetFloatArray("_luminance_from_radiance", luminanceFromRadiance);
            m_precomputeMaterial.SetTexture("transmittance_texture", m_transmittanceLUT);
            m_precomputeMaterial.SetTexture("single_rayleigh_scattering_texture", DeltaRayleighScatteringTexture);
            m_precomputeMaterial.SetTexture("single_mie_scattering_texture", DeltaMieScatteringTexture);
            m_precomputeMaterial.SetTexture("multiple_scattering_texture", DeltaMultipleScatteringTexture);
            m_precomputeMaterial.SetTexture("irradiance_texture", DeltaIrradianceTexture);
            m_precomputeMaterial.SetTexture("scattering_density_texture", DeltaScatteringDensityTexture);

            // Compute Transmittance LUT.
            Utilities.GraphicsHelpers.Blit(m_transmittanceLUT, m_precomputeMaterial, (int)PrecomputePass.ComputeTransmittance);

            // Compute Direct Irradiance into DeltaIrradianceTexture and either initialize m_irradianceLUT with 0 or leave it unchanged. (We don't want the direct irradiance in m_irradianceLUT, but only the irradiance from the sky)
            if (accumulate)
                Utilities.GraphicsHelpers.Blit3D(DeltaIrradianceTexture, m_precomputeMaterial, (int)PrecomputePass.ComputeDirectIrradiance);
            else
                Utilities.GraphicsHelpers.Blit3D(new RenderTexture[2] { DeltaIrradianceTexture, m_irradianceLUT }, m_precomputeMaterial, (int)PrecomputePass.ComputeDirectIrradiance);

            // Compute Raylie and Mie Single Scattering, and store them in DeltaRayleighScatteringTexture, DeltaMieScatteringTexture as well as scatteringLUT
            if (accumulate)
                Utilities.GraphicsHelpers.Blit3D(new RenderTexture[2] { DeltaRayleighScatteringTexture, DeltaMieScatteringTexture }, m_precomputeMaterial, (int)PrecomputePass.ComputeSingleScattering);
            else
                Utilities.GraphicsHelpers.Blit3D(new RenderTexture[3] { DeltaRayleighScatteringTexture, DeltaMieScatteringTexture, m_scatteringLUT }, m_precomputeMaterial, (int)PrecomputePass.ComputeSingleScattering);
            Utilities.GraphicsHelpers.Blit3D(m_scatteringLUT, m_precomputeMaterial, (int)PrecomputePass.AccumulateSingleScattering);

            // Compute to the nth order of scattering, in sequence
            for (int scatteringOrder = 2; scatteringOrder <= m_scatteringOrders; ++scatteringOrder)
            {
                m_precomputeMaterial.SetInt("scattering_order", scatteringOrder);

                // Compute the scattering density, and store it in DeltaScatteringDensityTexture.
                Utilities.GraphicsHelpers.Blit3D(DeltaScatteringDensityTexture, m_precomputeMaterial, (int)PrecomputePass.ComputeScatteringDensity);

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
            SunAngularRadius = Mathf.Clamp(SunAngularRadius, 0.001f, 0.1f);
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