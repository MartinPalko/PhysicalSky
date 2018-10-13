using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace PhysicalSky
{
    /// <summary>
    /// Handles the majority of the internal functionality of the AtmosphereModel
    /// </summary>
    public partial class AtmosphereModel
    {
        private static class Internal
        {
            public class TextureFactory
            {
                // Same as defined in PhysicalSkyCommon
                const int TRANSMITTANCE_TEXTURE_WIDTH = 256;
                const int TRANSMITTANCE_TEXTURE_HEIGHT = 64;
                const int SCATTERING_TEXTURE_R_SIZE = 32;
                const int SCATTERING_TEXTURE_MU_SIZE = 128;
                const int SCATTERING_TEXTURE_MU_S_SIZE = 32;
                const int SCATTERING_TEXTURE_NU_SIZE = 8;
                const int SCATTERING_TEXTURE_WIDTH = SCATTERING_TEXTURE_NU_SIZE * SCATTERING_TEXTURE_MU_S_SIZE;
                static readonly int SCATTERING_TEXTURE_HEIGHT = Supports3DTextures() ? SCATTERING_TEXTURE_MU_SIZE : SCATTERING_TEXTURE_R_SIZE * SCATTERING_TEXTURE_MU_SIZE;
                static readonly int SCATTERING_TEXTURE_DEPTH = Supports3DTextures() ? SCATTERING_TEXTURE_R_SIZE : 1;
                const int IRRADIANCE_TEXTURE_WIDTH = 64;
                const int IRRADIANCE_TEXTURE_HEIGHT = 16;

                public enum Preset
                {
                    Transmittance,
                    Irradiance,
                    Scattering
                }

                public static bool Supports3DTextures()
                {
                    return true;
                }

                private static FilterMode GetFilterMode(Preset preset)
                {
                    if (preset == Preset.Scattering)
                        return FilterMode.Trilinear;
                    else
                        return FilterMode.Bilinear;
                }

                private static RenderTextureDescriptor GetTextureDescriptor(Preset preset)
                {
                    RenderTextureDescriptor desc;
                    switch (preset)
                    {
                        case Preset.Transmittance:
                            desc = new RenderTextureDescriptor(TRANSMITTANCE_TEXTURE_WIDTH, TRANSMITTANCE_TEXTURE_HEIGHT);
                            break;
                        case Preset.Irradiance:
                            desc = new RenderTextureDescriptor(IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT);
                            break;
                        case Preset.Scattering:
                            desc = new RenderTextureDescriptor(SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT);
                            if (Supports3DTextures())
                            {
                                desc.dimension = UnityEngine.Rendering.TextureDimension.Tex3D;
                                desc.volumeDepth = SCATTERING_TEXTURE_DEPTH;
                            }
                            break;
                        default:
                            throw new System.NotImplementedException();
                    }

                    desc.useMipMap = false;
                    desc.sRGB = false;
                    desc.colorFormat = RenderTextureFormat.ARGBFloat; // TODO: Try changing to half or float
                    desc.depthBufferBits = 0;
                    return desc;
                }

                public static void CreateRenderTexture(ref RenderTexture tex, Preset preset)
                {
                    if (tex == null)
                        tex = new RenderTexture(GetTextureDescriptor(preset));

                    if (!tex.IsCreated())
                    {
                        tex.wrapMode = TextureWrapMode.Clamp;
                        tex.filterMode = GetFilterMode(preset);
                        tex.Create();
                    }
                }

                public static void ReleaseRenderTexture(ref RenderTexture tex)
                {
                    if (tex != null)
                    {
                        tex.Release();
                    }
                    tex = null;
                }

                public static RenderTexture CreateTempRenderTexture(Preset preset)
                {
                    RenderTexture tex = RenderTexture.GetTemporary(GetTextureDescriptor(preset));
                    tex.wrapMode = TextureWrapMode.Clamp;
                    tex.filterMode = GetFilterMode(preset);
                    return tex;
                }

                public static void ReleaseTempRenderTexture(ref RenderTexture tex)
                {
                    RenderTexture.ReleaseTemporary(tex);
                    tex = null;
                }
            }

            [Serializable]
            public struct ComputedValues
            {
                public float m_sunSolidAngle;
                public Vector3 m_sky_k;
                public Vector3 m_sun_k;
                public DensityProfile m_rayleighDensity;
                public DensityProfile m_mieDensity;
                public DensityProfile m_absorptionDensity;
                public Vector3 m_solarIrradiance;
                public Vector3 m_rayleighScattering;
                public Vector3 m_mieScattering;
                public Vector3 m_mieExtinction;
                public Vector3 m_absorptionExtinction;
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

            // Constants
            public const float LENGTH_UNIT_IN_METERS = 1000.0f;

            // The conversion factor between watts and lumens.
            const float MAX_LUMINOUS_EFFICACY = 683.0f;

            const float LAMBDA_R = 680.0f;
            const float LAMBDA_G = 550.0f;
            const float LAMBDA_B = 440.0f;
            readonly static Vector3 LAMBDA_V = new Vector3(LAMBDA_R, LAMBDA_G, LAMBDA_B);

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

            private static float Interpolate(float[] wavelengths, float[] wavelength_function, float wavelength)
            {
                Debug.Assert(wavelengths.Length > 0);
                Debug.Assert(wavelength_function.Length == wavelengths.Length);

                if (wavelength < wavelengths[0])
                {
                    return wavelength_function[0];
                }
                for (int i = 0; i < wavelengths.Length - 1; ++i)
                {
                    if (wavelength < wavelengths[i + 1])
                    {
                        float u = (wavelength - wavelengths[i]) / (wavelengths[i + 1] - wavelengths[i]);
                        return wavelength_function[i] * (1.0f - u) + wavelength_function[i + 1] * u;
                    }
                }
                return wavelength_function[wavelength_function.Length - 1];
            }

            private static float CieColorMatchingFunctionTableValue(float wavelength, int column)
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
            private static Vector3 ComputeSpectralRadianceToLuminanceFactors(
                float[] wavelengths,
                float[] solar_irradiance,
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

            private static Vector4 ScaleToWavelengths(float[] wavelengths, float[] v, Vector3 lambdas, float scale)
            {
                return new Vector4(
                    (Interpolate(wavelengths, v, lambdas.x) * scale),
                    (Interpolate(wavelengths, v, lambdas.y) * scale),
                    (Interpolate(wavelengths, v, lambdas.z) * scale),
                    1.0f);
            }

            private static void SetDensityProfileLayerShaderUniforms(DensityProfileLayer l, string name, Material m)
            {
                m.SetFloatArray(name, new float[5] {
                    l.width / LENGTH_UNIT_IN_METERS,
                    l.exp_term,
                    l.exp_scale * LENGTH_UNIT_IN_METERS,
                    l.linear_term * LENGTH_UNIT_IN_METERS,
                    l.constant_term });
            }

            private static void SetDensityProfileShaderUniforms(DensityProfile p, string name, Material m)
            {
                SetDensityProfileLayerShaderUniforms(p.layer0, name + "0", m);
                SetDensityProfileLayerShaderUniforms(p.layer1, name + "1", m);
            }

            public static void SetShaderUniforms(AtmosphereModel model, Material material)
            {
                material.SetTexture("_transmittance_texture", model.m_transmittanceLUT);
                material.SetTexture("_scattering_texture", model.m_scatteringLUT);
                material.SetTexture("_irradiance_texture", model.m_irradianceLUT);
                SetShaderUniforms(model.m_computedValues, model.m_computedParameters, material);
            }

            private static void SetShaderUniforms(ComputedValues values, AtmosphereParameters parameters, Material material)
            {
                material.SetInt("_use_luminance", parameters.luminance == AtmosphereParameters.LuminanceType.none ? 0 : 1);
                material.SetVector("_sky_spectral_radiance_to_luminance", parameters.luminance != AtmosphereParameters.LuminanceType.none ? values.m_sky_k : Vector3.one);
                material.SetVector("_sun_spectral_radiance_to_luminance", parameters.luminance != AtmosphereParameters.LuminanceType.none ? values.m_sun_k : Vector3.one);

                material.SetVector("_solar_irradiance", values.m_solarIrradiance);
                material.SetFloat("_sun_angular_radius", parameters.sunAngularRadius);
                material.SetFloat("_bottom_radius", parameters.planetaryRadius / LENGTH_UNIT_IN_METERS);
                material.SetFloat("_top_radius", (parameters.planetaryRadius + parameters.atmosphereThickness) / LENGTH_UNIT_IN_METERS);
                SetDensityProfileShaderUniforms(values.m_rayleighDensity, "_rayleigh_density", material);
                material.SetVector("_rayleigh_scattering", values.m_rayleighScattering);
                SetDensityProfileShaderUniforms(values.m_mieDensity, "_mie_density", material);
                material.SetVector("_mie_scattering", values.m_mieScattering);
                material.SetVector("_mie_extinction", values.m_mieExtinction);
                material.SetFloat("_mie_phase_function_g", parameters.miePhaseFunctionG);
                SetDensityProfileShaderUniforms(values.m_absorptionDensity, "_absorption_density", material);
                material.SetVector("_absorption_extinction", values.m_absorptionExtinction);
                material.SetVector("_ground_albedo", Vector3.one * parameters.groundAlbedo);
                material.SetFloat("_mu_s_min", Mathf.Cos(parameters.maxSunZenithAngle));

                material.SetVector("sun_size", new Vector3(Mathf.Tan(parameters.sunAngularRadius), Mathf.Cos(parameters.sunAngularRadius), parameters.sunAngularRadius));
            }

            public static void Compute(AtmosphereModel model, AtmosphereParameters parameters)
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

                ComputedValues values;

                values.m_sunSolidAngle = Mathf.PI * Mathf.Pow(parameters.sunAngularRadius, 2);
                const int wavelengthStepSize = 10;
                const int numWavelengthSteps = ((LAMBDA_MAX - LAMBDA_MIN) / wavelengthStepSize) + 1;

                float[] wavelengths = new float[numWavelengthSteps];
                float[] solarIrradiance = new float[numWavelengthSteps];
                float[] rayleighScattering = new float[numWavelengthSteps];
                float[] mieScattering = new float[numWavelengthSteps];
                float[] mieExtinction = new float[numWavelengthSteps];
                float[] absorptionExtinction = new float[numWavelengthSteps];

                {
                    int i = 0;
                    for (int l = LAMBDA_MIN; l <= LAMBDA_MAX; l += wavelengthStepSize)
                    {
                        float lambda = l * 1e-3f;  // micro-meters
                        float mie = parameters.mieAngstromBeta / parameters.mieScaleHeight * Mathf.Pow(lambda, -parameters.mieAngstromAlpha);
                        wavelengths[i] = l;

                        if (parameters.useConstantSolarSpectrum)
                            solarIrradiance[i] = parameters.constantSolarIrradiance;
                        else
                            solarIrradiance[i] = SOLAR_IRRADIANCE[(l - LAMBDA_MIN) / 10];

                        rayleighScattering[i] = parameters.rayleigh * Mathf.Pow(lambda, -4);
                        mieScattering[i] = mie * parameters.mieSingleScatteringAlbedo;
                        mieExtinction[i] = mie;
                        absorptionExtinction[i] = parameters.useOzone ? kMaxOzoneNumberDensity * OZONE_CROSS_SECTION[(l - LAMBDA_MIN) / 10] : 0.0f;

                        i++;
                    }
                }

                values.m_rayleighDensity = new DensityProfile(new DensityProfileLayer(0.0f, 1.0f, -1.0f / parameters.rayleighScaleHeight, 0.0f, 0.0f));
                values.m_mieDensity = new DensityProfile(new DensityProfileLayer(0.0f, 1.0f, -1.0f / parameters.mieScaleHeight, 0.0f, 0.0f));

                // Density profile increasing linearly from 0 to 1 between 10 and 25km, and
                // decreasing linearly from 1 to 0 between 25 and 40km. This is an approximate
                // profile from http://www.kln.ac.lk/science/Chemistry/Teaching_Resources/
                // Documents/Introduction%20to%20atmospheric%20chemistry.pdf (page 10).
                // (ozone density)
                values.m_absorptionDensity = new DensityProfile(
                    new DensityProfileLayer(25000.0f, 0.0f, 0.0f, 1.0f / 15000.0f, -2.0f / 3.0f),
                    new DensityProfileLayer(0.0f, 0.0f, 0.0f, -1.0f / 15000.0f, 8.0f / 3.0f));

                int numPrecomputedWavelengths = parameters.luminance == AtmosphereParameters.LuminanceType.precomputed ? 15 : 3;
                bool precomputeIlluminance = numPrecomputedWavelengths > 3;

                // Compute the values for sky_spectral_radiance_to_luminance. In theory
                // this should be 1 in precomputed illuminance mode (because the precomputed
                // textures already contain illuminance values). In practice, however, storing
                // true illuminance values in half precision textures yields artefacts
                // (because the values are too large), so we store illuminance values divided
                // by MAX_LUMINOUS_EFFICACY instead. This is why, in precomputed illuminance
                // mode, we set sky_spectral_radiance_to_luminance to MAX_LUMINOUS_EFFICACY.
                if (precomputeIlluminance)
                    values.m_sky_k = new Vector3(MAX_LUMINOUS_EFFICACY, MAX_LUMINOUS_EFFICACY, MAX_LUMINOUS_EFFICACY);
                else
                    values.m_sky_k = ComputeSpectralRadianceToLuminanceFactors(wavelengths, solarIrradiance, -3 /* lambda_power */);

                values.m_sun_k = ComputeSpectralRadianceToLuminanceFactors(wavelengths, solarIrradiance, 0 /* lambda_power */);

                Material precomputeMaterial = new Material(model.m_PrecomputeShader);

                if (numPrecomputedWavelengths <= 3)
                {
                    Vector3 lambdas = LAMBDA_V;
                    values.m_solarIrradiance = ScaleToWavelengths(wavelengths, solarIrradiance, lambdas, 1.0f);
                    values.m_rayleighScattering = ScaleToWavelengths(wavelengths, rayleighScattering, lambdas, LENGTH_UNIT_IN_METERS);
                    values.m_mieScattering = ScaleToWavelengths(wavelengths, mieScattering, lambdas, LENGTH_UNIT_IN_METERS);
                    values.m_mieExtinction = ScaleToWavelengths(wavelengths, mieExtinction, lambdas, LENGTH_UNIT_IN_METERS);
                    values.m_absorptionExtinction = ScaleToWavelengths(wavelengths, absorptionExtinction, lambdas, LENGTH_UNIT_IN_METERS);

                    float[] luminanceFromRadiance = new float[9]{
                    1.0f, 0.0f, 0.0f,
                    0.0f, 1.0f, 0.0f,
                    0.0f, 0.0f, 1.0f};

                    ComputeTextures(model, values, parameters, precomputeMaterial, luminanceFromRadiance, false);
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
                        values.m_solarIrradiance = ScaleToWavelengths(wavelengths, solarIrradiance, lambdas, 1.0f);
                        values.m_rayleighScattering = ScaleToWavelengths(wavelengths, rayleighScattering, lambdas, LENGTH_UNIT_IN_METERS);
                        values.m_mieScattering = ScaleToWavelengths(wavelengths, mieScattering, lambdas, LENGTH_UNIT_IN_METERS);
                        values.m_mieExtinction = ScaleToWavelengths(wavelengths, mieExtinction, lambdas, LENGTH_UNIT_IN_METERS);
                        values.m_absorptionExtinction = ScaleToWavelengths(wavelengths, absorptionExtinction, lambdas, LENGTH_UNIT_IN_METERS);

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

                        ComputeTextures(model, values, parameters, precomputeMaterial, luminanceFromRadiance, i == 0 ? false : true);
                    }

                    // After the above iterations, the transmittance texture contains the
                    // transmittance for the 3 wavelengths used at the last iteration. But we
                    // want the transmittance at kLambdaR, kLambdaG, kLambdaB instead, so we
                    // must recompute it here for these 3 wavelengths:

                    values.m_solarIrradiance = ScaleToWavelengths(wavelengths, solarIrradiance, LAMBDA_V, 1.0f);
                    values.m_rayleighScattering = ScaleToWavelengths(wavelengths, rayleighScattering, LAMBDA_V, LENGTH_UNIT_IN_METERS);
                    values.m_mieScattering = ScaleToWavelengths(wavelengths, mieScattering, LAMBDA_V, LENGTH_UNIT_IN_METERS);
                    values.m_mieExtinction = ScaleToWavelengths(wavelengths, mieExtinction, LAMBDA_V, LENGTH_UNIT_IN_METERS);
                    values.m_absorptionExtinction = ScaleToWavelengths(wavelengths, absorptionExtinction, LAMBDA_V, LENGTH_UNIT_IN_METERS);

                    SetShaderUniforms(values, parameters, precomputeMaterial);

                    Utilities.GraphicsHelpers.Blit(model.m_transmittanceLUT, precomputeMaterial, (int)PrecomputePass.ComputeTransmittance);

                    model.m_computedParameters = parameters;
                    model.m_computedValues = values;
                }

                float timerEndCompute = Time.realtimeSinceStartup;
#if PHYSICAL_SKY_DEBUG
                Debug.Log("Computed atmospheric lookup textures in " + (timerEndCompute - timerStartCompute) * 1000.0f + "ms");
#endif
            }

            private static void ComputeTextures(AtmosphereModel model, ComputedValues values, AtmosphereParameters parameters, Material precomputeMaterial, float[] luminanceFromRadiance, bool accumulate)
            {
                Debug.Assert(luminanceFromRadiance.Length == 9);

                RenderTexture DeltaIrradianceTexture = TextureFactory.CreateTempRenderTexture(TextureFactory.Preset.Irradiance);
                RenderTexture DeltaRayleighScatteringTexture = TextureFactory.CreateTempRenderTexture(TextureFactory.Preset.Scattering);
                RenderTexture DeltaMieScatteringTexture = TextureFactory.CreateTempRenderTexture(TextureFactory.Preset.Scattering);
                RenderTexture DeltaScatteringDensityTexture = TextureFactory.CreateTempRenderTexture(TextureFactory.Preset.Scattering);

                // DeltaMultipleScatteringTexture is only needed to compute scattering
                // order 3 or more, while DeltaRayleighScatteringTexture and
                // DeltaMieScatteringTexture are only needed to compute double scattering.
                // Therefore, to save memory, we can store DeltaRayleighScatteringTexture
                // and DeltaMultipleScatteringTexture in the same GPU texture.
                RenderTexture DeltaMultipleScatteringTexture = DeltaRayleighScatteringTexture;

                SetShaderUniforms(values, parameters, precomputeMaterial);
                precomputeMaterial.SetFloatArray("_luminance_from_radiance", luminanceFromRadiance);
                precomputeMaterial.SetTexture("transmittance_texture", model.m_transmittanceLUT);
                precomputeMaterial.SetTexture("single_rayleigh_scattering_texture", DeltaRayleighScatteringTexture);
                precomputeMaterial.SetTexture("single_mie_scattering_texture", DeltaMieScatteringTexture);
                precomputeMaterial.SetTexture("multiple_scattering_texture", DeltaMultipleScatteringTexture);
                precomputeMaterial.SetTexture("irradiance_texture", DeltaIrradianceTexture);
                precomputeMaterial.SetTexture("scattering_density_texture", DeltaScatteringDensityTexture);

                // Compute Transmittance LUT.
                Utilities.GraphicsHelpers.Blit(model.m_transmittanceLUT, precomputeMaterial, (int)PrecomputePass.ComputeTransmittance);

                // Compute Direct Irradiance into DeltaIrradianceTexture and either initialize m_irradianceLUT with 0 or leave it unchanged. (We don't want the direct irradiance in m_irradianceLUT, but only the irradiance from the sky)
                if (accumulate)
                    Utilities.GraphicsHelpers.Blit3D(DeltaIrradianceTexture, precomputeMaterial, (int)PrecomputePass.ComputeDirectIrradiance);
                else
                    Utilities.GraphicsHelpers.Blit3D(new RenderTexture[2] { DeltaIrradianceTexture, model.m_irradianceLUT }, precomputeMaterial, (int)PrecomputePass.ComputeDirectIrradiance);

                // Compute Raylie and Mie Single Scattering, and store them in DeltaRayleighScatteringTexture, DeltaMieScatteringTexture as well as scatteringLUT
                if (accumulate)
                    Utilities.GraphicsHelpers.Blit3D(new RenderTexture[2] { DeltaRayleighScatteringTexture, DeltaMieScatteringTexture }, precomputeMaterial, (int)PrecomputePass.ComputeSingleScattering);
                else
                    Utilities.GraphicsHelpers.Blit3D(new RenderTexture[3] { DeltaRayleighScatteringTexture, DeltaMieScatteringTexture, model.m_scatteringLUT }, precomputeMaterial, (int)PrecomputePass.ComputeSingleScattering);
                Utilities.GraphicsHelpers.Blit3D(model.m_scatteringLUT, precomputeMaterial, (int)PrecomputePass.AccumulateSingleScattering);
                // Compute to the nth order of scattering, in sequence
                for (int scatteringOrder = 2; scatteringOrder <= parameters.scatteringOrders; ++scatteringOrder)
                {
                    precomputeMaterial.SetInt("scattering_order", scatteringOrder);

                    // Compute the scattering density, and store it in DeltaScatteringDensityTexture.
                    Utilities.GraphicsHelpers.Blit3D(DeltaScatteringDensityTexture, precomputeMaterial, (int)PrecomputePass.ComputeScatteringDensity);

                    // Compute the indirect irradiance, store it in DeltaIrradianceTexture and accumulate it in irradianceLUT.
                    precomputeMaterial.SetInt("scattering_order", scatteringOrder - 1);
                    Utilities.GraphicsHelpers.Blit3D(DeltaIrradianceTexture, precomputeMaterial, (int)PrecomputePass.ComputeIndirectIrradiance);
                    Utilities.GraphicsHelpers.Blit3D(model.m_irradianceLUT, precomputeMaterial, (int)PrecomputePass.AccumulateIndirectIrradiance);

                    // Compute the multiple scattering, store it in DeltaMultipleScatteringTexture, and accumulate it in scatteringLUT.
                    Utilities.GraphicsHelpers.Blit3D(DeltaMultipleScatteringTexture, precomputeMaterial, (int)PrecomputePass.ComputeMultipleScattering);
                    Utilities.GraphicsHelpers.Blit3D(model.m_scatteringLUT, precomputeMaterial, (int)PrecomputePass.AccumulateMultipleScattering);
                }

                Graphics.ClearRandomWriteTargets();
                RenderTexture.active = null;

                // Release temporary textures.
                TextureFactory.ReleaseTempRenderTexture(ref DeltaIrradianceTexture);
                TextureFactory.ReleaseTempRenderTexture(ref DeltaRayleighScatteringTexture);
                TextureFactory.ReleaseTempRenderTexture(ref DeltaMieScatteringTexture);
                TextureFactory.ReleaseTempRenderTexture(ref DeltaScatteringDensityTexture);
            }
        }
    }
}