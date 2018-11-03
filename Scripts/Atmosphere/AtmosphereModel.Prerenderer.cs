/**
* Copyright (c) 2017 Eric Bruneton
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
* 1. Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in the
*    documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holders nor the names of its
*    contributors may be used to endorse or promote products derived from
*    this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED. in NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER in
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
* ARISING in ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
* THE POSSIBILITY OF SUCH DAMAGE.
*
* Precomputed Atmospheric Scattering
* Copyright (c) 2008 INRIA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
* 1. Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in the
*    documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holders nor the names of its
*    contributors may be used to endorse or promote products derived from
*    this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED. in NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER in
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
* ARISING in ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
* THE POSSIBILITY OF SUCH DAMAGE.
*/

using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace PhysicalSky
{
    public partial class AtmosphereModel
    {
        private abstract class Prerenderer
        {
            protected static Vector4[] GetDensityProfileShaderValues(AtmosphereLib.DensityProfile dp)
            {
                return new Vector4[3] {
                        new Vector4(
                            dp.layer0.width / LENGTH_UNIT_IN_METERS,
                            dp.layer0.exp_term,
                            dp.layer0.exp_scale * LENGTH_UNIT_IN_METERS,
                            dp.layer0.linear_term * LENGTH_UNIT_IN_METERS),
                        new Vector4(
                            dp.layer0.constant_term,
                            dp.layer1.width / LENGTH_UNIT_IN_METERS,
                            dp.layer1.exp_term,
                            dp.layer1.exp_scale * LENGTH_UNIT_IN_METERS),
                        new Vector4(
                            dp.layer1.linear_term * LENGTH_UNIT_IN_METERS,
                            dp.layer1.constant_term,
                            0,
                            0)};
            }

            protected AtmosphereParameters m_parameters;

            protected float m_sunSolidAngle;
            protected Vector3 m_sky_k;
            protected Vector3 m_sun_k;
            protected AtmosphereLib.DensityProfile m_rayleighDensity;
            protected AtmosphereLib.DensityProfile m_mieDensity;
            protected AtmosphereLib.DensityProfile m_absorptionDensity;

            protected const int WAVELENGTH_STEP_SIZE = 10;
            protected const int WAVELENGTH_STEP_COUNT = ((LAMBDA_MAX - LAMBDA_MIN) / WAVELENGTH_STEP_SIZE) + 1;
            protected float[] m_wavelengths = new float[WAVELENGTH_STEP_COUNT];
            protected float[] m_solarIrradiance = new float[WAVELENGTH_STEP_COUNT];
            protected float[] m_rayleighScattering = new float[WAVELENGTH_STEP_COUNT];
            protected float[] m_mieScattering = new float[WAVELENGTH_STEP_COUNT];
            protected float[] m_mieExtinction = new float[WAVELENGTH_STEP_COUNT];
            protected float[] m_absorptionExtinction = new float[WAVELENGTH_STEP_COUNT];

            protected int m_numPrecomputedWavelengths;

            protected RenderTexture m_transmittanceLUT = null;
            protected RenderTexture m_scatteringLUT = null;
            protected RenderTexture m_irradianceLUT = null;

            protected struct TemporaryTextures
            {
                public RenderTexture DeltaRayleighScatteringTexture;
                public RenderTexture DeltaMieScatteringTexture;
                public RenderTexture DeltaMultipleScatteringTexture;
                public RenderTexture DeltaIrradianceTexture;
                public RenderTexture DeltaScatteringDensityTexture;
            };
            protected TemporaryTextures m_tempTextures;

            protected class Matrix3x3
            {
                Vector3 a, b, c;

                public Matrix3x3(float[] f)
                {
                    a.x = f[0];
                    a.y = f[1];
                    a.z = f[2];
                    b.x = f[3];
                    b.y = f[4];
                    b.z = f[5];
                    c.x = f[6];
                    c.y = f[7];
                    c.z = f[8];
                }
                public static implicit operator float[] (Matrix3x3 m)
                {
                    return new float[9] {
                        m.a.x, m.a.y, m.a.z,
                        m.b.x, m.b.y, m.b.z,
                        m.c.x, m.c.y, m.c.z};
                }
                public static implicit operator Vector4[] (Matrix3x3 m)
                {
                    return new Vector4[3] {
                        new Vector4(m.a.x, m.a.y, m.a.z, 0),
                        new Vector4(m.b.x, m.b.y, m.b.z, 0),
                        new Vector4( m.c.x, m.c.y, m.c.z, 0)};
                }

                public static Matrix3x3 identity
                {
                    get
                    {
                        return new Matrix3x3(new float[9]{
                            1.0f, 0.0f, 0.0f,
                            0.0f, 1.0f, 0.0f,
                            0.0f, 0.0f, 1.0f});
                    }
                }
            }

            // Abstract functions
            protected abstract void RenderTransmittanceLUT();
            protected abstract void SetupPrerender(Vector3 lamdas, Matrix3x3 luminanceFromRadiance);
            protected abstract void CleanupPrerenderer();
            protected abstract void RenderSingleScattering(bool accumulate);
            protected abstract void RenderMultipleScattering(int scatteringOrder);

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

            protected Vector4 ScaleToWavelengths(float[] v, Vector3 lambdas, float scale)
            {
                return new Vector4(
                    (Interpolate(m_wavelengths, v, lambdas.x) * scale),
                    (Interpolate(m_wavelengths, v, lambdas.y) * scale),
                    (Interpolate(m_wavelengths, v, lambdas.z) * scale),
                    1.0f);
            }

            public bool Compute(AtmosphereParameters parameters, RenderTexture transmittanceLUT, RenderTexture scatteringLUT, RenderTexture irradianceLUT, ref AtmosphereLib.AtmosphereRenderParams renderValues)
            {
#if PHYSICAL_SKY_DEBUG
                Debug.Log("Computing Atmospheric Lookup Textures");
                float timerStartCompute = Time.realtimeSinceStartup;
#endif

                m_parameters = parameters;

                m_sunSolidAngle = Mathf.PI * Mathf.Pow(m_parameters.sunAngularRadius, 2);


                for (int l = LAMBDA_MIN, i = 0; l <= LAMBDA_MAX; l += WAVELENGTH_STEP_SIZE, i++)
                {
                    float lambda = l * 1e-3f;  // micro-meters
                    float mie = m_parameters.mieAngstromBeta / m_parameters.mieScaleHeight * Mathf.Pow(lambda, -m_parameters.mieAngstromAlpha);
                    m_wavelengths[i] = l;

                    if (m_parameters.useConstantSolarSpectrum)
                        m_solarIrradiance[i] = m_parameters.constantSolarIrradiance;
                    else
                        m_solarIrradiance[i] = SOLAR_IRRADIANCE[(l - LAMBDA_MIN) / 10];

                    m_rayleighScattering[i] = m_parameters.rayleigh * Mathf.Pow(lambda, -4);
                    m_mieScattering[i] = mie * m_parameters.mieSingleScatteringAlbedo;
                    m_mieExtinction[i] = mie;
                    m_absorptionExtinction[i] = m_parameters.useOzone ? kMaxOzoneNumberDensity * OZONE_CROSS_SECTION[(l - LAMBDA_MIN) / 10] : 0.0f;
                }

                m_rayleighDensity = new AtmosphereLib.DensityProfile(new AtmosphereLib.DensityProfileLayer(0.0f, 1.0f, -1.0f / m_parameters.rayleighScaleHeight, 0.0f, 0.0f));
                m_mieDensity = new AtmosphereLib.DensityProfile(new AtmosphereLib.DensityProfileLayer(0.0f, 1.0f, -1.0f / m_parameters.mieScaleHeight, 0.0f, 0.0f));

                // Density profile increasing linearly from 0 to 1 between 10 and 25km, and
                // decreasing linearly from 1 to 0 between 25 and 40km. This is an approximate
                // profile from http://www.kln.ac.lk/science/Chemistry/Teaching_Resources/
                // Documents/Introduction%20to%20atmospheric%20chemistry.pdf (page 10).
                // (ozone density)
                m_absorptionDensity = new AtmosphereLib.DensityProfile(
                    new AtmosphereLib.DensityProfileLayer(25000.0f, 0.0f, 0.0f, 1.0f / 15000.0f, -2.0f / 3.0f),
                    new AtmosphereLib.DensityProfileLayer(0.0f, 0.0f, 0.0f, -1.0f / 15000.0f, 8.0f / 3.0f));

                m_numPrecomputedWavelengths = m_parameters.luminance == AtmosphereParameters.LuminanceType.precomputed ? 15 : 3;

                // Compute the values for sky_spectral_radiance_to_luminance. In theory
                // this should be 1 in precomputed illuminance mode (because the precomputed
                // textures already contain illuminance values). In practice, however, storing
                // true illuminance values in half precision textures yields artefacts
                // (because the values are too large), so we store illuminance values divided
                // by MAX_LUMINOUS_EFFICACY instead. This is why, in precomputed illuminance
                // mode, we set sky_spectral_radiance_to_luminance to MAX_LUMINOUS_EFFICACY.
                if (m_numPrecomputedWavelengths > 3)
                    m_sky_k = new Vector3(MAX_LUMINOUS_EFFICACY, MAX_LUMINOUS_EFFICACY, MAX_LUMINOUS_EFFICACY);
                else
                    m_sky_k = ComputeSpectralRadianceToLuminanceFactors(m_wavelengths, m_solarIrradiance, -3 /* lambda_power */);

                m_sun_k = ComputeSpectralRadianceToLuminanceFactors(m_wavelengths, m_solarIrradiance, 0 /* lambda_power */);

                m_transmittanceLUT = transmittanceLUT;
                m_scatteringLUT = scatteringLUT;
                m_irradianceLUT = irradianceLUT;

                m_tempTextures.DeltaIrradianceTexture = TextureFactory.CreateTempRenderTexture(TextureFactory.Preset.Irradiance);
                m_tempTextures.DeltaRayleighScatteringTexture = TextureFactory.CreateTempRenderTexture(TextureFactory.Preset.Scattering);
                m_tempTextures.DeltaMieScatteringTexture = TextureFactory.CreateTempRenderTexture(TextureFactory.Preset.Scattering);
                m_tempTextures.DeltaScatteringDensityTexture = TextureFactory.CreateTempRenderTexture(TextureFactory.Preset.Scattering);

                // DeltaMultipleScatteringTexture is only needed to compute scattering
                // order 3 or more, while DeltaRayleighScatteringTexture and
                // DeltaMieScatteringTexture are only needed to compute double scattering.
                // Therefore, to save memory, we can store DeltaRayleighScatteringTexture
                // and DeltaMultipleScatteringTexture in the same GPU texture.
                m_tempTextures.DeltaMultipleScatteringTexture = m_tempTextures.DeltaRayleighScatteringTexture;

                if (m_numPrecomputedWavelengths <= 3)
                {
                    SetupPrerender(LAMBDA_V, Matrix3x3.identity);
                    RenderTransmittanceLUT();
                    RenderSingleScattering(false);

                    for (int scatteringOrder = 2; scatteringOrder <= m_parameters.scatteringOrders; ++scatteringOrder)
                    {
                        RenderMultipleScattering(scatteringOrder);
                    }
                }
                else
                {
                    int num_iterations = (m_numPrecomputedWavelengths + 2) / 3;
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

                        Matrix3x3 luminanceFromRadiance = new Matrix3x3(new float[9]{
                            coeff(lambdas[0], 0), coeff(lambdas[1], 0), coeff(lambdas[2], 0),
                            coeff(lambdas[0], 1), coeff(lambdas[1], 1), coeff(lambdas[2], 1),
                            coeff(lambdas[0], 2), coeff(lambdas[1], 2), coeff(lambdas[2], 2)});

                        SetupPrerender(lambdas, luminanceFromRadiance);
                        RenderTransmittanceLUT();
                        RenderSingleScattering(i != 0);

                        for (int scatteringOrder = 2; scatteringOrder <= m_parameters.scatteringOrders; ++scatteringOrder)
                        {
                            RenderMultipleScattering(scatteringOrder);
                        }
                    }

                    // After the above iterations, the transmittance texture contains the
                    // transmittance for the 3 wavelengths used at the last iteration. But we
                    // want the transmittance at kLambdaR, kLambdaG, kLambdaB instead, so we
                    // must recompute it here for these 3 wavelengths:
                    SetupPrerender(LAMBDA_V, Matrix3x3.identity);
                    RenderTransmittanceLUT();
                }

                CleanupPrerenderer();
                // Release temporary textures.
                TextureFactory.ReleaseTempRenderTexture(ref m_tempTextures.DeltaIrradianceTexture);
                TextureFactory.ReleaseTempRenderTexture(ref m_tempTextures.DeltaRayleighScatteringTexture);
                TextureFactory.ReleaseTempRenderTexture(ref m_tempTextures.DeltaMieScatteringTexture);
                TextureFactory.ReleaseTempRenderTexture(ref m_tempTextures.DeltaScatteringDensityTexture);

                m_transmittanceLUT = null;
                m_scatteringLUT = null;
                m_irradianceLUT = null;

                // Assign to computed values
                renderValues.sky_spectral_radiance_to_luminance = parameters.luminance != AtmosphereParameters.LuminanceType.none ? m_sky_k : Vector3.one;
                renderValues.sun_spectral_radiance_to_luminance = parameters.luminance != AtmosphereParameters.LuminanceType.none ? m_sun_k : Vector3.one;
                renderValues.solar_irradiance = ScaleToWavelengths(m_solarIrradiance, LAMBDA_V, 1.0f);
                renderValues.sun_angular_radius = parameters.sunAngularRadius;
                renderValues.bottom_radius = m_parameters.planetaryRadius / LENGTH_UNIT_IN_METERS;
                renderValues.top_radius = (m_parameters.planetaryRadius + m_parameters.atmosphereThickness) / LENGTH_UNIT_IN_METERS;
                renderValues.rayleigh_density = m_rayleighDensity;
                renderValues.rayleigh_scattering = ScaleToWavelengths(m_rayleighScattering, LAMBDA_V, LENGTH_UNIT_IN_METERS);
                renderValues.mie_density = m_mieDensity;
                renderValues.mie_scattering = ScaleToWavelengths(m_mieScattering, LAMBDA_V, LENGTH_UNIT_IN_METERS);
                renderValues.mie_extinction = ScaleToWavelengths(m_mieExtinction, LAMBDA_V, LENGTH_UNIT_IN_METERS);
                renderValues.mie_phase_function_g = parameters.miePhaseFunctionG;
                renderValues.absorption_density = m_absorptionDensity;
                renderValues.absorption_extinction = ScaleToWavelengths(m_absorptionExtinction, LAMBDA_V, LENGTH_UNIT_IN_METERS);
                renderValues.ground_albedo = new Vector3(parameters.groundAlbedo, parameters.groundAlbedo, parameters.groundAlbedo);
                renderValues.mu_s_min = parameters.maxSunZenithAngle;

#if PHYSICAL_SKY_DEBUG
                float timerEndCompute = Time.realtimeSinceStartup;
                Debug.Log("Computed atmospheric lookup textures in " + (timerEndCompute - timerStartCompute) * 1000.0f + "ms");
#endif

                return true;
            }
        }
    }
}