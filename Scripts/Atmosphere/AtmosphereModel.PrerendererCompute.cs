using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Rendering;

namespace PhysicalSky
{
    public partial class AtmosphereModel
    {
        private sealed class PrerendererCompute : Prerenderer
        {
            private enum KernelID : int
            {
                ComputeTransmittance = 0,
                ComputeDirectIrradiance,
                ComputeDirectIrradiance_Accumulate,
                ComputeSingleScattering,
                ComputeSingleScattering_Accumulate,
                ComputeScatteringDensity,
                ComputeIndirectIrradiance,
                ComputeMultipleScattering,
                count
            }

            ComputeShader m_computeShader = null;

            public PrerendererCompute() : base()
            {
                m_computeShader = PhysicalSkyProjectSettings.PrerenderComputeShader;
            }

            public static bool Supported()
            {
                return SystemInfo.supportsComputeShaders;
            }

            protected override void SetupPrerender(Vector3 lambdas, Matrix3x3 luminanceFromRadiance)
            {
                for (int id = 0; id < (int)KernelID.count; id++)
                {
                    m_computeShader.SetTexture(id, "transmittance_texture", m_transmittanceLUT);
                    m_computeShader.SetTexture(id, "single_rayleigh_scattering_texture", m_tempTextures.DeltaRayleighScatteringTexture);
                    m_computeShader.SetTexture(id, "single_mie_scattering_texture", m_tempTextures.DeltaMieScatteringTexture);
                    m_computeShader.SetTexture(id, "multiple_scattering_texture", m_tempTextures.DeltaMultipleScatteringTexture);
                    m_computeShader.SetTexture(id, "irradiance_texture", m_tempTextures.DeltaIrradianceTexture);
                    m_computeShader.SetTexture(id, "scattering_density_texture", m_tempTextures.DeltaScatteringDensityTexture);
                }

                m_computeShader.SetVectorArray("_luminance_from_radiance", luminanceFromRadiance);

                m_computeShader.SetVector("_sky_spectral_radiance_to_luminance", m_parameters.luminance != AtmosphereParameters.LuminanceType.none ? m_sky_k : Vector3.one);
                m_computeShader.SetVector("_sun_spectral_radiance_to_luminance", m_parameters.luminance != AtmosphereParameters.LuminanceType.none ? m_sun_k : Vector3.one);

                m_computeShader.SetVector("_solar_irradiance", ScaleToWavelengths(m_solarIrradiance, lambdas, 1.0f));
                m_computeShader.SetFloat("_sun_angular_radius", m_parameters.sunAngularRadius);
                m_computeShader.SetFloat("_bottom_radius", m_parameters.planetaryRadius / LENGTH_UNIT_IN_METERS);
                m_computeShader.SetFloat("_top_radius", (m_parameters.planetaryRadius + m_parameters.atmosphereThickness) / LENGTH_UNIT_IN_METERS);
                m_computeShader.SetVectorArray("_rayleigh_density", m_rayleighDensity.GetShaderValues());
                m_computeShader.SetVector("_rayleigh_scattering", ScaleToWavelengths(m_rayleighScattering, lambdas, LENGTH_UNIT_IN_METERS));
                m_computeShader.SetVectorArray("_mie_density", m_mieDensity.GetShaderValues());
                m_computeShader.SetVector("_mie_scattering", ScaleToWavelengths(m_mieScattering, lambdas, LENGTH_UNIT_IN_METERS));
                m_computeShader.SetVector("_mie_extinction", ScaleToWavelengths(m_mieExtinction, lambdas, LENGTH_UNIT_IN_METERS));
                m_computeShader.SetFloat("_mie_phase_function_g", m_parameters.miePhaseFunctionG);
                m_computeShader.SetVectorArray("_absorption_density", m_absorptionDensity.GetShaderValues());
                m_computeShader.SetVector("_absorption_extinction", ScaleToWavelengths(m_absorptionExtinction, lambdas, LENGTH_UNIT_IN_METERS));
                m_computeShader.SetVector("_ground_albedo", Vector3.one * m_parameters.groundAlbedo);
                m_computeShader.SetFloat("_mu_s_min", Mathf.Cos(m_parameters.maxSunZenithAngle));
            }

            protected override void CleanupPrerenderer()
            {
            }

            protected override void RenderTransmittanceLUT()
            {
                // Compute Transmittance LUT.
                DispatchKernel(KernelID.ComputeTransmittance, m_transmittanceLUT);
            }

            protected override void RenderSingleScattering(bool accumulate)
            {
                // Compute Direct Irradiance into DeltaIrradianceTexture and either initialize m_irradianceLUT with 0 or leave it unchanged. (We don't want the direct irradiance in m_irradianceLUT, but only the irradiance from the sky)
                DispatchKernel(accumulate ? KernelID.ComputeDirectIrradiance_Accumulate : KernelID.ComputeDirectIrradiance, new RenderTexture[2] { m_tempTextures.DeltaIrradianceTexture, m_irradianceLUT });

                // Compute Raylie and Mie Single Scattering, and store them in DeltaRayleighScatteringTexture, DeltaMieScatteringTexture as well as scatteringLUT
                DispatchKernel(accumulate ? KernelID.ComputeSingleScattering_Accumulate : KernelID.ComputeSingleScattering, new RenderTexture[3] { m_tempTextures.DeltaRayleighScatteringTexture, m_tempTextures.DeltaMieScatteringTexture, m_scatteringLUT });
            }

            protected override void RenderMultipleScattering(int scatteringOrder)
            {
                // Compute the scattering density, and store it in DeltaScatteringDensityTexture.
                m_computeShader.SetInt("scattering_order", scatteringOrder);
                DispatchKernel(KernelID.ComputeScatteringDensity, m_tempTextures.DeltaScatteringDensityTexture);

                // Compute the indirect irradiance, store it in DeltaIrradianceTexture and accumulate it in irradianceLUT.
                m_computeShader.SetInt("scattering_order", scatteringOrder - 1);
                DispatchKernel(KernelID.ComputeIndirectIrradiance, new RenderTexture[2] { m_tempTextures.DeltaIrradianceTexture, m_irradianceLUT });

                // Compute the multiple scattering, store it in DeltaMultipleScatteringTexture, and accumulate it in scatteringLUT.
                DispatchKernel(KernelID.ComputeMultipleScattering, new RenderTexture[2] { m_tempTextures.DeltaMultipleScatteringTexture, m_scatteringLUT });
            }

            private Vector3Int GetThreadGroupsForRenderTexture(KernelID id, RenderTexture rt)
            {
                uint groupSizeX, groupSizeY, groupSizeZ;
                m_computeShader.GetKernelThreadGroupSizes((int)id, out groupSizeX, out groupSizeY, out groupSizeZ);
                return new Vector3Int(
                    rt.width / (int)groupSizeX + 1,
                    rt.height / (int)groupSizeY + 1,
                    rt.volumeDepth / (int)groupSizeZ + 1);
            }

            private void DispatchKernel(KernelID id, RenderTexture[] writeTargets)
            {
                for (int i = 0; i < writeTargets.Length; i++)
                {
                    m_computeShader.SetTexture((int)id, "Result" + i, writeTargets[i]);
                }

                Vector3Int threadGroups = GetThreadGroupsForRenderTexture(id, writeTargets[0]);
                m_computeShader.Dispatch((int)id, threadGroups.x, threadGroups.y, threadGroups.z);
                Graphics.ClearRandomWriteTargets();
            }

            private void DispatchKernel(KernelID id, RenderTexture writeTarget)
            {
                DispatchKernel(id, new RenderTexture[] { writeTarget });
            }
        }
    }
}