using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace PhysicalSky
{
    public partial class AtmosphereModel
    {
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

        private sealed class PrerendererBlit : Prerenderer
        {
            Material m_material;

            public PrerendererBlit(Shader shader) : base()
            {
                m_material = new Material(shader);
            }

            public override bool Supported()
            {
                return SystemInfo.graphicsShaderLevel >= 50;
            }

            protected override void SetupPrerender(Vector3 lambdas, Matrix3x3 luminanceFromRadiance)
            {
                m_material.SetFloatArray("_luminance_from_radiance", luminanceFromRadiance);

                m_material.SetTexture("transmittance_texture", m_transmittanceLUT);
                m_material.SetTexture("single_rayleigh_scattering_texture", m_tempTextures.DeltaRayleighScatteringTexture);
                m_material.SetTexture("single_mie_scattering_texture", m_tempTextures.DeltaMieScatteringTexture);
                m_material.SetTexture("multiple_scattering_texture", m_tempTextures.DeltaMultipleScatteringTexture);
                m_material.SetTexture("irradiance_texture", m_tempTextures.DeltaIrradianceTexture);
                m_material.SetTexture("scattering_density_texture", m_tempTextures.DeltaScatteringDensityTexture);

                m_material.SetVector("_sky_spectral_radiance_to_luminance", m_parameters.luminance != AtmosphereParameters.LuminanceType.none ? m_sky_k : Vector3.one);
                m_material.SetVector("_sun_spectral_radiance_to_luminance", m_parameters.luminance != AtmosphereParameters.LuminanceType.none ? m_sun_k : Vector3.one);

                m_material.SetVector("_solar_irradiance", ScaleToWavelengths(m_solarIrradiance, lambdas, 1.0f));
                m_material.SetFloat("_sun_angular_radius", m_parameters.sunAngularRadius);
                m_material.SetFloat("_bottom_radius", m_parameters.planetaryRadius / LENGTH_UNIT_IN_METERS);
                m_material.SetFloat("_top_radius", (m_parameters.planetaryRadius + m_parameters.atmosphereThickness) / LENGTH_UNIT_IN_METERS);
                m_material.SetFloatArray("_rayleigh_density0", m_rayleighDensity.layer0.GetShaderValues());
                m_material.SetFloatArray("_rayleigh_density1", m_rayleighDensity.layer1.GetShaderValues());
                m_material.SetVector("_rayleigh_scattering", ScaleToWavelengths(m_rayleighScattering, lambdas, LENGTH_UNIT_IN_METERS));
                m_material.SetFloatArray("_mie_density0", m_mieDensity.layer0.GetShaderValues());
                m_material.SetFloatArray("_mie_density1", m_mieDensity.layer1.GetShaderValues());
                m_material.SetVector("_mie_scattering", ScaleToWavelengths(m_mieScattering, lambdas, LENGTH_UNIT_IN_METERS));
                m_material.SetVector("_mie_extinction", ScaleToWavelengths(m_mieExtinction, lambdas, LENGTH_UNIT_IN_METERS));
                m_material.SetFloat("_mie_phase_function_g", m_parameters.miePhaseFunctionG);
                m_material.SetFloatArray("_absorption_density0", m_absorptionDensity.layer0.GetShaderValues());
                m_material.SetFloatArray("_absorption_density1", m_absorptionDensity.layer1.GetShaderValues());
                m_material.SetVector("_absorption_extinction", ScaleToWavelengths(m_absorptionExtinction, lambdas, LENGTH_UNIT_IN_METERS));
                m_material.SetVector("_ground_albedo", Vector3.one * m_parameters.groundAlbedo);
                m_material.SetFloat("_mu_s_min", Mathf.Cos(m_parameters.maxSunZenithAngle));
            }

            protected override void CleanupPrerenderer()
            {
                RenderTexture.active = null;
            }

            protected override void RenderTransmittanceLUT()
            {
                // Compute Transmittance LUT.
                Utilities.GraphicsHelpers.Blit3D(m_transmittanceLUT, m_material, (int)PrecomputePass.ComputeTransmittance);
                Utilities.GraphicsHelpers.DumpTexture(m_transmittanceLUT, "C:\\debug\\m_transmittanceLUT");
            }

            protected override void RenderSingleScattering(bool accumulate)
            {
                // Compute Direct Irradiance into DeltaIrradianceTexture and either initialize m_irradianceLUT with 0 or leave it unchanged. (We don't want the direct irradiance in m_irradianceLUT, but only the irradiance from the sky)
                if (accumulate)
                    Utilities.GraphicsHelpers.Blit3D(m_tempTextures.DeltaIrradianceTexture, m_material, (int)PrecomputePass.ComputeDirectIrradiance);
                else
                    Utilities.GraphicsHelpers.Blit3D(new RenderTexture[2] { m_tempTextures.DeltaIrradianceTexture, m_irradianceLUT }, m_material, (int)PrecomputePass.ComputeDirectIrradiance);
                
                // Compute Raylie and Mie Single Scattering, and store them in DeltaRayleighScatteringTexture, DeltaMieScatteringTexture as well as scatteringLUT
                if (accumulate)
                    Utilities.GraphicsHelpers.Blit3D(new RenderTexture[2] { m_tempTextures.DeltaRayleighScatteringTexture, m_tempTextures.DeltaMieScatteringTexture }, m_material, (int)PrecomputePass.ComputeSingleScattering);
                else
                    Utilities.GraphicsHelpers.Blit3D(new RenderTexture[3] { m_tempTextures.DeltaRayleighScatteringTexture, m_tempTextures.DeltaMieScatteringTexture, m_scatteringLUT }, m_material, (int)PrecomputePass.ComputeSingleScattering);
                
                Utilities.GraphicsHelpers.Blit3D(m_scatteringLUT, m_material, (int)PrecomputePass.AccumulateSingleScattering);
            }

            protected override void RenderMultipleScattering(int scatteringOrder)
            {
                // Compute the scattering density, and store it in DeltaScatteringDensityTexture.
                m_material.SetInt("scattering_order", scatteringOrder);
                Utilities.GraphicsHelpers.Blit3D(m_tempTextures.DeltaScatteringDensityTexture, m_material, (int)PrecomputePass.ComputeScatteringDensity);

                // Compute the indirect irradiance, store it in DeltaIrradianceTexture and accumulate it in irradianceLUT.
                m_material.SetInt("scattering_order", scatteringOrder - 1);
                Utilities.GraphicsHelpers.Blit3D(m_tempTextures.DeltaIrradianceTexture, m_material, (int)PrecomputePass.ComputeIndirectIrradiance);
                Utilities.GraphicsHelpers.Blit3D(m_irradianceLUT, m_material, (int)PrecomputePass.AccumulateIndirectIrradiance);

                // Compute the multiple scattering, store it in DeltaMultipleScatteringTexture, and accumulate it in scatteringLUT.
                Utilities.GraphicsHelpers.Blit3D(m_tempTextures.DeltaMultipleScatteringTexture, m_material, (int)PrecomputePass.ComputeMultipleScattering);
                Utilities.GraphicsHelpers.Blit3D(m_scatteringLUT, m_material, (int)PrecomputePass.AccumulateMultipleScattering);
            }
        }
    }
}