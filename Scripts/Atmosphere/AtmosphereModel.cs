using System.Collections.Generic;
using UnityEngine;
using System;

namespace PhysicalSky
{
    [CreateAssetMenu(fileName = "NewAtmosphereModel", menuName = "PhysicalSky/AtmosphereModel")]
    public partial class AtmosphereModel : ScriptableObject, IAtmosphereModel
    {
        [Serializable]
        private struct ComputedRuntimeValues
        {
            public float m_sunSolidAngle;
            public Vector3 m_sky_k;
            public Vector3 m_sun_k;
            public Vector3 m_solarIrradiance;
            public Vector3 m_rayleighScattering;
            public Vector3 m_mieScattering;
            public Vector3 m_mieExtinction;
            public Vector3 m_absorptionExtinction;
        }

        [SerializeField]
        private AtmosphereParameters m_parameters = AtmosphereParameters.defaultEarth;
        public AtmosphereParameters Parameters
        {
            get { return m_parameters; }
            set
            {
                if (m_parameters != value)
                {
                    m_parameters = value;
                }
            }
        }

        // Computed values
        private AtmosphereParameters m_computedParameters;
        private ComputedRuntimeValues m_computedValues;
        // Computed textures
        private RenderTexture m_transmittanceLUT = null;
        private RenderTexture m_scatteringLUT = null;
        private RenderTexture m_irradianceLUT = null;

        public bool NeedsRecompute()
        {
            return (m_computedParameters != m_parameters) || TexturesInvalid();
        }

        public void SetShaderUniforms(Material material)
        {
            // TODO: Set only run-time required values here
            material.SetTexture("_transmittance_texture", m_transmittanceLUT);
            material.SetTexture("_scattering_texture", m_scatteringLUT);
            material.SetTexture("_irradiance_texture", m_irradianceLUT);
            
            material.SetVector("_sky_spectral_radiance_to_luminance", m_computedParameters.luminance != AtmosphereParameters.LuminanceType.none ? m_computedValues.m_sky_k : Vector3.one);
            material.SetVector("_sun_spectral_radiance_to_luminance", m_computedParameters.luminance != AtmosphereParameters.LuminanceType.none ? m_computedValues.m_sun_k : Vector3.one);

            material.SetVector("_solar_irradiance", m_computedValues.m_solarIrradiance);
            material.SetFloat("_sun_angular_radius", m_computedParameters.sunAngularRadius);
            material.SetFloat("_bottom_radius", m_computedParameters.planetaryRadius / Prerenderer.LENGTH_UNIT_IN_METERS);
            material.SetFloat("_top_radius", (m_computedParameters.planetaryRadius + m_computedParameters.atmosphereThickness) / Prerenderer.LENGTH_UNIT_IN_METERS);
            material.SetVector("_rayleigh_scattering", m_computedValues.m_rayleighScattering);
            material.SetVector("_mie_scattering", m_computedValues.m_mieScattering);
            material.SetVector("_mie_extinction", m_computedValues.m_mieExtinction);
            material.SetFloat("_mie_phase_function_g", m_computedParameters.miePhaseFunctionG);
            material.SetVector("_absorption_extinction", m_computedValues.m_absorptionExtinction);
            material.SetFloat("_mu_s_min", Mathf.Cos(m_computedParameters.maxSunZenithAngle));

            material.SetVector("sun_size", new Vector3(Mathf.Tan(m_computedParameters.sunAngularRadius), Mathf.Cos(m_computedParameters.sunAngularRadius), m_computedParameters.sunAngularRadius));
        }

        public bool Compute(bool force = false)
        {
            if (NeedsRecompute() || force)
            {
                AllocateLookupTextures();

                Prerenderer preRenderer;

                if (PrerendererCompute.Supported())
                    preRenderer = new PrerendererCompute();
                else if (PrerendererBlit.Supported())
                    preRenderer = new PrerendererBlit();
                else
                {
                    Debug.LogError("No atmosphere model pre-renderer is supported on the current platform. Cannot compute lookup textures!");
                    return false;
                }

                if (preRenderer.Compute(m_parameters, m_transmittanceLUT, m_scatteringLUT, m_irradianceLUT, ref m_computedValues))
                {
                    m_computedParameters = m_parameters;
                    return true;
                }
                else
                {
                    return false;
                }
            }
            return false;
        }

        public void ReleaseResources()
        {
#if PHYSICAL_SKY_DEBUG
            Debug.Log("Released Atmosphere Resources");
#endif
            ReleaseLookupTextures();
            m_computedParameters = new AtmosphereParameters();
        }

        private void AllocateLookupTextures()
        {
            TextureFactory.CreateRenderTexture(ref m_transmittanceLUT, TextureFactory.Preset.Transmittance);
            TextureFactory.CreateRenderTexture(ref m_irradianceLUT, TextureFactory.Preset.Irradiance);
            TextureFactory.CreateRenderTexture(ref m_scatteringLUT, TextureFactory.Preset.Scattering);
        }

        private void ReleaseLookupTextures()
        {
            TextureFactory.ReleaseRenderTexture(ref m_transmittanceLUT);
            TextureFactory.ReleaseRenderTexture(ref m_irradianceLUT);
            TextureFactory.ReleaseRenderTexture(ref m_scatteringLUT);
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

        private void OnDestroy()
        {
            ReleaseResources();
        }
    }
}