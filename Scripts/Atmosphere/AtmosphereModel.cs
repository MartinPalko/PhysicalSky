using System.Collections.Generic;
using UnityEngine;
using System;

namespace PhysicalSky
{
    [CreateAssetMenu(fileName = "NewAtmosphereModel", menuName = "PhysicalSky/AtmosphereModel")]
    public partial class AtmosphereModel : ScriptableObject, IAtmosphereModel
    {
        // Matches that in PhysicalSkyCommon.cginc
        public struct DensityProfileLayer
        {
            public float width;
            public float exp_term;
            public float exp_scale;
            public float linear_term;
            public float constant_term;

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
        };

        public struct DensityProfile
        {
            public DensityProfileLayer layer0;
            public DensityProfileLayer layer1;

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
        }

        // Computed values
        private float m_sunSolidAngle;
        Vector3 m_sky_k;
        Vector3 m_sun_k;
        private List<float> m_wavelengths = new List<float>();
        private List<float> m_solarIrradiance = new List<float>();
        private List<float> m_rayleighScattering = new List<float>();
        private List<float> m_mieScattering = new List<float>();
        private List<float> m_mieExtinction = new List<float>();
        private List<float> m_absorptionExtinction = new List<float>();
        private DensityProfile m_rayleighDensity;
        private DensityProfile m_mieDensity;
        private DensityProfile m_absorptionDensity;

        // Computed textures
        private RenderTexture m_transmittanceLUT = null;
        private RenderTexture m_scatteringLUT = null;
        private RenderTexture m_irradianceLUT = null;
        
        public bool NeedsRecompute { get { return (m_computedParameters != m_parameters) || TexturesInvalid(); } }

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
        private AtmosphereParameters m_computedParameters;

        // Shader and material used for precompute
        [SerializeField]
        private Shader m_PrecomputeShader = null;
        public Shader PrecomputeShader { get { return m_PrecomputeShader; } set { m_PrecomputeShader = value; } }

        public void SetShaderUniforms(Material m)
        {
            Internal.SetShaderUniforms(this, m_computedParameters, m);
        }

        public bool Compute(bool force = false)
        {
            if (NeedsRecompute || force)
            {
                AllocateLookupTextures();
                Internal.Compute(this, m_parameters);
                m_computedParameters = m_parameters;
                return true;
            }
            return false;
        }

        private void AllocateLookupTextures()
        {
            Internal.TextureFactory.CreateRenderTexture(ref m_transmittanceLUT, Internal.TextureFactory.Preset.Transmittance);
            Internal.TextureFactory.CreateRenderTexture(ref m_irradianceLUT, Internal.TextureFactory.Preset.Irradiance);
            Internal.TextureFactory.CreateRenderTexture(ref m_scatteringLUT, Internal.TextureFactory.Preset.Scattering);
        }

        private void ReleaseLookupTextures()
        {
            Internal.TextureFactory.ReleaseRenderTexture(ref m_transmittanceLUT);
            Internal.TextureFactory.ReleaseRenderTexture(ref m_irradianceLUT);
            Internal.TextureFactory.ReleaseRenderTexture(ref m_scatteringLUT);
        }

        public void ReleaseResources()
        {
#if PHYSICAL_SKY_DEBUG
            Debug.Log("Released Atmosphere Resources");
#endif
            ReleaseLookupTextures();
            m_computedParameters = new AtmosphereParameters();
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