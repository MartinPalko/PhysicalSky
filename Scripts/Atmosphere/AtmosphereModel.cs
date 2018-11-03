using System.Collections.Generic;
using UnityEngine;
#if UNITY_EDITOR
using UnityEditor;
#endif
using System;

using PhysicalSky.Utilities;

namespace PhysicalSky
{
    [CreateAssetMenu(fileName = "NewAtmosphereModel", menuName = "PhysicalSky/AtmosphereModel")]
    public partial class AtmosphereModel : ScriptableObject, IAtmosphereModel
    {
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
        [SerializeField] [HideInInspector]
        private AtmosphereParameters m_computedParameters;
        public AtmosphereParameters ComputedParameters
        {
            get { return m_computedParameters; }
        }

        [SerializeField] [HideInInspector]
        private AtmosphereLib.AtmosphereRenderParams m_renderParams;

        // Temporary lookup render textures output from the computation
        private RenderTexture m_gpuTransmittanceLUT = null;
        private RenderTexture m_gpuScatteringLUT = null;
        private RenderTexture m_gpuIrradianceLUT = null;

        // CPU accessible transmittance lookup texture
        private Texture2D m_cpuTransmittanceLUT = null;
        private Texture2D GetCPUTransmittanceLUT()
        {
            if (m_savedTransmittanceLUT)
            {
                return m_savedTransmittanceLUT;
            }
            else if (m_cpuTransmittanceLUT)
            {
                return m_cpuTransmittanceLUT;
            }
            else if (m_gpuTransmittanceLUT)
            {
                m_cpuTransmittanceLUT = Utilities.GraphicsHelpers.Texture2DFromRenderTexture(m_gpuTransmittanceLUT);
                return m_cpuTransmittanceLUT;
            }
            else
            {
                return null;
            }
        }

        private enum SerializeTextureMethod
        {
            SubAsset = 0,
            SeperateAsset = 1,
            DontSerialize = 2,
        }

        [SerializeField]
        private SerializeTextureMethod m_serializeTextures = SerializeTextureMethod.SubAsset;

        [SerializeField] [HideInInspector]
        private Texture2D m_savedTransmittanceLUT = null;
        [SerializeField] [HideInInspector]
        private Texture2D m_savedScatteringLUT = null;
        [SerializeField] [HideInInspector]
        private Texture2D m_savedIrradianceLUT = null;

#if UNITY_EDITOR
        private void DeleteTexture(ref Texture2D texture)
        {
            if (m_savedTransmittanceLUT)
            {
                string transmittanceAssetPath = AssetDatabase.GetAssetPath(m_savedTransmittanceLUT);

                if (transmittanceAssetPath != AssetDatabase.GetAssetPath(this))
                    AssetDatabase.DeleteAsset(transmittanceAssetPath);

                DestroyImmediate(m_savedTransmittanceLUT, true);
                m_savedTransmittanceLUT = null;
            }
        }

        private Texture2D SaveTexture(RenderTexture texture, string name)
        {
            Texture2D savedTexture = GraphicsHelpers.Texture2DFromRenderTexture(texture);
            savedTexture.name = name;

            string atmosphereAssetPath = AssetDatabase.GetAssetPath(this);
            string atmosphereAssetPathNoExtension = System.IO.Path.ChangeExtension(atmosphereAssetPath, null);

            if (m_serializeTextures == SerializeTextureMethod.SubAsset)
            {
                AssetDatabase.AddObjectToAsset(savedTexture, atmosphereAssetPath);
            }
            else if (m_serializeTextures == SerializeTextureMethod.SeperateAsset)
            {
                AssetDatabase.CreateAsset(savedTexture, atmosphereAssetPathNoExtension + "." + name + ".asset");
            }
            else
            {
                throw new System.NotImplementedException();
            }

            return savedTexture;
        }

        private void UpdateSavedTextures()
        {
            DeleteTexture(ref m_savedTransmittanceLUT);
            DeleteTexture(ref m_savedScatteringLUT);
            DeleteTexture(ref m_savedIrradianceLUT);

            // Clean up just incase any sub-assets were left over
            foreach (UnityEngine.Object subAsset in AssetDatabase.LoadAllAssetsAtPath(AssetDatabase.GetAssetPath(this)))
            {
                if (subAsset != this)
                    DestroyImmediate(subAsset, true);
            }

            if (m_serializeTextures != SerializeTextureMethod.DontSerialize)
            {
                m_savedTransmittanceLUT = SaveTexture(m_gpuTransmittanceLUT, "TransmittanceLUT");
                m_savedScatteringLUT = SaveTexture(m_gpuScatteringLUT, "ScatteringLUT");
                m_savedIrradianceLUT = SaveTexture(m_gpuIrradianceLUT, "IrradianceLUT");
            }

            AssetDatabase.SaveAssets();
            AssetDatabase.Refresh();
        }
#endif

        public bool HasComputedData()
        {
            return !TexturesInvalid() && m_computedParameters != new AtmosphereParameters();
        }

        public bool NeedsRecompute()
        {
            return (m_computedParameters != m_parameters) || !HasComputedData();
        }

        private void SetShaderUniforms(GraphicsHelpers.IMaterialProperties m)
        {
            // TODO: Set only run-time required values here
            // TODO: Set all from RenderValues
            m.SetTexture("_transmittance_texture", m_gpuTransmittanceLUT != null ? (Texture)m_gpuTransmittanceLUT : m_savedTransmittanceLUT);
            m.SetTexture("_scattering_texture", m_gpuScatteringLUT != null ? (Texture)m_gpuScatteringLUT : m_savedScatteringLUT);
            m.SetTexture("_irradiance_texture", m_gpuIrradianceLUT != null ? (Texture)m_gpuIrradianceLUT : m_savedIrradianceLUT);
            m.SetVector("_sky_spectral_radiance_to_luminance", m_renderParams.sky_spectral_radiance_to_luminance);
            m.SetVector("_sun_spectral_radiance_to_luminance", m_renderParams.sun_spectral_radiance_to_luminance);
            m.SetVector("_solar_irradiance", m_renderParams.solar_irradiance);
            m.SetFloat("_sun_angular_radius", m_computedParameters.sunAngularRadius);
            m.SetFloat("_bottom_radius", m_computedParameters.planetaryRadius / Prerenderer.LENGTH_UNIT_IN_METERS);
            m.SetFloat("_top_radius", (m_computedParameters.planetaryRadius + m_computedParameters.atmosphereThickness) / Prerenderer.LENGTH_UNIT_IN_METERS);
            m.SetVector("_rayleigh_scattering", m_renderParams.rayleigh_scattering);
            m.SetVector("_mie_scattering", m_renderParams.mie_scattering);
            m.SetVector("_mie_extinction", m_renderParams.mie_extinction);
            m.SetFloat("_mie_phase_function_g", m_computedParameters.miePhaseFunctionG);
            m.SetVector("_absorption_extinction", m_renderParams.absorption_extinction);
            m.SetFloat("_mu_s_min", Mathf.Cos(m_computedParameters.maxSunZenithAngle));
            m.SetVector("sun_size", new Vector3(Mathf.Tan(m_computedParameters.sunAngularRadius), Mathf.Cos(m_computedParameters.sunAngularRadius), m_computedParameters.sunAngularRadius));
        }

        public void SetShaderUniforms(Material material)
        {
            SetShaderUniforms(material.ToMaterialPropertyInterface());
        }

        public void SetShaderUniforms(MaterialPropertyBlock propertyBlock)
        {
            SetShaderUniforms(propertyBlock.ToMaterialPropertyInterface());
        }

        public bool Compute(bool force = false)
        {
            if (NeedsRecompute() || force)
            {
                AllocateRenderTextures();

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

                if (preRenderer.Compute(m_parameters, m_gpuTransmittanceLUT, m_gpuScatteringLUT, m_gpuIrradianceLUT, ref m_renderParams))
                {
                    m_computedParameters = m_parameters;
#if UNITY_EDITOR
                    UpdateSavedTextures();

                    if (m_serializeTextures != SerializeTextureMethod.DontSerialize)
                    {
                        ReleaseTextures();
                    }

                    EditorUtility.SetDirty(this);
#endif
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
            ReleaseTextures();
            m_computedParameters = new AtmosphereParameters();
        }

        private void AllocateRenderTextures()
        {
            TextureFactory.CreateRenderTexture(ref m_gpuTransmittanceLUT, TextureFactory.Preset.Transmittance);
            TextureFactory.CreateRenderTexture(ref m_gpuIrradianceLUT, TextureFactory.Preset.Irradiance);
            TextureFactory.CreateRenderTexture(ref m_gpuScatteringLUT, TextureFactory.Preset.Scattering);
        }

        private void ReleaseTextures()
        {
            TextureFactory.ReleaseRenderTexture(ref m_gpuTransmittanceLUT);
            TextureFactory.ReleaseRenderTexture(ref m_gpuIrradianceLUT);
            TextureFactory.ReleaseRenderTexture(ref m_gpuScatteringLUT);

            if (m_cpuTransmittanceLUT)
            {
                Destroy(m_cpuTransmittanceLUT);
                m_cpuTransmittanceLUT = null;
            }
        }

        public bool TexturesInvalid()
        {
            if (m_savedIrradianceLUT == null && (m_gpuTransmittanceLUT == null || !m_gpuTransmittanceLUT.IsCreated()))
                return true;
            else if (m_savedScatteringLUT == null && (m_gpuScatteringLUT == null || !m_gpuScatteringLUT.IsCreated()))
                return true;
            else if (m_savedIrradianceLUT == null && (m_gpuIrradianceLUT == null || !m_gpuIrradianceLUT.IsCreated()))
                return true;
            else
                return false;
        }

        private void OnDestroy()
        {
            ReleaseResources();
        }

        public Color GetSunIlluminance(Vector3 positionRelativePlanetCenter, Vector3 sunDirection)
        {
            if (HasComputedData())
                return AtmosphereLib.GetSunIlluminance(m_renderParams, GetCPUTransmittanceLUT(), positionRelativePlanetCenter, sunDirection);
            else
                return Color.black;
        }
    }
}