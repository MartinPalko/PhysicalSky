﻿using UnityEngine;

using PhysicalSky.Utilities;

namespace PhysicalSky
{
    [ExecuteInEditMode]
    [RequireComponent(typeof(Light))]
    public class PhysicalSky : MonoBehaviour, IPhysicalSky
    {
        public Vector3 SunDirection { get { return -transform.forward; } set { transform.rotation = Quaternion.LookRotation(-value); } }
        public Quaternion StarRotation { get { return m_StarMeshObject ? m_StarMeshObject.transform.rotation : Quaternion.identity; } set { if (m_StarMeshObject) m_StarMeshObject.transform.rotation = value; } }
        
        [SerializeField]
        private float m_SunLightBrightnessMultiplier = 1.0f;
        public float SunLightBrightnessMultiplier { get { return m_SunLightBrightnessMultiplier; } set { m_SunLightBrightnessMultiplier = Mathf.Max(value, 0.0f); } }

        [SerializeField]
        private float m_SunBrightnessMultiplier = 1.0f;
        public float SunBrightnessMultiplier { get { return m_SunBrightnessMultiplier; } set { m_SunBrightnessMultiplier = Mathf.Max(value, 0.0f); } }

        [SerializeField]
        private float m_SunSaturationMultiplier = 1.0f;
        public float SunSaturationMultiplier { get { return m_SunSaturationMultiplier; } set { m_SunSaturationMultiplier = Mathf.Max(value, 0.0f); } }

        [SerializeField]
        private float m_SkyExposure = 20.0f;
        public float SkyExposure { get { return m_SkyExposure; } set { m_SkyExposure = Mathf.Max(value, 0.0f); } }

        [SerializeField]
        private float m_Altitude = 0.1f;
        public float Altitude { get { return m_Altitude; } set { m_Altitude = Mathf.Max(value, 0.0f); DynamicGI.UpdateEnvironment(); } }

        [SerializeField]
        private AtmosphereModel m_Atmosphere;
        public IAtmosphereModel Atmosphere { get { return m_Atmosphere; } set { m_Atmosphere = value as AtmosphereModel; } }

        [SerializeField]
        private bool m_AutoRecomputeAtmosphere = true;

        [SerializeField]
        private float m_StarBrightnessMultiplier = 1.0f;
        public float StarBrightnessMultiplier { get { return m_StarBrightnessMultiplier; } set { m_StarBrightnessMultiplier = Mathf.Max(value, 0.0f); } }

        [SerializeField]
        private float m_StarBrightnessPower = 0.454545f;
        public float StarBrightnessPower { get { return m_StarBrightnessPower; } set { m_StarBrightnessPower = Mathf.Max(value, 0.0f); } }

        [SerializeField]
        private StarMap m_StarMap;
        public StarMap StarMap { get { return m_StarMap; } set { m_StarMap = value; } }

        [SerializeField]
        private Shader m_SkyShader = null;
        public Shader SkyShader { get { return m_SkyShader; } set { m_SkyShader = value; } }

        [SerializeField]
        private Shader m_SunRadianceShader = null;
        public Shader SunRadianceShader { get { return m_SunRadianceShader; } set { m_SunRadianceShader = value; } }

        [SerializeField]
        private Shader m_StarMeshShader = null;
        public Shader StarMeshShader { get { return m_StarMeshShader; } set { m_StarMeshShader = value; } }

        private Material m_SkyMaterial = null;
        private Material m_SunRadianceMaterial = null;

        private Material m_StarMeshMaterial = null;
        private GameObject m_StarMeshObject = null;
        private MeshFilter m_StarMeshFilter = null;

        private const int SUN_RADIANCE_TEXTURE_SIZE = 8;
        private Light m_SunLight = null;
        private Texture2D m_SunRadianceTexture = null;

        private void Start()
        {
            m_SkyMaterial = new Material(m_SkyShader);
            m_SunRadianceMaterial = new Material(m_SunRadianceShader);

            if (StarMap)
            {
                m_StarMeshMaterial = new Material(m_StarMeshShader);

                m_StarMeshObject = new GameObject("StarMesh");
#if PHYSICAL_SKY_DEBUG
                m_StarMeshObject.hideFlags = HideFlags.DontSave;
#else
                m_StarMeshObject.hideFlags = HideFlags.HideAndDontSave;
#endif
                m_StarMeshFilter = m_StarMeshObject.AddComponent<MeshFilter>();
                m_StarMeshFilter.mesh = StarMap.CreateStarMesh();
                MeshRenderer starMeshRenderer = m_StarMeshObject.AddComponent<MeshRenderer>();
                starMeshRenderer.material = m_StarMeshMaterial;
            }
            
            m_SunLight = GetComponent<Light>();
            m_SunLight.type = LightType.Directional;
        }

        private void OnDestroy()
        {
#if UNITY_EDITOR
            if (!Application.isPlaying)
            {
                if (m_StarMeshObject)
                    DestroyImmediate(m_StarMeshObject);
            }
            else
#endif
            {
                if (m_StarMeshObject)
                    Destroy(m_StarMeshObject);
            }
        }

        private void SetShaderUniforms(GraphicsHelpers.IMaterialProperties m)
        {
            m.SetVector("sun_direction", -transform.forward);

            m.SetFloat("sky_exposure", SkyExposure);

            // TEMP/HACK: When using luminance, values are much brighter than radiance.
            // Need a better solution to handle that.
            float luminanceCompensation = m_Atmosphere.ComputedParameters.luminance != AtmosphereParameters.LuminanceType.none ? 1e-05f : 1.0f;
            luminanceCompensation /= Mathf.PI; // TODO: Nicer conversion from illuminance to luminance;
            m.SetFloat("sun_brightness", SunBrightnessMultiplier * luminanceCompensation);

            m.SetVector("camera", new Vector3(0, (m_Atmosphere.ComputedParameters.planetaryRadius / 1000) + m_Altitude, 0));
        }

        public void SetShaderUniforms(Material material)
        {
            m_Atmosphere.SetShaderUniforms(material);
            SetShaderUniforms(material.ToMaterialPropertyInterface());
        }

        public void SetShaderUniforms(MaterialPropertyBlock materialPropertyBlock)
        {
            m_Atmosphere.SetShaderUniforms(materialPropertyBlock);
            SetShaderUniforms(materialPropertyBlock.ToMaterialPropertyInterface());
        }

        private void LateUpdate()
        {
#if UNITY_EDITOR
            if (!Application.isPlaying)
            {
                m_SunLight = GetComponent<Light>();
            }
#endif

            if (!m_Atmosphere)
            {
                Debug.LogWarning("No atmosphere");
                return;
            }

            if (m_AutoRecomputeAtmosphere && m_Atmosphere.NeedsRecompute())
                m_Atmosphere.Compute();

            if (!m_SunRadianceTexture)
            {
                m_SunRadianceTexture = new Texture2D(SUN_RADIANCE_TEXTURE_SIZE, SUN_RADIANCE_TEXTURE_SIZE, TextureFormat.RGBAHalf, false, false);
            }
            
            {
                float radius = m_Atmosphere.ComputedParameters.planetaryRadius / 1000;
                float altitude = m_Altitude;
                Color sunLuminance = Atmosphere.GetSunIlluminance(new Vector3(0, radius + altitude, 0), -m_SunLight.transform.forward);
                float hue, saturation, value;
                Color.RGBToHSV(sunLuminance, out hue, out saturation, out value);
                m_SunLight.intensity = value * (m_Atmosphere.ComputedParameters.luminance != AtmosphereParameters.LuminanceType.none ? 1e-05f : 1.0f);
                m_SunLight.color = Color.HSVToRGB(hue, saturation * SunSaturationMultiplier, 1.0f).gamma; // TODO: Why does this need the .gamma to seem right?
            }

            // Assign sky material as skybox.
            if (m_SkyMaterial)
            {
                SetShaderUniforms(m_SkyMaterial);
                

                if (m_StarMap)
                {
                    m_SkyMaterial.SetTexture("star_cubemap", m_StarMap.BackgroundCube);
                    m_SkyMaterial.SetFloat("star_brightness", m_StarMap.BackgroundCubeBrightness * m_StarBrightnessMultiplier);
                    m_SkyMaterial.SetMatrix("star_rotation", Matrix4x4.TRS(Vector3.zero, StarRotation * Quaternion.Euler(m_StarMap.BackgroundRotation), Vector3.one));
                }
                else
                {
                    m_SkyMaterial.SetTexture("star_cubemap", null);
                    m_SkyMaterial.SetFloat("star_brightness", 0.0f);
                    m_SkyMaterial.SetMatrix("star_rotation", Matrix4x4.identity);
                }

                RenderSettings.skybox = m_SkyMaterial;
                RenderSettings.sun = m_SunLight;
            }

            if (m_StarMeshMaterial)
            {
                SetShaderUniforms(m_StarMeshMaterial);

                m_StarMeshMaterial.SetFloat("star_intensity_multiplier", m_StarBrightnessMultiplier * SkyExposure);
                m_StarMeshMaterial.SetFloat("star_intensity_power", m_StarBrightnessPower);
                
                float D = m_Atmosphere.ComputedParameters.planetaryRadius + m_Altitude * 1000;
                float R = m_Atmosphere.ComputedParameters.planetaryRadius;
                float planetAngularRad = 2 * Mathf.Acos(Mathf.Sqrt(Mathf.Pow(D, 2) - Mathf.Pow(R, 2)) / D);
                m_StarMeshMaterial.SetFloat("planet_size", planetAngularRad);              
            }

#if UNITY_EDITOR
            // Force update of realtime reflection probes while in edit mode.
            if (!Application.isPlaying)
            {
                ReflectionProbe[] reflectionProbes = FindObjectsOfType<ReflectionProbe>();
                for (int i = 0; i < reflectionProbes.Length; i++)
                {
                    ReflectionProbe probe = reflectionProbes[i];
                    if (probe.mode == UnityEngine.Rendering.ReflectionProbeMode.Realtime)
                    {
                        probe.RenderProbe();
                    }
                }
            }
#endif
        }
    }
}