using UnityEngine;

namespace PhysicalSky
{
    [ExecuteInEditMode]
    [RequireComponent(typeof(Light))]
    public class PhysicalSky : MonoBehaviour, IPhysicalSky
    {
        public Vector3 SunDirection { get { return -transform.forward; } set { transform.rotation = Quaternion.LookRotation(-value); } }
        public Quaternion StarRotation { get { return m_StarMeshObject.transform.rotation; } set { m_StarMeshObject.transform.rotation = value; } }
        
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
            m_StarMeshMaterial = new Material(m_StarMeshShader);

            m_StarMeshObject = new GameObject("StarMesh");
#if PHYSICAL_SKY_DEBUG
            m_StarMeshObject.hideFlags = HideFlags.DontSave;
#else
            starMeshObject.hideFlags = HideFlags.HideAndDontSave;
#endif
            m_StarMeshFilter = m_StarMeshObject.AddComponent<MeshFilter>();
            m_StarMeshFilter.mesh = StarMap.CreateStarMesh();     
            MeshRenderer starMeshRenderer = m_StarMeshObject.AddComponent<MeshRenderer>();
            starMeshRenderer.material = m_StarMeshMaterial;
            
            m_SunLight = GetComponent<Light>();
            m_SunLight.type = LightType.Directional;
        }

        private void OnDestroy()
        {
#if UNITY_EDITOR
            if (!Application.isPlaying)
            {
                DestroyImmediate(m_StarMeshObject);
            }
            else
#endif
            {
                Destroy(m_StarMeshObject);
            }
        }

        private void ConfigureMaterial(Material m)
        {
            m_Atmosphere.SetShaderUniforms(m);
            m.SetTexture("transmittance_texture", m_Atmosphere.TransmittanceLUT);
            m.SetTexture("scattering_texture", m_Atmosphere.ScatteringLUT);
            m.SetTexture("irradiance_texture", m_Atmosphere.IrradianceLUT);

            m.SetVector("camera", new Vector3(0, (m_Atmosphere.PlanetaryRadius / 1000) + m_Altitude, 0));
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

            if (m_Atmosphere.NeedsRecompute)
                m_Atmosphere.Compute();

            if (!m_SunRadianceTexture)
            {
                m_SunRadianceTexture = new Texture2D(SUN_RADIANCE_TEXTURE_SIZE, SUN_RADIANCE_TEXTURE_SIZE, TextureFormat.RGBAHalf, false, false);
            }

            // Compute sun radiance, and apply to light
            if (m_SunRadianceMaterial)
            {
                ConfigureMaterial(m_SunRadianceMaterial);
                m_SunRadianceMaterial.SetMatrix("sun_rotation_matrix", transform.localToWorldMatrix);
                m_SunRadianceMaterial.SetVector("sun_direction", -transform.forward);

                RenderTexture temporary = RenderTexture.GetTemporary(SUN_RADIANCE_TEXTURE_SIZE, SUN_RADIANCE_TEXTURE_SIZE, 0, RenderTextureFormat.ARGBHalf);
                RenderTexture.active = temporary;
                Graphics.Blit(null, temporary, m_SunRadianceMaterial);
                m_SunRadianceTexture.ReadPixels(new Rect(0, 0, SUN_RADIANCE_TEXTURE_SIZE, SUN_RADIANCE_TEXTURE_SIZE), 0, 0);
                RenderTexture.active = null;
                RenderTexture.ReleaseTemporary(temporary);

                Color[] colors = m_SunRadianceTexture.GetPixels();

                Color radianceSum = Color.black;
                int nonZero = 0;

                for (int i = 0; i < colors.Length; i++)
                {
                    radianceSum += colors[i];
                    nonZero += colors[i] == Color.black ? 0 : 1;
                }
                Color averageRadiance = radianceSum / nonZero;
                averageRadiance /= 3000; // Arbitrary value scales down to unity's lighting.
                averageRadiance *= SunLightBrightnessMultiplier;

                // Set directional sunlight based on average radiance.
                {
                    float hue, saturation, value;
                    Color.RGBToHSV(averageRadiance, out hue, out saturation, out value);
                    m_SunLight.intensity = value;
                    m_SunLight.color = Color.HSVToRGB(hue, saturation * SunSaturationMultiplier, 1.0f).gamma;
                }
            }

            // Assign sky material as skybox.
            if (m_SkyMaterial)
            {
                ConfigureMaterial(m_SkyMaterial);
                m_SkyMaterial.SetFloat("sky_exposure", SkyExposure);

                // TEMP/HACK: When using luminance, values are much brighter than radiance.
                // Need a better solution to handle that.
                float luminanceCompensation = m_Atmosphere.Luminance != AtmosphereModel.LuminanceType.none ? 1e-05f : 1.0f;
                m_SkyMaterial.SetFloat("sun_brightness", SunBrightnessMultiplier * luminanceCompensation);

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
                ConfigureMaterial(m_StarMeshMaterial);

                m_StarMeshMaterial.SetFloat("star_intensity_multiplier", m_StarBrightnessMultiplier * SkyExposure);
                m_StarMeshMaterial.SetFloat("star_intensity_power", m_StarBrightnessPower);
                
                float D = m_Atmosphere.PlanetaryRadius + m_Altitude * 1000;
                float R = m_Atmosphere.PlanetaryRadius;
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