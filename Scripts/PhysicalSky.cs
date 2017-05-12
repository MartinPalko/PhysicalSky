using UnityEngine;

namespace PhysicalSky
{
    [ExecuteInEditMode]
    [RequireComponent(typeof(Light))]
    public class PhysicalSky : MonoBehaviour, IPhysicalSky
    {
        public Vector3 SunDirection { get { return -transform.forward; } set { transform.rotation = Quaternion.LookRotation(-value); } }
        public Quaternion StarRotation { get { return starMeshObject.transform.rotation; } set { starMeshObject.transform.rotation = value; } }

        [SerializeField]
        float sunBrightnessMultiplier = 1.0f;
        public float SunBrightnessMultiplier { get { return sunBrightnessMultiplier; } set { sunBrightnessMultiplier = Mathf.Max(value, 0.0f); } }

        [SerializeField]
        float sunSaturationMultiplier = 1.0f;
        public float SunSaturationMultiplier { get { return sunSaturationMultiplier; } set { sunSaturationMultiplier = Mathf.Max(value, 0.0f); } }

        [SerializeField]
        float skyExposure = 20.0f;
        public float SkyExposure { get { return skyExposure; } set { skyExposure = Mathf.Max(value, 0.0f); } }

        [SerializeField]
        private float altitude = 0.1f;
        public float Altitude { get { return altitude; } set { altitude = Mathf.Max(value, 0.0f); DynamicGI.UpdateEnvironment(); } }

        [SerializeField]
        private AtmosphereModel atmosphere;
        public IAtmosphereModel Atmosphere { get { return atmosphere; } set { atmosphere = value as AtmosphereModel; } }

        [SerializeField]
        float starBrightnessMultiplier = 1.0f;
        public float StarBrightnessMultiplier { get { return starBrightnessMultiplier; } set { starBrightnessMultiplier = Mathf.Max(value, 0.0f); } }

        [SerializeField]
        float starBrightnessPower = 0.454545f;
        public float StarBrightnessPower { get { return starBrightnessPower; } set { starBrightnessPower = Mathf.Max(value, 0.0f); } }

        [SerializeField]
        private StarMap starMap;
        public StarMap StarMap { get { return starMap; } set { starMap = value; } }

        [SerializeField]
        private Shader skyShader = null;
        [SerializeField]
        private Shader sunRadianceShader = null;
        [SerializeField]
        private Shader starMeshShader = null;

        private Material skyMaterial = null;
        private Material sunRadianceMaterial = null;

        private Material starMeshMaterial = null;
        private GameObject starMeshObject = null;
        private MeshFilter starMeshFilter = null;

        private const int SUN_RADIANCE_TEXTURE_SIZE = 8;
        private Light sunLight = null;
        private Texture2D sunRadianceTexture = null;

        private void Start()
        {
            skyMaterial = new Material(skyShader);
            sunRadianceMaterial = new Material(sunRadianceShader);
            starMeshMaterial = new Material(starMeshShader);

            starMeshObject = new GameObject("StarMesh");
#if PHYSICAL_SKY_DEBUG
            starMeshObject.hideFlags = HideFlags.DontSave;
#else
            starMeshObject.hideFlags = HideFlags.HideAndDontSave;
#endif
            starMeshFilter = starMeshObject.AddComponent<MeshFilter>();
            starMeshFilter.mesh = StarMap.CreateStarMesh();     
            MeshRenderer starMeshRenderer = starMeshObject.AddComponent<MeshRenderer>();
            starMeshRenderer.material = starMeshMaterial;
            
            sunLight = GetComponent<Light>();
            sunLight.type = LightType.Directional;
        }

        private void OnDestroy()
        {
#if UNITY_EDITOR
            if (!Application.isPlaying)
            {
                DestroyImmediate(starMeshObject);
            }
            else
#endif
            {
                Destroy(starMeshObject);
            }
        }

        private void ConfigureMaterial(Material m)
        {
            atmosphere.SetShaderUniforms(m);
            m.SetTexture("transmittance_texture", atmosphere.TransmittanceLUT);
            m.SetTexture("scattering_texture", atmosphere.ScatteringLUT);
            m.SetTexture("irradiance_texture", atmosphere.IrradianceLUT);

            m.SetVector("camera", new Vector3(0, (atmosphere.PlanetaryRadius / 1000) + altitude, 0));
        }

        private void Update()
        {
#if UNITY_EDITOR
            if (!Application.isPlaying)
            {
                sunLight = GetComponent<Light>();                
            }
#endif

            if (!atmosphere)
            {
                Debug.LogWarning("No atmosphere");
                return;
            }

            if (atmosphere.NeedsRecompute)
                atmosphere.Compute();

            if (!sunRadianceTexture)
            {
                sunRadianceTexture = new Texture2D(SUN_RADIANCE_TEXTURE_SIZE, SUN_RADIANCE_TEXTURE_SIZE, TextureFormat.RGBAHalf, false, false);
            }

            // Compute sun radiance, and apply to light
            if (sunRadianceMaterial)
            {
                ConfigureMaterial(sunRadianceMaterial);
                sunRadianceMaterial.SetMatrix("sun_rotation_matrix", transform.localToWorldMatrix);
                sunRadianceMaterial.SetVector("sun_direction", -transform.forward);

                RenderTexture temporary = RenderTexture.GetTemporary(SUN_RADIANCE_TEXTURE_SIZE, SUN_RADIANCE_TEXTURE_SIZE, 0, RenderTextureFormat.ARGBHalf);
                RenderTexture.active = temporary;
                Graphics.Blit(null, temporary, sunRadianceMaterial);
                sunRadianceTexture.ReadPixels(new Rect(0, 0, SUN_RADIANCE_TEXTURE_SIZE, SUN_RADIANCE_TEXTURE_SIZE), 0, 0);
                RenderTexture.active = null;
                RenderTexture.ReleaseTemporary(temporary);

                Color[] colors = sunRadianceTexture.GetPixels();

                Color radianceSum = Color.black;
                int nonZero = 0;

                for (int i = 0; i < colors.Length; i++)
                {
                    radianceSum += colors[i];
                    nonZero += colors[i] == Color.black ? 0 : 1;
                }
                Color averageRadiance = radianceSum / nonZero;
                averageRadiance /= 3000; // Arbitrary value scales down to unity's lighting.
                averageRadiance *= SunBrightnessMultiplier;

                Utilities.HSBColor hsbLightColor = new Utilities.HSBColor(averageRadiance);
                sunLight.intensity = hsbLightColor.b;
                hsbLightColor.b = 1.0f;
                hsbLightColor.s = Mathf.Clamp01(hsbLightColor.s * SunSaturationMultiplier);
                sunLight.color = hsbLightColor.ToColor().gamma;
            }

            // Assign sky material as skybox.
            if (skyMaterial)
            {
                ConfigureMaterial(skyMaterial);
                skyMaterial.SetFloat("sky_exposure", SkyExposure);

                if (starMap)
                {
                    skyMaterial.SetTexture("star_cubemap", starMap.BackgroundCube);
                    skyMaterial.SetFloat("star_brightness", starMap.BackgroundCubeBrightness * starBrightnessMultiplier);
                    skyMaterial.SetMatrix("star_rotation", Matrix4x4.TRS(Vector3.zero, Quaternion.Euler(starMap.BackgroundRotation), Vector3.one));
                }
                else
                {
                    skyMaterial.SetTexture("star_cubemap", null);
                    skyMaterial.SetFloat("star_brightness", 0.0f);
                    skyMaterial.SetMatrix("star_rotation", Matrix4x4.identity);
                }

                RenderSettings.skybox = skyMaterial;
                RenderSettings.sun = sunLight;
            }

            if (starMeshMaterial)
            {
                ConfigureMaterial(starMeshMaterial);

                starMeshMaterial.SetFloat("star_intensity_multiplier", starBrightnessMultiplier * SkyExposure);
                starMeshMaterial.SetFloat("star_intensity_power", starBrightnessPower);
                
                float D = atmosphere.PlanetaryRadius + altitude * 1000;
                float R = atmosphere.PlanetaryRadius;
                float planetAngularRad = 2 * Mathf.Acos(Mathf.Sqrt(Mathf.Pow(D, 2) - Mathf.Pow(R, 2)) / D);
                starMeshMaterial.SetFloat("planet_size", planetAngularRad);              
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