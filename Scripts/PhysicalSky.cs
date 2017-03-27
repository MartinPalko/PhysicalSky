using UnityEngine;
using System;

namespace PhysicalSky
{
    [ExecuteInEditMode]
    [RequireComponent(typeof(Light))]
    public class PhysicalSky : MonoBehaviour, IPhysicalSky
    {
        public Vector3 SunDirection { get { return -transform.forward; } set { transform.rotation = Quaternion.LookRotation(-value); } }

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
        public float Altitude { get { return altitude; } set { altitude = Mathf.Max(value, 0.0f); } }

        [SerializeField]
        private AtmosphereModel atmosphere;
        public IAtmosphereModel Atmosphere { get { return atmosphere; } set { atmosphere = value as AtmosphereModel; } }

        [SerializeField]
        private Shader skyShader;
        [SerializeField]
        private Shader sunRadianceShader;

        private Material skyMaterial;
        private Material sunRadianceMaterial;

        private const int SUN_RADIANCE_TEXTURE_SIZE = 8;
        private Light sunLight;
        private Texture2D sunRadianceTexture;

        private void Start()
        {
            skyMaterial = new Material(skyShader);
            sunRadianceMaterial = new Material(sunRadianceShader);
            sunLight = GetComponent<Light>();
            sunLight.type = LightType.Directional;
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
                RenderSettings.skybox = skyMaterial;
                RenderSettings.sun = sunLight;
            }
        }
    }
}