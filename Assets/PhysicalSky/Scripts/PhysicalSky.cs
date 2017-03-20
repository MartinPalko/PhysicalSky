using UnityEngine;
using System;

namespace PhysicalSky
{
    public class PhysicalSky : MonoBehaviour
    {
        public AtmosphereModel atmosphereModel;

        public Light sunLight;

        public Shader skyShader;
        private Material skyMaterial;

        // Serialized Fields.
        [SerializeField]
        private float altitude = 0.1f;
        public float Altitude
        {
            get { return altitude; }
            set { altitude = value; }
        }

        private void Start()
        {
            skyMaterial = new Material(skyShader);

            atmosphereModel.ComputeLookupTextures();
        }

        DateTime currentTime = DateTime.Now - TimeSpan.FromHours(5); // Time in EST
        double latitude = 45.5017 / (Math.PI * 2.0);
        double longitude = 73.5673 / (Math.PI * 2.0);
        Vector3 sunDirection;

        private void Update()
        {
            if (atmosphereModel.NeedsRecompute)
                atmosphereModel.ComputeLookupTextures();

            //// Slows down at sunrise/sunset, and speeds up overnight.
            //float sunSpeed = Mathf.Abs(sunDirection.y + 0.1f) * (sunDirection.y < -0.025f ? 10.0f : 1.0f) + 0.2f;
            //sunSpeed = 0;
            //currentTime = currentTime.AddHours(Time.deltaTime * sunSpeed);
            //double sunAltitude;
            //double sunAzimuth;
            //SunPosition.CalculateSunPosition(currentTime, latitude, longitude, out sunAltitude, out sunAzimuth);
            //CartesianCoords.SphericalToCartesian(1, (float)sunAzimuth, (float)sunAltitude, out sunDirection);
            //sunLight.transform.rotation = Quaternion.LookRotation(-sunDirection);


            Vector3 sunDirection = -sunLight.transform.forward;

            atmosphereModel.SetAtmosphereUniforms(skyMaterial);
            skyMaterial.SetTexture("transmittance_texture", atmosphereModel.TransmittanceLUT);
            skyMaterial.SetTexture("scattering_texture", atmosphereModel.ScatteringLUT);
            skyMaterial.SetTexture("irradiance_texture", atmosphereModel.IrradianceLUT);

            Color sunRadiance = new Color(1.0f, 1.0f, 1.0f);
            Vector3 sunSize = new Vector3(1.0f, 1.0f, 1.0f);

            skyMaterial.SetVector("camera", new Vector3(0, (float)(atmosphereModel.PlanetaryRadius / 1000) + altitude, 0));
            skyMaterial.SetVector("sun_direction", sunDirection.normalized);

            RenderSettings.skybox = skyMaterial;
        }
    }
}