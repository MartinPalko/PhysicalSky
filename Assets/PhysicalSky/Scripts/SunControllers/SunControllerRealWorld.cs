using System;
using UnityEngine;
using PhysicalSky.Utilities;

namespace PhysicalSky
{
    public class SunControllerRealWorld : SunController
    {
        [SerializeField]
        protected float latitude = 45.5017f;
        [SerializeField]
        protected float longitude = -73.5673f;
        [SerializeField]
        protected DateTime time = DateTime.Now;
        [SerializeField]
        protected TimeSpan timezone = TimeSpan.FromHours(5); // EST
        [SerializeField]
        protected float minutesPerHour = 60.0f;
        [SerializeField]
        protected float sunsetSpeedMultiplier = 0.2f;
        [SerializeField]
        protected float nightSpeedMultiplier = 10.0f;

        private void Update()
        {
            float sunsetY = -0.025f; // TODO: Find actual horizon angle based on altitude & planetary radius.
            float sunSpeed = 1.0f;

            sunSpeed *= sky.SunDirection.y < sunsetY ? nightSpeedMultiplier : 1.0f;
            float distanceToSunset = Mathf.Clamp01(Mathf.Abs(sky.SunDirection.y - sunsetY));
            sunSpeed *= Mathf.Lerp(sunsetSpeedMultiplier, 1.0f, distanceToSunset);

            time = time.AddHours((Time.deltaTime / 3600.0f) * sunSpeed * (60.0f / minutesPerHour));

            double sunAltitude;
            double sunAzimuth;            
            SunPosition.CalculateSunPosition(time, latitude, longitude, out sunAltitude, out sunAzimuth);

            sky.SunDirection = CartesianCoords.SphericalToCartesian((float)sunAzimuth, (float)sunAltitude);
        }

    }
}