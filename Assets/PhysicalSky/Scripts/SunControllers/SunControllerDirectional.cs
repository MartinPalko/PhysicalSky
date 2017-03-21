using UnityEngine;
using PhysicalSky.Utilities;

namespace PhysicalSky
{
    public class SunControllerDirectional : SunController
    {
        [SerializeField][Range(-360.0f, 360.0f)]
        float sunAzimuth = 0.0f;
        [SerializeField][Range(-90.0f, 90.0f)]
        float sunAltitude = 0.0f;       

        protected void Update()
        {
            sky.SunDirection = CartesianCoords.SphericalToCartesian(sunAzimuth * (Mathf.PI / 180.0f), sunAltitude * (Mathf.PI / 180.0f));
        }
    }
}