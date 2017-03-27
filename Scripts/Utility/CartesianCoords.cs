using UnityEngine;

namespace PhysicalSky.Utilities
{
    public struct CartesianCoords
    {
        public static Vector3 SphericalToCartesian(float polar, float elevation)
        {
            float a = Mathf.Cos(elevation);
            return new Vector3(a * Mathf.Cos(polar), Mathf.Sin(elevation), a * Mathf.Sin(polar));
        }

        public static void CartesianToSpherical(Vector3 cartCoords, out float outRadius, out float outPolar, out float outElevation)
        {
            if (cartCoords.x == 0)
                cartCoords.x = Mathf.Epsilon;
            outRadius = Mathf.Sqrt((cartCoords.x * cartCoords.x)
                            + (cartCoords.y * cartCoords.y)
                            + (cartCoords.z * cartCoords.z));
            outPolar = Mathf.Atan(cartCoords.z / cartCoords.x);
            if (cartCoords.x < 0)
                outPolar += Mathf.PI;
            outElevation = Mathf.Asin(cartCoords.y / outRadius);
        }
    }
}