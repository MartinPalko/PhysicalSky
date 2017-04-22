using UnityEngine;

namespace PhysicalSky
{
    namespace CelestialCoordinates
    {
        [System.Serializable]
        public struct CartesianCoords
        {
            public float x;
            public float y;
            public float z;

            public CartesianCoords(float x, float y, float z)
            {
                this.x = x;
                this.y = y;
                this.z = z;
            }

            public override string ToString()
            {
                return "(x" + x + ", y" + y + ", z" + z + ")";
            }

            public static implicit operator CartesianCoords(Vector3 v) { return new CartesianCoords(v.x, v.z, v.y); }
            public static implicit operator CartesianCoords(SphericalCoords s) { return Utility.SphericalToCartesian(s); }
            public static implicit operator CartesianCoords(GeographicCoords g) { return Utility.SphericalToCartesian(g); }
            public static implicit operator CartesianCoords(EquitorialCoordinates e) { return Utility.EquitorialToCartesian(e); }

            public Vector3 ToVector3()
            {
                // Cartesian Coords are Z-up, whereas unity is Y-up
                return new Vector3(x, z, y);
            }
        }

        [System.Serializable]
        public struct SphericalCoords
        {
            // Zenith and Azimuth in radians
            public float zenith;
            public float azimuth;
            public float radius;

            public SphericalCoords(float zenith, float azimuth, float radius = 1.0f)
            {
                this.radius = radius;
                this.zenith = zenith;
                this.azimuth = azimuth;
            }

            public override string ToString()
            {
                return "(zenith" + zenith + ", azimuth" + azimuth + ", radius" + radius + ")";
            }

            public static implicit operator SphericalCoords(CartesianCoords c) { return Utility.CartesianToSpherical(c); }
            public static implicit operator SphericalCoords(GeographicCoords g) { return new SphericalCoords((Mathf.Deg2Rad * 90.0f) - g.latitude, g.longitude); }
        }

        [System.Serializable]
        public struct GeographicCoords
        {
            // Latitude and Longitude in radians
            public float latitude;
            public float longitude;

            public GeographicCoords(float latitude, float longitude)
            {
                this.latitude = latitude;
                this.longitude = longitude;
            }

            public override string ToString()
            {
                return "(latitude" + latitude + ", longitude" + longitude + ")";
            }

            public static implicit operator GeographicCoords(SphericalCoords s) { return new GeographicCoords(s.zenith + (Mathf.Deg2Rad * 90.0f), s.azimuth); }
            public static implicit operator GeographicCoords(CartesianCoords c) { return (SphericalCoords)c; }
        }

        [System.Serializable]
        public struct EquitorialCoordinates
        {
            public float ra; // Right Ascension stored in hours, properties provided for minutes, and seconds.            
            public float dec; // Declination stored in degrees, properties provided for arc minutes, and arc seconds.
            public float d; // Distance

            public EquitorialCoordinates(float ra, float dec, float d = 1.0f)
            {
                this.ra = ra;
                this.dec = dec;
                this.d = d;
            }

            public EquitorialCoordinates(int raHours, int raMinutes, float raSeconds, int decDegrees, int decArcMinutes, float decArcSeconds, float distance = 1.0f)
            {
                ra = 0.0f;
                dec = 0.0f;
                d = distance;

                SetRA(raHours, raMinutes, raSeconds);
                SetDec(decDegrees, decArcMinutes, decArcSeconds);
            }

            public override string ToString()
            {
                return "(ra" + ra + ", dec" + dec + ", d" + d + ")";
            }

            public int RAHours { get { return (int)(ra); } }
            public int RAMinutes { get { return (int)((ra * 60.0f) % 60.0f); } }
            public float RASeconds { get { return (float)((ra * 60.0f * 60.0f) % 60.0f); } }

            public int DecDegrees { get { return (int)(dec); } }
            public int DecArcMin { get { return (int)((dec * 60.0f) % 60.0f); } }
            public float DecArcSeconds { get { return (float)((dec * 60.0f * 60.0f) % 60.0f); } }

            public void SetRA(int hours, int minutes, float seconds)
            {
                ra = (hours + minutes / 60.0f + seconds / (60.0f * 60.0f));
            }

            public void SetDec(int degrees, int arcMinutes, float arcSeconds)
            {
                dec = (degrees + arcMinutes / 60.0f + arcSeconds / (60.0f * 60.0f));
            }
        }

        public static class Utility
        {
            public static SphericalCoords CartesianToSpherical(CartesianCoords cartesian)
            {
                float radius = Mathf.Sqrt(cartesian.x * cartesian.x + cartesian.y * cartesian.y + cartesian.z * cartesian.z);
                float azimuth = Mathf.Atan2(cartesian.y, cartesian.x);
                float zenith = Mathf.Acos(cartesian.z / radius);

                return new SphericalCoords(zenith, azimuth, radius);
            }

            public static CartesianCoords SphericalToCartesian(SphericalCoords spherical)
            {
                float sinZenith = Mathf.Sin(spherical.zenith);

                float x = spherical.radius * sinZenith * Mathf.Cos(spherical.azimuth);
                float y = spherical.radius * sinZenith * Mathf.Sin(spherical.azimuth);
                float z = spherical.radius * Mathf.Cos(spherical.zenith);

                return new CartesianCoords(x, y, z);
            }

            public static CartesianCoords EquitorialToCartesian(EquitorialCoordinates equitorial)
            {
                float ra = (equitorial.ra * 15.0f * Mathf.PI) / 180.0f;
                float dec = (equitorial.dec * Mathf.PI) / 180.0f;
                float d = equitorial.d;
                float cosDec = Mathf.Cos(dec);

                float x = d * cosDec * Mathf.Cos(ra);
                float y = d * cosDec * Mathf.Sin(ra);
                float z = d * Mathf.Sin(dec);

                return new CartesianCoords((float)x, (float)y, (float)z);
            }

            public static Matrix4x4 GalacitcToJ2000Transform()
            {
                Vector3 xAxis = ((CartesianCoords)GalacticToJ2000(new CartesianCoords(-1.0f, 0.0f, 0.0f))).ToVector3();
                Vector3 yAxis = ((CartesianCoords)GalacticToJ2000(new CartesianCoords(0.0f, 1.0f, 0.0f))).ToVector3();
                Vector3 zAxis = ((CartesianCoords)GalacticToJ2000(new CartesianCoords(0.0f, 0.0f, 1.0f))).ToVector3();

                //Debug.DrawRay(Vector3.zero, xAxis, Color.red);
                //Debug.DrawRay(Vector3.zero, yAxis, Color.green);
                //Debug.DrawRay(Vector3.zero, zAxis, Color.blue);

                Matrix4x4 returnValue = Matrix4x4.identity;

                returnValue.m00 = xAxis.x;
                returnValue.m01 = xAxis.y;
                returnValue.m02 = xAxis.z;

                returnValue.m10 = yAxis.x;
                returnValue.m11 = yAxis.y;
                returnValue.m12 = yAxis.z;

                returnValue.m20 = zAxis.x;
                returnValue.m21 = zAxis.y;
                returnValue.m22 = zAxis.z;

                return returnValue;
            }

            public static EquitorialCoordinates GalacticToJ2000(GeographicCoords galactic)
            {
                //https://gist.github.com/barentsen/2367839#file-ga2equ-py-L15
                float l = galactic.latitude;
                float b = galactic.longitude;

                // North galactic pole (J2000)
                float pole_ra = 192.859508f * Mathf.Deg2Rad;
                float pole_dec = 27.128336f * Mathf.Deg2Rad;
                float posangle = (122.932f - 90.0f) * Mathf.Deg2Rad;

                float ra = Mathf.Atan2((Mathf.Cos(b) * Mathf.Cos(l - posangle)), (Mathf.Sin(b) * Mathf.Cos(pole_dec) - Mathf.Cos(b) * Mathf.Sin(pole_dec) * Mathf.Sin(l - posangle))) + pole_ra;
                float dec = Mathf.Asin(Mathf.Cos(b) * Mathf.Cos(pole_dec) * Mathf.Sin(l - posangle) + Mathf.Sin(b) * Mathf.Sin(pole_dec));

                return new EquitorialCoordinates((ra * Mathf.Rad2Deg / 360.0f) * 24.0f, dec * Mathf.Rad2Deg);
            }
        }
    }
}