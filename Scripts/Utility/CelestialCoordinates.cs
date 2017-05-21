using UnityEngine;
using Math = System.Math;
using DateTime = System.DateTime;

namespace PhysicalSky
{
    namespace CelestialCoordinates
    {
        [System.Serializable]
        public struct CartesianCoords
        {
            public double x;
            public double y;
            public double z;

            public CartesianCoords(double x, double y, double z)
            {
                this.x = x;
                this.y = y;
                this.z = z;
            }

            public CartesianCoords(Vector3 v)
            {
                x = v.x;
                y = v.y;
                z = v.z;
            }

            public override string ToString()
            {
                return "(x" + x + ", y" + y + ", z" + z + ")";
            }

            public static implicit operator CartesianCoords(Vector3 v) { return new CartesianCoords(v.x, v.y, v.z); }
            public static implicit operator CartesianCoords(GeographicCoords g) { return Utility.GeographicToCartesian(g); }
            public static implicit operator CartesianCoords(HorizontalCoords h) { return Utility.HorizontalToCartesian(h); }
            public static implicit operator CartesianCoords(EquitorialCoords e) { return Utility.EquitorialToCartesian(e); }

            public Vector3 ToVector3()
            {
                // Cartesian Coords are Z-up, whereas unity is Y-up
                return new Vector3((float)x, (float)y, (float)z);
            }
        }

        [System.Serializable]
        public struct GeographicCoords
        {
            // Latitude and Longitude in radians
            public double latitude;
            public double longitude;

            public GeographicCoords(double latitude, double longitude)
            {
                this.latitude = latitude;
                this.longitude = longitude;
            }

            public override string ToString()
            {
                return "(latitude" + latitude + ", longitude" + longitude + ")";
            }
            
            public static implicit operator GeographicCoords(CartesianCoords c) { return Utility.CartesianToGeographic(c); }
        }

        [System.Serializable]
        public struct HorizontalCoords
        {
            // Altitude and Azimuth in radians
            public double altitude;
            public double azimuth;

            public HorizontalCoords(double altitude, double azimuth)
            {
                this.altitude = altitude;
                this.azimuth = azimuth;
            }

            public override string ToString()
            {
                return "(altitude " + altitude + ", azimuth " + azimuth + ")";
            }

            public static implicit operator HorizontalCoords(CartesianCoords c) { return Utility.CartesianToHorizontal(c); }
        }

        [System.Serializable]
        public struct EquitorialCoords
        {
            public double ra; // Right Ascension stored in hours, properties provided for minutes, and seconds.            
            public double dec; // Declination stored in degrees, properties provided for arc minutes, and arc seconds.
            public double d; // Distance

            // ra and dec in Radians
            public double raRad
            {
                get { return ra * (Math.PI / 12.0); }
                set { ra = value / (Mathf.PI / 12.0); }
            }

            public double decRad
            {
                get { return dec * Utility.Deg2Rad; }
                set { dec = value * Utility.Rad2Deg; }
            }

            public EquitorialCoords(double ra, double dec, double d = 1.0)
            {
                this.ra = ra;
                this.dec = dec;
                this.d = d;
            }

            public EquitorialCoords(int raHours, int raMinutes, double raSeconds, int decDegrees, int decArcMinutes, double decArcSeconds, double distance = 1.0f)
            {
                ra = 0.0;
                dec = 0.0;
                d = distance;

                SetRA(raHours, raMinutes, raSeconds);
                SetDec(decDegrees, decArcMinutes, decArcSeconds);
            }

            public override string ToString()
            {
                return "(ra" + ra + ", dec" + dec + ", d" + d + ")";
            }

            public int RAHours { get { return (int)(ra); } }
            public int RAMinutes { get { return (int)((ra * 60.0) % 60.0); } }
            public double RASeconds { get { return (double)((ra * 60.0 * 60.0) % 60.0); } }

            public int DecDegrees { get { return (int)(dec); } }
            public int DecArcMin { get { return (int)((dec * 60.0) % 60.0); } }
            public double DecArcSeconds { get { return (double)((dec * 60.0 * 60.0) % 60.0); } }

            public void SetRA(int hours, int minutes, double seconds)
            {
                ra = (hours + minutes / 60.0 + seconds / (60.0 * 60.0));
            }

            public void SetDec(int degrees, int arcMinutes, double arcSeconds)
            {
                dec = (degrees + arcMinutes / 60.0 + arcSeconds / (60.0 * 60.0));
            }

            public static implicit operator EquitorialCoords(CartesianCoords c) { return Utility.CartesianToEquitorial(c); }
        }

        public static class Utility
        {
            public const double twoPI = Math.PI * 2.0;
            public const double twoPIf = Mathf.PI * 2.0;
            public const double Deg2Rad = Math.PI / 180.0;
            public const double Rad2Deg = 180.0 / Math.PI;

            public static GeographicCoords CartesianToGeographic(CartesianCoords cartesian)
            {
                return new GeographicCoords(Math.Asin(cartesian.z), Math.Atan2(cartesian.y, cartesian.x));
            }

            public static CartesianCoords GeographicToCartesian(GeographicCoords geographic)
            {
                double cosAltitude = Math.Cos(geographic.latitude);

                double x = cosAltitude * Math.Cos(geographic.longitude);
                double y = cosAltitude * Math.Sin(geographic.longitude);
                double z = Math.Sin(geographic.latitude);

                return new CartesianCoords(x, y, z);
            }

            public static HorizontalCoords CartesianToHorizontal(CartesianCoords cartesian)
            {
                return new HorizontalCoords(Math.Asin(cartesian.z) , Math.Atan2(cartesian.y, cartesian.x));
            }

            public static CartesianCoords HorizontalToCartesian(HorizontalCoords horizontal)
            {
                double cosAltitude = Math.Cos(horizontal.altitude);

                double x = cosAltitude * Math.Cos(horizontal.azimuth);
                double y = cosAltitude * Math.Sin(horizontal.azimuth);
                double z = Math.Sin(horizontal.altitude);

                return new CartesianCoords(x, y, z);
            }

            public static CartesianCoords EquitorialToCartesian(EquitorialCoords equitorial)
            {
                double ra = equitorial.raRad;
                double dec = equitorial.decRad;
                double d = equitorial.d;
                double cosDec = Math.Cos(dec);

                double x = d * cosDec * Math.Cos(ra);
                double y = d * cosDec * Math.Sin(ra);
                double z = d * Math.Sin(dec);

                return new CartesianCoords(x, y, z);
            }

            public static EquitorialCoords CartesianToEquitorial(CartesianCoords cartesian)
            {
                throw new System.NotImplementedException();
            }

            public static Matrix4x4 GalacitcToJ2000Transform()
            {
                Vector3 xAxis = ((CartesianCoords)GalacticToJ2000(new CartesianCoords(1.0, 0.0, 0.0))).ToVector3();
                Vector3 yAxis = ((CartesianCoords)GalacticToJ2000(new CartesianCoords(0.0, -1.0, 0.0))).ToVector3();
                Vector3 zAxis = ((CartesianCoords)GalacticToJ2000(new CartesianCoords(0.0, 0.0, 1.0))).ToVector3();

                Debug.DrawRay(Vector3.zero, xAxis, Color.red);
                Debug.DrawRay(Vector3.zero, yAxis, Color.green);
                Debug.DrawRay(Vector3.zero, zAxis, Color.blue);

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

            public static EquitorialCoords GalacticToJ2000(GeographicCoords galactic)
            {
                // https://gist.github.com/barentsen/2367839#file-ga2equ-py-L15
                double l = galactic.latitude;
                double b = galactic.longitude;

                // North galactic pole (J2000)
                double pole_ra = 192.859508 * Mathf.Deg2Rad;
                double pole_dec = 27.128336 * Mathf.Deg2Rad;
                double posangle = (122.932 - 90.0) * Mathf.Deg2Rad;

                double ra = Math.Atan2((Math.Cos(b) * Math.Cos(l - posangle)), (Math.Sin(b) * Math.Cos(pole_dec) - Math.Cos(b) * Math.Sin(pole_dec) * Math.Sin(l - posangle))) + pole_ra;
                double dec = Math.Asin(Math.Cos(b) * Math.Cos(pole_dec) * Math.Sin(l - posangle) + Math.Sin(b) * Math.Sin(pole_dec));

                return new EquitorialCoords((ra * Mathf.Rad2Deg / 360.0) * 24.0, dec * Mathf.Rad2Deg);
            }

            public static HorizontalCoords EquitorialToHorizontal(EquitorialCoords equitorial, GeographicCoords observer, DateTime time)
            {
                // http://www.convertalot.com/celestial_horizon_co-ordinates_calculator.html

                double ra = equitorial.raRad;
                double dec = equitorial.decRad;

                double lat = observer.latitude;
                double lon = observer.longitude;

                var ha = MeanSiderealTime(time, lon) - ra;
                if (ha < 0)
                    ha += twoPI;

                var sin_alt = Math.Sin(dec) * Math.Sin(lat) + Math.Cos(dec) * Math.Cos(lat) * Math.Cos(ha);
                var alt = Math.Asin(sin_alt);

                // compute azimuth in radians
                // divide by zero error at poles or if alt = 90 deg
                var cos_az = (Math.Sin(dec) - Math.Sin(alt) * Math.Sin(lat)) / (Math.Cos(alt) * Math.Cos(lat));
                var az = Math.Acos(cos_az);

                var hrz_altitude = alt;
                var hrz_azimuth = az;

                // choose hemisphere
                if (Math.Sin(ha) > 0)
                    hrz_azimuth = (twoPI) - hrz_azimuth;

                return new HorizontalCoords(hrz_altitude, hrz_azimuth);
            }

            public static double MeanSiderealTime(DateTime time, double longitude = 0.0)
            {
                double lon = longitude * (180.0 / Math.PI); // Rad to deg

                int year = time.Year;
                int month = time.Month;
                int day = time.Day;
                int hour = time.Hour;
                int minute = time.Minute;
                int second = time.Second;

                if ((month == 1) || (month == 2))
                {
                    year = year - 1;
                    month = month + 12;
                }

                var a = Math.Floor(year / 100.0);
                var b = 2 - a + Math.Floor(a / 4.0);
                var c = Math.Floor(365.25 * year);
                var d = Math.Floor(30.6001 * (month + 1));

                // Days since J2000.0
                var jd = b + c + d - 730550.5 + day + (hour + minute / 60.0 + second / 3600.0) / 24.0;

                // Julian centuries since J2000.0
                var jt = jd / 36525.0;

                // The mean sidereal time in degrees
                var mst = 280.46061837 + 360.98564736629 * jd + 0.000387933 * jt * jt - jt * jt * jt / 38710000 + lon;

                // In degrees modulo 360.0
                if (mst > 0.0)
                    while (mst > 360.0) mst = mst - 360.0;
                else
                    while (mst < 0.0) mst = mst + 360.0;

                return Math.PI * mst / 180.0; // Deg to rad

            }
        }
    }
}