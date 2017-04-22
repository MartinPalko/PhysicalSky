#if PHYSICAL_SKY_TESTS
using UnityEngine;

namespace PhysicalSky
{
    using NUnit.Framework;
    using CelestialCoordinates;

    [TestFixture]
    public class CelestialCoordinatesTests
    {
        static float roundFloat(float source, int decimalPoints)
        {
            float mult = Mathf.Pow(10, decimalPoints);
            return Mathf.Round((source * mult)) / mult;
        }

        [Test]
        public void CartesianToSpherical()
        {
            CartesianCoords sourceCoords = new CartesianCoords(3.0f, 4.0f, 5.0f);
            SphericalCoords castCoords = sourceCoords;
            SphericalCoords answer = new SphericalCoords(0.78539816339745f, 0.92729521800161f, 7.0710678118655f);
            Assert.AreEqual(answer.ToString(), castCoords.ToString()); // Silly, but compare string versions, because they have rounded floats the mitigate floating point precision inaccuracies.
        }

        [Test]
        public void SphericalToCartesian()
        {
            SphericalCoords sourceCoords = new SphericalCoords(1.236f, 0.35f, 5.0f);
            CartesianCoords castCoords = sourceCoords;
            CartesianCoords answer = new CartesianCoords(4.43608079f, 1.619295894f, 1.64288406f);
            Assert.AreEqual(answer.ToString(), castCoords.ToString());
        }

        [Test]
        public void EquitorialToCartesian()
        {
            EquitorialCoordinates sourceCoords = new EquitorialCoordinates(14.4966f, -62.681f, 1.29f);
            CartesianCoords castCoords = sourceCoords;

            // Round, because our test case result has limited precision.
            castCoords.x = roundFloat(castCoords.x, 3);
            castCoords.y = roundFloat(castCoords.y, 3);
            castCoords.z = roundFloat(castCoords.z, 3);

            CartesianCoords answer = new CartesianCoords(-0.470f, -0.360f, -1.146f);
            Assert.AreEqual(answer, castCoords);
        }

        [Test]
        public void GalacticToJ2000()
        {
            // Test case 1
            {
                GeographicCoords testCoord = new GeographicCoords(359.944f * Mathf.Deg2Rad, -0.0461f * Mathf.Deg2Rad);
                EquitorialCoordinates converted = Utility.GalacticToJ2000(testCoord);
                converted.ra = roundFloat(converted.ra, 3);
                converted.dec = roundFloat(converted.dec, 3);
                converted.d = roundFloat(converted.d, 3);

                EquitorialCoordinates answer = new EquitorialCoordinates((266.417f / 360.0f) * 24.0f, -29.008f);
                answer.ra = roundFloat(answer.ra, 3);
                answer.dec = roundFloat(answer.dec, 3);
                answer.d = roundFloat(answer.d, 3);

                Assert.AreEqual(answer, converted);
            }

            // Test case 2
            {
                GeographicCoords testCoord = new GeographicCoords(0.0f * Mathf.Deg2Rad, 0.0f * Mathf.Deg2Rad);
                EquitorialCoordinates converted = Utility.GalacticToJ2000(testCoord);
                converted.ra = roundFloat(converted.ra, 3);
                converted.dec = roundFloat(converted.dec, 3);
                converted.d = roundFloat(converted.d, 3);

                EquitorialCoordinates answer = new EquitorialCoordinates((266.405f / 360.0f) * 24.0f, -28.936f);
                answer.ra = roundFloat(answer.ra, 3);
                answer.dec = roundFloat(answer.dec, 3);
                answer.d = roundFloat(answer.d, 3);

                Assert.AreEqual(answer, converted);
            }
        }

        [Test]
        public void EquitorialCoordinates()
        {
            int testRAHours = 12;
            int testRAMinutes = 23;
            float testRASeconds = 26.20f;
            int testDecDegrees = 84;
            int testDecArcMinutes = 31;
            float testDecArcSeconds = 58.64f;

            // Test converting to/from Hours/minutes/seconds and degrees/arcmin/arcsec.
            EquitorialCoordinates testCoordinate = new EquitorialCoordinates(testRAHours, testRAMinutes, testRASeconds, testDecDegrees, testDecArcMinutes, testDecArcSeconds);
            testRASeconds = roundFloat(testRASeconds, 0);
            testDecArcSeconds = roundFloat(testDecArcSeconds, 0);
            float answerRASeconds = roundFloat(testCoordinate.RASeconds, 0);
            float answerDecArcSeconds = roundFloat(testCoordinate.DecArcSeconds, 0);
            Assert.AreEqual(testRAHours, testCoordinate.RAHours, "Ra Hours");
            Assert.AreEqual(testRAMinutes, testCoordinate.RAMinutes, "Ra Minutes");
            Assert.AreEqual(testRASeconds, answerRASeconds, "Ra Seconds");
            Assert.AreEqual(testDecDegrees, testCoordinate.DecDegrees, "Dec Degrees");
            Assert.AreEqual(testDecArcMinutes, testCoordinate.DecArcMin, "Dec Arc Minutes");
            Assert.AreEqual(testDecArcSeconds, answerDecArcSeconds, "Dec Arc Seconds");
        }

    }
}
#endif