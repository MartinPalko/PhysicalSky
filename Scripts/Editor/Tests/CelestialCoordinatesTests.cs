#if PHYSICAL_SKY_TESTS
using UnityEngine;
using Math = System.Math;

namespace PhysicalSky
{
    using NUnit.Framework;
    using CelestialCoordinates;

    [TestFixture]
    public class CelestialCoordinatesTests
    {
        [Test]
        public void DegRadConversions()
        {
            Assert.AreEqual(Math.Round(Utility.Deg2Rad, 4), Math.Round(Mathf.Deg2Rad, 4));
            Assert.AreEqual(Math.Round(Utility.Rad2Deg, 4), Math.Round(Mathf.Rad2Deg, 4));
        }

        [Test]
        public void EquitorialToCartesian()
        {
            EquitorialCoords sourceCoords = new EquitorialCoords(14.4966, -62.681, 1.29);
            CartesianCoords castCoords = sourceCoords;

            // Round, because our test case result has limited precision.
            castCoords.x = Math.Round(castCoords.x, 3);
            castCoords.y = Math.Round(castCoords.y, 3);
            castCoords.z = Math.Round(castCoords.z, 3);

            CartesianCoords answer = new CartesianCoords(-0.470, -0.360, -1.146);
            Assert.AreEqual(answer, castCoords);
        }

        [Test]
        public void GalacticToJ2000()
        {
            // Test case 1
            {
                GeographicCoords testCoord = new GeographicCoords(359.944 * Utility.Deg2Rad, -0.0461 * Utility.Deg2Rad);
                EquitorialCoords converted = Utility.GalacticToJ2000(testCoord);
                converted.ra = Math.Round(converted.ra, 2);
                converted.dec = Math.Round(converted.dec, 2);
                converted.d = Math.Round(converted.d, 2);

                EquitorialCoords answer = new EquitorialCoords((266.417 / 360.0) * 24.0, -29.008);
                answer.ra = Math.Round(answer.ra, 2);
                answer.dec = Math.Round(answer.dec, 2);
                answer.d = Math.Round(answer.d, 2);

                Assert.AreEqual(answer, converted);
            }

            // Test case 2
            {
                GeographicCoords testCoord = new GeographicCoords(0.0 * Utility.Deg2Rad, 0.0 * Utility.Deg2Rad);
                EquitorialCoords converted = Utility.GalacticToJ2000(testCoord);
                converted.ra = Math.Round(converted.ra, 2);
                converted.dec = Math.Round(converted.dec, 2);
                converted.d = Math.Round(converted.d, 2);

                EquitorialCoords answer = new EquitorialCoords((266.405 / 360.0) * 24.0, -28.936);
                answer.ra = Math.Round(answer.ra, 2);
                answer.dec = Math.Round(answer.dec, 2);
                answer.d = Math.Round(answer.d, 2);

                Assert.AreEqual(answer, converted);
            }
        }

        [Test]
        public void EquitorialCoordinates()
        {
            int testRAHours = 12;
            int testRAMinutes = 23;
            double testRASeconds = 26.20;
            int testDecDegrees = 84;
            int testDecArcMinutes = 31;
            double testDecArcSeconds = 58.64;

            // Test converting to/from Hours/minutes/seconds and degrees/arcmin/arcsec.
            EquitorialCoords testCoordinate = new EquitorialCoords(testRAHours, testRAMinutes, testRASeconds, testDecDegrees, testDecArcMinutes, testDecArcSeconds);
            testRASeconds = Math.Round(testRASeconds, 1);
            testDecArcSeconds = Math.Round(testDecArcSeconds, 1);
            double answerRASeconds = Math.Round(testCoordinate.RASeconds, 1);
            double answerDecArcSeconds = Math.Round(testCoordinate.DecArcSeconds, 1);
            Assert.AreEqual(testRAHours, testCoordinate.RAHours, "Ra Hours");
            Assert.AreEqual(testRAMinutes, testCoordinate.RAMinutes, "Ra Minutes");
            Assert.AreEqual(testRASeconds, answerRASeconds, "Ra Seconds");
            Assert.AreEqual(testDecDegrees, testCoordinate.DecDegrees, "Dec Degrees");
            Assert.AreEqual(testDecArcMinutes, testCoordinate.DecArcMin, "Dec Arc Minutes");
            Assert.AreEqual(testDecArcSeconds, answerDecArcSeconds, "Dec Arc Seconds");
        }

        [Test]
        public void EquitorialToHorizontal()
        {
            // Tested at: http://www.convertalot.com/celestial_horizon_co-ordinates_calculator.html
            EquitorialCoords testPos = new EquitorialCoords(18, 10, 0, 4, 46, 0);
            GeographicCoords observerPos = new GeographicCoords(10 * Utility.Deg2Rad, 45 * Utility.Deg2Rad);
            System.DateTime testTime = new System.DateTime(2017, 08, 12, 08, 52, 16);

            HorizontalCoords expectedResult = new HorizontalCoords(-41.00848542273538 * Utility.Deg2Rad, 73.23004573869906 * Utility.Deg2Rad);
            HorizontalCoords actualResult = Utility.EquitorialToHorizontal(testPos, observerPos, testTime);

            Assert.AreEqual(expectedResult, actualResult);
        }

    }
}
#endif