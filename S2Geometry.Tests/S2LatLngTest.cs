using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Google.Common.Geometry;
using NUnit.Framework;

namespace S2Geometry.Tests
{
    public class S2LatLngTest : GeometryTestCase
    {
        [Test]
        public void testBasic()
        {
            var llRad = S2LatLng.FromRadians(S2.PiOver4, S2.PiOver2);
            assertTrue(llRad.Lat.Radians == S2.PiOver4);
            assertTrue(llRad.Lng.Radians == S2.PiOver2);
            assertTrue(llRad.IsValid);
            var llDeg = S2LatLng.FromDegrees(45, 90);
            assertEquals(llDeg, llRad);
            assertTrue(llDeg.IsValid);
            assertTrue(!S2LatLng.FromDegrees(-91, 0).IsValid);
            assertTrue(!S2LatLng.FromDegrees(0, 181).IsValid);

            var bad = S2LatLng.FromDegrees(120, 200);
            assertTrue(!bad.IsValid);
            var better = bad.Normalized;
            assertTrue(better.IsValid);
            assertEquals(better.Lat, S1Angle.FromDegrees(90));
            assertDoubleNear(better.Lng.Radians, S1Angle.FromDegrees(-160).Radians);

            bad = S2LatLng.FromDegrees(-100, -360);
            assertTrue(!bad.IsValid);
            better = bad.Normalized;
            assertTrue(better.IsValid);
            assertEquals(better.Lat, S1Angle.FromDegrees(-90));
            assertDoubleNear(better.Lng.Radians, 0);

            assertTrue((S2LatLng.FromDegrees(10, 20) + S2LatLng.FromDegrees(20, 30)).ApproxEquals(
                S2LatLng.FromDegrees(30, 50)));
            assertTrue((S2LatLng.FromDegrees(10, 20) - S2LatLng.FromDegrees(20, 30)).ApproxEquals(
                S2LatLng.FromDegrees(-10, -10)));
            assertTrue((S2LatLng.FromDegrees(10, 20)*0.5).ApproxEquals(S2LatLng.FromDegrees(5, 10)));
        }

        [Test]
        public void testConversion()
        {
            // Test special cases: poles, "date line"
            assertDoubleNear(
                new S2LatLng(S2LatLng.FromDegrees(90.0, 65.0).ToPoint()).Lat.Degrees, 90.0);
            assertEquals(
                new S2LatLng(S2LatLng.FromRadians(-S2.PiOver2, 1).ToPoint()).Lat.Radians, -S2.PiOver2);
            assertDoubleNear(
                Math.Abs(new S2LatLng(S2LatLng.FromDegrees(12.2, 180.0).ToPoint()).Lng.Degrees), 180.0);
            assertEquals(
                Math.Abs(new S2LatLng(S2LatLng.FromRadians(0.1, -S2.Pi).ToPoint()).Lng.Radians),
                S2.Pi);

            // Test a bunch of random points.
            for (var i = 0; i < 100000; ++i)
            {
                var p = randomPoint();
                assertTrue(S2.ApproxEquals(p, new S2LatLng(p).ToPoint()));
            }

            // Test generation from E5
            var test = S2LatLng.FromE5(123456, 98765);
            assertDoubleNear(test.Lat.Degrees, 1.23456);
            assertDoubleNear(test.Lng.Degrees, 0.98765);
        }

        [Test]
        public void testDistance()
        {
            assertEquals(
                S2LatLng.FromDegrees(90, 0).GetDistance(S2LatLng.FromDegrees(90, 0)).Radians, 0.0);
            assertDoubleNear(
                S2LatLng.FromDegrees(-37, 25).GetDistance(S2LatLng.FromDegrees(-66, -155)).Degrees, 77,
                1e-13);
            assertDoubleNear(
                S2LatLng.FromDegrees(0, 165).GetDistance(S2LatLng.FromDegrees(0, -80)).Degrees, 115,
                1e-13);
            assertDoubleNear(
                S2LatLng.FromDegrees(47, -127).GetDistance(S2LatLng.FromDegrees(-47, 53)).Degrees, 180,
                2e-6);
        }
    }
}