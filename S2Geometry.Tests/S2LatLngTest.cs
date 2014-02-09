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
            var llRad = S2LatLng.fromRadians(S2.PiOver4, S2.PiOver2);
            assertTrue(llRad.lat().Radians == S2.PiOver4);
            assertTrue(llRad.lng().Radians == S2.PiOver2);
            assertTrue(llRad.isValid());
            var llDeg = S2LatLng.fromDegrees(45, 90);
            assertEquals(llDeg, llRad);
            assertTrue(llDeg.isValid());
            assertTrue(!S2LatLng.fromDegrees(-91, 0).isValid());
            assertTrue(!S2LatLng.fromDegrees(0, 181).isValid());

            var bad = S2LatLng.fromDegrees(120, 200);
            assertTrue(!bad.isValid());
            var better = bad.normalized();
            assertTrue(better.isValid());
            assertEquals(better.lat(), S1Angle.FromDegrees(90));
            assertDoubleNear(better.lng().Radians, S1Angle.FromDegrees(-160).Radians);

            bad = S2LatLng.fromDegrees(-100, -360);
            assertTrue(!bad.isValid());
            better = bad.normalized();
            assertTrue(better.isValid());
            assertEquals(better.lat(), S1Angle.FromDegrees(-90));
            assertDoubleNear(better.lng().Radians, 0);

            assertTrue((S2LatLng.fromDegrees(10, 20).add(S2LatLng.fromDegrees(20, 30))).approxEquals(
                S2LatLng.fromDegrees(30, 50)));
            assertTrue((S2LatLng.fromDegrees(10, 20).sub(S2LatLng.fromDegrees(20, 30))).approxEquals(
                S2LatLng.fromDegrees(-10, -10)));
            assertTrue((S2LatLng.fromDegrees(10, 20).mul(0.5)).approxEquals(S2LatLng.fromDegrees(5, 10)));
        }

        [Test]
        public void testConversion()
        {
            // Test special cases: poles, "date line"
            assertDoubleNear(
                new S2LatLng(S2LatLng.fromDegrees(90.0, 65.0).toPoint()).lat().Degrees, 90.0);
            assertEquals(
                new S2LatLng(S2LatLng.fromRadians(-S2.PiOver2, 1).toPoint()).lat().Radians, -S2.PiOver2);
            assertDoubleNear(
                Math.Abs(new S2LatLng(S2LatLng.fromDegrees(12.2, 180.0).toPoint()).lng().Degrees), 180.0);
            assertEquals(
                Math.Abs(new S2LatLng(S2LatLng.fromRadians(0.1, -S2.Pi).toPoint()).lng().Radians),
                S2.Pi);

            // Test a bunch of random points.
            for (var i = 0; i < 100000; ++i)
            {
                var p = randomPoint();
                assertTrue(S2.ApproxEquals(p, new S2LatLng(p).toPoint()));
            }

            // Test generation from E5
            var test = S2LatLng.fromE5(123456, 98765);
            assertDoubleNear(test.lat().Degrees, 1.23456);
            assertDoubleNear(test.lng().Degrees, 0.98765);
        }

        [Test]
        public void testDistance()
        {
            assertEquals(
                S2LatLng.fromDegrees(90, 0).getDistance(S2LatLng.fromDegrees(90, 0)).Radians, 0.0);
            assertDoubleNear(
                S2LatLng.fromDegrees(-37, 25).getDistance(S2LatLng.fromDegrees(-66, -155)).Degrees, 77,
                1e-13);
            assertDoubleNear(
                S2LatLng.fromDegrees(0, 165).getDistance(S2LatLng.fromDegrees(0, -80)).Degrees, 115,
                1e-13);
            assertDoubleNear(
                S2LatLng.fromDegrees(47, -127).getDistance(S2LatLng.fromDegrees(-47, 53)).Degrees, 180,
                2e-6);
        }
    }
}