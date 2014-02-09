using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Google.Common.Geometry;
using NUnit.Framework;

namespace S2Geometry.Tests
{
    [TestFixture]
    public class S1AngleTest
    {
        [Test]
        public void S1AngleBasicTest()
        {
            // Check that the conversion between Pi radians and 180 degrees is exact.
            JavaAssert.Equal(S1Angle.FromRadians(Math.PI).Radians, Math.PI);
            JavaAssert.Equal(S1Angle.FromRadians(Math.PI).Degrees, 180.0);
            JavaAssert.Equal(S1Angle.FromDegrees(180).Radians, Math.PI);
            JavaAssert.Equal(S1Angle.FromDegrees(180).Degrees, 180.0);

            JavaAssert.Equal(S1Angle.FromRadians(Math.PI/2).Degrees, 90.0);

            // Check negative angles.
            JavaAssert.Equal(S1Angle.FromRadians(-Math.PI/2).Degrees, -90.0);
            JavaAssert.Equal(S1Angle.FromDegrees(-45).Radians, -Math.PI/4);

            // Check that E5/E6/E7 representations work as expected.
            JavaAssert.Equal(S1Angle.E5(2000000), S1Angle.FromDegrees(20));
            JavaAssert.Equal(S1Angle.E6(-60000000), S1Angle.FromDegrees(-60));
            JavaAssert.Equal(S1Angle.E7(750000000), S1Angle.FromDegrees(75));
            JavaAssert.Equal(S1Angle.FromDegrees(12.34567).E5(), 1234567);
            JavaAssert.Equal(S1Angle.FromDegrees(12.345678).E6(), 12345678);
            JavaAssert.Equal(S1Angle.FromDegrees(-12.3456789).E7(), -123456789);
        }
    }
}