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
            JavaAssert.Equal(S1Angle.radians(Math.PI).radians(), Math.PI);
            JavaAssert.Equal(S1Angle.radians(Math.PI).degrees(), 180.0);
            JavaAssert.Equal(S1Angle.degrees(180).radians(), Math.PI);
            JavaAssert.Equal(S1Angle.degrees(180).degrees(), 180.0);

            JavaAssert.Equal(S1Angle.radians(Math.PI/2).degrees(), 90.0);

            // Check negative angles.
            JavaAssert.Equal(S1Angle.radians(-Math.PI/2).degrees(), -90.0);
            JavaAssert.Equal(S1Angle.degrees(-45).radians(), -Math.PI/4);

            // Check that E5/E6/E7 representations work as expected.
            JavaAssert.Equal(S1Angle.e5(2000000), S1Angle.degrees(20));
            JavaAssert.Equal(S1Angle.e6(-60000000), S1Angle.degrees(-60));
            JavaAssert.Equal(S1Angle.e7(750000000), S1Angle.degrees(75));
            JavaAssert.Equal(S1Angle.degrees(12.34567).e5(), 1234567);
            JavaAssert.Equal(S1Angle.degrees(12.345678).e6(), 12345678);
            JavaAssert.Equal(S1Angle.degrees(-12.3456789).e7(), -123456789);
        }
    }
}