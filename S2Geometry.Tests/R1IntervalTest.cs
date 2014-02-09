using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Google.Common.Geometry;
using NUnit.Framework;

namespace S2Geometry.Tests
{
    public class R1IntervalTest : GeometryTestCase
    {
        /**
  * Test all of the interval operations on the given pair of intervals.
  * "expected_relation" is a sequence of "T" and "F" characters corresponding
  * to the expected results of contains(), interiorContains(), Intersects(),
  * and InteriorIntersects() respectively.
  */

        private void testIntervalOps(R1Interval x, R1Interval y, String expectedRelation)
        {
            JavaAssert.Equal(x.Contains(y), expectedRelation[0] == 'T');
            JavaAssert.Equal(x.InteriorContains(y), expectedRelation[1] == 'T');
            JavaAssert.Equal(x.Intersects(y), expectedRelation[2] == 'T');
            JavaAssert.Equal(x.InteriorIntersects(y), expectedRelation[3] == 'T');

            JavaAssert.Equal(x.Contains(y), x.Union(y).Equals(x));
            JavaAssert.Equal(x.Intersects(y), !x.Intersection(y).IsEmpty);
        }

        [Test]
        public void R1IntervalBasicTest()
        {
            // Constructors and accessors.
            var unit = new R1Interval(0, 1);
            var negunit = new R1Interval(-1, 0);
            JavaAssert.Equal(unit.Lo, 0.0);
            JavaAssert.Equal(unit.Hi, 1.0);
            JavaAssert.Equal(negunit.Lo, -1.0);
            JavaAssert.Equal(negunit.Hi, 0.0);

            // is_empty()
            var half = new R1Interval(0.5, 0.5);
            Assert.True(!unit.IsEmpty);
            Assert.True(!half.IsEmpty);
            var empty = R1Interval.Empty;
            Assert.True(empty.IsEmpty);

            // GetCenter(), GetLength()
            JavaAssert.Equal(unit.Center, 0.5);
            JavaAssert.Equal(half.Center, 0.5);
            JavaAssert.Equal(negunit.Length, 1.0);
            JavaAssert.Equal(half.Length, 0.0);
            Assert.True(empty.Length < 0);

            // contains(double), interiorContains(double)
            Assert.True(unit.Contains(0.5));
            Assert.True(unit.InteriorContains(0.5));
            Assert.True(unit.Contains(0));
            Assert.True(!unit.InteriorContains(0));
            Assert.True(unit.Contains(1));
            Assert.True(!unit.InteriorContains(1));

            // contains(R1Interval), interiorContains(R1Interval)
            // Intersects(R1Interval), InteriorIntersects(R1Interval)
            testIntervalOps(empty, empty, "TTFF");
            testIntervalOps(empty, unit, "FFFF");
            testIntervalOps(unit, half, "TTTT");
            testIntervalOps(unit, unit, "TFTT");
            testIntervalOps(unit, empty, "TTFF");
            testIntervalOps(unit, negunit, "FFTF");
            testIntervalOps(unit, new R1Interval(0, 0.5), "TFTT");
            testIntervalOps(half, new R1Interval(0, 0.5), "FFTF");

            // addPoint()
            R1Interval r;
            r = empty.AddPoint(5);
            Assert.True(r.Lo == 5.0 && r.Hi == 5.0);
            r = r.AddPoint(-1);
            Assert.True(r.Lo == -1.0 && r.Hi == 5.0);
            r = r.AddPoint(0);
            Assert.True(r.Lo == -1.0 && r.Hi == 5.0);

            // fromPointPair()
            JavaAssert.Equal(R1Interval.FromPointPair(4, 4), new R1Interval(4, 4));
            JavaAssert.Equal(R1Interval.FromPointPair(-1, -2), new R1Interval(-2, -1));
            JavaAssert.Equal(R1Interval.FromPointPair(-5, 3), new R1Interval(-5, 3));

            // expanded()
            JavaAssert.Equal(empty.Expanded(0.45), empty);
            JavaAssert.Equal(unit.Expanded(0.5), new R1Interval(-0.5, 1.5));

            // union(), intersection()
            Assert.True(new R1Interval(99, 100).Union(empty).Equals(new R1Interval(99, 100)));
            Assert.True(empty.Union(new R1Interval(99, 100)).Equals(new R1Interval(99, 100)));
            Assert.True(new R1Interval(5, 3).Union(new R1Interval(0, -2)).IsEmpty);
            Assert.True(new R1Interval(0, -2).Union(new R1Interval(5, 3)).IsEmpty);
            Assert.True(unit.Union(unit).Equals(unit));
            Assert.True(unit.Union(negunit).Equals(new R1Interval(-1, 1)));
            Assert.True(negunit.Union(unit).Equals(new R1Interval(-1, 1)));
            Assert.True(half.Union(unit).Equals(unit));
            Assert.True(unit.Intersection(half).Equals(half));
            Assert.True(unit.Intersection(negunit).Equals(new R1Interval(0, 0)));
            Assert.True(negunit.Intersection(half).IsEmpty);
            Assert.True(unit.Intersection(empty).IsEmpty);
            Assert.True(empty.Intersection(unit).IsEmpty);
        }
    }
}