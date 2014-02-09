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
            JavaAssert.Equal(x.contains(y), expectedRelation[0] == 'T');
            JavaAssert.Equal(x.interiorContains(y), expectedRelation[1] == 'T');
            JavaAssert.Equal(x.intersects(y), expectedRelation[2] == 'T');
            JavaAssert.Equal(x.interiorIntersects(y), expectedRelation[3] == 'T');

            JavaAssert.Equal(x.contains(y), x.union(y).Equals(x));
            JavaAssert.Equal(x.intersects(y), !x.intersection(y).isEmpty());
        }

        [Test]
        public void R1IntervalBasicTest()
        {
            // Constructors and accessors.
            var unit = new R1Interval(0, 1);
            var negunit = new R1Interval(-1, 0);
            JavaAssert.Equal(unit.lo(), 0.0);
            JavaAssert.Equal(unit.hi(), 1.0);
            JavaAssert.Equal(negunit.lo(), -1.0);
            JavaAssert.Equal(negunit.hi(), 0.0);

            // is_empty()
            var half = new R1Interval(0.5, 0.5);
            Assert.True(!unit.isEmpty());
            Assert.True(!half.isEmpty());
            var empty = R1Interval.empty();
            Assert.True(empty.isEmpty());

            // GetCenter(), GetLength()
            JavaAssert.Equal(unit.getCenter(), 0.5);
            JavaAssert.Equal(half.getCenter(), 0.5);
            JavaAssert.Equal(negunit.getLength(), 1.0);
            JavaAssert.Equal(half.getLength(), 0.0);
            Assert.True(empty.getLength() < 0);

            // contains(double), interiorContains(double)
            Assert.True(unit.contains(0.5));
            Assert.True(unit.interiorContains(0.5));
            Assert.True(unit.contains(0));
            Assert.True(!unit.interiorContains(0));
            Assert.True(unit.contains(1));
            Assert.True(!unit.interiorContains(1));

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
            r = empty.addPoint(5);
            Assert.True(r.lo() == 5.0 && r.hi() == 5.0);
            r = r.addPoint(-1);
            Assert.True(r.lo() == -1.0 && r.hi() == 5.0);
            r = r.addPoint(0);
            Assert.True(r.lo() == -1.0 && r.hi() == 5.0);

            // fromPointPair()
            JavaAssert.Equal(R1Interval.fromPointPair(4, 4), new R1Interval(4, 4));
            JavaAssert.Equal(R1Interval.fromPointPair(-1, -2), new R1Interval(-2, -1));
            JavaAssert.Equal(R1Interval.fromPointPair(-5, 3), new R1Interval(-5, 3));

            // expanded()
            JavaAssert.Equal(empty.expanded(0.45), empty);
            JavaAssert.Equal(unit.expanded(0.5), new R1Interval(-0.5, 1.5));

            // union(), intersection()
            Assert.True(new R1Interval(99, 100).union(empty).Equals(new R1Interval(99, 100)));
            Assert.True(empty.union(new R1Interval(99, 100)).Equals(new R1Interval(99, 100)));
            Assert.True(new R1Interval(5, 3).union(new R1Interval(0, -2)).isEmpty());
            Assert.True(new R1Interval(0, -2).union(new R1Interval(5, 3)).isEmpty());
            Assert.True(unit.union(unit).Equals(unit));
            Assert.True(unit.union(negunit).Equals(new R1Interval(-1, 1)));
            Assert.True(negunit.union(unit).Equals(new R1Interval(-1, 1)));
            Assert.True(half.union(unit).Equals(unit));
            Assert.True(unit.intersection(half).Equals(half));
            Assert.True(unit.intersection(negunit).Equals(new R1Interval(0, 0)));
            Assert.True(negunit.intersection(half).isEmpty());
            Assert.True(unit.intersection(empty).isEmpty());
            Assert.True(empty.intersection(unit).isEmpty());
        }
    }
}