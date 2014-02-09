using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Google.Common.Geometry;
using NUnit.Framework;

namespace S2Geometry.Tests
{
    public class S2CellUnionTest : GeometryTestCase
    {
        // Test helper methods for testing the traversal order.
        private static int swapAxes(int ij)
        {
            return ((ij >> 1) & 1) + ((ij & 1) << 1);
        }

        private static int invertBits(int ij)
        {
            return ij ^ 3;
        }

        // Note: obviously, I could have defined a bundle of metrics like this in the
        // S2 class itself rather than just for testing. However, it's not clear that
        // this is useful other than for testing purposes, and I find
        // S2.kMinWidth.GetMaxLevel(width) to be slightly more readable than
        // than S2.kWidth.Min().GetMaxLevel(width). Also, there is no fundamental
        // reason that we need to analyze the minimum, maximum, and average values of
        // every metric; it would be perfectly reasonable to just define one of these.

        public class MetricBundle
        {
            public S2.Metric avg_;
            public S2.Metric max_;
            public S2.Metric min_;

            public MetricBundle(S2.Metric Min, S2.Metric Max, S2.Metric avg)
            {
                min_ = Min;
                max_ = Max;
                avg_ = avg;
            }
        }

        public void testMinMaxAvg(MetricBundle bundle)
        {
            assertTrue(bundle.min_.deriv() < bundle.avg_.deriv());
            assertTrue(bundle.avg_.deriv() < bundle.max_.deriv());
        }

        public void testLessOrEqual(MetricBundle a, MetricBundle b)
        {
            assertTrue(a.min_.deriv() <= b.min_.deriv());
            assertTrue(a.max_.deriv() <= b.max_.deriv());
            assertTrue(a.avg_.deriv() <= b.avg_.deriv());
        }

        [Test]
        public void testAngleArea()
        {
            var pz = new S2Point(0, 0, 1);
            var p000 = new S2Point(1, 0, 0);
            var p045 = new S2Point(1, 1, 0);
            var p090 = new S2Point(0, 1, 0);
            var p180 = new S2Point(-1, 0, 0);
            assertDoubleNear(S2.angle(p000, pz, p045), S2.M_PI_4);
            assertDoubleNear(S2.angle(p045, pz, p180), 3*S2.M_PI_4);
            assertDoubleNear(S2.angle(p000, pz, p180), S2.M_PI);
            assertDoubleNear(S2.angle(pz, p000, pz), 0);
            assertDoubleNear(S2.angle(pz, p000, p045), S2.M_PI_2);

            assertDoubleNear(S2.area(p000, p090, pz), S2.M_PI_2);
            assertDoubleNear(S2.area(p045, pz, p180), 3*S2.M_PI_4);

            // Make sure that area() has good *relative* accuracy even for
            // very small areas.
            var eps = 1e-10;
            var pepsx = new S2Point(eps, 0, 1);
            var pepsy = new S2Point(0, eps, 1);
            var expected1 = 0.5*eps*eps;
            assertDoubleNear(S2.area(pepsx, pepsy, pz), expected1, 1e-14*expected1);

            // Make sure that it can handle degenerate triangles.
            var pr = new S2Point(0.257, -0.5723, 0.112);
            var pq = new S2Point(-0.747, 0.401, 0.2235);
            assertEquals(S2.area(pr, pr, pr), 0.0);
            // TODO: The following test is not exact in optimized mode because the
            // compiler chooses to mix 64-bit and 80-bit intermediate results.
            assertDoubleNear(S2.area(pr, pq, pr), 0);
            assertEquals(S2.area(p000, p045, p090), 0.0);

            double maxGirard = 0;
            for (var i = 0; i < 10000; ++i)
            {
                var p0 = randomPoint();
                var d1 = randomPoint();
                var d2 = randomPoint();
                var p1 = S2Point.add(p0, S2Point.mul(d1, 1e-15));
                var p2 = S2Point.add(p0, S2Point.mul(d2, 1e-15));
                // The actual displacement can be as much as 1.2e-15 due to roundoff.
                // This yields a maximum triangle area of about 0.7e-30.
                assertTrue(S2.area(p0, p1, p2) < 0.7e-30);
                maxGirard = Math.Max(maxGirard, S2.girardArea(p0, p1, p2));
            }
            Console.WriteLine("Worst case Girard for triangle area 1e-30: " + maxGirard);

            // Try a very long and skinny triangle.
            var p045eps = new S2Point(1, 1, eps);
            var expected2 = 5.8578643762690495119753e-11; // Mathematica.
            assertDoubleNear(S2.area(p000, p045eps, p090), expected2, 1e-9*expected2);

            // Triangles with near-180 degree edges that sum to a quarter-sphere.
            var eps2 = 1e-10;
            var p000eps2 = new S2Point(1, 0.1*eps2, eps2);
            var quarterArea1 =
                S2.area(p000eps2, p000, p090) + S2.area(p000eps2, p090, p180) + S2.area(p000eps2, p180, pz)
                + S2.area(p000eps2, pz, p000);
            assertDoubleNear(quarterArea1, S2.M_PI);

            // Four other triangles that sum to a quarter-sphere.
            var p045eps2 = new S2Point(1, 1, eps2);
            var quarterArea2 =
                S2.area(p045eps2, p000, p090) + S2.area(p045eps2, p090, p180) + S2.area(p045eps2, p180, pz)
                + S2.area(p045eps2, pz, p000);
            assertDoubleNear(quarterArea2, S2.M_PI);
        }

        [Test]
        public void testCCW()
        {
            var a = new S2Point(0.72571927877036835, 0.46058825605889098, 0.51106749730504852);
            var b = new S2Point(0.7257192746638208, 0.46058826573818168, 0.51106749441312738);
            var c = new S2Point(0.72571927671709457, 0.46058826089853633, 0.51106749585908795);
            assertTrue(S2.robustCCW(a, b, c) != 0);
        }

        [Test]
        public void testExp()
        {
            for (var i = 0; i < 10; ++i)
            {
                assertEquals(i + 1, S2.exp(Math.Pow(2, i)));
            }

            for (var i = 0; i < 10; ++i)
            {
                assertEquals(i + 1, S2.exp(-Math.Pow(2, i)));
            }

            assertEquals(0, S2.exp(0));
            assertEquals(2, S2.exp(3));
            assertEquals(3, S2.exp(5));
        }

        [Test]
        public void testFaceUVtoXYZ()
        {
            // Check that each face appears exactly once.
            var sum = new S2Point();
            for (var face = 0; face < 6; ++face)
            {
                var center = S2Projections.faceUvToXyz(face, 0, 0);
                assertEquals(S2Projections.getNorm(face), center);
                assertEquals(Math.Abs(center.get(center.largestAbsComponent())), 1.0);
                sum = S2Point.add(sum, S2Point.fabs(center));
            }
            assertEquals(sum, new S2Point(2, 2, 2));

            // Check that each face has a right-handed coordinate system.
            for (var face = 0; face < 6; ++face)
            {
                assertEquals(
                    S2Point.crossProd(S2Projections.getUAxis(face), S2Projections.getVAxis(face)).dotProd(
                        S2Projections.faceUvToXyz(face, 0, 0)), 1.0);
            }

            // Check that the Hilbert curves on each face combine to form a
            // continuous curve over the entire cube.
            for (var face = 0; face < 6; ++face)
            {
                // The Hilbert curve on each face starts at (-1,-1) and terminates
                // at either (1,-1) (if axes not swapped) or (-1,1) (if swapped).
                var sign = ((face & S2.SWAP_MASK) != 0) ? -1 : 1;
                assertEquals(S2Projections.faceUvToXyz(face, sign, -sign),
                             S2Projections.faceUvToXyz((face + 1)%6, -1, -1));
            }
        }

        [Test]
        public void testMetrics()
        {
            var angleSpan = new MetricBundle(
                S2Projections.MIN_ANGLE_SPAN, S2Projections.MAX_ANGLE_SPAN, S2Projections.AVG_ANGLE_SPAN);
            var width =
                new MetricBundle(S2Projections.MIN_WIDTH, S2Projections.MAX_WIDTH, S2Projections.AVG_WIDTH);
            var edge =
                new MetricBundle(S2Projections.MIN_EDGE, S2Projections.MAX_EDGE, S2Projections.AVG_EDGE);
            var diag =
                new MetricBundle(S2Projections.MIN_DIAG, S2Projections.MAX_DIAG, S2Projections.AVG_DIAG);
            var area =
                new MetricBundle(S2Projections.MIN_AREA, S2Projections.MAX_AREA, S2Projections.AVG_AREA);

            // First, check that Min <= avg <= Max for each metric.
            testMinMaxAvg(angleSpan);
            testMinMaxAvg(width);
            testMinMaxAvg(edge);
            testMinMaxAvg(diag);
            testMinMaxAvg(area);

            // Check that the maximum aspect ratio of an individual cell is consistent
            // with the global minimums and maximums.
            assertTrue(S2Projections.MAX_EDGE_ASPECT >= 1.0);
            assertTrue(S2Projections.MAX_EDGE_ASPECT
                       < S2Projections.MAX_EDGE.deriv()/S2Projections.MIN_EDGE.deriv());
            assertTrue(S2Projections.MAX_DIAG_ASPECT >= 1);
            assertTrue(S2Projections.MAX_DIAG_ASPECT
                       < S2Projections.MAX_DIAG.deriv()/S2Projections.MIN_DIAG.deriv());

            // Check various conditions that are provable mathematically.
            testLessOrEqual(width, angleSpan);
            testLessOrEqual(width, edge);
            testLessOrEqual(edge, diag);

            assertTrue(S2Projections.MIN_AREA.deriv()
                       >= S2Projections.MIN_WIDTH.deriv()*S2Projections.MIN_EDGE.deriv() - 1e-15);
            assertTrue(S2Projections.MAX_AREA.deriv()
                       < S2Projections.MAX_WIDTH.deriv()*S2Projections.MAX_EDGE.deriv() + 1e-15);

            // GetMinLevelForLength() and friends have built-in assertions, we just need
            // to call these functions to test them.
            //
            // We don't actually check that the metrics are correct here, e.g. that
            // GetMinWidth(10) is a lower bound on the width of cells at level 10.
            // It is easier to check these properties in s2cell_unittest, since
            // S2Cell has methods to compute the cell vertices, etc.

            for (var level = -2; level <= S2CellId.MAX_LEVEL + 3; ++level)
            {
                var dWidth = (2*S2Projections.MIN_WIDTH.deriv())*Math.Pow(2, -level);
                if (level >= S2CellId.MAX_LEVEL + 3)
                {
                    dWidth = 0;
                }

                // Check boundary cases (exactly equal to a threshold value).
                var expectedLevel = Math.Max(0, Math.Min(S2CellId.MAX_LEVEL, level));
                assertEquals(S2Projections.MIN_WIDTH.getMinLevel(dWidth), expectedLevel);
                assertEquals(S2Projections.MIN_WIDTH.getMaxLevel(dWidth), expectedLevel);
                assertEquals(S2Projections.MIN_WIDTH.getClosestLevel(dWidth), expectedLevel);

                // Also check non-boundary cases.
                assertEquals(S2Projections.MIN_WIDTH.getMinLevel(1.2*dWidth), expectedLevel);
                assertEquals(S2Projections.MIN_WIDTH.getMaxLevel(0.8*dWidth), expectedLevel);
                assertEquals(S2Projections.MIN_WIDTH.getClosestLevel(1.2*dWidth), expectedLevel);
                assertEquals(S2Projections.MIN_WIDTH.getClosestLevel(0.8*dWidth), expectedLevel);

                // Same thing for area1.
                var area1 = (4*S2Projections.MIN_AREA.deriv())*Math.Pow(4, -level);
                if (level <= -3)
                {
                    area1 = 0;
                }
                assertEquals(S2Projections.MIN_AREA.getMinLevel(area1), expectedLevel);
                assertEquals(S2Projections.MIN_AREA.getMaxLevel(area1), expectedLevel);
                assertEquals(S2Projections.MIN_AREA.getClosestLevel(area1), expectedLevel);
                assertEquals(S2Projections.MIN_AREA.getMinLevel(1.2*area1), expectedLevel);
                assertEquals(S2Projections.MIN_AREA.getMaxLevel(0.8*area1), expectedLevel);
                assertEquals(S2Projections.MIN_AREA.getClosestLevel(1.2*area1), expectedLevel);
                assertEquals(S2Projections.MIN_AREA.getClosestLevel(0.8*area1), expectedLevel);
            }
        }

        [Test]
        public void testSTUV()
        {
            // Check boundary conditions.
            for (double x = -1; x <= 1; ++x)
            {
                assertEquals(S2Projections.stToUV(x), x);
                assertEquals(S2Projections.uvToST(x), x);
            }
            // Check that UVtoST and STtoUV are inverses.
            for (double x = -1; x <= 1; x += 0.0001)
            {
                assertDoubleNear(S2Projections.uvToST(S2Projections.stToUV(x)), x);
                assertDoubleNear(S2Projections.stToUV(S2Projections.uvToST(x)), x);
            }
        }

        [Test]
        public void testTraversalOrder()
        {
            for (var r = 0; r < 4; ++r)
            {
                for (var i = 0; i < 4; ++i)
                {
                    // Check consistency with respect to swapping axes.
                    assertEquals(S2.ijToPos(r, i),
                                 S2.ijToPos(r ^ S2.SWAP_MASK, swapAxes(i)));
                    assertEquals(S2.posToIJ(r, i),
                                 swapAxes(S2.posToIJ(r ^ S2.SWAP_MASK, i)));

                    // Check consistency with respect to reversing axis directions.
                    assertEquals(S2.ijToPos(r, i),
                                 S2.ijToPos(r ^ S2.INVERT_MASK, invertBits(i)));
                    assertEquals(S2.posToIJ(r, i),
                                 invertBits(S2.posToIJ(r ^ S2.INVERT_MASK, i)));

                    // Check that the two tables are inverses of each other.
                    assertEquals(S2.ijToPos(r, S2.posToIJ(r, i)), i);
                    assertEquals(S2.posToIJ(r, S2.ijToPos(r, i)), i);
                }
            }
        }

        [Test]
        public void testUVAxes()
        {
            // Check that axes are consistent with FaceUVtoXYZ.
            for (var face = 0; face < 6; ++face)
            {
                assertEquals(S2Projections.getUAxis(face), S2Point.sub(
                    S2Projections.faceUvToXyz(face, 1, 0), S2Projections.faceUvToXyz(face, 0, 0)));
                assertEquals(S2Projections.getVAxis(face), S2Point.sub(
                    S2Projections.faceUvToXyz(face, 0, 1), S2Projections.faceUvToXyz(face, 0, 0)));
            }
        }

        [Test]
        public void testUVNorms()
        {
            // Check that GetUNorm and GetVNorm compute right-handed normals for
            // an edge in the increasing U or V direction.
            for (var face = 0; face < 6; ++face)
            {
                for (double x = -1; x <= 1; x += 1/1024.0)
                {
                    assertDoubleNear(
                        S2Point.crossProd(
                            S2Projections.faceUvToXyz(face, x, -1), S2Projections.faceUvToXyz(face, x, 1))
                               .angle(S2Projections.getUNorm(face, x)), 0);
                    assertDoubleNear(
                        S2Point.crossProd(
                            S2Projections.faceUvToXyz(face, -1, x), S2Projections.faceUvToXyz(face, 1, x))
                               .angle(S2Projections.getVNorm(face, x)), 0);
                }
            }
        }
    }
}