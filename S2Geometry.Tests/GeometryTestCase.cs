using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Google.Common.Geometry;
using NUnit.Framework;

namespace S2Geometry.Tests
{
    [TestFixture]
    public abstract class GeometryTestCase
    {
        public Random rand { get; private set; }

        [TestFixtureSetUp]
        protected virtual void SetUp()
        {
            rand = new Random(123456);
        }

        protected long LongRandom()
        {
            // This is how Java's nextLong works
            var bytes = new byte[4];
            rand.NextBytes(bytes);
            var bytes1 = new byte[4];
            rand.NextBytes(bytes1);

            var i1 = BitConverter.ToInt32(bytes, 0);
            var i2 = BitConverter.ToInt32(bytes1, 0);
            return ((long)(i1 << 32)) + i2;
        }

        public void assertDoubleNear(double a, double b)
        {
            assertDoubleNear(a, b, 1e-9);
        }

        public void assertDoubleNear(double a, double b, double error)
        {
            Assert.True(a + error > b);
            Assert.True(a < b + error);
        }

        // maybe these should be put in a special testing util class
        /** Return a random unit-length vector. */

        public S2Point randomPoint()
        {
            return S2Point.Normalize(new S2Point(
                                         2*rand.NextDouble() - 1,
                                         2*rand.NextDouble() - 1,
                                         2*rand.NextDouble() - 1));
        }


        /**
         * Return a right-handed coordinate frame (three orthonormal vectors). Returns
         * an array of three points: x,y,z
         */

        public IReadOnlyList<S2Point> getRandomFrame()
        {
            var p0 = randomPoint();
            var p1 = S2Point.Normalize(S2Point.CrossProd(p0, randomPoint()));
            var p2 = S2Point.Normalize(S2Point.CrossProd(p0, p1));
            return new List<S2Point>(new[] {p0, p1, p2});
        }

        /**
         * Return a random cell id at the given level or at a randomly chosen level.
         * The distribution is uniform over the space of cell ids, but only
         * approximately uniform over the surface of the sphere.
         */

        public S2CellId getRandomCellId(int level)
        {
            var face = random(S2CellId.NUM_FACES);

            var pos = (ulong)LongRandom() & ((1L << (2*S2CellId.MAX_LEVEL)) - 1);
            return S2CellId.fromFacePosLevel(face, pos, level);
        }

        public S2CellId getRandomCellId()
        {
            return getRandomCellId(random(S2CellId.MAX_LEVEL + 1));
        }

        protected int random(int n)
        {
            if (n == 0)
            {
                return 0;
            }
            return rand.Next(n);
        }


        // Pick "base" uniformly from range [0,maxLog] and then return
        // "base" random bits. The effect is to pick a number in the range
        // [0,2^maxLog-1] with bias towards smaller numbers.
        protected int skewed(int maxLog)
        {
            var @base = Math.Abs(rand.Next())%(maxLog + 1);
            // if (!base) return 0; // if 0==base, we & with 0 below.
            //
            // this distribution differs slightly from ACMRandom's Skewed,
            // since 0 occurs approximately 3 times more than 1 here, and
            // ACMRandom's Skewed never outputs 0.
            return rand.Next() & ((1 << @base) - 1);
        }

        /**
         * Checks that "covering" completely covers the given region. If "check_tight"
         * is true, also checks that it does not contain any cells that do not
         * intersect the given region. ("id" is only used internally.)
         */

        protected void checkCovering(IS2Region region, S2CellUnion covering, bool checkTight, S2CellId id)
        {
            if (!id.isValid())
            {
                for (var face = 0; face < 6; ++face)
                {
                    checkCovering(region, covering, checkTight, S2CellId.fromFacePosLevel(face, 0, 0));
                }
                return;
            }

            if (!region.MayIntersect(new S2Cell(id)))
            {
                // If region does not intersect id, then neither should the covering.
                if (checkTight)
                {
                    Assert.True(!covering.intersects(id));
                }
            }
            else if (!covering.contains(id))
            {
                // The region may intersect id, but we can't assert that the covering
                // intersects id because we may discover that the region does not actually
                // intersect upon further subdivision. (MayIntersect is not exact.)
                Assert.True(!region.Contains(new S2Cell(id)));
                var result = !id.isLeaf();
                Assert.True(result);
                var end = id.childEnd();
                for (var child = id.childBegin(); !child.Equals(end); child = child.next())
                {
                    checkCovering(region, covering, checkTight, child);
                }
            }
        }

        protected S2Cap getRandomCap(double minArea, double maxArea)
        {
            var capArea = maxArea
                          *Math.Pow(minArea/maxArea, rand.NextDouble());
            Assert.True(capArea >= minArea && capArea <= maxArea);

            // The surface area of a cap is 2*Pi times its height.
            return S2Cap.fromAxisArea(randomPoint(), capArea);
        }

        protected S2Point samplePoint(S2Cap cap)
        {
            // We consider the cap axis to be the "z" axis. We choose two other axes to
            // complete the coordinate frame.

            var z = cap.axis();
            var x = z.Ortho;
            var y = S2Point.CrossProd(z, x);

            // The surface area of a spherical cap is directly proportional to its
            // height. First we choose a random height, and then we choose a random
            // point along the circle at that height.

            var h = rand.NextDouble()*cap.height();
            var theta = 2*S2.Pi*rand.NextDouble();
            var r = Math.Sqrt(h*(2 - h)); // Radius of circle.

            // (cos(theta)*r*x + sin(theta)*r*y + (1-h)*z).Normalize()
            return S2Point.Normalize(((x * Math.Cos(theta)*r) + (y * Math.Sin(theta)*r)) + (z * (1 - h)));
        }

        private static void parseVertices(String str, List<S2Point> vertices)
        {
            if (str == null)
            {
                return;
            }

            foreach (var token in str.Split(','))
            {
                var colon = token.IndexOf(':');
                if (colon == -1)
                {
                    throw new ArgumentException(
                        "Illegal string:" + token + ". Should look like '35:20'");
                }
                var lat = Double.Parse(token.Substring(0, colon));
                var lng = Double.Parse(token.Substring(colon + 1));
                vertices.Add(S2LatLng.fromDegrees(lat, lng).toPoint());
            }
        }

        protected static S2Point makePoint(String str)
        {
            var vertices = new List<S2Point>();
            parseVertices(str, vertices);
            return vertices.Single();
        }

        protected static S2Loop makeLoop(String str)
        {
            var vertices = new List<S2Point>();
            parseVertices(str, vertices);
            return new S2Loop(vertices);
        }

        protected static S2Polygon makePolygon(String str)
        {
            var loops = new List<S2Loop>();

            foreach (var token in str.Split(new[] {';'}, StringSplitOptions.RemoveEmptyEntries))
            {
//Splitter.on(';').omitEmptyStrings().split(str)) {
                var loop = makeLoop(token);
                loop.normalize();
                loops.Add(loop);
            }

            return new S2Polygon(loops);
        }

        [DebuggerNonUserCode]
        [DebuggerStepThrough]
        protected static void assertEquals(object actual, object expected)
        {
            JavaAssert.Equal(actual, expected);
        }

        [DebuggerNonUserCode]
        [DebuggerStepThrough]
        protected static void assertEquals(double actual, double expected, double delta)
        {
            Assert.AreEqual(expected, actual, delta);
        }

        [DebuggerNonUserCode]
        [DebuggerStepThrough]
        protected static void assertTrue(bool value)
        {
            Assert.True(value);
        }

        [DebuggerNonUserCode]
        [DebuggerStepThrough]
        protected static void assertFalse(bool value)
        {
            Assert.False(value);
        }

        [DebuggerNonUserCode]
        [DebuggerStepThrough]
        protected static void assertFalse(string message, bool value)
        {
            Assert.False(value, message);
        }

        [DebuggerNonUserCode]
        [DebuggerStepThrough]
        protected static void assertTrue(string message, bool value)
        {
            Assert.True(value, message);
        }

        protected static S2Polyline makePolyline(String str)
        {
            var vertices = new List<S2Point>();
            parseVertices(str, vertices);
            return new S2Polyline(vertices);
        }
    }
}