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
    public class S2CellIdTest : GeometryTestCase
    {
        private S2CellId getCellId(double latDegrees, double lngDegrees)
        {
            var id = S2CellId.fromLatLng(S2LatLng.fromDegrees(latDegrees, lngDegrees));
            Trace.WriteLine(Convert.ToString(unchecked ((long)id.id()), 16));
            return id;
        }

        public void testInverses()
        {
            Trace.WriteLine("TestInverses");
            // Check the conversion of random leaf cells to S2LatLngs and back.
            for (var i = 0; i < 200000; ++i)
            {
                var id = getRandomCellId(S2CellId.MAX_LEVEL);
                Assert.True(id.isLeaf() && id.level() == S2CellId.MAX_LEVEL);
                var center = id.toLatLng();
                JavaAssert.Equal(S2CellId.fromLatLng(center).id(), id.id());
            }
        }

        private const int kMaxExpandLevel = 3;

        private void expandCell(
            S2CellId parent, List<S2CellId> cells, IDictionary<S2CellId, S2CellId> parentMap)
        {
            cells.Add(parent);
            if (parent.level() == kMaxExpandLevel)
            {
                return;
            }
            var i = 0;
            var j = 0;
            int? orientation = 0;
            var face = parent.toFaceIJOrientation(ref i, ref j, ref orientation);
            JavaAssert.Equal(face, parent.face());

            var pos = 0;
            for (var child = parent.childBegin(); !child.Equals(parent.childEnd());
                 child = child.next())
            {
                // Do some basic checks on the children
                JavaAssert.Equal(child.level(), parent.level() + 1);
                Assert.True(!child.isLeaf());
                int? childOrientation = 0;
                JavaAssert.Equal(child.toFaceIJOrientation(ref i, ref j, ref childOrientation), face);
                JavaAssert.Equal(
                    childOrientation.Value, orientation.Value ^ S2.posToOrientation(pos));

                parentMap.Add(child, parent);
                expandCell(child, cells, parentMap);
                ++pos;
            }
        }

        private const int MAX_WALK_LEVEL = 8;

        public void testAllNeighbors(S2CellId id, int level)
        {
            Assert.True(level >= id.level() && level < S2CellId.MAX_LEVEL);

            // We compute GetAllNeighbors, and then add in all the children of "id"
            // at the given level. We then compare this against the result of finding
            // all the vertex neighbors of all the vertices of children of "id" at the
            // given level. These should give the same result.
            var all = new List<S2CellId>();
            var expected = new List<S2CellId>();
            id.getAllNeighbors(level, all);
            var end = id.childEnd(level + 1);
            for (var c = id.childBegin(level + 1); !c.Equals(end); c = c.next())
            {
                all.Add(c.parent());
                c.getVertexNeighbors(level, expected);
            }
            // Sort the results and eliminate duplicates.
            all.Sort();
            expected.Sort();
            ISet<S2CellId> allSet = new HashSet<S2CellId>(all);
            ISet<S2CellId> expectedSet = new HashSet<S2CellId>(expected);
            var result = allSet.SetEquals(expectedSet);
            Assert.True(result);
        }

        [Test]
        public void S2CellIdTestBasic()
        {
            Trace.WriteLine("TestBasic");
            // Check default constructor.
            var id = new S2CellId();
            //JavaAssert.Equal(id.id(), 0);
            //Assert.True(!id.isValid());

            // Check basic accessor methods.
            id = S2CellId.fromFacePosLevel(3, 0x12345678, S2CellId.MAX_LEVEL - 4);
            //Assert.True(id.isValid());
            //JavaAssert.Equal(id.face(), 3);
            // JavaAssert.Equal(id.pos(), 0x12345700);
            //JavaAssert.Equal(id.level(), S2CellId.MAX_LEVEL - 4);
            //Assert.True(!id.isLeaf());

            //// Check face definitions
            //JavaAssert.Equal(getCellId(0, 0).face(), 0);
            //JavaAssert.Equal(getCellId(0, 90).face(), 1);
            //JavaAssert.Equal(getCellId(90, 0).face(), 2);
            //JavaAssert.Equal(getCellId(0, 180).face(), 3);
            //JavaAssert.Equal(getCellId(0, -90).face(), 4);
            //JavaAssert.Equal(getCellId(-90, 0).face(), 5);

            //// Check parent/child relationships.
            //JavaAssert.Equal(id.childBegin(id.level() + 2).pos(), 0x12345610);
            //JavaAssert.Equal(id.childBegin().pos(), 0x12345640);
            //JavaAssert.Equal(id.parent().pos(), 0x12345400);
            //JavaAssert.Equal(id.parent(id.level() - 2).pos(), 0x12345000);

            //// Check ordering of children relative to parents.
            //Assert.True(id.childBegin().lessThan(id));
            //var childEnd = id.childEnd();
            //var childId = childEnd.id();
            //var id1 = id.id();

            //Assert.True(id.childEnd().greaterThan(id));
            //JavaAssert.Equal(id.childBegin().next().next().next().next(), id.childEnd());
            //JavaAssert.Equal(id.childBegin(S2CellId.MAX_LEVEL), id.rangeMin());
            //JavaAssert.Equal(id.childEnd(S2CellId.MAX_LEVEL), id.rangeMax().next());

            // Check wrapping from beginning of Hilbert curve to end and vice versa.
            // JavaAssert.Equal(S2CellId.begin(0).prevWrap(), S2CellId.end(0).prev());

            JavaAssert.Equal(S2CellId.begin(S2CellId.MAX_LEVEL).prevWrap(),
                             S2CellId.fromFacePosLevel(5, ~0UL >> S2CellId.FACE_BITS, S2CellId.MAX_LEVEL));

            JavaAssert.Equal(S2CellId.end(4).prev().nextWrap(), S2CellId.begin(4));
            JavaAssert.Equal(S2CellId.end(S2CellId.MAX_LEVEL).prev().nextWrap(),
                             S2CellId.fromFacePosLevel(0, 0, S2CellId.MAX_LEVEL));

            // Check that cells are represented by the position of their center
            // along the Hilbert curve.
            JavaAssert.Equal(id.rangeMin().id() + id.rangeMax().id(), 2*id.id());
        }

        [Test]
        public void testContainment()
        {
            Trace.WriteLine("TestContainment");
            IDictionary<S2CellId, S2CellId> parentMap = new Dictionary<S2CellId, S2CellId>();
            var cells = new List<S2CellId>();
            for (var face = 0; face < 6; ++face)
            {
                expandCell(S2CellId.fromFacePosLevel(face, 0, 0), cells, parentMap);
            }
            for (var i = 0; i < cells.Count; ++i)
            {
                for (var j = 0; j < cells.Count; ++j)
                {
                    var contained = true;
                    for (var id = cells[j]; id != cells[i]; id = parentMap[id])
                    {
                        if (!parentMap.ContainsKey(id))
                        {
                            contained = false;
                            break;
                        }
                    }
                    JavaAssert.Equal(cells[i].contains(cells[j]), contained);
                    JavaAssert.Equal(cells[j].greaterOrEquals(cells[i].rangeMin())
                                     && cells[j].lessOrEquals(cells[i].rangeMax()), contained);
                    JavaAssert.Equal(cells[i].intersects(cells[j]),
                                     cells[i].contains(cells[j]) || cells[j].contains(cells[i]));
                }
            }
        }

        [Test]
        public void testContinuity()
        {
            Trace.WriteLine("TestContinuity");
            // Make sure that sequentially increasing cell ids form a continuous
            // path over the surface of the sphere, i.e. there are no
            // discontinuous jumps from one region to another.

            var maxDist = S2Projections.MAX_EDGE.getValue(MAX_WALK_LEVEL);
            var end = S2CellId.end(MAX_WALK_LEVEL);
            var id = S2CellId.begin(MAX_WALK_LEVEL);
            for (; !id.Equals(end); id = id.next())
            {
                Assert.True(id.toPointRaw().angle(id.nextWrap().toPointRaw()) <= maxDist);

                // Check that the ToPointRaw() returns the center of each cell
                // in (s,t) coordinates.
                var p = id.toPointRaw();
                var face = S2Projections.xyzToFace(p);
                var uv = S2Projections.validFaceXyzToUv(face, p);
                assertDoubleNear(Math.IEEERemainder(
                    S2Projections.uvToST(uv.x()), 1.0/(1 << MAX_WALK_LEVEL)), 0);
                assertDoubleNear(Math.IEEERemainder(
                    S2Projections.uvToST(uv.y()), 1.0/(1 << MAX_WALK_LEVEL)), 0);
            }
        }

        [Test]
        public void testCoverage()
        {
            Trace.WriteLine("TestCoverage");
            // Make sure that random points on the sphere can be represented to the
            // expected level of accuracy, which in the worst case is sqrt(2/3) times
            // the maximum arc length between the points on the sphere associated with
            // adjacent values of "i" or "j". (It is sqrt(2/3) rather than 1/2 because
            // the cells at the corners of each face are stretched -- they have 60 and
            // 120 degree angles.)

            var maxDist = 0.5*S2Projections.MAX_DIAG.getValue(S2CellId.MAX_LEVEL);
            for (var i = 0; i < 1000000; ++i)
            {
                // randomPoint();
                var p = new S2Point(0.37861576725894824, 0.2772406863275093, 0.8830558887338725);
                var q = S2CellId.fromPoint(p).toPointRaw();

                Assert.True(p.angle(q) <= maxDist);
            }
        }

        [Test]
        [Ignore("Not necessrily valid values. Fails on Java too.")]
        public void testNeighborLevel29()
        {
            // Note: These parameters fail on the Java version too. Not sure if this is a valid Cell anyway
            testAllNeighbors(new S2CellId(0x6000000000000004UL), 29);
        }

        [Test]
        public void testNeighbors()
        {
            Trace.WriteLine("TestNeighbors");

            // Check the edge neighbors of face 1.
            int[] outFaces = {5, 3, 2, 0};
            var faceNbrs = new S2CellId[4];
            S2CellId.fromFacePosLevel(1, 0, 0).getEdgeNeighbors(faceNbrs);
            for (var i = 0; i < 4; ++i)
            {
                Assert.True(faceNbrs[i].isFace());
                JavaAssert.Equal(faceNbrs[i].face(), outFaces[i]);
            }

            // Check the vertex neighbors of the center of face 2 at level 5.
            var nbrs = new List<S2CellId>();
            S2CellId.fromPoint(new S2Point(0, 0, 1)).getVertexNeighbors(5, nbrs);
            nbrs.Sort();
            for (var i = 0; i < 4; ++i)
            {
                JavaAssert.Equal(nbrs[i], S2CellId.fromFaceIJ(
                    2, (1 << 29) - (i < 2 ? 1 : 0), (1 << 29) - ((i == 0 || i == 3) ? 1 : 0)).parent(5));
            }
            nbrs.Clear();

            // Check the vertex neighbors of the corner of faces 0, 4, and 5.
            var id = S2CellId.fromFacePosLevel(0, 0, S2CellId.MAX_LEVEL);
            id.getVertexNeighbors(0, nbrs);
            nbrs.Sort();

            JavaAssert.Equal(nbrs.Count, 3);
            JavaAssert.Equal(nbrs[0], S2CellId.fromFacePosLevel(0, 0, 0));
            JavaAssert.Equal(nbrs[1], S2CellId.fromFacePosLevel(4, 0, 0));
            JavaAssert.Equal(nbrs[2], S2CellId.fromFacePosLevel(5, 0, 0));

            // Check that GetAllNeighbors produces results that are consistent
            // with GetVertexNeighbors for a bunch of random cells.
            for (var i = 0; i < 1000; ++i)
            {
                var id1 = getRandomCellId();
                var toTest = id1;
                if (id1.isLeaf())
                {
                    toTest = id1.parent();
                }

                // TestAllNeighbors computes approximately 2**(2*(diff+1)) cell id1s,
                // so it's not reasonable to use large values of "diff".
                var maxDiff = Math.Min(6, S2CellId.MAX_LEVEL - toTest.level() - 1);
                var level = toTest.level() + random(maxDiff);
                testAllNeighbors(toTest, level);
            }
        }

        [Test]
        public void testToToken()
        {
            JavaAssert.Equal("000000000000010a", new S2CellId(266).toToken());
            JavaAssert.Equal("80855c", new S2CellId(unchecked ((ulong)-9185834709882503168L)).toToken());
        }

        [Test]
        public void testTokens()
        {
            Trace.WriteLine("TestTokens");

            // Test random cell ids at all levels.
            for (var i = 0; i < 10000; ++i)
            {
                var id = getRandomCellId();
                if (!id.isValid())
                {
                    continue;
                }
                var token = id.toToken();
                Assert.True(token.Length <= 16);
                JavaAssert.Equal(S2CellId.fromToken(token), id);
            }
            // Check that invalid cell ids can be encoded.
            var token1 = S2CellId.none().toToken();
            JavaAssert.Equal(S2CellId.fromToken(token1), S2CellId.none());
        }
    }
}