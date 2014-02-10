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
            var id = S2CellId.FromLatLng(S2LatLng.FromDegrees(latDegrees, lngDegrees));
            Trace.WriteLine(Convert.ToString(unchecked ((long)id.Id), 16));
            return id;
        }

        public void testInverses()
        {
            Trace.WriteLine("TestInverses");
            // Check the conversion of random leaf cells to S2LatLngs and back.
            for (var i = 0; i < 200000; ++i)
            {
                var id = getRandomCellId(S2CellId.MaxLevel);
                Assert.True(id.IsLeaf && id.Level == S2CellId.MaxLevel);
                var center = id.ToLatLng();
                JavaAssert.Equal(S2CellId.FromLatLng(center).Id, id.Id);
            }
        }

        private const int kMaxExpandLevel = 3;

        private void expandCell(
            S2CellId parent, List<S2CellId> cells, IDictionary<S2CellId, S2CellId> parentMap)
        {
            cells.Add(parent);
            if (parent.Level == kMaxExpandLevel)
            {
                return;
            }
            var i = 0;
            var j = 0;
            int? orientation = 0;
            var face = parent.ToFaceIjOrientation(ref i, ref j, ref orientation);
            JavaAssert.Equal(face, parent.Face);

            var pos = 0;
            for (var child = parent.ChildBegin; !child.Equals(parent.ChildEnd);
                 child = child.Next)
            {
                // Do some basic checks on the children
                JavaAssert.Equal(child.Level, parent.Level + 1);
                Assert.True(!child.IsLeaf);
                int? childOrientation = 0;
                JavaAssert.Equal(child.ToFaceIjOrientation(ref i, ref j, ref childOrientation), face);
                JavaAssert.Equal(
                    childOrientation.Value, orientation.Value ^ S2.PosToOrientation(pos));

                parentMap.Add(child, parent);
                expandCell(child, cells, parentMap);
                ++pos;
            }
        }

        private const int MAX_WALK_LEVEL = 8;

        public void testAllNeighbors(S2CellId id, int level)
        {
            Assert.True(level >= id.Level && level < S2CellId.MaxLevel);

            // We compute GetAllNeighbors, and then add in all the children of "id"
            // at the given level. We then compare this against the result of finding
            // all the vertex neighbors of all the vertices of children of "id" at the
            // given level. These should give the same result.
            var all = new List<S2CellId>();
            var expected = new List<S2CellId>();
            id.GetAllNeighbors(level, all);
            var end = id.ChildEndForLevel(level + 1);
            for (var c = id.ChildBeginForLevel(level + 1); !c.Equals(end); c = c.Next)
            {
                all.Add(c.Parent);
                c.GetVertexNeighbors(level, expected);
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
            id = S2CellId.FromFacePosLevel(3, 0x12345678, S2CellId.MaxLevel - 4);
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

            JavaAssert.Equal(S2CellId.Begin(S2CellId.MaxLevel).PreviousWithWrap,
                             S2CellId.FromFacePosLevel(5, ~0UL >> S2CellId.FaceBits, S2CellId.MaxLevel));

            JavaAssert.Equal(S2CellId.End(4).Previous.NextWithWrap, S2CellId.Begin(4));
            JavaAssert.Equal(S2CellId.End(S2CellId.MaxLevel).Previous.NextWithWrap,
                             S2CellId.FromFacePosLevel(0, 0, S2CellId.MaxLevel));

            // Check that cells are represented by the position of their center
            // along the Hilbert curve.
            JavaAssert.Equal(id.RangeMin.Id + id.RangeMax.Id, 2*id.Id);
        }

        [Test]
        public void testContainment()
        {
            Trace.WriteLine("TestContainment");
            IDictionary<S2CellId, S2CellId> parentMap = new Dictionary<S2CellId, S2CellId>();
            var cells = new List<S2CellId>();
            for (var face = 0; face < 6; ++face)
            {
                expandCell(S2CellId.FromFacePosLevel(face, 0, 0), cells, parentMap);
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
                    JavaAssert.Equal(cells[i].Contains(cells[j]), contained);
                    JavaAssert.Equal(cells[j] >= cells[i].RangeMin
                                     && cells[j] <= cells[i].RangeMax, contained);
                    JavaAssert.Equal(cells[i].Intersects(cells[j]),
                                     cells[i].Contains(cells[j]) || cells[j].Contains(cells[i]));
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

            var maxDist = S2Projections.MaxEdge.GetValue(MAX_WALK_LEVEL);
            var end = S2CellId.End(MAX_WALK_LEVEL);
            var id = S2CellId.Begin(MAX_WALK_LEVEL);
            for (; !id.Equals(end); id = id.Next)
            {
                Assert.True(id.ToPointRaw().Angle(id.NextWithWrap.ToPointRaw()) <= maxDist);

                // Check that the ToPointRaw() returns the center of each cell
                // in (s,t) coordinates.
                var p = id.ToPointRaw();
                var face = S2Projections.XyzToFace(p);
                var uv = S2Projections.ValidFaceXyzToUv(face, p);
                assertDoubleNear(Math.IEEERemainder(
                    S2Projections.UvToSt(uv.X), 1.0/(1 << MAX_WALK_LEVEL)), 0);
                assertDoubleNear(Math.IEEERemainder(
                    S2Projections.UvToSt(uv.Y), 1.0/(1 << MAX_WALK_LEVEL)), 0);
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

            var maxDist = 0.5*S2Projections.MaxDiag.GetValue(S2CellId.MaxLevel);
            for (var i = 0; i < 1000000; ++i)
            {
                // randomPoint();
                var p = new S2Point(0.37861576725894824, 0.2772406863275093, 0.8830558887338725);
                var q = S2CellId.FromPoint(p).ToPointRaw();

                Assert.True(p.Angle(q) <= maxDist);
            }
        }

        //[Test]
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
            
            var faceNbrs = S2CellId.FromFacePosLevel(1, 0, 0).GetEdgeNeighbors();
            for (var i = 0; i < 4; ++i)
            {
                Assert.True(faceNbrs[i].IsFace);
                JavaAssert.Equal(faceNbrs[i].Face, outFaces[i]);
            }

            // Check the vertex neighbors of the center of face 2 at level 5.
            var nbrs = new List<S2CellId>();
            S2CellId.FromPoint(new S2Point(0, 0, 1)).GetVertexNeighbors(5, nbrs);
            nbrs.Sort();
            for (var i = 0; i < 4; ++i)
            {
                JavaAssert.Equal(nbrs[i], S2CellId.FromFaceIj(
                    2, (1 << 29) - (i < 2 ? 1 : 0), (1 << 29) - ((i == 0 || i == 3) ? 1 : 0)).ParentForLevel(5));
            }
            nbrs.Clear();

            // Check the vertex neighbors of the corner of faces 0, 4, and 5.
            var id = S2CellId.FromFacePosLevel(0, 0, S2CellId.MaxLevel);
            id.GetVertexNeighbors(0, nbrs);
            nbrs.Sort();

            JavaAssert.Equal(nbrs.Count, 3);
            JavaAssert.Equal(nbrs[0], S2CellId.FromFacePosLevel(0, 0, 0));
            JavaAssert.Equal(nbrs[1], S2CellId.FromFacePosLevel(4, 0, 0));
            JavaAssert.Equal(nbrs[2], S2CellId.FromFacePosLevel(5, 0, 0));

            // Check that GetAllNeighbors produces results that are consistent
            // with GetVertexNeighbors for a bunch of random cells.
            for (var i = 0; i < 1000; ++i)
            {
                var id1 = getRandomCellId();
                var toTest = id1;
                if (id1.IsLeaf)
                {
                    toTest = id1.Parent;
                }

                // TestAllNeighbors computes approximately 2**(2*(diff+1)) cell id1s,
                // so it's not reasonable to use large values of "diff".
                var maxDiff = Math.Min(6, S2CellId.MaxLevel - toTest.Level - 1);
                var level = toTest.Level + random(maxDiff);
                testAllNeighbors(toTest, level);
            }
        }

        [Test]
        public void testToToken()
        {
            JavaAssert.Equal("000000000000010a", new S2CellId(266).ToToken());
            JavaAssert.Equal("80855c", new S2CellId(unchecked ((ulong)-9185834709882503168L)).ToToken());
        }

        [Test]
        public void testTokens()
        {
            Trace.WriteLine("TestTokens");

            // Test random cell ids at all levels.
            for (var i = 0; i < 10000; ++i)
            {
                var id = getRandomCellId();
                if (!id.IsValid)
                {
                    continue;
                }
                var token = id.ToToken();
                Assert.True(token.Length <= 16);
                JavaAssert.Equal(S2CellId.FromToken(token), id);
            }
            // Check that invalid cell ids can be encoded.
            var token1 = S2CellId.None.ToToken();
            JavaAssert.Equal(S2CellId.FromToken(token1), S2CellId.None);
        }
    }
}