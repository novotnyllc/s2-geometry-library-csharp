using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using C5;

namespace Google.Common.Geometry
{
    /**
 * An S2Polygon is an SI2Region object that represents a polygon. A polygon
 * consists of zero or more {@link S2Loop loops} representing "shells" and
 * "holes". All loops should be oriented CCW, i.e. the shell or hole is on the
 * left side of the loop. Loops may be specified in any order. A point is
 * defined to be inside the polygon if it is contained by an odd number of
 * loops.
 *
 *  Polygons have the following restrictions:
 *
 *  - Loops may not cross, i.e. the boundary of a loop may not intersect both
 * the interior and exterior of any other loop.
 *
 *  - Loops may not share edges, i.e. if a loop contains an edge AB, then no
 * other loop may contain AB or BA.
 *
 *  - No loop may cover more than half the area of the sphere. This ensures that
 * no loop properly contains the complement of any other loop, even if the loops
 * are from different polygons. (Loops that represent exact hemispheres are
 * allowed.)
 *
 *  Loops may share vertices, however no vertex may appear twice in a single
 * loop.
 *
 */

    public sealed class S2Polygon : IS2Region, IComparable<S2Polygon>
    {
        private readonly List<S2Loop> _loops;

        private S2LatLngRect _bound;
        private bool _hasHoles;
        private int _numVertices;

        /**
   * Creates an empty polygon that should be initialized by calling Init().
   */

        public S2Polygon()
        {
            _loops = new List<S2Loop>();
            _bound = S2LatLngRect.Empty;
            _hasHoles = false;
            _numVertices = 0;
        }

        /**
   * Convenience constructor that calls Init() with the given loops. Clears the
   * given list.
   */

        public S2Polygon(System.Collections.Generic.IList<S2Loop> loops)
        {
            _loops = new List<S2Loop>();
            _bound = S2LatLngRect.Empty;

            Init(loops);
        }

        /**
   * Copy constructor.
   */

        public S2Polygon(S2Loop loop)
        {
            _loops = new List<S2Loop>();
            _bound = loop.RectBound;
            _hasHoles = false;
            _numVertices = loop.NumVertices;

            _loops.Add(loop);
        }

        /**
   * Copy constructor.
   */

        public S2Polygon(S2Polygon src)
        {
            _loops = new List<S2Loop>();
            _bound = src.RectBound;
            _hasHoles = src._hasHoles;
            _numVertices = src._numVertices;

            for (var i = 0; i < src.NumLoops; ++i)
            {
                _loops.Add(new S2Loop(src.Loop(i)));
            }
        }

        public int NumLoops
        {
            get { return _loops.Count; }
        }

        public S2AreaCentroid AreaAndCentroid
        {
            get { return GetAreaCentroid(true); }
        }

        /**
   * Return the area of the polygon interior, i.e. the region on the left side
   * of an odd number of loops. The return value is between 0 and 4*Pi.
   */

        public double Area
        {
            get { return GetAreaCentroid(false).Area; }
        }

        /**
   * Return the true centroid of the polygon multiplied by the area of the
   * polygon (see s2.h for details on centroids). Note that the centroid may not
   * be contained by the polygon.
   */

        public S2Point? Centroid
        {
            get { return GetAreaCentroid(true).Centroid; }
        }

        public int NumVertices
        {
            get { return _numVertices; }
        }

        public bool IsNormalized
        {
            get
            {
                var vertices = new HashBag<S2Point>();


                S2Loop lastParent = null;
                for (var i = 0; i < NumLoops; ++i)
                {
                    var child = Loop(i);
                    if (child.Depth == 0)
                    {
                        continue;
                    }
                    var parent = Loop(GetParent(i));
                    if (parent != lastParent)
                    {
                        vertices.Clear();
                        for (var j = 0; j < parent.NumVertices; ++j)
                        {
                            vertices.Add(parent.Vertex(j));
                        }
                        lastParent = parent;
                    }
                    var count = 0;
                    for (var j = 0; j < child.NumVertices; ++j)
                    {
                        var item = child.Vertex(j);
                        if (vertices.Any(p => p.Equals(item)))
                        {
                            ++count;
                        }
                    }
                    if (count > 1)
                    {
                        return false;
                    }
                }
                return true;
            }
        }

        /**
   * Comparator (needed by Comparable interface). For two polygons to be
   * compared as equal: - the must have the same number of loops; - the loops
   * must be ordered in the same way (this is guaranteed by the total ordering
   * imposed by sortValueLoops). - loops must be logically equivalent (even if
   * ordered with a different starting point, e.g. ABCD and BCDA).
   */

        public int CompareTo(S2Polygon other)
        {
            // If number of loops differ, use that.
            if (NumLoops != other.NumLoops)
            {
                return NumLoops - other.NumLoops;
            }
            for (var i = 0; i < NumLoops; ++i)
            {
                var compare = _loops[i].CompareTo(other._loops[i]);
                if (compare != 0)
                {
                    return compare;
                }
            }
            return 0;
        }

        public S2Cap CapBound
        {
            get { return _bound.CapBound; }
        }


        /** Return a bounding latitude-longitude rectangle. */

        public S2LatLngRect RectBound
        {
            get { return _bound; }
        }

        /**
   * If this method returns true, the region completely contains the given cell.
   * Otherwise, either the region does not contain the cell or the containment
   * relationship could not be determined.
   */

        public bool Contains(S2Cell cell)
        {
            if (NumLoops == 1)
            {
                return Loop(0).Contains(cell);
            }
            var cellBound = cell.RectBound;
            if (!_bound.Contains(cellBound))
            {
                return false;
            }

            var cellLoop = new S2Loop(cell, cellBound);
            var cellPoly = new S2Polygon(cellLoop);
            return Contains(cellPoly);
        }

        /**
   * If this method returns false, the region does not intersect the given cell.
   * Otherwise, either region intersects the cell, or the intersection
   * relationship could not be determined.
   */

        public bool MayIntersect(S2Cell cell)
        {
            if (NumLoops == 1)
            {
                return Loop(0).MayIntersect(cell);
            }
            var cellBound = cell.RectBound;
            if (!_bound.Intersects(cellBound))
            {
                return false;
            }

            var cellLoop = new S2Loop(cell, cellBound);
            var cellPoly = new S2Polygon(cellLoop);
            return Intersects(cellPoly);
        }

        /**
   * Initialize a polygon by taking ownership of the given loops and clearing
   * the given list. This method figures out the loop nesting hierarchy and then
   * reorders the loops by following a preorder traversal. This implies that
   * each loop is immediately followed by its descendants in the nesting
   * hierarchy. (See also getParent and getLastDescendant.)
   */

        public void Init(System.Collections.Generic.IList<S2Loop> loops)
        {
            // assert isValid(loops);
            // assert (this.loops.isEmpty());

            //Dictionary<S2Loop, List<S2Loop>> loopMap =new Dictionary<S2Loop, List<S2Loop>>();
            // Note: We're using C5's HashDictionary because SCG's Dictionary<,> does not allow
            // NULL keys
            var loopMap = new HashDictionary<S2Loop, List<S2Loop>>();
            // Yes, a null key is valid. It is used here to refer to the root of the
            // loopMap
            loopMap[null] = new List<S2Loop>();

            foreach (var loop in loops)
            {
                InsertLoop(loop, null, loopMap);
                _numVertices += loop.NumVertices;
            }
            loops.Clear();

            // Sort all of the lists of loops; in this way we guarantee a total ordering
            // on loops in the polygon. Loops will be sorted by their natural ordering,
            // while also preserving the requirement that each loop is immediately
            // followed by its descendants in the nesting hierarchy.
            //
            // TODO(andriy): as per kirilll in CL 18750833 code review comments:
            // This should work for now, but I think it's possible to guarantee the
            // correct order inside insertLoop by searching for the correct position in
            // the children list before inserting.
            SortValueLoops(loopMap);

            // Reorder the loops in depth-first traversal order.
            // Starting at null == starting at the root
            InitLoop(null, -1, loopMap);

            // TODO(dbeaumont): Add tests or preconditions for these asserts (here and elesewhere).
            // forall i != j : containsChild(loop(i), loop(j), loopMap) == loop(i).containsNested(loop(j)));

            // Compute the bounding rectangle of the entire polygon.
            _hasHoles = false;
            _bound = S2LatLngRect.Empty;
            for (var i = 0; i < NumLoops; ++i)
            {
                if (Loop(i).Sign < 0)
                {
                    _hasHoles = true;
                }
                else
                {
                    _bound = _bound.Union(Loop(i).RectBound);
                }
            }
        }

        /**
   * Release ownership of the loops of this polygon by appending them to the
   * given list. Resets the polygon to be empty.
   */

        public void Release(System.Collections.Generic.IList<S2Loop> loops)
        {
            foreach (var item in _loops)
                loops.Add(item);

            _loops.Clear();
            _bound = S2LatLngRect.Empty;
            _hasHoles = false;
            _numVertices = 0;
        }

        /**
   * Return true if the given loops form a valid polygon. Assumes that that all
   * of the given loops have already been validated.
   */

        public static bool IsValidPolygon(IReadOnlyList<S2Loop> loops)
        {
            // If a loop contains an edge AB, then no other loop may contain AB or BA.
            // We only need this test if there are at least two loops, assuming that
            // each loop has already been validated.
            if (loops.Count > 1)
            {
                System.Collections.Generic.IDictionary<UndirectedEdge, LoopVertexIndexPair> edges = new Dictionary<UndirectedEdge, LoopVertexIndexPair>();
                for (var i = 0; i < loops.Count; ++i)
                {
                    var lp = loops[i];
                    for (var j = 0; j < lp.NumVertices; ++j)
                    {
                        var key = new UndirectedEdge(lp.Vertex(j), lp.Vertex(j + 1));
                        var value = new LoopVertexIndexPair(i, j);
                        if (edges.ContainsKey(key))
                        {
                            var other = edges[key];
                            Debug.WriteLine(
                                "Duplicate edge: loop " + i + ", edge " + j + " and loop " + other.LoopIndex
                                + ", edge " + other.VertexIndex);
                            return false;
                        }
                        else
                        {
                            edges[key] = value;
                        }
                    }
                }
            }

            // Verify that no loop covers more than half of the sphere, and that no
            // two loops cross.
            for (var i = 0; i < loops.Count; ++i)
            {
                if (!loops[i].IsNormalized)
                {
                    Debug.WriteLine("Loop " + i + " encloses more than half the sphere");
                    return false;
                }
                for (var j = i + 1; j < loops.Count; ++j)
                {
                    // This test not only checks for edge crossings, it also detects
                    // cases where the two boundaries cross at a shared vertex.
                    if (loops[i].ContainsOrCrosses(loops[j]) < 0)
                    {
                        Debug.WriteLine("Loop " + i + " crosses loop " + j);
                        return false;
                    }
                }
            }
            return true;
        }

        public S2Loop Loop(int k)
        {
            return _loops[k];
        }

        /**
   * Return the index of the parent of loop k, or -1 if it has no parent.
   */

        public int GetParent(int k)
        {
            var depth = Loop(k).Depth;
            if (depth == 0)
            {
                return -1; // Optimization.
            }
            while (--k >= 0 && Loop(k).Depth >= depth)
            {
                // spin
            }
            return k;
        }

        /**
   * Return the index of the last loop that is contained within loop k. Returns
   * num_loops() - 1 if k < 0. Note that loops are indexed according to a
   * preorder traversal of the nesting hierarchy, so the immediate children of
   * loop k can be found by iterating over loops (k+1)..getLastDescendant(k) and
   * selecting those whose depth is equal to (loop(k).depth() + 1).
   */

        public int GetLastDescendant(int k)
        {
            if (k < 0)
            {
                return NumLoops - 1;
            }
            var depth = Loop(k).Depth;
            while (++k < NumLoops && Loop(k).Depth > depth)
            {
                // spin
            }
            return k - 1;
        }

        private S2AreaCentroid GetAreaCentroid(bool doCentroid)
        {
            double areaSum = 0;
            var centroidSum = new S2Point(0, 0, 0);
            for (var i = 0; i < NumLoops; ++i)
            {
                var areaCentroid = doCentroid ? (S2AreaCentroid?)Loop(i).AreaAndCentroid : null;
                var loopArea = doCentroid ? areaCentroid.Value.Area : Loop(i).Area;

                var loopSign = Loop(i).Sign;
                areaSum += loopSign*loopArea;
                if (doCentroid)
                {
                    var currentCentroid = areaCentroid.Value.Centroid.Value;
                    centroidSum =
                        new S2Point(centroidSum.X + loopSign*currentCentroid.X,
                                    centroidSum.Y + loopSign*currentCentroid.Y,
                                    centroidSum.Z + loopSign*currentCentroid.Z);
                }
            }

            return new S2AreaCentroid(areaSum, doCentroid ? (S2Point?)centroidSum : null);
        }

        /**
   * Return the area of the polygon interior, i.e. the region on the left side
   * of an odd number of loops (this value return value is between 0 and 4*Pi)
   * and the true centroid of the polygon multiplied by the area of the polygon
   * (see s2.h for details on centroids). Note that the centroid may not be
   * contained by the polygon.
   */

        /**
   * Returns the shortest distance from a point P to this polygon, given as the
   * angle formed between P, the origin and the nearest point on the polygon to
   * P. This angle in radians is equivalent to the arclength along the unit
   * sphere.
   *
   * If the point is contained inside the polygon, the distance returned is 0.
   */

        public S1Angle GetDistance(S2Point p)
        {
            if (Contains(p))
            {
                return S1Angle.FromRadians(0);
            }

            // The furthest point from p on the sphere is its antipode, which is an
            // angle of PI radians. This is an upper bound on the angle.
            var minDistance = S1Angle.FromRadians(Math.PI);
            for (var i = 0; i < NumLoops; i++)
            {
                minDistance = S1Angle.Min(minDistance, Loop(i).GetDistance(p));
            }

            return minDistance;
        }


        /**
   * Return true if this polygon contains the given other polygon, i.e. if
   * polygon A contains all points contained by polygon B.
   */

        public bool Contains(S2Polygon b)
        {
            // If both polygons have one loop, use the more efficient S2Loop method.
            // Note that S2Loop.contains does its own bounding rectangle check.
            if (NumLoops == 1 && b.NumLoops == 1)
            {
                return Loop(0).Contains(b.Loop(0));
            }

            // Otherwise if neither polygon has holes, we can still use the more
            // efficient S2Loop::Contains method (rather than ContainsOrCrosses),
            // but it's worthwhile to do our own bounds check first.
            if (!_bound.Contains(b.RectBound))
            {
                // If the union of the bounding boxes spans the full longitude range,
                // it is still possible that polygon A contains B. (This is only
                // possible if at least one polygon has multiple shells.)
                if (!_bound.Lng.Union(b.RectBound.Lng).IsFull)
                {
                    return false;
                }
            }
            if (!_hasHoles && !b._hasHoles)
            {
                for (var j = 0; j < b.NumLoops; ++j)
                {
                    if (!AnyLoopContains(b.Loop(j)))
                    {
                        return false;
                    }
                }
                return true;
            }

            // This could be implemented more efficiently for polygons with lots of
            // holes by keeping a copy of the LoopMap computed during initialization.
            // However, in practice most polygons are one loop, and multiloop polygons
            // tend to consist of many shells rather than holes. In any case, the real
            // way to get more efficiency is to implement a sub-quadratic algorithm
            // such as building a trapezoidal map.

            // Every shell of B must be contained by an odd number of loops of A,
            // and every hole of A must be contained by an even number of loops of B.
            return ContainsAllShells(b) && b.ExcludesAllHoles(this);
        }

        /**
   * Return true if this polygon intersects the given other polygon, i.e. if
   * there is a point that is contained by both polygons.
   */

        public bool Intersects(S2Polygon b)
        {
            // A.intersects(B) if and only if !complement(A).contains(B). However,
            // implementing a complement() operation is trickier than it sounds,
            // and in any case it's more efficient to test for intersection directly.

            // If both polygons have one loop, use the more efficient S2Loop method.
            // Note that S2Loop.intersects does its own bounding rectangle check.
            if (NumLoops == 1 && b.NumLoops == 1)
            {
                return Loop(0).Intersects(b.Loop(0));
            }

            // Otherwise if neither polygon has holes, we can still use the more
            // efficient S2Loop.intersects method. The polygons intersect if and
            // only if some pair of loop regions intersect.
            if (!_bound.Intersects(b.RectBound))
            {
                return false;
            }
            if (!_hasHoles && !b._hasHoles)
            {
                for (var i = 0; i < NumLoops; ++i)
                {
                    for (var j = 0; j < b.NumLoops; ++j)
                    {
                        if (Loop(i).Intersects(b.Loop(j)))
                        {
                            return true;
                        }
                    }
                }
                return false;
            }

            // Otherwise if any shell of B is contained by an odd number of loops of A,
            // or any shell of A is contained by an odd number of loops of B, there is
            // an intersection.
            return IntersectsAnyShell(b) || b.IntersectsAnyShell(this);
        }

        /**
   *  Indexing structure to efficiently clipEdge() of a polygon. This is an
   * abstract class because we need to use if for both polygons (for
   * initToIntersection() and friends) and for sets of lists of points (for
   * initToSimplified()).
   *
   *  Usage -- in your subclass, create an array of vertex counts for each loop
   * in the loop sequence and pass it to this constructor. Overwrite
   * edgeFromTo(), calling decodeIndex() and use the resulting two indices to
   * access your accessing vertices.
   */

        private static void AddIntersection(S2Point a0,
                                            S2Point a1,
                                            S2Point b0,
                                            S2Point b1,
                                            bool addSharedEdges,
                                            int crossing, System.Collections.Generic.ICollection<ParametrizedS2Point> intersections)
        {
            if (crossing > 0)
            {
                // There is a proper edge crossing.
                var x = S2EdgeUtil.GetIntersection(a0, a1, b0, b1);
                var t = S2EdgeUtil.GetDistanceFraction(x, a0, a1);
                intersections.Add(new ParametrizedS2Point(t, x));
            }
            else if (S2EdgeUtil.VertexCrossing(a0, a1, b0, b1))
            {
                // There is a crossing at one of the vertices. The basic rule is simple:
                // if a0 equals one of the "b" vertices, the crossing occurs at t=0;
                // otherwise, it occurs at t=1.
                //
                // This has the effect that when two symmetric edges are encountered (an
                // edge an its reverse), neither one is included in the output. When two
                // duplicate edges are encountered, both are included in the output. The
                // "addSharedEdges" flag allows one of these two copies to be removed by
                // changing its intersection parameter from 0 to 1.
                double t = (a0 == b0 || a0 == b1) ? 0 : 1;
                if (!addSharedEdges && a1 == b1)
                {
                    t = 1;
                }
                intersections.Add(new ParametrizedS2Point(t, t == 0 ? a0 : a1));
            }
        }

        /**
   * Find all points where the polygon B intersects the edge (a0,a1), and add
   * the corresponding parameter values (in the range [0,1]) to "intersections".
   */

        private static void ClipEdge(S2Point a0, S2Point a1, S2LoopSequenceIndex bIndex,
                                     bool addSharedEdges, System.Collections.Generic.ICollection<ParametrizedS2Point> intersections)
        {
            var it = new S2EdgeIndex.DataEdgeIterator(bIndex);
            it.GetCandidates(a0, a1);
            var crosser = new EdgeCrosser(a0, a1, a0);
            //S2Point from;
            var to = default (S2Point);

            foreach (var index in it)
            {
                var previousTo = to;
                var fromTo = bIndex.EdgeFromTo(index);
                var from = fromTo.Start;
                to = fromTo.End;
                if (previousTo != from)
                {
                    crosser.RestartAt(from);
                }
                var crossing = crosser.RobustCrossing(to);
                if (crossing < 0)
                {
                    continue;
                }
                AddIntersection(a0, a1, from, to, addSharedEdges, crossing, intersections);
            }
        }

        /**
   * Clip the boundary of A to the interior of B, and add the resulting edges to
   * "builder". Shells are directed CCW and holes are directed clockwise, unless
   * "reverseA" or "reverseB" is true in which case these directions in the
   * corresponding polygon are reversed. If "invertB" is true, the boundary of A
   * is clipped to the exterior rather than the interior of B. If
   * "adSharedEdges" is true, then the output will include any edges that are
   * shared between A and B (both edges must be in the same direction after any
   * edge reversals are taken into account).
   */

        private static void ClipBoundary(S2Polygon a,
                                         bool reverseA,
                                         S2Polygon b,
                                         bool reverseB,
                                         bool invertB,
                                         bool addSharedEdges,
                                         S2PolygonBuilder builder)
        {
            var bIndex = new S2PolygonIndex(b, reverseB);
            bIndex.PredictAdditionalCalls(a.NumVertices);

            var intersections = new List<ParametrizedS2Point>();
            foreach (var aLoop in a._loops)
            {
                var n = aLoop.NumVertices;
                var dir = (aLoop.IsHole ^ reverseA) ? -1 : 1;
                var inside = b.Contains(aLoop.Vertex(0)) ^ invertB;
                for (var j = (dir > 0) ? 0 : n; n > 0; --n, j += dir)
                {
                    var a0 = aLoop.Vertex(j);
                    var a1 = aLoop.Vertex(j + dir);
                    intersections.Clear();
                    ClipEdge(a0, a1, bIndex, addSharedEdges, intersections);

                    if (inside)
                    {
                        intersections.Add(new ParametrizedS2Point(0.0, a0));
                    }
                    inside = ((intersections.Count & 0x1) == 0x1);
                    // assert ((b.contains(a1) ^ invertB) == inside);
                    if (inside)
                    {
                        intersections.Add(new ParametrizedS2Point(1.0, a1));
                    }

                    // Remove duplicates and produce a list of unique intersections.
                    intersections.Sort();
                    for (int size = intersections.Count, i = 1; i < size; i += 2)
                    {
                        builder.AddEdge(intersections[i - 1].Point, intersections[i].Point);
                    }
                }
            }
        }

        /**
   * Returns total number of vertices in all loops.
   */

        /**
   * Initialize this polygon to the intersection, union, or difference (A - B)
   * of the given two polygons. The "vertexMergeRadius" determines how close two
   * vertices must be to be merged together and how close a vertex must be to an
   * edge in order to be spliced into it (see S2PolygonBuilder for details). By
   * default, the merge radius is just large enough to compensate for errors
   * that occur when computing intersection points between edges
   * (S2EdgeUtil.DEFAULT_INTERSECTION_TOLERANCE).
   *
   *  If you are going to convert the resulting polygon to a lower-precision
   * format, it is necessary to increase the merge radius in order to get a
   * valid result after rounding (i.e. no duplicate vertices, etc). For example,
   * if you are going to convert them to geostore.PolygonProto format, then
   * S1Angle.e7(1) is a good value for "vertex_merge_radius".
   */

        public void InitToIntersection(S2Polygon a, S2Polygon b)
        {
            InitToIntersectionSloppy(a, b, S2EdgeUtil.DefaultIntersectionTolerance);
        }

        public void InitToIntersectionSloppy(
            S2Polygon a, S2Polygon b, S1Angle vertexMergeRadius)
        {
            Preconditions.CheckState(NumLoops == 0);
            if (!a._bound.Intersects(b._bound))
            {
                return;
            }

            // We want the boundary of A clipped to the interior of B,
            // plus the boundary of B clipped to the interior of A,
            // plus one copy of any directed edges that are in both boundaries.

            var options = S2PolygonBuilderOptions.DirectedXor;
            options.MergeDistance = vertexMergeRadius;
            var builder = new S2PolygonBuilder(options);
            ClipBoundary(a, false, b, false, false, true, builder);
            ClipBoundary(b, false, a, false, false, false, builder);
            if (!builder.AssemblePolygon(this, null))
            {
                // TODO (andriy): do something more meaningful here.
                Debug.WriteLine("Bad directed edges");
            }
        }

        public void InitToUnion(S2Polygon a, S2Polygon b)
        {
            InitToUnionSloppy(a, b, S2EdgeUtil.DefaultIntersectionTolerance);
        }

        public void InitToUnionSloppy(S2Polygon a, S2Polygon b, S1Angle vertexMergeRadius)
        {
            Preconditions.CheckState(NumLoops == 0);

            // We want the boundary of A clipped to the exterior of B,
            // plus the boundary of B clipped to the exterior of A,
            // plus one copy of any directed edges that are in both boundaries.

            var options = S2PolygonBuilderOptions.DirectedXor;
            options.MergeDistance = vertexMergeRadius;
            var builder = new S2PolygonBuilder(options);
            ClipBoundary(a, false, b, false, true, true, builder);
            ClipBoundary(b, false, a, false, true, false, builder);
            if (!builder.AssemblePolygon(this, null))
            {
                // TODO(andriy): do something more meaningful here.
                Debug.WriteLine("Bad directed edges");
            }
        }

        /**
   * Return a polygon which is the union of the given polygons. Note: clears the
   * List!
   */

        public static S2Polygon DestructiveUnion(System.Collections.Generic.ICollection<S2Polygon> polygons)
        {
            return DestructiveUnionSloppy(polygons, S2EdgeUtil.DefaultIntersectionTolerance);
        }

        /**
   * Return a polygon which is the union of the given polygons; combines
   * vertices that form edges that are almost identical, as defined by
   * vertexMergeRadius. Note: clears the List!
   */

        public static S2Polygon DestructiveUnionSloppy(System.Collections.Generic.ICollection<S2Polygon> polygons, S1Angle vertexMergeRadius)
        {
            // Effectively create a priority queue of polygons in order of number of
            // vertices. Repeatedly union the two smallest polygons and add the result
            // to the queue until we have a single polygon to return.

            // map: # of vertices -> polygon
            //TreeMultimap<Integer, S2Polygon> queue = TreeMultimap.create();

            var queue = new MultiMap<int, S2Polygon>();

            foreach (var polygon in polygons)
            {
                queue.Add(polygon.NumVertices, polygon);
            }
            polygons.Clear();

            // Java uses a live-view that maps to the underlying structure
            //Set<Map.Entry<Integer, S2Polygon>> queueSet = queue.entries();

            var enumer = queue.SortedValues;

            while (queue.CountIsAtLeast(2))
            {
                // Pop two simplest polygons from queue.
//      queueSet = queue.entries();
                //Iterator<Map.Entry<Integer, S2Polygon>> smallestIter = queueSet.iterator();

                var smallestTwo = enumer.Take(2).ToList();

                //Map.Entry<Integer, S2Polygon> smallest = smallestIter.next();
                //int aSize = smallest.getKey().intValue();
                //S2Polygon aPolygon = smallest.getValue();
                //smallestIter.remove();

                //smallest = smallestIter.next();
                //int bSize = smallest.getKey().intValue();
                //S2Polygon bPolygon = smallest.getValue();
                //smallestIter.remove();

                foreach (var item in smallestTwo)
                    queue.Remove(item);


                // Union and add result back to queue.
                var unionPolygon = new S2Polygon();
                unionPolygon.InitToUnionSloppy(smallestTwo[0].Value, smallestTwo[1].Value, vertexMergeRadius);
                var unionSize = smallestTwo[0].Key + smallestTwo[1].Key;
                queue.Add(unionSize, unionPolygon);
                // We assume that the number of vertices in the union polygon is the
                // sum of the number of vertices in the original polygons, which is not
                // always true, but will almost always be a decent approximation, and
                // faster than recomputing.
            }

            if (queue.Count == 0)
            {
                return new S2Polygon();
            }
            else
            {
                //return queue.get(queue.asMap().firstKey()).first();
                return queue.SortedValues.First().Value;
            }
        }

        /**
   * Return true if two polygons have the same boundary except for vertex
   * perturbations. Both polygons must have loops with the same cyclic vertex
   * order and the same nesting hierarchy, but the vertex locations are allowed
   * to differ by up to "max_error". Note: This method mostly useful only for
   * testing purposes.
   */

        internal bool BoundaryApproxEquals(S2Polygon b, double maxError)
        {
            if (NumLoops != b.NumLoops)
            {
                Debug.WriteLine(
                    "!= loops: " + NumLoops + " vs. " + b.NumLoops);
                return false;
            }

            // For now, we assume that there is at most one candidate match for each
            // loop. (So far this method is just used for testing.)
            for (var i = 0; i < NumLoops; ++i)
            {
                var aLoop = Loop(i);
                var success = false;
                for (var j = 0; j < NumLoops; ++j)
                {
                    var bLoop = b.Loop(j);
                    if (bLoop.Depth == aLoop.Depth && bLoop.BoundaryApproxEquals(aLoop, maxError))
                    {
                        success = true;
                        break;
                    }
                }
                if (!success)
                {
                    return false;
                }
            }
            return true;
        }

        // S2Region interface (see S2Region.java for details):

        /** Return a bounding spherical cap. */

        /**
   * The point 'p' does not need to be normalized.
   */

        public bool Contains(S2Point p)
        {
            if (NumLoops == 1)
            {
                return Loop(0).Contains(p); // Optimization.
            }
            if (!_bound.Contains(p))
            {
                return false;
            }
            var inside = false;
            for (var i = 0; i < NumLoops; ++i)
            {
                inside ^= Loop(i).Contains(p);
                if (inside && !_hasHoles)
                {
                    break; // Shells are disjoint.
                }
            }
            return inside;
        }

        // For each map entry, sorts the value list.
        private static void SortValueLoops(C5.IDictionary<S2Loop, List<S2Loop>> loopMap)
        {
            foreach (var key in loopMap.Keys)
            {
                loopMap[key].Sort();
            }
        }

        private static void InsertLoop(S2Loop newLoop, S2Loop parent, C5.IDictionary<S2Loop, List<S2Loop>> loopMap)
        {
            List<S2Loop> children = null;
            if (loopMap.Contains(parent))
                children = loopMap[parent];

            if (children == null)
            {
                children = new List<S2Loop>();
                loopMap[parent] = children;
            }

            foreach (var child in children)
            {
                if (child.ContainsNested(newLoop))
                {
                    InsertLoop(newLoop, child, loopMap);
                    return;
                }
            }

            // No loop may contain the complement of another loop. (Handling this case
            // is significantly more complicated.)
            // assert (parent == null || !newLoop.containsNested(parent));

            // Some of the children of the parent loop may now be children of
            // the new loop.
            List<S2Loop> newChildren = null;
            if (loopMap.Contains(newLoop))
                newChildren = loopMap[newLoop];
            for (var i = 0; i < children.Count;)
            {
                var child = children[i];
                if (newLoop.ContainsNested(child))
                {
                    if (newChildren == null)
                    {
                        newChildren = new List<S2Loop>();
                        loopMap[newLoop] = newChildren;
                    }
                    newChildren.Add(child);
                    children.RemoveAt(i);
                }
                else
                {
                    ++i;
                }
            }
            children.Add(newLoop);
        }

        private void InitLoop(S2Loop loop, int depth, C5.IDictionary<S2Loop, List<S2Loop>> loopMap)
        {
            if (loop != null)
            {
                loop.Depth = depth;
                _loops.Add(loop);
            }
            List<S2Loop> children = null;
            if (loopMap.Contains(loop))
                children = loopMap[loop];
            if (children != null)
            {
                foreach (var child in children)
                {
                    InitLoop(child, depth + 1, loopMap);
                }
            }
        }

        private int ContainsOrCrosses(S2Loop b)
        {
            var inside = false;
            for (var i = 0; i < NumLoops; ++i)
            {
                var result = Loop(i).ContainsOrCrosses(b);
                if (result < 0)
                {
                    return -1; // The loop boundaries intersect.
                }
                if (result > 0)
                {
                    inside ^= true;
                }
            }
            return inside ? 1 : 0; // True if loop B is contained by the polygon.
        }

        /** Return true if any loop contains the given loop. */

        private bool AnyLoopContains(S2Loop b)
        {
            for (var i = 0; i < NumLoops; ++i)
            {
                if (Loop(i).Contains(b))
                {
                    return true;
                }
            }
            return false;
        }

        /** Return true if this polygon (A) contains all the shells of B. */

        private bool ContainsAllShells(S2Polygon b)
        {
            for (var j = 0; j < b.NumLoops; ++j)
            {
                if (b.Loop(j).Sign < 0)
                {
                    continue;
                }
                if (ContainsOrCrosses(b.Loop(j)) <= 0)
                {
                    // Shell of B is not contained by A, or the boundaries intersect.
                    return false;
                }
            }
            return true;
        }

        /**
   * Return true if this polygon (A) excludes (i.e. does not intersect) all
   * holes of B.
   */

        private bool ExcludesAllHoles(S2Polygon b)
        {
            for (var j = 0; j < b.NumLoops; ++j)
            {
                if (b.Loop(j).Sign > 0)
                {
                    continue;
                }
                if (ContainsOrCrosses(b.Loop(j)) != 0)
                {
                    // Hole of B is contained by A, or the boundaries intersect.
                    return false;
                }
            }
            return true;
        }

        /** Return true if this polygon (A) intersects any shell of B. */

        private bool IntersectsAnyShell(S2Polygon b)
        {
            for (var j = 0; j < b.NumLoops; ++j)
            {
                if (b.Loop(j).Sign < 0)
                {
                    continue;
                }
                if (ContainsOrCrosses(b.Loop(j)) != 0)
                {
                    // Shell of B is contained by A, or the boundaries intersect.
                    return true;
                }
            }
            return false;
        }

        /**
   * A human readable representation of the polygon
   */

        public override String ToString()
        {
            var sb = new StringBuilder();
            sb.Append("Polygon: (").Append(NumLoops).Append(") loops:\n");
            for (var i = 0; i < NumLoops; ++i)
            {
                var s2Loop = Loop(i);
                sb.Append("loop <\n");
                for (var v = 0; v < s2Loop.NumVertices; ++v)
                {
                    var s2Point = s2Loop.Vertex(v);
                    sb.Append(s2Point.ToDegreesString());
                    sb.Append("\n"); // end of vertex
                }
                sb.Append(">\n"); // end of loop
            }
            return sb.ToString();
        }

        private struct LoopVertexIndexPair
        {
            private readonly int _loopIndex;
            private readonly int _vertexIndex;

            public LoopVertexIndexPair(int loopIndex, int vertexIndex)
            {
                _loopIndex = loopIndex;
                _vertexIndex = vertexIndex;
            }

            public int LoopIndex
            {
                get { return _loopIndex; }
            }

            public int VertexIndex
            {
                get { return _vertexIndex; }
            }
        }

        /**
   * An S2Point that also has a parameter associated with it, which corresponds
   * to a time-like order on the points.
   */

        private struct ParametrizedS2Point : IComparable<ParametrizedS2Point>
        {
            private readonly S2Point _point;
            private readonly double _time;

            public ParametrizedS2Point(double time, S2Point point)
            {
                _time = time;
                _point = point;
            }

            public double Time
            {
                get { return _time; }
            }

            public S2Point Point
            {
                get { return _point; }
            }

            public int CompareTo(ParametrizedS2Point o)
            {
                var compareTime = _time.CompareTo(o._time);
                if (compareTime != 0)
                {
                    return compareTime;
                }
                return _point.CompareTo(o._point);
            }
        }

        private abstract class S2LoopSequenceIndex : S2EdgeIndex
        {
            /** Map from the unidimensional edge index to the loop this edge belongs to. */
            private readonly int[] _indexToLoop;

            /**
     * Reverse of indexToLoop: maps a loop index to the unidimensional index
     * of the first edge in the loop.
     */
            private readonly int[] _loopToFirstIndex;

            /**
     * Must be called by each subclass with the array of vertices per loop. The
     * length of the array is the number of loops, and the <code>i</code>
     * <sup>th</sup> loop's vertex count is in the <code>i</code>
     * <sup>th</sup> index of the array.
     */

            protected S2LoopSequenceIndex(System.Collections.Generic.IList<int> numVertices)
            {
                var totalEdges = 0;
                foreach (var edges in numVertices)
                {
                    totalEdges += edges;
                }
                _indexToLoop = new int[totalEdges];
                _loopToFirstIndex = new int[numVertices.Count];

                totalEdges = 0;
                for (var j = 0; j < numVertices.Count; j++)
                {
                    _loopToFirstIndex[j] = totalEdges;
                    for (var i = 0; i < numVertices[j]; i++)
                    {
                        _indexToLoop[totalEdges] = j;
                        totalEdges++;
                    }
                }
            }

            protected override int NumEdges
            {
                get { return _indexToLoop.Length; }
            }

            protected LoopVertexIndexPair DecodeIndex(int index)
            {
                var loopIndex = _indexToLoop[index];
                var vertexInLoop = index - _loopToFirstIndex[loopIndex];
                return new LoopVertexIndexPair(loopIndex, vertexInLoop);
            }

            // It is faster to return both vertices at once. It makes a difference
            // for small polygons.
            public abstract S2Edge EdgeFromTo(int index);


            protected override S2Point EdgeFrom(int index)
            {
                var fromTo = EdgeFromTo(index);
                var from = fromTo.Start;
                return from;
            }

            protected override S2Point EdgeTo(int index)
            {
                var fromTo = EdgeFromTo(index);
                var to = fromTo.End;
                return to;
            }
        }

        // Indexing structure for an S2Polygon.
        private sealed class S2PolygonIndex : S2LoopSequenceIndex
        {
            private readonly S2Polygon _poly;
            private readonly bool _reverse;

            public S2PolygonIndex(S2Polygon poly, bool reverse) : base(GetVertices(poly))
            {
                _poly = poly;
                _reverse = reverse;
            }

            private static int[] GetVertices(S2Polygon poly)
            {
                var vertices = new int[poly.NumLoops];
                for (var i = 0; i < vertices.Length; i++)
                {
                    vertices[i] = poly.Loop(i).NumVertices;
                }
                return vertices;
            }

            public override S2Edge EdgeFromTo(int index)
            {
                var indices = DecodeIndex(index);
                var loopIndex = indices.LoopIndex;
                var vertexInLoop = indices.VertexIndex;
                var loop = _poly.Loop(loopIndex);
                int fromIndex;
                int toIndex;
                if (loop.IsHole ^ _reverse)
                {
                    fromIndex = loop.NumVertices - 1 - vertexInLoop;
                    toIndex = 2*loop.NumVertices - 2 - vertexInLoop;
                }
                else
                {
                    fromIndex = vertexInLoop;
                    toIndex = vertexInLoop + 1;
                }
                var from = loop.Vertex(fromIndex);
                var to = loop.Vertex(toIndex);
                return new S2Edge(from, to);
            }
        }

        private struct UndirectedEdge : IEquatable<UndirectedEdge>
        {
            private readonly S2Point _a;
            private readonly S2Point _b;

            public UndirectedEdge(S2Point start, S2Point end)
            {
                _a = start;
                _b = end;
            }

            public S2Point Start
            {
                get { return _a; }
            }

            public S2Point End
            {
                get { return _b; }
            }

            public bool Equals(UndirectedEdge other)
            {
                return ((Start.Equals(other.Start) && End.Equals(other.End))
                        || (Start.Equals(other.End) && End.Equals(other.Start)));
            }

            public override bool Equals(object obj)
            {
                if (ReferenceEquals(null, obj)) return false;
                return obj is UndirectedEdge && Equals((UndirectedEdge)obj);
            }

            public override int GetHashCode()
            {
                unchecked
                {
                    return (_a.GetHashCode()*397) ^ _b.GetHashCode();
                }
            }

            public static bool operator ==(UndirectedEdge left, UndirectedEdge right)
            {
                return Equals(left, right);
            }

            public static bool operator !=(UndirectedEdge left, UndirectedEdge right)
            {
                return !Equals(left, right);
            }

            // Note: An UndirectedEdge and an S2Edge can never be considered equal (in
            // terms of the equals() method) and hence they re not be related types.
            // If you need to convert between the types then separate conversion
            // methods should be introduced.

            public override String ToString()
            {
                return String.Format("Edge: ({0} <-> {1})\n   or [{2} <-> {3}]",
                                     _a.ToDegreesString(), _b.ToDegreesString(), _a, _b);
            }
        }
    }
}