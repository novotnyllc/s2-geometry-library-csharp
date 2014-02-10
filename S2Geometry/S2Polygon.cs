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
 * An S2Polygon is an S2Region object that represents a polygon. A polygon
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

    public class S2Polygon : IS2Region, IComparable<S2Polygon>
    {
        private readonly List<S2Loop> loops;

        private S2LatLngRect bound;
        private bool hasHoles;
        private int numVertices;

        /**
   * Creates an empty polygon that should be initialized by calling Init().
   */

        public S2Polygon()
        {
            loops = new List<S2Loop>();
            bound = S2LatLngRect.Empty;
            hasHoles = false;
            numVertices = 0;
        }

        /**
   * Convenience constructor that calls Init() with the given loops. Clears the
   * given list.
   */

        public S2Polygon(List<S2Loop> loops)
        {
            this.loops = new List<S2Loop>();
            bound = S2LatLngRect.Empty;

            init(loops);
        }

        /**
   * Copy constructor.
   */

        public S2Polygon(S2Loop loop)
        {
            loops = new List<S2Loop>();
            bound = loop.RectBound;
            hasHoles = false;
            numVertices = loop.NumVertices;

            loops.Add(loop);
        }

        /**
   * Copy constructor.
   */

        public S2Polygon(S2Polygon src)
        {
            loops = new List<S2Loop>();
            bound = src.RectBound;
            hasHoles = src.hasHoles;
            numVertices = src.numVertices;

            for (var i = 0; i < src.numLoops(); ++i)
            {
                loops.Add(new S2Loop(src.loop(i)));
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
            if (numLoops() != other.numLoops())
            {
                return numLoops() - other.numLoops();
            }
            for (var i = 0; i < numLoops(); ++i)
            {
                var compare = loops[i].CompareTo(other.loops[i]);
                if (compare != 0)
                {
                    return compare;
                }
            }
            return 0;
        }

        public S2Cap CapBound
        {
            get { return bound.CapBound; }
        }


        /** Return a bounding latitude-longitude rectangle. */

        public S2LatLngRect RectBound
        {
            get { return bound; }
        }

        /**
   * If this method returns true, the region completely contains the given cell.
   * Otherwise, either the region does not contain the cell or the containment
   * relationship could not be determined.
   */

        public bool Contains(S2Cell cell)
        {
            if (numLoops() == 1)
            {
                return loop(0).Contains(cell);
            }
            var cellBound = cell.RectBound;
            if (!bound.Contains(cellBound))
            {
                return false;
            }

            var cellLoop = new S2Loop(cell, cellBound);
            var cellPoly = new S2Polygon(cellLoop);
            return contains(cellPoly);
        }

        /**
   * If this method returns false, the region does not intersect the given cell.
   * Otherwise, either region intersects the cell, or the intersection
   * relationship could not be determined.
   */

        public bool MayIntersect(S2Cell cell)
        {
            if (numLoops() == 1)
            {
                return loop(0).MayIntersect(cell);
            }
            var cellBound = cell.RectBound;
            if (!bound.Intersects(cellBound))
            {
                return false;
            }

            var cellLoop = new S2Loop(cell, cellBound);
            var cellPoly = new S2Polygon(cellLoop);
            return intersects(cellPoly);
        }

        /**
   * Initialize a polygon by taking ownership of the given loops and clearing
   * the given list. This method figures out the loop nesting hierarchy and then
   * reorders the loops by following a preorder traversal. This implies that
   * each loop is immediately followed by its descendants in the nesting
   * hierarchy. (See also getParent and getLastDescendant.)
   */

        public void init(List<S2Loop> loops)
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
                insertLoop(loop, null, loopMap);
                numVertices += loop.NumVertices;
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
            sortValueLoops(loopMap);

            // Reorder the loops in depth-first traversal order.
            // Starting at null == starting at the root
            initLoop(null, -1, loopMap);

            // TODO(dbeaumont): Add tests or preconditions for these asserts (here and elesewhere).
            // forall i != j : containsChild(loop(i), loop(j), loopMap) == loop(i).containsNested(loop(j)));

            // Compute the bounding rectangle of the entire polygon.
            hasHoles = false;
            bound = S2LatLngRect.Empty;
            for (var i = 0; i < numLoops(); ++i)
            {
                if (loop(i).Sign < 0)
                {
                    hasHoles = true;
                }
                else
                {
                    bound = bound.Union(loop(i).RectBound);
                }
            }
        }

        /**
   * Release ownership of the loops of this polygon by appending them to the
   * given list. Resets the polygon to be empty.
   */

        public void release(List<S2Loop> loops)
        {
            loops.AddRange(this.loops);
            this.loops.Clear();
            bound = S2LatLngRect.Empty;
            hasHoles = false;
            numVertices = 0;
        }

        /**
   * Return true if the given loops form a valid polygon. Assumes that that all
   * of the given loops have already been validated.
   */

        public static bool isValid(List<S2Loop> loops)
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
                                "Duplicate edge: loop " + i + ", edge " + j + " and loop " + other.getLoopIndex()
                                + ", edge " + other.getVertexIndex());
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

        public int numLoops()
        {
            return loops.Count;
        }

        public S2Loop loop(int k)
        {
            return loops[k];
        }

        /**
   * Return the index of the parent of loop k, or -1 if it has no parent.
   */

        public int getParent(int k)
        {
            var depth = loop(k).Depth;
            if (depth == 0)
            {
                return -1; // Optimization.
            }
            while (--k >= 0 && loop(k).Depth >= depth)
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

        public int getLastDescendant(int k)
        {
            if (k < 0)
            {
                return numLoops() - 1;
            }
            var depth = loop(k).Depth;
            while (++k < numLoops() && loop(k).Depth > depth)
            {
                // spin
            }
            return k - 1;
        }

        private S2AreaCentroid getAreaCentroid(bool doCentroid)
        {
            double areaSum = 0;
            var centroidSum = new S2Point(0, 0, 0);
            for (var i = 0; i < numLoops(); ++i)
            {
                var areaCentroid = doCentroid ? (S2AreaCentroid?) loop(i).AreaAndCentroid : null;
                var loopArea = doCentroid ? areaCentroid.Value.Area : loop(i).Area;

                var loopSign = loop(i).Sign;
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

        public S2AreaCentroid getAreaAndCentroid()
        {
            return getAreaCentroid(true);
        }

        /**
   * Return the area of the polygon interior, i.e. the region on the left side
   * of an odd number of loops. The return value is between 0 and 4*Pi.
   */

        public double getArea()
        {
            return getAreaCentroid(false).Area;
        }

        /**
   * Return the true centroid of the polygon multiplied by the area of the
   * polygon (see s2.h for details on centroids). Note that the centroid may not
   * be contained by the polygon.
   */

        public S2Point? getCentroid()
        {
            return getAreaCentroid(true).Centroid;
        }

        /**
   * Returns the shortest distance from a point P to this polygon, given as the
   * angle formed between P, the origin and the nearest point on the polygon to
   * P. This angle in radians is equivalent to the arclength along the unit
   * sphere.
   *
   * If the point is contained inside the polygon, the distance returned is 0.
   */

        public S1Angle getDistance(S2Point p)
        {
            if (contains(p))
            {
                return S1Angle.FromRadians(0);
            }

            // The furthest point from p on the sphere is its antipode, which is an
            // angle of PI radians. This is an upper bound on the angle.
            var minDistance = S1Angle.FromRadians(Math.PI);
            for (var i = 0; i < numLoops(); i++)
            {
                minDistance = S1Angle.Min(minDistance, loop(i).GetDistance(p));
            }

            return minDistance;
        }


        /**
   * Return true if this polygon contains the given other polygon, i.e. if
   * polygon A contains all points contained by polygon B.
   */

        public bool contains(S2Polygon b)
        {
            // If both polygons have one loop, use the more efficient S2Loop method.
            // Note that S2Loop.contains does its own bounding rectangle check.
            if (numLoops() == 1 && b.numLoops() == 1)
            {
                return loop(0).Contains(b.loop(0));
            }

            // Otherwise if neither polygon has holes, we can still use the more
            // efficient S2Loop::Contains method (rather than ContainsOrCrosses),
            // but it's worthwhile to do our own bounds check first.
            if (!bound.Contains(b.RectBound))
            {
                // If the union of the bounding boxes spans the full longitude range,
                // it is still possible that polygon A contains B. (This is only
                // possible if at least one polygon has multiple shells.)
                if (!bound.Lng.Union(b.RectBound.Lng).IsFull)
                {
                    return false;
                }
            }
            if (!hasHoles && !b.hasHoles)
            {
                for (var j = 0; j < b.numLoops(); ++j)
                {
                    if (!anyLoopContains(b.loop(j)))
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
            return containsAllShells(b) && b.excludesAllHoles(this);
        }

        /**
   * Return true if this polygon intersects the given other polygon, i.e. if
   * there is a point that is contained by both polygons.
   */

        public bool intersects(S2Polygon b)
        {
            // A.intersects(B) if and only if !complement(A).contains(B). However,
            // implementing a complement() operation is trickier than it sounds,
            // and in any case it's more efficient to test for intersection directly.

            // If both polygons have one loop, use the more efficient S2Loop method.
            // Note that S2Loop.intersects does its own bounding rectangle check.
            if (numLoops() == 1 && b.numLoops() == 1)
            {
                return loop(0).Intersects(b.loop(0));
            }

            // Otherwise if neither polygon has holes, we can still use the more
            // efficient S2Loop.intersects method. The polygons intersect if and
            // only if some pair of loop regions intersect.
            if (!bound.Intersects(b.RectBound))
            {
                return false;
            }
            if (!hasHoles && !b.hasHoles)
            {
                for (var i = 0; i < numLoops(); ++i)
                {
                    for (var j = 0; j < b.numLoops(); ++j)
                    {
                        if (loop(i).Intersects(b.loop(j)))
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
            return intersectsAnyShell(b) || b.intersectsAnyShell(this);
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

        private static void addIntersection(S2Point a0,
                                            S2Point a1,
                                            S2Point b0,
                                            S2Point b1,
                                            bool addSharedEdges,
                                            int crossing,
                                            List<ParametrizedS2Point> intersections)
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

        private static void clipEdge(S2Point a0, S2Point a1, S2LoopSequenceIndex bIndex,
                                     bool addSharedEdges, List<ParametrizedS2Point> intersections)
        {
            var it = new S2EdgeIndex.DataEdgeIterator(bIndex);
            it.GetCandidates(a0, a1);
            var crosser = new EdgeCrosser(a0, a1, a0);
            //S2Point from;
            var to = default (S2Point);
            
            foreach (var index in it)
            {
                var previousTo = to;
                var fromTo = bIndex.edgeFromTo(index);
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
                addIntersection(a0, a1, from, to, addSharedEdges, crossing, intersections);
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

        private static void clipBoundary(S2Polygon a,
                                         bool reverseA,
                                         S2Polygon b,
                                         bool reverseB,
                                         bool invertB,
                                         bool addSharedEdges,
                                         S2PolygonBuilder builder)
        {
            var bIndex = new S2PolygonIndex(b, reverseB);
            bIndex.PredictAdditionalCalls(a.getNumVertices());

            var intersections = new List<ParametrizedS2Point>();
            foreach (var aLoop in a.loops)
            {
                var n = aLoop.NumVertices;
                var dir = (aLoop.IsHole ^ reverseA) ? -1 : 1;
                var inside = b.contains(aLoop.Vertex(0)) ^ invertB;
                for (var j = (dir > 0) ? 0 : n; n > 0; --n, j += dir)
                {
                    var a0 = aLoop.Vertex(j);
                    var a1 = aLoop.Vertex(j + dir);
                    intersections.Clear();
                    clipEdge(a0, a1, bIndex, addSharedEdges, intersections);

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
                        builder.addEdge(intersections[i - 1].getPoint(), intersections[i].getPoint());
                    }
                }
            }
        }

        /**
   * Returns total number of vertices in all loops.
   */

        public int getNumVertices()
        {
            return numVertices;
        }

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

        public void initToIntersection(S2Polygon a, S2Polygon b)
        {
            initToIntersectionSloppy(a, b, S2EdgeUtil.DefaultIntersectionTolerance);
        }

        public void initToIntersectionSloppy(
            S2Polygon a, S2Polygon b, S1Angle vertexMergeRadius)
        {
            Preconditions.CheckState(numLoops() == 0);
            if (!a.bound.Intersects(b.bound))
            {
                return;
            }

            // We want the boundary of A clipped to the interior of B,
            // plus the boundary of B clipped to the interior of A,
            // plus one copy of any directed edges that are in both boundaries.

            var options = S2PolygonBuilder.Options.DIRECTED_XOR;
            options.setMergeDistance(vertexMergeRadius);
            var builder = new S2PolygonBuilder(options);
            clipBoundary(a, false, b, false, false, true, builder);
            clipBoundary(b, false, a, false, false, false, builder);
            if (!builder.assemblePolygon(this, null))
            {
                // TODO (andriy): do something more meaningful here.
                Debug.WriteLine("Bad directed edges");
            }
        }

        public void initToUnion(S2Polygon a, S2Polygon b)
        {
            initToUnionSloppy(a, b, S2EdgeUtil.DefaultIntersectionTolerance);
        }

        public void initToUnionSloppy(S2Polygon a, S2Polygon b, S1Angle vertexMergeRadius)
        {
            Preconditions.CheckState(numLoops() == 0);

            // We want the boundary of A clipped to the exterior of B,
            // plus the boundary of B clipped to the exterior of A,
            // plus one copy of any directed edges that are in both boundaries.

            var options = S2PolygonBuilder.Options.DIRECTED_XOR;
            options.setMergeDistance(vertexMergeRadius);
            var builder = new S2PolygonBuilder(options);
            clipBoundary(a, false, b, false, true, true, builder);
            clipBoundary(b, false, a, false, true, false, builder);
            if (!builder.assemblePolygon(this, null))
            {
                // TODO(andriy): do something more meaningful here.
                Debug.WriteLine("Bad directed edges");
            }
        }

        /**
   * Return a polygon which is the union of the given polygons. Note: clears the
   * List!
   */

        public static S2Polygon destructiveUnion(List<S2Polygon> polygons)
        {
            return destructiveUnionSloppy(polygons, S2EdgeUtil.DefaultIntersectionTolerance);
        }

        /**
   * Return a polygon which is the union of the given polygons; combines
   * vertices that form edges that are almost identical, as defined by
   * vertexMergeRadius. Note: clears the List!
   */

        public static S2Polygon destructiveUnionSloppy(
            List<S2Polygon> polygons, S1Angle vertexMergeRadius)
        {
            // Effectively create a priority queue of polygons in order of number of
            // vertices. Repeatedly union the two smallest polygons and add the result
            // to the queue until we have a single polygon to return.

            // map: # of vertices -> polygon
            //TreeMultimap<Integer, S2Polygon> queue = TreeMultimap.create();

            var queue = new MultiMap<int, S2Polygon>();

            foreach (var polygon in polygons)
            {
                queue.Add(polygon.getNumVertices(), polygon);
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
                unionPolygon.initToUnionSloppy(smallestTwo[0].Value, smallestTwo[1].Value, vertexMergeRadius);
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

        public bool isNormalized()
        {
            var vertices = new HashBag<S2Point>();


            S2Loop lastParent = null;
            for (var i = 0; i < numLoops(); ++i)
            {
                var child = loop(i);
                if (child.Depth == 0)
                {
                    continue;
                }
                var parent = loop(getParent(i));
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

        /**
   * Return true if two polygons have the same boundary except for vertex
   * perturbations. Both polygons must have loops with the same cyclic vertex
   * order and the same nesting hierarchy, but the vertex locations are allowed
   * to differ by up to "max_error". Note: This method mostly useful only for
   * testing purposes.
   */

        internal bool boundaryApproxEquals(S2Polygon b, double maxError)
        {
            if (numLoops() != b.numLoops())
            {
                Debug.WriteLine(
                    "!= loops: " + numLoops() + " vs. " + b.numLoops());
                return false;
            }

            // For now, we assume that there is at most one candidate match for each
            // loop. (So far this method is just used for testing.)
            for (var i = 0; i < numLoops(); ++i)
            {
                var aLoop = loop(i);
                var success = false;
                for (var j = 0; j < numLoops(); ++j)
                {
                    var bLoop = b.loop(j);
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

        public bool contains(S2Point p)
        {
            if (numLoops() == 1)
            {
                return loop(0).Contains(p); // Optimization.
            }
            if (!bound.Contains(p))
            {
                return false;
            }
            var inside = false;
            for (var i = 0; i < numLoops(); ++i)
            {
                inside ^= loop(i).Contains(p);
                if (inside && !hasHoles)
                {
                    break; // Shells are disjoint.
                }
            }
            return inside;
        }

        // For each map entry, sorts the value list.
        private static void sortValueLoops(C5.IDictionary<S2Loop, List<S2Loop>> loopMap)
        {
            foreach (var key in loopMap.Keys)
            {
                loopMap[key].Sort();
            }
        }

        private static void insertLoop(S2Loop newLoop, S2Loop parent, C5.IDictionary<S2Loop, List<S2Loop>> loopMap)
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
                    insertLoop(newLoop, child, loopMap);
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

        private void initLoop(S2Loop loop, int depth, C5.IDictionary<S2Loop, List<S2Loop>> loopMap)
        {
            if (loop != null)
            {
                loop.Depth = depth;
                loops.Add(loop);
            }
            List<S2Loop> children = null;
            if (loopMap.Contains(loop))
                children = loopMap[loop];
            if (children != null)
            {
                foreach (var child in children)
                {
                    initLoop(child, depth + 1, loopMap);
                }
            }
        }

        private int containsOrCrosses(S2Loop b)
        {
            var inside = false;
            for (var i = 0; i < numLoops(); ++i)
            {
                var result = loop(i).ContainsOrCrosses(b);
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

        private bool anyLoopContains(S2Loop b)
        {
            for (var i = 0; i < numLoops(); ++i)
            {
                if (loop(i).Contains(b))
                {
                    return true;
                }
            }
            return false;
        }

        /** Return true if this polygon (A) contains all the shells of B. */

        private bool containsAllShells(S2Polygon b)
        {
            for (var j = 0; j < b.numLoops(); ++j)
            {
                if (b.loop(j).Sign < 0)
                {
                    continue;
                }
                if (containsOrCrosses(b.loop(j)) <= 0)
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

        private bool excludesAllHoles(S2Polygon b)
        {
            for (var j = 0; j < b.numLoops(); ++j)
            {
                if (b.loop(j).Sign > 0)
                {
                    continue;
                }
                if (containsOrCrosses(b.loop(j)) != 0)
                {
                    // Hole of B is contained by A, or the boundaries intersect.
                    return false;
                }
            }
            return true;
        }

        /** Return true if this polygon (A) intersects any shell of B. */

        private bool intersectsAnyShell(S2Polygon b)
        {
            for (var j = 0; j < b.numLoops(); ++j)
            {
                if (b.loop(j).Sign < 0)
                {
                    continue;
                }
                if (containsOrCrosses(b.loop(j)) != 0)
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
            sb.Append("Polygon: (").Append(numLoops()).Append(") loops:\n");
            for (var i = 0; i < numLoops(); ++i)
            {
                var s2Loop = loop(i);
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

        private sealed class LoopVertexIndexPair
        {
            private readonly int loopIndex;
            private readonly int vertexIndex;

            public LoopVertexIndexPair(int loopIndex, int vertexIndex)
            {
                this.loopIndex = loopIndex;
                this.vertexIndex = vertexIndex;
            }

            public int getLoopIndex()
            {
                return loopIndex;
            }

            public int getVertexIndex()
            {
                return vertexIndex;
            }
        }

        /**
   * An S2Point that also has a parameter associated with it, which corresponds
   * to a time-like order on the points.
   */

        private sealed class ParametrizedS2Point : IComparable<ParametrizedS2Point>
        {
            private readonly S2Point point;
            private readonly double time;

            public ParametrizedS2Point(double time, S2Point point)
            {
                this.time = time;
                this.point = point;
            }

            public int CompareTo(ParametrizedS2Point o)
            {
                var compareTime = time.CompareTo(o.time);
                if (compareTime != 0)
                {
                    return compareTime;
                }
                return point.CompareTo(o.point);
            }

            public double getTime()
            {
                return time;
            }

            public S2Point getPoint()
            {
                return point;
            }
        }

        private abstract class S2LoopSequenceIndex : S2EdgeIndex
        {
            /** Map from the unidimensional edge index to the loop this edge belongs to. */
            private readonly int[] indexToLoop;

            /**
     * Reverse of indexToLoop: maps a loop index to the unidimensional index
     * of the first edge in the loop.
     */
            private readonly int[] loopToFirstIndex;

            /**
     * Must be called by each subclass with the array of vertices per loop. The
     * length of the array is the number of loops, and the <code>i</code>
     * <sup>th</sup> loop's vertex count is in the <code>i</code>
     * <sup>th</sup> index of the array.
     */

            public S2LoopSequenceIndex(int[] numVertices)
            {
                var totalEdges = 0;
                foreach (var edges in numVertices)
                {
                    totalEdges += edges;
                }
                indexToLoop = new int[totalEdges];
                loopToFirstIndex = new int[numVertices.Length];

                totalEdges = 0;
                for (var j = 0; j < numVertices.Length; j++)
                {
                    loopToFirstIndex[j] = totalEdges;
                    for (var i = 0; i < numVertices[j]; i++)
                    {
                        indexToLoop[totalEdges] = j;
                        totalEdges++;
                    }
                }
            }

            public LoopVertexIndexPair decodeIndex(int index)
            {
                var loopIndex = indexToLoop[index];
                var vertexInLoop = index - loopToFirstIndex[loopIndex];
                return new LoopVertexIndexPair(loopIndex, vertexInLoop);
            }

            // It is faster to return both vertices at once. It makes a difference
            // for small polygons.
            public abstract S2Edge edgeFromTo(int index);


            protected override int NumEdges
            {
                get { return indexToLoop.Length; }
            }

            protected override S2Point EdgeFrom(int index)
            {
                var fromTo = edgeFromTo(index);
                var from = fromTo.Start;
                return from;
            }

            protected override S2Point EdgeTo(int index)
            {
                var fromTo = edgeFromTo(index);
                var to = fromTo.End;
                return to;
            }
        }

        // Indexing structure for an S2Polygon.
        private sealed class S2PolygonIndex : S2LoopSequenceIndex
        {
            private readonly S2Polygon poly;
            private readonly bool reverse;

            public S2PolygonIndex(S2Polygon poly, bool reverse) : base(getVertices(poly))
            {
                this.poly = poly;
                this.reverse = reverse;
            }

            private static int[] getVertices(S2Polygon poly)
            {
                var vertices = new int[poly.numLoops()];
                for (var i = 0; i < vertices.Length; i++)
                {
                    vertices[i] = poly.loop(i).NumVertices;
                }
                return vertices;
            }

            public override S2Edge edgeFromTo(int index)
            {
                var indices = decodeIndex(index);
                var loopIndex = indices.getLoopIndex();
                var vertexInLoop = indices.getVertexIndex();
                var loop = poly.loop(loopIndex);
                int fromIndex;
                int toIndex;
                if (loop.IsHole ^ reverse)
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

        private sealed class UndirectedEdge : IEquatable<UndirectedEdge>
        {
            private readonly S2Point a;
            private readonly S2Point b;

            public UndirectedEdge(S2Point start, S2Point end)
            {
                a = start;
                b = end;
            }

            public bool Equals(UndirectedEdge other)
            {
                if (ReferenceEquals(null, other)) return false;
                if (ReferenceEquals(this, other)) return true;
                return ((getStart().Equals(other.getStart()) && getEnd().Equals(other.getEnd()))
                        || (getStart().Equals(other.getEnd()) && getEnd().Equals(other.getStart())));
            }

            public override bool Equals(object obj)
            {
                if (ReferenceEquals(null, obj)) return false;
                if (ReferenceEquals(this, obj)) return true;
                return obj is UndirectedEdge && Equals((UndirectedEdge)obj);
            }

            public override int GetHashCode()
            {
                unchecked
                {
                    return (a.GetHashCode()*397) ^ b.GetHashCode();
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

            public S2Point getStart()
            {
                return a;
            }

            public S2Point getEnd()
            {
                return b;
            }

            public override String ToString()
            {
                return String.Format("Edge: ({0} <-> {1})\n   or [{2} <-> {3}]",
                                     a.ToDegreesString(), b.ToDegreesString(), a, b);
            }
        }
    }
}