using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Google.Common.Geometry
{
    /**
     *
     * An S2Loop represents a simple spherical polygon. It consists of a single
     * chain of vertices where the first vertex is implicitly connected to the last.
     * All loops are defined to have a CCW orientation, i.e. the interior of the
     * polygon is on the left side of the edges. This implies that a clockwise loop
     * enclosing a small area is interpreted to be a CCW loop enclosing a very large
     * area.
     *
     *  Loops are not allowed to have any duplicate vertices (whether adjacent or
     * not), and non-adjacent edges are not allowed to intersect. Loops must have at
     * least 3 vertices. Although these restrictions are not enforced in optimized
     * code, you may get unexpected results if they are violated.
     *
     *  Point containment is defined such that if the sphere is subdivided into
     * faces (loops), every point is contained by exactly one face. This implies
     * that loops do not necessarily contain all (or any) of their vertices An
     * S2LatLngRect represents a latitude-longitude rectangle. It is capable of
     * representing the empty and full rectangles as well as single points.
     *
     */

    public sealed class S2Loop : IS2Region, IComparable<S2Loop>
    {
        /**
   * Max angle that intersections can be off by and yet still be considered
   * colinear.
   */
        private const double MaxIntersectionError = 1e-15;

        /**
   * Edge index used for performance-critical operations. For example,
   * contains() can determine whether a point is inside a loop in nearly
   * constant time, whereas without an edge index it is forced to compare the
   * query point against every edge in the loop.
   */
        private readonly int _numVertices;
        private readonly S2Point[] _vertices;

        /**
   * The index (into "vertices") of the vertex that comes first in the total
   * ordering of all vertices in this loop.
   */

        private S2LatLngRect _bound;
        private int _depth;
        private int _firstLogicalVertex;
        private S2EdgeIndex _index;
        private bool _originInside;
        private Dictionary<S2Point, int> _vertexToIndex;

        /**
   * Initialize a loop connecting the given vertices. The last vertex is
   * implicitly connected to the first. All points should be unit length. Loops
   * must have at least 3 vertices.
   *
   * @param vertices
   */

        public S2Loop(IEnumerable<S2Point> vertices)
        {
            _vertices = vertices.ToArray();
            _numVertices = _vertices.Length;
            _bound = S2LatLngRect.Full;
            _depth = 0;

            // if (debugMode) {
            //  assert (isValid(vertices, DEFAULT_MAX_ADJACENT));
            // }


            // initOrigin() must be called before InitBound() because the latter
            // function expects Contains() to work properly.
            InitOrigin();
            InitBound();
            InitFirstLogicalVertex();
        }

        /**
   * Initialize a loop corresponding to the given cell.
   */

        public S2Loop(S2Cell cell) : this(cell, cell.RectBound)
        {
        }

        /**
   * Like the constructor above, but assumes that the cell's bounding rectangle
   * has been precomputed.
   *
   * @param cell
   * @param bound
   */

        public S2Loop(S2Cell cell, S2LatLngRect bound)
        {
            _bound = bound;
            _numVertices = 4;
            _vertices = new S2Point[_numVertices];
            _vertexToIndex = null;
            _index = null;
            _depth = 0;
            for (var i = 0; i < 4; ++i)
            {
                _vertices[i] = cell.GetVertex(i);
            }
            InitOrigin();
            InitFirstLogicalVertex();
        }

        /**
   * Copy constructor.
   */

        public S2Loop(S2Loop src)
        {
            _numVertices = src._numVertices;
            _vertices = (S2Point[])src._vertices.Clone();
            _vertexToIndex = src._vertexToIndex;
            _index = src._index;
            _firstLogicalVertex = src._firstLogicalVertex;
            _bound = src.RectBound;
            _originInside = src._originInside;
            _depth = src._depth;
        }

        public int Depth
        {
            get { return _depth; }
            set { _depth = value; }
        }

        /**
   * Return true if this loop represents a hole in its containing polygon.
   */

        public bool IsHole
        {
            get { return (_depth & 1) != 0; }
        }

        /**
   * The sign of a loop is -1 if the loop represents a hole in its containing
   * polygon, and +1 otherwise.
   */

        public int Sign
        {
            get { return IsHole ? -1 : 1; }
        }

        public int NumVertices
        {
            get { return _numVertices; }
        }

        public bool IsNormalized
        {
            get
            {
                // We allow a bit of error so that exact hemispheres are
                // considered normalized.
                return Area <= 2*S2.Pi + 1e-14;
            }
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
   * polygon (see {@link S2} for details on centroids). Note that the centroid
   * may not be contained by the polygon.
   */

        public S2Point? Centroid
        {
            get { return GetAreaCentroid(true).Centroid; }
        }

        public bool IsValid
        {
            get
            {
                if (_numVertices < 3)
                {
                    Debug.WriteLine("Degenerate loop");
                    return false;
                }

                // All vertices must be unit length.
                for (var i = 0; i < _numVertices; ++i)
                {
                    if (!S2.IsUnitLength(Vertex(i)))
                    {
                        Debug.WriteLine("Vertex " + i + " is not unit length");
                        return false;
                    }
                }

                // Loops are not allowed to have any duplicate vertices.
                var vmap = new Dictionary<S2Point, int>();
                for (var i = 0; i < _numVertices; ++i)
                {
                    var key = Vertex(i);
                    var contains = vmap.ContainsKey(key);
                    if (contains)
                    {
                        var prevIndex = vmap[key];
                        Debug.WriteLine("Duplicate vertices: " + prevIndex + " and " + i);
                    }
                    // update always
                    vmap[key] = i;
                    if (contains)
                        return false;
                }


                // Non-adjacent edges are not allowed to intersect.
                // var crosses = false;
                var it = GetEdgeIterator(_numVertices);
                for (var a1 = 0; a1 < _numVertices; a1++)
                {
                    var a2 = (a1 + 1)%_numVertices;
                    var crosser = new EdgeCrosser(Vertex(a1), Vertex(a2), Vertex(0));
                    var previousIndex = -2;
                    it.GetCandidates(Vertex(a1), Vertex(a2));
                    foreach (var b1 in it) // it.GetCandidates(vertex(a1), vertex(a2)); it.HasNext; it.Next())
                    {
                        //var b1 = it.Index;
                        var b2 = (b1 + 1)%_numVertices;
                        // If either 'a' index equals either 'b' index, then these two edges
                        // share a vertex. If a1==b1 then it must be the case that a2==b2, e.g.
                        // the two edges are the same. In that case, we skip the test, since we
                        // don't want to test an edge against itself. If a1==b2 or b1==a2 then
                        // we have one edge ending at the start of the other, or in other words,
                        // the edges share a vertex -- and in S2 space, where edges are always
                        // great circle segments on a sphere, edges can only intersect at most
                        // once, so we don't need to do further checks in that case either.
                        if (a1 != b2 && a2 != b1 && a1 != b1)
                        {
                            // WORKAROUND(shakusa, ericv): S2.robustCCW() currently
                            // requires arbitrary-precision arithmetic to be truly robust. That
                            // means it can give the wrong answers in cases where we are trying
                            // to determine edge intersections. The workaround is to ignore
                            // intersections between edge pairs where all four points are
                            // nearly colinear.
                            var abc = S2.Angle(Vertex(a1), Vertex(a2), Vertex(b1));
                            var abcNearlyLinear = S2.ApproxEquals(abc, 0D, MaxIntersectionError) ||
                                                  S2.ApproxEquals(abc, S2.Pi, MaxIntersectionError);
                            var abd = S2.Angle(Vertex(a1), Vertex(a2), Vertex(b2));
                            var abdNearlyLinear = S2.ApproxEquals(abd, 0D, MaxIntersectionError) ||
                                                  S2.ApproxEquals(abd, S2.Pi, MaxIntersectionError);
                            if (abcNearlyLinear && abdNearlyLinear)
                            {
                                continue;
                            }

                            if (previousIndex != b1)
                            {
                                crosser.RestartAt(Vertex(b1));
                            }

                            // Beware, this may return the loop is valid if there is a
                            // "vertex crossing".
                            // TODO(user): Fix that.
                            var crosses = crosser.RobustCrossing(Vertex(b2)) > 0;
                            previousIndex = b2;
                            if (crosses)
                            {
                                Debug.WriteLine("Edges " + a1 + " and " + b1 + " cross");
                                Debug.WriteLine("Edge locations in degrees: " + "{0}-{1} and {2}-{3}",
                                                new S2LatLng(Vertex(a1)).ToStringDegrees(),
                                                new S2LatLng(Vertex(a2)).ToStringDegrees(),
                                                new S2LatLng(Vertex(b1)).ToStringDegrees(),
                                                new S2LatLng(Vertex(b2)).ToStringDegrees());
                                return false;
                            }
                        }
                    }
                }

                return true;
            }
        }

        public int CompareTo(S2Loop other)
        {
            if (NumVertices != other.NumVertices)
            {
                return NumVertices - other.NumVertices;
            }
            // Compare the two loops' vertices, starting with each loop's
            // firstLogicalVertex. This allows us to always catch cases where logically
            // identical loops have different vertex orderings (e.g. ABCD and BCDA).
            var maxVertices = NumVertices;
            var iThis = _firstLogicalVertex;
            var iOther = other._firstLogicalVertex;
            for (var i = 0; i < maxVertices; ++i, ++iThis, ++iOther)
            {
                var compare = Vertex(iThis).CompareTo(other.Vertex(iOther));
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
            // It is faster to construct a bounding rectangle for an S2Cell than for
            // a general polygon. A future optimization could also take advantage of
            // the fact than an S2Cell is convex.

            var cellBound = cell.RectBound;
            if (!_bound.Contains(cellBound))
            {
                return false;
            }
            var cellLoop = new S2Loop(cell, cellBound);
            return Contains(cellLoop);
        }

        /**
   * If this method returns false, the region does not intersect the given cell.
   * Otherwise, either region intersects the cell, or the intersection
   * relationship could not be determined.
   */

        public bool MayIntersect(S2Cell cell)
        {
            // It is faster to construct a bounding rectangle for an S2Cell than for
            // a general polygon. A future optimization could also take advantage of
            // the fact than an S2Cell is convex.

            var cellBound = cell.RectBound;
            if (!_bound.Intersects(cellBound))
            {
                return false;
            }
            return new S2Loop(cell, cellBound).Intersects(this);
        }

        /**
* The depth of a loop is defined as its nesting level within its containing
* polygon. "Outer shell" loops have depth 0, holes within those loops have
* depth 1, shells within those holes have depth 2, etc. This field is only
* used by the S2Polygon implementation.
*
* @param depth
*/

        /**
   * For convenience, we make two entire copies of the vertex list available:
   * vertex(n..2*n-1) is mapped to vertex(0..n-1), where n == numVertices().
   */

        public S2Point Vertex(int i)
        {
            try
            {
                return _vertices[i >= _vertices.Length ? i - _vertices.Length : i];
            }
            catch (IndexOutOfRangeException)
            {
                throw new InvalidOperationException("Invalid vertex index");
            }
        }

        /**
   * Comparator (needed by Comparable interface)
   */

        /**
   * Calculates firstLogicalVertex, the vertex in this loop that comes first in
   * a total ordering of all vertices (by way of S2Point's compareTo function).
   */

        private void InitFirstLogicalVertex()
        {
            var first = 0;
            for (var i = 1; i < _numVertices; ++i)
            {
                if (Vertex(i).CompareTo(Vertex(first)) < 0)
                {
                    first = i;
                }
            }
            _firstLogicalVertex = first;
        }

        /**
   * Return true if the loop area is at most 2*Pi.
   */

        /**
   * Invert the loop if necessary so that the area enclosed by the loop is at
   * most 2*Pi.
   */

        public void Normalize()
        {
            if (!IsNormalized)
            {
                Invert();
            }
        }

        /**
   * Reverse the order of the loop vertices, effectively complementing the
   * region represented by the loop.
   */

        public void Invert()
        {
            var last = NumVertices - 1;
            for (var i = (last - 1)/2; i >= 0; --i)
            {
                var t = _vertices[i];
                _vertices[i] = _vertices[last - i];
                _vertices[last - i] = t;
            }
            _vertexToIndex = null;
            _index = null;
            _originInside ^= true;
            if (_bound.Lat.Lo > -S2.PiOver2 && _bound.Lat.Hi < S2.PiOver2)
            {
                // The complement of this loop contains both poles.
                _bound = S2LatLngRect.Full;
            }
            else
            {
                InitBound();
            }
            InitFirstLogicalVertex();
        }

        /**
   * Helper method to get area and optionally centroid.
   */

        private S2AreaCentroid GetAreaCentroid(bool doCentroid)
        {
            // Don't crash even if loop is not well-defined.
            if (NumVertices < 3)
            {
                return new S2AreaCentroid(0D);
            }

            // The triangle area calculation becomes numerically unstable as the length
            // of any edge approaches 180 degrees. However, a loop may contain vertices
            // that are 180 degrees apart and still be valid, e.g. a loop that defines
            // the northern hemisphere using four points. We handle this case by using
            // triangles centered around an origin that is slightly displaced from the
            // first vertex. The amount of displacement is enough to get plenty of
            // accuracy for antipodal points, but small enough so that we still get
            // accurate areas for very tiny triangles.
            //
            // Of course, if the loop contains a point that is exactly antipodal from
            // our slightly displaced vertex, the area will still be unstable, but we
            // expect this case to be very unlikely (i.e. a polygon with two vertices on
            // opposite sides of the Earth with one of them displaced by about 2mm in
            // exactly the right direction). Note that the approximate point resolution
            // using the E7 or S2CellId representation is only about 1cm.

            var origin = Vertex(0);
            var axis = (origin.LargestAbsComponent + 1)%3;
            var slightlyDisplaced = origin[axis] + S2.E*1e-10;
            origin =
                new S2Point((axis == 0) ? slightlyDisplaced : origin.X,
                            (axis == 1) ? slightlyDisplaced : origin.Y, (axis == 2) ? slightlyDisplaced : origin.Z);
            origin = S2Point.Normalize(origin);

            double areaSum = 0;
            var centroidSum = new S2Point(0, 0, 0);
            for (var i = 1; i <= NumVertices; ++i)
            {
                areaSum += S2.SignedArea(origin, Vertex(i - 1), Vertex(i));
                if (doCentroid)
                {
                    // The true centroid is already premultiplied by the triangle area.
                    var trueCentroid = S2.TrueCentroid(origin, Vertex(i - 1), Vertex(i));
                    centroidSum = centroidSum + trueCentroid;
                }
            }
            // The calculated area at this point should be between -4*Pi and 4*Pi,
            // although it may be slightly larger or smaller than this due to
            // numerical errors.
            // assert (Math.abs(areaSum) <= 4 * S2.M_PI + 1e-12);

            if (areaSum < 0)
            {
                // If the area is negative, we have computed the area to the right of the
                // loop. The area to the left is 4*Pi - (-area). Amazingly, the centroid
                // does not need to be changed, since it is the negative of the integral
                // of position over the region to the right of the loop. This is the same
                // as the integral of position over the region to the left of the loop,
                // since the integral of position over the entire sphere is (0, 0, 0).
                areaSum += 4*S2.Pi;
            }
            // The loop's sign() does not affect the return result and should be taken
            // into account by the caller.
            S2Point? centroid = null;
            if (doCentroid)
            {
                centroid = centroidSum;
            }
            return new S2AreaCentroid(areaSum, centroid);
        }

        /**
   * Return the area of the loop interior, i.e. the region on the left side of
   * the loop. The return value is between 0 and 4*Pi and the true centroid of
   * the loop multiplied by the area of the loop (see S2.java for details on
   * centroids). Note that the centroid may not be contained by the loop.
   */

        // The following are the possible relationships between two loops A and B:
        //
        // (1) A and B do not intersect.
        // (2) A contains B.
        // (3) B contains A.
        // (4) The boundaries of A and B cross (i.e. the boundary of A
        // intersects the interior and exterior of B and vice versa).
        // (5) (A union B) is the entire sphere (i.e. A contains the
        // complement of B and vice versa).
        //
        // More than one of these may be true at the same time, for example if
        // A == B or A == Complement(B).

        /**
   * Return true if the region contained by this loop is a superset of the
   * region contained by the given other loop.
   */

        public bool Contains(S2Loop b)
        {
            // For this loop A to contains the given loop B, all of the following must
            // be true:
            //
            // (1) There are no edge crossings between A and B except at vertices.
            //
            // (2) At every vertex that is shared between A and B, the local edge
            // ordering implies that A contains B.
            //
            // (3) If there are no shared vertices, then A must contain a vertex of B
            // and B must not contain a vertex of A. (An arbitrary vertex may be
            // chosen in each case.)
            //
            // The second part of (3) is necessary to detect the case of two loops whose
            // union is the entire sphere, i.e. two loops that contains each other's
            // boundaries but not each other's interiors.

            if (!_bound.Contains(b.RectBound))
            {
                return false;
            }

            // Unless there are shared vertices, we need to check whether A contains a
            // vertex of B. Since shared vertices are rare, it is more efficient to do
            // this test up front as a quick rejection test.
            if (!Contains(b.Vertex(0)) && FindVertex(b.Vertex(0)) < 0)
            {
                return false;
            }

            // Now check whether there are any edge crossings, and also check the loop
            // relationship at any shared vertices.
            if (CheckEdgeCrossings(b, new WedgeContains()) <= 0)
            {
                return false;
            }

            // At this point we know that the boundaries of A and B do not intersect,
            // and that A contains a vertex of B. However we still need to check for
            // the case mentioned above, where (A union B) is the entire sphere.
            // Normally this check is very cheap due to the bounding box precondition.
            if (_bound.Union(b.RectBound).IsFull)
            {
                if (b.Contains(Vertex(0)) && b.FindVertex(Vertex(0)) < 0)
                {
                    return false;
                }
            }
            return true;
        }

        /**
   * Return true if the region contained by this loop intersects the region
   * contained by the given other loop.
   */

        public bool Intersects(S2Loop b)
        {
            // a->Intersects(b) if and only if !a->Complement()->Contains(b).
            // This code is similar to Contains(), but is optimized for the case
            // where both loops enclose less than half of the sphere.

            if (!_bound.Intersects(b.RectBound))
            {
                return false;
            }

            // Normalize the arguments so that B has a smaller longitude span than A.
            // This makes intersection tests much more efficient in the case where
            // longitude pruning is used (see CheckEdgeCrossings).
            if (b.RectBound.Lng.Length > _bound.Lng.Length)
            {
                return b.Intersects(this);
            }

            // Unless there are shared vertices, we need to check whether A contains a
            // vertex of B. Since shared vertices are rare, it is more efficient to do
            // this test up front as a quick acceptance test.
            if (Contains(b.Vertex(0)) && FindVertex(b.Vertex(0)) < 0)
            {
                return true;
            }

            // Now check whether there are any edge crossings, and also check the loop
            // relationship at any shared vertices.
            if (CheckEdgeCrossings(b, new WedgeIntersects()) < 0)
            {
                return true;
            }

            // We know that A does not contain a vertex of B, and that there are no edge
            // crossings. Therefore the only way that A can intersect B is if B
            // entirely contains A. We can check this by testing whether B contains an
            // arbitrary non-shared vertex of A. Note that this check is cheap because
            // of the bounding box precondition and the fact that we normalized the
            // arguments so that A's longitude span is at least as long as B's.
            if (b.RectBound.Contains(_bound))
            {
                if (b.Contains(Vertex(0)) && b.FindVertex(Vertex(0)) < 0)
                {
                    return true;
                }
            }

            return false;
        }

        /**
   * Given two loops of a polygon, return true if A contains B. This version of
   * contains() is much cheaper since it does not need to check whether the
   * boundaries of the two loops cross.
   */

        public bool ContainsNested(S2Loop b)
        {
            if (!_bound.Contains(b.RectBound))
            {
                return false;
            }

            // We are given that A and B do not share any edges, and that either one
            // loop contains the other or they do not intersect.
            var m = FindVertex(b.Vertex(1));
            if (m < 0)
            {
                // Since b->vertex(1) is not shared, we can check whether A contains it.
                return Contains(b.Vertex(1));
            }
            // Check whether the edge order around b->vertex(1) is compatible with
            // A containin B.
            return (new WedgeContains()).Test(
                Vertex(m - 1), Vertex(m), Vertex(m + 1), b.Vertex(0), b.Vertex(2)) > 0;
        }

        /**
   * Return +1 if A contains B (i.e. the interior of B is a subset of the
   * interior of A), -1 if the boundaries of A and B cross, and 0 otherwise.
   * Requires that A does not properly contain the complement of B, i.e. A and B
   * do not contain each other's boundaries. This method is used for testing
   * whether multi-loop polygons contain each other.
   */

        public int ContainsOrCrosses(S2Loop b)
        {
            // There can be containment or crossing only if the bounds intersect.
            if (!_bound.Intersects(b.RectBound))
            {
                return 0;
            }

            // Now check whether there are any edge crossings, and also check the loop
            // relationship at any shared vertices. Note that unlike Contains() or
            // Intersects(), we can't do a point containment test as a shortcut because
            // we need to detect whether there are any edge crossings.
            var result = CheckEdgeCrossings(b, new WedgeContainsOrCrosses());

            // If there was an edge crossing or a shared vertex, we know the result
            // already. (This is true even if the result is 1, but since we don't
            // bother keeping track of whether a shared vertex was seen, we handle this
            // case below.)
            if (result <= 0)
            {
                return result;
            }

            // At this point we know that the boundaries do not intersect, and we are
            // given that (A union B) is a proper subset of the sphere. Furthermore
            // either A contains B, or there are no shared vertices (due to the check
            // above). So now we just need to distinguish the case where A contains B
            // from the case where B contains A or the two loops are disjoint.
            if (!_bound.Contains(b.RectBound))
            {
                return 0;
            }
            if (!Contains(b.Vertex(0)) && FindVertex(b.Vertex(0)) < 0)
            {
                return 0;
            }

            return 1;
        }

        /**
   * Returns true if two loops have the same boundary except for vertex
   * perturbations. More precisely, the vertices in the two loops must be in the
   * same cyclic order, and corresponding vertex pairs must be separated by no
   * more than maxError. Note: This method mostly useful only for testing
   * purposes.
   */

        internal bool BoundaryApproxEquals(S2Loop b, double maxError)
        {
            if (NumVertices != b.NumVertices)
            {
                return false;
            }
            var maxVertices = NumVertices;
            var iThis = _firstLogicalVertex;
            var iOther = b._firstLogicalVertex;
            for (var i = 0; i < maxVertices; ++i, ++iThis, ++iOther)
            {
                if (!S2.ApproxEquals(Vertex(iThis), b.Vertex(iOther), maxError))
                {
                    return false;
                }
            }
            return true;
        }

        // S2Region interface (see {@code S2Region} for details):

        /** Return a bounding spherical cap. */

        /**
   * The point 'p' does not need to be normalized.
   */

        public bool Contains(S2Point p)
        {
            if (!_bound.Contains(p))
            {
                return false;
            }

            var inside = _originInside;
            var origin = S2.Origin;
            var crosser = new EdgeCrosser(origin, p,
                                          _vertices[_numVertices - 1]);

            // The s2edgeindex library is not optimized yet for long edges,
            // so the tradeoff to using it comes with larger loops.
            if (_numVertices < 2000)
            {
                for (var i = 0; i < _numVertices; i++)
                {
                    inside ^= crosser.EdgeOrVertexCrossing(_vertices[i]);
                }
            }
            else
            {
                var it = GetEdgeIterator(_numVertices);
                it.GetCandidates(origin, p);
                var previousIndex = -2;
                foreach (var ai in it) // it.GetCandidates(origin, p); it.HasNext; it.Next())
                {
                    //var ai = it.Index;
                    if (previousIndex != ai - 1)
                    {
                        crosser.RestartAt(_vertices[ai]);
                    }
                    previousIndex = ai;
                    inside ^= crosser.EdgeOrVertexCrossing(Vertex(ai + 1));
                }
            }

            return inside;
        }

        /**
   * Returns the shortest distance from a point P to this loop, given as the
   * angle formed between P, the origin and the nearest point on the loop to P.
   * This angle in radians is equivalent to the arclength along the unit sphere.
   */

        public S1Angle GetDistance(S2Point p)
        {
            var normalized = S2Point.Normalize(p);

            // The furthest point from p on the sphere is its antipode, which is an
            // angle of PI radians. This is an upper bound on the angle.
            var minDistance = S1Angle.FromRadians(Math.PI);
            for (var i = 0; i < NumVertices; i++)
            {
                minDistance =
                    S1Angle.Min(minDistance, S2EdgeUtil.GetDistance(normalized, Vertex(i), Vertex(i + 1)));
            }
            return minDistance;
        }

        /**
   * Creates an edge index over the vertices, which by itself takes no time.
   * Then the expected number of queries is used to determine whether brute
   * force lookups are likely to be slower than really creating an index, and if
   * so, we do so. Finally an iterator is returned that can be used to perform
   * edge lookups.
   */

        private S2EdgeIndex.DataEdgeIterator GetEdgeIterator(int expectedQueries)
        {
            if (_index == null)
            {
                _index = new AnonS2EdgeIndex(this);
            }
            _index.PredictAdditionalCalls(expectedQueries);

            return new S2EdgeIndex.DataEdgeIterator(_index);
        }


        /** Return true if this loop is valid. */

        /**
   * Static version of isValid(), to be used only when an S2Loop instance is not
   * available, but validity of the points must be checked.
   *
   * @return true if the given loop is valid. Creates an instance of S2Loop and
   *         defers this call to {@link #isValid()}.
   */

        public static bool IsValidLoop(IEnumerable<S2Point> vertices)
        {
            return new S2Loop(vertices).IsValid;
        }

        public override String ToString()
        {
            var builder = new StringBuilder("S2Loop, ");

            builder.Append(_vertices.Length).Append(" points. [");

            foreach (var v in _vertices)
            {
                builder.Append(v.ToString()).Append(" ");
            }
            builder.Append("]");

            return builder.ToString();
        }

        private void InitOrigin()
        {
            // The bounding box does not need to be correct before calling this
            // function, but it must at least contain vertex(1) since we need to
            // do a Contains() test on this point below.
            Preconditions.CheckState(_bound.Contains(Vertex(1)));

            // To ensure that every point is contained in exactly one face of a
            // subdivision of the sphere, all containment tests are done by counting the
            // edge crossings starting at a fixed point on the sphere (S2::Origin()).
            // We need to know whether this point is inside or outside of the loop.
            // We do this by first guessing that it is outside, and then seeing whether
            // we get the correct containment result for vertex 1. If the result is
            // incorrect, the origin must be inside the loop.
            //
            // A loop with consecutive vertices A,B,C contains vertex B if and only if
            // the fixed vector R = S2::Ortho(B) is on the left side of the wedge ABC.
            // The test below is written so that B is inside if C=R but not if A=R.

            _originInside = false; // Initialize before calling Contains().
            var v1Inside = S2.OrderedCcw(S2.Ortho(Vertex(1)), Vertex(0), Vertex(2), Vertex(1));
            if (v1Inside != Contains(Vertex(1)))
            {
                _originInside = true;
            }
        }

        private void InitBound()
        {
            // The bounding rectangle of a loop is not necessarily the same as the
            // bounding rectangle of its vertices. First, the loop may wrap entirely
            // around the sphere (e.g. a loop that defines two revolutions of a
            // candy-cane stripe). Second, the loop may include one or both poles.
            // Note that a small clockwise loop near the equator contains both poles.

            var bounder = new RectBounder();
            for (var i = 0; i <= NumVertices; ++i)
            {
                bounder.AddPoint(Vertex(i));
            }
            var b = bounder.Bound;
            // Note that we need to initialize bound with a temporary value since
            // contains() does a bounding rectangle check before doing anything else.
            _bound = S2LatLngRect.Full;
            if (Contains(new S2Point(0, 0, 1)))
            {
                b = new S2LatLngRect(new R1Interval(b.Lat.Lo, S2.PiOver2), S1Interval.Full);
            }
            // If a loop contains the south pole, then either it wraps entirely
            // around the sphere (full longitude range), or it also contains the
            // north pole in which case b.lng().isFull() due to the test above.

            if (b.Lng.IsFull && Contains(new S2Point(0, 0, -1)))
            {
                b = new S2LatLngRect(new R1Interval(-S2.PiOver2, b.Lat.Hi), b.Lng);
            }
            _bound = b;
        }

        /**
   * Return the index of a vertex at point "p", or -1 if not found. The return
   * value is in the range 1..num_vertices_ if found.
   */

        private int FindVertex(S2Point p)
        {
            if (_vertexToIndex == null)
            {
                _vertexToIndex = new Dictionary<S2Point, int>();
                for (var i = 1; i <= _numVertices; i++)
                {
                    _vertexToIndex[Vertex(i)] = i;
                }
            }

            if (!_vertexToIndex.ContainsKey(p))
            {
                return -1;
            }

            return _vertexToIndex[p];
        }

        /**
   * This method encapsulates the common code for loop containment and
   * intersection tests. It is used in three slightly different variations to
   * implement contains(), intersects(), and containsOrCrosses().
   *
   *  In a nutshell, this method checks all the edges of this loop (A) for
   * intersection with all the edges of B. It returns -1 immediately if any edge
   * intersections are found. Otherwise, if there are any shared vertices, it
   * returns the minimum value of the given WedgeRelation for all such vertices
   * (returning immediately if any wedge returns -1). Returns +1 if there are no
   * intersections and no shared vertices.
   */

        private int CheckEdgeCrossings(S2Loop b, IWedgeRelation relation)
        {
            var it = GetEdgeIterator(b._numVertices);
            var result = 1;
            // since 'this' usually has many more vertices than 'b', use the index on
            // 'this' and loop over 'b'
            for (var j = 0; j < b.NumVertices; ++j)
            {
                var crosser =
                    new EdgeCrosser(b.Vertex(j), b.Vertex(j + 1), Vertex(0));
                var previousIndex = -2;

                it.GetCandidates(b.Vertex(j), b.Vertex(j + 1));
                foreach (var i in it) // it.GetCandidates(b.vertex(j), b.vertex(j + 1)); it.HasNext; it.Next())
                {
                    //    var i = it.Index;
                    if (previousIndex != i - 1)
                    {
                        crosser.RestartAt(Vertex(i));
                    }
                    previousIndex = i;
                    var crossing = crosser.RobustCrossing(Vertex(i + 1));
                    if (crossing < 0)
                    {
                        continue;
                    }
                    if (crossing > 0)
                    {
                        return -1; // There is a proper edge crossing.
                    }
                    if (Vertex(i + 1).Equals(b.Vertex(j + 1)))
                    {
                        result = Math.Min(result, relation.Test(
                            Vertex(i), Vertex(i + 1), Vertex(i + 2), b.Vertex(j), b.Vertex(j + 2)));
                        if (result < 0)
                        {
                            return result;
                        }
                    }
                }
            }
            return result;
        }

        private sealed class AnonS2EdgeIndex : S2EdgeIndex
        {
            private readonly S2Loop _this;

            public AnonS2EdgeIndex(S2Loop This)
            {
                _this = This;
            }

            protected override int NumEdges
            {
                get { return _this._numVertices; }
            }

            protected override S2Point EdgeFrom(int index)
            {
                return _this.Vertex(index);
            }

            protected override S2Point EdgeTo(int index)
            {
                return _this.Vertex(index + 1);
            }
        }
    }
}