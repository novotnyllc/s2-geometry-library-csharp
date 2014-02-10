using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Google.Common.Geometry
{
    /**
     * An S2Cell is an S2Region object that represents a cell. Unlike S2CellIds, it
     * supports efficient containment and intersection tests. However, it is also a
     * more expensive representation.
     *
     */

    public sealed class S2Cell : IS2Region, IEquatable<S2Cell>
    {
        private const int MaxCellSize = 1 << S2CellId.MaxLevel;
        private const double MaxError = 1.0/(1L << 51);

        // The 4 cells around the equator extend to +/-45 degrees latitude at the
        // midpoints of their top and bottom edges. The two cells covering the
        // poles extend down to +/-35.26 degrees at their vertices.
        // adding kMaxError (as opposed to the C version) because of asin and atan2
        // roundoff errors
        private static readonly double PoleMinLat = Math.Asin(Math.Sqrt(1.0/3.0)) - MaxError;
        private S2CellId _cellId;

        private byte _face;
        private byte _level;
        private byte _orientation;

        private double[][] _uv = new double[2][]
        {
            new double[2],
            new double[2],
        };

        /**
   * Default constructor used only internally.
   */

        internal S2Cell()
        {
        }

        /**
   * An S2Cell always corresponds to a particular S2CellId. The other
   * constructors are just convenience methods.
   */

        public S2Cell(S2CellId id)
        {
            Init(id);
        }

        // This is a static method in order to provide named parameters.

        // Convenience methods.
        public S2Cell(S2Point p)
        {
            Init(S2CellId.FromPoint(p));
        }

        public S2Cell(S2LatLng ll)
        {
            Init(S2CellId.FromLatLng(ll));
        }

        public S2CellId Id
        {
            get { return _cellId; }
        }

        public int Face
        {
            get { return _face; }
        }

        public byte Level
        {
            get { return _level; }
        }

        public byte Orientation
        {
            get { return _orientation; }
        }

        public bool IsLeaf
        {
            get { return _level == S2CellId.MaxLevel; }
        }

        public S2Point Center
        {
            get { return S2Point.Normalize(CenterRaw); }
        }

        public S2Point CenterRaw
        {
            get { return _cellId.ToPointRaw(); }
        }

        /**
   * Return the center of the cell in (u,v) coordinates (see {@code
   * S2Projections}). Note that the center of the cell is defined as the point
   * at which it is recursively subdivided into four children; in general, it is
   * not at the midpoint of the (u,v) rectangle covered by the cell
   */

        public R2Vector CenterUv
        {
            get
            {
                var i = 0;
                var j = 0;
                int? notUsed = null;
                _cellId.ToFaceIjOrientation(ref i, ref j, ref notUsed);
                var cellSize = 1 << (S2CellId.MaxLevel - _level);

                // TODO(dbeaumont): Figure out a better naming of the variables here (and elsewhere).
                var si = (i & -cellSize)*2 + cellSize - MaxCellSize;
                var x = S2Projections.stToUV((1.0/MaxCellSize)*si);

                var sj = (j & -cellSize)*2 + cellSize - MaxCellSize;
                var y = S2Projections.stToUV((1.0/MaxCellSize)*sj);

                return new R2Vector(x, y);
            }
        }

        public bool Equals(S2Cell other)
        {
            if (ReferenceEquals(null, other)) return false;
            if (ReferenceEquals(this, other)) return true;
            return _cellId.Equals(other._cellId) && _level == other._level && _face == other._face && _orientation == other._orientation;
        }

        public S2Cap CapBound
        {
            get
            {
                // Use the cell center in (u,v)-space as the cap axis. This vector is
                // very close to GetCenter() and faster to compute. Neither one of these
                // vectors yields the bounding cap with minimal surface area, but they
                // are both pretty close.
                //
                // It's possible to show that the two vertices that are furthest from
                // the (u,v)-origin never determine the maximum cap size (this is a
                // possible future optimization).

                var u = 0.5*(_uv[0][0] + _uv[0][1]);
                var v = 0.5*(_uv[1][0] + _uv[1][1]);
                var cap = S2Cap.FromAxisHeight(S2Point.Normalize(S2Projections.faceUvToXyz(_face, u, v)), 0);
                for (var k = 0; k < 4; ++k)
                {
                    cap = cap.AddPoint(GetVertex(k));
                }
                return cap;
            }
        }

        public S2LatLngRect RectBound
        {
            get
            {
                if (_level > 0)
                {
                    // Except for cells at level 0, the latitude and longitude extremes are
                    // attained at the vertices. Furthermore, the latitude range is
                    // determined by one pair of diagonally opposite vertices and the
                    // longitude range is determined by the other pair.
                    //
                    // We first determine which corner (i,j) of the cell has the largest
                    // absolute latitude. To maximize latitude, we want to find the point in
                    // the cell that has the largest absolute z-coordinate and the smallest
                    // absolute x- and y-coordinates. To do this we look at each coordinate
                    // (u and v), and determine whether we want to minimize or maximize that
                    // coordinate based on the axis direction and the cell's (u,v) quadrant.
                    var u = _uv[0][0] + _uv[0][1];
                    var v = _uv[1][0] + _uv[1][1];
                    var i = S2Projections.getUAxis(_face).Z == 0 ? (u < 0 ? 1 : 0) : (u > 0 ? 1 : 0);
                    var j = S2Projections.getVAxis(_face).Z == 0 ? (v < 0 ? 1 : 0) : (v > 0 ? 1 : 0);


                    var lat = R1Interval.FromPointPair(GetLatitude(i, j), GetLatitude(1 - i, 1 - j));
                    lat = lat.Expanded(MaxError).Intersection(S2LatLngRect.fullLat());
                    if (lat.Lo == -S2.PiOver2 || lat.Hi == S2.PiOver2)
                    {
                        return new S2LatLngRect(lat, S1Interval.Full);
                    }
                    var lng = S1Interval.FromPointPair(GetLongitude(i, 1 - j), GetLongitude(1 - i, j));
                    return new S2LatLngRect(lat, lng.Expanded(MaxError));
                }


                // The face centers are the +X, +Y, +Z, -X, -Y, -Z axes in that order.
                // assert (S2Projections.getNorm(face).get(face % 3) == ((face < 3) ? 1 : -1));
                switch (_face)
                {
                    case 0:
                        return new S2LatLngRect(
                            new R1Interval(-S2.PiOver4, S2.PiOver4), new S1Interval(-S2.PiOver4, S2.PiOver4));
                    case 1:
                        return new S2LatLngRect(
                            new R1Interval(-S2.PiOver4, S2.PiOver4), new S1Interval(S2.PiOver4, 3*S2.PiOver4));
                    case 2:
                        return new S2LatLngRect(
                            new R1Interval(PoleMinLat, S2.PiOver2), new S1Interval(-S2.Pi, S2.Pi));
                    case 3:
                        return new S2LatLngRect(
                            new R1Interval(-S2.PiOver4, S2.PiOver4), new S1Interval(3*S2.PiOver4, -3*S2.PiOver4));
                    case 4:
                        return new S2LatLngRect(
                            new R1Interval(-S2.PiOver4, S2.PiOver4), new S1Interval(-3*S2.PiOver4, -S2.PiOver4));
                    default:
                        return new S2LatLngRect(
                            new R1Interval(-S2.PiOver2, -PoleMinLat), new S1Interval(-S2.Pi, S2.Pi));
                }
            }
        }

        public bool MayIntersect(S2Cell cell)
        {
            return _cellId.Intersects(cell._cellId);
        }

        public bool Contains(S2Cell cell)
        {
            return _cellId.Contains(cell._cellId);
        }

        public override bool Equals(object obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            if (ReferenceEquals(this, obj)) return true;
            if (obj.GetType() != GetType()) return false;
            return Equals((S2Cell)obj);
        }

        public override int GetHashCode()
        {
            unchecked
            {
                var hashCode = _cellId.GetHashCode();
                hashCode = (hashCode*397) ^ _level.GetHashCode();
                hashCode = (hashCode*397) ^ _face.GetHashCode();
                hashCode = (hashCode*397) ^ _orientation.GetHashCode();
                return hashCode;
            }
        }

        public static bool operator ==(S2Cell left, S2Cell right)
        {
            return Equals(left, right);
        }

        public static bool operator !=(S2Cell left, S2Cell right)
        {
            return !Equals(left, right);
        }

        public static S2Cell FromFacePosLevel(int face, byte pos, int level)
        {
            return new S2Cell(S2CellId.FromFacePosLevel(face, pos, level));
        }


        public S2Point GetVertex(int k)
        {
            return S2Point.Normalize(GetVertexRaw(k));
        }

        /**
   * Return the k-th vertex of the cell (k = 0,1,2,3). Vertices are returned in
   * CCW order. The points returned by GetVertexRaw are not necessarily unit
   * length.
   */

        public S2Point GetVertexRaw(int k)
        {
            // Vertices are returned in the order SW, SE, NE, NW.
            return S2Projections.faceUvToXyz(_face, _uv[0][(k >> 1) ^ (k & 1)], _uv[1][k >> 1]);
        }

        public S2Point GetEdge(int k)
        {
            return S2Point.Normalize(GetEdgeRaw(k));
        }

        public S2Point GetEdgeRaw(int k)
        {
            switch (k)
            {
                case 0:
                    return S2Projections.getVNorm(_face, _uv[1][0]); // South
                case 1:
                    return S2Projections.getUNorm(_face, _uv[0][1]); // East
                case 2:
                    return -S2Projections.getVNorm(_face, _uv[1][1]); // North
                default:
                    return -S2Projections.getUNorm(_face, _uv[0][0]); // West
            }
        }

        /**
   * Return the inward-facing normal of the great circle passing through the
   * edge from vertex k to vertex k+1 (mod 4). The normals returned by
   * GetEdgeRaw are not necessarily unit length.
   *
   *  If this is not a leaf cell, set children[0..3] to the four children of
   * this cell (in traversal order) and return true. Otherwise returns false.
   * This method is equivalent to the following:
   *
   *  for (pos=0, id=child_begin(); id != child_end(); id = id.next(), ++pos)
   * children[i] = S2Cell(id);
   *
   * except that it is more than two times faster.
   */

        public bool Subdivide(IReadOnlyList<S2Cell> children)
        {
            // This function is equivalent to just iterating over the child cell ids
            // and calling the S2Cell constructor, but it is about 2.5 times faster.

            if (_cellId.IsLeaf)
            {
                return false;
            }

            // Compute the cell midpoint in uv-space.
            var uvMid = CenterUv;

            // Create four children with the appropriate bounds.
            var id = _cellId.ChildBegin;
            for (var pos = 0; pos < 4; ++pos, id = id.Next)
            {
                var child = children[pos];
                child._face = _face;
                child._level = (byte)(_level + 1);
                child._orientation = (byte)(_orientation ^ S2.PosToOrientation(pos));
                child._cellId = id;
                var ij = S2.PosToIj(_orientation, pos);
                for (var d = 0; d < 2; ++d)
                {
                    // The dimension 0 index (i/u) is in bit 1 of ij.
                    var m = 1 - ((ij >> (1 - d)) & 1);
                    child._uv[d][m] = uvMid[d];
                    child._uv[d][1 - m] = _uv[d][1 - m];
                }
            }
            return true;
        }

        /**
   * Return the direction vector corresponding to the center in (s,t)-space of
   * the given cell. This is the point at which the cell is divided into four
   * subcells; it is not necessarily the centroid of the cell in (u,v)-space or
   * (x,y,z)-space. The point returned by GetCenterRaw is not necessarily unit
   * length.
   */

        /**
   * Return the average area for cells at the given level.
   */

        public static double AverageArea(int level)
        {
            return S2Projections.AVG_AREA.GetValue(level);
        }

        /**
   * Return the average area of cells at this level. This is accurate to within
   * a factor of 1.7 (for S2_QUADRATIC_PROJECTION) and is extremely cheap to
   * compute.
   */

        public double AverageArea()
        {
            return AverageArea(_level);
        }

        /**
   * Return the approximate area of this cell. This method is accurate to within
   * 3% percent for all cell sizes and accurate to within 0.1% for cells at
   * level 5 or higher (i.e. 300km square or smaller). It is moderately cheap to
   * compute.
   */

        public double ApproxArea()
        {
            // All cells at the first two levels have the same area.
            if (_level < 2)
            {
                return AverageArea(_level);
            }

            // First, compute the approximate area of the cell when projected
            // perpendicular to its normal. The cross product of its diagonals gives
            // the normal, and the length of the normal is twice the projected area.
            var flatArea = 0.5*S2Point.CrossProd(
                GetVertex(2) - GetVertex(0), GetVertex(3) - GetVertex(1)).Norm;

            // Now, compensate for the curvature of the cell surface by pretending
            // that the cell is shaped like a spherical cap. The ratio of the
            // area of a spherical cap to the area of its projected disc turns out
            // to be 2 / (1 + sqrt(1 - r*r)) where "r" is the radius of the disc.
            // For example, when r=0 the ratio is 1, and when r=1 the ratio is 2.
            // Here we set Pi*r*r == flat_area to find the equivalent disc.
            return flatArea*2/(1 + Math.Sqrt(1 - Math.Min(S2.InversePi*flatArea, 1.0)));
        }

        /**
   * Return the area of this cell as accurately as possible. This method is more
   * expensive but it is accurate to 6 digits of precision even for leaf cells
   * (whose area is approximately 1e-18).
   */

        public double ExactArea()
        {
            var v0 = GetVertex(0);
            var v1 = GetVertex(1);
            var v2 = GetVertex(2);
            var v3 = GetVertex(3);
            return S2.Area(v0, v1, v2) + S2.Area(v0, v2, v3);
        }

        // //////////////////////////////////////////////////////////////////////
        // S2Region interface (see {@code S2Region} for details):

        public IS2Region Clone()
        {
            var clone = new S2Cell();
            clone._face = _face;
            clone._level = _level;
            clone._orientation = _orientation;

            clone._uv = (double[][])_uv.Clone();

            return clone;
        }

        public bool Contains(S2Point p)
        {
            // We can't just call XYZtoFaceUV, because for points that lie on the
            // boundary between two faces (i.e. u or v is +1/-1) we need to return
            // true for both adjacent cells.
            var uvPoint = S2Projections.faceXyzToUv(_face, p);
            if (uvPoint == null)
            {
                return false;
            }
            return (uvPoint.Value.X >= _uv[0][0] && uvPoint.Value.X <= _uv[0][1]
                    && uvPoint.Value.Y >= _uv[1][0] && uvPoint.Value.Y <= _uv[1][1]);
        }

        // The point 'p' does not need to be normalized.

        private void Init(S2CellId id)
        {
            _cellId = id;
            var ij = new int[2];
            int? mOrientation = 0;

            for (var d = 0; d < 2; ++d)
            {
                ij[d] = 0;
            }

            _face = (byte)id.ToFaceIjOrientation(ref ij[0], ref ij[1], ref mOrientation);
            _orientation = (byte)mOrientation.Value; // Compress int to a byte.
            _level = (byte)id.Level;
            var cellSize = 1 << (S2CellId.MaxLevel - _level);
            for (var d = 0; d < 2; ++d)
            {
                // Compute the cell bounds in scaled (i,j) coordinates.
                var sijLo = (ij[d] & -cellSize)*2 - MaxCellSize;
                var sijHi = sijLo + cellSize*2;
                _uv[d][0] = S2Projections.stToUV((1.0/MaxCellSize)*sijLo);
                _uv[d][1] = S2Projections.stToUV((1.0/MaxCellSize)*sijHi);
            }
        }


        // Internal method that does the actual work in the constructors.

        private double GetLatitude(int i, int j)
        {
            var p = S2Projections.faceUvToXyz(_face, _uv[0][i], _uv[1][j]);
            return Math.Atan2(p.Z, Math.Sqrt(p.X*p.X + p.Y*p.Y));
        }

        private double GetLongitude(int i, int j)
        {
            var p = S2Projections.faceUvToXyz(_face, _uv[0][i], _uv[1][j]);
            return Math.Atan2(p.Y, p.X);
        }

        // Return the latitude or longitude of the cell vertex given by (i,j),
        // where "i" and "j" are either 0 or 1.


        public override String ToString()
        {
            return "[" + _face + ", " + _level + ", " + _orientation + ", " + _cellId + "]";
        }
    }
}