using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Google.Common.Geometry
{
    /**
     * An S2LatLngRect represents a latitude-longitude rectangle. It is capable of
     * representing the empty and full rectangles as well as single points.
     *
     */

    public struct S2LatLngRect : IS2Region, IEquatable<S2LatLngRect>
    {
        public static readonly S2LatLngRect Empty = new S2LatLngRect(R1Interval.Empty, S1Interval.Empty);


        /** The full allowable range of latitudes. */

        public static readonly R1Interval FullLat = new R1Interval(-S2.PiOver2, S2.PiOver2);

        /**
   * The full allowable range of longitudes.
   */

        public static readonly S1Interval FullLng = S1Interval.Full;

        /** The canonical full rectangle. */

        public static readonly S2LatLngRect Full = new S2LatLngRect(FullLat, FullLng);
        private readonly R1Interval _lat;
        private readonly S1Interval _lng;

        /**
   * Construct a rectangle from minimum and maximum latitudes and longitudes. If
   * lo.Lng > hi.Lng, the rectangle spans the 180 degree longitude line.
   */

        public S2LatLngRect(S2LatLng lo, S2LatLng hi)
        {
            _lat = new R1Interval(lo.Lat.Radians, hi.Lat.Radians);
            _lng = new S1Interval(lo.Lng.Radians, hi.Lng.Radians);
            // assert (isValid());
        }

        /** Construct a rectangle from latitude and longitude intervals. */

        public S2LatLngRect(R1Interval lat, S1Interval lng)
        {
            _lat = lat;
            _lng = lng;
            // assert (isValid());
        }

        public R1Interval Lat
        {
            get { return _lat; }
        }

        public S1Interval Lng
        {
            get { return _lng; }
        }

        public bool IsValid
        {
            get
            {
                // The lat/lng ranges must either be both empty or both non-empty.
                return (Math.Abs(_lat.Lo) <= S2.PiOver2 && Math.Abs(_lat.Hi) <= S2.PiOver2
                        && _lng.IsValid && _lat.IsEmpty == _lng.IsEmpty);
            }
        }

        // Accessor methods.

        public S1Angle LatLo
        {
            get { return S1Angle.FromRadians(_lat.Lo); }
        }

        public S1Angle LatHi
        {
            get { return S1Angle.FromRadians(_lat.Hi); }
        }

        public S1Angle LngLo
        {
            get { return S1Angle.FromRadians(_lng.Lo); }
        }

        public S1Angle LngHi
        {
            get { return S1Angle.FromRadians(_lng.Hi); }
        }

        public S2LatLng Lo
        {
            get { return new S2LatLng(LatLo, LngLo); }
        }

        public S2LatLng Hi
        {
            get { return new S2LatLng(LatHi, LngHi); }
        }

        /**
   * Return true if the rectangle is empty, i.e. it contains no points at all.
   */

        public bool IsEmpty
        {
            get { return _lat.IsEmpty; }
        }

        // Return true if the rectangle is full, i.e. it contains all points.

        public bool IsFull
        {
            get { return _lat.Equals(FullLat) && _lng.IsFull; }
        }

        /**
   * Return true if lng_.lo() > lng_.hi(), i.e. the rectangle crosses the 180
   * degree latitude line.
   */

        public bool IsInverted
        {
            get { return _lng.IsInverted; }
        }

        public S2LatLng Center
        {
            get { return S2LatLng.FromRadians(_lat.Center, _lng.Center); }
        }

        public S2LatLng Size
        {
            get { return S2LatLng.FromRadians(_lat.Length, _lng.Length); }
        }

        public double Area
        {
            get
            {
                if (IsEmpty)
                {
                    return 0;
                }

                // This is the size difference of the two spherical caps, multiplied by
                // the longitude ratio.
                return Lng.Length*Math.Abs(Math.Sin(LatHi.Radians) - Math.Sin(LatLo.Radians));
            }
        }

        public bool Equals(S2LatLngRect other)
        {
            return Equals(_lat, other._lat) && Equals(_lng, other._lng);
        }

        public S2Cap CapBound
        {
            get
            {
                // We consider two possible bounding caps, one whose axis passes
                // through the center of the lat-long rectangle and one whose axis
                // is the north or south pole. We return the smaller of the two caps.

                if (IsEmpty)
                {
                    return S2Cap.Empty;
                }

                double poleZ, poleAngle;
                if (_lat.Lo + _lat.Hi < 0)
                {
                    // South pole axis yields smaller cap.
                    poleZ = -1;
                    poleAngle = S2.PiOver2 + _lat.Hi;
                }
                else
                {
                    poleZ = 1;
                    poleAngle = S2.PiOver2 - _lat.Lo;
                }
                var poleCap = S2Cap.FromAxisAngle(new S2Point(0, 0, poleZ), S1Angle
                                                                                .FromRadians(poleAngle));

                // For bounding rectangles that span 180 degrees or less in longitude, the
                // maximum cap size is achieved at one of the rectangle vertices. For
                // rectangles that are larger than 180 degrees, we punt and always return a
                // bounding cap centered at one of the two poles.
                var lngSpan = _lng.Hi - _lng.Lo;
                if (Math.IEEERemainder(lngSpan, 2*S2.Pi) >= 0)
                {
                    if (lngSpan < 2*S2.Pi)
                    {
                        var midCap = S2Cap.FromAxisAngle(Center.ToPoint(), S1Angle
                                                                               .FromRadians(0));
                        for (var k = 0; k < 4; ++k)
                        {
                            midCap = midCap.AddPoint(GetVertex(k).ToPoint());
                        }
                        if (midCap.Height < poleCap.Height)
                        {
                            return midCap;
                        }
                    }
                }
                return poleCap;
            }
        }

        public S2LatLngRect RectBound
        {
            get { return this; }
        }


        public bool Contains(S2Cell cell)
        {
            // A latitude-longitude rectangle contains a cell if and only if it contains
            // the cell's bounding rectangle. (This is an exact test.)
            return Contains(cell.RectBound);
        }

        /**
   * This test is cheap but is NOT exact. Use Intersects() if you want a more
   * accurate and more expensive test. Note that when this method is used by an
   * S2RegionCoverer, the accuracy isn't all that important since if a cell may
   * intersect the region then it is subdivided, and the accuracy of this method
   * goes up as the cells get smaller.
   */

        public bool MayIntersect(S2Cell cell)
        {
            // This test is cheap but is NOT exact (see s2latlngrect.h).
            return Intersects(cell.RectBound);
        }

        public override bool Equals(object obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            return obj is S2LatLngRect && Equals((S2LatLngRect)obj);
        }

        public override int GetHashCode()
        {
            unchecked
            {
                return ((_lat != null ? _lat.GetHashCode() : 0)*397) ^ (_lng != null ? _lng.GetHashCode() : 0);
            }
        }

        public static bool operator ==(S2LatLngRect left, S2LatLngRect right)
        {
            return left.Equals(right);
        }

        public static bool operator !=(S2LatLngRect left, S2LatLngRect right)
        {
            return !left.Equals(right);
        }

        /** The canonical empty rectangle */

        /**
   * Construct a rectangle from a center point (in lat-lng space) and size in
   * each dimension. If size.Lng is greater than 360 degrees it is clamped,
   * and latitudes greater than +/- 90 degrees are also clamped. So for example,
   * FromCenterSize((80,170),(20,20)) -> (lo=(60,150),hi=(90,-170)).
   */

        public static S2LatLngRect FromCenterSize(S2LatLng center, S2LatLng size)
        {
            return FromPoint(center).Expanded(size*0.5);
        }

        /** Convenience method to construct a rectangle containing a single point. */

        public static S2LatLngRect FromPoint(S2LatLng p)
        {
            // assert (p.isValid());
            return new S2LatLngRect(p, p);
        }

        /**
   * Convenience method to construct the minimal bounding rectangle containing
   * the two given points. This is equivalent to starting with an empty
   * rectangle and calling AddPoint() twice. Note that it is different than the
   * S2LatLngRect(lo, hi) constructor, where the first point is always used as
   * the lower-left corner of the resulting rectangle.
   */

        public static S2LatLngRect FromPointPair(S2LatLng p1, S2LatLng p2)
        {
            // assert (p1.isValid() && p2.isValid());
            return new S2LatLngRect(R1Interval.FromPointPair(p1.Lat.Radians, p2.Lat.Radians), S1Interval.FromPointPair(p1.Lng.Radians, p2.Lng.Radians));
        }

        /**
   * Return a latitude-longitude rectangle that contains the edge from "a" to
   * "b". Both points must be unit-length. Note that the bounding rectangle of
   * an edge can be larger than the bounding rectangle of its endpoints.
   */

        public static S2LatLngRect FromEdge(S2Point a, S2Point b)
        {
            // assert (S2.isUnitLength(a) && S2.isUnitLength(b));
            var r = FromPointPair(new S2LatLng(a), new S2LatLng(b));

            // Check whether the min/max latitude occurs in the edge interior.
            // We find the normal to the plane containing AB, and then a vector "dir" in
            // this plane that also passes through the equator. We use RobustCrossProd
            // to ensure that the edge normal is accurate even when the two points are
            // very close together.
            var ab = S2.RobustCrossProd(a, b);
            var dir = S2Point.CrossProd(ab, new S2Point(0, 0, 1));
            var da = dir.DotProd(a);
            var db = dir.DotProd(b);
            if (da*db >= 0)
            {
                // Minimum and maximum latitude are attained at the vertices.
                return r;
            }
            // Minimum/maximum latitude occurs in the edge interior. This affects the
            // latitude bounds but not the longitude bounds.
            var absLat = Math.Acos(Math.Abs(ab.Z/ab.Norm));
            if (da < 0)
            {
                return new S2LatLngRect(new R1Interval(r.Lat.Lo, absLat), r.Lng);
            }
            else
            {
                return new S2LatLngRect(new R1Interval(-absLat, r.Lat.Hi), r.Lng);
            }
        }

        /**
   * Return true if the rectangle is valid, which essentially just means that
   * the latitude bounds do not exceed Pi/2 in absolute value and the longitude
   * bounds do not exceed Pi in absolute value.
   *
   */

        /** Return the k-th vertex of the rectangle (k = 0,1,2,3) in CCW order. */

        public S2LatLng GetVertex(int k)
        {
            // Return the points in CCW order (SW, SE, NE, NW).
            switch (k)
            {
                case 0:
                    return S2LatLng.FromRadians(_lat.Lo, _lng.Lo);
                case 1:
                    return S2LatLng.FromRadians(_lat.Lo, _lng.Hi);
                case 2:
                    return S2LatLng.FromRadians(_lat.Hi, _lng.Hi);
                case 3:
                    return S2LatLng.FromRadians(_lat.Hi, _lng.Lo);
                default:
                    throw new ArgumentException("Invalid vertex index.");
            }
        }

        /**
   * Return the center of the rectangle in latitude-longitude space (in general
   * this is not the center of the region on the sphere).
   */

        /**
   * Return the minimum distance (measured along the surface of the sphere)
   * from a given point to the rectangle (both its boundary and its interior).
   * The latLng must be valid.
   */

        public S1Angle GetDistance(S2LatLng p)
        {
            // The algorithm here is the same as in getDistance(S2LagLngRect), only
            // with simplified calculations.
            var a = this;

            Preconditions.CheckState(!a.IsEmpty);
            Preconditions.CheckArgument(p.IsValid);

            if (a.Lng.Contains(p.Lng.Radians))
            {
                return S1Angle.FromRadians(Math.Max(0.0, Math.Max(p.Lat.Radians - a.Lat.Hi,
                                                                  a.Lat.Lo - p.Lat.Radians)));
            }

            var interval = new S1Interval(a.Lng.Hi, a.Lng.Complement.Center);
            var aLng = a.Lng.Lo;
            if (interval.Contains(p.Lng.Radians))
            {
                aLng = a.Lng.Hi;
            }

            var lo = S2LatLng.FromRadians(a.Lat.Lo, aLng).ToPoint();
            var hi = S2LatLng.FromRadians(a.Lat.Hi, aLng).ToPoint();
            var loCrossHi =
                S2LatLng.FromRadians(0, aLng - S2.PiOver2).Normalized.ToPoint();
            return S2EdgeUtil.GetDistance(p.ToPoint(), lo, hi, loCrossHi);
        }

        /**
   * Return the minimum distance (measured along the surface of the sphere) to
   * the given S2LatLngRect. Both S2LatLngRects must be non-empty.
   */

        public S1Angle GetDistance(S2LatLngRect other)
        {
            var a = this;
            var b = other;

            Preconditions.CheckState(!a.IsEmpty);
            Preconditions.CheckArgument(!b.IsEmpty);

            // First, handle the trivial cases where the longitude intervals overlap.
            if (a.Lng.Intersects(b.Lng))
            {
                if (a.Lat.Intersects(b.Lat))
                {
                    return S1Angle.FromRadians(0); // Intersection between a and b.
                }

                // We found an overlap in the longitude interval, but not in the latitude
                // interval. This means the shortest path travels along some line of
                // longitude connecting the high-latitude of the lower rect with the
                // low-latitude of the higher rect.
                S1Angle lo, hi;
                if (a.Lat.Lo > b.Lat.Hi)
                {
                    lo = b.LatHi;
                    hi = a.LatLo;
                }
                else
                {
                    lo = a.LatHi;
                    hi = b.LatLo;
                }
                return S1Angle.FromRadians(hi.Radians - lo.Radians);
            }

            // The longitude intervals don't overlap. In this case, the closest points
            // occur somewhere on the pair of longitudinal edges which are nearest in
            // longitude-space.
            S1Angle aLng, bLng;
            var loHi = S1Interval.FromPointPair(a.Lng.Lo, b.Lng.Hi);
            var hiLo = S1Interval.FromPointPair(a.Lng.Hi, b.Lng.Lo);
            if (loHi.Length < hiLo.Length)
            {
                aLng = a.LngLo;
                bLng = b.LngHi;
            }
            else
            {
                aLng = a.LngHi;
                bLng = b.LngLo;
            }

            // The shortest distance between the two longitudinal segments will include
            // at least one segment endpoint. We could probably narrow this down further
            // to a single point-edge distance by comparing the relative latitudes of the
            // endpoints, but for the sake of clarity, we'll do all four point-edge
            // distance tests.
            var aLo = new S2LatLng(a.LatLo, aLng).ToPoint();
            var aHi = new S2LatLng(a.LatHi, aLng).ToPoint();
            var aLoCrossHi =
                S2LatLng.FromRadians(0, aLng.Radians - S2.PiOver2).Normalized.ToPoint();
            var bLo = new S2LatLng(b.LatLo, bLng).ToPoint();
            var bHi = new S2LatLng(b.LatHi, bLng).ToPoint();
            var bLoCrossHi =
                S2LatLng.FromRadians(0, bLng.Radians - S2.PiOver2).Normalized.ToPoint();

            return S1Angle.Min(S2EdgeUtil.GetDistance(aLo, bLo, bHi, bLoCrossHi),
                               S1Angle.Min(S2EdgeUtil.GetDistance(aHi, bLo, bHi, bLoCrossHi),
                                           S1Angle.Min(S2EdgeUtil.GetDistance(bLo, aLo, aHi, aLoCrossHi),
                                                       S2EdgeUtil.GetDistance(bHi, aLo, aHi, aLoCrossHi))));
        }

        /**
   * Return the width and height of this rectangle in latitude-longitude space.
   * Empty rectangles have a negative width and height.
   */

        /**
   * More efficient version of Contains() that accepts a S2LatLng rather than an
   * S2Point.
   */

        public bool Contains(S2LatLng ll)
        {
            // assert (ll.isValid());
            return (_lat.Contains(ll.Lat.Radians) && _lng.Contains(ll.Lng.Radians));
        }

        /**
   * Return true if and only if the given point is contained in the interior of
   * the region (i.e. the region excluding its boundary). The point 'p' does not
   * need to be normalized.
   */

        public bool InteriorContains(S2Point p)
        {
            return InteriorContains(new S2LatLng(p));
        }

        /**
   * More efficient version of InteriorContains() that accepts a S2LatLng rather
   * than an S2Point.
   */

        public bool InteriorContains(S2LatLng ll)
        {
            // assert (ll.isValid());
            return (_lat.InteriorContains(ll.Lat.Radians) && _lng
                                                                 .InteriorContains(ll.Lng.Radians));
        }

        /**
   * Return true if and only if the rectangle contains the given other
   * rectangle.
   */

        public bool Contains(S2LatLngRect other)
        {
            return _lat.Contains(other._lat) && _lng.Contains(other._lng);
        }

        /**
   * Return true if and only if the interior of this rectangle contains all
   * points of the given other rectangle (including its boundary).
   */

        public bool InteriorContains(S2LatLngRect other)
        {
            return (_lat.InteriorContains(other._lat) && _lng
                                                             .InteriorContains(other._lng));
        }

        /** Return true if this rectangle and the given other rectangle have any
  points in common. */

        public bool Intersects(S2LatLngRect other)
        {
            return _lat.Intersects(other._lat) && _lng.Intersects(other._lng);
        }

        /**
   * Returns true if this rectangle intersects the given cell. (This is an exact
   * test and may be fairly expensive, see also MayIntersect below.)
   */

        public bool Intersects(S2Cell cell)
        {
            // First we eliminate the cases where one region completely contains the
            // other. Once these are disposed of, then the regions will intersect
            // if and only if their boundaries intersect.

            if (IsEmpty)
            {
                return false;
            }
            if (Contains(cell.Center))
            {
                return true;
            }
            if (cell.Contains(Center.ToPoint()))
            {
                return true;
            }

            // Quick rejection test (not required for correctness).
            if (!Intersects(cell.RectBound))
            {
                return false;
            }

            // Now check whether the boundaries intersect. Unfortunately, a
            // latitude-longitude rectangle does not have straight edges -- two edges
            // are curved, and at least one of them is concave.

            // Precompute the cell vertices as points and latitude-longitudes.
            var cellV = new S2Point[4];
            var cellLl = new S2LatLng[4];
            for (var i = 0; i < 4; ++i)
            {
                cellV[i] = cell.GetVertex(i); // Must be normalized.
                cellLl[i] = new S2LatLng(cellV[i]);
                if (Contains(cellLl[i]))
                {
                    return true; // Quick acceptance test.
                }
            }

            for (var i = 0; i < 4; ++i)
            {
                var edgeLng = S1Interval.FromPointPair(
                    cellLl[i].Lng.Radians, cellLl[(i + 1) & 3].Lng.Radians);
                if (!_lng.Intersects(edgeLng))
                {
                    continue;
                }

                var a = cellV[i];
                var b = cellV[(i + 1) & 3];
                if (edgeLng.Contains(_lng.Lo))
                {
                    if (IntersectsLngEdge(a, b, _lat, _lng.Lo))
                    {
                        return true;
                    }
                }
                if (edgeLng.Contains(_lng.Hi))
                {
                    if (IntersectsLngEdge(a, b, _lat, _lng.Hi))
                    {
                        return true;
                    }
                }
                if (IntersectsLatEdge(a, b, _lat.Lo, _lng))
                {
                    return true;
                }
                if (IntersectsLatEdge(a, b, _lat.Hi, _lng))
                {
                    return true;
                }
            }
            return false;
        }

        /**
   * Return true if and only if the interior of this rectangle intersects any
   * point (including the boundary) of the given other rectangle.
   */

        public bool InteriorIntersects(S2LatLngRect other)
        {
            return (_lat.InteriorIntersects(other._lat) && _lng
                                                               .InteriorIntersects(other._lng));
        }

        public S2LatLngRect AddPoint(S2Point p)
        {
            return AddPoint(new S2LatLng(p));
        }

        // Increase the size of the bounding rectangle to include the given point.
        // The rectangle is expanded by the minimum amount possible.
        public S2LatLngRect AddPoint(S2LatLng ll)
        {
            // assert (ll.isValid());
            var newLat = _lat.AddPoint(ll.Lat.Radians);
            var newLng = _lng.AddPoint(ll.Lng.Radians);
            return new S2LatLngRect(newLat, newLng);
        }

        /**
   * Return a rectangle that contains all points whose latitude distance from
   * this rectangle is at most margin.Lat, and whose longitude distance from
   * this rectangle is at most margin.Lng. In particular, latitudes are
   * clamped while longitudes are wrapped. Note that any expansion of an empty
   * interval remains empty, and both components of the given margin must be
   * non-negative.
   *
   * NOTE: If you are trying to grow a rectangle by a certain *distance* on the
   * sphere (e.g. 5km), use the ConvolveWithCap() method instead.
   */

        public S2LatLngRect Expanded(S2LatLng margin)
        {
            // assert (margin.Lat.radians() >= 0 && margin.Lng.radians() >= 0);
            if (IsEmpty)
            {
                return this;
            }
            return new S2LatLngRect(_lat.Expanded(margin.Lat.Radians).Intersection(
                FullLat), _lng.Expanded(margin.Lng.Radians));
        }

        /**
   * Return the smallest rectangle containing the union of this rectangle and
   * the given rectangle.
   */

        public S2LatLngRect Union(S2LatLngRect other)
        {
            return new S2LatLngRect(_lat.Union(other._lat), _lng.Union(other._lng));
        }

        /**
   * Return the smallest rectangle containing the intersection of this rectangle
   * and the given rectangle. Note that the region of intersection may consist
   * of two disjoint rectangles, in which case a single rectangle spanning both
   * of them is returned.
   */

        public S2LatLngRect Intersection(S2LatLngRect other)
        {
            var intersectLat = _lat.Intersection(other._lat);
            var intersectLng = _lng.Intersection(other._lng);
            if (intersectLat.IsEmpty || intersectLng.IsEmpty)
            {
                // The lat/lng ranges must either be both empty or both non-empty.
                return Empty;
            }
            return new S2LatLngRect(intersectLat, intersectLng);
        }

        /**
   * Return a rectangle that contains the convolution of this rectangle with a
   * cap of the given angle. This expands the rectangle by a fixed distance (as
   * opposed to growing the rectangle in latitude-longitude space). The returned
   * rectangle includes all points whose minimum distance to the original
   * rectangle is at most the given angle.
   */

        public S2LatLngRect ConvolveWithCap(S1Angle angle)
        {
            // The most straightforward approach is to build a cap centered on each
            // vertex and take the union of all the bounding rectangles (including the
            // original rectangle; this is necessary for very large rectangles).

            // Optimization: convert the angle to a height exactly once.
            var cap = S2Cap.FromAxisAngle(new S2Point(1, 0, 0), angle);

            var r = this;
            for (var k = 0; k < 4; ++k)
            {
                var vertexCap = S2Cap.FromAxisHeight(GetVertex(k).ToPoint(), cap.Height);
                r = r.Union(vertexCap.RectBound);
            }
            return r;
        }

        /** Return the surface area of this rectangle on the unit sphere. */


        /**
   * Return true if the latitude and longitude intervals of the two rectangles
   * are the same up to the given tolerance (see r1interval.h and s1interval.h
   * for details).
   */

        public bool ApproxEquals(S2LatLngRect other, double maxError)
        {
            return (_lat.ApproxEquals(other._lat, maxError) && _lng.ApproxEquals(
                other._lng, maxError));
        }

        public bool ApproxEquals(S2LatLngRect other)
        {
            return ApproxEquals(other, 1e-15);
        }

        // //////////////////////////////////////////////////////////////////////
        // S2Region interface (see {@code S2Region} for details):

        public IS2Region Clone()
        {
            return new S2LatLngRect(Lo, Hi);
        }

        /** The point 'p' does not need to be normalized. */

        public bool Contains(S2Point p)
        {
            return Contains(new S2LatLng(p));
        }

        /**
   * Return true if the edge AB intersects the given edge of constant longitude.
   */

        private static bool IntersectsLngEdge(S2Point a, S2Point b,
                                              R1Interval lat, double lng)
        {
            // Return true if the segment AB intersects the given edge of constant
            // longitude. The nice thing about edges of constant longitude is that
            // they are straight lines on the sphere (geodesics).

            return S2.SimpleCrossing(a, b, S2LatLng.FromRadians(lat.Lo, lng)
                                                   .ToPoint(), S2LatLng.FromRadians(lat.Hi, lng).ToPoint());
        }

        /**
   * Return true if the edge AB intersects the given edge of constant latitude.
   */

        private static bool IntersectsLatEdge(S2Point a, S2Point b, double lat,
                                              S1Interval lng)
        {
            // Return true if the segment AB intersects the given edge of constant
            // latitude. Unfortunately, lines of constant latitude are curves on
            // the sphere. They can intersect a straight edge in 0, 1, or 2 points.
            // assert (S2.isUnitLength(a) && S2.isUnitLength(b));

            // First, compute the normal to the plane AB that points vaguely north.
            var z = S2Point.Normalize(S2.RobustCrossProd(a, b));
            if (z.Z < 0)
            {
                z = -z;
            }

            // Extend this to an orthonormal frame (x,y,z) where x is the direction
            // where the great circle through AB achieves its maximium latitude.
            var y = S2Point.Normalize(S2.RobustCrossProd(z, new S2Point(0, 0, 1)));
            var x = S2Point.CrossProd(y, z);
            // assert (S2.isUnitLength(x) && x.z >= 0);

            // Compute the angle "theta" from the x-axis (in the x-y plane defined
            // above) where the great circle intersects the given line of latitude.
            var sinLat = Math.Sin(lat);
            if (Math.Abs(sinLat) >= x.Z)
            {
                return false; // The great circle does not reach the given latitude.
            }
            // assert (x.z > 0);
            var cosTheta = sinLat/x.Z;
            var sinTheta = Math.Sqrt(1 - cosTheta*cosTheta);
            var theta = Math.Atan2(sinTheta, cosTheta);

            // The candidate intersection points are located +/- theta in the x-y
            // plane. For an intersection to be valid, we need to check that the
            // intersection point is contained in the interior of the edge AB and
            // also that it is contained within the given longitude interval "lng".

            // Compute the range of theta values spanned by the edge AB.
            var abTheta = S1Interval.FromPointPair(Math.Atan2(
                a.DotProd(y), a.DotProd(x)), Math.Atan2(b.DotProd(y), b.DotProd(x)));

            if (abTheta.Contains(theta))
            {
                // Check if the intersection point is also in the given "lng" interval.
                var isect = (x*cosTheta) + (y*sinTheta);
                if (lng.Contains(Math.Atan2(isect.Y, isect.X)))
                {
                    return true;
                }
            }
            if (abTheta.Contains(-theta))
            {
                // Check if the intersection point is also in the given "lng" interval.
                var intersection = (x*cosTheta) - (y*sinTheta);
                if (lng.Contains(Math.Atan2(intersection.Y, intersection.X)))
                {
                    return true;
                }
            }
            return false;
        }

        public override String ToString()
        {
            return "[Lo=" + Lo + ", Hi=" + Hi + "]";
        }
    }
}