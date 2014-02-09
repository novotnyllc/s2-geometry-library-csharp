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
        private readonly R1Interval lat;
        private readonly S1Interval lng;

        /**
   * Construct a rectangle from minimum and maximum latitudes and longitudes. If
   * lo.Lng > hi.Lng, the rectangle spans the 180 degree longitude line.
   */

        public S2LatLngRect(S2LatLng lo, S2LatLng hi)
        {
            lat = new R1Interval(lo.lat().Radians, hi.lat().Radians);
            lng = new S1Interval(lo.lng().Radians, hi.lng().Radians);
            // assert (isValid());
        }

        /** Construct a rectangle from latitude and longitude intervals. */

        public S2LatLngRect(R1Interval lat, S1Interval lng)
        {
            this.lat = lat;
            this.lng = lng;
            // assert (isValid());
        }

        public R1Interval Lat
        {
            get { return lat; }
        }

        public S1Interval Lng
        {
            get { return lng; }
        }

        public bool Equals(S2LatLngRect other)
        {
            return Equals(lat, other.lat) && Equals(lng, other.lng);
        }

        public S2Cap CapBound
        {
            get
            {
                // We consider two possible bounding caps, one whose axis passes
                // through the center of the lat-long rectangle and one whose axis
                // is the north or south pole. We return the smaller of the two caps.

                if (isEmpty())
                {
                    return S2Cap.empty();
                }

                double poleZ, poleAngle;
                if (lat.Lo + lat.Hi < 0)
                {
                    // South pole axis yields smaller cap.
                    poleZ = -1;
                    poleAngle = S2.PiOver2 + lat.Hi;
                }
                else
                {
                    poleZ = 1;
                    poleAngle = S2.PiOver2 - lat.Lo;
                }
                var poleCap = S2Cap.fromAxisAngle(new S2Point(0, 0, poleZ), S1Angle
                                                                                .FromRadians(poleAngle));

                // For bounding rectangles that span 180 degrees or less in longitude, the
                // maximum cap size is achieved at one of the rectangle vertices. For
                // rectangles that are larger than 180 degrees, we punt and always return a
                // bounding cap centered at one of the two poles.
                var lngSpan = lng.Hi - lng.Lo;
                if (Math.IEEERemainder(lngSpan, 2*S2.Pi) >= 0)
                {
                    if (lngSpan < 2*S2.Pi)
                    {
                        var midCap = S2Cap.fromAxisAngle(getCenter().toPoint(), S1Angle
                                                                                    .FromRadians(0));
                        for (var k = 0; k < 4; ++k)
                        {
                            midCap = midCap.addPoint(getVertex(k).toPoint());
                        }
                        if (midCap.height() < poleCap.height())
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
            return contains(cell.RectBound);
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
            return intersects(cell.RectBound);
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
                return ((lat != null ? lat.GetHashCode() : 0)*397) ^ (lng != null ? lng.GetHashCode() : 0);
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

        public static S2LatLngRect empty()
        {
            return new S2LatLngRect(R1Interval.Empty, S1Interval.Empty);
        }

        /** The canonical full rectangle. */

        public static S2LatLngRect full()
        {
            return new S2LatLngRect(fullLat(), fullLng());
        }

        /** The full allowable range of latitudes. */

        public static R1Interval fullLat()
        {
            return new R1Interval(-S2.PiOver2, S2.PiOver2);
        }

        /**
   * The full allowable range of longitudes.
   */

        public static S1Interval fullLng()
        {
            return S1Interval.Full;
        }

        /**
   * Construct a rectangle from a center point (in lat-lng space) and size in
   * each dimension. If size.Lng is greater than 360 degrees it is clamped,
   * and latitudes greater than +/- 90 degrees are also clamped. So for example,
   * FromCenterSize((80,170),(20,20)) -> (lo=(60,150),hi=(90,-170)).
   */

        public static S2LatLngRect fromCenterSize(S2LatLng center, S2LatLng size)
        {
            return fromPoint(center).expanded(size.mul(0.5));
        }

        /** Convenience method to construct a rectangle containing a single point. */

        public static S2LatLngRect fromPoint(S2LatLng p)
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

        public static S2LatLngRect fromPointPair(S2LatLng p1, S2LatLng p2)
        {
            // assert (p1.isValid() && p2.isValid());
            return new S2LatLngRect(R1Interval.FromPointPair(p1.lat().Radians, p2
                                                                                   .lat().Radians), S1Interval.FromPointPair(p1.lng().Radians, p2.lng().Radians));
        }

        /**
   * Return a latitude-longitude rectangle that contains the edge from "a" to
   * "b". Both points must be unit-length. Note that the bounding rectangle of
   * an edge can be larger than the bounding rectangle of its endpoints.
   */

        public static S2LatLngRect fromEdge(S2Point a, S2Point b)
        {
            // assert (S2.isUnitLength(a) && S2.isUnitLength(b));
            var r = fromPointPair(new S2LatLng(a), new S2LatLng(b));

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

        public bool isValid()
        {
            // The lat/lng ranges must either be both empty or both non-empty.
            return (Math.Abs(lat.Lo) <= S2.PiOver2 && Math.Abs(lat.Hi) <= S2.PiOver2
                    && lng.IsValid && lat.IsEmpty == lng.IsEmpty);
        }

        // Accessor methods.
        public S1Angle latLo()
        {
            return S1Angle.FromRadians(lat.Lo);
        }

        public S1Angle latHi()
        {
            return S1Angle.FromRadians(lat.Hi);
        }

        public S1Angle lngLo()
        {
            return S1Angle.FromRadians(lng.Lo);
        }

        public S1Angle lngHi()
        {
            return S1Angle.FromRadians(lng.Hi);
        }

        public S2LatLng lo()
        {
            return new S2LatLng(latLo(), lngLo());
        }

        public S2LatLng hi()
        {
            return new S2LatLng(latHi(), lngHi());
        }

        /**
   * Return true if the rectangle is empty, i.e. it contains no points at all.
   */

        public bool isEmpty()
        {
            return lat.IsEmpty;
        }

        // Return true if the rectangle is full, i.e. it contains all points.
        public bool isFull()
        {
            return lat.Equals(fullLat()) && lng.IsFull;
        }

        /**
   * Return true if lng_.lo() > lng_.hi(), i.e. the rectangle crosses the 180
   * degree latitude line.
   */

        public bool isInverted()
        {
            return lng.IsInverted;
        }

        /** Return the k-th vertex of the rectangle (k = 0,1,2,3) in CCW order. */

        public S2LatLng getVertex(int k)
        {
            // Return the points in CCW order (SW, SE, NE, NW).
            switch (k)
            {
                case 0:
                    return S2LatLng.fromRadians(lat.Lo, lng.Lo);
                case 1:
                    return S2LatLng.fromRadians(lat.Lo, lng.Hi);
                case 2:
                    return S2LatLng.fromRadians(lat.Hi, lng.Hi);
                case 3:
                    return S2LatLng.fromRadians(lat.Hi, lng.Lo);
                default:
                    throw new ArgumentException("Invalid vertex index.");
            }
        }

        /**
   * Return the center of the rectangle in latitude-longitude space (in general
   * this is not the center of the region on the sphere).
   */

        public S2LatLng getCenter()
        {
            return S2LatLng.fromRadians(lat.Center, lng.Center);
        }

        /**
   * Return the minimum distance (measured along the surface of the sphere)
   * from a given point to the rectangle (both its boundary and its interior).
   * The latLng must be valid.
   */

        public S1Angle getDistance(S2LatLng p)
        {
            // The algorithm here is the same as in getDistance(S2LagLngRect), only
            // with simplified calculations.
            var a = this;

            Preconditions.CheckState(!a.isEmpty());
            Preconditions.CheckArgument(p.isValid());

            if (a.Lng.Contains(p.lng().Radians))
            {
                return S1Angle.FromRadians(Math.Max(0.0, Math.Max(p.lat().Radians - a.Lat.Hi,
                                                              a.Lat.Lo - p.lat().Radians)));
            }

            var interval = new S1Interval(a.Lng.Hi, a.Lng.Complement.Center);
            var aLng = a.Lng.Lo;
            if (interval.Contains(p.lng().Radians))
            {
                aLng = a.Lng.Hi;
            }

            var lo = S2LatLng.fromRadians(a.Lat.Lo, aLng).toPoint();
            var hi = S2LatLng.fromRadians(a.Lat.Hi, aLng).toPoint();
            var loCrossHi =
                S2LatLng.fromRadians(0, aLng - S2.PiOver2).normalized().toPoint();
            return S2EdgeUtil.getDistance(p.toPoint(), lo, hi, loCrossHi);
        }

        /**
   * Return the minimum distance (measured along the surface of the sphere) to
   * the given S2LatLngRect. Both S2LatLngRects must be non-empty.
   */

        public S1Angle getDistance(S2LatLngRect other)
        {
            var a = this;
            var b = other;

            Preconditions.CheckState(!a.isEmpty());
            Preconditions.CheckArgument(!b.isEmpty());

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
                    lo = b.latHi();
                    hi = a.latLo();
                }
                else
                {
                    lo = a.latHi();
                    hi = b.latLo();
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
                aLng = a.lngLo();
                bLng = b.lngHi();
            }
            else
            {
                aLng = a.lngHi();
                bLng = b.lngLo();
            }

            // The shortest distance between the two longitudinal segments will include
            // at least one segment endpoint. We could probably narrow this down further
            // to a single point-edge distance by comparing the relative latitudes of the
            // endpoints, but for the sake of clarity, we'll do all four point-edge
            // distance tests.
            var aLo = new S2LatLng(a.latLo(), aLng).toPoint();
            var aHi = new S2LatLng(a.latHi(), aLng).toPoint();
            var aLoCrossHi =
                S2LatLng.fromRadians(0, aLng.Radians - S2.PiOver2).normalized().toPoint();
            var bLo = new S2LatLng(b.latLo(), bLng).toPoint();
            var bHi = new S2LatLng(b.latHi(), bLng).toPoint();
            var bLoCrossHi =
                S2LatLng.fromRadians(0, bLng.Radians - S2.PiOver2).normalized().toPoint();

            return S1Angle.Min(S2EdgeUtil.getDistance(aLo, bLo, bHi, bLoCrossHi),
                               S1Angle.Min(S2EdgeUtil.getDistance(aHi, bLo, bHi, bLoCrossHi),
                                           S1Angle.Min(S2EdgeUtil.getDistance(bLo, aLo, aHi, aLoCrossHi),
                                                       S2EdgeUtil.getDistance(bHi, aLo, aHi, aLoCrossHi))));
        }

        /**
   * Return the width and height of this rectangle in latitude-longitude space.
   * Empty rectangles have a negative width and height.
   */

        public S2LatLng getSize()
        {
            return S2LatLng.fromRadians(lat.Length, lng.Length);
        }

        /**
   * More efficient version of Contains() that accepts a S2LatLng rather than an
   * S2Point.
   */

        public bool contains(S2LatLng ll)
        {
            // assert (ll.isValid());
            return (lat.Contains(ll.lat().Radians) && lng.Contains(ll.lng().Radians));
        }

        /**
   * Return true if and only if the given point is contained in the interior of
   * the region (i.e. the region excluding its boundary). The point 'p' does not
   * need to be normalized.
   */

        public bool interiorContains(S2Point p)
        {
            return interiorContains(new S2LatLng(p));
        }

        /**
   * More efficient version of InteriorContains() that accepts a S2LatLng rather
   * than an S2Point.
   */

        public bool interiorContains(S2LatLng ll)
        {
            // assert (ll.isValid());
            return (lat.InteriorContains(ll.lat().Radians) && lng
                                                                    .InteriorContains(ll.lng().Radians));
        }

        /**
   * Return true if and only if the rectangle contains the given other
   * rectangle.
   */

        public bool contains(S2LatLngRect other)
        {
            return lat.Contains(other.lat) && lng.Contains(other.lng);
        }

        /**
   * Return true if and only if the interior of this rectangle contains all
   * points of the given other rectangle (including its boundary).
   */

        public bool interiorContains(S2LatLngRect other)
        {
            return (lat.InteriorContains(other.lat) && lng
                                                           .InteriorContains(other.lng));
        }

        /** Return true if this rectangle and the given other rectangle have any
  points in common. */

        public bool intersects(S2LatLngRect other)
        {
            return lat.Intersects(other.lat) && lng.Intersects(other.lng);
        }

        /**
   * Returns true if this rectangle intersects the given cell. (This is an exact
   * test and may be fairly expensive, see also MayIntersect below.)
   */

        public bool intersects(S2Cell cell)
        {
            // First we eliminate the cases where one region completely contains the
            // other. Once these are disposed of, then the regions will intersect
            // if and only if their boundaries intersect.

            if (isEmpty())
            {
                return false;
            }
            if (contains(cell.getCenter()))
            {
                return true;
            }
            if (cell.contains(getCenter().toPoint()))
            {
                return true;
            }

            // Quick rejection test (not required for correctness).
            if (!intersects(cell.RectBound))
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
                cellV[i] = cell.getVertex(i); // Must be normalized.
                cellLl[i] = new S2LatLng(cellV[i]);
                if (contains(cellLl[i]))
                {
                    return true; // Quick acceptance test.
                }
            }

            for (var i = 0; i < 4; ++i)
            {
                var edgeLng = S1Interval.FromPointPair(
                    cellLl[i].lng().Radians, cellLl[(i + 1) & 3].lng().Radians);
                if (!lng.Intersects(edgeLng))
                {
                    continue;
                }

                var a = cellV[i];
                var b = cellV[(i + 1) & 3];
                if (edgeLng.Contains(lng.Lo))
                {
                    if (intersectsLngEdge(a, b, lat, lng.Lo))
                    {
                        return true;
                    }
                }
                if (edgeLng.Contains(lng.Hi))
                {
                    if (intersectsLngEdge(a, b, lat, lng.Hi))
                    {
                        return true;
                    }
                }
                if (intersectsLatEdge(a, b, lat.Lo, lng))
                {
                    return true;
                }
                if (intersectsLatEdge(a, b, lat.Hi, lng))
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

        public bool interiorIntersects(S2LatLngRect other)
        {
            return (lat.InteriorIntersects(other.lat) && lng
                                                             .InteriorIntersects(other.lng));
        }

        public S2LatLngRect addPoint(S2Point p)
        {
            return addPoint(new S2LatLng(p));
        }

        // Increase the size of the bounding rectangle to include the given point.
        // The rectangle is expanded by the minimum amount possible.
        public S2LatLngRect addPoint(S2LatLng ll)
        {
            // assert (ll.isValid());
            var newLat = lat.AddPoint(ll.lat().Radians);
            var newLng = lng.AddPoint(ll.lng().Radians);
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

        public S2LatLngRect expanded(S2LatLng margin)
        {
            // assert (margin.Lat.radians() >= 0 && margin.Lng.radians() >= 0);
            if (isEmpty())
            {
                return this;
            }
            return new S2LatLngRect(lat.Expanded(margin.lat().Radians).Intersection(
                fullLat()), lng.Expanded(margin.lng().Radians));
        }

        /**
   * Return the smallest rectangle containing the union of this rectangle and
   * the given rectangle.
   */

        public S2LatLngRect union(S2LatLngRect other)
        {
            return new S2LatLngRect(lat.Union(other.lat), lng.Union(other.lng));
        }

        /**
   * Return the smallest rectangle containing the intersection of this rectangle
   * and the given rectangle. Note that the region of intersection may consist
   * of two disjoint rectangles, in which case a single rectangle spanning both
   * of them is returned.
   */

        public S2LatLngRect intersection(S2LatLngRect other)
        {
            var intersectLat = lat.Intersection(other.lat);
            var intersectLng = lng.Intersection(other.lng);
            if (intersectLat.IsEmpty || intersectLng.IsEmpty)
            {
                // The lat/lng ranges must either be both empty or both non-empty.
                return empty();
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

        public S2LatLngRect convolveWithCap(S1Angle angle)
        {
            // The most straightforward approach is to build a cap centered on each
            // vertex and take the union of all the bounding rectangles (including the
            // original rectangle; this is necessary for very large rectangles).

            // Optimization: convert the angle to a height exactly once.
            var cap = S2Cap.fromAxisAngle(new S2Point(1, 0, 0), angle);

            var r = this;
            for (var k = 0; k < 4; ++k)
            {
                var vertexCap = S2Cap.fromAxisHeight(getVertex(k).toPoint(), cap
                                                                                 .height());
                r = r.union(vertexCap.RectBound);
            }
            return r;
        }

        /** Return the surface area of this rectangle on the unit sphere. */

        public double area()
        {
            if (isEmpty())
            {
                return 0;
            }

            // This is the size difference of the two spherical caps, multiplied by
            // the longitude ratio.
            return Lng.Length*Math.Abs(Math.Sin(latHi().Radians) - Math.Sin(latLo().Radians));
        }


        /**
   * Return true if the latitude and longitude intervals of the two rectangles
   * are the same up to the given tolerance (see r1interval.h and s1interval.h
   * for details).
   */

        public bool approxEquals(S2LatLngRect other, double maxError)
        {
            return (lat.ApproxEquals(other.lat, maxError) && lng.ApproxEquals(
                other.lng, maxError));
        }

        public bool approxEquals(S2LatLngRect other)
        {
            return approxEquals(other, 1e-15);
        }

        // //////////////////////////////////////////////////////////////////////
        // S2Region interface (see {@code S2Region} for details):

        public IS2Region clone()
        {
            return new S2LatLngRect(lo(), hi());
        }

        /** The point 'p' does not need to be normalized. */

        public bool contains(S2Point p)
        {
            return contains(new S2LatLng(p));
        }

        /**
   * Return true if the edge AB intersects the given edge of constant longitude.
   */

        private static bool intersectsLngEdge(S2Point a, S2Point b,
                                              R1Interval lat, double lng)
        {
            // Return true if the segment AB intersects the given edge of constant
            // longitude. The nice thing about edges of constant longitude is that
            // they are straight lines on the sphere (geodesics).

            return S2.SimpleCrossing(a, b, S2LatLng.fromRadians(lat.Lo, lng)
                                                   .toPoint(), S2LatLng.fromRadians(lat.Hi, lng).toPoint());
        }

        /**
   * Return true if the edge AB intersects the given edge of constant latitude.
   */

        private static bool intersectsLatEdge(S2Point a, S2Point b, double lat,
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
                var isect = (x * cosTheta) + (y * sinTheta);
                if (lng.Contains(Math.Atan2(isect.Y, isect.X)))
                {
                    return true;
                }
            }
            if (abTheta.Contains(-theta))
            {
                // Check if the intersection point is also in the given "lng" interval.
                var intersection = (x * cosTheta) - (y * sinTheta);
                if (lng.Contains(Math.Atan2(intersection.Y, intersection.X)))
                {
                    return true;
                }
            }
            return false;
        }

        public override String ToString()
        {
            return "[Lo=" + lo() + ", Hi=" + hi() + "]";
        }
    }
}