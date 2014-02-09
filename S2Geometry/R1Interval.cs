using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Google.Common.Geometry
{
    /// <summary>
    ///     An R1Interval represents a closed, bounded interval on the real line. It is
    ///     capable of representing the empty interval (containing no points) and
    ///     zero-length intervals (containing a single point).
    /// </summary>
    public sealed class R1Interval : IEquatable<R1Interval>
    {
        private readonly double _hi;
        private readonly double _lo;

        public R1Interval(double lo, double hi)
        {
            _lo = lo;
            _hi = hi;
        }

        public bool Equals(R1Interval other)
        {
            if (ReferenceEquals(null, other)) return false;
            if (ReferenceEquals(this, other)) return true;
            return (_hi.Equals(other._hi) && _lo.Equals(other._lo)) || (isEmpty() && other.isEmpty());
        }

        public override bool Equals(object obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            if (ReferenceEquals(this, obj)) return true;
            if (obj.GetType() != GetType()) return false;
            return Equals((R1Interval)obj);
        }

        public override int GetHashCode()
        {
            unchecked
            {
                return (_hi.GetHashCode()*397) ^ _lo.GetHashCode();
            }
        }

        public static bool operator ==(R1Interval left, R1Interval right)
        {
            return Equals(left, right);
        }

        public static bool operator !=(R1Interval left, R1Interval right)
        {
            return !Equals(left, right);
        }

        /**
   * Returns an empty interval. (Any interval where lo > hi is considered
   * empty.)
   */

        public static R1Interval empty()
        {
            return new R1Interval(1, 0);
        }

        /**
   * Convenience method to construct an interval containing a single point.
   */

        public static R1Interval fromPoint(double p)
        {
            return new R1Interval(p, p);
        }

        /**
   * Convenience method to construct the minimal interval containing the two
   * given points. This is equivalent to starting with an empty interval and
   * calling AddPoint() twice, but it is more efficient.
   */

        public static R1Interval fromPointPair(double p1, double p2)
        {
            if (p1 <= p2)
            {
                return new R1Interval(p1, p2);
            }
            else
            {
                return new R1Interval(p2, p1);
            }
        }

        public double lo()
        {
            return _lo;
        }

        public double hi()
        {
            return _hi;
        }

        /**
   * Return true if the interval is empty, i.e. it contains no points.
   */

        public bool isEmpty()
        {
            return lo() > hi();
        }

        /**
   * Return the center of the interval. For empty intervals, the result is
   * arbitrary.
   */

        public double getCenter()
        {
            return 0.5*(lo() + hi());
        }

        /**
   * Return the length of the interval. The length of an empty interval is
   * negative.
   */

        public double getLength()
        {
            return hi() - lo();
        }

        public bool contains(double p)
        {
            return p >= lo() && p <= hi();
        }

        public bool interiorContains(double p)
        {
            return p > lo() && p < hi();
        }

        /** Return true if this interval contains the interval 'y'. */

        public bool contains(R1Interval y)
        {
            if (y.isEmpty())
            {
                return true;
            }
            return y.lo() >= lo() && y.hi() <= hi();
        }

        /**
   * Return true if the interior of this interval contains the entire interval
   * 'y' (including its boundary).
   */

        public bool interiorContains(R1Interval y)
        {
            if (y.isEmpty())
            {
                return true;
            }
            return y.lo() > lo() && y.hi() < hi();
        }

        /**
   * Return true if this interval intersects the given interval, i.e. if they
   * have any points in common.
   */

        public bool intersects(R1Interval y)
        {
            if (lo() <= y.lo())
            {
                return y.lo() <= hi() && y.lo() <= y.hi();
            }
            else
            {
                return lo() <= y.hi() && lo() <= hi();
            }
        }

        /**
   * Return true if the interior of this interval intersects any point of the
   * given interval (including its boundary).
   */

        public bool interiorIntersects(R1Interval y)
        {
            return y.lo() < hi() && lo() < y.hi() && lo() < hi() && y.lo() <= y.hi();
        }

        /** Expand the interval so that it contains the given point "p". */

        public R1Interval addPoint(double p)
        {
            if (isEmpty())
            {
                return fromPoint(p);
            }
            else if (p < lo())
            {
                return new R1Interval(p, hi());
            }
            else if (p > hi())
            {
                return new R1Interval(lo(), p);
            }
            else
            {
                return new R1Interval(lo(), hi());
            }
        }

        /**
   * Return an interval that contains all points with a distance "radius" of a
   * point in this interval. Note that the expansion of an empty interval is
   * always empty.
   */

        public R1Interval expanded(double radius)
        {
            // assert (radius >= 0);
            if (isEmpty())
            {
                return this;
            }
            return new R1Interval(lo() - radius, hi() + radius);
        }

        /**
   * Return the smallest interval that contains this interval and the given
   * interval "y".
   */

        public R1Interval union(R1Interval y)
        {
            if (isEmpty())
            {
                return y;
            }
            if (y.isEmpty())
            {
                return this;
            }
            return new R1Interval(Math.Min(lo(), y.lo()), Math.Max(hi(), y.hi()));
        }

        /**
   * Return the intersection of this interval with the given interval. Empty
   * intervals do not need to be special-cased.
   */

        public R1Interval intersection(R1Interval y)
        {
            return new R1Interval(Math.Max(lo(), y.lo()), Math.Min(hi(), y.hi()));
        }

        public bool approxEquals(R1Interval y)
        {
            return approxEquals(y, 1e-15);
        }

        /**
   * Return true if length of the symmetric difference between the two intervals
   * is at most the given tolerance.
   *
   */

        public bool approxEquals(R1Interval y, double maxError)
        {
            if (isEmpty())
            {
                return y.getLength() <= maxError;
            }
            if (y.isEmpty())
            {
                return getLength() <= maxError;
            }
            return Math.Abs(y.lo() - lo()) + Math.Abs(y.hi() - hi()) <= maxError;
        }

        public override string ToString()
        {
            return "[" + lo() + ", " + hi() + "]";
        }
    }
}