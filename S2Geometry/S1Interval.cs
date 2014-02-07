using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Google.Common.Geometry
{
    public sealed class S1Interval : IEquatable<S1Interval>
    {
        public bool Equals(S1Interval other)
        {
            if (ReferenceEquals(null, other)) return false;
            if (ReferenceEquals(this, other)) return true;
            return _lo.Equals(other._lo) && _hi.Equals(other._hi);
        }

        public override bool Equals(object obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            if (ReferenceEquals(this, obj)) return true;
            return obj is S1Interval && Equals((S1Interval)obj);
        }

        public override int GetHashCode()
        {
            unchecked
            {
                return (_lo.GetHashCode()*397) ^ _hi.GetHashCode();
            }
        }

        public static bool operator ==(S1Interval left, S1Interval right)
        {
            return Equals(left, right);
        }

        public static bool operator !=(S1Interval left, S1Interval right)
        {
            return !Equals(left, right);
        }

        private readonly double _lo;
    private readonly double _hi;

  /**
   * Both endpoints must be in the range -Pi to Pi inclusive. The value -Pi is
   * converted internally to Pi except for the Full() and Empty() intervals.
   */
  public S1Interval(double lo, double hi) : this(lo, hi, false){
    
  }

  /**
   * Copy constructor. Assumes that the given interval is valid.
   *
   * TODO(dbeaumont): Make this class immutable and remove this method.
   */
  public S1Interval(S1Interval interval) {
    this._lo = interval._lo;
    this._hi = interval._hi;
  }

  /**
   * Internal constructor that assumes that both arguments are in the correct
   * range, i.e. normalization from -Pi to Pi is already done.
   */
  private S1Interval(double lo, double hi, bool @checked) {
    double newLo = lo;
    double newHi = hi;
    if (!@checked) {
      if (lo == -S2.M_PI && hi != S2.M_PI) {
        newLo = S2.M_PI;
      }
      if (hi == -S2.M_PI && lo != S2.M_PI) {
        newHi = S2.M_PI;
      }
    }
    this._lo = newLo;
    this._hi = newHi;
  }

  public static S1Interval empty() {
    return new S1Interval(S2.M_PI, -S2.M_PI, true);
  }

  public static S1Interval full() {
    return new S1Interval(-S2.M_PI, S2.M_PI, true);
  }

  /** Convenience method to construct an interval containing a single point. */
  public static S1Interval fromPoint(double p) {
    if (p == -S2.M_PI) {
      p = S2.M_PI;
    }
    return new S1Interval(p, p, true);
  }

  /**
   * Convenience method to construct the minimal interval containing the two
   * given points. This is equivalent to starting with an empty interval and
   * calling AddPoint() twice, but it is more efficient.
   */
  public static S1Interval fromPointPair(double p1, double p2) {
    // assert (Math.Abs(p1) <= S2.M_PI && Math.Abs(p2) <= S2.M_PI);
    if (p1 == -S2.M_PI) {
      p1 = S2.M_PI;
    }
    if (p2 == -S2.M_PI) {
      p2 = S2.M_PI;
    }
    if (positiveDistance(p1, p2) <= S2.M_PI) {
      return new S1Interval(p1, p2, true);
    } else {
      return new S1Interval(p2, p1, true);
    }
  }

  public double lo() {
    return _lo;
  }

  public double hi() {
    return _hi;
  }

  /**
   * An interval is valid if neither bound exceeds Pi in absolute value, and the
   * value -Pi appears only in the Empty() and Full() intervals.
   */
  public bool isValid() {
    return (Math.Abs(lo()) <= S2.M_PI && Math.Abs(hi()) <= S2.M_PI
        && !(lo() == -S2.M_PI && hi() != S2.M_PI) && !(hi() == -S2.M_PI && lo() != S2.M_PI));
  }

  /** Return true if the interval contains all points on the unit circle. */
  public bool isFull() {
    return hi() - lo() == 2 * S2.M_PI;
  }


  /** Return true if the interval is empty, i.e. it contains no points. */
  public bool isEmpty() {
    return lo() - hi() == 2 * S2.M_PI;
  }


  /* Return true if lo() > hi(). (This is true for empty intervals.) */
  public bool isInverted() {
    return lo() > hi();
  }

  /**
   * Return the midpoint of the interval. For full and empty intervals, the
   * result is arbitrary.
   */
  public double getCenter() {
    double center = 0.5 * (lo() + hi());
    if (!isInverted()) {
      return center;
    }
    // Return the center in the range (-Pi, Pi].
    return (center <= 0) ? (center + S2.M_PI) : (center - S2.M_PI);
  }

  /**
   * Return the length of the interval. The length of an empty interval is
   * negative.
   */
  public double getLength() {
    double length = hi() - lo();
    if (length >= 0) {
      return length;
    }
    length += 2 * S2.M_PI;
    // Empty intervals have a negative length.
    return (length > 0) ? length : -1;
  }

  /**
   * Return the complement of the interior of the interval. An interval and its
   * complement have the same boundary but do not share any interior values. The
   * complement operator is not a bijection, since the complement of a singleton
   * interval (containing a single value) is the same as the complement of an
   * empty interval.
   */
  public S1Interval complement() {
    if (lo() == hi()) {
      return full(); // Singleton.
    }
    return new S1Interval(hi(), lo(), true); // Handles
    // empty and
    // full.
  }

  /** Return true if the interval (which is closed) contains the point 'p'. */
  public bool contains(double p) {
    // Works for empty, full, and singleton intervals.
    // assert (Math.Abs(p) <= S2.M_PI);
    if (p == -S2.M_PI) {
      p = S2.M_PI;
    }
    return fastContains(p);
  }

  /**
   * Return true if the interval (which is closed) contains the point 'p'. Skips
   * the normalization of 'p' from -Pi to Pi.
   *
   */
  public bool fastContains(double p) {
    if (isInverted()) {
      return (p >= lo() || p <= hi()) && !isEmpty();
    } else {
      return p >= lo() && p <= hi();
    }
  }

  /** Return true if the interior of the interval contains the point 'p'. */
  public bool interiorContains(double p) {
    // Works for empty, full, and singleton intervals.
    // assert (Math.Abs(p) <= S2.M_PI);
    if (p == -S2.M_PI) {
      p = S2.M_PI;
    }

    if (isInverted()) {
      return p > lo() || p < hi();
    } else {
      return (p > lo() && p < hi()) || isFull();
    }
  }

  /**
   * Return true if the interval contains the given interval 'y'. Works for
   * empty, full, and singleton intervals.
   */
  public bool contains(S1Interval y) {
    // It might be helpful to compare the structure of these tests to
    // the simpler Contains(double) method above.

    if (isInverted()) {
      if (y.isInverted()) {
        return y.lo() >= lo() && y.hi() <= hi();
      }
      return (y.lo() >= lo() || y.hi() <= hi()) && !isEmpty();
    } else {
      if (y.isInverted()) {
        return isFull() || y.isEmpty();
      }
      return y.lo() >= lo() && y.hi() <= hi();
    }
  }

  /**
   * Returns true if the interior of this interval contains the entire interval
   * 'y'. Note that x.InteriorContains(x) is true only when x is the empty or
   * full interval, and x.InteriorContains(S1Interval(p,p)) is equivalent to
   * x.InteriorContains(p).
   */
  public bool interiorContains(S1Interval y) {
    if (isInverted()) {
      if (!y.isInverted()) {
        return y.lo() > lo() || y.hi() < hi();
      }
      return (y.lo() > lo() && y.hi() < hi()) || y.isEmpty();
    } else {
      if (y.isInverted()) {
        return isFull() || y.isEmpty();
      }
      return (y.lo() > lo() && y.hi() < hi()) || isFull();
    }
  }

  /**
   * Return true if the two intervals contain any points in common. Note that
   * the point +/-Pi has two representations, so the intervals [-Pi,-3] and
   * [2,Pi] intersect, for example.
   */
  public bool intersects(S1Interval y) {
    if (isEmpty() || y.isEmpty()) {
      return false;
    }
    if (isInverted()) {
      // Every non-empty inverted interval contains Pi.
      return y.isInverted() || y.lo() <= hi() || y.hi() >= lo();
    } else {
      if (y.isInverted()) {
        return y.lo() <= hi() || y.hi() >= lo();
      }
      return y.lo() <= hi() && y.hi() >= lo();
    }
  }

  /**
   * Return true if the interior of this interval contains any point of the
   * interval 'y' (including its boundary). Works for empty, full, and singleton
   * intervals.
   */
  public bool interiorIntersects(S1Interval y) {
    if (isEmpty() || y.isEmpty() || lo() == hi()) {
      return false;
    }
    if (isInverted()) {
      return y.isInverted() || y.lo() < hi() || y.hi() > lo();
    } else {
      if (y.isInverted()) {
        return y.lo() < hi() || y.hi() > lo();
      }
      return (y.lo() < hi() && y.hi() > lo()) || isFull();
    }
  }

  /**
   * Expand the interval by the minimum amount necessary so that it contains the
   * given point "p" (an angle in the range [-Pi, Pi]).
   */
  public S1Interval addPoint(double p) {
    // assert (Math.Abs(p) <= S2.M_PI);
    if (p == -S2.M_PI) {
      p = S2.M_PI;
    }

    if (fastContains(p)) {
      return new S1Interval(this);
    }

    if (isEmpty()) {
      return S1Interval.fromPoint(p);
    } else {
      // Compute distance from p to each endpoint.
      double dlo = positiveDistance(p, lo());
      double dhi = positiveDistance(hi(), p);
      if (dlo < dhi) {
        return new S1Interval(p, hi());
      } else {
        return new S1Interval(lo(), p);
      }
      // Adding a point can never turn a non-full interval into a full one.
    }
  }

  /**
   * Return an interval that contains all points within a distance "radius" of
   * a point in this interval. Note that the expansion of an empty interval is
   * always empty. The radius must be non-negative.
   */
  public S1Interval expanded(double radius) {
    // assert (radius >= 0);
    if (isEmpty()) {
      return this;
    }

    // Check whether this interval will be full after expansion, allowing
    // for a 1-bit rounding error when computing each endpoint.
    if (getLength() + 2 * radius >= 2 * S2.M_PI - 1e-15) {
      return full();
    }

    // NOTE(dbeaumont): Should this remainder be 2 * M_PI or just M_PI ??
    double lo = Math.IEEERemainder(this.lo() - radius, 2 * S2.M_PI);
    double hi = Math.IEEERemainder(this.hi() + radius, 2 * S2.M_PI);
    if (lo == -S2.M_PI) {
      lo = S2.M_PI;
    }
    return new S1Interval(lo, hi);
  }

  /**
   * Return the smallest interval that contains this interval and the given
   * interval "y".
   */
  public S1Interval union(S1Interval y) {
    // The y.is_full() case is handled correctly in all cases by the code
    // below, but can follow three separate code paths depending on whether
    // this interval is inverted, is non-inverted but contains Pi, or neither.

    if (y.isEmpty()) {
      return this;
    }
    if (fastContains(y.lo())) {
      if (fastContains(y.hi())) {
        // Either this interval contains y, or the union of the two
        // intervals is the Full() interval.
        if (contains(y)) {
          return this; // is_full() code path
        }
        return full();
      }
      return new S1Interval(lo(), y.hi(), true);
    }
    if (fastContains(y.hi())) {
      return new S1Interval(y.lo(), hi(), true);
    }

    // This interval contains neither endpoint of y. This means that either y
    // contains all of this interval, or the two intervals are disjoint.
    if (isEmpty() || y.fastContains(lo())) {
      return y;
    }

    // Check which pair of endpoints are closer together.
    double dlo = positiveDistance(y.hi(), lo());
    double dhi = positiveDistance(hi(), y.lo());
    if (dlo < dhi) {
      return new S1Interval(y.lo(), hi(), true);
    } else {
      return new S1Interval(lo(), y.hi(), true);
    }
  }

  /**
   * Return the smallest interval that contains the intersection of this
   * interval with "y". Note that the region of intersection may consist of two
   * disjoint intervals.
   */
  public S1Interval intersection(S1Interval y) {
    // The y.is_full() case is handled correctly in all cases by the code
    // below, but can follow three separate code paths depending on whether
    // this interval is inverted, is non-inverted but contains Pi, or neither.

    if (y.isEmpty()) {
      return empty();
    }
    if (fastContains(y.lo())) {
      if (fastContains(y.hi())) {
        // Either this interval contains y, or the region of intersection
        // consists of two disjoint subintervals. In either case, we want
        // to return the shorter of the two original intervals.
        if (y.getLength() < getLength()) {
          return y; // is_full() code path
        }
        return this;
      }
      return new S1Interval(y.lo(), hi(), true);
    }
    if (fastContains(y.hi())) {
      return new S1Interval(lo(), y.hi(), true);
    }

    // This interval contains neither endpoint of y. This means that either y
    // contains all of this interval, or the two intervals are disjoint.

    if (y.fastContains(lo())) {
      return this; // is_empty() okay here
    }
    // assert (!intersects(y));
    return empty();
  }

  /**
   * Return true if the length of the symmetric difference between the two
   * intervals is at most the given tolerance.
   */
  public bool approxEquals(S1Interval y, double maxError) {
    if (isEmpty()) {
      return y.getLength() <= maxError;
    }
    if (y.isEmpty()) {
      return getLength() <= maxError;
    }
    return (Math.Abs(Math.IEEERemainder(y.lo() - lo(), 2 * S2.M_PI))
        + Math.Abs(Math.IEEERemainder(y.hi() - hi(), 2 * S2.M_PI))) <= maxError;
  }

  public bool approxEquals(S1Interval y) {
    return approxEquals(y, 1e-9);
  }

 
        public override string ToString()
        {
    return "[" + this.lo() + ", " + this.hi() + "]";
        }

   

  /**
   * Compute the distance from "a" to "b" in the range [0, 2*Pi). This is
   * equivalent to (drem(b - a - S2.M_PI, 2 * S2.M_PI) + S2.M_PI), except that
   * it is more numerically stable (it does not lose precision for very small
   * positive distances).
   */
  public static double positiveDistance(double a, double b) {
    double d = b - a;
    if (d >= 0) {
      return d;
    }
    // We want to ensure that if b == Pi and a == (-Pi + eps),
    // the return result is approximately 2*Pi and not zero.
    return (b + S2.M_PI) - (a - S2.M_PI);
  }

    }
}
