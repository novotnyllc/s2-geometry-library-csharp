using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Google.Common.Geometry
{
    public sealed class S1Angle : IEquatable<S1Angle>, IComparable<S1Angle>
    {
        public bool Equals(S1Angle other)
        {
            if (ReferenceEquals(null, other)) return false;
            if (ReferenceEquals(this, other)) return true;
            return _radians.Equals(other._radians);
        }

        public int CompareTo(S1Angle other)
        {
            return this._radians < other._radians ? -1 : this._radians > other._radians ? 1 : 0;
        }

        public override bool Equals(object obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            if (ReferenceEquals(this, obj)) return true;
            if (obj.GetType() != this.GetType()) return false;
            return Equals((S1Angle)obj);
        }

        public override int GetHashCode()
        {
            return _radians.GetHashCode();
        }

        public static bool operator ==(S1Angle left, S1Angle right)
        {
            return Equals(left, right);
        }

        public static bool operator !=(S1Angle left, S1Angle right)
        {
            return !Equals(left, right);
        }

        private readonly double _radians;

  public double radians() {
    return _radians;
  }

  public double degrees() {
    return _radians * (180 / Math.PI);
  }

  public long e5() {
    return (long)Math.Round(degrees() * 1e5);
  }

  public long e6() {
    return (long)Math.Round(degrees() * 1e6);
  }

  public long e7() {
    return (long)Math.Round(degrees() * 1e7);
  }

  /**
   * The default constructor yields a zero angle.
   */
  public S1Angle() {
    this._radians = 0;
  }

  private S1Angle(double radians) {
    this._radians = radians;
  }

  /**
   * Return the angle between two points, which is also equal to the distance
   * between these points on the unit sphere. The points do not need to be
   * normalized.
   */
  public S1Angle(S2Point x, S2Point y) {
    this._radians = x.angle(y);
  }

 
  public bool lessThan(S1Angle that) {
    return this.radians() < that.radians();
  }

  public bool greaterThan(S1Angle that) {
    return this.radians() > that.radians();
  }

  public bool lessOrEquals(S1Angle that) {
    return this.radians() <= that.radians();
  }

  public bool greaterOrEquals(S1Angle that) {
    return this.radians() >= that.radians();
  }

  public static S1Angle max(S1Angle left, S1Angle right) {
    return right.greaterThan(left) ? right : left;
  }

  public static S1Angle min(S1Angle left, S1Angle right) {
    return right.greaterThan(left) ? left : right;
  }

  public static S1Angle radians(double radians) {
    return new S1Angle(radians);
  }

  public static S1Angle degrees(double degrees) {
    return new S1Angle(degrees * (Math.PI / 180));
  }

  public static S1Angle e5(long e5) {
    return degrees(e5 * 1e-5);
  }

  public static S1Angle e6(long e6) {
    // Multiplying by 1e-6 isn't quite as accurate as dividing by 1e6,
    // but it's about 10 times faster and more than accurate enough.
    return degrees(e6 * 1e-6);
  }

  public static S1Angle e7(long e7) {
    return degrees(e7 * 1e-7);
  }

  /**
   * Writes the angle in degrees with a "d" suffix, e.g. "17.3745d". By default
   * 6 digits are printed; this can be changed using setprecision(). Up to 17
   * digits are required to distinguish one angle from another.
   */
  public override string ToString() {
    return degrees() + "d";
  }


    }
}
