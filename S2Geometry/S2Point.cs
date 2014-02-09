using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Google.Common.Geometry
{
    public sealed class S2Point : IEquatable<S2Point>, IComparable<S2Point>
    {
// coordinates of the points
        public readonly double x;
        public readonly double y;
        public readonly double z;

        public S2Point()
        {
        }

        public S2Point(double x, double y, double z)
        {
            this.x = x;
            this.y = y;
            this.z = z;
        }

        public int CompareTo(S2Point other)
        {
            return (lessThan(other) ? -1 : (Equals(other) ? 0 : 1));
        }

        public bool Equals(S2Point other)
        {
            if (ReferenceEquals(null, other)) return false;
            if (ReferenceEquals(this, other)) return true;
            return x.Equals(other.x) && y.Equals(other.y) && z.Equals(other.z);
        }

        public override bool Equals(object obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            if (ReferenceEquals(this, obj)) return true;
            if (obj.GetType() != GetType()) return false;
            return Equals((S2Point)obj);
        }

        /**
   * Calcualates hashcode based on stored coordinates. Since we want +0.0 and
   * -0.0 to be treated the same, we ignore the sign of the coordinates.
   */

        public override int GetHashCode()
        {
            unchecked
            {
                var hashCode = Math.Abs(x).GetHashCode();
                hashCode = (hashCode*397) ^ Math.Abs(y).GetHashCode();
                hashCode = (hashCode*397) ^ Math.Abs(z).GetHashCode();
                return hashCode;
            }
        }

        public static bool operator ==(S2Point left, S2Point right)
        {
            return Equals(left, right);
        }

        public static bool operator !=(S2Point left, S2Point right)
        {
            return !Equals(left, right);
        }

        public static S2Point minus(S2Point p1, S2Point p2)
        {
            return sub(p1, p2);
        }

        public static S2Point neg(S2Point p)
        {
            return new S2Point(-p.x, -p.y, -p.z);
        }

        public double norm2()
        {
            return x*x + y*y + z*z;
        }

        public double norm()
        {
            return Math.Sqrt(norm2());
        }

        public static S2Point crossProd(S2Point p1, S2Point p2)
        {
            return new S2Point(
                p1.y*p2.z - p1.z*p2.y, p1.z*p2.x - p1.x*p2.z, p1.x*p2.y - p1.y*p2.x);
        }

        public static S2Point add(S2Point p1, S2Point p2)
        {
            return new S2Point(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
        }

        public static S2Point sub(S2Point p1, S2Point p2)
        {
            return new S2Point(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
        }

        public double dotProd(S2Point that)
        {
            return x*that.x + y*that.y + z*that.z;
        }

        public static S2Point mul(S2Point p, double m)
        {
            return new S2Point(m*p.x, m*p.y, m*p.z);
        }

        public static S2Point div(S2Point p, double m)
        {
            return new S2Point(p.x/m, p.y/m, p.z/m);
        }

        /** return a vector orthogonal to this one */

        public S2Point ortho()
        {
            var k = largestAbsComponent();
            S2Point temp;
            if (k == 1)
            {
                temp = new S2Point(1, 0, 0);
            }
            else if (k == 2)
            {
                temp = new S2Point(0, 1, 0);
            }
            else
            {
                temp = new S2Point(0, 0, 1);
            }
            return normalize(crossProd(this, temp));
        }

        /** Return the index of the largest component fabs */

        public int largestAbsComponent()
        {
            var temp = fabs(this);
            if (temp.x > temp.y)
            {
                if (temp.x > temp.z)
                {
                    return 0;
                }
                else
                {
                    return 2;
                }
            }
            else
            {
                if (temp.y > temp.z)
                {
                    return 1;
                }
                else
                {
                    return 2;
                }
            }
        }

        public static S2Point fabs(S2Point p)
        {
            return new S2Point(Math.Abs(p.x), Math.Abs(p.y), Math.Abs(p.z));
        }

        public static S2Point normalize(S2Point p)
        {
            var norm = p.norm();
            if (norm != 0)
            {
                norm = 1.0/norm;
            }
            return mul(p, norm);
        }

        public double get(int axis)
        {
            return (axis == 0) ? x : (axis == 1) ? y : z;
        }

        /** Return the angle between two vectors in radians */

        public double angle(S2Point va)
        {
            return Math.Atan2(crossProd(this, va).norm(), dotProd(va));
        }

        /**
   * Compare two vectors, return true if all their components are within a
   * difference of margin.
   */

        public bool aequal(S2Point that, double margin)
        {
            return (Math.Abs(x - that.x) < margin) && (Math.Abs(y - that.y) < margin)
                   && (Math.Abs(z - that.z) < margin);
        }


        public bool lessThan(S2Point vb)
        {
            if (x < vb.x)
            {
                return true;
            }
            if (vb.x < x)
            {
                return false;
            }
            if (y < vb.y)
            {
                return true;
            }
            if (vb.y < y)
            {
                return false;
            }
            if (z < vb.z)
            {
                return true;
            }
            return false;
        }


        public override string ToString()
        {
            return "(" + x + ", " + y + ", " + z + ")";
        }

        public string toDegreesString()
        {
            var s2LatLng = new S2LatLng(this);
            return "(" + s2LatLng.latDegrees() + ", "
                   + s2LatLng.lngDegrees() + ")";
        }
    }
}