using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Google.Common.Geometry
{
    public sealed class R2Vector : IEquatable<R2Vector>
    {
        private readonly double _x;
        private readonly double _y;

        public R2Vector() : this(0, 0)
        {
        }

        public R2Vector(double x, double y)
        {
            _x = x;
            _y = y;
        }

        public R2Vector(double[] coord)
        {
            if (coord.Length != 2)
            {
                throw new ArgumentException("Points must have exactly 2 coordinates", "coord");
            }
            _x = coord[0];
            _y = coord[1];
        }

        public bool Equals(R2Vector other)
        {
            if (ReferenceEquals(null, other)) return false;
            if (ReferenceEquals(this, other)) return true;
            return _y.Equals(other._y) && _x.Equals(other._x);
        }

        public override bool Equals(object obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            if (ReferenceEquals(this, obj)) return true;
            return obj is R2Vector && Equals((R2Vector)obj);
        }


        /**
     * Calcualates hashcode based on stored coordinates. Since we want +0.0 and
     * -0.0 to be treated the same, we ignore the sign of the coordinates.
     */

        public override int GetHashCode()
        {
            unchecked
            {
                return (Math.Abs(_y).GetHashCode()*397) ^ Math.Abs(_x).GetHashCode();
            }
        }

        public static bool operator ==(R2Vector left, R2Vector right)
        {
            return Equals(left, right);
        }

        public static bool operator !=(R2Vector left, R2Vector right)
        {
            return !Equals(left, right);
        }

        public double x()
        {
            return _x;
        }

        public double y()
        {
            return _y;
        }

        public double get(int index)
        {
            if (index > 1)
            {
                throw new ArgumentOutOfRangeException("index");
            }
            return index == 0 ? _x : _y;
        }

        public static R2Vector add(R2Vector p1, R2Vector p2)
        {
            return new R2Vector(p1._x + p2._x, p1._y + p2._y);
        }

        public static R2Vector mul(R2Vector p, double m)
        {
            return new R2Vector(m*p._x, m*p._y);
        }

        public double norm2()
        {
            return (_x*_x) + (_y*_y);
        }

        public static double dotProd(R2Vector p1, R2Vector p2)
        {
            return (p1._x*p2._x) + (p1._y*p2._y);
        }

        public double dotProd(R2Vector that)
        {
            return dotProd(this, that);
        }

        public double crossProd(R2Vector that)
        {
            return _x*that._y - _y*that._x;
        }

        public bool lessThan(R2Vector vb)
        {
            if (_x < vb._x)
            {
                return true;
            }
            if (vb._x < _x)
            {
                return false;
            }
            if (_y < vb._y)
            {
                return true;
            }
            return false;
        }


        public override string ToString()
        {
            return "(" + _x + ", " + _y + ")";
        }
    }
}