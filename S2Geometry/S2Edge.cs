using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Google.Common.Geometry
{
    /**
 * An abstract directed edge from one S2Point to another S2Point.
 *
 * @author kirilll@google.com (Kirill Levin)
 */

    public struct S2Edge : IEquatable<S2Edge>
    {
        private readonly S2Point _end;
        private readonly S2Point _start;

        public S2Edge(S2Point start, S2Point end)
        {
            _start = start;
            _end = end;
        }

        public S2Point Start
        {
            get { return _start; }
        }

        public S2Point End
        {
            get { return _end; }
        }

        public bool Equals(S2Edge other)
        {
            return _end.Equals(other._end) && _start.Equals(other._start);
        }

        public override bool Equals(object obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            return obj is S2Edge && Equals((S2Edge)obj);
        }

        public override int GetHashCode()
        {
            unchecked
            {
                return (_end.GetHashCode()*397) ^ _start.GetHashCode();
            }
        }


        public static bool operator ==(S2Edge left, S2Edge right)
        {
            return Equals(left, right);
        }

        public static bool operator !=(S2Edge left, S2Edge right)
        {
            return !Equals(left, right);
        }

        public override string ToString()
        {
            return string.Format("Edge: ({0} -> {1})\n   or [{2} -> {3}]",
                                 _start.ToDegreesString(), _end.ToDegreesString(), _start, _end);
        }
    }
}