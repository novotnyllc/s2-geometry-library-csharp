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

    public sealed class S2Edge : IEquatable<S2Edge>
    {
        private readonly S2Point end;
        private readonly S2Point start;

        public S2Edge(S2Point start, S2Point end)
        {
            this.start = start;
            this.end = end;
        }

        public bool Equals(S2Edge other)
        {
            if (ReferenceEquals(null, other)) return false;
            if (ReferenceEquals(this, other)) return true;
            return start.Equals(other.start) && end.Equals(other.end);
        }

        public override bool Equals(object obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            if (ReferenceEquals(this, obj)) return true;
            return obj is S2Edge && Equals((S2Edge)obj);
        }

        public override int GetHashCode()
        {
            unchecked
            {
                return (start.GetHashCode()*397) ^ end.GetHashCode();
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

        public S2Point getStart()
        {
            return start;
        }

        public S2Point getEnd()
        {
            return end;
        }

        public override String ToString()
        {
            return String.Format("Edge: ({0} -> {1})\n   or [{2} -> {3}]",
                                 start.toDegreesString(), end.toDegreesString(), start, end);
        }
    }
}