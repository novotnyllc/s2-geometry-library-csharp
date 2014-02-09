using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Google.Common.Geometry
{
    internal class Preconditions
    {
        public static void checkArgument(bool expression, string message = null)
        {
            if (!expression)
                throw new ArgumentException(message ?? string.Empty);
        }

        public static void checkState(bool expression, string message = null)
        {
            if (!expression)
                throw new InvalidOperationException(message ?? "bad state");
        }
    }
}