using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Google.Common.Geometry
{

    /**
     * An S2Region represents a two-dimensional region over the unit sphere. It is
     * an abstract interface with various concrete subtypes.
     *
     *  The main purpose of this interface is to allow complex regions to be
     * approximated as simpler regions. So rather than having a wide variety of
     * virtual methods that are implemented by all subtypes, the interface is
     * restricted to methods that are useful for computing approximations.
     *
     *
     */
    public interface IS2Region
    {
        /** Return a bounding spherical cap. */
        S2Cap getCapBound();


        /** Return a bounding latitude-longitude rectangle. */
        S2LatLngRect getRectBound();

        /**
         * If this method returns true, the region completely contains the given cell.
         * Otherwise, either the region does not contain the cell or the containment
         * relationship could not be determined.
         */
        bool contains(S2Cell cell);

        /**
         * If this method returns false, the region does not intersect the given cell.
         * Otherwise, either region intersects the cell, or the intersection
         * relationship could not be determined.
         */
        bool mayIntersect(S2Cell cell);
    }
}
