using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Google.Common.Geometry.DataStructures;

namespace Google.Common.Geometry
{
    /**
 * An S2RegionCoverer is a class that allows arbitrary regions to be
 * approximated as unions of cells (S2CellUnion). This is useful for
 * implementing various sorts of search and precomputation operations.
 *
 * Typical usage: {@code S2RegionCoverer coverer; coverer.setMaxCells(5); S2Cap
 * cap = S2Cap.fromAxisAngle(...); S2CellUnion covering;
 * coverer.getCovering(cap, covering); * }
 *
 * This yields a cell union of at most 5 cells that is guaranteed to cover the
 * given cap (a disc-shaped region on the sphere).
 *
 *  The approximation algorithm is not optimal but does a pretty good job in
 * practice. The output does not always use the maximum number of cells allowed,
 * both because this would not always yield a better approximation, and because
 * max_cells() is a limit on how much work is done exploring the possible
 * covering as well as a limit on the final output size.
 *
 *  One can also generate interior coverings, which are sets of cells which are
 * entirely contained within a region. Interior coverings can be empty, even for
 * non-empty regions, if there are no cells that satisfy the provided
 * constraints and are contained by the region. Note that for performance
 * reasons, it is wise to specify a max_level when computing interior coverings
 * - otherwise for regions with small or zero area, the algorithm may spend a
 * lot of time subdividing cells all the way to leaf level to try to find
 * contained cells.
 *
 *  This class is thread-unsafe. Simultaneous calls to any of the getCovering
 * methods will conflict and produce unpredictable results.
 *
 */

    public sealed class S2RegionCoverer
    {
        /**
   * By default, the covering uses at most 8 cells at any level. This gives a
   * reasonable tradeoff between the number of cells used and the accuracy of
   * the approximation (see table below).
   */
        public const int DefaultMaxCells = 8;

        private static readonly S2Cell[] FaceCells = new S2Cell[6];
        private readonly PriorityQueue<QueueEntry> _candidateQueue;
        private readonly List<S2CellId> _result;
        private int _candidatesCreatedCounter;
        private bool _interiorCovering;

        private int _levelMod;
        private int _maxCells;
        private int _maxLevel;
        private int _minLevel;

        // True if we're computing an interior covering.

        // Counter of number of candidates created, for performance evaluation.

        /**
   * We save a temporary copy of the pointer passed to GetCovering() in order to
   * avoid passing this parameter around internally. It is only used (and only
   * valid) for the duration of a single GetCovering() call.
   */
        private IS2Region _region;

        static S2RegionCoverer()
        {
            for (var face = 0; face < 6; ++face)
            {
                FaceCells[face] = S2Cell.FromFacePosLevel(face, (byte)0, 0);
            }
        }

        /**
   * A temporary variable used by GetCovering() that holds the cell ids that
   * have been added to the covering so far.
   */

        /**
   * Default constructor, sets all fields to default values.
   */

        public S2RegionCoverer()
        {
            _minLevel = 0;
            _maxLevel = S2CellId.MaxLevel;
            _levelMod = 1;
            _maxCells = DefaultMaxCells;
            _region = null;
            _result = new List<S2CellId>();
            // TODO(kirilll?): 10 is a completely random number, work out a better
            // estimate
            _candidateQueue = new PriorityQueue<QueueEntry>();
        }

        // Set the minimum and maximum cell level to be used. The default is to use
        // all cell levels. Requires: max_level() >= min_level().
        //
        // To find the cell level corresponding to a given physical distance, use
        // the S2Cell metrics defined in s2.h. For example, to find the cell
        // level that corresponds to an average edge length of 10km, use:
        //
        // int level = S2::kAvgEdge.GetClosestLevel(
        // geostore::S2Earth::KmToRadians(length_km));
        //
        // Note: min_level() takes priority over max_cells(), i.e. cells below the
        // given level will never be used even if this causes a large number of
        // cells to be returned.

        /**
   * Sets the minimum level to be used.
   */

        public int MinLevel
        {
            set
            {
                // assert (minLevel >= 0 && minLevel <= S2CellId.MAX_LEVEL);
                _minLevel = Math.Max(0, Math.Min(S2CellId.MaxLevel, value));
            }
            get { return _minLevel; }
        }

        /**
   * Sets the maximum level to be used.
   */

        public int MaxLevel
        {
            set
            {
                // assert (maxLevel >= 0 && maxLevel <= S2CellId.MAX_LEVEL);
                _maxLevel = Math.Max(0, Math.Min(S2CellId.MaxLevel, value));
            }
            get { return _maxLevel; }
        }

        /**
  * Sets the maximum desired number of cells in the approximation (defaults to
  * kDefaultMaxCells). Note the following:
  *
  * <ul>
  * <li>For any setting of max_cells(), up to 6 cells may be returned if that
  * is the minimum number of cells required (e.g. if the region intersects all
  * six face cells). Up to 3 cells may be returned even for very tiny convex
  * regions if they happen to be located at the intersection of three cube
  * faces.
  *
  * <li>For any setting of max_cells(), an arbitrary number of cells may be
  * returned if min_level() is too high for the region being approximated.
  *
  * <li>If max_cells() is less than 4, the area of the covering may be
  * arbitrarily large compared to the area of the original region even if the
  * region is convex (e.g. an S2Cap or S2LatLngRect).
  * </ul>
  *
  * Accuracy is measured by dividing the area of the covering by the area of
  * the original region. The following table shows the median and worst case
  * values for this area ratio on a test case consisting of 100,000 spherical
  * caps of random size (generated using s2regioncoverer_unittest):
  *
  * <pre>
  * max_cells: 3 4 5 6 8 12 20 100 1000
  * median ratio: 5.33 3.32 2.73 2.34 1.98 1.66 1.42 1.11 1.01
  * worst case: 215518 14.41 9.72 5.26 3.91 2.75 1.92 1.20 1.02
  * </pre>
  */

        public int MaxCells
        {
            get { return _maxCells; }
            set { _maxCells = value; }
        }

        /**
   * If specified, then only cells where (level - min_level) is a multiple of
   * "level_mod" will be used (default 1). This effectively allows the branching
   * factor of the S2CellId hierarchy to be increased. Currently the only
   * parameter values allowed are 1, 2, or 3, corresponding to branching factors
   * of 4, 16, and 64 respectively.
   */

        public int LevelMod
        {
            set
            {
                Debug.Assert(value >= 1 && value <= 3);
                _levelMod = Math.Max(1, Math.Min(3, value));
            }
            get { return _levelMod; }
        }

        private int MaxChildrenShift
        {
            get { return 2*_levelMod; }
        }


        /**
   * Computes a list of cell ids that covers the given region and satisfies the
   * various restrictions specified above.
   *
   * @param region The region to cover
   * @param covering The list filled in by this method
   */

        public void GetCovering(IS2Region region, ICollection<S2CellId> covering)
        {
            // Rather than just returning the raw list of cell ids generated by
            // GetCoveringInternal(), we construct a cell union and then denormalize it.
            // This has the effect of replacing four child cells with their parent
            // whenever this does not violate the covering parameters specified
            // (min_level, level_mod, etc). This strategy significantly reduces the
            // number of cells returned in many cases, and it is cheap compared to
            // computing the covering in the first place.

            var tmp = GetCovering(region);
            tmp.Denormalize(MinLevel, LevelMod, covering);
        }

        /**
   * Computes a list of cell ids that is contained within the given region and
   * satisfies the various restrictions specified above.
   *
   * @param region The region to fill
   * @param interior The list filled in by this method
   */

        public void GetInteriorCovering(IS2Region region, List<S2CellId> interior)
        {
            var tmp = GetInteriorCovering(region);
            tmp.Denormalize(MinLevel, LevelMod, interior);
        }

        /**
   * Return a normalized cell union that covers the given region and satisfies
   * the restrictions *EXCEPT* for min_level() and level_mod(). These criteria
   * cannot be satisfied using a cell union because cell unions are
   * automatically normalized by replacing four child cells with their parent
   * whenever possible. (Note that the list of cell ids passed to the cell union
   * constructor does in fact satisfy all the given restrictions.)
   */

        public S2CellUnion GetCovering(IS2Region region)
        {
            var covering = new S2CellUnion();
            GetCovering(region, covering);
            return covering;
        }

        public void GetCovering(IS2Region region, S2CellUnion covering)
        {
            _interiorCovering = false;
            GetCoveringInternal(region);
            covering.InitSwap(_result);
        }

        /**
   * Return a normalized cell union that is contained within the given region
   * and satisfies the restrictions *EXCEPT* for min_level() and level_mod().
   */

        public S2CellUnion GetInteriorCovering(IS2Region region)
        {
            var covering = new S2CellUnion();
            GetInteriorCovering(region, covering);
            return covering;
        }

        public void GetInteriorCovering(IS2Region region, S2CellUnion covering)
        {
            _interiorCovering = true;
            GetCoveringInternal(region);
            covering.InitSwap(_result);
        }

        /**
   * Given a connected region and a starting point, return a set of cells at the
   * given level that cover the region.
   */

        public static void GetSimpleCovering(
            IS2Region region, S2Point start, int level, List<S2CellId> output)
        {
            FloodFill(region, S2CellId.FromPoint(start).ParentForLevel(level), output);
        }

        /**
   * If the cell intersects the given region, return a new candidate with no
   * children, otherwise return null. Also marks the candidate as "terminal" if
   * it should not be expanded further.
   */

        private Candidate NewCandidate(S2Cell cell)
        {
            if (!_region.MayIntersect(cell))
            {
                return null;
            }

            var isTerminal = false;
            if (cell.Level >= _minLevel)
            {
                if (_interiorCovering)
                {
                    if (_region.Contains(cell))
                    {
                        isTerminal = true;
                    }
                    else if (cell.Level + _levelMod > _maxLevel)
                    {
                        return null;
                    }
                }
                else
                {
                    if (cell.Level + _levelMod > _maxLevel || _region.Contains(cell))
                    {
                        isTerminal = true;
                    }
                }
            }
            var candidate = new Candidate();
            candidate.Cell = cell;
            candidate.IsTerminal = isTerminal;
            if (!isTerminal)
            {
                candidate.Children = new Candidate[1 << MaxChildrenShift];
            }
            _candidatesCreatedCounter++;
            return candidate;
        }

        /** Return the log base 2 of the maximum number of children of a candidate. */

        /**
   * Process a candidate by either adding it to the result list or expanding its
   * children and inserting it into the priority queue. Passing an argument of
   * NULL does nothing.
   */

        private void AddCandidate(Candidate candidate)
        {
            if (candidate == null)
            {
                return;
            }

            if (candidate.IsTerminal)
            {
                _result.Add(candidate.Cell.Id);
                return;
            }
            // assert (candidate.numChildren == 0);

            // Expand one level at a time until we hit min_level_ to ensure that
            // we don't skip over it.
            var numLevels = (candidate.Cell.Level < _minLevel) ? 1 : _levelMod;
            var numTerminals = ExpandChildren(candidate, candidate.Cell, numLevels);

            if (candidate.NumChildren == 0)
            {
                // Do nothing
            }
            else if (!_interiorCovering && numTerminals == 1 << MaxChildrenShift
                     && candidate.Cell.Level >= _minLevel)
            {
                // Optimization: add the parent cell rather than all of its children.
                // We can't do this for interior coverings, since the children just
                // intersect the region, but may not be contained by it - we need to
                // subdivide them further.
                candidate.IsTerminal = true;
                AddCandidate(candidate);
            }
            else
            {
                // We negate the priority so that smaller absolute priorities are returned
                // first. The heuristic is designed to refine the largest cells first,
                // since those are where we have the largest potential gain. Among cells
                // at the same level, we prefer the cells with the smallest number of
                // intersecting children. Finally, we prefer cells that have the smallest
                // number of children that cannot be refined any further.
                var priority = -((((candidate.Cell.Level << MaxChildrenShift) + candidate.NumChildren)
                                  << MaxChildrenShift) + numTerminals);
                var entry = new QueueEntry(priority, candidate);
                _candidateQueue.Enqueue(entry);
                // logger.info("Push: " + candidate.cell.id() + " (" + priority + ") ");
            }
        }

        /**
   * Populate the children of "candidate" by expanding the given number of
   * levels from the given cell. Returns the number of children that were marked
   * "terminal".
   */

        private int ExpandChildren(Candidate candidate, S2Cell cell, int numLevels)
        {
            numLevels--;
            var childCells = new S2Cell[4];
            for (var i = 0; i < 4; ++i)
            {
                childCells[i] = new S2Cell();
            }
            cell.Subdivide(childCells);
            var numTerminals = 0;
            for (var i = 0; i < 4; ++i)
            {
                if (numLevels > 0)
                {
                    if (_region.MayIntersect(childCells[i]))
                    {
                        numTerminals += ExpandChildren(candidate, childCells[i], numLevels);
                    }
                    continue;
                }
                var child = NewCandidate(childCells[i]);
                if (child != null)
                {
                    candidate.Children[candidate.NumChildren++] = child;
                    if (child.IsTerminal)
                    {
                        ++numTerminals;
                    }
                }
            }
            return numTerminals;
        }

        /** Computes a set of initial candidates that cover the given region. */

        private void GetInitialCandidates()
        {
            // Optimization: if at least 4 cells are desired (the normal case),
            // start with a 4-cell covering of the region's bounding cap. This
            // lets us skip quite a few levels of refinement when the region to
            // be covered is relatively small.
            if (_maxCells >= 4)
            {
                // Find the maximum level such that the bounding cap contains at most one
                // cell vertex at that level.
                var cap = _region.CapBound;
                var level = Math.Min(S2Projections.MinWidth.GetMaxLevel(2*cap.Angle.Radians),
                                     Math.Min(MaxLevel, S2CellId.MaxLevel - 1));
                if (LevelMod > 1 && level > MinLevel)
                {
                    level -= (level - MinLevel)%LevelMod;
                }
                // We don't bother trying to optimize the level == 0 case, since more than
                // four face cells may be required.
                if (level > 0)
                {
                    // Find the leaf cell containing the cap axis, and determine which
                    // subcell of the parent cell contains it.
                    var @base = new List<S2CellId>(4);
                    var id = S2CellId.FromPoint(cap.Axis);
                    id.GetVertexNeighbors(level, @base);
                    for (var i = 0; i < @base.Count; ++i)
                    {
                        AddCandidate(NewCandidate(new S2Cell(@base[i])));
                    }
                    return;
                }
            }
            // Default: start with all six cube faces.
            for (var face = 0; face < 6; ++face)
            {
                AddCandidate(NewCandidate(FaceCells[face]));
            }
        }

        /** Generates a covering and stores it in result. */

        private void GetCoveringInternal(IS2Region region)
        {
            // Strategy: Start with the 6 faces of the cube. Discard any
            // that do not intersect the shape. Then repeatedly choose the
            // largest cell that intersects the shape and subdivide it.
            //
            // result contains the cells that will be part of the output, while the
            // priority queue contains cells that we may still subdivide further. Cells
            // that are entirely contained within the region are immediately added to
            // the output, while cells that do not intersect the region are immediately
            // discarded.
            // Therefore pq_ only contains cells that partially intersect the region.
            // Candidates are prioritized first according to cell size (larger cells
            // first), then by the number of intersecting children they have (fewest
            // children first), and then by the number of fully contained children
            // (fewest children first).

            Preconditions.CheckState(_candidateQueue.Count == 0 && _result.Count == 0);

            _region = region;
            _candidatesCreatedCounter = 0;

            GetInitialCandidates();
            while (_candidateQueue.Count != 0 && (!_interiorCovering || _result.Count < _maxCells))
            {
                var qe = _candidateQueue.Dequeue();
                var candidate = qe.Candidate;
                // logger.info("Pop: " + candidate.cell.id());
                if (candidate.Cell.Level < _minLevel || candidate.NumChildren == 1
                    || _result.Count + (_interiorCovering ? 0 : _candidateQueue.Count) + candidate.NumChildren
                    <= _maxCells)
                {
                    // Expand this candidate into its children.
                    for (var i = 0; i < candidate.NumChildren; ++i)
                    {
                        AddCandidate(candidate.Children[i]);
                    }
                }
                else if (_interiorCovering)
                {
                    // Do nothing
                }
                else
                {
                    candidate.IsTerminal = true;
                    AddCandidate(candidate);
                }
            }

            _candidateQueue.Clear();
            _region = null;
        }

        /**
   * Given a region and a starting cell, return the set of all the
   * edge-connected cells at the same level that intersect "region". The output
   * cells are returned in arbitrary order.
   */

        private static void FloodFill(IS2Region region, S2CellId start, List<S2CellId> output)
        {
            var all = new HashSet<S2CellId>();
            var frontier = new List<S2CellId>();
            output.Clear();
            all.Add(start);
            frontier.Add(start);
            while (frontier.Any())
            {
                var id = frontier[frontier.Count - 1];
                frontier.RemoveAt(frontier.Count - 1);
                if (!region.MayIntersect(new S2Cell(id)))
                {
                    continue;
                }
                output.Add(id);

                var neighbors = id.GetEdgeNeighbors();
                for (var edge = 0; edge < 4; ++edge)
                {
                    var nbr = neighbors[edge];
                    var hasNbr = all.Contains(nbr);
                    if (!all.Contains(nbr))
                    {
                        frontier.Add(nbr);
                        all.Add(nbr);
                    }
                }
            }
        }

        private sealed class Candidate
        {
            public S2Cell Cell;
            public Candidate[] Children; // Actual size may be 0, 4, 16, or 64
            public bool IsTerminal; // Cell should not be expanded further.
            public int NumChildren; // Number of children that intersect the region.
            // elements.
        }

        private struct QueueEntry : IComparable<QueueEntry>
        {
            public readonly Candidate Candidate;
            private readonly int _id;

            public QueueEntry(int id, Candidate candidate)
            {
                _id = id;
                Candidate = candidate;
            }

            public int CompareTo(QueueEntry other)
            {
                return _id < other._id ? 1 : (_id > other._id ? -1 : 0);
            }
        }
    }
}