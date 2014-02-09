using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Google.Common.Geometry
{
    public abstract class S2EdgeIndex
    {
        /**
   * Thicken the edge in all directions by roughly 1% of the edge length when
   * thickenEdge is true.
   */
        private const double THICKENING = 0.01;

        /**
   * Threshold for small angles, that help lenientCrossing to determine whether
   * two edges are likely to intersect.
   */
        private const double MAX_DET_ERROR = 1e-14;

        /**
   * The cell containing each edge, as given in the parallel array
   * <code>edges</code>.
   */
        private ulong[] cells;

        /**
   * The edge contained by each cell, as given in the parallel array
   * <code>cells</code>.
   */
        private int[] edges;

        /**
   * No cell strictly below this level appears in mapping. Initially leaf level,
   * that's the minimum level at which we will ever look for test edges.
   */

        /**
   * Has the index been computed already?
   */
        private bool indexComputed;
        private int minimumS2LevelUsed;

        /**
   * Number of queries so far
   */
        private int queryCount;

        /**
   * Empties the index in case it already contained something.
   */

        public void reset()
        {
            minimumS2LevelUsed = S2CellId.MAX_LEVEL;
            indexComputed = false;
            queryCount = 0;
            cells = null;
            edges = null;
        }

        /**
   * Compares [cell1, edge1] to [cell2, edge2], by cell first and edge second.
   *
   * @return -1 if [cell1, edge1] is less than [cell2, edge2], 1 if [cell1,
   *         edge1] is greater than [cell2, edge2], 0 otherwise.
   */

        private static int compare(ulong cell1, int edge1, ulong cell2, int edge2)
        {
            if (cell1 < cell2)
            {
                return -1;
            }
            else if (cell1 > cell2)
            {
                return 1;
            }
            else if (edge1 < edge2)
            {
                return -1;
            }
            else if (edge1 > edge2)
            {
                return 1;
            }
            else
            {
                return 0;
            }
        }

        /** Computes the index (if it has not been previously done). */

        public void computeIndex()
        {
            if (indexComputed)
            {
                return;
            }
            var cellList = new List<ulong>();
            var edgeList = new List<int>();
            for (var i = 0; i < getNumEdges(); ++i)
            {
                var from = edgeFrom(i);
                var to = edgeTo(i);
                var cover = new List<S2CellId>();
                var level = getCovering(from, to, true, cover);
                minimumS2LevelUsed = Math.Min(minimumS2LevelUsed, level);
                foreach (var cellId in cover)
                {
                    cellList.Add(cellId.id());
                    edgeList.Add(i);
                }
            }
            cells = new ulong[cellList.Count];
            edges = new int[edgeList.Count];
            for (var i = 0; i < cells.Length; i++)
            {
                cells[i] = cellList[i];
                edges[i] = edgeList[i];
            }
            sortIndex();
            indexComputed = true;
        }

        /** Sorts the parallel <code>cells</code> and <code>edges</code> arrays. */

        private void sortIndex()
        {
            // create an array of indices and sort based on the values in the parallel
            // arrays at each index
            var indices = new int[cells.Length];
            for (var i = 0; i < indices.Length; i++)
            {
                indices[i] = i;
            }

            Array.Sort(indices, (index1, index2) => compare(cells[index1], edges[index1], cells[index2], edges[index2]));


            // copy the cells and edges in the order given by the sorted list of indices
            var newCells = new ulong[cells.Length];
            var newEdges = new int[edges.Length];
            for (var i = 0; i < indices.Length; i++)
            {
                newCells[i] = cells[indices[i]];
                newEdges[i] = edges[indices[i]];
            }
            // replace the cells and edges with the sorted arrays
            cells = newCells;
            edges = newEdges;
        }

        public bool isIndexComputed()
        {
            return indexComputed;
        }

        /**
   * Tell the index that we just received a new request for candidates. Useful
   * to compute when to switch to quad tree.
   */

        protected void incrementQueryCount()
        {
            ++queryCount;
        }

        /**
   * If the index hasn't been computed yet, looks at how much work has gone into
   * iterating using the brute force method, and how much more work is planned
   * as defined by 'cost'. If it were to have been cheaper to use a quad tree
   * from the beginning, then compute it now. This guarantees that we will never
   * use more than twice the time we would have used had we known in advance
   * exactly how many edges we would have wanted to test. It is the theoretical
   * best.
   *
   *  The value 'n' is the number of iterators we expect to request from this
   * edge index.
   *
   *  If we have m data edges and n query edges, then the brute force cost is m
   * * n * testCost where testCost is taken to be the cost of
   * EdgeCrosser.robustCrossing, measured to be about 30ns at the time of this
   * writing.
   *
   *  If we compute the index, the cost becomes: m * costInsert + n *
   * costFind(m)
   *
   *  - costInsert can be expected to be reasonably stable, and was measured at
   * 1200ns with the BM_QuadEdgeInsertionCost benchmark.
   *
   *  - costFind depends on the length of the edge . For m=1000 edges, we got
   * timings ranging from 1ms (edge the length of the polygon) to 40ms. The
   * latter is for very long query edges, and needs to be optimized. We will
   * assume for the rest of the discussion that costFind is roughly 3ms.
   *
   *  When doing one additional query, the differential cost is m * testCost -
   * costFind(m) With the numbers above, it is better to use the quad tree (if
   * we have it) if m >= 100.
   *
   *  If m = 100, 30 queries will give m*n*testCost = m * costInsert = 100ms,
   * while the marginal cost to find is 3ms. Thus, this is a reasonable thing to
   * do.
   */

        public void predictAdditionalCalls(int n)
        {
            if (indexComputed)
            {
                return;
            }
            if (getNumEdges() > 100 && (queryCount + n) > 30)
            {
                computeIndex();
            }
        }

        /**
   * Overwrite these functions to give access to the underlying data. The
   * function getNumEdges() returns the number of edges in the index, while
   * edgeFrom(index) and edgeTo(index) return the "from" and "to" endpoints of
   * the edge at the given index.
   */
        protected abstract int getNumEdges();

        protected abstract S2Point edgeFrom(int index);

        protected abstract S2Point edgeTo(int index);

        /**
   * Appends to "candidateCrossings" all edge references which may cross the
   * given edge. This is done by covering the edge and then finding all
   * references of edges whose coverings overlap this covering. Parent cells are
   * checked level by level. Child cells are checked all at once by taking
   * advantage of the natural ordering of S2CellIds.
   */

        protected void findCandidateCrossings(S2Point a, S2Point b, List<int> candidateCrossings)
        {
            Preconditions.CheckState(indexComputed);
            var cover = new List<S2CellId>();
            getCovering(a, b, false, cover);

            // Edge references are inserted into the map once for each covering cell, so
            // absorb duplicates here
            var uniqueSet = new HashSet<int>();
            getEdgesInParentCells(cover, uniqueSet);

            // TODO(user): An important optimization for long query
            // edges (Contains queries): keep a bounding cap and clip the query
            // edge to the cap before starting the descent.
            getEdgesInChildrenCells(a, b, cover, uniqueSet);

            candidateCrossings.Clear();
            candidateCrossings.AddRange(uniqueSet);
        }

        /**
   * Returns the smallest cell containing all four points, or
   * {@link S2CellId#sentinel()} if they are not all on the same face. The
   * points don't need to be normalized.
   */

        private static S2CellId containingCell(S2Point pa, S2Point pb, S2Point pc, S2Point pd)
        {
            var a = S2CellId.fromPoint(pa);
            var b = S2CellId.fromPoint(pb);
            var c = S2CellId.fromPoint(pc);
            var d = S2CellId.fromPoint(pd);

            if (a.face() != b.face() || a.face() != c.face() || a.face() != d.face())
            {
                return S2CellId.sentinel();
            }

            while (!a.Equals(b) || !a.Equals(c) || !a.Equals(d))
            {
                a = a.parent();
                b = b.parent();
                c = c.parent();
                d = d.parent();
            }
            return a;
        }

        /**
   * Returns the smallest cell containing both points, or Sentinel if they are
   * not all on the same face. The points don't need to be normalized.
   */

        private static S2CellId containingCell(S2Point pa, S2Point pb)
        {
            var a = S2CellId.fromPoint(pa);
            var b = S2CellId.fromPoint(pb);

            if (a.face() != b.face())
            {
                return S2CellId.sentinel();
            }

            while (!a.Equals(b))
            {
                a = a.parent();
                b = b.parent();
            }
            return a;
        }

        /**
   * Computes a cell covering of an edge. Clears edgeCovering and returns the
   * level of the s2 cells used in the covering (only one level is ever used for
   * each call).
   *
   *  If thickenEdge is true, the edge is thickened and extended by 1% of its
   * length.
   *
   *  It is guaranteed that no child of a covering cell will fully contain the
   * covered edge.
   */

        private int getCovering(
            S2Point a, S2Point b, bool thickenEdge, List<S2CellId> edgeCovering)
        {
            edgeCovering.Clear();

            // Selects the ideal s2 level at which to cover the edge, this will be the
            // level whose S2 cells have a width roughly commensurate to the length of
            // the edge. We multiply the edge length by 2*THICKENING to guarantee the
            // thickening is honored (it's not a big deal if we honor it when we don't
            // request it) when doing the covering-by-cap trick.
            var edgeLength = a.Angle(b);
            var idealLevel = S2Projections.MIN_WIDTH.GetMaxLevel(edgeLength*(1 + 2*THICKENING));

            S2CellId containingCellId;
            if (!thickenEdge)
            {
                containingCellId = containingCell(a, b);
            }
            else
            {
                if (idealLevel == S2CellId.MAX_LEVEL)
                {
                    // If the edge is tiny, instabilities are more likely, so we
                    // want to limit the number of operations.
                    // We pretend we are in a cell much larger so as to trigger the
                    // 'needs covering' case, so we won't try to thicken the edge.
                    containingCellId = (new S2CellId(0xFFF0)).parent(3);
                }
                else
                {
                    var pq = (b - a)*THICKENING;
                    var ortho = (S2Point.Normalize(S2Point.CrossProd(pq, a)))*edgeLength*THICKENING;
                    var p = a - pq;
                    var q = b + pq;
                    // If p and q were antipodal, the edge wouldn't be lengthened,
                    // and it could even flip! This is not a problem because
                    // idealLevel != 0 here. The farther p and q can be is roughly
                    // a quarter Earth away from each other, so we remain
                    // Theta(THICKENING).
                    containingCellId = containingCell(p - ortho, p + ortho, q - ortho, q + ortho);
                }
            }

            // Best case: edge is fully contained in a cell that's not too big.
            if (!containingCellId.Equals(S2CellId.sentinel())
                && containingCellId.level() >= idealLevel - 2)
            {
                edgeCovering.Add(containingCellId);
                return containingCellId.level();
            }

            if (idealLevel == 0)
            {
                // Edge is very long, maybe even longer than a face width, so the
                // trick below doesn't work. For now, we will add the whole S2 sphere.
                // TODO(user): Do something a tad smarter (and beware of the
                // antipodal case).
                for (var cellid = S2CellId.begin(0); !cellid.Equals(S2CellId.end(0));
                     cellid = cellid.next())
                {
                    edgeCovering.Add(cellid);
                }
                return 0;
            }
            // TODO(user): Check trick below works even when vertex is at
            // interface
            // between three faces.

            // Use trick as in S2PolygonBuilder.PointIndex.findNearbyPoint:
            // Cover the edge by a cap centered at the edge midpoint, then cover
            // the cap by four big-enough cells around the cell vertex closest to the
            // cap center.
            var middle = S2Point.Normalize((a + b) / 2);
            var actualLevel = Math.Min(idealLevel, S2CellId.MAX_LEVEL - 1);
            S2CellId.fromPoint(middle).getVertexNeighbors(actualLevel, edgeCovering);
            return actualLevel;
        }

        /**
   * Filters a list of entries down to the inclusive range defined by the given
   * cells, in <code>O(log N)</code> time.
   *
   * @param cell1 One side of the inclusive query range.
   * @param cell2 The other side of the inclusive query range.
   * @return An array of length 2, containing the start/end indices.
   */

        private int[] getEdges(ulong cell1, ulong cell2)
        {
            // ensure cell1 <= cell2
            if (cell1 > cell2)
            {
                var temp = cell1;
                cell1 = cell2;
                cell2 = temp;
            }
            // The binary search returns -N-1 to indicate an insertion point at index N,
            // if an exact match cannot be found. Since the edge indices queried for are
            // not valid edge indices, we will always get -N-1, so we immediately
            // convert to N.
            return new int[]
            {
                -1 - binarySearch(cell1, int.MinValue),
                -1 - binarySearch(cell2, int.MaxValue)
            };
        }

        private int binarySearch(ulong cell, int edge)
        {
            var low = 0;
            var high = cells.Length - 1;
            while (low <= high)
            {
                var mid = unchecked ((low + high) >> 1);
                var cmp = compare(cells[mid], edges[mid], cell, edge);
                if (cmp < 0)
                {
                    low = mid + 1;
                }
                else if (cmp > 0)
                {
                    high = mid - 1;
                }
                else
                {
                    return mid;
                }
            }
            return -(low + 1);
        }

        /**
   * Adds to candidateCrossings all the edges present in any ancestor of any
   * cell of cover, down to minimumS2LevelUsed. The cell->edge map is in the
   * variable mapping.
   */

        private void getEdgesInParentCells(List<S2CellId> cover, HashSet<int> candidateCrossings)
        {
            // Find all parent cells of covering cells.
            var parentCells = new HashSet<S2CellId>();
            foreach (var coverCell in cover)
            {
                for (var parentLevel = coverCell.level() - 1; parentLevel >= minimumS2LevelUsed;
                     --parentLevel)
                {
                    if (!parentCells.Add(coverCell.parent(parentLevel)))
                    {
                        break; // cell is already in => parents are too.
                    }
                }
            }

            // Put parent cell edge references into result.
            foreach (var parentCell in parentCells)
            {
                var bounds = getEdges(parentCell.id(), parentCell.id());
                for (var i = bounds[0]; i < bounds[1]; i++)
                {
                    candidateCrossings.Add(edges[i]);
                }
            }
        }

        /**
   * Returns true if ab possibly crosses cd, by clipping tiny angles to zero.
   */

        private static bool lenientCrossing(S2Point a, S2Point b, S2Point c, S2Point d)
        {
            // assert (S2.isUnitLength(a));
            // assert (S2.isUnitLength(b));
            // assert (S2.isUnitLength(c));

            var acb = S2Point.CrossProd(a, c).DotProd(b);
            var bda = S2Point.CrossProd(b, d).DotProd(a);
            if (Math.Abs(acb) < MAX_DET_ERROR || Math.Abs(bda) < MAX_DET_ERROR)
            {
                return true;
            }
            if (acb*bda < 0)
            {
                return false;
            }
            var cbd = S2Point.CrossProd(c, b).DotProd(d);
            var dac = S2Point.CrossProd(c, a).DotProd(c);
            if (Math.Abs(cbd) < MAX_DET_ERROR || Math.Abs(dac) < MAX_DET_ERROR)
            {
                return true;
            }
            return (acb*cbd >= 0) && (acb*dac >= 0);
        }

        /**
   * Returns true if the edge and the cell (including boundary) intersect.
   */

        private static bool edgeIntersectsCellBoundary(S2Point a, S2Point b, S2Cell cell)
        {
            var vertices = new S2Point[4];
            for (var i = 0; i < 4; ++i)
            {
                vertices[i] = cell.getVertex(i);
            }
            for (var i = 0; i < 4; ++i)
            {
                var fromPoint = vertices[i];
                var toPoint = vertices[(i + 1)%4];
                if (lenientCrossing(a, b, fromPoint, toPoint))
                {
                    return true;
                }
            }
            return false;
        }

        /**
   * Appends to candidateCrossings the edges that are fully contained in an S2
   * covering of edge. The covering of edge used is initially cover, but is
   * refined to eliminate quickly subcells that contain many edges but do not
   * intersect with edge.
   */

        private void getEdgesInChildrenCells(S2Point a, S2Point b, List<S2CellId> cover,
                                             HashSet<int> candidateCrossings)
        {
            // Put all edge references of (covering cells + descendant cells) into
            // result.
            // This relies on the natural ordering of S2CellIds.
            S2Cell[] children = null;
            while (cover.Any())
            {
                var cell = cover[cover.Count - 1];
                cover.RemoveAt(cover.Count - 1);
                var bounds = getEdges(cell.rangeMin().id(), cell.rangeMax().id());
                if (bounds[1] - bounds[0] <= 16)
                {
                    for (var i = bounds[0]; i < bounds[1]; i++)
                    {
                        candidateCrossings.Add(edges[i]);
                    }
                }
                else
                {
                    // Add cells at this level
                    bounds = getEdges(cell.id(), cell.id());
                    for (var i = bounds[0]; i < bounds[1]; i++)
                    {
                        candidateCrossings.Add(edges[i]);
                    }
                    // Recurse on the children -- hopefully some will be empty.
                    if (children == null)
                    {
                        children = new S2Cell[4];
                        for (var i = 0; i < 4; ++i)
                        {
                            children[i] = new S2Cell();
                        }
                    }
                    new S2Cell(cell).subdivide(children);
                    foreach (var child in children)
                    {
                        // TODO(user): Do the check for the four cells at once,
                        // as it is enough to check the four edges between the cells. At
                        // this time, we are checking 16 edges, 4 times too many.
                        //
                        // Note that given the guarantee of AppendCovering, it is enough
                        // to check that the edge intersect with the cell boundary as it
                        // cannot be fully contained in a cell.
                        if (edgeIntersectsCellBoundary(a, b, child))
                        {
                            cover.Add(child.id());
                        }
                    }
                }
            }
        }

        /*
   * An iterator on data edges that may cross a query edge (a,b). Create the
   * iterator, call getCandidates(), then hasNext()/next() repeatedly.
   *
   * The current edge in the iteration has index index(), goes between from()
   * and to().
   */

        public class DataEdgeIterator
        {
            /**
     * The structure containing the data edges.
     */
            private readonly List<int> candidates;
            private readonly S2EdgeIndex edgeIndex;

            /**
     * Tells whether getCandidates() obtained the candidates through brute force
     * iteration or using the quad tree structure.
     */

            /**
     * Index of the current edge and of the edge before the last next() call.
     */
            private int currentIndex;

            /**
     * Cache of edgeIndex.getNumEdges() so that hasNext() doesn't make an extra
     * call
     */

            /**
     * Index within array above. We have: currentIndex =
     * candidates.get(currentIndexInCandidates).
     */
            private int currentIndexInCandidates;
            private bool isBruteForce;
            private int numEdges;

            public DataEdgeIterator(S2EdgeIndex edgeIndex)
            {
                this.edgeIndex = edgeIndex;
                candidates = new List<int>();
            }

            /**
     * Initializes the iterator to iterate over a set of candidates that may
     * cross the edge (a,b).
     */

            public void getCandidates(S2Point a, S2Point b)
            {
                edgeIndex.predictAdditionalCalls(1);
                isBruteForce = !edgeIndex.isIndexComputed();
                if (isBruteForce)
                {
                    edgeIndex.incrementQueryCount();
                    currentIndex = 0;
                    numEdges = edgeIndex.getNumEdges();
                }
                else
                {
                    candidates.Clear();
                    edgeIndex.findCandidateCrossings(a, b, candidates);
                    currentIndexInCandidates = 0;
                    if (candidates.Any())
                    {
                        currentIndex = candidates[0];
                    }
                }
            }

            /**
     * Index of the current edge in the iteration.
     */

            public int index()
            {
                Preconditions.CheckState(hasNext());
                return currentIndex;
            }

            /**
     * False if there are no more candidates; true otherwise.
     */

            public bool hasNext()
            {
                if (isBruteForce)
                {
                    return (currentIndex < numEdges);
                }
                else
                {
                    return currentIndexInCandidates < candidates.Count;
                }
            }

            /**
     * Iterate to the next available candidate.
     */

            public void next()
            {
                Preconditions.CheckState(hasNext());
                if (isBruteForce)
                {
                    ++currentIndex;
                }
                else
                {
                    ++currentIndexInCandidates;
                    if (currentIndexInCandidates < candidates.Count)
                    {
                        currentIndex = candidates[currentIndexInCandidates];
                    }
                }
            }
        }
    }
}