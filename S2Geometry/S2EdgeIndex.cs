using System;
using System.Collections;
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
        private const double Thickening = 0.01;

        /**
   * Threshold for small angles, that help lenientCrossing to determine whether
   * two edges are likely to intersect.
   */
        private const double MaxDetError = 1e-14;

        /**
   * The cell containing each edge, as given in the parallel array
   * <code>edges</code>.
   */
        private ulong[] _cells;

        /**
   * The edge contained by each cell, as given in the parallel array
   * <code>cells</code>.
   */
        private int[] _edges;

        /**
   * No cell strictly below this level appears in mapping. Initially leaf level,
   * that's the minimum level at which we will ever look for test edges.
   */

        /**
   * Has the index been computed already?
   */
        private bool _indexComputed;
        private int _minimumS2LevelUsed;

        /**
   * Number of queries so far
   */
        private int _queryCount;

        /**
   * Empties the index in case it already contained something.
   */

        public void Reset()
        {
            _minimumS2LevelUsed = S2CellId.MaxLevel;
            _indexComputed = false;
            _queryCount = 0;
            _cells = null;
            _edges = null;
        }

        /**
   * Compares [cell1, edge1] to [cell2, edge2], by cell first and edge second.
   *
   * @return -1 if [cell1, edge1] is less than [cell2, edge2], 1 if [cell1,
   *         edge1] is greater than [cell2, edge2], 0 otherwise.
   */

        private static int Compare(ulong cell1, int edge1, ulong cell2, int edge2)
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

        public void ComputeIndex()
        {
            if (_indexComputed)
            {
                return;
            }
            var cellList = new List<ulong>();
            var edgeList = new List<int>();
            for (var i = 0; i < NumEdges; ++i)
            {
                var from = EdgeFrom(i);
                var to = EdgeTo(i);
                var cover = new List<S2CellId>();
                var level = GetCovering(from, to, true, cover);
                _minimumS2LevelUsed = Math.Min(_minimumS2LevelUsed, level);
                foreach (var cellId in cover)
                {
                    cellList.Add(cellId.Id);
                    edgeList.Add(i);
                }
            }
            _cells = new ulong[cellList.Count];
            _edges = new int[edgeList.Count];
            for (var i = 0; i < _cells.Length; i++)
            {
                _cells[i] = cellList[i];
                _edges[i] = edgeList[i];
            }
            SortIndex();
            _indexComputed = true;
        }

        /** Sorts the parallel <code>cells</code> and <code>edges</code> arrays. */

        private void SortIndex()
        {
            // create an array of indices and sort based on the values in the parallel
            // arrays at each index
            var indices = new int[_cells.Length];
            for (var i = 0; i < indices.Length; i++)
            {
                indices[i] = i;
            }

            Array.Sort(indices, (index1, index2) => Compare(_cells[index1], _edges[index1], _cells[index2], _edges[index2]));


            // copy the cells and edges in the order given by the sorted list of indices
            var newCells = new ulong[_cells.Length];
            var newEdges = new int[_edges.Length];
            for (var i = 0; i < indices.Length; i++)
            {
                newCells[i] = _cells[indices[i]];
                newEdges[i] = _edges[indices[i]];
            }
            // replace the cells and edges with the sorted arrays
            _cells = newCells;
            _edges = newEdges;
        }

        public bool IsIndexComputed
        {
            get { return _indexComputed; }
        }

        /**
   * Tell the index that we just received a new request for candidates. Useful
   * to compute when to switch to quad tree.
   */

        protected void IncrementQueryCount()
        {
            ++_queryCount;
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

        public void PredictAdditionalCalls(int n)
        {
            if (_indexComputed)
            {
                return;
            }
            if (NumEdges > 100 && (_queryCount + n) > 30)
            {
                ComputeIndex();
            }
        }

        /**
   * Overwrite these functions to give access to the underlying data. The
   * function getNumEdges() returns the number of edges in the index, while
   * edgeFrom(index) and edgeTo(index) return the "from" and "to" endpoints of
   * the edge at the given index.
   */
        protected abstract int NumEdges { get; }

        protected abstract S2Point EdgeFrom(int index);

        protected abstract S2Point EdgeTo(int index);

        /**
   * Appends to "candidateCrossings" all edge references which may cross the
   * given edge. This is done by covering the edge and then finding all
   * references of edges whose coverings overlap this covering. Parent cells are
   * checked level by level. Child cells are checked all at once by taking
   * advantage of the natural ordering of S2CellIds.
   */

        protected void FindCandidateCrossings(S2Point a, S2Point b, IList<int> candidateCrossings)
        {
            Preconditions.CheckState(_indexComputed);
            var cover = new List<S2CellId>();
            GetCovering(a, b, false, cover);

            // Edge references are inserted into the map once for each covering cell, so
            // absorb duplicates here
            var uniqueSet = new HashSet<int>();
            GetEdgesInParentCells(cover, uniqueSet);

            // TODO(user): An important optimization for long query
            // edges (Contains queries): keep a bounding cap and clip the query
            // edge to the cap before starting the descent.
            GetEdgesInChildrenCells(a, b, cover, uniqueSet);

            candidateCrossings.Clear();

            foreach (var item in uniqueSet)
                candidateCrossings.Add(item);
        }

        /**
   * Returns the smallest cell containing all four points, or
   * {@link S2CellId#sentinel()} if they are not all on the same face. The
   * points don't need to be normalized.
   */

        private static S2CellId ContainingCell(S2Point pa, S2Point pb, S2Point pc, S2Point pd)
        {
            var a = S2CellId.FromPoint(pa);
            var b = S2CellId.FromPoint(pb);
            var c = S2CellId.FromPoint(pc);
            var d = S2CellId.FromPoint(pd);

            if (a.Face != b.Face || a.Face != c.Face || a.Face != d.Face)
            {
                return S2CellId.Sentinel;
            }

            while (!a.Equals(b) || !a.Equals(c) || !a.Equals(d))
            {
                a = a.Parent;
                b = b.Parent;
                c = c.Parent;
                d = d.Parent;
            }
            return a;
        }

        /**
   * Returns the smallest cell containing both points, or Sentinel if they are
   * not all on the same face. The points don't need to be normalized.
   */

        private static S2CellId ContainingCell(S2Point pa, S2Point pb)
        {
            var a = S2CellId.FromPoint(pa);
            var b = S2CellId.FromPoint(pb);

            if (a.Face != b.Face)
            {
                return S2CellId.Sentinel;
            }

            while (!a.Equals(b))
            {
                a = a.Parent;
                b = b.Parent;
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

        private int GetCovering(
            S2Point a, S2Point b, bool thickenEdge, List<S2CellId> edgeCovering)
        {
            edgeCovering.Clear();

            // Selects the ideal s2 level at which to cover the edge, this will be the
            // level whose S2 cells have a width roughly commensurate to the length of
            // the edge. We multiply the edge length by 2*THICKENING to guarantee the
            // thickening is honored (it's not a big deal if we honor it when we don't
            // request it) when doing the covering-by-cap trick.
            var edgeLength = a.Angle(b);
            var idealLevel = S2Projections.MinWidth.GetMaxLevel(edgeLength*(1 + 2*Thickening));

            S2CellId containingCellId;
            if (!thickenEdge)
            {
                containingCellId = ContainingCell(a, b);
            }
            else
            {
                if (idealLevel == S2CellId.MaxLevel)
                {
                    // If the edge is tiny, instabilities are more likely, so we
                    // want to limit the number of operations.
                    // We pretend we are in a cell much larger so as to trigger the
                    // 'needs covering' case, so we won't try to thicken the edge.
                    containingCellId = (new S2CellId(0xFFF0)).ParentForLevel(3);
                }
                else
                {
                    var pq = (b - a)*Thickening;
                    var ortho = (S2Point.Normalize(S2Point.CrossProd(pq, a)))*edgeLength*Thickening;
                    var p = a - pq;
                    var q = b + pq;
                    // If p and q were antipodal, the edge wouldn't be lengthened,
                    // and it could even flip! This is not a problem because
                    // idealLevel != 0 here. The farther p and q can be is roughly
                    // a quarter Earth away from each other, so we remain
                    // Theta(THICKENING).
                    containingCellId = ContainingCell(p - ortho, p + ortho, q - ortho, q + ortho);
                }
            }

            // Best case: edge is fully contained in a cell that's not too big.
            if (!containingCellId.Equals(S2CellId.Sentinel)
                && containingCellId.Level >= idealLevel - 2)
            {
                edgeCovering.Add(containingCellId);
                return containingCellId.Level;
            }

            if (idealLevel == 0)
            {
                // Edge is very long, maybe even longer than a face width, so the
                // trick below doesn't work. For now, we will add the whole S2 sphere.
                // TODO(user): Do something a tad smarter (and beware of the
                // antipodal case).
                for (var cellid = S2CellId.Begin(0); !cellid.Equals(S2CellId.End(0));
                     cellid = cellid.Next)
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
            var actualLevel = Math.Min(idealLevel, S2CellId.MaxLevel - 1);
            S2CellId.FromPoint(middle).GetVertexNeighbors(actualLevel, edgeCovering);
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

        private int[] GetEdges(ulong cell1, ulong cell2)
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
                -1 - BinarySearch(cell1, int.MinValue),
                -1 - BinarySearch(cell2, int.MaxValue)
            };
        }

        private int BinarySearch(ulong cell, int edge)
        {
            var low = 0;
            var high = _cells.Length - 1;
            while (low <= high)
            {
                var mid = unchecked ((low + high) >> 1);
                var cmp = Compare(_cells[mid], _edges[mid], cell, edge);
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

        private void GetEdgesInParentCells(IEnumerable<S2CellId> cover, ISet<int> candidateCrossings)
        {
            // Find all parent cells of covering cells.
            var parentCells = new HashSet<S2CellId>();
            foreach (var coverCell in cover)
            {
                for (var parentLevel = coverCell.Level - 1; parentLevel >= _minimumS2LevelUsed;
                     --parentLevel)
                {
                    if (!parentCells.Add(coverCell.ParentForLevel(parentLevel)))
                    {
                        break; // cell is already in => parents are too.
                    }
                }
            }

            // Put parent cell edge references into result.
            foreach (var parentCell in parentCells)
            {
                var bounds = GetEdges(parentCell.Id, parentCell.Id);
                for (var i = bounds[0]; i < bounds[1]; i++)
                {
                    candidateCrossings.Add(_edges[i]);
                }
            }
        }

        /**
   * Returns true if ab possibly crosses cd, by clipping tiny angles to zero.
   */

        private static bool LenientCrossing(S2Point a, S2Point b, S2Point c, S2Point d)
        {
            // assert (S2.isUnitLength(a));
            // assert (S2.isUnitLength(b));
            // assert (S2.isUnitLength(c));

            var acb = S2Point.CrossProd(a, c).DotProd(b);
            var bda = S2Point.CrossProd(b, d).DotProd(a);
            if (Math.Abs(acb) < MaxDetError || Math.Abs(bda) < MaxDetError)
            {
                return true;
            }
            if (acb*bda < 0)
            {
                return false;
            }
            var cbd = S2Point.CrossProd(c, b).DotProd(d);
            var dac = S2Point.CrossProd(c, a).DotProd(c);
            if (Math.Abs(cbd) < MaxDetError || Math.Abs(dac) < MaxDetError)
            {
                return true;
            }
            return (acb*cbd >= 0) && (acb*dac >= 0);
        }

        /**
   * Returns true if the edge and the cell (including boundary) intersect.
   */

        private static bool EdgeIntersectsCellBoundary(S2Point a, S2Point b, S2Cell cell)
        {
            var vertices = new S2Point[4];
            for (var i = 0; i < 4; ++i)
            {
                vertices[i] = cell.GetVertex(i);
            }
            for (var i = 0; i < 4; ++i)
            {
                var fromPoint = vertices[i];
                var toPoint = vertices[(i + 1)%4];
                if (LenientCrossing(a, b, fromPoint, toPoint))
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

        private void GetEdgesInChildrenCells(S2Point a, S2Point b, IList<S2CellId> cover,
                                             ISet<int> candidateCrossings)
        {
            // Put all edge references of (covering cells + descendant cells) into
            // result.
            // This relies on the natural ordering of S2CellIds.
            S2Cell[] children = null;
            while (cover.Any())
            {
                var cell = cover[cover.Count - 1];
                cover.RemoveAt(cover.Count - 1);
                var bounds = GetEdges(cell.RangeMin.Id, cell.RangeMax.Id);
                if (bounds[1] - bounds[0] <= 16)
                {
                    for (var i = bounds[0]; i < bounds[1]; i++)
                    {
                        candidateCrossings.Add(_edges[i]);
                    }
                }
                else
                {
                    // Add cells at this level
                    bounds = GetEdges(cell.Id, cell.Id);
                    for (var i = bounds[0]; i < bounds[1]; i++)
                    {
                        candidateCrossings.Add(_edges[i]);
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
                    new S2Cell(cell).Subdivide(children);
                    foreach (var child in children)
                    {
                        // TODO(user): Do the check for the four cells at once,
                        // as it is enough to check the four edges between the cells. At
                        // this time, we are checking 16 edges, 4 times too many.
                        //
                        // Note that given the guarantee of AppendCovering, it is enough
                        // to check that the edge intersect with the cell boundary as it
                        // cannot be fully contained in a cell.
                        if (EdgeIntersectsCellBoundary(a, b, child))
                        {
                            cover.Add(child.Id);
                        }
                    }
                }
            }
        }

        /// <summary>
        /// An iterator on data edges that may cross a query edge (a,b). Create the
        /// iterator, call getCandidates(), then enumerating.
         ///
       /// The current edge in the iteration has index, goes between from()
       /// and to().
        /// </summary>
        /// <returns></returns>
        public DataEdgeIterator GetIterator()
        {
            return new DataEdgeIterator(this);
        }

        

        /*
   * An iterator on data edges that may cross a query edge (a,b). Create the
   * iterator, call getCandidates(), then hasNext()/next() repeatedly.
   *
   * The current edge in the iteration has index index(), goes between from()
   * and to().
   */

        public sealed class DataEdgeIterator : IEnumerable<int>
        {
            /**
     * The structure containing the data edges.
     */
            private readonly List<int> _candidates;
            private readonly S2EdgeIndex _edgeIndex;

            /**
     * Tells whether getCandidates() obtained the candidates through brute force
     * iteration or using the quad tree structure.
     */

            /**
     * Index of the current edge and of the edge before the last next() call.
     */
            private int _currentIndex;

            /**
     * Cache of edgeIndex.getNumEdges() so that hasNext() doesn't make an extra
     * call
     */

            /**
     * Index within array above. We have: currentIndex =
     * candidates.get(currentIndexInCandidates).
     */
            private int _currentIndexInCandidates;
            private bool _isBruteForce;
            private int _numEdges;

            public DataEdgeIterator(S2EdgeIndex edgeIndex)
            {
                this._edgeIndex = edgeIndex;
                _candidates = new List<int>();
            }

            /**
     * Initializes the iterator to iterate over a set of candidates that may
     * cross the edge (a,b).
     */

            public void GetCandidates(S2Point a, S2Point b)
            {
                _edgeIndex.PredictAdditionalCalls(1);
                _isBruteForce = !_edgeIndex.IsIndexComputed;
                if (_isBruteForce)
                {
                    _edgeIndex.IncrementQueryCount();
                    _currentIndex = 0;
                    _numEdges = _edgeIndex.NumEdges;
                }
                else
                {
                    _candidates.Clear();
                    _edgeIndex.FindCandidateCrossings(a, b, _candidates);
                    _currentIndexInCandidates = 0;
                    if (_candidates.Any())
                    {
                        _currentIndex = _candidates[0];
                    }
                }

            }

            public IEnumerator<int> GetEnumerator()
            {
                while (_isBruteForce ? _currentIndex < _numEdges : _currentIndexInCandidates < _candidates.Count)
                {
                    yield return _currentIndex;

                    if (_isBruteForce)
                    {
                        ++_currentIndex;
                    }
                    else
                    {
                        ++_currentIndexInCandidates;
                        if (_currentIndexInCandidates < _candidates.Count)
                        {
                            _currentIndex = _candidates[_currentIndexInCandidates];
                        }
                    }
                }
            }

            IEnumerator IEnumerable.GetEnumerator()
            {
                return GetEnumerator();
            }
        }
    }
}