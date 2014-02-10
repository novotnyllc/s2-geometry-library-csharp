using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Google.Common.Geometry
{
    /**
  * Normalizes the cell union by discarding cells that are contained by other
  * cells, replacing groups of 4 child cells by their parent cell whenever
  * possible, and sorting all the cell ids in increasing order. Returns true if
  * the number of cells was reduced.
  *
  *  This method *must* be called before doing any calculations on the cell
  * union, such as Intersects() or Contains().
  *
  * @return true if the normalize operation had any effect on the cell union,
  *         false if the union was already normalized
  */

    public sealed class S2CellUnion : IS2Region, IEnumerable<S2CellId>, IEquatable<S2CellUnion>
    {
        private List<S2CellId> _cellIds = new List<S2CellId>();

        public S2CellUnion()
        {
        }

        public int Count
        {
            get { return _cellIds.Count; }
        }

        public IList<S2CellId> CellIds
        {
            get { return _cellIds; }
        }

        public ulong LeafCellsCovered
        {
            get
            {
                ulong numLeaves = 0;
                foreach (var cellId in _cellIds)
                {
                    var invertedLevel = S2CellId.MaxLevel - cellId.Level;
                    numLeaves += (ulong)(1L << (invertedLevel << 1));
                }
                return numLeaves;
            }
        }


        /**
   * Approximate this cell union's area by summing the average area of
   * each contained cell's average area, using {@link S2Cell#averageArea()}.
   * This is equivalent to the number of leaves covered, multiplied by
   * the average area of a leaf.
   * Note that {@link S2Cell#averageArea()} does not take into account
   * distortion of cell, and thus may be off by up to a factor of 1.7.
   * NOTE: Since this is proportional to LeafCellsCovered(), it is
   * always better to use the other function if all you care about is
   * the relative average area between objects.
   *
   * @return the sum of the average area of each contained cell's average area
   */

        public double AverageBasedArea
        {
            get { return S2Cell.AverageArea(S2CellId.MaxLevel)*LeafCellsCovered; }
        }

        /**
   * Calculates this cell union's area by summing the approximate area for each
   * contained cell, using {@link S2Cell#approxArea()}.
   *
   * @return approximate area of the cell union
   */

        public double ApproxArea
        {
            get
            {
                double area = 0;
                foreach (var cellId in _cellIds)
                {
                    area += new S2Cell(cellId).ApproxArea();
                }
                return area;
            }
        }

        /**
   * Calculates this cell union's area by summing the exact area for each
   * contained cell, using the {@link S2Cell#exactArea()}.
   *
   * @return the exact area of the cell union
   */

        public double ExactArea
        {
            get
            {
                double area = 0;
                foreach (var cellId in _cellIds)
                {
                    area += new S2Cell(cellId).ExactArea();
                }
                return area;
            }
        }

        public IEnumerator<S2CellId> GetEnumerator()
        {
            return _cellIds.GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }

        public bool Equals(S2CellUnion other)
        {
            if (ReferenceEquals(null, other)) return false;
            if (ReferenceEquals(this, other)) return true;
            return _cellIds.SequenceEqual(other._cellIds);
        }

        public bool Contains(S2Cell cell)
        {
            return Contains(cell.Id);
        }

        public S2Cap CapBound
        {
            get
            {
                // Compute the approximate centroid of the region. This won't produce the
                // bounding cap of minimal area, but it should be close enough.
                if (_cellIds.Count == 0)
                {
                    return S2Cap.Empty;
                }
                var centroid = new S2Point(0, 0, 0);
                foreach (var id in this)
                {
                    var area = S2Cell.AverageArea(id.Level);
                    centroid = centroid + (id.ToPoint()*area);
                }
                if (centroid.Equals(new S2Point(0, 0, 0)))
                {
                    centroid = new S2Point(1, 0, 0);
                }
                else
                {
                    centroid = S2Point.Normalize(centroid);
                }

                // Use the centroid as the cap axis, and expand the cap angle so that it
                // contains the bounding caps of all the individual cells. Note that it is
                // *not* sufficient to just bound all the cell vertices because the bounding
                // cap may be concave (i.e. cover more than one hemisphere).
                var cap = S2Cap.FromAxisHeight(centroid, 0);
                foreach (var id in this)
                {
                    cap = cap.AddCap(new S2Cell(id).CapBound);
                }
                return cap;
            }
        }

        public S2LatLngRect RectBound
        {
            get
            {
                var bound = S2LatLngRect.Empty;
                foreach (var id in this)
                {
                    bound = bound.Union(new S2Cell(id).RectBound);
                }
                return bound;
            }
        }


        /** This is a fast operation (logarithmic in the size of the cell union). */

        public bool MayIntersect(S2Cell cell)
        {
            return Intersects(cell.Id);
        }

        public override bool Equals(object obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            if (ReferenceEquals(this, obj)) return true;
            if (obj.GetType() != GetType()) return false;
            return Equals((S2CellUnion)obj);
        }

        public override int GetHashCode()
        {
            var value = 17;
            foreach (var id in this)
            {
                value = 37*value + id.GetHashCode();
            }
            return value;
        }

        public static bool operator ==(S2CellUnion left, S2CellUnion right)
        {
            return Equals(left, right);
        }

        public static bool operator !=(S2CellUnion left, S2CellUnion right)
        {
            return !Equals(left, right);
        }

        /** The CellIds that form the Union */

        public void InitFromCellIds(IEnumerable<S2CellId> cellIds)
        {
            InitRawCellIds(cellIds);
            Normalize();
        }

        /**
   * Populates a cell union with the given S2CellIds or 64-bit cells ids, and
   * then calls Normalize(). The InitSwap() version takes ownership of the
   * vector data without copying and clears the given vector. These methods may
   * be called multiple times.
   */

        public void InitFromIds(IEnumerable<ulong> cellIds)
        {
            InitRawIds(cellIds);
            Normalize();
        }

        public void InitSwap(ICollection<S2CellId> cellIds)
        {
            InitRawSwap(cellIds);
            Normalize();
        }

        public void InitRawCellIds(IEnumerable<S2CellId> cellIds)
        {
            _cellIds = new List<S2CellId>(cellIds);
        }

        public void InitRawIds(IEnumerable<ulong> cellIds)
        {
            _cellIds = cellIds
                .Select(id => new S2CellId(id))
                .ToList();
        }

        /**
   * Like Init(), but does not call Normalize(). The cell union *must* be
   * normalized before doing any calculations with it, so it is the caller's
   * responsibility to make sure that the input is normalized. This method is
   * useful when converting cell unions to another representation and back.
   * These methods may be called multiple times.
   */

        public void InitRawSwap(ICollection<S2CellId> cellIds)
        {
            _cellIds = new List<S2CellId>(cellIds);
            cellIds.Clear();
        }

        /** Convenience methods for accessing the individual cell ids. */

        public S2CellId CellId(int i)
        {
            return _cellIds[i];
        }

        /** Direct access to the underlying vector for iteration . */

        /**
   * Replaces "output" with an expanded version of the cell union where any
   * cells whose level is less than "min_level" or where (level - min_level) is
   * not a multiple of "level_mod" are replaced by their children, until either
   * both of these conditions are satisfied or the maximum level is reached.
   *
   *  This method allows a covering generated by S2RegionCoverer using
   * min_level() or level_mod() constraints to be stored as a normalized cell
   * union (which allows various geometric computations to be done) and then
   * converted back to the original list of cell ids that satisfies the desired
   * constraints.
   */

        public void Denormalize(int minLevel, int levelMod, ICollection<S2CellId> output)
        {
            // assert (minLevel >= 0 && minLevel <= S2CellId.MAX_LEVEL);
            // assert (levelMod >= 1 && levelMod <= 3);

            output.Clear();
            foreach (var id in this)
            {
                var level = id.Level;
                var newLevel = Math.Max(minLevel, level);
                if (levelMod > 1)
                {
                    // Round up so that (new_level - min_level) is a multiple of level_mod.
                    // (Note that S2CellId::kMaxLevel is a multiple of 1, 2, and 3.)
                    newLevel += (S2CellId.MaxLevel - (newLevel - minLevel))%levelMod;
                    newLevel = Math.Min(S2CellId.MaxLevel, newLevel);
                }
                if (newLevel == level)
                {
                    output.Add(id);
                }
                else
                {
                    var end = id.ChildEndForLevel(newLevel);
                    for (var idInner = id.ChildBeginForLevel(newLevel); !idInner.Equals(end); idInner = idInner.Next)
                    {
                        output.Add(idInner);
                    }
                }
            }
        }

        /**
   * If there are more than "excess" elements of the cell_ids() vector that are
   * allocated but unused, reallocate the array to eliminate the excess space.
   * This reduces memory usage when many cell unions need to be held in memory
   * at once.
   */

        public void Pack()
        {
            _cellIds.TrimExcess();
        }


        /**
   * Return true if the cell union contains the given cell id. Containment is
   * defined with respect to regions, e.g. a cell contains its 4 children. This
   * is a fast operation (logarithmic in the size of the cell union).
   */

        public bool Contains(S2CellId id)
        {
            // This function requires that Normalize has been called first.
            //
            // This is an exact test. Each cell occupies a linear span of the S2
            // space-filling curve, and the cell id is simply the position at the center
            // of this span. The cell union ids are sorted in increasing order along
            // the space-filling curve. So we simply find the pair of cell ids that
            // surround the given cell id (using binary search). There is containment
            // if and only if one of these two cell ids contains this cell.

            var pos = _cellIds.BinarySearch(id);
            if (pos < 0)
            {
                pos = -pos - 1;
            }
            if (pos < _cellIds.Count && _cellIds[pos].RangeMin <= id)
            {
                return true;
            }
            return pos != 0 && _cellIds[pos - 1].RangeMax >= id;
        }

        /**
   * Return true if the cell union intersects the given cell id. This is a fast
   * operation (logarithmic in the size of the cell union).
   */

        public bool Intersects(S2CellId id)
        {
            // This function requires that Normalize has been called first.
            // This is an exact test; see the comments for Contains() above.
            var pos = _cellIds.BinarySearch(id);

            if (pos < 0)
            {
                pos = -pos - 1;
            }


            if (pos < _cellIds.Count && _cellIds[pos].RangeMin <= id.RangeMax)
            {
                return true;
            }
            return pos != 0 && _cellIds[pos - 1].RangeMax >= id.RangeMin;
        }

        public bool Contains(S2CellUnion that)
        {
            // TODO(kirilll?): A divide-and-conquer or alternating-skip-search approach
            // may be significantly faster in both the average and worst case.
            foreach (var id in that)
            {
                if (!Contains(id))
                {
                    return false;
                }
            }
            return true;
        }

        /** This is a fast operation (logarithmic in the size of the cell union). */

        /**
   * Return true if this cell union contain/intersects the given other cell
   * union.
   */

        public bool Intersects(S2CellUnion union)
        {
            // TODO(kirilll?): A divide-and-conquer or alternating-skip-search approach
            // may be significantly faster in both the average and worst case.
            foreach (var id in union)
            {
                if (Intersects(id))
                {
                    return true;
                }
            }
            return false;
        }

        public void GetUnion(S2CellUnion x, S2CellUnion y)
        {
            // assert (x != this && y != this);
            _cellIds.Clear();

            _cellIds.AddRange(x._cellIds);
            _cellIds.AddRange(y._cellIds);
            Normalize();
        }

        /**
   * Specialized version of GetIntersection() that gets the intersection of a
   * cell union with the given cell id. This can be useful for "splitting" a
   * cell union into chunks.
   */

        public void GetIntersection(S2CellUnion x, S2CellId id)
        {
            // assert (x != this);
            _cellIds.Clear();
            if (x.Contains(id))
            {
                _cellIds.Add(id);
            }
            else
            {
                var pos = x._cellIds.BinarySearch(id.RangeMin);

                if (pos < 0)
                {
                    pos = -pos - 1;
                }

                var idmax = id.RangeMax;
                var size = x._cellIds.Count;
                while (pos < size && x._cellIds[pos] <= idmax)
                {
                    _cellIds.Add(x._cellIds[pos++]);
                }
            }
        }

        /**
   * Initialize this cell union to the union or intersection of the two given
   * cell unions. Requires: x != this and y != this.
   */

        public void GetIntersection(S2CellUnion x, S2CellUnion y)
        {
            // assert (x != this && y != this);

            // This is a fairly efficient calculation that uses binary search to skip
            // over sections of both input vectors. It takes constant time if all the
            // cells of "x" come before or after all the cells of "y" in S2CellId order.

            _cellIds.Clear();

            var i = 0;
            var j = 0;

            while (i < x._cellIds.Count && j < y._cellIds.Count)
            {
                var imin = x.CellId(i).RangeMin;
                var jmin = y.CellId(j).RangeMin;
                if (imin > jmin)
                {
                    // Either j->contains(*i) or the two cells are disjoint.
                    if (x.CellId(i) <= y.CellId(j).RangeMax)
                    {
                        _cellIds.Add(x.CellId(i++));
                    }
                    else
                    {
                        // Advance "j" to the first cell possibly contained by *i.
                        j = IndexedBinarySearch(y._cellIds, imin, j + 1);
                        // The previous cell *(j-1) may now contain *i.
                        if (x.CellId(i) <= y.CellId(j - 1).RangeMax)
                        {
                            --j;
                        }
                    }
                }
                else if (jmin >= imin)
                {
                    // Identical to the code above with "i" and "j" reversed.
                    if (y.CellId(j) <= x.CellId(i).RangeMax)
                    {
                        _cellIds.Add(y.CellId(j++));
                    }
                    else
                    {
                        i = IndexedBinarySearch(x._cellIds, jmin, i + 1);
                        if (y.CellId(j) <= x.CellId(i - 1).RangeMax)
                        {
                            --i;
                        }
                    }
                }
                else
                {
                    // "i" and "j" have the same range_min(), so one contains the other.
                    if (x.CellId(i) < y.CellId(j))
                    {
                        _cellIds.Add(x.CellId(i++));
                    }
                    else
                    {
                        _cellIds.Add(y.CellId(j++));
                    }
                }
            }
            // The output is generated in sorted order, and there should not be any
            // cells that can be merged (provided that both inputs were normalized).
            // assert (!normalize());
        }

        /**
   * Just as normal binary search, except that it allows specifying the starting
   * value for the lower bound.
   *
   * @return The position of the searched element in the list (if found), or the
   *         position where the element could be inserted without violating the
   *         order.
   */

        private static int IndexedBinarySearch(IReadOnlyList<S2CellId> list, S2CellId key, int low)
        {
            var high = list.Count - 1;

            while (low <= high)
            {
                var mid = (low + high) >> 1;
                var midVal = list[mid];
                var cmp = midVal.CompareTo(key);

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
                    return mid; // key found
                }
            }
            return low; // key not found
        }

        /**
   * Expands the cell union such that it contains all cells of the given level
   * that are adjacent to any cell of the original union. Two cells are defined
   * as adjacent if their boundaries have any points in common, i.e. most cells
   * have 8 adjacent cells (not counting the cell itself).
   *
   *  Note that the size of the output is exponential in "level". For example,
   * if level == 20 and the input has a cell at level 10, there will be on the
   * order of 4000 adjacent cells in the output. For most applications the
   * Expand(min_fraction, min_distance) method below is easier to use.
   */

        public void Expand(int level)
        {
            var output = new List<S2CellId>();
            var levelLsb = S2CellId.LowestOnBitForLevel(level);
            var i = Count - 1;
            do
            {
                var id = CellId(i);
                if (id.LowestOnBit < levelLsb)
                {
                    id = id.ParentForLevel(level);
                    // Optimization: skip over any cells contained by this one. This is
                    // especially important when very small regions are being expanded.
                    while (i > 0 && id.Contains(CellId(i - 1)))
                    {
                        --i;
                    }
                }
                output.Add(id);
                id.GetAllNeighbors(level, output);
            } while (--i >= 0);
            InitSwap(output);
        }

        /**
   * Expand the cell union such that it contains all points whose distance to
   * the cell union is at most minRadius, but do not use cells that are more
   * than maxLevelDiff levels higher than the largest cell in the input. The
   * second parameter controls the tradeoff between accuracy and output size
   * when a large region is being expanded by a small amount (e.g. expanding
   * Canada by 1km).
   *
   *  For example, if maxLevelDiff == 4, the region will always be expanded by
   * approximately 1/16 the width of its largest cell. Note that in the worst
   * case, the number of cells in the output can be up to 4 * (1 + 2 **
   * maxLevelDiff) times larger than the number of cells in the input.
   */

        public void Expand(S1Angle minRadius, int maxLevelDiff)
        {
            var minLevel = S2CellId.MaxLevel;
            foreach (var id in this)
            {
                minLevel = Math.Min(minLevel, id.Level);
            }
            // Find the maximum level such that all cells are at least "min_radius"
            // wide.
            var radiusLevel = S2Projections.MIN_WIDTH.GetMaxLevel(minRadius.Radians);
            if (radiusLevel == 0 && minRadius.Radians > S2Projections.MIN_WIDTH.GetValue(0))
            {
                // The requested expansion is greater than the width of a face cell.
                // The easiest way to handle this is to expand twice.
                Expand(0);
            }
            Expand(Math.Min(minLevel + maxLevelDiff, radiusLevel));
        }


        public IS2Region Clone()
        {
            var copy = new S2CellUnion();
            copy.InitRawCellIds(_cellIds);
            return copy;
        }

        /**
   * The point 'p' does not need to be normalized. This is a fast operation
   * (logarithmic in the size of the cell union).
   */

        public bool Contains(S2Point p)
        {
            return Contains(S2CellId.FromPoint(p));
        }

        /**
   * The number of leaf cells covered by the union.
   * This will be no more than 6*2^60 for the whole sphere.
   *
   * @return the number of leaf cells covered by the union
   */


        /**
   * Normalizes the cell union by discarding cells that are contained by other
   * cells, replacing groups of 4 child cells by their parent cell whenever
   * possible, and sorting all the cell ids in increasing order. Returns true if
   * the number of cells was reduced.
   *
   *  This method *must* be called before doing any calculations on the cell
   * union, such as Intersects() or Contains().
   *
   * @return true if the normalize operation had any effect on the cell union,
   *         false if the union was already normalized
   */

        public bool Normalize()
        {
            // Optimize the representation by looking for cases where all subcells
            // of a parent cell are present.

            var output = new List<S2CellId>(_cellIds.Count);
            _cellIds.Sort();


            foreach (var idLoop in this)
            {
                var id = idLoop;
                var sze = output.Count;
                // Check whether this cell is contained by the previous cell.
                if (output.Any() && output[sze - 1].Contains(id))
                {
                    continue;
                }

                // Discard any previous cells contained by this cell.
                while (output.Any() && id.Contains(output[output.Count - 1]))
                {
                    output.RemoveAt(output.Count - 1);
                }

                // Check whether the last 3 elements of "output" plus "id" can be
                // collapsed into a single parent cell.
                while (output.Count >= 3)
                {
                    sze = output.Count;
                    // A necessary (but not sufficient) condition is that the XOR of the
                    // four cells must be zero. This is also very fast to test.
                    if ((output[sze - 3].Id ^ output[sze - 2].Id ^ output[sze - 1].Id)
                        != id.Id)
                    {
                        break;
                    }

                    // Now we do a slightly more expensive but exact test. First, compute a
                    // mask that blocks out the two bits that encode the child position of
                    // "id" with respect to its parent, then check that the other three
                    // children all agree with "mask.
                    var mask = id.LowestOnBit << 1;
                    mask = ~(mask + (mask << 1));
                    var idMasked = (id.Id & mask);
                    if ((output[sze - 3].Id & mask) != idMasked
                        || (output[sze - 2].Id & mask) != idMasked
                        || (output[sze - 1].Id & mask) != idMasked || id.IsFace)
                    {
                        break;
                    }

                    // Replace four children by their parent cell.
                    output.RemoveAt(sze - 1);
                    output.RemoveAt(sze - 2);
                    output.RemoveAt(sze - 3);
                    id = id.Parent;
                }
                output.Add(id);
            }
            if (output.Count < Count)
            {
                InitRawSwap(output);
                return true;
            }
            return false;
        }
    }
}