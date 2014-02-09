using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Google.Common.Geometry;
using NUnit.Framework;

namespace S2Geometry.Tests
{
    public class S2RegionCovererTest : GeometryTestCase
    {
        public void checkCovering(
            S2RegionCoverer coverer, IS2Region region, List<S2CellId> covering, bool interior)
        {
            // Keep track of how many cells have the same coverer.min_level() ancestor.
            IDictionary<S2CellId, int> minLevelCells = new Dictionary<S2CellId, int>();
            for (var i = 0; i < covering.Count; ++i)
            {
                var level = covering[i].level();
                assertTrue(level >= coverer.minLevel());
                assertTrue(level <= coverer.maxLevel());
                assertEquals((level - coverer.minLevel())%coverer.levelMod(), 0);
                var key = covering[i].parent(coverer.minLevel());
                if (!minLevelCells.ContainsKey(key))
                {
                    minLevelCells.Add(key, 1);
                }
                else
                {
                    minLevelCells[key] = minLevelCells[key] + 1;
                }
            }
            if (covering.Count > coverer.maxCells())
            {
                // If the covering has more than the requested number of cells, then check
                // that the cell count cannot be reduced by using the parent of some cell.
                foreach (var i in minLevelCells.Values)
                {
                    assertEquals(i, 1);
                }
            }

            if (interior)
            {
                for (var i = 0; i < covering.Count; ++i)
                {
                    assertTrue(region.Contains(new S2Cell(covering[i])));
                }
            }
            else
            {
                var cellUnion = new S2CellUnion();
                cellUnion.initFromCellIds(covering);
                checkCovering(region, cellUnion, true, new S2CellId());
            }
        }

        [Test]
        public void testRandomCaps()
        {
            Console.WriteLine("TestRandomCaps");

            var kMaxLevel = S2CellId.MAX_LEVEL;
            var coverer = new S2RegionCoverer();
            for (var i = 0; i < 1000; ++i)
            {
                do
                {
                    coverer.setMinLevel(random(kMaxLevel + 1));
                    coverer.setMaxLevel(random(kMaxLevel + 1));
                } while (coverer.minLevel() > coverer.maxLevel());
                coverer.setMaxCells(skewed(10));
                coverer.setLevelMod(1 + random(3));
                var maxArea = Math.Min(
                    4*S2.Pi, (3*coverer.maxCells() + 1)*S2Cell.AverageArea(coverer.minLevel()));
                var cap = getRandomCap(0.1*S2Cell.AverageArea(kMaxLevel), maxArea);
                var covering = new List<S2CellId>();
                var interior = new List<S2CellId>();

                coverer.getCovering(cap, covering);
                checkCovering(coverer, cap, covering, false);

                coverer.getInteriorCovering(cap, interior);
                checkCovering(coverer, cap, interior, true);


                // Check that GetCovering is deterministic.
                var covering2 = new List<S2CellId>();
                coverer.getCovering(cap, covering2);
                assertTrue(covering.SequenceEqual(covering2));

                // Also check S2CellUnion.denormalize(). The denormalized covering
                // may still be different and smaller than "covering" because
                // S2RegionCoverer does not guarantee that it will not output all four
                // children of the same parent.
                var cells = new S2CellUnion();
                cells.initFromCellIds(covering);
                var denormalized = new List<S2CellId>();
                cells.denormalize(coverer.minLevel(), coverer.levelMod(), denormalized);
                checkCovering(coverer, cap, denormalized, false);
            }
        }

        [Test]
        public void testRandomCells()
        {
            Console.WriteLine("TestRandomCells");

            var coverer = new S2RegionCoverer();
            coverer.setMaxCells(1);

            // Test random cell ids at all levels.
            for (var i = 0; i < 10000; ++i)
            {
                var id = getRandomCellId();
                var covering = new S2CellUnion();
                coverer.getCovering(new S2Cell(id), covering.cellIds());
                assertEquals(covering.size(), 1);
                assertEquals(covering.cellId(0), id);
            }
        }

        [Test]
        public void testSimpleCoverings()
        {
            Console.WriteLine("TestSimpleCoverings");

            var kMaxLevel = S2CellId.MAX_LEVEL;
            var coverer = new S2RegionCoverer();
            coverer.setMaxCells(int.MaxValue);
            for (var i = 0; i < 1000; ++i)
            {
                var level = random(kMaxLevel + 1);
                coverer.setMinLevel(level);
                coverer.setMaxLevel(level);
                var maxArea = Math.Min(4*S2.Pi, 1000*S2Cell.AverageArea(level));
                var cap = getRandomCap(0.1*S2Cell.AverageArea(kMaxLevel), maxArea);
                var covering = new List<S2CellId>();
                S2RegionCoverer.getSimpleCovering(cap, cap.Axis, level, covering);
                checkCovering(coverer, cap, covering, false);
            }
        }
    }
}