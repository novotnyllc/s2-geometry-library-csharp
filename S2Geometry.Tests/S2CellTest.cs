using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Google.Common.Geometry;
using NUnit.Framework;

namespace S2Geometry.Tests
{
    public class S2CellTest : GeometryTestCase
    {
        public const bool DEBUG_MODE = true;

        [Test]
  public void testFaces() {
    IDictionary<S2Point, int> edgeCounts = new Dictionary<S2Point, int>();
    IDictionary<S2Point, int> vertexCounts = new Dictionary<S2Point, int>();
    for (int face = 0; face < 6; ++face) {
      S2CellId id = S2CellId.fromFacePosLevel(face, 0, 0);
      S2Cell cell = new S2Cell(id);
      JavaAssert.Equal(cell.id(), id);
      JavaAssert.Equal(cell.face(), face);
      JavaAssert.Equal(cell.level(), 0);
      // Top-level faces have alternating orientations to get RHS coordinates.
      JavaAssert.Equal(cell.orientation(), face & S2.SWAP_MASK);
      Assert.True(!cell.isLeaf());
      for (int k = 0; k < 4; ++k) {
        if (edgeCounts.ContainsKey(cell.getEdgeRaw(k))) {
          edgeCounts[cell.getEdgeRaw(k)] = edgeCounts[cell
            .getEdgeRaw(k)] + 1;
        } else {
          edgeCounts[cell.getEdgeRaw(k)]= 1;
        }

        if (vertexCounts.ContainsKey(cell.getVertexRaw(k))) {
          vertexCounts[cell.getVertexRaw(k)] = vertexCounts[cell
            .getVertexRaw(k)] + 1;
        } else {
          vertexCounts[cell.getVertexRaw(k)]= 1;
        }
        assertDoubleNear(cell.getVertexRaw(k).dotProd(cell.getEdgeRaw(k)), 0);
        assertDoubleNear(cell.getVertexRaw((k + 1) & 3).dotProd(
          cell.getEdgeRaw(k)), 0);
        assertDoubleNear(S2Point.normalize(
          S2Point.crossProd(cell.getVertexRaw(k), cell
            .getVertexRaw((k + 1) & 3))).dotProd(cell.getEdge(k)), 1.0);
      }
    }
    // Check that edges have multiplicity 2 and vertices have multiplicity 3.
    foreach (int i in edgeCounts.Values) {
      JavaAssert.Equal(i, 2);
    }
    foreach (int i in vertexCounts.Values) {
      JavaAssert.Equal(i, 3);
    }
  }

  class LevelStats {
    public double count;
    public double minArea, maxArea, avgArea;
    public double minWidth, maxWidth, avgWidth;
    public double minEdge, maxEdge, avgEdge, maxEdgeAspect;
    public double minDiag, maxDiag, avgDiag, maxDiagAspect;
    public double minAngleSpan, maxAngleSpan, avgAngleSpan;
    public double minApproxRatio, maxApproxRatio;

    public LevelStats() {
      count = 0;
      minArea = 100;
      maxArea = 0;
      avgArea = 0;
      minWidth = 100;
      maxWidth = 0;
      avgWidth = 0;
      minEdge = 100;
      maxEdge = 0;
      avgEdge = 0;
      maxEdgeAspect = 0;
      minDiag = 100;
      maxDiag = 0;
      avgDiag = 0;
      maxDiagAspect = 0;
      minAngleSpan = 100;
      maxAngleSpan = 0;
      avgAngleSpan = 0;
      minApproxRatio = 100;
      maxApproxRatio = 0;
    }
  }

  static List<LevelStats> levelStats = new List<LevelStats>(
    S2CellId.MAX_LEVEL + 1);

  static S2CellTest() {
    for (int i = 0; i < S2CellId.MAX_LEVEL + 1; ++i) {
      levelStats.Add(new LevelStats());
    }
  }

  static void gatherStats(S2Cell cell) {
    LevelStats s = levelStats[cell.level()];
    double exactArea = cell.exactArea();
    double approxArea = cell.approxArea();
    double minEdge = 100, maxEdge = 0, avgEdge = 0;
    double minDiag = 100, maxDiag = 0;
    double minWidth = 100, maxWidth = 0;
    double minAngleSpan = 100, maxAngleSpan = 0;
    for (int i = 0; i < 4; ++i) {
      double edge = cell.getVertexRaw(i).angle(cell.getVertexRaw((i + 1) & 3));
      minEdge = Math.Min(edge, minEdge);
      maxEdge = Math.Max(edge, maxEdge);
      avgEdge += 0.25 * edge;
      S2Point mid = S2Point.add(cell.getVertexRaw(i), cell
        .getVertexRaw((i + 1) & 3));
      double width = S2.M_PI_2 - mid.angle(cell.getEdgeRaw(i ^ 2));
      minWidth = Math.Min(width, minWidth);
      maxWidth = Math.Max(width, maxWidth);
      if (i < 2) {
        double diag = cell.getVertexRaw(i).angle(cell.getVertexRaw(i ^ 2));
        minDiag = Math.Min(diag, minDiag);
        maxDiag = Math.Max(diag, maxDiag);
        double angleSpan = cell.getEdgeRaw(i).angle(
          S2Point.neg(cell.getEdgeRaw(i ^ 2)));
        minAngleSpan = Math.Min(angleSpan, minAngleSpan);
        maxAngleSpan = Math.Max(angleSpan, maxAngleSpan);
      }
    }
    s.count += 1;
    s.minArea = Math.Min(exactArea, s.minArea);
    s.maxArea = Math.Max(exactArea, s.maxArea);
    s.avgArea += exactArea;
    s.minWidth = Math.Min(minWidth, s.minWidth);
    s.maxWidth = Math.Max(maxWidth, s.maxWidth);
    s.avgWidth += 0.5 * (minWidth + maxWidth);
    s.minEdge = Math.Min(minEdge, s.minEdge);
    s.maxEdge = Math.Max(maxEdge, s.maxEdge);
    s.avgEdge += avgEdge;
    s.maxEdgeAspect = Math.Max(maxEdge / minEdge, s.maxEdgeAspect);
    s.minDiag = Math.Min(minDiag, s.minDiag);
    s.maxDiag = Math.Max(maxDiag, s.maxDiag);
    s.avgDiag += 0.5 * (minDiag + maxDiag);
    s.maxDiagAspect = Math.Max(maxDiag / minDiag, s.maxDiagAspect);
    s.minAngleSpan = Math.Min(minAngleSpan, s.minAngleSpan);
    s.maxAngleSpan = Math.Max(maxAngleSpan, s.maxAngleSpan);
    s.avgAngleSpan += 0.5 * (minAngleSpan + maxAngleSpan);
    double approxRatio = approxArea / exactArea;
    s.minApproxRatio = Math.Min(approxRatio, s.minApproxRatio);
    s.maxApproxRatio = Math.Max(approxRatio, s.maxApproxRatio);
  }

  public void testSubdivide(S2Cell cell) {
    gatherStats(cell);
    if (cell.isLeaf()) {
      return;
    }

    S2Cell[] children = new S2Cell[4];
    for (int i = 0; i < children.Length; ++i) {
      children[i] = new S2Cell();
    }
    Assert.True(cell.subdivide(children));
    S2CellId childId = cell.id().childBegin();
    double exactArea = 0;
    double approxArea = 0;
    double averageArea = 0;
    for (int i = 0; i < 4; ++i, childId = childId.next()) {
      exactArea += children[i].exactArea();
      approxArea += children[i].approxArea();
      averageArea += children[i].averageArea();

      // Check that the child geometry is consistent with its cell id.
      JavaAssert.Equal(children[i].id(), childId);
      Assert.True(children[i].getCenter().aequal(childId.toPoint(), 1e-15));
      S2Cell direct = new S2Cell(childId);
      JavaAssert.Equal(children[i].face(), direct.face());
      JavaAssert.Equal(children[i].level(), direct.level());
      JavaAssert.Equal(children[i].orientation(), direct.orientation());
      JavaAssert.Equal(children[i].getCenterRaw(), direct.getCenterRaw());
      for (int k = 0; k < 4; ++k) {
        JavaAssert.Equal(children[i].getVertexRaw(k), direct.getVertexRaw(k));
        JavaAssert.Equal(children[i].getEdgeRaw(k), direct.getEdgeRaw(k));
      }

      // Test Contains() and MayIntersect().
      Assert.True(cell.contains(children[i]));
      Assert.True(cell.mayIntersect(children[i]));
      Assert.True(!children[i].contains(cell));
      Assert.True(cell.contains(children[i].getCenterRaw()));
      for (int j = 0; j < 4; ++j) {
        Assert.True(cell.contains(children[i].getVertexRaw(j)));
        if (j != i) {
          Assert.True(!children[i].contains(children[j].getCenterRaw()));
          Assert.True(!children[i].mayIntersect(children[j]));
        }
      }

      // Test GetCapBound and GetRectBound.
      S2Cap parentCap = cell.getCapBound();
      S2LatLngRect parentRect = cell.getRectBound();
      if (cell.contains(new S2Point(0, 0, 1))
        || cell.contains(new S2Point(0, 0, -1))) {
        Assert.True(parentRect.Lng.isFull());
      }
      S2Cap childCap = children[i].getCapBound();
      S2LatLngRect childRect = children[i].getRectBound();
      Assert.True(childCap.contains(children[i].getCenter()));
      Assert.True(childRect.contains(children[i].getCenterRaw()));
      Assert.True(parentCap.contains(children[i].getCenter()));
      Assert.True(parentRect.contains(children[i].getCenterRaw()));
      for (int j = 0; j < 4; ++j) {
        Assert.True(childCap.contains(children[i].getVertex(j)));
        Assert.True(childRect.contains(children[i].getVertex(j)));
        Assert.True(childRect.contains(children[i].getVertexRaw(j)));
        Assert.True(parentCap.contains(children[i].getVertex(j)));
        if (!parentRect.contains(children[i].getVertex(j))) {
          Console.WriteLine("cell: " + cell + " i: " + i + " j: " + j);
          Console.WriteLine("Children " + i + ": " + children[i]);
          Console.WriteLine("Parent rect: " + parentRect);
          Console.WriteLine("Vertex raw(j) " + children[i].getVertex(j));
          Console.WriteLine("Latlng of vertex: " + new S2LatLng(children[i].getVertex(j)));
          cell.getRectBound();
        }
        Assert.True(parentRect.contains(children[i].getVertex(j)));
        if (!parentRect.contains(children[i].getVertexRaw(j))) {
          Console.WriteLine("cell: " + cell + " i: " + i + " j: " + j);
          Console.WriteLine("Children " + i + ": " + children[i]);
          Console.WriteLine("Parent rect: " + parentRect);
          Console.WriteLine("Vertex raw(j) " + children[i].getVertexRaw(j));
          Console.WriteLine("Latlng of vertex: " + new S2LatLng(children[i].getVertexRaw(j)));
          cell.getRectBound();
        }
        Assert.True(parentRect.contains(children[i].getVertexRaw(j)));
        if (j != i) {
          // The bounding caps and rectangles should be tight enough so that
          // they exclude at least two vertices of each adjacent cell.
          int capCount = 0;
          int rectCount = 0;
          for (int k = 0; k < 4; ++k) {
            if (childCap.contains(children[j].getVertex(k))) {
              ++capCount;
            }
            if (childRect.contains(children[j].getVertexRaw(k))) {
              ++rectCount;
            }
          }
          Assert.True(capCount <= 2);
          if (childRect.latLo().radians() > -S2.M_PI_2
            && childRect.latHi().radians() < S2.M_PI_2) {
            // Bounding rectangles may be too large at the poles because the
            // pole itself has an arbitrary fixed longitude.
            Assert.True(rectCount <= 2);
          }
        }
      }

      // Check all children for the first few levels, and then sample randomly.
      // Also subdivide one corner cell, one edge cell, and one center cell
      // so that we have a better chance of sample the minimum metric values.
      bool forceSubdivide = false;
      S2Point center = S2Projections.getNorm(children[i].face());
      S2Point edge = S2Point.add(center, S2Projections.getUAxis(children[i].face()));
      S2Point corner = S2Point.add(edge, S2Projections.getVAxis(children[i].face()));
      for (int j = 0; j < 4; ++j) {
        S2Point p = children[i].getVertexRaw(j);
        if (p.Equals(center) || p.Equals(edge) || p.Equals(corner))
        {
          forceSubdivide = true;
        }
      }
      if (forceSubdivide || cell.level() < (DEBUG_MODE ? 5 : 6)
        || random(DEBUG_MODE ? 10 : 4) == 0) {
        testSubdivide(children[i]);
      }
    }

    // Check sum of child areas equals parent area.
    //
    // For ExactArea(), the best relative error we can expect is about 1e-6
    // because the precision of the unit vector coordinates is only about 1e-15
    // and the edge length of a leaf cell is about 1e-9.
    //
    // For ApproxArea(), the areas are accurate to within a few percent.
    //
    // For AverageArea(), the areas themselves are not very accurate, but
    // the average area of a parent is exactly 4 times the area of a child.

    Assert.True(Math.Abs(Math.Log(exactArea / cell.exactArea())) <= Math
      .Abs(Math.Log(1 + 1e-6)));
    Assert.True(Math.Abs(Math.Log(approxArea / cell.approxArea())) <= Math
      .Abs(Math.Log(1.03)));
    Assert.True(Math.Abs(Math.Log(averageArea / cell.averageArea())) <= Math
      .Abs(Math.Log(1 + 1e-15)));
  }

  public void testMinMaxAvg(String label, int level, double count,
      double absError, double minValue, double maxValue, double avgValue,
      S2.Metric minMetric, S2.Metric maxMetric, S2.Metric avgMetric) {

    // All metrics are minimums, maximums, or averages of differential
    // quantities, and therefore will not be exact for cells at any finite
    // level. The differential minimum is always a lower bound, and the maximum
    // is always an upper bound, but these minimums and maximums may not be
    // achieved for two different reasons. First, the cells at each level are
    // sampled and we may miss the most extreme examples. Second, the actual
    // metric for a cell is obtained by integrating the differential quantity,
    // which is not constant across the cell. Therefore cells at low levels
    // (bigger cells) have smaller variations.
    //
    // The "tolerance" below is an attempt to model both of these effects.
    // At low levels, error is dominated by the variation of differential
    // quantities across the cells, while at high levels error is dominated by
    // the effects of random sampling.
    double tolerance = (maxMetric.getValue(level) - minMetric.getValue(level))
      / Math.Sqrt(Math.Min(count, 0.5 * (1L << level))) * 10;
    if (tolerance == 0) {
      tolerance = absError;
    }

    double minError = minValue - minMetric.getValue(level);
    double maxError = maxMetric.getValue(level) - maxValue;
    double avgError = Math.Abs(avgMetric.getValue(level) - avgValue);
    Console.WriteLine(
      "%-10s (%6.0f samples, tolerance %8.3g) - Min (%9.3g : %9.3g) "
        + "Max (%9.3g : %9.3g), avg (%9.3g : %9.3g)\n", label, count,
      tolerance, minError / minValue, minError / tolerance, maxError
        / maxValue, maxError / tolerance, avgError / avgValue, avgError
        / tolerance);

    Assert.True(minMetric.getValue(level) <= minValue + absError);
    Assert.True(minMetric.getValue(level) >= minValue - tolerance);
    Console.WriteLine("Level: " + maxMetric.getValue(level) + " Max " +  (maxValue + tolerance));
    Assert.True(maxMetric.getValue(level) <= maxValue + tolerance);
    Assert.True(maxMetric.getValue(level) >= maxValue - absError);
    assertDoubleNear(avgMetric.getValue(level), avgValue, 10 * tolerance);
  }

        [Test]
  public void testSubdivide() {
    for (int face = 0; face < 6; ++face) {
      testSubdivide(S2Cell.fromFacePosLevel(face, (byte) 0, 0));
    }

    // The maximum edge *ratio* is the ratio of the longest edge of any cell to
    // the shortest edge of any cell at the same level (and similarly for the
    // maximum diagonal ratio).
    //
    // The maximum edge *aspect* is the maximum ratio of the longest edge of a
    // cell to the shortest edge of that same cell (and similarly for the
    // maximum diagonal aspect).

    Console
      .WriteLine("Level    Area      Edge          Diag          Approx       Average\n");
    Console
      .WriteLine("        Ratio  Ratio Aspect  Ratio Aspect    Min    Max    Min    Max\n");
    for (int i = 0; i <= S2CellId.MAX_LEVEL; ++i) {
      LevelStats s = levelStats[i];
      if (s.count > 0) {
        s.avgArea /= s.count;
        s.avgWidth /= s.count;
        s.avgEdge /= s.count;
        s.avgDiag /= s.count;
        s.avgAngleSpan /= s.count;
      }
      Console.WriteLine(
        "%5d  %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n", i,
        s.maxArea / s.minArea, s.maxEdge / s.minEdge, s.maxEdgeAspect,
        s.maxDiag / s.minDiag, s.maxDiagAspect, s.minApproxRatio,
        s.maxApproxRatio, S2Cell.averageArea(i) / s.maxArea, S2Cell
          .averageArea(i)
          / s.minArea);
    }

    // Now check the validity of the S2 length and area metrics.
    for (int i = 0; i <= S2CellId.MAX_LEVEL; ++i) {
      LevelStats s = levelStats[i];
      if (s.count == 0) {
        continue;
      }

      Console.WriteLine(
        "Level {0} - metric (error/actual : error/tolerance)\n", i);

      // The various length calculations are only accurate to 1e-15 or so,
      // so we need to allow for this amount of discrepancy with the theoretical
      // minimums and maximums. The area calculation is accurate to about 1e-15
      // times the cell width.
      testMinMaxAvg("area", i, s.count, 1e-15 * s.minWidth, s.minArea,
        s.maxArea, s.avgArea, S2Projections.MIN_AREA, S2Projections.MAX_AREA,
        S2Projections.AVG_AREA);
      testMinMaxAvg("width", i, s.count, 1e-15, s.minWidth, s.maxWidth,
        s.avgWidth, S2Projections.MIN_WIDTH, S2Projections.MAX_WIDTH,
        S2Projections.AVG_WIDTH);
      testMinMaxAvg("edge", i, s.count, 1e-15, s.minEdge, s.maxEdge,
        s.avgEdge, S2Projections.MIN_EDGE, S2Projections.MAX_EDGE,
        S2Projections.AVG_EDGE);
      testMinMaxAvg("diagonal", i, s.count, 1e-15, s.minDiag, s.maxDiag,
        s.avgDiag, S2Projections.MIN_DIAG, S2Projections.MAX_DIAG,
        S2Projections.AVG_DIAG);
      testMinMaxAvg("angle span", i, s.count, 1e-15, s.minAngleSpan,
        s.maxAngleSpan, s.avgAngleSpan, S2Projections.MIN_ANGLE_SPAN,
        S2Projections.MAX_ANGLE_SPAN, S2Projections.AVG_ANGLE_SPAN);

      // The aspect ratio calculations are ratios of lengths and are therefore
      // less accurate at higher subdivision levels.
      Assert.True(s.maxEdgeAspect <= S2Projections.MAX_EDGE_ASPECT + 1e-15
        * (1 << i));
      Assert.True(s.maxDiagAspect <= S2Projections.MAX_DIAG_ASPECT + 1e-15
        * (1 << i));
    }
  }

  const int MAX_LEVEL = DEBUG_MODE ? 6 : 10;

  public void expandChildren1(S2Cell cell) {
    S2Cell[] children = new S2Cell[4];
    Assert.True(cell.subdivide(children));
    if (children[0].level() < MAX_LEVEL) {
      for (int pos = 0; pos < 4; ++pos) {
        expandChildren1(children[pos]);
      }
    }
  }

  public void expandChildren2(S2Cell cell) {
    S2CellId id = cell.id().childBegin();
    for (int pos = 0; pos < 4; ++pos, id = id.next()) {
      S2Cell child = new S2Cell(id);
      if (child.level() < MAX_LEVEL) {
        expandChildren2(child);
      }
    }
  }
    }
}
