using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Google.Common.Geometry;
using NUnit.Framework;

namespace S2Geometry.Tests
{
    public class S2PolygonBuilderTest : GeometryTestCase
    {
        // A chain represents either a polyline or a loop, depending
        // on whether "closed" is true.
        private class Chain
        {
            public readonly bool closed;
            public readonly String str;

            public Chain(String str, bool closed)
            {
                this.str = str;
                this.closed = closed;
            }
        }

        private class TestCase
        {
            // +1 = undirected, -1 = directed, 0 = either one

            // Each test case consists of a set of input loops and polylines.
            public readonly Chain[] chainsIn;

            // The expected set of output loops, directed appropriately.
            public readonly String[] loopsOut;
            public readonly double maxMerge;
            public readonly double minMerge;

            // The expected number of unused edges.
            public readonly int numUnusedEdges;
            public readonly int undirectedEdges;

            // +1 = XOR, -1 = don't XOR, 0 = either one
            public readonly int xorEdges;

            public TestCase(int undirectedEdges,
                            int xorEdges,
                            double minMerge,
                            double maxMerge,
                            Chain[] chainsIn,
                            String[] loopsOut,
                            int numUnusedEdges)
            {
                this.undirectedEdges = undirectedEdges;
                this.xorEdges = xorEdges;
                this.minMerge = minMerge;
                this.maxMerge = maxMerge;
                this.chainsIn = chainsIn;
                this.loopsOut = loopsOut;
                this.numUnusedEdges = numUnusedEdges;
            }
        }

        private readonly TestCase[] testCases = new TestCase[]
        {
            // 0: No loops.
            new TestCase(0, 0, 0.0, 10.0, new Chain[] {new Chain(null, false)}, new String[] {}, 0),

            // 1: One loop with some extra edges.
            new TestCase(0,
                         0,
                         0.0,
                         4.0,
                         new Chain[]
                         {
                             new Chain("0:0, 0:10, 10:5", true), new Chain("0:0, 5:5", false),
                             new Chain("10:5, 20:7, 30:10, 40:15, 50:3, 60:-20", false)
                         },
                         new String[] {"0:0, 0:10, 10:5"},
                         6),

            // 2: One loop that has an edge removed by XORing, plus lots of
            // extra edges.
            new TestCase(0, 1, 0.0, 1.0, // XOR
                         new Chain[]
                         {
                             new Chain("0:0, 0:10, 5:15, 10:10, 10:0", true),
                             new Chain("10:10, 12:12, 14:14, 16:16, 18:18", false),
                             new Chain("14:14, 14:16, 14:18, 14:20", false),
                             new Chain("14:18, 16:20, 18:22", false),
                             new Chain("18:12, 16:12, 14:12, 12:12", false),
                             new Chain("20:18, 18:16, 16:14, 14:12", false),
                             new Chain("20:14, 18:14, 16:14", false),
                             new Chain("5:15, 0:10", false)
                         },
                         new String[] {},
                         21),

            // 3: Three loops (two shells and one hole) that combine into one.
            new TestCase(0, 1, 0.0, 4.0, // XOR
                         new Chain[]
                         {
                             new Chain("0:0, 0:10, 5:10, 10:10, 10:5, 10:0", true),
                             new Chain("0:10, 0:15, 5:15, 5:10", true),
                             new Chain("10:10, 5:10, 5:5, 10:5", true),
                         },
                         new String[] {"0:0, 0:10, 0:15, 5:15, 5:10, 5:5, 10:5, 10:0"},
                         0),

            // 4: A big CCW triangle contained 3 CW triangular holes. The whole thing
            // looks like a pyramid of nine small triangles (with two extra edges).
            new TestCase(-1, 0, 0.0, 0.9, // Directed edges required for unique result.
                         new Chain[]
                         {
                             new Chain("0:0, 0:2, 0:4, 0:6, 1:5, 2:4, 3:3, 2:2, 1:1", true),
                             new Chain("0:2, 1:1, 1:3", true),
                             new Chain("0:4, 1:3, 1:5", true),
                             new Chain("1:3, 2:2, 2:4", true),
                             new Chain("0:0, 0:1", false),
                             new Chain("1:3, 5:7", false)
                         },
                         new String[]
                         {
                             "0:0, 0:2, 1:1",
                             "0:2, 0:4, 1:3",
                             "0:4, 0:6, 1:5",
                             "1:1, 1:3, 2:2",
                             "1:3, 1:5, 2:4",
                             "2:2, 2:4, 3:3"
                         },
                         2),

            // 5: A square divided into four subsquares. In this case we want
            // to extract the four loops rather than taking their union.
            // There are four extra edges as well.
            new TestCase(0, -1, 0.0, 4.0, // Don't XOR
                         new Chain[]
                         {
                             new Chain("0:0, 0:5, 5:5, 5:0", true),
                             new Chain("0:5, 0:10, 5:10, 5:5", true),
                             new Chain("5:0, 5:5, 10:5, 10:0", true),
                             new Chain("5:5, 5:10, 10:10, 10:5", true),
                             new Chain("0:10, 0:15, 0:20", false),
                             new Chain("20:0, 15:0, 10:0", false)
                         },
                         new String[]
                         {
                             "0:0, 0:5, 5:5, 5:0", "0:5, 0:10, 5:10, 5:5", "5:0, 5:5, 10:5, 10:0",
                             "5:5, 5:10, 10:10, 10:5"
                         },
                         4),

            // 6: Five nested loops that touch at a point.
            new TestCase(0,
                         0,
                         0.0,
                         0.8,
                         new Chain[]
                         {
                             new Chain("0:0, 0:10, 10:10, 10:0", true),
                             new Chain("0:0, 1:9, 9:9, 9:1", true), new Chain("0:0, 2:8, 8:8, 8:2", true),
                             new Chain("0:0, 3:7, 7:7, 7:3", true), new Chain("0:0, 4:6, 6:6, 6:4", true)
                         },
                         new String[]
                         {
                             "0:0, 0:10, 10:10, 10:0", "0:0, 1:9, 9:9, 9:1", "0:0, 2:8, 8:8, 8:2",
                             "0:0, 3:7, 7:7, 7:3", "0:0, 4:6, 6:6, 6:4"
                         },
                         0),


            // 7: Four diamonds nested within each other touching at two points.
            new TestCase(-1, 0, 0.0, 4.0, // Directed edges required for unique result.
                         new Chain[]
                         {
                             new Chain("0:-20, -10:0, 0:20, 10:0", true),
                             new Chain("0:10, -10:0, 0:-10, 10:0", true),
                             new Chain("0:-10, -5:0, 0:10, 5:0", true), new Chain("0:5, -5:0, 0:-5, 5:0", true)
                         },
                         new String[]
                         {
                             "0:-20, -10:0, 0:-10, 10:0", "0:-10, -5:0, 0:-5, 5:0",
                             "0:5, -5:0, 0:10, 5:0", "0:10, -10:0, 0:20, 10:0"
                         },
                         0),

            // 8: Seven diamonds nested within each other touching at one
            // point between each nested pair.
            new TestCase(0,
                         0,
                         0.0,
                         9.0,
                         new Chain[]
                         {
                             new Chain("0:-70, -70:0, 0:70, 70:0", true),
                             new Chain("0:-70, -60:0, 0:60, 60:0", true),
                             new Chain("0:-50, -60:0, 0:50, 50:0", true),
                             new Chain("0:-40, -40:0, 0:50, 40:0", true),
                             new Chain("0:-30, -30:0, 0:30, 40:0", true),
                             new Chain("0:-20, -20:0, 0:30, 20:0", true),
                             new Chain("0:-10, -20:0, 0:10, 10:0", true)
                         },
                         new String[]
                         {
                             "0:-70, -70:0, 0:70, 70:0",
                             "0:-70, -60:0, 0:60, 60:0",
                             "0:-50, -60:0, 0:50, 50:0",
                             "0:-40, -40:0, 0:50, 40:0",
                             "0:-30, -30:0, 0:30, 40:0",
                             "0:-20, -20:0, 0:30, 20:0",
                             "0:-10, -20:0, 0:10, 10:0"
                         },
                         0),

            // 9: A triangle and a self-intersecting bowtie.
            new TestCase(0,
                         0,
                         0.0,
                         4.0,
                         new Chain[]
                         {
                             new Chain("0:0, 0:10, 5:5", true), new Chain("0:20, 0:30, 10:20", false),
                             new Chain("10:20, 10:30, 0:20", false)
                         },
                         new String[] {"0:0, 0:10, 5:5"},
                         4),

            // 10: Two triangles that intersect each other.
            new TestCase(0,
                         0,
                         0.0,
                         2.0,
                         new Chain[] {new Chain("0:0, 0:10, 5:5", true), new Chain("2:2, 2:12, 7:7", true)},
                         new String[] {},
                         6),

            // 11: Four squares that combine to make a big square. The nominal
            // edges of the square are at +/-8.5 degrees in latitude and longitude.
            // All vertices except the center vertex are perturbed by up to 0.5
            // degrees in latitude and/or longitude. The various copies of the
            // center vertex are misaligned by more than this (i.e. they are
            // structured as a tree where adjacent vertices are separated by at
            // most 1 degree in latitude and/or longitude) so that the clustering
            // algorithm needs more than one iteration to find them all. Note that
            // the merged position of this vertex doesn't matter because it is XORed
            // away in the output.
            new TestCase(0, 1, 1.5, 5.8, // XOR, min_merge > sqrt(2), max_merge < 6.
                         new Chain[]
                         {
                             new Chain("-8:-8, -8:0", false),
                             new Chain("-8:1, -8:8", false),
                             new Chain("0:-9, -2:0", false),
                             new Chain("-1:1, 1:9", false),
                             new Chain("0:8, 2:2", false),
                             new Chain("0:-2, 1:-8", false),
                             new Chain("8:9, 9:1", false),
                             new Chain("9:0, 8:-9", false),
                             new Chain("9:-9, 0:-8", false),
                             new Chain("1:-9, -9:-9", false),
                             new Chain("8:0, 1:0", false),
                             new Chain("1:2, -8:0", false),
                             new Chain("-8:1, 1:-1", false),
                             new Chain("0:1, 8:1", false),
                             new Chain("-9:8, 1:8", false),
                             new Chain("0:9, 8:8", false)
                         },
                         new String[]
                         {
                             "8.5:8.5, 8.5:0.5, 8.5:-8.5, 0.5:-8.5, "
                             + "-8.5:-8.5, -8.5:0.5, -8.5:8.5, 0.5:8.5"
                         },
                         0)
        };

        private void getVertices(String str,
                                 S2Point x,
                                 S2Point y,
                                 S2Point z,
                                 double maxPerturbation,
                                 List<S2Point> vertices)
        {
            // Parse the vertices, perturb them if desired, and transform them into the
            // given frame.
            var line = makePolyline(str);

            for (var i = 0; i < line.numVertices(); ++i)
            {
                var p = line.vertex(i);
                // (p[0]*x + p[1]*y + p[2]*z).Normalize()
                var axis = S2Point.Normalize(((x * p.X) + (y * p.Y)) + (z * p.Z));
                var cap = S2Cap.FromAxisAngle(axis, S1Angle.FromRadians(maxPerturbation));
                vertices.Add(samplePoint(cap));
            }
        }

        private bool loopsEqual(S2Loop a, S2Loop b, double maxError)
        {
            // Return true if two loops have the same cyclic vertex sequence.

            if (a.NumVertices != b.NumVertices)
            {
                return false;
            }
            for (var offset = 0; offset < a.NumVertices; ++offset)
            {
                if (S2.ApproxEquals(a.Vertex(offset), b.Vertex(0), maxError))
                {
                    var success = true;
                    for (var i = 0; i < a.NumVertices; ++i)
                    {
                        if (!S2.ApproxEquals(a.Vertex(i + offset), b.Vertex(i), maxError))
                        {
                            success = false;
                            break;
                        }
                    }
                    if (success)
                    {
                        return true;
                    }
                    // Otherwise continue looping. There may be more than one candidate
                    // starting offset since vertices are only matched approximately.
                }
            }
            return false;
        }

        private bool findLoop(S2Loop loop, List<S2Loop> candidates, double maxError)
        {
            for (var i = 0; i < candidates.Count; ++i)
            {
                if (loopsEqual(loop, candidates[i], maxError))
                {
                    return true;
                }
            }
            return false;
        }

        private bool findMissingLoops(
            List<S2Loop> actual, List<S2Loop> expected, double maxError, String label)
        {
            // Dump any loops from "actual" that are not present in "expected".
            var found = false;
            for (var i = 0; i < actual.Count; ++i)
            {
                if (findLoop(actual[i], expected, maxError))
                {
                    continue;
                }
                Console.Error.WriteLine(label + " loop " + i + ":\n");
                var loop = actual[i];
                for (var j = 0; j < loop.NumVertices; ++j)
                {
                    var p = loop.Vertex(j);
                    Console.Error.WriteLine("   [" + p.X + ", " + p.Y + ", " + p.Z + "]\n");
                }
                found = true;
            }
            return found;
        }

        private void addChain(Chain chain,
                              S2Point x,
                              S2Point y,
                              S2Point z,
                              double maxPerturbation,
                              S2PolygonBuilder builder)
        {
            // Transform the given edge chain to the frame (x,y,z), perturb each vertex
            // up to the given distance, and add it to the builder.

            var vertices = new List<S2Point>();
            getVertices(chain.str, x, y, z, maxPerturbation, vertices);
            if (chain.closed)
            {
                vertices.Add(vertices[0]);
            }
            for (var i = 1; i < vertices.Count; ++i)
            {
                builder.AddEdge(vertices[i - 1], vertices[i]);
            }
        }

        private bool evalTristate(int state)
        {
            return (state > 0) ? true : (state < 0) ? false : (rand.NextDouble() > 0.5);
        }

        private bool testBuilder(TestCase test)
        {
            for (var iter = 0; iter < 200; ++iter)
            {
                // Initialize to the default options, which are changed below
                var options = S2PolygonBuilderOptions.DirectedXor;

                options.UndirectedEdges = evalTristate(test.undirectedEdges);
                options.XorEdges = evalTristate(test.xorEdges);

                // Each test has a minimum and a maximum merge distance. The merge
                // distance must be at least the given minimum to ensure that all expected
                // merging will take place, and it must be at most the given maximum to
                // ensure that no unexpected merging takes place.
                //
                // If the minimum and maximum values are different, we have some latitude
                // to perturb the vertices as long as the merge distance is adjusted
                // appropriately. If "p" is the maximum perturbation distance, "min" and
                // "max" are the min/max merge distances, and "m" is the actual merge
                // distance for this test, we require that
                //
                // x >= min + 2*p and x <= max - 2*p .
                //
                // This implies that p <= 0.25 * (max - min). We choose "p" so that it is
                // zero half of the time, and otherwise chosen randomly up to this limit.

                var minMerge = S1Angle.FromDegrees(test.minMerge).Radians;
                var maxMerge = S1Angle.FromDegrees(test.maxMerge).Radians;
                var r = Math.Max(0.0, 2*rand.NextDouble() - 1);
                var maxPerturbation = r*0.25*(maxMerge - minMerge);

                // Now we set the merge distance chosen randomly within the limits above
                // (min + 2*p and max - 2*p). Half of the time we set the merge distance
                // to the minimum value.

                r = Math.Max(0.0, 2*rand.NextDouble() - 1);
                options.MergeDistance = S1Angle.FromRadians(
                    minMerge + 2*maxPerturbation + r*(maxMerge - minMerge - 4*maxPerturbation));

                options.Validate = true;
                var builder = new S2PolygonBuilder(options);

                // On each iteration we randomly rotate the test case around the sphere.
                // This causes the S2PolygonBuilder to choose different first edges when
                // trying to build loops.
                var x = randomPoint();
                var y = S2Point.Normalize(S2Point.CrossProd(x, randomPoint()));
                var z = S2Point.Normalize(S2Point.CrossProd(x, y));

                foreach (var chain in test.chainsIn)
                {
                    addChain(chain, x, y, z, maxPerturbation, builder);
                }
                var loops = new List<S2Loop>();
                var unusedEdges = new List<S2Edge>();
                if (test.xorEdges < 0)
                {
                    builder.AssembleLoops(loops, unusedEdges);
                }
                else
                {
                    var polygon = new S2Polygon();
                    builder.AssemblePolygon(polygon, unusedEdges);
                    polygon.Release(loops);
                }
                var expected = new List<S2Loop>();
                foreach (var loop in test.loopsOut)
                {
                    var vertices = new List<S2Point>();
                    getVertices(loop, x, y, z, 0, vertices);
                    expected.Add(new S2Loop(vertices));
                }
                // We assume that the vertex locations in the expected output polygon
                // are separated from the corresponding vertex locations in the input
                // edges by at most half of the minimum merge distance. Essentially
                // this means that the expected output vertices should be near the
                // centroid of the various input vertices.
                var maxError = 0.5*minMerge + maxPerturbation;

                // Note single "|" below so that we print both sets of loops.
                if (findMissingLoops(loops, expected, maxError, "Actual")
                    | findMissingLoops(expected, loops, maxError, "Expected"))
                {
                    Console.Error.WriteLine(
                        "During iteration " + iter + ", undirected: " + options.UndirectedEdges + ", xor: "
                        + options.XorEdges + "\n\n");
                    return false;
                }
                if (unusedEdges.Count != test.numUnusedEdges)
                {
                    Console.Error.WriteLine("Wrong number of unused edges: " + unusedEdges.Count + " (should be "
                                            + test.numUnusedEdges + ")\n");
                    return false;
                }
            }
            return true;
        }

        [Test]
        public void testAssembleLoops()
        {
            var success = true;
            for (var i = 0; i < testCases.Length; ++i)
            {
                Console.WriteLine("Starting test case " + i);

                var caseSuccess = testBuilder(testCases[i]);

                Console.WriteLine("Test case " + i + " finished: " + ((caseSuccess) ? "SUCCESS" : "FAILED"));

                success &= caseSuccess;
            }
            assertTrue(success);
        }
    }
}