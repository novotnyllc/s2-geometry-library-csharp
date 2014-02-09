using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Google.Common.Geometry;
using NUnit.Framework;

namespace S2Geometry.Tests
{
    public class S2CapTest : GeometryTestCase
    {
        public S2Point getLatLngPoint(double latDegrees, double lngDegrees)
        {
            return S2LatLng.fromDegrees(latDegrees, lngDegrees).toPoint();
        }

        // About 9 times the double-precision roundoff relative error.
        public const double EPS = 1e-15;

        public void testRectBound()
        {
            // Empty and full caps.
            Assert.True(S2Cap.empty().RectBound.isEmpty());
            Assert.True(S2Cap.full().RectBound.isFull());

            var kDegreeEps = 1e-13;
            // Maximum allowable error for latitudes and longitudes measured in
            // degrees. (assertDoubleNear uses a fixed tolerance that is too small.)

            // Cap that includes the south pole.
            var rect =
                S2Cap.fromAxisAngle(getLatLngPoint(-45, 57), S1Angle.FromDegrees(50)).RectBound;
            assertDoubleNear(rect.latLo().Degrees, -90, kDegreeEps);
            assertDoubleNear(rect.latHi().Degrees, 5, kDegreeEps);
            Assert.True(rect.Lng.isFull());

            // Cap that is tangent to the north pole.
            rect = S2Cap.fromAxisAngle(S2Point.normalize(new S2Point(1, 0, 1)), S1Angle.FromRadians(S2.M_PI_4)).RectBound;
            assertDoubleNear(rect.Lat.Lo, 0);
            assertDoubleNear(rect.Lat.Hi, S2.M_PI_2);
            Assert.True(rect.Lng.isFull());

            rect = S2Cap
                .fromAxisAngle(S2Point.normalize(new S2Point(1, 0, 1)), S1Angle.FromDegrees(45)).RectBound;
            assertDoubleNear(rect.latLo().Degrees, 0, kDegreeEps);
            assertDoubleNear(rect.latHi().Degrees, 90, kDegreeEps);
            Assert.True(rect.Lng.isFull());

            // The eastern hemisphere.
            rect = S2Cap
                .fromAxisAngle(new S2Point(0, 1, 0), S1Angle.FromRadians(S2.M_PI_2 + 5e-16)).RectBound;
            assertDoubleNear(rect.latLo().Degrees, -90, kDegreeEps);
            assertDoubleNear(rect.latHi().Degrees, 90, kDegreeEps);
            Assert.True(rect.Lng.isFull());

            // A cap centered on the equator.
            rect = S2Cap.fromAxisAngle(getLatLngPoint(0, 50), S1Angle.FromDegrees(20)).RectBound;
            assertDoubleNear(rect.latLo().Degrees, -20, kDegreeEps);
            assertDoubleNear(rect.latHi().Degrees, 20, kDegreeEps);
            assertDoubleNear(rect.lngLo().Degrees, 30, kDegreeEps);
            assertDoubleNear(rect.lngHi().Degrees, 70, kDegreeEps);

            // A cap centered on the north pole.
            rect = S2Cap.fromAxisAngle(getLatLngPoint(90, 123), S1Angle.FromDegrees(10)).RectBound;
            assertDoubleNear(rect.latLo().Degrees, 80, kDegreeEps);
            assertDoubleNear(rect.latHi().Degrees, 90, kDegreeEps);
            Assert.True(rect.Lng.isFull());
        }

        public void testCells()
        {
            // For each cube face, we construct some cells on
            // that face and some caps whose positions are relative to that face,
            // and then check for the expected intersection/containment results.

            // The distance from the center of a face to one of its vertices.
            var kFaceRadius = Math.Atan(S2.M_SQRT2);

            for (var face = 0; face < 6; ++face)
            {
                // The cell consisting of the entire face.
                var rootCell = S2Cell.fromFacePosLevel(face, (byte)0, 0);

                // A leaf cell at the midpoint of the v=1 edge.
                var edgeCell = new S2Cell(S2Projections.faceUvToXyz(face, 0, 1 - EPS));

                // A leaf cell at the u=1, v=1 corner.
                var cornerCell = new S2Cell(S2Projections.faceUvToXyz(face, 1 - EPS, 1 - EPS));

                // Quick check for full and empty caps.
                Assert.True(S2Cap.full().Contains(rootCell));
                Assert.True(!S2Cap.empty().MayIntersect(rootCell));

                // Check intersections with the bounding caps of the leaf cells that are
                // adjacent to 'corner_cell' along the Hilbert curve. Because this corner
                // is at (u=1,v=1), the curve stays locally within the same cube face.
                var first = cornerCell.id().prev().prev().prev();
                var last = cornerCell.id().next().next().next().next();
                for (var id = first; id.lessThan(last); id = id.next())
                {
                    var cell = new S2Cell(id);
                    JavaAssert.Equal(cell.CapBound.Contains(cornerCell), id.Equals(cornerCell.id()));
                    JavaAssert.Equal(
                        cell.CapBound.MayIntersect(cornerCell), id.parent().contains(cornerCell.id()));
                }

                var antiFace = (face + 3)%6; // Opposite face.
                for (var capFace = 0; capFace < 6; ++capFace)
                {
                    // A cap that barely contains all of 'cap_face'.
                    var center = S2Projections.getNorm(capFace);
                    var covering = S2Cap.fromAxisAngle(center, S1Angle.FromRadians(kFaceRadius + EPS));
                    JavaAssert.Equal(covering.Contains(rootCell), capFace == face);
                    JavaAssert.Equal(covering.MayIntersect(rootCell), capFace != antiFace);
                    JavaAssert.Equal(covering.Contains(edgeCell), center.dotProd(edgeCell.getCenter()) > 0.1);
                    JavaAssert.Equal(covering.Contains(edgeCell), covering.MayIntersect(edgeCell));
                    JavaAssert.Equal(covering.Contains(cornerCell), capFace == face);
                    JavaAssert.Equal(
                        covering.MayIntersect(cornerCell), center.dotProd(cornerCell.getCenter()) > 0);

                    // A cap that barely intersects the edges of 'cap_face'.
                    var bulging = S2Cap.fromAxisAngle(center, S1Angle.FromRadians(S2.M_PI_4 + EPS));
                    Assert.True(!bulging.Contains(rootCell));
                    JavaAssert.Equal(bulging.MayIntersect(rootCell), capFace != antiFace);
                    JavaAssert.Equal(bulging.Contains(edgeCell), capFace == face);
                    JavaAssert.Equal(bulging.MayIntersect(edgeCell), center.dotProd(edgeCell.getCenter()) > 0.1);
                    Assert.True(!bulging.Contains(cornerCell));
                    Assert.True(!bulging.MayIntersect(cornerCell));

                    // A singleton cap.
                    var singleton = S2Cap.fromAxisAngle(center, S1Angle.FromRadians(0));
                    JavaAssert.Equal(singleton.MayIntersect(rootCell), capFace == face);
                    Assert.True(!singleton.MayIntersect(edgeCell));
                    Assert.True(!singleton.MayIntersect(cornerCell));
                }
            }
        }

        [Test]
        public void S2CapBasicTest()
        {
            // Test basic properties of empty and full caps.
            var empty = S2Cap.empty();
            var full = S2Cap.full();
            Assert.True(empty.isValid());
            Assert.True(empty.isEmpty());
            Assert.True(empty.complement().isFull());
            Assert.True(full.isValid());
            Assert.True(full.isFull());
            Assert.True(full.complement().isEmpty());
            JavaAssert.Equal(full.height(), 2.0);
            assertDoubleNear(full.angle().Degrees, 180);

            // Containment and intersection of empty and full caps.
            Assert.True(empty.contains(empty));
            Assert.True(full.contains(empty));
            Assert.True(full.contains(full));
            Assert.True(!empty.interiorIntersects(empty));
            Assert.True(full.interiorIntersects(full));
            Assert.True(!full.interiorIntersects(empty));

            // Singleton cap containing the x-axis.
            var xaxis = S2Cap.fromAxisHeight(new S2Point(1, 0, 0), 0);
            Assert.True(xaxis.contains(new S2Point(1, 0, 0)));
            Assert.True(!xaxis.contains(new S2Point(1, 1e-20, 0)));
            JavaAssert.Equal(xaxis.angle().Radians, 0.0);

            // Singleton cap containing the y-axis.
            var yaxis = S2Cap.fromAxisAngle(new S2Point(0, 1, 0), S1Angle.FromRadians(0));
            Assert.True(!yaxis.contains(xaxis.axis()));
            JavaAssert.Equal(xaxis.height(), 0.0);

            // Check that the complement of a singleton cap is the full cap.
            var xcomp = xaxis.complement();
            Assert.True(xcomp.isValid());
            Assert.True(xcomp.isFull());
            Assert.True(xcomp.contains(xaxis.axis()));

            // Check that the complement of the complement is *not* the original.
            Assert.True(xcomp.complement().isValid());
            Assert.True(xcomp.complement().isEmpty());
            Assert.True(!xcomp.complement().contains(xaxis.axis()));

            // Check that very small caps can be represented accurately.
            // Here "kTinyRad" is small enough that unit vectors perturbed by this
            // amount along a tangent do not need to be renormalized.
            var kTinyRad = 1e-10;
            var tiny =
                S2Cap.fromAxisAngle(S2Point.normalize(new S2Point(1, 2, 3)), S1Angle.FromRadians(kTinyRad));
            var tangent = S2Point.normalize(S2Point.crossProd(tiny.axis(), new S2Point(3, 2, 1)));
            Assert.True(tiny.contains(S2Point.add(tiny.axis(), S2Point.mul(tangent, 0.99*kTinyRad))));
            Assert.True(!tiny.contains(S2Point.add(tiny.axis(), S2Point.mul(tangent, 1.01*kTinyRad))));

            // Basic tests on a hemispherical cap.
            var hemi = S2Cap.fromAxisHeight(S2Point.normalize(new S2Point(1, 0, 1)), 1);
            JavaAssert.Equal(hemi.complement().axis(), S2Point.neg(hemi.axis()));
            JavaAssert.Equal(hemi.complement().height(), 1.0);
            Assert.True(hemi.contains(new S2Point(1, 0, 0)));
            Assert.True(!hemi.complement().contains(new S2Point(1, 0, 0)));
            Assert.True(hemi.contains(S2Point.normalize(new S2Point(1, 0, -(1 - EPS)))));
            Assert.True(!hemi.interiorContains(S2Point.normalize(new S2Point(1, 0, -(1 + EPS)))));

            // A concave cap.
            var concave = S2Cap.fromAxisAngle(getLatLngPoint(80, 10), S1Angle.FromDegrees(150));
            Assert.True(concave.contains(getLatLngPoint(-70*(1 - EPS), 10)));
            Assert.True(!concave.contains(getLatLngPoint(-70*(1 + EPS), 10)));
            Assert.True(concave.contains(getLatLngPoint(-50*(1 - EPS), -170)));
            Assert.True(!concave.contains(getLatLngPoint(-50*(1 + EPS), -170)));

            // Cap containment tests.
            Assert.True(!empty.contains(xaxis));
            Assert.True(!empty.interiorIntersects(xaxis));
            Assert.True(full.contains(xaxis));
            Assert.True(full.interiorIntersects(xaxis));
            Assert.True(!xaxis.contains(full));
            Assert.True(!xaxis.interiorIntersects(full));
            Assert.True(xaxis.contains(xaxis));
            Assert.True(!xaxis.interiorIntersects(xaxis));
            Assert.True(xaxis.contains(empty));
            Assert.True(!xaxis.interiorIntersects(empty));
            Assert.True(hemi.contains(tiny));
            Assert.True(hemi.contains(
                S2Cap.fromAxisAngle(new S2Point(1, 0, 0), S1Angle.FromRadians(S2.M_PI_4 - EPS))));
            Assert.True(!hemi.contains(
                S2Cap.fromAxisAngle(new S2Point(1, 0, 0), S1Angle.FromRadians(S2.M_PI_4 + EPS))));
            Assert.True(concave.contains(hemi));
            Assert.True(concave.interiorIntersects(hemi.complement()));
            Assert.True(!concave.contains(S2Cap.fromAxisHeight(S2Point.neg(concave.axis()), 0.1)));
        }
    }
}