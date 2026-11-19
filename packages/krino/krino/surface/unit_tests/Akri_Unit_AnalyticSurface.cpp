#include <gtest/gtest.h>

#include <stk_math/StkVector.hpp>
#include <Akri_UnitTestUtils.hpp>
#include <Akri_AnalyticSurf.hpp>
#include <Akri_CubitJournal.hpp>
#include <Akri_UnitCubitJournalUtils.hpp>
#include <Akri_ML_Models.hpp>
#include <Akri_OutputUtils.hpp>

#include <array>
#include <cmath>
#include <vector>

namespace krino {

TEST(sphere, intersectWithEdge_locationCorrect)
{
  Sphere sphere{stk::math::Vector3d::ZERO, 0.35};
  const stk::math::Vector3d pt0 = stk::math::Vector3d::ZERO;
  const stk::math::Vector3d pt1(0.5,0.0,0.0);
  const double edgeCrossingTol = 1.e-5;
  const auto & [crossing, loc] = sphere.compute_intersection_with_segment(pt0, pt1, edgeCrossingTol);
  EXPECT_EQ(+1, crossing);
  EXPECT_NEAR(0.7, loc, edgeCrossingTol);
}

TEST(sphere, addToJournal_contourCloseToTheSphereInTheRegion)
{
  const double radius = 0.35;
  Sphere sphere{stk::math::Vector3d::ZERO, radius};
  CubitJournal journal;
  sphere.add_to_journal(journal, BoundingBox(stk::math::Vector3d(0.,0.,0.), stk::math::Vector3d(0.5,0.5,0.5)));
  const ParsedCubitJournal parsed = parse_cubit_journal(journal.str());

  ASSERT_FALSE(parsed.tris.empty());
  for (auto && x : parsed.nodes)
  {
    EXPECT_NEAR(radius, x.length(), 1.e-3);
    EXPECT_GE(x[0], 0.);
  }
}

TEST(cuboid, unrotatedCuboid_correctSignedDistanceForEveryFaceEdgeAndVertex)
{
  stk::math::Vector3d center = stk::math::Vector3d::ZERO;
  stk::math::Vector3d dimensions(1.,2.,0.5);
  Cuboid cuboid{center, dimensions};

  const double tol = 1.e-5;
  EXPECT_NEAR(-0.25, cuboid.point_signed_distance(center), tol);
  // face distances
  EXPECT_NEAR(0.2, cuboid.point_signed_distance(stk::math::Vector3d(0.7,0.,0.)), tol);
  EXPECT_NEAR(-0.1, cuboid.point_signed_distance(stk::math::Vector3d(0.4,0.,0.)), tol);
  EXPECT_NEAR(0.2, cuboid.point_signed_distance(stk::math::Vector3d(-0.7,0.,0.)), tol);
  EXPECT_NEAR(-0.1, cuboid.point_signed_distance(stk::math::Vector3d(-0.4,0.,0.)), tol);
  EXPECT_NEAR(0.2, cuboid.point_signed_distance(stk::math::Vector3d(0.,1.2,0.)), tol);
  EXPECT_NEAR(-0.1, cuboid.point_signed_distance(stk::math::Vector3d(0.,0.9,0.)), tol);
  EXPECT_NEAR(0.2, cuboid.point_signed_distance(stk::math::Vector3d(0.,0.,0.45)), tol);
  EXPECT_NEAR(-0.1, cuboid.point_signed_distance(stk::math::Vector3d(0.,0.,0.15)), tol);
  EXPECT_NEAR(0.2, cuboid.point_signed_distance(stk::math::Vector3d(0.,0.,-0.45)), tol);
  EXPECT_NEAR(-0.1, cuboid.point_signed_distance(stk::math::Vector3d(0.,0.,-0.15)), tol);
  // edge distances
  EXPECT_NEAR(std::sqrt(2.), cuboid.point_signed_distance(0.5*(cuboid.vertex_location(0)+cuboid.vertex_location(1))+stk::math::Vector3d(0.,-1.,-1.)), tol);
  EXPECT_NEAR(std::sqrt(2.), cuboid.point_signed_distance(0.5*(cuboid.vertex_location(1)+cuboid.vertex_location(2))+stk::math::Vector3d(+1.,0.,-1.)), tol);
  EXPECT_NEAR(std::sqrt(2.), cuboid.point_signed_distance(0.5*(cuboid.vertex_location(2)+cuboid.vertex_location(3))+stk::math::Vector3d(0.,+1.,-1.)), tol);
  EXPECT_NEAR(std::sqrt(2.), cuboid.point_signed_distance(0.5*(cuboid.vertex_location(3)+cuboid.vertex_location(0))+stk::math::Vector3d(-1.,0.,-1.)), tol);
  EXPECT_NEAR(std::sqrt(2.), cuboid.point_signed_distance(0.5*(cuboid.vertex_location(4)+cuboid.vertex_location(5))+stk::math::Vector3d(0.,-1.,+1.)), tol);
  EXPECT_NEAR(std::sqrt(2.), cuboid.point_signed_distance(0.5*(cuboid.vertex_location(5)+cuboid.vertex_location(6))+stk::math::Vector3d(+1.,0.,+1.)), tol);
  EXPECT_NEAR(std::sqrt(2.), cuboid.point_signed_distance(0.5*(cuboid.vertex_location(6)+cuboid.vertex_location(7))+stk::math::Vector3d(0.,+1.,+1.)), tol);
  EXPECT_NEAR(std::sqrt(2.), cuboid.point_signed_distance(0.5*(cuboid.vertex_location(7)+cuboid.vertex_location(4))+stk::math::Vector3d(-1.,0.,+1.)), tol);
  EXPECT_NEAR(std::sqrt(2.), cuboid.point_signed_distance(0.5*(cuboid.vertex_location(0)+cuboid.vertex_location(4))+stk::math::Vector3d(-1.,-1.,0.)), tol);
  EXPECT_NEAR(std::sqrt(2.), cuboid.point_signed_distance(0.5*(cuboid.vertex_location(1)+cuboid.vertex_location(5))+stk::math::Vector3d(+1.,-1.,0.)), tol);
  EXPECT_NEAR(std::sqrt(2.), cuboid.point_signed_distance(0.5*(cuboid.vertex_location(2)+cuboid.vertex_location(6))+stk::math::Vector3d(+1.,+1.,0.)), tol);
  EXPECT_NEAR(std::sqrt(2.), cuboid.point_signed_distance(0.5*(cuboid.vertex_location(3)+cuboid.vertex_location(7))+stk::math::Vector3d(-1.,+1.,0.)), tol);
  //corner distances
  EXPECT_NEAR(std::sqrt(3.), cuboid.point_signed_distance(cuboid.vertex_location(0)+stk::math::Vector3d(-1.,-1.,-1.)), tol);
  EXPECT_NEAR(std::sqrt(3.), cuboid.point_signed_distance(cuboid.vertex_location(1)+stk::math::Vector3d(+1.,-1.,-1.)), tol);
  EXPECT_NEAR(std::sqrt(3.), cuboid.point_signed_distance(cuboid.vertex_location(2)+stk::math::Vector3d(+1.,+1.,-1.)), tol);
  EXPECT_NEAR(std::sqrt(3.), cuboid.point_signed_distance(cuboid.vertex_location(3)+stk::math::Vector3d(-1.,+1.,-1.)), tol);
  EXPECT_NEAR(std::sqrt(3.), cuboid.point_signed_distance(cuboid.vertex_location(4)+stk::math::Vector3d(-1.,-1.,+1.)), tol);
  EXPECT_NEAR(std::sqrt(3.), cuboid.point_signed_distance(cuboid.vertex_location(5)+stk::math::Vector3d(+1.,-1.,+1.)), tol);
  EXPECT_NEAR(std::sqrt(3.), cuboid.point_signed_distance(cuboid.vertex_location(6)+stk::math::Vector3d(+1.,+1.,+1.)), tol);
  EXPECT_NEAR(std::sqrt(3.), cuboid.point_signed_distance(cuboid.vertex_location(7)+stk::math::Vector3d(-1.,+1.,+1.)), tol);
}

TEST(cuboid, rotatedAndTranslatedCuboid_correctSignedDistanceFromCenter)
{
  stk::math::Vector3d center(1.,-2.,3.);
  stk::math::Vector3d dimensions(1.,2.,0.5);
  stk::math::Vector3d rotation(0.5,0.25,0.1);
  Cuboid cuboid{center, dimensions, rotation};

  const double tol = 1.e-5;
  EXPECT_NEAR(-0.25, cuboid.point_signed_distance(center), tol);
}

// The vertices of a rotated cuboid off the origin, and the midpoints of its edges, are on its
// surface.
TEST(cuboid, rotatedAndTranslatedCuboid_verticesAndEdgeMidpointsOnSurface)
{
  const stk::math::Vector3d center(1.,-2.,3.);
  const stk::math::Vector3d dimensions(1.,2.,0.5);
  const Cuboid cuboid{center, dimensions, stk::math::Vector3d(0.5,0.25,0.1)};

  const double tol = 1.e-12;
  const double halfDiagonal = 0.5*dimensions.length();
  for (unsigned i=0; i<8; ++i)
  {
    EXPECT_NEAR(0., cuboid.point_signed_distance(cuboid.vertex_location(i)), tol) << "vertex " << i;
    EXPECT_NEAR(halfDiagonal, (cuboid.vertex_location(i)-center).length(), tol) << "vertex " << i;
  }
  const std::array<std::array<unsigned,2>,12> edges{{{{0,1}}, {{1,2}}, {{2,3}}, {{3,0}}, {{4,5}}, {{5,6}}, {{6,7}}, {{7,4}}, {{0,4}}, {{1,5}}, {{2,6}}, {{3,7}}}};
  for (auto & edge : edges)
  {
    const stk::math::Vector3d midpoint = 0.5*(cuboid.vertex_location(edge[0])+cuboid.vertex_location(edge[1]));
    EXPECT_NEAR(0., cuboid.point_signed_distance(midpoint), tol) << "edge " << edge[0] << "-" << edge[1];
  }
}

namespace
{
// A quarter turn about z of a cuboid off the origin swaps its x and y extents about its center.
const stk::math::Vector3d quarterTurnCenter(1.,-2.,3.);
const stk::math::Vector3d quarterTurnHalfExtents(1.,0.5,0.25);
Cuboid quarter_turned_cuboid(const int sign = 1)
{
  return Cuboid{quarterTurnCenter, stk::math::Vector3d(1.,2.,0.5), stk::math::Vector3d(0.,0.,0.5*M_PI), sign};
}
} // namespace

TEST(cuboid, quarterTurnAboutZ_distancesAndBoundingBoxOfTheTurnedCuboid)
{
  const Cuboid cuboid = quarter_turned_cuboid();
  const double tol = 1.e-12;
  EXPECT_NEAR(0.2, cuboid.point_signed_distance(quarterTurnCenter+stk::math::Vector3d(1.2,0.,0.)), tol);
  EXPECT_NEAR(-0.1, cuboid.point_signed_distance(quarterTurnCenter+stk::math::Vector3d(0.,-0.4,0.)), tol);
  EXPECT_NEAR(0.2, cuboid.point_signed_distance(quarterTurnCenter+stk::math::Vector3d(0.,0.,0.45)), tol);

  BoundingBox bbox;
  cuboid.insert_into(bbox);
  const double bboxTol = 1.e-5; // BoundingBox is single precision
  for (int i=0; i<3; ++i)
  {
    EXPECT_NEAR((quarterTurnCenter-quarterTurnHalfExtents)[i], bbox.get_min()[i], bboxTol);
    EXPECT_NEAR((quarterTurnCenter+quarterTurnHalfExtents)[i], bbox.get_max()[i], bboxTol);
  }
}

TEST(cuboid, quarterTurnAboutZ_cornerInTetAndEdgeCrossingTriangle)
{
  const Cuboid cuboid = quarter_turned_cuboid();
  const double tol = 1.e-12;
  const stk::math::Vector3d corner = quarterTurnCenter+quarterTurnHalfExtents;
  std::vector<stk::math::Vector3d> intParamCoords;

  const std::array<stk::math::Vector3d,4> tetAroundCorner{{corner+stk::math::Vector3d(-0.05,-0.05,-0.05),
      corner+stk::math::Vector3d(0.15,-0.05,-0.05), corner+stk::math::Vector3d(-0.05,0.15,-0.05), corner+stk::math::Vector3d(-0.05,-0.05,0.15)}};
  cuboid.fill_tetrahedon_intersection_parametric_coordinates(tetAroundCorner, intParamCoords);
  ASSERT_EQ(1u, intParamCoords.size());
  const stk::math::Vector3d & p = intParamCoords[0];
  const stk::math::Vector3d cornerFromTet = tetAroundCorner[0] + p[0]*(tetAroundCorner[1]-tetAroundCorner[0]) +
      p[1]*(tetAroundCorner[2]-tetAroundCorner[0]) + p[2]*(tetAroundCorner[3]-tetAroundCorner[0]);
  EXPECT_NEAR(0., (cornerFromTet-corner).length(), tol);

  // The edge along z through this corner crosses a triangle at the height of the center.
  const stk::math::Vector3d edgePoint(corner[0], corner[1], quarterTurnCenter[2]);
  const std::array<stk::math::Vector3d,3> triAroundEdge{{edgePoint+stk::math::Vector3d(-0.05,-0.05,0.),
      edgePoint+stk::math::Vector3d(0.15,-0.05,0.), edgePoint+stk::math::Vector3d(-0.05,0.15,0.)}};
  cuboid.fill_triangle_intersection_parametric_coordinates(triAroundEdge, intParamCoords);
  ASSERT_EQ(1u, intParamCoords.size());
  const stk::math::Vector3d & q = intParamCoords[0];
  const stk::math::Vector3d edgePointFromTri = triAroundEdge[0] + q[0]*(triAroundEdge[1]-triAroundEdge[0]) + q[1]*(triAroundEdge[2]-triAroundEdge[0]);
  EXPECT_NEAR(0., (edgePointFromTri-edgePoint).length(), tol);

  // Neither a corner nor an edge is near the middle of a face.
  const stk::math::Vector3d faceMiddle = quarterTurnCenter+stk::math::Vector3d(quarterTurnHalfExtents[0],0.,0.);
  const std::array<stk::math::Vector3d,4> tetAtFaceMiddle{{faceMiddle+stk::math::Vector3d(-0.05,-0.05,-0.05),
      faceMiddle+stk::math::Vector3d(0.15,-0.05,-0.05), faceMiddle+stk::math::Vector3d(-0.05,0.15,-0.05), faceMiddle+stk::math::Vector3d(-0.05,-0.05,0.15)}};
  cuboid.fill_tetrahedon_intersection_parametric_coordinates(tetAtFaceMiddle, intParamCoords);
  EXPECT_TRUE(intParamCoords.empty());
  cuboid.fill_triangle_intersection_parametric_coordinates({{tetAtFaceMiddle[0], tetAtFaceMiddle[1], tetAtFaceMiddle[2]}}, intParamCoords);
  EXPECT_TRUE(intParamCoords.empty());
}

namespace
{
void expect_segment_crossings(const Surface & surface,
    const stk::math::Vector3d & pt0,
    const stk::math::Vector3d & pt1,
    const std::vector<double> & goldPositions)
{
  std::vector<double> positions;
  surface.fill_segment_crossings(pt0, pt1, 1.e-10, positions);
  ASSERT_EQ(goldPositions.size(), positions.size()) << pt0 << " to " << pt1;
  for (size_t i = 0; i < positions.size(); ++i)
  {
    EXPECT_NEAR(goldPositions[i], positions[i], 1.e-12) << pt0 << " to " << pt1;
  }
}
} // namespace

namespace
{
// A surface whose signed distance is 0.1 x (x - 0.5) + distanceAtOrigin: on the x axis it is on the
// surface at the origin, to within distanceAtOrigin, below it up to x = 0.5, and above it beyond.
class SurfaceDippingBesideOriginForSegmentTest : public SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign
{
public:
  explicit SurfaceDippingBesideOriginForSegmentTest(const double distanceAtOrigin) : myDistanceAtOrigin(distanceAtOrigin) {}
  Surface_Type type() const override { return STRING_FUNCTION; }
  size_t storage_size() const override { return sizeof(SurfaceDippingBesideOriginForSegmentTest); }
  double point_signed_distance(const stk::math::Vector3d & x) const override
  {
    return 0.1 * x[0] * (x[0] - 0.5) + myDistanceAtOrigin;
  }

private:
  double myDistanceAtOrigin;
};
} // namespace

// A segment from a point on the surface, such as a node snapped onto it, runs below the surface and
// crosses back at x = 0.5.  The distance at that end is roundoff, of either sign, and the crossing is
// found whatever its sign, walking the segment either way.
TEST(surface, segmentCrossings_crossingBesideEndOnTheSurfaceFoundForEitherSignOfRoundoffThere)
{
  const stk::math::Vector3d origin = stk::math::Vector3d::ZERO;
  const stk::math::Vector3d end(1., 0., 0.);
  for (const double distanceAtOrigin : {-1.e-17, 0., 1.e-17})
  {
    const SurfaceDippingBesideOriginForSegmentTest surface(distanceAtOrigin);
    expect_segment_crossings(surface, origin, end, {0.5});
    expect_segment_crossings(surface, end, origin, {0.5});
  }
}

// Every crossing of the faces x = 0.5 and y = 0.5 of a long cuboid, also where the signs of the
// distances at the ends of the segment do not change.
TEST(cuboid, segmentCrossings_alsoThoseTheEndSignsDoNotShow)
{
  const Cuboid cuboid{stk::math::Vector3d::ZERO, stk::math::Vector3d(1., 1., 10.)};
  // From a node on the +x face, inside, then out through y = 0.5 (the crossing at the start is left
  // out).
  expect_segment_crossings(cuboid, {0.5, 0.3, 0.}, {0.45, 0.7, 0.}, {0.5});
  // Both ends outside, clipping the region near the edge x = y = 0.5.
  expect_segment_crossings(cuboid, {0.6, 0.3, 0.}, {0.3, 0.8, 0.}, {1. / 3., 0.4});
  // From inside to outside, through x = 0.5.
  expect_segment_crossings(cuboid, {0., 0., 0.}, {1., 0., 0.}, {0.5});
  // Missing the cuboid, and in the plane of its +x face.
  expect_segment_crossings(cuboid, {0.6, 0., 0.}, {0.8, 0., 0.}, {});
  expect_segment_crossings(cuboid, {0.5, -0.2, 0.}, {0.5, 0.2, 0.}, {});
}

TEST(cuboid, quarterTurnAboutZ_segmentCrossingsOfTheTurnedCuboid)
{
  const Cuboid cuboid = quarter_turned_cuboid();
  // Along x through the center: the turned cuboid extends 1 either side.
  expect_segment_crossings(cuboid,
      quarterTurnCenter - stk::math::Vector3d(2., 0., 0.),
      quarterTurnCenter + stk::math::Vector3d(2., 0., 0.),
      {0.25, 0.75});
}

TEST(cuboid, negativeSign_flipsDistance)
{
  const Cuboid cuboid = quarter_turned_cuboid();
  const Cuboid flipped = quarter_turned_cuboid(-1);
  for (const stk::math::Vector3d & x : {quarterTurnCenter, quarterTurnCenter+stk::math::Vector3d(1.2,0.,0.), quarterTurnCenter+quarterTurnHalfExtents+stk::math::Vector3d(0.1,0.2,0.3)})
  {
    EXPECT_DOUBLE_EQ(-cuboid.point_signed_distance(x), flipped.point_signed_distance(x));
  }
}

TEST(INR, sphere_cone)
{
  if (!have_ML_Model_support())
  {
    GTEST_SKIP() << "Skipping because ML Model Support is not enabled.";
  }

  const std::string modelBaseName = "test";
  if (!does_file_exist(modelBaseName + ".model"))
  {
    GTEST_SKIP() << "Skipping because file " << modelBaseName << ".model was not found.";
  }

  Implicit_Neural_Representation inr(modelBaseName);
  const double tol = 1.e-6;

  EXPECT_NEAR(-0.1795189, inr.point_signed_distance(stk::math::Vector3d(0., 0., 0.)), tol);
  EXPECT_NEAR(0.52966934, inr.point_signed_distance(stk::math::Vector3d(1., 0., 0.)), tol);
  EXPECT_NEAR(0.6886717, inr.point_signed_distance(stk::math::Vector3d(0., 1., 0.)), tol);
  EXPECT_NEAR(0.6959222, inr.point_signed_distance(stk::math::Vector3d(0., 0., 1.)), tol);
}

}
