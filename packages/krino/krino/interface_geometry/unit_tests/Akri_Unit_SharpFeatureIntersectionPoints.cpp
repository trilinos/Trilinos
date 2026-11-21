// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include <gtest/gtest.h>

#include <Akri_AnalyticSurf.hpp>
#include <Akri_AnalyticSurfaceInterfaceGeometry.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_StkMeshFixture.hpp>
#include <Akri_Unit_SharpFeatureIntersectionPointChecks.hpp>

namespace krino {

namespace {

constexpr double edgeTol = 1.e-3;

template <typename MESHSPEC>
class SharpFeatureIntersectionPointsFixture : public SharpFeaturePointChecks<StkMeshTetFixture>
{
public:
  SharpFeatureIntersectionPointsFixture()
  {
    set_valid_proc_sizes_for_test({1});
    build_mesh_from_spec(meshSpec);
  }

  std::unique_ptr<InterfaceGeometry> build_geometry(const Surface & surface)
  {
    auto geometry = std::make_unique<AnalyticSurfaceInterfaceGeometry>(get_aux_meta().active_part(),
        CDFEM_Support::get(mMesh.mesh_meta_data()),
        Phase_Support::get(mMesh.mesh_meta_data()));
    geometry->add_surface(Surface_Identifier(0), surface, get_aux_meta().active_part());
    return geometry;
  }

protected:
  MESHSPEC meshSpec;
};

typedef SharpFeatureIntersectionPointsFixture<TwoRightTets> TwoRightTetsSharpFeaturePoints;
typedef SharpFeatureIntersectionPointsFixture<FourRightTets> FourRightTetsSharpFeaturePoints;

constexpr int cornerCodimension = 3;
constexpr int sharpEdgeCodimension = 2;
}

// A sharp edge along x at y=z=0.3 crosses the face x=0 shared by the two tets, and an outer face of
// each: both tets have the same point on the shared face.
TEST_F(TwoRightTetsSharpFeaturePoints, sharpEdgeCrossingSharedFace_samePointInBothTets)
{
  if (is_valid_proc_size_for_test())
  {
    const Cuboid cuboid(stk::math::Vector3d(0., 0.8, 0.8), stk::math::Vector3d(4., 1., 1.));
    const auto geometry = build_geometry(cuboid);

    EXPECT_NEAR(0.4,
        expect_same_point_on_in_elements(
            *geometry, {0, 1}, {0, 2, 4}, sharpEdgeCodimension, edgeTol),
        1.e-12);
    expect_same_point_on_in_elements(*geometry, {0}, {1, 2, 4}, sharpEdgeCodimension, edgeTol);
    expect_same_point_on_in_elements(*geometry, {1}, {2, 3, 4}, sharpEdgeCodimension, edgeTol);
    for (size_t elementIndex : {0u, 1u})
    {
      EXPECT_EQ(2u, element_points(*geometry, elementIndex, sharpEdgeCodimension, edgeTol).size())
          << "element " << elementIndex;
      EXPECT_TRUE(element_points(*geometry, elementIndex, cornerCodimension, edgeTol).empty())
          << "element " << elementIndex;
    }
  }
}

// A sharp edge along y at x=0.0003, z=0.4 crosses the face y=0 shared by the first and last tets
// within the tolerance of their edge from node 0 to node 5: both have that point on the edge.  Its
// crossings of their outer faces through node 1 are as near the edges opposite node 1.  The other
// tets, at x<0, have none.
TEST_F(FourRightTetsSharpFeaturePoints, sharpEdgeCrossingFaceNearTetEdge_edgePoint)
{
  if (is_valid_proc_size_for_test())
  {
    const Cuboid cuboid(stk::math::Vector3d(1.00015, 0., 0.9), stk::math::Vector3d(1.9997, 4., 1.));
    const auto geometry = build_geometry(cuboid);

    EXPECT_NEAR(0.5997 / 0.9997,
        expect_same_point_on_in_elements(*geometry, {0, 3}, {0, 5}, sharpEdgeCodimension, edgeTol),
        1.e-12);
    expect_same_point_on_in_elements(*geometry, {0}, {2, 5}, sharpEdgeCodimension, edgeTol);
    expect_same_point_on_in_elements(*geometry, {3}, {4, 5}, sharpEdgeCodimension, edgeTol);
    for (size_t elementIndex : {0u, 3u})
    {
      EXPECT_EQ(2u, element_points(*geometry, elementIndex, sharpEdgeCodimension, edgeTol).size())
          << "element " << elementIndex;
    }
    for (size_t elementIndex : {1u, 2u})
    {
      EXPECT_TRUE(element_points(*geometry, elementIndex, sharpEdgeCodimension, edgeTol).empty())
          << "element " << elementIndex;
    }
  }
}

// A corner at (0.0002, 0.3, 0.3), inside the first tet within the tolerance of the face it shares
// with the second: a point on that face.
TEST_F(TwoRightTetsSharpFeaturePoints, cornerNearSharedFace_facePoint)
{
  if (is_valid_proc_size_for_test())
  {
    const Cuboid cuboid(stk::math::Vector3d(0.5002, 0.8, 0.8), stk::math::Vector3d(1., 1., 1.));
    const auto geometry = build_geometry(cuboid);

    const std::vector<IntersectionPoint> corners =
        element_points(*geometry, 0, cornerCodimension, edgeTol);
    EXPECT_EQ(1u, corners.size());
    EXPECT_EQ(1u, points_on(corners, {0, 2, 4}).size());
  }
}

// A corner at (0.0002, 0.0002, 0.0002), within the tolerance of node 0: a point at that node.
TEST_F(TwoRightTetsSharpFeaturePoints, cornerNearNode_nodePoint)
{
  if (is_valid_proc_size_for_test())
  {
    const Cuboid cuboid(stk::math::Vector3d(0.5002, 0.5002, 0.5002), stk::math::Vector3d(1., 1., 1.));
    const auto geometry = build_geometry(cuboid);

    const std::vector<IntersectionPoint> corners =
        element_points(*geometry, 0, cornerCodimension, edgeTol);
    EXPECT_EQ(1u, corners.size());
    EXPECT_EQ(1u, points_on(corners, {0}).size());
  }
}

typedef SharpFeatureIntersectionPointsFixture<RightTet> RightTetSharpFeaturePoints;

// Of a mix of points, only those of the rank asked for, each reduced, in canonical order (more
// nodes first), without repeats: the same point listed with its nodes in another order is removed,
// but a point within the tolerance of another is kept.
TEST_F(RightTetSharpFeaturePoints, reducedPointsOfCodimension_onlyThatRankReducedSortedUnique)
{
  if (is_valid_proc_size_for_test())
  {
    const stk::mesh::Entity n0 = get_assigned_node_for_index(0);
    const stk::mesh::Entity n1 = get_assigned_node_for_index(1);
    const stk::mesh::Entity n2 = get_assigned_node_for_index(2);
    const stk::mesh::Entity n3 = get_assigned_node_for_index(3);
    const std::vector<IntersectionPoint> points{
        IntersectionPoint(true, {n1, n2, n3}, {0.3, 0.3, 0.4}, {0}),
        IntersectionPoint(true, {n2, n1, n3}, {0.4995, 0.5, 0.0005}, {0, 0}),
        IntersectionPoint(true, {n0, n1, n2, n3}, {0.0002, 0.3, 0.3, 0.3998}, {0, 0, 0}),
        IntersectionPoint(true, {n0, n1, n2}, {0.2, 0.3, 0.5}, {0, 0, 0}),
        IntersectionPoint(true, {n0, n1, n2}, {0.2, 0.3, 0.5}, {0, 1, 2}),
        IntersectionPoint(true, {n0, n1, n2}, {0.2, 0.3, 0.5}, {0, 0}),
        IntersectionPoint(true, {n0, n2, n1}, {0.2, 0.5, 0.3}, {0, 0}),
        IntersectionPoint(true, {n0, n2, n1}, {0.2003, 0.4997, 0.3}, {0, 0})};

    const std::vector<IntersectionPoint> edgeCrossings =
        reduced_intersection_points_of_codimension(mMesh, points, sharpEdgeCodimension, edgeTol);
    ASSERT_EQ(4u, edgeCrossings.size());
    EXPECT_EQ(std::vector<int>({0, 1, 2}), edgeCrossings[0].get_sorted_domains());
    for (size_t i : {1u, 2u})
    {
      EXPECT_EQ(3u, edgeCrossings[i].get_nodes().size());
      EXPECT_EQ(std::vector<int>({0, 0}), edgeCrossings[i].get_sorted_domains());
    }
    EXPECT_NEAR(0.2, weight_of_node(edgeCrossings[1], 0), 1.e-12);
    EXPECT_NEAR(0.2003, weight_of_node(edgeCrossings[2], 0), 1.e-12);
    EXPECT_EQ(1u, points_on(edgeCrossings, {1, 2}).size());
    EXPECT_NEAR(0.5, weight_of_node(edgeCrossings[3], 1), 1.e-3);

    const std::vector<IntersectionPoint> corners =
        reduced_intersection_points_of_codimension(mMesh, points, cornerCodimension, edgeTol);
    ASSERT_EQ(2u, corners.size());
    EXPECT_EQ(1u, points_on(corners, {1, 2, 3}).size());
    EXPECT_EQ(1u, points_on(corners, {0, 1, 2}).size());
    for (auto & corner : corners)
    {
      for (size_t i = 1; i < corner.get_nodes().size(); ++i)
      {
        EXPECT_LT(
            mMesh.identifier(corner.get_nodes()[i - 1]), mMesh.identifier(corner.get_nodes()[i]));
      }
    }
  }
}

// A cuboid whose corner is eta inside the face x+y+z=1 of the right tet (the face of nodes 1, 2,
// 3): the corner's weight for node 0 is eta, and its three edges cross that face eta from it.
static Cuboid cuboid_with_corner_inside_face_of_right_tet(const double eta)
{
  const double a = (1. - eta) / 3.;
  return Cuboid(stk::math::Vector3d(a + 0.5, a + 0.5, a + 0.5), stk::math::Vector3d(1., 1., 1.));
}

// Within the tolerance of the face, the corner is on it.  Its edges' crossings of the face are
// three points, though within the tolerance of each other and of the corner: that the corner
// captures them is decided when cutting, by the domains of the corner's node.
TEST_F(RightTetSharpFeaturePoints, cornerJustInsideFaceWithinTolerance_cornerOnFace)
{
  if (is_valid_proc_size_for_test())
  {
    for (const double eta : {1.e-4, 5.e-4})
    {
      const Cuboid cuboid = cuboid_with_corner_inside_face_of_right_tet(eta);
      const auto geometry = build_geometry(cuboid);

      const std::vector<IntersectionPoint> corners =
          element_points(*geometry, 0, cornerCodimension, edgeTol);
      ASSERT_EQ(1u, corners.size()) << "eta " << eta;
      ASSERT_EQ(1u, points_on(corners, {1, 2, 3}).size()) << "eta " << eta;
      for (const double weight : corners[0].get_weights())
      {
        EXPECT_NEAR(1. / 3., weight, 1.e-12) << "eta " << eta;
      }

      const std::vector<IntersectionPoint> crossings =
          element_points(*geometry, 0, sharpEdgeCodimension, edgeTol);
      EXPECT_EQ(3u, crossings.size()) << "eta " << eta;
      EXPECT_EQ(3u, points_on(crossings, {1, 2, 3}).size()) << "eta " << eta;
    }
  }
}

// Beyond the tolerance, the corner is inside the tet, and its three edges cross the face at three
// separate points, each a weight eta nearer a different node of the face.
TEST_F(
    RightTetSharpFeaturePoints, cornerJustInsideFaceBeyondTolerance_cornerAndItsThreeEdgeCrossings)
{
  if (is_valid_proc_size_for_test())
  {
    for (const double eta : {1.5e-3, 1.e-2, 1.e-1})
    {
      const double a = (1. - eta) / 3.;
      const Cuboid cuboid = cuboid_with_corner_inside_face_of_right_tet(eta);
      const auto geometry = build_geometry(cuboid);

      const std::vector<IntersectionPoint> corners =
          element_points(*geometry, 0, cornerCodimension, edgeTol);
      ASSERT_EQ(1u, corners.size()) << "eta " << eta;
      ASSERT_EQ(1u, points_on(corners, {0, 1, 2, 3}).size()) << "eta " << eta;
      EXPECT_NEAR(eta, weight_of_node(corners[0], 0), 1.e-12) << "eta " << eta;

      const std::vector<IntersectionPoint> crossings =
          element_points(*geometry, 0, sharpEdgeCodimension, edgeTol);
      ASSERT_EQ(3u, crossings.size()) << "eta " << eta;
      ASSERT_EQ(3u, points_on(crossings, {1, 2, 3}).size()) << "eta " << eta;
      for (size_t nodeIndex : {1u, 2u, 3u})
      {
        const size_t numNearerThisNode = std::count_if(crossings.begin(),
            crossings.end(),
            [&](const IntersectionPoint & crossing)
            { return std::abs(weight_of_node(crossing, nodeIndex) - (a + eta)) < 1.e-12; });
        EXPECT_EQ(1u, numNearerThisNode) << "eta " << eta << " node " << nodeIndex;
      }
    }
  }
}

} // namespace krino
