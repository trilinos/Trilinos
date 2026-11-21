/*
 * Akri_Unit_InterfaceGeometry.cpp
 *
 *  Created on: May 4, 2023
 *      Author: drnoble
 */
#include <algorithm>
#include <memory>

#include <gtest/gtest.h>

#include <Akri_AnalyticSurf.hpp>
#include <Akri_AnalyticSurfaceInterfaceGeometry.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_Edge.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_Intersection_Points.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_SecondOrderMeshSpecs.hpp>
#include <Akri_Sign.hpp>
#include <Akri_StkMeshFixture.hpp>
#include <Akri_Transformation.hpp>
#include <Akri_TransformedSurface.hpp>
#include <Akri_UnitTestUtils.hpp>
#include <Akri_Unit_InterfaceGeometry.hpp>

namespace krino {

// One right tet, linear or quadratic (MESHSPEC is RightTet or RightTet10).
template <typename MESHSPEC>
class RightTetIntersectionFixture : public StkMeshFixture<MESHSPEC::TOPOLOGY>
{
public:
  RightTetIntersectionFixture()
  {
    this->set_valid_proc_sizes_for_test({1});
    this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1});
  }
  stk::mesh::Entity get_element()
  {
    const std::vector<stk::mesh::Entity> ownedElements = this->get_owned_elements();
    return ownedElements[0];
  }
  CDFEM_Support & get_cdfem_support() { return CDFEM_Support::get(this->mMesh.mesh_meta_data()); }
  Phase_Support & get_phase_support() { return Phase_Support::get(this->mMesh.mesh_meta_data()); }
  Edge get_edge(unsigned edgeOrdinal) { std::vector<Edge> elemEdges; fill_entity_edges(this->mMesh, get_element(), elemEdges); return elemEdges[edgeOrdinal]; }

  std::vector<IntersectionPoint> find_intersections_with_cuboid(const stk::math::Vector3d & center, const stk::math::Vector3d & dimensions)
  {
    Cuboid cuboid{center, dimensions};
    return find_intersections_with_surface(cuboid, 0.);
  }

  std::vector<IntersectionPoint> find_intersections_with_surface(const Surface & surface, const double time)
  {
    AnalyticSurfaceInterfaceGeometry geomIn(this->get_aux_meta().active_part(), get_cdfem_support(), get_phase_support());
    geomIn.add_surface(Surface_Identifier(0),  surface,  this->get_aux_meta().active_part());
    geomIn.set_time(time);

    const NodeToCapturedDomainsMap nodesToCapturedDomains;
    return build_all_intersection_points(this->mMesh, this->mMesh.mesh_meta_data().universal_part(), geomIn, nodesToCapturedDomains);
  }

  // Every point is in the stencil of the tet's vertices, with weights that interpolate them, and is
  // on the surface.
  void expect_points_on_vertices_and_surface(const std::vector<IntersectionPoint> & intPts, const Surface & surface)
  {
    const stk::mesh::Entity * elemNodes = this->mMesh.begin_nodes(get_element());
    const std::vector<stk::mesh::Entity> vertices(elemNodes, elemNodes + 4);
    const FieldRef coordsField(this->mMesh.mesh_meta_data().coordinate_field());
    for (auto & intPt : intPts)
    {
      double sumWeights = 0.;
      stk::math::Vector3d location = stk::math::Vector3d::ZERO;
      for (size_t i = 0; i < intPt.get_nodes().size(); ++i)
      {
        EXPECT_NE(vertices.end(), std::find(vertices.begin(), vertices.end(), intPt.get_nodes()[i]));
        EXPECT_LE(0., intPt.get_weights()[i]);
        EXPECT_GE(1., intPt.get_weights()[i]);
        sumWeights += intPt.get_weights()[i];
        location += intPt.get_weights()[i] * get_vector_field(this->mMesh, coordsField, intPt.get_nodes()[i], 3);
      }
      EXPECT_NEAR(1., sumWeights, 1.e-12);
      EXPECT_NEAR(0., surface.point_signed_distance(location), 1.e-12);
    }
  }

  void expect_num_intersections(const std::vector<IntersectionPoint> & intPts, const size_t goldNumEdgeIntPts, const size_t goldNumFaceIntPts, const size_t goldNumTetIntPts )
  {
    size_t numEdgeIntPts = 0;
    size_t numFaceIntPts = 0;
    size_t numTetIntPts = 0;

    for (auto & intPt : intPts)
    {
      const unsigned numIntPtNodes = intPt.get_nodes().size();
      if (2 == numIntPtNodes) ++numEdgeIntPts;
      else if (3 == numIntPtNodes) ++numFaceIntPts;
      else if (4 == numIntPtNodes) ++numTetIntPts;
    }
    EXPECT_EQ(goldNumEdgeIntPts, numEdgeIntPts);
    EXPECT_EQ(goldNumFaceIntPts, numFaceIntPts);
    EXPECT_EQ(goldNumTetIntPts, numTetIntPts);
  }

protected:
  MESHSPEC meshSpec;
  std::unique_ptr<AnalyticSurfaceInterfaceGeometry> geom;
};

using TetIntersectionFixture = RightTetIntersectionFixture<RightTet>;
using Tet10IntersectionFixture = RightTetIntersectionFixture<RightTet10>;

TEST_F(TetIntersectionFixture, cuboidThatIntersects3Edges)
{
  if(is_valid_proc_size_for_test())
  {
    stk::math::Vector3d dimensions(1.,2.,2.);
    stk::math::Vector3d center(1.,0.,0.);

    const std::vector<IntersectionPoint> intPts = find_intersections_with_cuboid(center, dimensions);
    expect_num_intersections(intPts, 3, 0, 0);
  }
}

// A moving cuboid is where its motion puts it at the time of the geometry: away from the tet at time
// 0, and where cuboidThatIntersects3Edges has it at time 1.
TEST_F(TetIntersectionFixture, movingCuboid_intersectionsAtTheTimeOfTheGeometry)
{
  if(is_valid_proc_size_for_test())
  {
    const stk::math::Vector3d velocity(5.,0.,0.);
    Transformation motion;
    motion.set_translational_velocity(velocity);
    TransformedSurface movingCuboid(std::make_shared<Cuboid>(stk::math::Vector3d(1.,0.,0.)-velocity, stk::math::Vector3d(1.,2.,2.)), motion);

    expect_num_intersections(find_intersections_with_surface(movingCuboid, 0.), 0, 0, 0);
    expect_num_intersections(find_intersections_with_surface(movingCuboid, 1.), 3, 0, 0);
  }
}

TEST_F(TetIntersectionFixture, cuboidThatIntersects3EdgesAnd2Faces)
{
  if(is_valid_proc_size_for_test())
  {
    stk::math::Vector3d dimensions(1.,1.,2.);
    stk::math::Vector3d center(-0.3,-0.3,0.);

    const std::vector<IntersectionPoint> intPts = find_intersections_with_cuboid(center, dimensions);
    expect_num_intersections(intPts, 3, 2, 0);
  }
}

TEST_F(TetIntersectionFixture, cuboidThatIntersects3EdgesAnd2FacesAndVol)
{
  if(is_valid_proc_size_for_test())
  {
    stk::math::Vector3d dimensions(2.,2.,2.);
    stk::math::Vector3d center(1.1,-0.9,-0.9);

    const std::vector<IntersectionPoint> intPts = find_intersections_with_cuboid(center, dimensions);
    expect_num_intersections(intPts, 3, 3, 1);
  }
}

// The filter applies to the points on faces and within the tet as to those on edges: the domains of
// a point repeat the surface for each of its sharp features that meet there, so filtering by their
// number drops the corner inside the tet (3 domains), then also the crossings of faces (2).
TEST_F(TetIntersectionFixture, cuboidThatIntersects3EdgesAnd2FacesAndVol_filteredPointsAreLeftOut)
{
  if (is_valid_proc_size_for_test())
  {
    const Cuboid cuboid(stk::math::Vector3d(1.1, -0.9, -0.9), stk::math::Vector3d(2., 2., 2.));
    AnalyticSurfaceInterfaceGeometry geometry(
        get_aux_meta().active_part(), get_cdfem_support(), get_phase_support());
    geometry.add_surface(Surface_Identifier(0), cuboid, get_aux_meta().active_part());

    auto filtered_points = [&](const size_t maxNumDomains)
    {
      const IntersectionPointFilter filter = [maxNumDomains](
                                                 const std::vector<stk::mesh::Entity> & /*nodes*/,
                                                 const std::vector<int> & sortedDomains)
      { return sortedDomains.size() <= maxNumDomains; };
      std::vector<IntersectionPoint> intPts;
      geometry.append_element_intersection_points(
          mMesh, NodeToCapturedDomainsMap{}, {get_element()}, filter, intPts);
      return intPts;
    };

    expect_num_intersections(filtered_points(3), 3, 3, 1);
    expect_num_intersections(filtered_points(2), 3, 3, 0);
    expect_num_intersections(filtered_points(1), 3, 0, 0);
  }
}

// On a tet10, the points on faces and within the tet are in the stencils of its vertices, as on a
// tet4: only vertices are snapped.
TEST_F(Tet10IntersectionFixture, cuboidThatIntersects3EdgesAnd2Faces_pointsOnTheVertices)
{
  if(is_valid_proc_size_for_test())
  {
    const Cuboid cuboid(stk::math::Vector3d(-0.3,-0.3,0.), stk::math::Vector3d(1.,1.,2.));

    const std::vector<IntersectionPoint> intPts = find_intersections_with_surface(cuboid, 0.);
    expect_num_intersections(intPts, 3, 2, 0);
    expect_points_on_vertices_and_surface(intPts, cuboid);
  }
}

TEST_F(Tet10IntersectionFixture, cuboidThatIntersects3EdgesAnd2FacesAndVol_pointsOnTheVertices)
{
  if(is_valid_proc_size_for_test())
  {
    const Cuboid cuboid(stk::math::Vector3d(1.1,-0.9,-0.9), stk::math::Vector3d(2.,2.,2.));

    const std::vector<IntersectionPoint> intPts = find_intersections_with_surface(cuboid, 0.);
    expect_num_intersections(intPts, 3, 3, 1);
    expect_points_on_vertices_and_surface(intPts, cuboid);
  }
}

// A cuboid clipping the tet's edge from (1,0,0) to (0,1,0), with every node outside it: the edge's
// two crossings are not edge intersection points, as the signs at its ends do not show them, but
// are among all of them.
TEST_F(
    TetIntersectionFixture, cuboidClippingEdge_multipleCrossingsOnlyAmongAllEdgeIntersectionPoints)
{
  if (is_valid_proc_size_for_test())
  {
    const Cuboid cuboid(stk::math::Vector3d(0.5, 0.5, 0.), stk::math::Vector3d(0.4, 0.4, 0.2));
    AnalyticSurfaceInterfaceGeometry geometry(
        get_aux_meta().active_part(), get_cdfem_support(), get_phase_support());
    geometry.add_surface(Surface_Identifier(0), cuboid, get_aux_meta().active_part());

    EXPECT_TRUE(geometry.get_edge_intersection_points(mMesh, NodeToCapturedDomainsMap{}).empty());

    const std::vector<IntersectionPoint> allIntPts =
        geometry.get_all_edge_intersection_points(mMesh);
    ASSERT_EQ(2u, allIntPts.size());
    const FieldRef coordsField(mMesh.mesh_meta_data().coordinate_field());
    std::vector<stk::math::Vector3d> locations;
    for (auto & intPt : allIntPts)
    {
      ASSERT_EQ(2u, intPt.get_nodes().size());
      locations.push_back(
          intPt.get_weights()[0] * get_vector_field(mMesh, coordsField, intPt.get_nodes()[0], 3) +
          intPt.get_weights()[1] * get_vector_field(mMesh, coordsField, intPt.get_nodes()[1], 3));
    }
    std::sort(locations.begin(),
        locations.end(),
        [](const stk::math::Vector3d & a, const stk::math::Vector3d & b) { return a[0] < b[0]; });
    expect_eq_absolute(stk::math::Vector3d(0.3, 0.7, 0.), locations[0], 1.e-12);
    expect_eq_absolute(stk::math::Vector3d(0.7, 0.3, 0.), locations[1], 1.e-12);
  }
}

static std::vector<double> multiple_crossings(
    const Surface & surface, const stk::math::Vector3d & x0, const stk::math::Vector3d & x1)
{
  std::vector<double> positions;
  fill_multiple_segment_crossings(surface, x0, x1, 1.e-10, positions);
  return positions;
}

static void
expect_positions(const std::vector<double> & expected, const std::vector<double> & positions)
{
  ASSERT_EQ(expected.size(), positions.size());
  for (size_t i = 0; i < expected.size(); ++i)
  {
    EXPECT_NEAR(expected[i], positions[i], 1.e-12);
  }
}

// The multiple crossings of an edge are those that the signs of its ends do not show, the same
// whichever end it starts from.
TEST(MultipleSegmentCrossings, cuboid_onlyCrossingsNotShownByEndSigns)
{
  const Cuboid cube(stk::math::Vector3d::ZERO, stk::math::Vector3d(1., 1., 1.));

  // From a node on a face, inside and out through the adjacent face.
  expect_positions({0.5},
      multiple_crossings(
          cube, stk::math::Vector3d(0.5, 0.3, 0.), stk::math::Vector3d(0.3, 0.7, 0.)));
  expect_positions({0.5},
      multiple_crossings(
          cube, stk::math::Vector3d(0.3, 0.7, 0.), stk::math::Vector3d(0.5, 0.3, 0.)));

  // Both ends outside, clipping a sharp edge.
  expect_positions({1. / 3., 2. / 3.},
      multiple_crossings(
          cube, stk::math::Vector3d(0.6, 0.3, 0.), stk::math::Vector3d(0.3, 0.6, 0.)));
  expect_positions({1. / 3., 2. / 3.},
      multiple_crossings(
          cube, stk::math::Vector3d(0.3, 0.6, 0.), stk::math::Vector3d(0.6, 0.3, 0.)));

  // A single crossing, which the end signs show.
  expect_positions({},
      multiple_crossings(cube, stk::math::Vector3d(0., 0., 0.), stk::math::Vector3d(1., 0., 0.)));
}

// An edge in the plane of a face, with its ends on the face to roundoff, runs along the surface
// rather than crossing it.
TEST(MultipleSegmentCrossings, cuboidEdgeInPlaneOfFace_none)
{
  const Cuboid cube(stk::math::Vector3d::ZERO, stk::math::Vector3d(1., 1., 1.));
  const stk::math::Vector3d x0(0.5 + 1.e-14, -0.2, 0.1), x1(0.5 - 1.e-14, 0.3, -0.1);
  std::vector<double> crossings;
  cube.fill_segment_crossings(x0, x1, 1.e-10, crossings);
  EXPECT_EQ(1u, crossings.size());

  expect_positions({}, multiple_crossings(cube, x0, x1));
  expect_positions({}, multiple_crossings(cube, x1, x0));
}

} // namespace krino
