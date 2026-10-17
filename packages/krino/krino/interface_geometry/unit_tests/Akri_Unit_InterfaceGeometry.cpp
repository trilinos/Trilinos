/*
 * Akri_Unit_InterfaceGeometry.cpp
 *
 *  Created on: May 4, 2023
 *      Author: drnoble
 */
#include <Akri_AnalyticSurf.hpp>
#include <Akri_AnalyticSurfaceInterfaceGeometry.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_Unit_InterfaceGeometry.hpp>

#include <Akri_Edge.hpp>
#include <Akri_Intersection_Points.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_Sign.hpp>
#include <Akri_StkMeshFixture.hpp>
#include <gtest/gtest.h>

namespace krino {

class TetIntersectionFixture : public StkMeshTetFixture
{
public:
  TetIntersectionFixture()
  {
    set_valid_proc_sizes_for_test({1});
    StkMeshTetFixture::build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1});
  }
  stk::mesh::Entity get_element()
  {
    const std::vector<stk::mesh::Entity> ownedElements = get_owned_elements();
    return ownedElements[0];
  }
  CDFEM_Support & get_cdfem_support() { return CDFEM_Support::get(mMesh.mesh_meta_data()); }
  Phase_Support & get_phase_support() { return Phase_Support::get(mMesh.mesh_meta_data()); }
  Edge get_edge(unsigned edgeOrdinal) { std::vector<Edge> elemEdges; fill_entity_edges(mMesh, get_element(), elemEdges); return elemEdges[edgeOrdinal]; }

  std::vector<IntersectionPoint> find_intersections_with_cuboid(const stk::math::Vector3d & center, const stk::math::Vector3d & dimensions)
  {
    Cuboid cuboid{center, dimensions};
    AnalyticSurfaceInterfaceGeometry geomIn(get_aux_meta().active_part(), get_cdfem_support(), get_phase_support());
    geomIn.add_surface(Surface_Identifier(0),  cuboid,  get_aux_meta().active_part());

    const NodeToCapturedDomainsMap nodesToCapturedDomains;
    return build_all_intersection_points(mMesh, mMesh.mesh_meta_data().universal_part(), geomIn, nodesToCapturedDomains);
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
  RightTet meshSpec;
  std::unique_ptr<AnalyticSurfaceInterfaceGeometry> geom;
};

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



}


