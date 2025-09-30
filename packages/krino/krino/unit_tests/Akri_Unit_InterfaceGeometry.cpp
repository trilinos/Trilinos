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

static std::vector<Edge> get_element_edges(const stk::mesh::BulkData & mesh, const std::vector<stk::mesh::Entity> & elementsToIntersect)
{
  std::vector<Edge> elemEdges;
  std::vector<Edge> uniqueEdges;
    for (auto && elem : elementsToIntersect)
  {
    fill_entity_edges(mesh, elem, elemEdges);
    uniqueEdges.insert(uniqueEdges.end(), elemEdges.begin(), elemEdges.end());
  }
  stk::util::sort_and_unique(uniqueEdges);
  return uniqueEdges;
}

void IntersectionPointFromNodalLevelsetInterfaceGeometry::append_element_intersection_points(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & /*nodesToCapturedDomains*/,
    const std::vector<stk::mesh::Entity> & elementsToIntersect,
    const IntersectionPointFilter & intersectionPointFilter,
    std::vector<IntersectionPoint> & intersectionPoints) const
{
  std::vector<Edge> edges = get_element_edges(mesh, elementsToIntersect);
  const bool intersectionPointIsOwned = true;
  for (auto edge : edges)
  {
    const std::array<stk::mesh::Entity,2> & edgeNodes = get_edge_nodes(edge);

    const double ls0 = nodeLSValues.at(edgeNodes[0]);
    const double ls1 = nodeLSValues.at(edgeNodes[1]);
    if (sign_change(ls0, ls1))
    {
      const std::vector<stk::mesh::Entity> intersectionPointNodes{edgeNodes[0], edgeNodes[1]};
      const std::vector<int> intersectionPointSortedDomains{0};
      const double loc = ls0 / (ls0-ls1);
      if (intersectionPointFilter(intersectionPointNodes, intersectionPointSortedDomains))
        intersectionPoints.emplace_back(intersectionPointIsOwned, intersectionPointNodes, std::vector<double>{1.-loc, loc}, intersectionPointSortedDomains);
    }
  }
}

void IntersectionPointFromNodalLevelsetInterfaceGeometry::set_nodal_levelset(const stk::mesh::BulkData & mesh, const std::vector<stk::mesh::EntityId> & nodeIds, const std::vector<double> & nodeLs)
{
  STK_ThrowRequire(nodeIds.size() == nodeLs.size());
  for (size_t n=0; n<nodeIds.size(); ++n)
  {
    stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodeIds[n]);
    if (mesh.is_valid(node))
      nodeLSValues[node] = nodeLs[n];
  }
}

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


