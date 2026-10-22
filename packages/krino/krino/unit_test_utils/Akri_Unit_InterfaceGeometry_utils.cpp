/*
 * Akri_Unit_InterfaceGeometry.cpp
 *
 * Out-of-line implementations for the shared test-support interface-geometry
 * fixtures declared in Akri_Unit_InterfaceGeometry.hpp.  These live in
 * unit_test_utils (built into krino_unit_test_utils_lib) so every per-library
 * test executable that uses the fixtures resolves them, rather than one test
 * translation unit.
 */
#include <Akri_Unit_InterfaceGeometry.hpp>

#include <Akri_Edge.hpp>
#include <Akri_Intersection_Points.hpp>
#include <Akri_Sign.hpp>

#include <stk_util/util/SortAndUnique.hpp>

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

void PhasePerElementInterfaceGeometry::set_element_phases(const stk::mesh::BulkData & mesh, const std::vector<stk::mesh::EntityId> & elemIds, const std::vector<int> & elemPhases)
{
  STK_ThrowRequire(elemIds.size() == elemPhases.size());
  for (size_t n=0; n<elemIds.size(); ++n)
  {
    stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, elemIds[n]);
    if (mesh.is_valid(elem))
      myUncutElementPhases[elem] = elemPhases[n];
  }
}

}
