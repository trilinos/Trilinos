#ifndef KRINO_KRINO_REFINEMENT_AKRI_EDGEMARKERUTILS_HPP_
#define KRINO_KRINO_REFINEMENT_AKRI_EDGEMARKERUTILS_HPP_

#include <vector>

#include <stk_mesh/base/Entity.hpp>
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>

namespace krino {

class Refinement;
struct Edge;
class NodeRefiner;
struct StkMeshEntities;

class ElementEdgeCaseIds
{
public:
  void set(const int preAdaptCaseId, const int postAdaptCaseId) { myPreAdaptCaseId=preAdaptCaseId; myPostAdaptCaseId=postAdaptCaseId; }
  bool has_changed() const { return myPreAdaptCaseId != myPostAdaptCaseId; }
  void clear() { myPreAdaptCaseId=0; myPostAdaptCaseId=0; }
  int pre_adapt_case_id() const { STK_ThrowAssert(has_changed()); return myPreAdaptCaseId; }
  int post_adapt_case_id() const { STK_ThrowAssert(has_changed()); return myPostAdaptCaseId; }
private:
  int myPreAdaptCaseId{0};
  int myPostAdaptCaseId{0};
};

namespace EdgeMarkerUtils {

bool node_is_parent_node(const StkMeshEntities & parentNodes, const stk::mesh::Entity node);

bool is_any_child_invalid_or_a_parent(const stk::mesh::BulkData & mesh,
    const Refinement & refinement,
    const std::vector<stk::mesh::Entity> & childElements);

std::vector<stk::mesh::Entity> get_child_nodes_that_are_not_parent_nodes(const stk::mesh::BulkData & mesh,
    const stk::mesh::Entity parentElem,
    const std::vector<stk::mesh::Entity> & childElems);

std::vector<stk::mesh::Entity> get_child_nodes_that_are_not_parent_nodes(const stk::mesh::BulkData & mesh,
    const Refinement & refinement,
    const stk::mesh::Entity parentElem);

bool is_element_a_parent_but_not_a_grandparent(const stk::mesh::BulkData & mesh,
    const Refinement & refinement,
    const stk::mesh::Entity elem);

void fill_refinement_caseIds_and_refined_edge_nodes_if_changed(const stk::mesh::BulkData & mesh,
    const Refinement & refinement,
    const NodeRefiner & nodeRefiner,
    const stk::mesh::Entity elem,
    const bool isElemCandidateForRefinement,
    const stk::topology & elemTopology,
    ElementEdgeCaseIds & elementEdgeCaseIds,
    std::vector<stk::mesh::Entity> & elemEdgeChildNodes);

void fill_unrefinement_caseIds_and_refined_edge_nodes_if_changed(const stk::mesh::BulkData & mesh,
    const Refinement & refinement,
    const NodeRefiner & nodeRefiner,
    const stk::mesh::Entity elem,
    const bool isElemCandidateForUnrefinement,
    const stk::topology & elemTopology,
    ElementEdgeCaseIds & elementEdgeCaseIds,
    std::vector<stk::mesh::Entity> & elemEdgeChildNodes);
}

}

#endif /* KRINO_KRINO_REFINEMENT_AKRI_EDGEMARKERUTILS_HPP_ */
