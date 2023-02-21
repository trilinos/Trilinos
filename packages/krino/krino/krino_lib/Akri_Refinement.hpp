#ifndef KRINO_KRINO_KRINO_LIB_AKRI_REFINEMENT_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_REFINEMENT_HPP_
#include <string>
#include <vector>
#include <Akri_Edge.hpp>
#include <stk_mesh/base/Types.hpp>
#include "Akri_EntityIdPool.hpp"
#include <stk_mesh/base/Entity.hpp>
#include <stk_math/StkVector.hpp>
#include "Akri_NodeRefiner.hpp"

namespace stk { namespace mesh { class MetaData; } }
namespace stk { class topology; }

namespace krino {

struct SideDescription;
class EdgeMarkerInterface;

class Refinement
{
public:
  enum RefinementMarker
  {
    COARSEN = -1,
    NOTHING = 0,
    REFINE = 1
  };
  Refinement(stk::mesh::MetaData & meta, stk::mesh::Part * activePart, const bool force64Bit, const bool assert32Bit);
  Refinement(stk::mesh::MetaData & meta, stk::mesh::Part * activePart);
  Refinement(stk::mesh::MetaData & meta);
  Refinement ( const Refinement & ) = delete;
  Refinement & operator= ( const Refinement & ) = delete;

  static unsigned get_num_children_when_fully_refined(const stk::topology elementTopology);

  stk::mesh::Part & child_part() const;
  stk::mesh::Part & parent_part() const;
  stk::mesh::Part & refined_edge_node_part() const;
  bool is_parent(const stk::mesh::Bucket & bucket) const;
  bool is_parent(const stk::mesh::Entity elem) const;
  bool is_this_parent_element_partially_refined(const stk::mesh::Entity parentElem) const;
  bool is_child(const stk::mesh::Bucket & bucket) const;
  bool is_child(const stk::mesh::Entity elem) const;
  int refinement_level(const stk::mesh::Entity elem) const;
  stk::mesh::Entity get_parent(const stk::mesh::Entity elem) const;
  std::pair<stk::mesh::EntityId,int> get_parent_id_and_parallel_owner_rank(const stk::mesh::Entity child) const;
  bool is_refined_edge_node(const stk::mesh::Entity node) const;
  std::array<stk::mesh::Entity,2> get_edge_parent_nodes(const stk::mesh::Entity edgeNode) const;
  std::tuple<const uint64_t *,unsigned> get_child_ids_and_num_children_when_fully_refined(const stk::mesh::Entity elem) const;
  unsigned get_num_children(const stk::mesh::Entity elem) const;
  unsigned get_num_children_when_fully_refined(const stk::mesh::Entity elem) const;
  void fill_children(const stk::mesh::Entity elem, std::vector<stk::mesh::Entity> & children) const;
  void fill_child_element_ids(const stk::mesh::Entity elem, std::vector<stk::mesh::EntityId> & childElemIds) const;
  std::vector<stk::mesh::Entity> get_children(const stk::mesh::Entity elem) const;
  stk::mesh::Entity get_edge_child_node(const Edge edge) const { return myNodeRefiner.get_edge_child_node(edge); }
  size_t get_num_edges_to_refine() const { return myNodeRefiner.get_num_edges_to_refine(); }

  void do_refinement(const EdgeMarkerInterface & edgeMarker);
  void do_uniform_refinement(const int numUniformRefinementLevels);

  void restore_after_restart();
  void parallel_sync_child_element_ids_fields();

  // public for unit testing
  void find_edges_to_refine(const EdgeMarkerInterface & edgeMarker);
  bool have_any_hanging_refined_nodes() const;

  // Currently only for unit testing
  void fully_unrefine_mesh();

private:
  typedef std::tuple<stk::topology,stk::mesh::PartVector,stk::mesh::EntityVector> BucketData;
  typedef std::pair<stk::mesh::Entity, stk::mesh::EntityId> ParentAndChildId;

  size_t count_new_child_elements(const EdgeMarkerInterface & edgeMarker, const std::vector<BucketData> & bucketsData) const;
  void refine_elements_with_refined_edges_and_store_sides_to_create(const EdgeMarkerInterface & edgeMarker, const std::vector<BucketData> & bucketsData, std::vector<SideDescription> & sideRequests, std::vector<stk::mesh::Entity> & elementsToDelete);
  void declare_refinement_parts();
  void declare_refinement_fields();

  stk::mesh::PartVector get_parts_for_child_elements(const stk::mesh::Bucket & parentBucket) const;
  std::vector<BucketData> get_buckets_data_for_candidate_elements_to_refine(const EdgeMarkerInterface & edgeMarker) const;
  stk::mesh::Field<uint64_t> & get_child_element_ids_field(const unsigned numChildWhenFullyRefined) const;

  void set_parent_id(const stk::mesh::Entity elem, const stk::mesh::EntityId parentElemId) const;
  stk::mesh::EntityId * get_child_element_ids(const unsigned numChildWhenFullyRefined, const stk::mesh::Entity parent) const;
  stk::mesh::EntityId * get_child_element_ids(const unsigned numChildWhenFullyRefined, const stk::mesh::Bucket & bucket) const;
  stk::math::Vector3d get_coordinates(const stk::mesh::Entity node, const int dim=3) const;

  bool locally_have_any_hanging_refined_nodes() const;
  void create_refined_nodes_elements_and_sides(const EdgeMarkerInterface & edgeMarker);
  void create_another_layer_of_refined_elements_and_sides_to_eliminate_hanging_nodes(const EdgeMarkerInterface & edgeMarker);
  bool do_unrefinement(const EdgeMarkerInterface & edgeMarker);
  void mark_already_refined_edges();

  stk::mesh::EntityId get_parent_id(const stk::mesh::Entity elem) const;
  int get_originating_processor_for_parent_element(const stk::mesh::Entity elem) const;
  void set_originating_processor_for_parent_element(const stk::mesh::Entity elem, const int originatingProc) const;
  void set_refinement_level(const stk::mesh::Entity elem, const int refinementLevel) const;
  void set_parent_parts_and_parent_child_relation_fields(const stk::mesh::Entity parentElement, const std::vector<stk::mesh::Entity> & childElements, const unsigned numChildWhenFullyRefined);
  void refine_element_if_it_has_refined_edges_and_append_sides_to_create(const stk::topology & elemTopology,
      const stk::mesh::PartVector & childParts,
      const stk::mesh::Entity elem,
      const std::vector<stk::mesh::Entity> & elemChildEdgeNodes,
      std::vector<SideDescription> & sideRequests,
      std::vector<stk::mesh::Entity> & elementsToDelete);
  void refine_tri_3_and_append_sides_to_create(const stk::mesh::PartVector & childParts, const stk::mesh::Entity parentElem, const std::vector<stk::mesh::Entity> & elemChildEdgeNodes, const int caseId, std::vector<SideDescription> & sideRequests);
  void refine_tet_4_and_append_sides_to_create(const stk::mesh::PartVector & childParts, const stk::mesh::Entity parentElem, const std::vector<stk::mesh::Entity> & elemChildEdgeNodes, const int caseId, std::vector<SideDescription> & sideRequests);
  stk::mesh::PartVector get_parts_for_new_refined_edge_nodes() const;
  void remove_parent_parts(const std::vector<stk::mesh::Entity> & elements);
  void fill_parents_and_children_and_parents_with_off_proc_child(std::vector<stk::mesh::Entity> & parents, std::vector<stk::mesh::Entity> & children, std::vector<ParentAndChildId> & parentsAndOffProcChildId) const;
  void restore_parent_and_child_element_parts(const std::vector<stk::mesh::Entity> & parents, const std::vector<stk::mesh::Entity> & children);
  void restore_child_element_ids_field(const std::vector<ParentAndChildId> & parentsAndOffProcChildId);
  void add_child_to_parent(const stk::mesh::EntityId childId, const stk::mesh::Entity parent);
  std::vector<int> get_originating_procs_for_elements(const std::vector<stk::mesh::Entity> & elements) const;
  void respect_originating_proc_for_parents_modified_by_unrefinement(const std::vector<stk::mesh::Entity> & parentsModifiedByUnrefinement, const std::vector<int> & originatingProcForParentsModifiedByUnrefinement);

  stk::mesh::MetaData & myMeta;
  bool myForce64Bit;
  bool myAssert32Bit;
  NodeRefiner myNodeRefiner;
  EntityIdPool myEntityIdPool;

  stk::mesh::Part * myActivePart {nullptr};
  stk::mesh::Part * myParentPart {nullptr};
  stk::mesh::Part * myChildPart {nullptr};
  stk::mesh::Part * myRefinedEdgeNodePart {nullptr};

  const stk::mesh::Field<double> * myCoordsField{nullptr};
  stk::mesh::Field<int> * myRefinementLevelField{nullptr};
  stk::mesh::Field<uint64_t> * myParentElementIdField{nullptr};
  stk::mesh::Field<uint64_t> * myChildElementIds4Field{nullptr};
  stk::mesh::Field<uint64_t> * myChildElementIds8Field{nullptr};
  stk::mesh::Field<uint64_t> * myRefinedEdgeNodeParentIdsField{nullptr};
  stk::mesh::Field<int> * myOriginatingProcForParentElementField{nullptr};
};

} // namespace krino
#endif /* KRINO_KRINO_KRINO_LIB_AKRI_REFINEMENT_HPP_ */
