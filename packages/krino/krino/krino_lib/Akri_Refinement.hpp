#ifndef KRINO_KRINO_KRINO_LIB_AKRI_REFINEMENT_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_REFINEMENT_HPP_
#include <string>
#include <vector>
#include <Akri_Edge.hpp>
#include <stk_mesh/base/Types.hpp>
#include "Akri_EntityIdPool.hpp"
#include <stk_mesh/base/Entity.hpp>
#include "Akri_FieldRef.hpp"
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

  static Refinement & get(const stk::mesh::MetaData & meta);
  static Refinement & create(stk::mesh::MetaData & meta, stk::mesh::Part * activePart=nullptr); // must be called before calling get

  Refinement(stk::mesh::MetaData & meta, stk::mesh::Part * activePart);
  Refinement ( const Refinement & ) = delete;
  Refinement & operator= ( const Refinement & ) = delete;

  stk::mesh::Part & child_part() const;
  stk::mesh::Part & parent_part() const;
  stk::mesh::Part & refined_edge_node_part() const;
  bool is_parent(const stk::mesh::Bucket & bucket) const;
  bool is_parent(const stk::mesh::Entity elem) const;
  bool is_a_partially_refined_parent_element(const stk::mesh::Entity parentElem) const;
  bool is_child(const stk::mesh::Bucket & bucket) const;
  bool is_child(const stk::mesh::Entity elem) const;
  stk::mesh::Entity get_parent(const stk::mesh::Entity elem) const;
  bool is_refined_edge_node(const stk::mesh::Entity node) const;
  std::array<stk::mesh::Entity,2> get_edge_parent_nodes(const stk::mesh::Entity edgeNode) const;
  std::tuple<const uint64_t *,unsigned> get_child_ids_and_num_children_when_fully_refined(const stk::mesh::Entity elem) const;
  unsigned get_num_children(const stk::mesh::Entity elem) const;
  unsigned get_num_children_when_fully_refined(const stk::mesh::Entity elem) const;
  void fill_children(const stk::mesh::Entity elem, std::vector<stk::mesh::Entity> & children) const;
  std::vector<stk::mesh::Entity> get_children(const stk::mesh::Entity elem) const;
  stk::mesh::Entity get_edge_child_node(const Edge edge) const { return myNodeRefiner.get_edge_child_node(edge); }
  size_t get_num_edges_to_refine() const { return myNodeRefiner.get_num_edges_to_refine(); }
  void mark_edge_for_refinement(const Edge & edge)  { return myNodeRefiner.mark_edge_for_refinement(edge); }

  void do_refinement(const EdgeMarkerInterface & edgeMarker);

  // public for unit testing
  void find_edges_to_refine(const EdgeMarkerInterface & edgeMarker);

  // Currently only for unit testing
  void fully_unrefine_mesh();

private:
  typedef std::tuple<stk::topology,stk::mesh::PartVector,stk::mesh::EntityVector> BucketData;

  size_t count_new_child_elements(const EdgeMarkerInterface & edgeMarker, const std::vector<BucketData> & bucketsData) const;
  void refine_elements_with_refined_edges_and_store_sides_to_create(const EdgeMarkerInterface & edgeMarker, const std::vector<BucketData> & bucketsData, std::vector<SideDescription> & sideRequests, std::vector<stk::mesh::Entity> & elementsToDelete);
  void declare_refinement_parts();
  void declare_refinement_fields();

  std::vector<BucketData> get_buckets_data_for_candidate_elements_to_refine(const EdgeMarkerInterface & edgeMarker) const;
  FieldRef get_child_element_ids_field(const unsigned numChildWhenFullyRefined) const;

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

  stk::mesh::MetaData & myMeta;
  NodeRefiner myNodeRefiner;
  EntityIdPool myEntityIdPool;

  stk::mesh::Part * myActivePart {nullptr};
  stk::mesh::Part * myParentPart {nullptr};
  stk::mesh::Part * myChildPart {nullptr};
  stk::mesh::Part * myRefinedEdgeNodePart {nullptr};

  FieldRef myCoordsField;
  FieldRef myParentElementIdField;
  FieldRef myChildElementIds4Field;
  FieldRef myChildElementIds8Field;
  FieldRef myRefinedEdgeNodeParentIdsField;
};

} // namespace krino
#endif /* KRINO_KRINO_KRINO_LIB_AKRI_REFINEMENT_HPP_ */
