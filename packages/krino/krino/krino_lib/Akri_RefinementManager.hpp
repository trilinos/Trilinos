#ifndef KRINO_KRINO_LIB_AKRI_REFINEMENT_MANAGER_HPP_
#define KRINO_KRINO_LIB_AKRI_REFINEMENT_MANAGER_HPP_
#include <vector>

#include <Akri_FieldRef.hpp>
#include "Akri_Refinement.hpp"
#include "stk_util/diag/Timer.hpp"

namespace stk { namespace mesh { struct Entity; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace diag { class Timer; } }

namespace krino {

class AuxMetaData;
class TransitionElementEdgeMarker;
class RefinementManager;

void clear_refinement_marker(const RefinementManager & refinement);
void mark_selected_elements_for_refinement(const RefinementManager & refinement, const stk::mesh::Selector & selector);
void mark_selected_elements_for_refinement(const RefinementManager & refinement, const int current_refinement_level, const int max_refinement_levels, const stk::mesh::Selector & selector);
void mark_elements_for_refinement(const RefinementManager & refinement, const std::vector<stk::mesh::Entity> & elemsToRefine);
void mark_elements_with_given_ids_for_refinement(const stk::mesh::BulkData & mesh, const RefinementManager & refinement, const std::vector<stk::mesh::EntityId> & idsOfElemsToRefine);
void mark_based_on_indicator_field(const stk::mesh::BulkData & mesh,
    const RefinementManager & refinement,
     const std::string & indicatorFieldName,
    const int maxRefinementLevel,
    const int currentRefinementLevel,
    const uint64_t targetElemCount);

void fill_leaf_children(const RefinementManager & refinement, const stk::mesh::Entity entity, std::vector<stk::mesh::Entity> & leaf_children);
void fill_leaf_children(const RefinementManager & refinement, const std::vector<stk::mesh::Entity> & children, std::vector<stk::mesh::Entity> & leaf_children);
void fill_all_children(const RefinementManager & refinement, const std::vector<stk::mesh::Entity> & children, std::vector<stk::mesh::Entity> & all_children);
void fill_all_children(const RefinementManager & refinement, const stk::mesh::Entity entity, std::vector<stk::mesh::Entity> & all_children);

class RefinementManager
{
public:
  static RefinementManager & get(const stk::mesh::MetaData & meta);
  static RefinementManager & create(stk::mesh::MetaData & meta);
  static RefinementManager & create(stk::mesh::MetaData & meta, stk::diag::Timer & timer);
  static RefinementManager & get_or_create(stk::mesh::MetaData & meta);
  static RefinementManager & get_or_create(stk::mesh::MetaData & meta, stk::diag::Timer & timer);
  static bool is_created(const stk::mesh::MetaData & meta);
  static void register_parts_and_fields_via_aux_meta_for_fmwk(stk::mesh::MetaData & meta);

  bool is_child(const stk::mesh::Entity entity) const;
  bool is_parent(const stk::mesh::Entity parent) const;
  bool is_parent_side(const stk::mesh::Entity side) const;
  bool is_transition(const stk::mesh::Entity parent) const;
  stk::mesh::Part & parent_part() const;
  stk::mesh::Part & child_part() const;
  stk::mesh::Entity get_parent(const stk::mesh::Entity elem) const;
  std::pair<stk::mesh::EntityId,int> get_parent_id_and_parallel_owner_rank(const stk::mesh::Entity child) const;
  unsigned get_num_children(const stk::mesh::Entity elem) const;
  void fill_children(const stk::mesh::Entity parent, std::vector<stk::mesh::Entity> & children) const;
  void fill_child_element_ids(const stk::mesh::Entity parent, std::vector<stk::mesh::EntityId> & childElemIds) const;
  void fill_dependents(const stk::mesh::Entity parent, std::vector<stk::mesh::Entity> & dependents) const;
  bool has_rebalance_constraint(const stk::mesh::Entity entity) const;
  void update_element_rebalance_weights_incorporating_parallel_owner_constraints(stk::mesh::Field<double> & elemWtField) const;
  int fully_refined_level(const stk::mesh::Entity elem) const;
  int partially_refined_level(const stk::mesh::Entity elem) const;
  std::string locally_check_leaf_children_have_parents_on_same_proc() const;

  FieldRef get_marker_field_and_sync_to_host() const;
  bool require_post_refinement_fixups() const { return false; };

  bool do_refinement(const int debugLevel = 0);

  bool do_uniform_refinement(const int numUniformRefinementLevels);

  void restore_after_restart();
  void set_marker_field(const std::string & markerFieldName);
  void delete_parent_elements() {myRefinement.delete_parent_elements();};

private:
  RefinementManager(stk::mesh::MetaData & meta, stk::mesh::Part * activePart, const bool force64Bit, const bool assert32Bit, stk::diag::Timer & parent_timer);
  std::array<unsigned, 2>  get_marked_element_counts() const;
  bool is_supported_uniform_refinement_element() const;
  TransitionElementEdgeMarker & setup_marker() const;
  TransitionElementEdgeMarker & get_marker() const;
  stk::mesh::MetaData & myMeta;
  FieldRef myElementMarkerField;
  mutable std::unique_ptr<TransitionElementEdgeMarker> myMarker;
  mutable stk::diag::Timer myRefinementTimer;
  Refinement myRefinement;
};

void check_leaf_children_have_parents_on_same_proc(const stk::ParallelMachine comm, const RefinementManager & refinement);

} // namespace krino

#endif /* KRINO_KRINO_LIB_AKRI_REFINEMENT_MANAGER_HPP_ */
