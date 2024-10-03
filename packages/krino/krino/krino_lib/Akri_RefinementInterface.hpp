/*
 * Akri_RefinementInterface.hpp
 *
 *  Created on: Oct 31, 2022
 *      Author: drnoble
 */

#ifndef KRINO_KRINO_ADAPTIVITY_INTERFACE_AKRI_REFINEMENTINTERFACE_HPP_
#define KRINO_KRINO_ADAPTIVITY_INTERFACE_AKRI_REFINEMENTINTERFACE_HPP_
#include <vector>

#include <Akri_FieldRef.hpp>
#include "Akri_Refinement.hpp"
#include "stk_util/diag/Timer.hpp"

namespace stk { namespace mesh { class Entity; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace diag { class Timer; } }

namespace krino {

class AuxMetaData;
class TransitionElementEdgeMarker;
class RefinementInterface;

void clear_refinement_marker(const RefinementInterface & refinement);
void mark_selected_elements_for_refinement(const RefinementInterface & refinement, const stk::mesh::Selector & selector);
void mark_selected_elements_for_refinement(const RefinementInterface & refinement, const int current_refinement_level, const int max_refinement_levels, const stk::mesh::Selector & selector);
void mark_elements_for_refinement(const RefinementInterface & refinement, const std::vector<stk::mesh::Entity> & elemsToRefine);
void mark_elements_with_given_ids_for_refinement(const stk::mesh::BulkData & mesh, const RefinementInterface & refinement, const std::vector<stk::mesh::EntityId> & idsOfElemsToRefine);
void mark_based_on_indicator_field(const stk::mesh::BulkData & mesh,
    const RefinementInterface & refinement,
     const std::string & indicatorFieldName,
    const int maxRefinementLevel,
    const int currentRefinementLevel,
    const uint64_t targetElemCount);

void fill_leaf_children(const RefinementInterface & refinement, const stk::mesh::Entity entity, std::vector<stk::mesh::Entity> & leaf_children);
void fill_leaf_children(const RefinementInterface & refinement, const std::vector<stk::mesh::Entity> & children, std::vector<stk::mesh::Entity> & leaf_children);
void fill_all_children(const RefinementInterface & refinement, const std::vector<stk::mesh::Entity> & children, std::vector<stk::mesh::Entity> & all_children);
void fill_all_children(const RefinementInterface & refinement, const stk::mesh::Entity entity, std::vector<stk::mesh::Entity> & all_children);

class RefinementInterface
{
public:
  virtual ~RefinementInterface() {}
  virtual bool is_child(const stk::mesh::Entity entity) const = 0;
  virtual bool is_parent(const stk::mesh::Entity parent) const = 0;
  virtual bool is_parent_side(const stk::mesh::Entity side) const = 0;
  virtual bool is_transition(const stk::mesh::Entity parent) const = 0;
  virtual stk::mesh::Part & parent_part() const = 0;
  virtual stk::mesh::Part & child_part() const = 0;
  virtual stk::mesh::Entity get_parent(const stk::mesh::Entity child) const = 0;
  virtual std::pair<stk::mesh::EntityId,int> get_parent_id_and_parallel_owner_rank(const stk::mesh::Entity child) const = 0;
  virtual void fill_children(const stk::mesh::Entity parent, std::vector<stk::mesh::Entity> & children) const = 0;
  virtual void fill_child_element_ids(const stk::mesh::Entity parent, std::vector<stk::mesh::EntityId> & childElemIds) const = 0;
  //Dependents must stay on the same proc as the entity they are dependent on. 
  //How many levels of children count as dependents is implementation specific
  virtual void fill_dependents(const stk::mesh::Entity parent, std::vector<stk::mesh::Entity> & dependents) const = 0;
  //Whether an entity is constrained or free to move in a rebalance
  virtual bool has_rebalance_constraint(const stk::mesh::Entity entity) const = 0;
  virtual void update_element_rebalance_weights_incorporating_parallel_owner_constraints(stk::mesh::Field<double> & elemWtField) const = 0;
  virtual unsigned get_num_children(const stk::mesh::Entity elem) const = 0;
  virtual int fully_refined_level(const stk::mesh::Entity elem) const = 0;
  virtual FieldRef get_marker_field_and_sync_to_host() const = 0;
  virtual bool require_post_refinement_fixups() const = 0;
  virtual std::string locally_check_leaf_children_have_parents_on_same_proc() const = 0;

  virtual bool do_refinement(const int debugLevel = 0) = 0;
  virtual bool do_uniform_refinement(const int numUniformRefinementLevels) = 0;
  virtual void delete_parent_elements() = 0;
};

class KrinoRefinement : public RefinementInterface
{
public:
  static KrinoRefinement & get(const stk::mesh::MetaData & meta);
  static KrinoRefinement & create(stk::mesh::MetaData & meta);
  static KrinoRefinement & create(stk::mesh::MetaData & meta, stk::diag::Timer & timer);
  static KrinoRefinement & get_or_create(stk::mesh::MetaData & meta);
  static KrinoRefinement & get_or_create(stk::mesh::MetaData & meta, stk::diag::Timer & timer);
  static bool is_created(const stk::mesh::MetaData & meta);
  static void register_parts_and_fields_via_aux_meta_for_fmwk(stk::mesh::MetaData & meta);

  virtual ~KrinoRefinement();

  virtual bool is_child(const stk::mesh::Entity entity) const override;
  virtual bool is_parent(const stk::mesh::Entity parent) const override;
  virtual bool is_parent_side(const stk::mesh::Entity side) const override;
  virtual bool is_transition(const stk::mesh::Entity parent) const override;
  virtual stk::mesh::Part & parent_part() const override;
  virtual stk::mesh::Part & child_part() const override;
  virtual stk::mesh::Entity get_parent(const stk::mesh::Entity elem) const override;
  virtual std::pair<stk::mesh::EntityId,int> get_parent_id_and_parallel_owner_rank(const stk::mesh::Entity child) const override;
  virtual unsigned get_num_children(const stk::mesh::Entity elem) const override;
  virtual void fill_children(const stk::mesh::Entity parent, std::vector<stk::mesh::Entity> & children) const override;
  virtual void fill_child_element_ids(const stk::mesh::Entity parent, std::vector<stk::mesh::EntityId> & childElemIds) const override;
  virtual void fill_dependents(const stk::mesh::Entity parent, std::vector<stk::mesh::Entity> & dependents) const override;
  virtual bool has_rebalance_constraint(const stk::mesh::Entity entity) const override;
  virtual void update_element_rebalance_weights_incorporating_parallel_owner_constraints(stk::mesh::Field<double> & elemWtField) const override;
  virtual int fully_refined_level(const stk::mesh::Entity elem) const override;
  int partially_refined_level(const stk::mesh::Entity elem) const;
  virtual std::string locally_check_leaf_children_have_parents_on_same_proc() const override;

  virtual FieldRef get_marker_field_and_sync_to_host() const override;
  virtual bool require_post_refinement_fixups() const override { return false; };

  virtual bool do_refinement(const int debugLevel = 0) override;

  virtual bool do_uniform_refinement(const int numUniformRefinementLevels) override;

  void restore_after_restart();
  void set_marker_field(const std::string & markerFieldName);
  virtual void delete_parent_elements() override {myRefinement.delete_parent_elements();};

private:
  KrinoRefinement(stk::mesh::MetaData & meta, stk::mesh::Part * activePart, const bool force64Bit, const bool assert32Bit, stk::diag::Timer & parent_timer);
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

void check_leaf_children_have_parents_on_same_proc(const stk::ParallelMachine comm, const RefinementInterface & refinement);

} // namespace krino

#endif /* KRINO_KRINO_ADAPTIVITY_INTERFACE_AKRI_REFINEMENTINTERFACE_HPP_ */
