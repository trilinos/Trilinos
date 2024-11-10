/*
 * Akri_RefinementSupport.hpp
 *
 *  Created on: Feb 1, 2023
 *      Author: drnoble
 */

#ifndef KRINO_KRINO_KRINO_LIB_AKRI_REFINEMENTSUPPORT_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_REFINEMENTSUPPORT_HPP_
#include <stk_util/diag/Timer.hpp>
#include <Akri_FieldRef.hpp>

namespace stk { namespace mesh { class MetaData; } }

namespace krino {

class RefinementInterface;

class RefinementSupport {
public:
  static RefinementSupport & get(stk::mesh::MetaData & meta);
  static RefinementSupport & get(const stk::mesh::MetaData & meta);
  static bool use_nonconformal_adaptivity(stk::mesh::MetaData & meta) { RefinementSupport & refinementSupport = get(meta); return refinementSupport.get_interface_maximum_refinement_level() > 0; }
  static stk::mesh::Selector do_not_refine_or_unrefine_selector(const stk::mesh::MetaData & meta);

  void set_initial_refinement_levels(int levels) { my_initial_refinement_levels = levels; }
  int get_initial_refinement_levels() const { return my_initial_refinement_levels; }

  void activate_interface_refinement(int minimum_level, int maximum_level);
  void activate_nonconformal_adaptivity(const int num_levels);

  void activate_nonconformal_adapt_target_count(uint64_t val);
  uint64_t get_nonconformal_adapt_target_count() const { return my_nonconformal_adapt_target_element_count; }
  int get_interface_minimum_refinement_level() const { return my_interface_minimum_refinement_level; }
  int get_interface_maximum_refinement_level() const { return my_interface_maximum_refinement_level; }
  void set_post_adapt_refinement_levels(int levels) { my_post_adapt_uniform_refinement_levels = levels; }
  int get_post_adapt_refinement_levels() const { return my_post_adapt_uniform_refinement_levels; }

  FieldRef get_nonconforming_refinement_node_marker_field() const { return myNonInterfaceConformingRefinementNodeMarkerField; }
  const std::string & get_nonconformal_adapt_indicator_name() const { return my_nonconformal_adapt_indicator_name; }

  void set_non_interface_conforming_refinement(RefinementInterface & refinement) { myNonInterfaceConformingRefinement = &refinement; }
  bool has_non_interface_conforming_refinement() const { return nullptr != myNonInterfaceConformingRefinement; }
  RefinementInterface & get_non_interface_conforming_refinement() const { return *myNonInterfaceConformingRefinement; }
  stk::mesh::Selector get_do_not_refine_or_unrefine_selector() const;

  void do_nearby_refinement_before_interface_refinement(bool flag) { myFlagDoNearbyRefinementBeforeInterfaceRefinement = flag; }
  bool do_nearby_refinement_before_interface_refinement() const { return myFlagDoNearbyRefinementBeforeInterfaceRefinement; }

  static bool is_valid_interval(const std::array<double,2> & refinementInterval) { return refinementInterval[0] <= refinementInterval[1]; }
  void set_refinement_interval(const std::array<double,2> & refinementInterval) { myRefinementInterval = refinementInterval; }
  bool has_refinement_interval() const { return is_valid_interval(myRefinementInterval); }
  const std::array<double,2> get_refinement_interval() const { return myRefinementInterval; }

  stk::diag::Timer & get_timer() const { return myTimer; }

private:
  void setup_refinement_node_marker();

private:
  RefinementSupport(stk::mesh::MetaData & meta);
  RefinementSupport(RefinementSupport const&) = delete;
  RefinementSupport& operator=(RefinementSupport const&) = delete;

  stk::mesh::MetaData & myMeta;
  int my_initial_refinement_levels{0};
  int my_interface_minimum_refinement_level{0};
  int my_interface_maximum_refinement_level{0};
  int my_post_adapt_uniform_refinement_levels{0};
  uint64_t my_nonconformal_adapt_target_element_count{0};
  std::array<double,2> myRefinementInterval{1.,-1.}; // bad initial interval
  FieldRef myNonInterfaceConformingRefinementNodeMarkerField;
  std::string my_nonconformal_adapt_indicator_name;
  RefinementInterface * myNonInterfaceConformingRefinement{nullptr};
  bool myFlagDoNearbyRefinementBeforeInterfaceRefinement{false};
  mutable stk::diag::Timer myTimer;
};

}


#endif /* KRINO_KRINO_KRINO_LIB_AKRI_REFINEMENTSUPPORT_HPP_ */
