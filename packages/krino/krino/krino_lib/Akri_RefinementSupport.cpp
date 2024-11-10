/*
 * Akri_RefinementSupport.cpp
 *
 *  Created on: Feb 1, 2023
 *      Author: drnoble
 */
#include <Akri_RefinementSupport.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_Phase_Support.hpp>

namespace krino {

RefinementSupport &
RefinementSupport::get(stk::mesh::MetaData & meta)
{
  RefinementSupport * support = const_cast<RefinementSupport *>(meta.get_attribute<RefinementSupport>());
  if (NULL == support)
  {
    support = new RefinementSupport(meta);
    meta.declare_attribute_with_delete<RefinementSupport>(support);
  }
  return *support;
}

RefinementSupport &
RefinementSupport::get(const stk::mesh::MetaData & meta)
{
  RefinementSupport * support = const_cast<RefinementSupport *>(meta.get_attribute<RefinementSupport>());
  STK_ThrowRequireMsg(nullptr != support, "Could not find RefinementSupport attribute on MetaData.");
  return *support;
}

stk::mesh::Selector RefinementSupport::do_not_refine_or_unrefine_selector(const stk::mesh::MetaData & meta)
{
  if (!CDFEM_Support::is_active(meta))
  {
    stk::mesh::Selector emptySelector;
    return emptySelector;
  }

  const CDFEM_Support & cdfemSupport = CDFEM_Support::get(meta);
  const stk::mesh::Selector parent_or_child_selector =
      cdfemSupport.get_child_part() | cdfemSupport.get_parent_part();
  const stk::mesh::Selector decomposed_blocks_selector =
      krino::Phase_Support::get(cdfemSupport.get_mesh_meta()).get_all_decomposed_blocks_selector();
  const stk::mesh::Selector do_not_refine_selector = (!decomposed_blocks_selector) | parent_or_child_selector;
  return do_not_refine_selector;
}

RefinementSupport::RefinementSupport(stk::mesh::MetaData & meta)
  : myMeta(meta),
    myTimer("Noninterface Conforming Adapt", sierra::Diag::sierraTimer())
{

}

stk::mesh::Selector RefinementSupport::get_do_not_refine_or_unrefine_selector() const
{
  return do_not_refine_or_unrefine_selector(myMeta);
}

void
RefinementSupport::activate_interface_refinement(int minimumLevel, int maximumLevel)
{
  STK_ThrowRequireMsg(my_interface_minimum_refinement_level == 0 && my_interface_maximum_refinement_level == 0,
      "Interface refinement levels should only be specified once.");
  STK_ThrowRequireMsg(maximumLevel >= minimumLevel || maximumLevel == 0,
      "Maximum interface refinement level must be greater than or equal to the minimum interface refinement level or left unspecified.");
  if (maximumLevel == 0) maximumLevel = minimumLevel;

  my_interface_minimum_refinement_level = minimumLevel;
  my_interface_maximum_refinement_level = maximumLevel;

  setup_refinement_node_marker();

  if (maximumLevel > 0)
    CDFEM_Support::get(myMeta).set_global_ids_are_NOT_parallel_consistent();
}

void
RefinementSupport::activate_nonconformal_adaptivity(const int numLevels)
{
  if (numLevels < my_interface_maximum_refinement_level)
  {
    krinolog << "Ignoring request to activate " << numLevels << " of CDFEM nonconformal adaptivity because a maximum of " << my_interface_maximum_refinement_level << " have already been activated." << stk::diag::dendl;
    return;
  }

  my_interface_minimum_refinement_level = numLevels;
  my_interface_maximum_refinement_level = numLevels;

  setup_refinement_node_marker();

  if (numLevels > 0)
    CDFEM_Support::get(myMeta).set_global_ids_are_NOT_parallel_consistent();
}

void
RefinementSupport::setup_refinement_node_marker()
{
  myNonInterfaceConformingRefinementNodeMarkerField = AuxMetaData::get(myMeta).register_field("REFINEMENT_NODE_MARKER", FieldType::INTEGER, stk::topology::NODE_RANK, 1, 1, myMeta.universal_part());
}

void
RefinementSupport::activate_nonconformal_adapt_target_count(const uint64_t target_count)
{
  my_nonconformal_adapt_target_element_count = target_count;
  my_nonconformal_adapt_indicator_name = "CDFEM_ADAPTIVITY_ERROR_INDICATOR";

  AuxMetaData::get(myMeta).register_field(my_nonconformal_adapt_indicator_name,
      FieldType::REAL,
      stk::topology::ELEMENT_RANK,
      1,
      1,
      myMeta.universal_part());
}


}


