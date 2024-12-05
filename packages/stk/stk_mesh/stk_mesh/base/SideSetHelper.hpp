// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions are
 // met:
 //
 //     * Redistributions of source code must retain the above copyright
 //       notice, this list of conditions and the following disclaimer.
 //
 //     * Redistributions in binary form must reproduce the above
 //       copyright notice, this list of conditions and the following
 //       disclaimer in the documentation and/or other materials provided
 //       with the distribution.
 //
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef STK_STK_MESH_STK_MESH_BASE_SIDESETHELPER_HPP_
#define STK_STK_MESH_STK_MESH_BASE_SIDESETHELPER_HPP_

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/SideSetEntry.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_util/util/DebugStream.hpp>
#include <stk_mesh/baseImpl/elementGraph/ParallelInfoForGraph.hpp>
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraphImpl.hpp"

namespace stk
{
namespace mesh
{
class BulkData;

struct PartChangeAccumulator
{
  Entity entity;
  OrdinalVector partOrdinals;

  PartChangeAccumulator()
  : entity(Entity()) {}

  PartChangeAccumulator(Entity entity_)
  : entity(entity_) {}

  PartChangeAccumulator(Entity entity_, const OrdinalVector& partOrdinals_)
  : entity(entity_), partOrdinals(partOrdinals_) {}

  bool operator!=(const Entity& entity_) const
  {
    return entity != entity_;
  }

  bool operator<(const Entity& entity_) const
  {
    return entity < entity_;
  }

  bool operator!=(const PartChangeAccumulator &rhs) const
  {
    return entity != rhs.entity;
  }

  bool operator<(const PartChangeAccumulator &rhs) const
  {
    return entity < rhs.entity;
  }
};

using PartChangeAccumulatorVector = std::vector<PartChangeAccumulator>;

struct SideSetHelper
{
  SideSetHelper(BulkData& mesh_, const Selector& activeSelector_, std::ostream *outputStream_ = nullptr)
  : mesh(mesh_)
  , activeSelector(activeSelector_)
  , outputStream(outputStream_)
  , m_modCycleWhenParallelPartInfoUpdated(0)
  { }

  void set_output_stream(std::ostream& ostrm) {
    outputStream = &ostrm;
  }

  void remove_side_entry_from_sideset(SideSet* sideset, const SideSetEntry& entry, std::set<const Part*> *touchedSidesetParts = nullptr);
  void remove_side_entry_from_sidesets(Entity elem, ConnectivityOrdinal ordinal, std::set<const Part*> *touchedSidesetParts = nullptr);
  void remove_side_entry_from_sidesets(SideSetVector& sidesets, Entity elem, ConnectivityOrdinal ordinal, std::set<const Part*> *touchedSidesetParts = nullptr);

  void remove_side_entries_from_sidesets(const Entity entity, std::set<const Part*> *touchedSidesetParts = nullptr);
  void remove_side_entries_from_sidesets(SideSetVector& sidesets, const Entity entity, std::set<const Part*> *touchedSidesetParts = nullptr);

  void remove_element_entries_from_sidesets(const Entity entity, std::set<const Part*> *touchedSidesetParts = nullptr);
  void remove_element_entries_from_sidesets(SideSetVector& sidesets, const Entity entity, std::set<const Part*> *touchedSidesetParts = nullptr);

  void remove_entity_entries_from_sidesets(const Entity entity, std::set<const Part*> *touchedSidesetParts = nullptr);
  void remove_entity_entries_from_sidesets(SideSetVector& sidesets, const Entity entity, std::set<const Part*> *touchedSidesetParts = nullptr);

  void add_element_side_entry_to_sidesets(SideSetSelectorVector& sidesets, Entity elem, ConnectivityOrdinal ordinal);
  void add_element_side_entry_to_sideset(SideSetSelector& sidesetSelector, Entity elem, ConnectivityOrdinal ordinal);
  void add_side_entries_to_sidesets(SideSetSelectorVector& sidesets, const Entity entity);
  void add_side_entries_to_sideset(SideSetSelector& sideset, const Entity entity);

  void add_coincident_side_entries_to_sideset(SideSetSelector& sideset, const Entity entity);
  void add_coincident_element_side_entry_to_sideset(SideSetSelector& sideset, const Entity side, const Entity element, const ConnectivityOrdinal ordinal);

  void add_sideset_entry_for_element_selected_by_sidesets(Entity entity, SideSetSelectorVector& sidesetsAndSelectors);

  bool element_side_has_coincidence(SideSetSelector& sideset, const Entity side, const Entity element, const ConnectivityOrdinal ordinal);
  bool element_side_has_coincidence_using_connectivity(const Entity side, const Entity element, const ConnectivityOrdinal ordinal);
  bool element_side_has_coincidence_using_elem_elem_graph(const Entity side, const Entity element, const ConnectivityOrdinal ordinal);

  void reset_internal_sideset_detection(bool elemOrSideRankedPartsChanged);
  void warn_internal_sideset_detection();

  void set_warn_about_internal_sideset(bool flag);

  const impl::ParallelPartInfo& get_parallel_part_info() const { return parallelPartInfo; }

private:
  BulkData& mesh;
  Selector activeSelector;
  std::ostream* outputStream = nullptr;

  bool internalSidesetWarningHasBeenIssued = false;
  impl::ParallelPartInfo parallelPartInfo;
  unsigned m_modCycleWhenParallelPartInfoUpdated;
  std::set<unsigned> internalSidesetOrdinals;
  bool warnAboutInternalSidesets = false;

  SideSetHelper();

  void warn_about_internal_sidesets();
  bool graph_edge_can_be_distinguished(const ElemElemGraph& eeGraph, const GraphEdge& graphEdge, const Selector& selector, bool selectorValue);

  void fill_coincident_sideset_entries_for_side(const Entity side, std::vector<SideSetEntry>& entries);
  void fill_coincident_sideset_entries_for_side_using_connectivity(const Entity side, std::vector<SideSetEntry>& entries);
  void fill_coincident_sideset_entries_for_side_using_elem_elem_graph(const Entity side, std::vector<SideSetEntry>& entries);

  bool element_side_can_be_distinguished(const Entity side, const Entity element, const ConnectivityOrdinal ordinal, const SideSetSelector& selector);
  bool element_side_can_be_distinguished_using_connectivity(const Entity side, const Entity element, const ConnectivityOrdinal ordinal, const SideSetSelector& selector);
  bool element_side_can_be_distinguished_using_elem_elem_graph(const Entity side, const Entity element, const ConnectivityOrdinal ordinal, const SideSetSelector& selector);

  bool element_side_has_remote_coincidence_using_elem_elem_graph(Entity element,  ConnectivityOrdinal ordinal);
  bool element_side_has_local_coincidence_using_elem_elem_graph(Entity element,  ConnectivityOrdinal ordinal);
};

inline void issue_internal_sideset_warning(const std::string& sidesetName, std::ostream& ostrm)
{
  ostrm <<"WARNING, Internal sideset ("<<sidesetName<<") detected. STK doesn't support internal sidesets\n"
      <<"(i.e., sidesets between elements where both elements are in the same block)\n"
      <<"Execution will continue but correct results are not guaranteed. Contact sierra-help@sandia.gov"<<std::endl;
}

}
}

#endif /* STK_STK_MESH_STK_MESH_BASE_SIDESETHELPER_HPP_ */
