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
//

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off

#include <stk_io/DynamicTopology.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace io {

unsigned DynamicTopologyObserver::count_output_field_variables(stk::mesh::EntityRank rank) const
{
  const std::vector<stk::io::FieldAndName>& definedFields = m_stkIoBroker.get_defined_output_fields(m_outputFileIndex);

  unsigned numFieldsOfRank = 0;
  for(const stk::io::FieldAndName& entry : definedFields) {
    if(entry.type() == rank) {
      numFieldsOfRank++;
    }
  }

  return numFieldsOfRank;
}

DynamicTopologyObserver::DynamicTopologyObserver(StkMeshIoBroker& stkIoBroker, size_t index,
                                                 const Ioss::FileControlOption fileControlOption_)
    : Ioss::DynamicTopologyObserver(nullptr)
    , m_stkIoBroker(stkIoBroker)
    , m_outputFileIndex(index)
    , m_fileControlOption(fileControlOption_)
{

}

void DynamicTopologyObserver::define_model()
{
  // TODO: Need to set is_automatic_restart() and is_restart_requested() on the observer here

  m_stkIoBroker.reset_output_mesh_definition(m_outputFileIndex);
  m_stkIoBroker.define_output_mesh(m_outputFileIndex);
}

void DynamicTopologyObserver::write_model()
{
  m_stkIoBroker.write_output_mesh(m_outputFileIndex);
}

void DynamicTopologyObserver::define_transient()
{
  // TODO: Need to set is_automatic_restart() and is_restart_requested() on the observer here
  m_stkIoBroker.define_output_fields(m_outputFileIndex);
}

Ioss::FileControlOption DynamicTopologyObserver::get_control_option() const
{
  return m_fileControlOption;
}

bool DynamicTopologyObserver::needs_new_output_file() const
{
  // See what type of topology modification has ocurred.
  stk::topology::rank_t sideRank = m_stkIoBroker.meta_data().side_rank();

  if (((get_topology_modification() & Ioss::TOPOLOGY_HADAPT) && get_output_refined_mesh()) ||
      ((get_topology_modification() & Ioss::TOPOLOGY_CREATEFACE) && (count_output_field_variables(sideRank) > 0)) ||
      ((get_topology_modification() & Ioss::TOPOLOGY_CREATEELEM) && (count_output_field_variables(stk::topology::ELEM_RANK) > 0)) ||
      ((get_topology_modification() & Ioss::TOPOLOGY_CREATENODE) && (count_output_field_variables(stk::topology::NODE_RANK) > 0)) ||
      (get_topology_modification() & Ioss::TOPOLOGY_UNKNOWN) ||
      (get_topology_modification() & Ioss::TOPOLOGY_SHUFFLE) ||
      (get_topology_modification() & Ioss::TOPOLOGY_REORDER)) {
    // See if database has been written to since being opened...
    if (get_cumulative_topology_modification() != 0) {
      return true;
    }
  }
  return false;
}

Ioss::FileControlOption get_ioss_file_control_option(FileOption fileOption)
{
  switch(fileOption) {
  default:
  case FileOption::NO_DYNAMIC_TOPOLOGY_FILE_CONTROL:
    return Ioss::FileControlOption::CONTROL_NONE;
    break;
  case FileOption::USE_DYNAMIC_TOPOLOGY_MULTI_FILE:
    return Ioss::FileControlOption::CONTROL_AUTO_MULTI_FILE;
    break;
  case FileOption::USE_DYNAMIC_TOPOLOGY_GROUP_FILE:
    return Ioss::FileControlOption::CONTROL_AUTO_GROUP_FILE;
    break;
  }

  return Ioss::FileControlOption::CONTROL_NONE;
}

void DynamicTopologyObserver::initialize_region()
{
  if (nullptr != m_region) {
    delete_selector_property(*m_region);
    m_region->reset_region();
  }
}

} // namespace io
} // namespace stk
