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

#include "BalanceIO.hpp"
#include "stk_balance/internal/privateDeclarations.hpp"
#include "stk_balance/fixSplitCoincidentElements.hpp"
#include "stk_tools/mesh_clone/MeshClone.hpp"
#include "stk_tools/transfer_utils/TransientFieldTransferById.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"

#include <sys/stat.h> // move us
#include <algorithm>
#include <cerrno>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>

namespace stk {
namespace balance {

void read_mesh_with_auto_decomp(stk::io::StkMeshIoBroker & stkIo,
                                const std::string& meshSpec,
                                stk::mesh::BulkData& bulkData,
                                const stk::balance::BalanceSettings & balanceSettings)
{
  stkIo.set_bulk_data(bulkData);
  stkIo.add_mesh_database(meshSpec, stk::io::READ_MESH);
  stkIo.create_input_mesh();
  stkIo.add_all_mesh_fields_as_input_fields();

  internal::register_internal_fields_and_parts(bulkData, balanceSettings);

  stkIo.populate_bulk_data();

  if (balanceSettings.getVertexWeightMethod() == stk::balance::VertexWeightMethod::FIELD) {
    const int numTimeSteps = stkIo.get_num_time_steps();
    stkIo.read_defined_input_fields_at_step(numTimeSteps, nullptr);
  }

  if(stkIo.check_integer_size_requirements() == 8) {
    bulkData.set_large_ids_flag(true);
  }
}

BalanceIO::BalanceIO(MPI_Comm comm, const BalanceSettings& settings)
  : m_comm(comm),
    m_settings(settings),
    m_inputBulk(stk::mesh::MeshBuilder(m_comm).create()),
    m_inputMeta(m_inputBulk->mesh_meta_data()),
    m_copyBulk(stk::mesh::MeshBuilder(m_comm).create()),
    m_copyMeta(m_copyBulk->mesh_meta_data()),
    m_mesh(nullptr)
{
  m_inputMeta.set_coordinate_field_name(m_settings.getCoordinateFieldName());

}

BalanceMesh& BalanceIO::initial_decomp()
{
  internal::logMessage(m_comm, "Reading mesh and performing initial decomposition");
  m_inputBroker.property_add(Ioss::Property("DECOMPOSITION_METHOD", m_settings.getInitialDecompMethod()));
  read_mesh_with_auto_decomp(m_inputBroker, m_settings.get_input_filename(), *m_inputBulk, m_settings);
  make_mesh_consistent_with_parallel_mesh_rule1(*m_inputBulk);

  if(stk::io::get_transient_fields(m_inputMeta).empty()) {
    m_mesh = std::unique_ptr<BalanceMesh>(new BalanceMesh(*m_inputBulk));
  }
  else {
    internal::logMessage(m_comm, "Copying input mesh to handle transient fields");
    stk::tools::copy_mesh(*m_inputBulk, m_inputMeta.universal_part(), *m_copyBulk);
    m_mesh = std::unique_ptr<BalanceMesh>(new BalanceMesh(*m_copyBulk));
  }

  return *m_mesh;
}

void BalanceIO::write(BalanceMesh& mesh)
{
  stk::io::StkMeshIoBroker outputBroker;
  outputBroker.set_bulk_data(mesh.get_bulk());
  outputBroker.set_attribute_field_ordering_stored_by_part_ordinal(m_inputBroker.get_attribute_field_ordering_stored_by_part_ordinal());
  m_inputBroker.cache_entity_list_for_transient_steps(true);

  stk::transfer_utils::TransientFieldTransferById transfer(m_inputBroker, outputBroker);
  transfer.transfer_and_write_transient_fields(m_settings.get_output_filename());

  internal::logMessage(m_comm, "Finished writing output mesh");
}

} }
