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
#include "M2NOutputMesh.hpp"
#include <stk_tools/mesh_clone/MeshClone.hpp>
#include <stk_balance/m2n/M2NInputMesh.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>
#include "stk_io/IossBridge.hpp"

namespace stk {
namespace balance {
namespace m2n {

OutputMesh::OutputMesh(const InputMesh& inputMesh,
                       const std::vector<unsigned>& targetSubdomains)
  : m_inputMesh(inputMesh),
    m_targetSubdomains(targetSubdomains),
    m_bulk(m_inputMesh.get_bulk().parallel()),
    m_meta(m_bulk.mesh_meta_data())
{
  clone_input_mesh();
  move_subdomain_to_owning_processor();
}

void
OutputMesh::clone_input_mesh()
{
  const stk::mesh::BulkData& inputBulk = m_inputMesh.get_bulk();

  stk::mesh::Selector subdomainSelector;
  for (unsigned subdomain : m_targetSubdomains) {
    if (is_valid_subdomain(subdomain)) {
      subdomainSelector |= *inputBulk.mesh_meta_data().get_part(m_inputMesh.get_subdomain_part_name(subdomain));
    }
  }

  stk::tools::copy_mesh(inputBulk, subdomainSelector, m_bulk);
}

void
OutputMesh::move_subdomain_to_owning_processor()
{
  stk::mesh::EntityProcVec subdomainDecomp;
  for (unsigned subdomain : m_targetSubdomains) {
    if (is_valid_subdomain(subdomain)) {
      const stk::mesh::Selector& locallyOwnedSubdomain = m_meta.locally_owned_part() &
                                                         *m_meta.get_part(m_inputMesh.get_subdomain_part_name(subdomain));
      for (const stk::mesh::Bucket* bucket : m_bulk.get_buckets(stk::topology::ELEM_RANK, locallyOwnedSubdomain)) {
        for (const stk::mesh::Entity& elem : *bucket) {
          subdomainDecomp.emplace_back(elem, subdomain);
        }
      }
    }
  }

  stk::balance::internal::rebalance(m_bulk, m_inputMesh.get_owner_for_each_final_subdomain(), subdomainDecomp);
}

void
OutputMesh::transfer_and_write()
{
  TransientFieldTransferById transientIo(m_inputMesh.get_io_broker(), m_inputMesh.get_num_output_processors());

  const unsigned mySubdomain = m_targetSubdomains[m_inputMesh.get_bulk().parallel_rank()];
  stk::io::EntitySharingInfo nodeSharingInfo = m_inputMesh.get_node_sharing_info(mySubdomain, m_targetSubdomains);
  m_bulk.switch_to_serial_mesh();

  transientIo.setup_subdomain(m_bulk, m_inputMesh.get_output_file_name(), mySubdomain, nodeSharingInfo,
                              m_inputMesh.get_global_num_nodes(), m_inputMesh.get_global_num_elements());
  transientIo.transfer_and_write_transient_data(mySubdomain);
}

}}}
