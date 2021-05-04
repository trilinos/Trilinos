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

#include "MtoNTransientFieldTransferById.hpp"
#include "stk_mesh/base/FieldParallel.hpp"
#include "stk_mesh/base/Field.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_io/StkIoUtils.hpp"
#include "Ioss_Region.h"

namespace stk {
namespace transfer_utils {

MtoNTransientFieldTransferById::MtoNTransientFieldTransferById(stk::io::StkMeshIoBroker& inputBroker, unsigned numSubDomain)
  : m_inputBroker(inputBroker),
    m_numSubDomain(numSubDomain)
{
  for (const std::string &name : inputBroker.meta_data().entity_rank_names()) {
    stk::mesh::EntityRank rank = inputBroker.meta_data().entity_rank(name);
    m_entityRanks.push_back(rank);
  }
}

MtoNTransientFieldTransferById::MtoNTransientFieldTransferById(stk::io::StkMeshIoBroker &inputBroker, unsigned numSubDomain,
                                                               const std::vector<stk::mesh::EntityRank> &entityRanks)
  : m_inputBroker(inputBroker),
    m_numSubDomain(numSubDomain),
    m_entityRanks(entityRanks)
{
}

MtoNTransientFieldTransferById::~MtoNTransientFieldTransferById()
{
  for (auto & writerPair : m_subdomainWriters) {
    SubdomainWriterBase * subdomainWriter = writerPair.second;
    delete subdomainWriter;
  }

  for (auto & transferPair : m_subdomainTransfers) {
    std::vector<TransientTransferByIdForRank *> & subdomainTransfers = transferPair.second;
    for (TransientTransferByIdForRank * transfer : subdomainTransfers) {
      delete transfer;
    }
  }
}

void MtoNTransientFieldTransferById::initialize_transfer(unsigned subdomain)
{
  for (stk::mesh::EntityRank rank : m_entityRanks) {
    if (stk::io::get_transient_fields(m_inputBroker.meta_data(), rank).size() > 0) {
      stk::mesh::MetaData& outputMeta = m_subdomainWriters[subdomain]->get_meta_data();
      TransientTransferByIdForRank *transfer = new TransientTransferByIdForRank(m_inputBroker.meta_data(),
                                                                                outputMeta,
                                                                                rank);
      transfer->initialize();
      m_subdomainTransfers[subdomain].push_back(transfer);
    }
  }
}

void MtoNTransientFieldTransferById::setup_subdomain(stk::mesh::BulkData& outputBulk, const std::string &filename,
                                                     unsigned subdomain, const stk::io::EntitySharingInfo& nodeSharingInfo,
                                                     int global_num_nodes, int global_num_elems)
{
  SubdomainWriterBase * subdomainWriter = nullptr;

  if (is_valid_subdomain(subdomain)) {
    subdomainWriter = new SubdomainWriter(m_inputBroker, &outputBulk, nodeSharingInfo);
  }
  else {
    subdomainWriter = new EmptySubdomainWriter(m_inputBroker, &outputBulk);
  }

  subdomainWriter->setup_output_file(filename, subdomain, m_numSubDomain, global_num_nodes, global_num_elems);
  m_subdomainWriters[subdomain] = subdomainWriter;

  initialize_transfer(subdomain);
}

void MtoNTransientFieldTransferById::transfer_transient_data(unsigned subdomain)
{
  stk::mesh::BulkData& bulk = m_subdomainWriters[subdomain]->get_bulk_data();
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  for (TransientTransferByIdForRank *transfer : m_subdomainTransfers[subdomain]) {
    transfer->do_transfer();
  }

  stk::mesh::FieldVector fieldVector = stk::io::get_transient_fields(meta);
  std::vector<const stk::mesh::FieldBase*> transientFields;
  transientFields.assign(fieldVector.begin(), fieldVector.end());

  stk::mesh::copy_owned_to_shared(bulk, transientFields); // FIXME: Do we need this?
}

void MtoNTransientFieldTransferById::transfer_and_write_transient_data(unsigned subdomain)
{
  SubdomainWriterBase & subdomainWriter = *m_subdomainWriters.at(subdomain);
  subdomainWriter.write_mesh();

  std::vector<double> inputTimeSteps = m_inputBroker.get_time_steps();

  for (double time : inputTimeSteps) {
    m_inputBroker.read_defined_input_fields(time);
    transfer_transient_data(subdomain);
    subdomainWriter.write_transient_data(time);
  }
}

SubdomainWriterBase &
MtoNTransientFieldTransferById::get_subdomain_writer(unsigned subdomain)
{
  return *m_subdomainWriters.at(subdomain);
}

}
}
