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
#include "M2NDecomposer.hpp"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/setup/M2NParser.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>

namespace stk {
namespace balance {
namespace internal {

M2NDecomposer::M2NDecomposer(stk::mesh::BulkData & bulkData,
                             const stk::balance::BalanceSettings & balanceSettings,
                             const stk::balance::M2NParsedOptions & parsedOptions)
  : m_bulkData(bulkData),
    m_balanceSettings(balanceSettings),
    m_parsedOptions(parsedOptions)
{
}

unsigned
M2NDecomposer::num_required_subdomains_for_each_proc()
{
  const unsigned numInitialSubdomains = m_bulkData.parallel_size();
  const unsigned numFinalSubdomains = m_parsedOptions.targetNumProcs;
  return (numFinalSubdomains / numInitialSubdomains) + (numFinalSubdomains % numInitialSubdomains > 0);
}

stk::mesh::EntityProcVec
M2NDecomposer::get_partition()
{
  stk::mesh::EntityProcVec decomp;
  std::vector<stk::mesh::Selector> selectors = { m_bulkData.mesh_meta_data().universal_part() };
  stk::balance::internal::calculateGeometricOrGraphBasedDecomp(m_bulkData, selectors,
                                                               m_bulkData.parallel(), m_parsedOptions.targetNumProcs,
                                                               m_balanceSettings, decomp);
  return decomp;
}

std::vector<unsigned>
M2NDecomposer::map_new_subdomains_to_original_processors()
{
  const unsigned numInitialSubdomains = m_bulkData.parallel_size();
  const unsigned numFinalSubdomains = m_parsedOptions.targetNumProcs;

  std::vector<unsigned> targetSubdomainsToProc(numFinalSubdomains, std::numeric_limits<unsigned>::max());
  for (unsigned i = 0; i < numFinalSubdomains; ++i) {
    targetSubdomainsToProc[i] = i % numInitialSubdomains;
  }
  return targetSubdomainsToProc;
}

M2NDecomposerNested::M2NDecomposerNested(stk::mesh::BulkData & bulkData,
                                         const stk::balance::BalanceSettings & balanceSettings,
                                         const stk::balance::M2NParsedOptions & parsedOptions)
  : M2NDecomposer(bulkData, balanceSettings, parsedOptions)
{
  const int numInitialSubdomains = m_bulkData.parallel_size();
  const int numFinalSubdomains = m_parsedOptions.targetNumProcs;
  ThrowRequireMsg((numFinalSubdomains % numInitialSubdomains) == 0,
                  "Final subdomains (" << numFinalSubdomains << ") must be an integer multiple of initial subdomains ("
                  << numInitialSubdomains << ")");

  m_numFinalSubdomainsPerProc = numFinalSubdomains / numInitialSubdomains;
}

stk::mesh::EntityProcVec
M2NDecomposerNested::get_partition()
{
  declare_all_initial_subdomain_parts();
  move_entities_into_initial_subdomain_part();

  stk::mesh::EntityProcVec decomp;
  std::vector<stk::mesh::Selector> selectors = { *m_bulkData.mesh_meta_data().get_part(get_initial_subdomain_part_name(m_bulkData.parallel_rank())) };
  const unsigned numLocalElems = stk::mesh::count_entities(m_bulkData, stk::topology::ELEM_RANK, selectors[0]);

  if (numLocalElems > 0) {
    stk::balance::internal::calculateGeometricOrGraphBasedDecomp(m_bulkData, selectors,
                                                                 MPI_COMM_SELF, m_numFinalSubdomainsPerProc,
                                                                 m_balanceSettings, decomp);
  }

  for (stk::mesh::EntityProc & entityProc : decomp) {
    int & targetProc = entityProc.second;
    targetProc += m_bulkData.parallel_rank()*m_numFinalSubdomainsPerProc;
  }

  return decomp;
}

stk::mesh::EntityProcVec
M2NDecomposerNested::get_partition_for_subdomain(int subdomainId)
{
  stk::mesh::EntityProcVec decomp;
  std::vector<stk::mesh::Selector> selectors = { *m_bulkData.mesh_meta_data().get_part(get_initial_subdomain_part_name(subdomainId)) };
  stk::balance::internal::calculateGeometricOrGraphBasedDecomp(m_bulkData, selectors,
                                                               m_bulkData.parallel(), m_numFinalSubdomainsPerProc,
                                                               m_balanceSettings, decomp);
  return decomp;
}

std::vector<unsigned>
M2NDecomposerNested::map_new_subdomains_to_original_processors()
{
  const unsigned numFinalSubdomains = m_parsedOptions.targetNumProcs;

  std::vector<unsigned> targetSubdomainsToProc(numFinalSubdomains, std::numeric_limits<unsigned>::max());
  for (unsigned i = 0; i < numFinalSubdomains; ++i) {
    targetSubdomainsToProc[i] = i / m_numFinalSubdomainsPerProc;
  }
  return targetSubdomainsToProc;
}

void
M2NDecomposerNested::declare_all_initial_subdomain_parts()
{
  for (int i = 0; i < m_bulkData.parallel_size(); ++i) {
    std::string partNameForSubdomain = get_initial_subdomain_part_name(i);
    m_bulkData.mesh_meta_data().declare_part(partNameForSubdomain, stk::topology::ELEMENT_RANK);
  }
}

void
M2NDecomposerNested::move_entities_into_initial_subdomain_part()
{
  stk::mesh::Part * subdomainPart = m_bulkData.mesh_meta_data().get_part(get_initial_subdomain_part_name(m_bulkData.parallel_rank()));
  stk::mesh::EntityVector localEntities;
  m_bulkData.get_entities(stk::topology::ELEM_RANK, m_bulkData.mesh_meta_data().locally_owned_part(), localEntities);
  m_bulkData.batch_change_entity_parts(localEntities, {subdomainPart}, {});
}

std::string
M2NDecomposerNested::get_initial_subdomain_part_name(int subdomainId)
{
  std::ostringstream out;
  out << "initial_subdomain_" << subdomainId;
  return out.str();
}

}
}
}
