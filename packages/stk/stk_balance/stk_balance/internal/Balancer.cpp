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

#include "Balancer.hpp"
#include "stk_balance/balance.hpp"
#include "stk_balance/internal/privateDeclarations.hpp"
#include "stk_balance/internal/balanceCoincidentElements.hpp"
#include "stk_balance/internal/DetectAndFixMechanisms.hpp"
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <limits>

namespace stk {
namespace balance {

void verify_mesh_ids_are_not_too_large(const stk::mesh::BulkData& bulk)
{
  stk::mesh::EntityId maxAllowableId = std::numeric_limits<BalanceGlobalNumber>::max();
  stk::mesh::EntityId maxElemId = stk::mesh::get_max_id_on_local_proc(bulk, stk::topology::ELEM_RANK);

  stk::mesh::EntityId maxNodeId = stk::mesh::get_max_id_on_local_proc(bulk, stk::topology::NODE_RANK);
  stk::mesh::EntityId maxMeshId = std::max(maxElemId, maxNodeId);

  bool localMeshIdTooLarge = maxMeshId > maxAllowableId;
  bool globalMeshIdTooLarge = stk::is_true_on_any_proc(bulk.parallel(), localMeshIdTooLarge);

  ThrowRequireMsg(!globalMeshIdTooLarge, "Mesh has global-ids too large to represent in the index type that was configured for Trilinos ("<<sizeof(BalanceGlobalNumber)<<" bytes).");
}

bool loadBalance(const BalanceSettings& balanceSettings, stk::mesh::BulkData& stkMeshBulkData, unsigned numSubdomainsToCreate, const std::vector<stk::mesh::Selector>& selectors)
{
    internal::logMessage(stkMeshBulkData.parallel(), "Computing new decomposition");

    if (sizeof(BalanceGlobalNumber) < sizeof(stk::mesh::EntityId)) {
      verify_mesh_ids_are_not_too_large(stkMeshBulkData);
    }

    stk::mesh::EntityProcVec decomp;
    internal::calculateGeometricOrGraphBasedDecomp(balanceSettings, numSubdomainsToCreate, decomp, stkMeshBulkData, selectors);

    DecompositionChangeList changeList(stkMeshBulkData, decomp);
    balanceSettings.modifyDecomposition(changeList);

    internal::logMessage(stkMeshBulkData.parallel(), "Moving coincident elements to the same processor");
    keep_coincident_elements_together(stkMeshBulkData, changeList);

    if (balanceSettings.shouldFixSpiders()) {
        internal::logMessage(stkMeshBulkData.parallel(), "Preventing unnecessary movement of spider elements");
        internal::keep_spiders_on_original_proc(stkMeshBulkData, balanceSettings, changeList);
    }

    const size_t num_global_entity_migrations = changeList.get_num_global_entity_migrations();
    const size_t max_global_entity_migrations = changeList.get_max_global_entity_migrations();

    if (num_global_entity_migrations > 0)
    {
        internal::logMessage(stkMeshBulkData.parallel(), "Moving elements to new processors");
        internal::rebalance(changeList);

        if (balanceSettings.shouldFixMechanisms())
        {
            internal::logMessage(stkMeshBulkData.parallel(), "Fixing mechanisms found during decomposition");
            stk::balance::internal::detectAndFixMechanisms(balanceSettings, stkMeshBulkData);
        }

        if (balanceSettings.shouldFixSpiders())
        {
            internal::logMessage(stkMeshBulkData.parallel(), "Fixing spider elements");
            stk::balance::internal::fix_spider_elements(balanceSettings, stkMeshBulkData);
        }

        if (balanceSettings.shouldPrintMetrics())
            internal::print_rebalance_metrics(num_global_entity_migrations, max_global_entity_migrations, stkMeshBulkData);
    }

    internal::logMessage(stkMeshBulkData.parallel(), "Finished rebalance");


    return (num_global_entity_migrations > 0);
}

Balancer::Balancer(const BalanceSettings& settings)
  : m_settings(settings)
{ }

bool Balancer::balance(BalanceMesh& mesh) const
{
  std::vector<stk::mesh::Selector> selectors = {mesh.get_bulk().mesh_meta_data().locally_owned_part()};
  return balance(mesh, selectors);
}

bool Balancer::balance(BalanceMesh& mesh, const std::vector<stk::mesh::Selector>& selectors) const
{
  if( m_settings.getGraphOption() == BalanceSettings::LOAD_BALANCE ) {
    return loadBalance(m_settings, mesh.get_bulk(), mesh.get_bulk().parallel_size(), selectors);
  }
  return false;
}

} }

