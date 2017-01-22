
#include <vector>

#include "balance.hpp"
#include "balanceUtils.hpp"               // for BalanceSettings, etc
#include "internal/privateDeclarations.hpp"  // for callZoltan1, etc
#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "stk_util/environment/ReportHandler.hpp"  // for ThrowRequireMsg
#include "stk_mesh/base/MeshDiagnostics.hpp"
#include "stk_util/parallel/ParallelReduceBool.hpp"

namespace stk
{
namespace balance
{

bool loadBalance(const BalanceSettings& balanceSettings, stk::mesh::BulkData& stkMeshBulkData, unsigned numSubdomainsToCreate, const std::vector<stk::mesh::Selector>& selectors)
{
    internal::logMessage(stkMeshBulkData.parallel(), "Starting rebalance.");

    stk::mesh::EntityProcVec decomp;
    internal::callZoltan2(balanceSettings, numSubdomainsToCreate, decomp, stkMeshBulkData, selectors);

    DecompositionChangeList changeList(stkMeshBulkData, decomp);
    balanceSettings.modifyDecomposition(changeList);
    const size_t num_global_entity_migrations = changeList.get_num_global_entity_migrations();
    const size_t max_global_entity_migrations = changeList.get_max_global_entity_migrations();

    if (num_global_entity_migrations > 0)
    {
        internal::rebalance(changeList);
        if(balanceSettings.shouldPrintMetrics())
            internal::print_rebalance_metrics(num_global_entity_migrations, max_global_entity_migrations, stkMeshBulkData);
    }

    internal::logMessage(stkMeshBulkData.parallel(), "Finished rebalance.");

    return (num_global_entity_migrations > 0);
}

bool balanceStkMesh(const BalanceSettings& balanceSettings, stk::mesh::BulkData& stkMeshBulkData)
{
    std::vector<stk::mesh::Selector> selectors = {stkMeshBulkData.mesh_meta_data().locally_owned_part()};
    return balanceStkMesh(balanceSettings, stkMeshBulkData, selectors);
}

bool balanceStkMesh(const BalanceSettings& balanceSettings, stk::mesh::BulkData& stkMeshBulkData, const std::vector<stk::mesh::Selector>& selectors)
{
    switch(balanceSettings.getGraphOption())
    {
        case BalanceSettings::LOADBALANCE:
            return loadBalance(balanceSettings, stkMeshBulkData, stkMeshBulkData.parallel_size(), selectors);
            break;
        case BalanceSettings::COLORING:
            ThrowRequireMsg(false, "Coloring not implemented yet.");
            break;
    }
    return false;
}

stk::mesh::EntityProcVec get_rebalance_by_moving_split_coincident_elements(stk::mesh::BulkData& bulkData, const stk::mesh::SplitCoincidentInfo &splitCoincidentElements)
{
    stk::mesh::EntityProcVec elementsToMigrate;
    elementsToMigrate.reserve(splitCoincidentElements.size());
    for(const stk::mesh::SplitCoincidentInfo::value_type& item : splitCoincidentElements)
    {
        int min_proc = bulkData.parallel_size();
        stk::mesh::Entity element = bulkData.get_entity(stk::topology::ELEM_RANK,item.first);
        for(size_t i=0;i<item.second.size();++i)
            min_proc=std::min(min_proc, item.second[i].second);
        if(min_proc<bulkData.parallel_rank())
            elementsToMigrate.push_back(stk::mesh::EntityProc(element, min_proc));
    }
    return elementsToMigrate;
}

void rebalance_elements(stk::mesh::BulkData& bulkData, const stk::mesh::EntityProcVec &elementsToMove)
{
    stk::balance::DecompositionChangeList changeList(bulkData, elementsToMove);
    if(bulkData.parallel_rank()==0)
        std::cerr << "Fixing up mesh due to detection of violation of parallel mesh rule #1.\n";
    stk::balance::internal::rebalance(changeList);
}

stk::mesh::EntityIdProcMap rebalance_mesh_to_avoid_split_coincident_elements(stk::mesh::BulkData& bulkData, const stk::mesh::SplitCoincidentInfo &splitCoincidentElements)
{
    stk::mesh::EntityProcVec elementsToMove = get_rebalance_by_moving_split_coincident_elements(bulkData, splitCoincidentElements);
    stk::mesh::EntityIdProcMap entityIdProcMap;
    for(size_t i=0;i<elementsToMove.size();++i)
        entityIdProcMap[bulkData.identifier(elementsToMove[i].first)] = elementsToMove[i].second;
    rebalance_elements(bulkData, elementsToMove);
    return entityIdProcMap;
}

stk::mesh::EntityIdProcMap make_mesh_consistent_with_parallel_mesh_rule1(stk::mesh::BulkData& bulkData)
{
    stk::mesh::SplitCoincidentInfo splitCoincidentElements = stk::mesh::get_split_coincident_elements(bulkData);
    bool allOkThisProc = splitCoincidentElements.empty();
    bool allOkEverywhere = stk::is_true_on_all_procs(bulkData.parallel(), allOkThisProc);
    if(!allOkEverywhere)
        return rebalance_mesh_to_avoid_split_coincident_elements(bulkData, splitCoincidentElements);
    return stk::mesh::EntityIdProcMap();
}

}
}
