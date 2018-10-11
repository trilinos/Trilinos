#include <vector>
#include "balance.hpp"
#include "balanceUtils.hpp"               // for BalanceSettings, etc
#include "internal/privateDeclarations.hpp"  // for callZoltan1, etc
#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "stk_util/util/ReportHandler.hpp"  // for ThrowRequireMsg
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include "stk_mesh/base/MeshDiagnostics.hpp"

namespace stk
{
namespace balance
{

void rebalance_elements(stk::mesh::BulkData& bulkData, const stk::mesh::EntityProcVec &elementsToMove)
{
    stk::balance::DecompositionChangeList changeList(bulkData, elementsToMove);
    if(bulkData.parallel_rank()==0)
        std::cerr << "Fixing up mesh due to detection of coincident elements split across multiple processors.\n";
    stk::balance::internal::rebalance(changeList);
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
    internal::logMessage(bulkData.parallel(), "Fixing split coincident elements");
    stk::mesh::SplitCoincidentInfo splitCoincidentElements = stk::mesh::get_split_coincident_elements(bulkData);
    bool allOkThisProc = splitCoincidentElements.empty();
    bool allOkEverywhere = stk::is_true_on_all_procs(bulkData.parallel(), allOkThisProc);
    if(!allOkEverywhere)
        return rebalance_mesh_to_avoid_split_coincident_elements(bulkData, splitCoincidentElements);
    return stk::mesh::EntityIdProcMap();
}

}}
