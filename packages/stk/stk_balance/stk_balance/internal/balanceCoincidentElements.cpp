#include <vector>
#include "balance.hpp"
#include "balanceUtils.hpp"               // for BalanceSettings, etc
#include "internal/privateDeclarations.hpp"  // for callZoltan1, etc
#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "stk_util/util/ReportHandler.hpp"  // for ThrowRequireMsg
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>

namespace stk
{
namespace balance
{

stk::mesh::impl::LocalId get_min_coincident(stk::mesh::impl::LocalId elemLocalId,
                                            const std::vector<stk::mesh::impl::LocalId>& coincidents)
{
    stk::mesh::impl::LocalId minCoin = elemLocalId;
    for(stk::mesh::impl::LocalId coin : coincidents)
        if(coin < minCoin)
            minCoin = coin;
    return minCoin;
}

void modify_list_to_keep_on_same_proc(stk::mesh::Entity coinElem,
                                      stk::mesh::Entity minCoinElem,
                                      int minCoinDestination,
                                      DecompositionChangeList& changeList)
{
    bool minCoinElemIsMoving = changeList.has_entity(minCoinElem);
    if(!minCoinElemIsMoving)
        changeList.delete_entity(coinElem);
    else
        changeList.set_entity_destination(coinElem, minCoinDestination);
}

void set_all_coincidents_to_min_destination(const stk::mesh::ElemElemGraph& graph,
                                            const std::vector<stk::mesh::impl::LocalId>& coincidents,
                                            stk::mesh::Entity elem,
                                            stk::mesh::impl::LocalId minCoinLocalId,
                                            DecompositionChangeList& changeList)
{
    stk::mesh::Entity minCoinElem = graph.get_entity_from_local_id(minCoinLocalId);
    int minCoinDestination = changeList.get_entity_destination(minCoinElem);
    modify_list_to_keep_on_same_proc(elem, minCoinElem, minCoinDestination, changeList);
    for(stk::mesh::impl::LocalId coin : coincidents)
    {
        stk::mesh::Entity coinElem = graph.get_entity_from_local_id(coin);
        modify_list_to_keep_on_same_proc(coinElem, minCoinElem, minCoinDestination, changeList);
    }
}

void move_coincidents_to_same_proc(const stk::mesh::ElemElemGraph& graph,
                                   stk::mesh::Entity elem,
                                   DecompositionChangeList& changeList)
{
    stk::mesh::impl::LocalId elemLocalId = graph.get_local_element_id(elem);
    std::vector<stk::mesh::impl::LocalId> coincidents = graph.get_coincident_elements(elemLocalId);
    if (!coincidents.empty()) {
        stk::mesh::impl::LocalId minCoinLocalId = get_min_coincident(elemLocalId, coincidents);
        set_all_coincidents_to_min_destination(graph, coincidents, elem, minCoinLocalId, changeList);
    }
}

void keep_coincident_elements_together(stk::mesh::BulkData &bulk, DecompositionChangeList &changeList)
{
    bulk.initialize_face_adjacent_element_graph();
    const stk::mesh::ElemElemGraph &graph = bulk.get_face_adjacent_element_graph();
    stk::mesh::EntityProcVec entityProcs = changeList.get_all_partition_changes();
    for(const stk::mesh::EntityProc &entityProc : entityProcs)
        if(bulk.entity_rank(entityProc.first) == stk::topology::ELEM_RANK)
            move_coincidents_to_same_proc(graph, entityProc.first, changeList);
}

}}
