#ifndef STK_BALANCE_MECHANISMS_HPP
#define STK_BALANCE_MECHANISMS_HPP

#include "stk_mesh/baseImpl/elementGraph/BulkDataIdMapper.hpp"
#include "stk_mesh/base/Types.hpp"

namespace stk { namespace balance { class BalanceSettings; } }
namespace stk { namespace mesh { class LocalIdMapper; } }
class Zoltan2ParallelGraph;

namespace stk { namespace balance { namespace internal {

bool detectAndFixMechanisms(const stk::balance::BalanceSettings& graphSettings, stk::mesh::BulkData &bulk);

void move_components(const Zoltan2ParallelGraph &zoltan2Graph,
                     const stk::mesh::impl::LocalIdMapper &localIds,
                     stk::mesh::BulkData& bulk,
                     const std::vector<stk::mesh::EntityVector>& elementsToMove,
                     const std::vector<int>& componentsToMove);

void fill_zoltan2_parallel_graph(stk::mesh::BulkData& bulk, const stk::balance::BalanceSettings& graphSettings, const stk::mesh::impl::LocalIdMapper& localIds, Zoltan2ParallelGraph &zoltan2Graph);

}}} // end of namespaces

#endif
