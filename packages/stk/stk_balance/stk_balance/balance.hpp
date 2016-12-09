#ifndef STKBALANCE_HPP
#define STKBALANCE_HPP

namespace stk { namespace balance { class BalanceSettings; } }
namespace stk { namespace mesh { class BulkData; } }

#include <stk_mesh/base/Selector.hpp>
#include <vector>

namespace stk
{
namespace balance
{

bool balanceStkMesh(const BalanceSettings& balanceSettings, stk::mesh::BulkData& stkMeshBulkData);
bool balanceStkMesh(const BalanceSettings& balanceSettings, stk::mesh::BulkData& stkMeshBulkData, const std::vector<stk::mesh::Selector>& selectors);
stk::mesh::EntityIdProcMap make_mesh_consistent_with_parallel_mesh_rule1(stk::mesh::BulkData& bulkData);

}
}
#endif
