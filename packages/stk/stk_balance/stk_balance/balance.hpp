#ifndef STKBALANCE_HPP
#define STKBALANCE_HPP

namespace stk { namespace balance { class BalanceSettings; } }
namespace stk { namespace mesh { class BulkData; } }

#include <stk_mesh/base/Selector.hpp>
#include <vector>
#include <string>
#include <mpi.h>

namespace stk
{
namespace balance
{

bool balanceStkMesh(const BalanceSettings& balanceSettings, stk::mesh::BulkData& stkMeshBulkData);
bool balanceStkMesh(const BalanceSettings& balanceSettings, stk::mesh::BulkData& stkMeshBulkData, const std::vector<stk::mesh::Selector>& selectors);
stk::mesh::EntityIdProcMap make_mesh_consistent_with_parallel_mesh_rule1(stk::mesh::BulkData& bulkData);
void run_stk_rebalance(const std::string& outputDirectory, const std::string& exodusFilename, MPI_Comm comm);
void run_stk_balance_with_settings(const std::string& outputDirectory, const std::string& exodusFilename, MPI_Comm comm, stk::balance::BalanceSettings& graphOptions);

}
}
#endif
