#include "MxNutils.hpp"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>

namespace stk {
namespace balance {
namespace internal {

std::vector<unsigned> assign_target_subdomains_roundrobin_to_procs(unsigned num_procs_M, unsigned num_procs_N)
{
    std::vector<unsigned> targetSubdomainsToProc(num_procs_N, std::numeric_limits<unsigned>::max());
    for(unsigned i = 0; i < num_procs_N; i++)
        targetSubdomainsToProc[i] = i % num_procs_M;
    return targetSubdomainsToProc;
}

void fill_decomp(const int num_partitions, stk::mesh::BulkData& bulk, const stk::balance::BalanceSettings &graphSettings, stk::mesh::EntityProcVec &decomp)
{
    std::vector<stk::mesh::Selector> selectors = { bulk.mesh_meta_data().locally_owned_part() };
    stk::balance::internal::calculateGeometricOrGraphBasedDecomp(graphSettings, num_partitions, decomp, bulk, selectors);
}

stk::mesh::EntityProcVec get_element_decomp(const int num_partitions, stk::mesh::BulkData& bulk, const stk::balance::BalanceSettings &graphSettings)
{
    stk::mesh::EntityProcVec decomp;
    fill_decomp(num_partitions, bulk, graphSettings, decomp);
    return decomp;
}

}}}
