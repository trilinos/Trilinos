#include <cmath>
#include <vector>
#include <Zoltan2_MeshAdapter.hpp>      // for MeshEntityType, etc
#include "Vertices.hpp"
#include "ZoltanGeometricAdapter.hpp"
#include <stk_balance/internal/privateDeclarations.hpp>

namespace stk {
namespace balance {

template<typename ZoltanAdapter>
std::vector<unsigned> get_decomp_per_entity(size_t num_entities, const typename ZoltanAdapter::part_t* decomposition)
{
    std::vector<unsigned> results(num_entities,0);
    for(size_t i=0;i<num_entities;i++)
        results[i] = decomposition[i];
    return results;
}

void setParamsGeometricMethod(Teuchos::ParameterList& params, const BalanceSettings& balanceSettings, int numParts)
{
    params.set("debug_level", "basic_status");
    params.set("imbalance_tolerance", 1.10);
    params.set("num_global_parts", numParts);
    params.set("algorithm", balanceSettings.getDecompMethod());
//    params.set("partitioning_objective", "balance_object_weight"); // RCB1
//    params.set("partitioning_objective", "multicriteria_minimize_total_weight"); // RCB2
//    params.set("partitioning_objective", "multicriteria_minimize_maximum_weight"); // RCB3
//    params.set("partitioning_objective", "multicriteria_balance_total_maximum"); // RCB4
    //    params.set("mj_parts", "1,4,4"); // multijagged special for x = one element thick
//    params.set("rectilinear", "1"); // This causes app test failures
    Teuchos::ParameterList &zparams = params.sublist("zoltan_parameters", false);
    zparams.set("RCB_RECOMPUTE_BOX", "1");
    zparams.set("REMAP", "0");
    if (balanceSettings.getDecompMethod() != "rcb" && balanceSettings.isIncrementalRebalance())
    {
        zparams.set("REMAP", "1");
    }
    zparams.set("debug_level", "0");
}


template<typename ZoltanAdapter>
std::vector<unsigned> get_solution_from_zoltan(ZoltanAdapter& adapter, Teuchos::ParameterList& params, MPI_Comm comm)
{
    Zoltan2::PartitioningProblem<ZoltanAdapter> problem(&adapter, &params, comm);
    internal::print_statistics(adapter, comm, stk::parallel_machine_rank(comm));
    std::srand(stk::parallel_machine_rank(comm)); // KHP: Temporary until an API is added to Zoltan2 for random seeds.
    problem.solve();
    internal::print_solution_statistics(adapter, problem.getSolution(), comm, stk::parallel_machine_rank(comm));
    std::vector<unsigned> proc_decomp = get_decomp_per_entity<ZoltanAdapter>(adapter.getLocalNumOf(Zoltan2::MESH_REGION), problem.getSolution().getPartListView());
    return proc_decomp;
}

std::vector<unsigned> get_decomposition(const stk::balance::internal::GeometricVertices& vertexInfo, const BalanceSettings& balanceSettings, int numParts, MPI_Comm comm)
{
    ZoltanGeometricAdapter adapter(vertexInfo);
    Teuchos::ParameterList params("test params");
    setParamsGeometricMethod(params, balanceSettings, numParts);
    return get_solution_from_zoltan(adapter, params, comm);
}

}
}
