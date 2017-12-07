#include <vector>

#include "balance.hpp"
#include "balanceUtils.hpp"               // for BalanceSettings, etc
#include "internal/privateDeclarations.hpp"  // for callZoltan1, etc
#include "internal/balanceCoincidentElements.hpp"
#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "stk_util/environment/ReportHandler.hpp"  // for ThrowRequireMsg
#include <stk_io/FillMesh.hpp>
#include <stk_io/WriteMesh.hpp>
#include "stk_io/StkIoUtils.hpp"
#include <stk_tools/mesh_clone/MeshClone.hpp>
#include "internal/LastStepFieldWriter.hpp"
#include "stk_balance/internal/TransientFieldTransferById.hpp"
#include "stk_balance/internal/DetectAndFixMechanisms.hpp"
#include "fixSplitCoincidentElements.hpp"

namespace stk
{
namespace balance
{

bool loadBalance(const BalanceSettings& balanceSettings, stk::mesh::BulkData& stkMeshBulkData, unsigned numSubdomainsToCreate, const std::vector<stk::mesh::Selector>& selectors)
{
    internal::logMessage(stkMeshBulkData.parallel(), "Starting rebalance.");

    stk::mesh::EntityProcVec decomp;
    internal::calculateGeometricOrGraphBasedDecomp(balanceSettings, numSubdomainsToCreate, decomp, stkMeshBulkData, selectors);

    DecompositionChangeList changeList(stkMeshBulkData, decomp);
    balanceSettings.modifyDecomposition(changeList);
    keep_coincident_elements_together(stkMeshBulkData, changeList);
    const size_t num_global_entity_migrations = changeList.get_num_global_entity_migrations();
    const size_t max_global_entity_migrations = changeList.get_max_global_entity_migrations();

    if (num_global_entity_migrations > 0)
    {
        internal::rebalance(changeList);
        if(balanceSettings.shouldFixMechanisms())
        {
            bool mechanismsFound = stk::balance::internal::detectAndFixMechanisms(balanceSettings, stkMeshBulkData);
            if(mechanismsFound)
                internal::logMessage(stkMeshBulkData.parallel(), "Mechanisms were found and fixed during decomposition\n");
        }
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

void run_static_stk_balance_with_settings(stk::io::StkMeshIoBroker &stkInput, stk::mesh::BulkData &inputBulk, const std::string& outputFilename, MPI_Comm comm, stk::balance::BalanceSettings& graphOptions)
{
    stk::mesh::MetaData metaB;
    stk::mesh::BulkData bulkB(metaB, comm);

    stk::mesh::BulkData *balancedBulk = nullptr;
    if(stk::io::get_transient_fields(inputBulk.mesh_meta_data()).empty())
    {
        balancedBulk = &inputBulk;
    }
    else
    {
        stk::tools::copy_mesh(inputBulk, inputBulk.mesh_meta_data().universal_part(), bulkB);
        balancedBulk = &bulkB;
    }

    stk::balance::balanceStkMesh(graphOptions, *balancedBulk);

    stk::io::StkMeshIoBroker stkOutput;
    stkOutput.set_bulk_data(*balancedBulk);
    stkOutput.set_attribute_field_ordering_stored_by_part_ordinal(stkInput.get_attribute_field_ordering_stored_by_part_ordinal());

    stk::balance::internal::TransientFieldTransferById transfer(stkInput, stkOutput);
    transfer.transfer_and_write_transient_fields(outputFilename);
}

void initial_decomp_and_balance(stk::mesh::BulkData &bulk,
                      stk::balance::BalanceSettings& graphOptions,
                      const std::string& exodusFilename,
                      const std::string& outputFilename)
{
    stk::io::StkMeshIoBroker stkInput;
    stkInput.property_add(Ioss::Property("DECOMPOSITION_METHOD", "LINEAR"));
    stk::io::fill_mesh_preexisting(stkInput, exodusFilename, bulk);
    make_mesh_consistent_with_parallel_mesh_rule1(bulk);
    run_static_stk_balance_with_settings(stkInput, bulk, outputFilename, bulk.parallel(), graphOptions);
}

void run_stk_balance_with_settings(const std::string& outputFilename, const std::string& exodusFilename, MPI_Comm comm, stk::balance::BalanceSettings& graphOptions)
{
    stk::mesh::MetaData meta;
    stk::mesh::BulkData bulk(meta, comm);
    initial_decomp_and_balance(bulk, graphOptions, exodusFilename, outputFilename);
}

void run_stk_rebalance(const std::string& outputDirectory, const std::string& exodusFilename, MPI_Comm comm)
{
    stk::balance::GraphCreationSettings graphOptions;
    std::string outputFilename = outputDirectory + "/" + exodusFilename;
    run_stk_balance_with_settings(outputFilename, exodusFilename, comm, graphOptions);
}

}
}
