#include <vector>

#include "balance.hpp"
#include "balanceUtils.hpp"               // for BalanceSettings, etc
#include "internal/privateDeclarations.hpp"  // for callZoltan1, etc
#include "internal/balanceCoincidentElements.hpp"
#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "stk_util/util/ReportHandler.hpp"  // for ThrowRequireMsg
#include <stk_io/FillMesh.hpp>
#include <stk_io/WriteMesh.hpp>
#include "stk_io/StkIoUtils.hpp"
#include <stk_tools/mesh_clone/MeshClone.hpp>
#include "internal/LastStepFieldWriter.hpp"
#include "stk_balance/internal/TransientFieldTransferById.hpp"
#include "stk_balance/internal/DetectAndFixMechanisms.hpp"
#include <stk_balance/internal/balanceDefaults.hpp>
#include <stk_balance/search_tolerance_algs/SecondShortestEdgeFaceSearchTolerance.hpp>
#include "fixSplitCoincidentElements.hpp"

namespace stk
{
namespace balance
{

bool loadBalance(const BalanceSettings& balanceSettings, stk::mesh::BulkData& stkMeshBulkData, unsigned numSubdomainsToCreate, const std::vector<stk::mesh::Selector>& selectors)
{
    internal::logMessage(stkMeshBulkData.parallel(), "Computing new decomposition");

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

bool balanceStkMesh(const BalanceSettings& balanceSettings, stk::mesh::BulkData& stkMeshBulkData)
{
    std::vector<stk::mesh::Selector> selectors = {stkMeshBulkData.mesh_meta_data().locally_owned_part()};
    return balanceStkMesh(balanceSettings, stkMeshBulkData, selectors);
}

bool balanceStkMesh(const BalanceSettings& balanceSettings, stk::mesh::BulkData& stkMeshBulkData, const std::vector<stk::mesh::Selector>& selectors)
{
  if( balanceSettings.getGraphOption() == BalanceSettings::LOADBALANCE )
  {
    return loadBalance(balanceSettings, stkMeshBulkData, stkMeshBulkData.parallel_size(), selectors);
  }
  ThrowRequireMsg(balanceSettings.getGraphOption() != BalanceSettings::COLORING, "Coloring not implemented yet.");
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
        internal::logMessage(inputBulk.parallel(), "Copying input mesh to handle transient fields");
        stk::tools::copy_mesh(inputBulk, inputBulk.mesh_meta_data().universal_part(), bulkB);
        balancedBulk = &bulkB;
    }

    stk::balance::balanceStkMesh(graphOptions, *balancedBulk);

    stk::io::StkMeshIoBroker stkOutput;
    stkOutput.set_bulk_data(*balancedBulk);
    stkOutput.set_attribute_field_ordering_stored_by_part_ordinal(stkInput.get_attribute_field_ordering_stored_by_part_ordinal());

    stk::balance::internal::TransientFieldTransferById transfer(stkInput, stkOutput);
    transfer.transfer_and_write_transient_fields(outputFilename);

    internal::logMessage(inputBulk.parallel(), "Finished writing output mesh");
}

void register_internal_fields(stk::mesh::BulkData& bulkData, stk::balance::BalanceSettings& balanceSettings)
{
    if (balanceSettings.shouldFixSpiders()) {
        stk::mesh::MetaData& meta = bulkData.mesh_meta_data();
        stk::mesh::Field<int> & field = meta.declare_field<stk::mesh::Field<int>>(stk::topology::NODE_RANK,
                                                                                  balanceSettings.getSpiderConnectivityCountFieldName());
        const int initValue = 0;
        stk::mesh::put_field_on_mesh(field, meta.universal_part(), &initValue);
    }
}

void read_mesh_with_auto_decomp(stk::io::StkMeshIoBroker & stkIo,
                                const std::string& meshSpec,
                                stk::mesh::BulkData& bulkData,
                                stk::balance::BalanceSettings & balanceSettings)
{
    stkIo.set_bulk_data(bulkData);
    stkIo.add_mesh_database(meshSpec, stk::io::READ_MESH);
    stkIo.create_input_mesh();
    stkIo.add_all_mesh_fields_as_input_fields();

    register_internal_fields(bulkData, balanceSettings);

    stkIo.populate_bulk_data();

    if(stkIo.check_integer_size_requirements() == 8) {
        bulkData.set_large_ids_flag(true);
    }
}

void initial_decomp_and_balance(stk::mesh::BulkData &bulk,
                                stk::balance::BalanceSettings& graphOptions,
                                const std::string& exodusFilename,
                                const std::string& outputFilename)
{
    stk::io::StkMeshIoBroker stkInput;
    stkInput.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RIB"));

    internal::logMessage(bulk.parallel(), "Reading mesh and performing initial decomposition");
    read_mesh_with_auto_decomp(stkInput, exodusFilename, bulk, graphOptions);

    make_mesh_consistent_with_parallel_mesh_rule1(bulk);
    run_static_stk_balance_with_settings(stkInput, bulk, outputFilename, bulk.parallel(), graphOptions);
}

void run_stk_balance_with_settings(const std::string& outputFilename, const std::string& exodusFilename, MPI_Comm comm, stk::balance::BalanceSettings& graphOptions)
{
    const std::string trimmedInputName = (exodusFilename.substr(0,2) == "./") ? exodusFilename.substr(2) : exodusFilename;
    const std::string trimmedOutputName = (outputFilename.substr(0,2) == "./") ? outputFilename.substr(2) : outputFilename;
    const bool isSerial = (stk::parallel_machine_size(comm) == 1);
    const bool inputEqualsOutput = (trimmedOutputName == trimmedInputName);
    ThrowRequireMsg(!(isSerial && inputEqualsOutput),
                    "Running on 1 MPI rank and input-file ("<<exodusFilename
                     <<") == output-file, doing nothing. Specify outputDirectory if you "
                     <<"wish to copy the input-file to an output-file of the same name.");

    stk::mesh::MetaData meta;
    stk::mesh::BulkData bulk(meta, comm);
    initial_decomp_and_balance(bulk, graphOptions, exodusFilename, outputFilename);
}

void run_stk_rebalance(const std::string& outputDirectory, const std::string& exodusFilename, stk::balance::AppTypeDefaults appType, MPI_Comm comm)
{
    stk::balance::GraphCreationSettings graphOptions;

    if (appType == stk::balance::SD_DEFAULTS)
    {
        graphOptions.setShouldFixSpiders(true);
    }
    else if (appType == stk::balance::SM_DEFAULTS)
    {
        graphOptions.setEdgeWeightForSearch(3.0);
        graphOptions.setVertexWeightMultiplierForVertexInSearch(10.0);
        graphOptions.setToleranceFunctionForFaceSearch(std::make_shared<stk::balance::SecondShortestEdgeFaceSearchTolerance>());
    }

    std::string outputFilename = outputDirectory + "/" + exodusFilename;
    run_stk_balance_with_settings(outputFilename, exodusFilename, comm, graphOptions);
}

}
}
