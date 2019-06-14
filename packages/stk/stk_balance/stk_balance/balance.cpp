#include <vector>

#include "balance.hpp"
#include "balanceUtils.hpp"               // for BalanceSettings, etc
#include "fixSplitCoincidentElements.hpp"
#include "internal/LastStepFieldWriter.hpp"
#include "internal/balanceCoincidentElements.hpp"
#include "internal/privateDeclarations.hpp"  // for callZoltan1, etc
#include "stk_balance/internal/DetectAndFixMechanisms.hpp"
#include "stk_tools/transfer_utils/TransientFieldTransferById.hpp"
#include "stk_io/StkIoUtils.hpp"
#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "stk_mesh/base/Comm.hpp"
#include "stk_util/diag/StringUtil.hpp"
#include "stk_util/parallel/ParallelReduce.hpp"
#include "stk_util/util/ReportHandler.hpp"  // for ThrowRequireMsg
#include <stk_balance/internal/balanceDefaults.hpp>
#include <stk_balance/search_tolerance_algs/SecondShortestEdgeFaceSearchTolerance.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_tools/mesh_clone/MeshClone.hpp>

namespace stk
{
namespace balance
{

namespace {
    bool check_if_mesh_has_coloring(const stk::mesh::MetaData& meta)
    {
        stk::mesh::PartVector coloringParts;
        fill_coloring_parts(meta, coloringParts);
        return !coloringParts.empty();
    }

    void move_entities_to_coloring_part(stk::mesh::BulkData& bulk,
                                        const stk::mesh::EntityRank rank,
                                        const stk::mesh::Part& rootTopologyPart,
                                        const stk::mesh::impl::LocalIdMapper& localIds,
                                        int* coloredGraphVertices)
    {

        stk::mesh::MetaData& meta = bulk.mesh_meta_data();
        stk::mesh::EntityVector entities;
        stk::mesh::Selector entitySelector = meta.locally_owned_part() & rootTopologyPart;
        stk::mesh::get_selected_entities(entitySelector, bulk.buckets(rank), entities);
        std::vector<stk::mesh::PartVector> addParts;
        std::vector<stk::mesh::PartVector> removeParts;
        for(stk::mesh::Entity entity : entities)
        {
            if(localIds.does_entity_have_local_id(entity))
            {
                unsigned localId = localIds.entity_to_local(entity);
                int color = coloredGraphVertices[localId];
                std::string partName = construct_coloring_part_name(color, rootTopologyPart);
                stk::mesh::Part* colorPart = meta.get_part(partName);
                ThrowRequireMsg(nullptr != colorPart, "Color Part for " << bulk.entity_key(entity) << " cannot be null!");
                addParts.push_back({colorPart});
                removeParts.push_back({});
            }
        }
        bulk.batch_change_entity_parts(entities, addParts, removeParts);
    }

    stk::mesh::FieldBase* get_coloring_field(const stk::mesh::MetaData& meta, const stk::mesh::Part& rootTopologyPart)
    {
        stk::mesh::FieldBase* colorField = nullptr;
        colorField = meta.get_field(stk::topology::ELEMENT_RANK, rootTopologyPart.topology().name() + "coloring");
        return colorField;
    }

    void update_color_fields(stk::mesh::BulkData& bulk,
                             const stk::mesh::EntityRank rank,
                             const stk::mesh::Part& rootTopologyPart,
                             const stk::mesh::impl::LocalIdMapper& localIds,
                             int* coloredGraphVertices)
    {
        stk::mesh::MetaData& meta = bulk.mesh_meta_data();
        stk::mesh::FieldBase* colorField = get_coloring_field(meta, rootTopologyPart);
        ThrowRequireMsg(colorField != nullptr, "Root topology part not supported, created after I/O for topology " << rootTopologyPart.topology().name());

        stk::mesh::EntityVector entities;
        stk::mesh::get_selected_entities(rootTopologyPart, bulk.buckets(rank), entities);
        for(stk::mesh::Entity entity : entities)
        {
            if(localIds.does_entity_have_local_id(entity))
            {
                unsigned localId = localIds.entity_to_local(entity);
                int* colorData = static_cast<int*>(stk::mesh::field_data(*colorField, entity));
                *colorData = coloredGraphVertices[localId];
            }
        }
    }

    void fill_coloring_set(int* coloredGraphVertices, size_t numColoredGraphVertices, std::set<int>& colors)
    {
        for(size_t vertex = 0; vertex < numColoredGraphVertices; ++vertex)
        {
            colors.insert(coloredGraphVertices[vertex]);
        }
    }

    int get_max_color(stk::mesh::MetaData& meta, const std::set<int>& colors)
    {
        int localMaxColor = -1;
        if(colors.size() > 0)
        {
            localMaxColor = *colors.rbegin();
        }
        int globalMaxColor;
        stk::all_reduce_max<int>(meta.mesh_bulk_data().parallel(), &localMaxColor, &globalMaxColor, 1);
        return globalMaxColor;
    }

    void create_coloring_parts(stk::mesh::MetaData& meta, const stk::mesh::Part& rootTopologyPart, const std::set<int>& colors)
    {
        int globalMaxColor = get_max_color(meta, colors);
        for(int color = 1; color <= globalMaxColor; ++color)
        {
            std::string partName = construct_coloring_part_name(color, rootTopologyPart);
            meta.declare_part(partName);
        }
    }
}

std::string get_selected_part_name(const stk::mesh::Part& part)
{
  std::string partName = "";
  partName += part.name() + "_";
  return partName;
}

std::string construct_coloring_part_name(const int color, const stk::mesh::Part& part)
{
    std::ostringstream oss;
    oss << get_coloring_part_base_name() << "_" << get_selected_part_name(part) << color;
    return oss.str();
}

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

bool colorMesh(const BalanceSettings& balanceSettings, stk::mesh::BulkData& bulk, const stk::mesh::PartVector& parts)
{
    ThrowRequireMsg(balanceSettings.getGraphOption() == BalanceSettings::COLOR_MESH ||
                    balanceSettings.getGraphOption() == BalanceSettings::COLOR_MESH_BY_TOPOLOGY ||
                    balanceSettings.getGraphOption() == BalanceSettings::COLOR_MESH_AND_OUTPUT_COLOR_FIELDS,
                    "colorMesh must be called with COLOR_MESH or COLOR_MESH_BY_TOPOLOGY Setting");

    internal::logMessage(bulk.parallel(), "Start Coloring Mesh");

    stk::mesh::MetaData& meta = bulk.mesh_meta_data();
    ThrowRequireMsg(!check_if_mesh_has_coloring(meta), "Mesh has already been colored!");

    const stk::mesh::EntityRank rank = stk::topology::ELEM_RANK;

    Teuchos::ParameterList params("stk_balance coloring mesh");

    int totalNumColors = 0;
    for (stk::mesh::Part* part : parts)
    {
      internal::logMessage(bulk.parallel(), "Coloring Part: " + part->name());
      std::vector<int> adjacencyProcs;
      stk::mesh::Selector selector = *part;
      stk::mesh::impl::LocalIdMapper localIds(bulk, rank, selector);

      Zoltan2ParallelGraph zoltan2Graph;
      zoltan2Graph.fillZoltan2AdapterDataFromStkMesh(bulk,
                                                     balanceSettings,
                                                     adjacencyProcs,
                                                     selector,
                                                     localIds);

      std::vector<size_t> counts;
      stk::mesh::comm_mesh_counts(bulk, counts, &selector);
      zoltan2Graph.set_num_global_elements(counts[rank]);

      zoltan2Graph.set_spatial_dim(meta.spatial_dimension());
      StkMeshZoltanAdapter stkMeshAdapter(zoltan2Graph);

      Zoltan2::ColoringProblem<StkMeshZoltanAdapter> problem(&stkMeshAdapter, &params);

      std::srand(bulk.parallel_rank());
      internal::logMessage(bulk.parallel(), "Solving graph for colors");
      problem.solve();

      Zoltan2::ColoringSolution<StkMeshZoltanAdapter> *soln = problem.getSolution();

      size_t numColoredGraphVertices = soln->getColorsSize();
      int* coloredGraphVertices = soln->getColors();

      std::set<int> colors;
      fill_coloring_set(coloredGraphVertices, numColoredGraphVertices, colors);

      create_coloring_parts(meta, *part, colors);

      internal::logMessage(bulk.parallel(), "Moving entities to coloring part");
      move_entities_to_coloring_part(bulk, rank, *part, localIds, coloredGraphVertices);

      if (balanceSettings.getGraphOption() == BalanceSettings::COLOR_MESH_AND_OUTPUT_COLOR_FIELDS)
      {
          update_color_fields(bulk, rank, *part, localIds, coloredGraphVertices);
      }
      totalNumColors += colors.size();
    }
    internal::logMessage(bulk.parallel(), "Finish coloring");
    return totalNumColors > 0;
}

bool balanceStkMesh(const BalanceSettings& balanceSettings, stk::mesh::BulkData& stkMeshBulkData)
{
    std::vector<stk::mesh::Selector> selectors = {stkMeshBulkData.mesh_meta_data().locally_owned_part()};
    return balanceStkMesh(balanceSettings, stkMeshBulkData, selectors);
}

bool balanceStkMesh(const BalanceSettings& balanceSettings, stk::mesh::BulkData& stkMeshBulkData, const std::vector<stk::mesh::Selector>& selectors)
{
    if( balanceSettings.getGraphOption() == BalanceSettings::LOAD_BALANCE )
    {
        return loadBalance(balanceSettings, stkMeshBulkData, stkMeshBulkData.parallel_size(), selectors);
    }
    return false;
}

bool colorStkMesh(const BalanceSettings& colorSettings, stk::mesh::BulkData& stkMeshBulkData)
{
    if(colorSettings.getGraphOption() == BalanceSettings::COLOR_MESH )
    {
      return colorMesh(colorSettings, stkMeshBulkData, {&(stkMeshBulkData.mesh_meta_data().locally_owned_part())});
    }
    else if (colorSettings.getGraphOption() == BalanceSettings::COLOR_MESH_BY_TOPOLOGY || colorSettings.getGraphOption() == BalanceSettings::COLOR_MESH_AND_OUTPUT_COLOR_FIELDS)
    {
      stk::mesh::PartVector rootTopologyParts = get_root_topology_parts_for_rank(stkMeshBulkData, stk::topology::ELEMENT_RANK);
      return colorMesh(colorSettings, stkMeshBulkData, rootTopologyParts);
    }
    return false;
}

void fill_coloring_parts(const stk::mesh::MetaData& meta, stk::mesh::PartVector& coloringParts)
{
    coloringParts.clear();
    const stk::mesh::PartVector& parts = meta.get_parts();
    const std::string& coloringPartBaseName = get_coloring_part_base_name();
    const unsigned length = coloringPartBaseName.length();
    for (stk::mesh::Part* part : parts)
    {
        std::string partSubName = part->name().substr(0, length);
        if (!sierra::case_strcmp(partSubName, coloringPartBaseName))
        {
            coloringParts.push_back(part);
        }
    }
}

void fill_coloring_parts_with_topology(const stk::mesh::MetaData& meta, const stk::topology topo, stk::mesh::PartVector& coloringParts)
{
    coloringParts.clear();
    const stk::mesh::PartVector& parts = meta.get_parts();
    const std::string& coloringPartBaseName = get_coloring_part_base_name();
    const std::string coloringTopoPartBaseName = coloringPartBaseName + "_" + meta.get_topology_root_part(topo).name();
    const unsigned length = coloringTopoPartBaseName.length();
    for (stk::mesh::Part* part : parts)
    {
        std::string partSubName = part->name().substr(0, length);
        if (!sierra::case_strcmp(partSubName, coloringTopoPartBaseName))
        {
            coloringParts.push_back(part);
        }
    }
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

    stk::transfer_utils::TransientFieldTransferById transfer(stkInput, stkOutput);
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
