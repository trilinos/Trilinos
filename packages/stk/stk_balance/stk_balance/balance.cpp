// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <vector>

#include "mesh/BalanceMesh.hpp"

#include "balance.hpp"
#include "balanceUtils.hpp"               // for BalanceSettings, etc
#include "fixSplitCoincidentElements.hpp"
#include "stk_balance/internal/Balancer.hpp"
#include "internal/DetectAndFixMechanisms.hpp"
#include "internal/LastStepFieldWriter.hpp"
#include "internal/balanceCoincidentElements.hpp"
#include "internal/privateDeclarations.hpp"  // for callZoltan1, etc
#include "internal/NodeBalancer.hpp"
#include "stk_io/StkIoUtils.hpp"
#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "stk_mesh/base/Comm.hpp"
#include "stk_tools/transfer_utils/TransientFieldTransferById.hpp"
#include "stk_util/diag/StringUtil.hpp"
#include "stk_util/parallel/ParallelReduce.hpp"
#include "stk_util/util/ReportHandler.hpp"  // for ThrowRequireMsg
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

bool balanceStkMesh(const BalanceSettings& balanceSettings, stk::mesh::BulkData& stkMeshBulkData, const std::vector<stk::mesh::Selector>& selectors)
{
  const Balancer balancer(balanceSettings);
  BalanceMesh mesh(stkMeshBulkData);
  return balancer.balance(mesh, selectors);
}

bool balanceStkMesh(const BalanceSettings& balanceSettings, stk::mesh::BulkData& stkMeshBulkData)
{
  const Balancer balancer(balanceSettings);
  BalanceMesh mesh(stkMeshBulkData);
  return balancer.balance(mesh);
}

bool balanceStkMeshNodes(const BalanceSettings& balanceSettings, stk::mesh::BulkData& stkMeshBulkData)
{
    if ((balanceSettings.getGraphOption() == BalanceSettings::LOAD_BALANCE) && balanceSettings.useNodeBalancer())
    {
      internal::NodeBalancer nodeBalancer(stkMeshBulkData);
      return nodeBalancer.balance_node_entities(balanceSettings.getNodeBalancerTargetLoadBalance(),
                                                balanceSettings.getNodeBalancerMaxIterations());
    }
    return false;
}

bool colorStkMesh(const BalanceSettings& colorSettings, stk::mesh::BulkData& stkMeshBulkData)
{
    if (colorSettings.getGraphOption() == BalanceSettings::COLOR_MESH )
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

}
}
