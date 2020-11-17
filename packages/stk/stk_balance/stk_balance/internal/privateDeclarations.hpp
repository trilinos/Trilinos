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

#ifndef BALANCE_PRIVATEDECLARATIONS_HPP
#define BALANCE_PRIVATEDECLARATIONS_HPP

#include <stk_util/util/ReportHandler.hpp>

#include <stk_mesh/base/GetEntities.hpp>
#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/CreateEdges.hpp>
#include <stk_mesh/base/CreateFaces.hpp>
#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/baseImpl/elementGraph/BulkDataIdMapper.hpp>
#include <stk_util/environment/Env.hpp>

#include <Teuchos_ParameterList.hpp>
#include <stk_balance/internal/StkMeshAdapterForZoltan2.hpp>
#include <stk_balance/internal/StkBalanceUtils.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_ColoringProblem.hpp>

#include <stk_balance/balanceUtils.hpp>

#include <zoltan.h>

#include <algorithm>

namespace stk
{
namespace balance
{
namespace internal
{

inline unsigned get_index(const int second_dim, const int third_dim, const int first_index, const int second_index, const int third_index)
{
    return first_index*third_dim*second_dim + second_index*third_dim + third_index;
}

void fillEntityCentroid(const stk::mesh::BulkData &stkMeshBulkData,  const stk::mesh::FieldBase* coord, stk::mesh::Entity entityOfConcern, double *elementCentroid);

void addBoxForFace(stk::mesh::BulkData &stkMeshBulkData, stk::mesh::Entity face, const double eps, SearchBoxIdentProcs &faceBoxes, const stk::mesh::FieldBase* coord);

void addGraphEdgesUsingBBSearch(stk::mesh::BulkData& stkMeshBulkData,
                                const BalanceSettings &balanceSettings,
                                std::vector<GraphEdge>& graphEdges,
                                const stk::mesh::Selector& searchSelector);

void calculateGeometricOrGraphBasedDecomp(const BalanceSettings& balanceSettings, const int num_procs_decomp, stk::mesh::EntityProcVec &decomp, stk::mesh::BulkData& stkMeshBulkData, const std::vector<stk::mesh::Selector>& selectors);

void rebalance(stk::mesh::BulkData& stkMeshBulkData, const stk::mesh::EntityProcVec& mockDecomposition);
void rebalance(stk::mesh::BulkData& stkMeshBulkData, const std::vector<unsigned>& mappings, const stk::mesh::EntityProcVec& mockDecomposition);
void rebalance(DecompositionChangeList &changeList);

void logMessage(MPI_Comm communicator, const std::string &message);

bool isElementPartOfSpider(const stk::mesh::BulkData & stkMeshBulkData,
                           const stk::mesh::Field<int> & beamConnectivityCountField,
                           stk::mesh::Entity element);

bool shouldOmitSpiderElement(const stk::mesh::BulkData & stkMeshBulkData,
                             const stk::balance::BalanceSettings & balanceSettings,
                             stk::mesh::Entity elem);

void fill_spider_connectivity_count_fields(stk::mesh::BulkData & bulk, const BalanceSettings & balanceSettings);
void keep_spiders_on_original_proc(stk::mesh::BulkData &bulk, const stk::balance::BalanceSettings & balanceSettings, DecompositionChangeList &changeList);
void fix_spider_elements(const BalanceSettings & balanceSettings, stk::mesh::BulkData & stkMeshBulkData);

void createZoltanParallelGraph(const BalanceSettings& balanceSettings, stk::mesh::BulkData& stkMeshBulkData,
                               const std::vector<stk::mesh::Selector>& selectors, const stk::mesh::impl::LocalIdMapper& localIds,
                               Zoltan2ParallelGraph& zoltan2Graph);

void add_connected_entities_of_rank(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Entity element, int newOwningProc, stk::mesh::EntityRank rank, std::vector<std::pair<stk::mesh::Entity, int> > &entityProcPairs);

unsigned get_local_id(const stk::mesh::impl::LocalIdMapper& localIds, stk::mesh::Entity entity);

bool is_geometric_method(const std::string& method);

stk::mesh::EntityVector get_entities_to_balance(stk::mesh::Selector selector, stk::mesh::EntityRank primaryRank, const stk::mesh::BulkData& bulkData);

void print_rebalance_metrics(const size_t num_global_entity_migrations, const size_t max_global_entity_migrations, stk::mesh::BulkData & stkMeshBulkData);

template <class ZoltanAdapter>
void print_zoltan_statistics(const ZoltanAdapter& stkMeshAdapter, Zoltan2::EvaluatePartition<ZoltanAdapter> &zoltanEvaluateParition, int parallel_rank)
{
    const double *kdd = NULL;
    int kddstride;
    stkMeshAdapter.getWeightsView(kdd, kddstride, 0);
    double localsum = 0.;
    for (size_t i = 0; i< stkMeshAdapter.getLocalNumIDs(); i++) localsum += kdd[i];

    sierra::Env::outputP0() << parallel_rank
            << " PRINTING METRICS nObj " << stkMeshAdapter.getLocalNumIDs()
            << " nwgts " << stkMeshAdapter.getNumWeightsPerID()
            << " sumwgt " << localsum << std::endl;
    zoltanEvaluateParition.printMetrics(sierra::Env::outputP0());
    //sierra::Env::outputP0() << parallel_rank
    //        << " PRINTING GRAPH METRICS" << std::endl;
    //zoltanEvaluateParition.printGraphMetrics(sierra::Env::outputP0());
}

template <class ZoltanAdapter>
void print_statistics(const ZoltanAdapter& stkMeshAdapter, MPI_Comm comm, int parallel_rank )
{
    //KDD Expect this interface to change -- I promise!!!
    auto teuchos_comm = Teuchos::rcp(new Teuchos::MpiComm<int>(comm));
    Teuchos::ParameterList pl;
    Zoltan2::EvaluatePartition<ZoltanAdapter> ep(&stkMeshAdapter, &pl, teuchos_comm, nullptr);
    print_zoltan_statistics(stkMeshAdapter, ep, parallel_rank);
}

template <class ZoltanAdapter>
void print_solution_statistics(const ZoltanAdapter& stkMeshAdapter, const Zoltan2::PartitioningSolution<ZoltanAdapter>& solution, MPI_Comm comm, int parallel_rank )
{
    //KDD Expect this interface to change -- I promise!!!
    auto teuchos_comm = Teuchos::rcp(new Teuchos::MpiComm<int>(comm));
    Teuchos::ParameterList pl;
    Zoltan2::EvaluatePartition<ZoltanAdapter> ep(&stkMeshAdapter, &pl, teuchos_comm, &solution);
    print_zoltan_statistics(stkMeshAdapter, ep, parallel_rank);
}

}
}
}
#endif /* PRIVATEDECLARATIONS_HPP_ */
