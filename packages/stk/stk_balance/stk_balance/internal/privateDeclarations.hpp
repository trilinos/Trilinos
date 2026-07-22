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
#include <stk_mesh/baseImpl/elementGraph/BulkDataIdMapper.hpp>
#include <stk_util/parallel/OutputStreams.hpp>

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

class ElementCountDiagnostic;
class TotalElementWeightDiagnostic;

namespace internal
{

inline unsigned get_index(const int numSelectors, const int numCriteria, const int entityIndex,
                          const int selectorIndex, const int criterionIndex)
{
  return entityIndex*numCriteria*numSelectors + selectorIndex*numCriteria + criterionIndex;
}

void fillEntityCentroid(const stk::mesh::BulkData &stkMeshBulkData,  const stk::mesh::FieldBase* coord, stk::mesh::Entity entityOfConcern, double *elementCentroid);

void addBoxForFace(stk::mesh::BulkData &stkMeshBulkData, stk::mesh::Entity face, const double eps, SearchBoxIdentProcs &faceBoxes, const stk::mesh::FieldBase* coord);

void addGraphEdgesUsingBBSearch(stk::mesh::BulkData& stkMeshBulkData,
                                const BalanceSettings &balanceSettings,
                                std::vector<GraphEdge>& graphEdges,
                                const stk::mesh::Selector& searchSelector);

void calculateGeometricOrGraphBasedDecomp(stk::mesh::BulkData & stkMeshBulkData,
                                          const std::vector<stk::mesh::Selector> & selectors,
                                          const stk::ParallelMachine & decompCommunicator,
                                          const int numSubdomainsToCreate,
                                          const BalanceSettings & balanceSettings,
                                          stk::mesh::EntityProcVec & decomp);

void rebalance(stk::mesh::BulkData& stkMeshBulkData, const stk::mesh::EntityProcVec& mockDecomposition);
void rebalance(stk::mesh::BulkData& stkMeshBulkData, const std::vector<unsigned>& mappings, const stk::mesh::EntityProcVec& mockDecomposition);
void rebalance(DecompositionChangeList &changeList);

void logMessage(MPI_Comm communicator, const std::string &message);

bool is_element_part_of_spider(const stk::mesh::BulkData & stkMeshBulkData,
                               const stk::mesh::Field<int> & beamConnectivityCountField,
                               stk::mesh::Entity element);

bool should_omit_spider_element(const stk::mesh::BulkData & stkMeshBulkData,
                                const stk::balance::BalanceSettings & balanceSettings,
                                stk::mesh::Entity elem);

void register_internal_fields_and_parts(stk::mesh::BulkData & bulk, const stk::balance::BalanceSettings & balanceSettings);
void fill_spider_connectivity_count_fields_and_parts(stk::mesh::BulkData & bulk, const BalanceSettings & balanceSettings);
void fill_output_subdomain_field(const stk::mesh::BulkData & bulk, const BalanceSettings & balanceSettings,
                                 stk::mesh::EntityProcVec & decomp);
void fix_spider_elements(const BalanceSettings & balanceSettings, stk::mesh::BulkData & stkMeshBulkData,
                         DecompositionChangeList & changeList);

void createZoltanParallelGraph(stk::mesh::BulkData & stkMeshBulkData,
                               const stk::mesh::Selector & selector,
                               const stk::ParallelMachine & decompCommunicator,
                               const BalanceSettings & balanceSettings,
                               Zoltan2ParallelGraph & zoltan2Graph);

void add_connected_entities_of_rank(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Entity element, int newOwningProc, stk::mesh::EntityRank rank, std::vector<std::pair<stk::mesh::Entity, int> > &entityProcPairs);

unsigned get_local_id(const stk::mesh::impl::LocalIdMapper& localIds, stk::mesh::Entity entity);

bool is_geometric_method(const std::string& method);

stk::mesh::EntityVector get_entities_to_balance(stk::mesh::Selector selector,
                                                stk::mesh::EntityRank primaryRank,
                                                const stk::mesh::BulkData& bulkData);

void print_rebalance_metrics(const size_t num_global_entity_migrations, const size_t max_global_entity_migrations, stk::mesh::BulkData & stkMeshBulkData);

template <class ZoltanAdapter>
void print_zoltan_statistics(const ZoltanAdapter& stkMeshAdapter, Zoltan2::EvaluatePartition<ZoltanAdapter> &zoltanEvaluateParition, int parallel_rank)
{
  const double *kdd = NULL;
  int kddstride;
  stkMeshAdapter.getWeightsView(kdd, kddstride, 0);
  double localsum = 0.;
  for (size_t i = 0; i< stkMeshAdapter.getLocalNumIDs(); i++) localsum += kdd[i];

  stk::outputP0() << parallel_rank
                          << " PRINTING METRICS nObj " << stkMeshAdapter.getLocalNumIDs()
                          << " nwgts " << stkMeshAdapter.getNumWeightsPerID()
                          << " sumwgt " << localsum << std::endl;
  zoltanEvaluateParition.printMetrics(stk::outputP0());
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

void compute_connectivity_weight(const stk::mesh::BulkData & bulk,  const BalanceSettings & balanceSettings);

void compute_balance_diagnostics(const stk::mesh::BulkData & bulk, const stk::balance::BalanceSettings & balanceSettings);
void compute_element_count_diagnostic(ElementCountDiagnostic & diag, const stk::mesh::BulkData & bulk, int rank);
void compute_total_element_weight_diagnostic(TotalElementWeightDiagnostic & diag, const stk::mesh::BulkData & bulk,
                                             const stk::balance::BalanceSettings & balanceSettings,
                                             const stk::mesh::Field<double> & weightField, int rank);


class EnableAura
{
public:
  EnableAura(stk::mesh::BulkData & bulk);
  ~EnableAura();

private:
  stk::mesh::BulkData & m_bulk;
  bool m_weTurnedOnAura;
};

} // namespace internal
} // namespace balance
} // namespace stk
#endif /* PRIVATEDECLARATIONS_HPP_ */
