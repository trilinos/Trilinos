#ifndef BALANCE_PRIVATEDECLARATIONS_HPP
#define BALANCE_PRIVATEDECLARATIONS_HPP

#include <stk_util/environment/ReportHandler.hpp>
// stk search (order matters!)
#include <stk_search/IdentProc.hpp>
#include <stk_search/BoundingBox.hpp>

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

#include <Teuchos_ParameterList.hpp>
#include <stk_balance/internal/StkMeshAdapterForZoltan2.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_ColoringProblem.hpp>

#include <stk_balance/balanceUtils.hpp>

#include <zoltan.h>

#include <algorithm>

typedef stk::search::IdentProc<stk::mesh::EntityId, int> StkMeshIdent;
typedef std::vector<std::pair<StkMeshIdent, StkMeshIdent> > StkSearchResults;
typedef std::pair<stk::search::Box<float>, StkMeshIdent> BoxWithStkId;
typedef std::vector< BoxWithStkId > BoxVectorWithStkId;


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

void addBoxForFace(stk::mesh::BulkData &stkMeshBulkData, stk::mesh::Entity face, const double eps, BoxVectorWithStkId &faceBoxes, const stk::mesh::FieldBase* coord);

void fillFaceBoxesWithIds(stk::mesh::BulkData &stkMeshBulkData, const double eps, const stk::mesh::FieldBase* coord, BoxVectorWithStkId &faceBoxes, const stk::mesh::Selector& searchSelector);

void createGraphEdgesUsingBBSearch(stk::mesh::BulkData& stkMeshBulkData, const BalanceSettings &balanceSettings,
                                   std::vector<GraphEdge>& graphEdges, const stk::mesh::Selector& searchSelector);

void callZoltan2(const BalanceSettings& balanceSettings, const int num_procs_decomp, stk::mesh::EntityProcVec &decomp, stk::mesh::BulkData& stkMeshBulkData, const std::vector<stk::mesh::Selector>& selectors);

void rebalance(stk::mesh::BulkData& stkMeshBulkData, const stk::mesh::EntityProcVec& mockDecomposition);
void rebalance(stk::mesh::BulkData& stkMeshBulkData, const std::vector<unsigned>& mappings, const stk::mesh::EntityProcVec& mockDecomposition);
void rebalance(DecompositionChangeList &changeList);

void logMessage(MPI_Comm communicator, const std::string &message);

void fill_zoltan2_graph(const BalanceSettings& balanceSettings, stk::mesh::BulkData& stkMeshBulkData, Zoltan2ParallelGraph& zoltan2Graph,
                        const stk::mesh::Selector& searchSelector, const stk::mesh::impl::LocalIdMapper& localIds);

void add_connected_entities_of_rank(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Entity element, int newOwningProc, stk::mesh::EntityRank rank, std::vector<std::pair<stk::mesh::Entity, int> > &entityProcPairs);

void fill_list_of_entities_to_send_for_aura_like_ghosting(stk::mesh::BulkData& bulkData, stk::mesh::EntityProcVec &entitiesToGhost);

stk::mesh::Ghosting * create_custom_ghosting(stk::mesh::BulkData & stkMeshBulkData);
void destroy_custom_ghosting(stk::mesh::BulkData & stkMeshBulkData, stk::mesh::Ghosting * customAura);

unsigned get_local_id(const stk::mesh::impl::LocalIdMapper& localIds, stk::mesh::Entity entity);

bool is_geometric_method(const std::string method);

const stk::mesh::FieldBase * get_coordinate_field(const stk::mesh::MetaData& meta_data,
                                                  const std::string& coordinateFieldName);

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
    if (parallel_rank == 0)
    {
        std::cout << parallel_rank
                << " PRINTING METRICS nObj " << stkMeshAdapter.getLocalNumIDs()
                << " nwgts " << stkMeshAdapter.getNumWeightsPerID()
                << " sumwgt " << localsum << std::endl;
        zoltanEvaluateParition.printMetrics(std::cout);
        //std::cout << parallel_rank
        //        << " PRINTING GRAPH METRICS" << std::endl;
        //zoltanEvaluateParition.printGraphMetrics(std::cout);
    }
}

template <class ZoltanAdapter>
void print_statistics(const ZoltanAdapter& stkMeshAdapter, int parallel_rank )
{
    //KDD Expect this interface to change -- I promise!!!
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm(); // This should be our stkMeshBulkData.parallel()
    Teuchos::ParameterList pl;
    Zoltan2::EvaluatePartition<ZoltanAdapter> ep(&stkMeshAdapter, &pl, comm, nullptr);
    print_zoltan_statistics(stkMeshAdapter, ep, parallel_rank);
}

template <class ZoltanAdapter>
void print_solution_statistics(const ZoltanAdapter& stkMeshAdapter, const Zoltan2::PartitioningSolution<ZoltanAdapter>& solution, int parallel_rank )
{
    //KDD Expect this interface to change -- I promise!!!
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm(); // This should be our stkMeshBulkData.parallel()
    Teuchos::ParameterList pl;
    Zoltan2::EvaluatePartition<ZoltanAdapter> ep(&stkMeshAdapter, &pl, comm, &solution);
    print_zoltan_statistics(stkMeshAdapter, ep, parallel_rank);
}

}
}
}
#endif /* PRIVATEDECLARATIONS_HPP_ */
