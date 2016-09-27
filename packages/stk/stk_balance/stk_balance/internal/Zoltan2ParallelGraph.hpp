
#ifndef ZOLTAN2PARALLELGRAPH_HPP_
#define ZOLTAN2PARALLELGRAPH_HPP_

#include "Vertices.hpp"
#include <stk_mesh/baseImpl/elementGraph/BulkDataIdMapper.hpp>
#include "balanceTypes.hpp"
#include <balanceUtils.hpp>

namespace stk { namespace balance { class GraphEdge;}}

class Zoltan2ParallelGraph : public stk::balance::internal::Vertices
{
public:
    const std::vector<double> &get_edge_weights() const { return mEdgeWeights; }
    const std::vector<BalanceLocalNumber> &get_offsets() const { return mOffsets; }
    const std::vector<BalanceGlobalNumber> &get_adjacency() const { return mAdjacency; }
    size_t get_num_global_elements() const { return mNumGlobalElements; }
    void set_num_global_elements( size_t num ) { mNumGlobalElements = num; }

    void adjust_vertex_weights(const stk::balance::BalanceSettings& balanceSettings,
                               stk::mesh::BulkData& stkMeshBulkData,
                               const std::vector<stk::mesh::Selector>& selectors,
                               const stk::mesh::impl::LocalIdMapper& localIds);

    void createGraphEdgesUsingNodeConnectivity(stk::mesh::BulkData &stkMeshBulkData,
                                               const stk::balance::BalanceSettings &balanceSettings,
                                               size_t numElements,
                                               std::vector<stk::balance::GraphEdge> &graphEdges,
                                               const stk::mesh::impl::LocalIdMapper& localIds);
    void fillZoltan2AdapterDataFromStkMesh(stk::mesh::BulkData &stkMeshBulkData,
                                           const stk::balance::BalanceSettings &balanceSettings,
                                           std::vector<int>& adjacencyProcs,
                                           const stk::mesh::Selector& searchSelector,
                                           const stk::mesh::impl::LocalIdMapper& localIds);
private:
    void convertGraphEdgesToZoltanGraph(const stk::mesh::BulkData& stkMeshBulkData,
                                          const std::vector<stk::balance::GraphEdge> &graphEdges,
                                          const unsigned numElements,
                                          std::vector<int>& adjacencyProcs,
                                          const stk::mesh::impl::LocalIdMapper& localIds);

    std::vector<double> mEdgeWeights;
    std::vector<BalanceLocalNumber> mOffsets;
    std::vector<BalanceGlobalNumber> mAdjacency;
    size_t mNumGlobalElements = 0;
};



#endif /* ZOLTAN2PARALLELGRAPH_HPP_ */
