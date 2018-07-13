// Copyright (c) 2013, Sandia Corporation.
 // Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 // the U.S. Government retains certain rights in this software.
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
 //     * Neither the name of Sandia Corporation nor the names of its
 //       contributors may be used to endorse or promote products derived
 //       from this software without specific prior written permission.
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

#ifndef ZOLTAN2PARALLELGRAPH_HPP_
#define ZOLTAN2PARALLELGRAPH_HPP_

#include "Vertices.hpp"
#include <stk_mesh/baseImpl/elementGraph/BulkDataIdMapper.hpp>
#include "balanceTypes.hpp"
#include <stk_balance/balanceUtils.hpp>

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

    void adjust_weights_for_small_meshes();

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

    void setMechanismCheckFlag() { mCheckingMechanisms = true; }
    bool amCheckingForMechanisms() const { return mCheckingMechanisms; }

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
    bool mCheckingMechanisms = false;
};



#endif /* ZOLTAN2PARALLELGRAPH_HPP_ */
