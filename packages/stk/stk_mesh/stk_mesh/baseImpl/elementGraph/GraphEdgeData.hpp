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

#ifndef GRAPHEDGEDATA_HPP_
#define GRAPHEDGEDATA_HPP_

#include "ElemElemGraphImpl.hpp"
#include <stk_util/util/PairIter.hpp>
#include <stk_util/util/ReportHandler.hpp>

namespace stk
{
namespace mesh
{

struct GraphEdge;

using GraphEdgesForElement = PairIter<const GraphEdge*>;

class Graph
{
public:
    void set_num_local_elements(size_t n);
    void add_new_element();
    size_t get_num_elements_in_graph() const;
    size_t get_num_edges() const;
    size_t get_num_edges_for_element(impl::LocalId elem) const;
    const GraphEdge & get_edge_for_element(impl::LocalId elem1, size_t index) const;
    GraphEdgesForElement get_edges_for_element(impl::LocalId elem) const;
    std::vector<GraphEdge> get_edges_for_element_side(impl::LocalId elem, int side) const;

    void add_sorted_edges(const std::vector<GraphEdge>& graphEdge);
    void replace_sorted_edges(std::vector<GraphEdge>& graphEdge);
    void delete_sorted_edges(const std::vector<GraphEdge>& graphEdge);
    void clear();

private:
    void set_offsets();
    std::vector<GraphEdge> m_graphEdges;
    std::vector<unsigned> m_elemOffsets;
};

class ParallelInfoForGraphEdges
{
public:
    ParallelInfoForGraphEdges(int procRank) : m_procRank(procRank) {}
    impl::ParallelInfo& get_parallel_info_for_graph_edge(const GraphEdge& graphEdge);
    const impl::ParallelInfo& get_parallel_info_for_graph_edge(const GraphEdge& graphEdge) const;
    impl::ParallelGraphInfo &get_parallel_graph_info() { return m_parallel_graph_info; }
    const impl::ParallelGraphInfo &get_parallel_graph_info() const { return m_parallel_graph_info; }

    bool insert_parallel_info_for_graph_edge(const GraphEdge& graphEdge, const impl::ParallelInfo& p_info);
    void erase_parallel_info_for_graph_edge(const GraphEdge& graphEdge);

    impl::LocalId convert_remote_global_id_to_negative_local_id(stk::mesh::EntityId remoteElementId) const;
    stk::mesh::EntityId convert_negative_local_id_to_remote_global_id(impl::LocalId remoteElementId) const;
    void clear();
private:
    impl::ParallelGraphInfo::const_iterator get_parallel_info_iterator_for_graph_edge(const GraphEdge& graphEdge) const;
    impl::ParallelGraphInfo m_parallel_graph_info;
    int m_procRank;
};

GraphEdge create_symmetric_edge(const GraphEdge& graphEdge);

}
}
#endif /* GRAPHEDGEDATA_HPP_ */
