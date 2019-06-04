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

#ifndef GRAPHEDGEDATA_HPP_
#define GRAPHEDGEDATA_HPP_

#include "ElemElemGraphImpl.hpp"

namespace stk
{
namespace mesh
{

struct GraphEdge;

class GraphEdgesForElement
{
public:
    typedef std::vector<GraphEdge>::const_iterator const_iterator;

    const_iterator begin() const
    {
        return graphEdges.begin();
    }
    const_iterator end() const
    {
        return graphEdges.end();
    }

    size_t size() const
    {
        return graphEdges.size();
    }
    void reserve(size_t size)
    {
        graphEdges.reserve(size);
    }
    const GraphEdge & get_edge_at_index(impl::LocalId elem) const
    {
        return graphEdges[elem];
    }
    void emplace_back(const GraphEdge &graphEdge)
    {
        graphEdges.push_back(graphEdge);
    }
    void erase_at_index(impl::LocalId elem)
    {
        graphEdges.erase(graphEdges.begin() + elem);
    }
    void clear()
    {
        graphEdges.clear();
    }
private:
    std::vector<GraphEdge> graphEdges;
};

class Graph
{
public:
    void set_num_local_elements(size_t n);
    void add_new_element();
    size_t get_num_elements_in_graph() const;
    size_t get_num_edges() const;
    size_t get_num_edges_for_element(impl::LocalId elem) const;
    const GraphEdge & get_edge_for_element(impl::LocalId elem1, size_t index) const;
    const GraphEdgesForElement& get_edges_for_element(impl::LocalId elem) const;
    std::vector<GraphEdge> get_edges_for_element_side(impl::LocalId elem, int side) const;

    void reserve_edges(impl::LocalId localElemId, size_t numEdges);
    void add_edge(const GraphEdge &graphEdge);
    void delete_edge_from_graph(impl::LocalId local_elem_id, int offset);
    void delete_edge(const GraphEdge &graphEdge);
    void delete_all_edges(impl::LocalId elem);
    void clear();
    void delete_vertex(stk::mesh::impl::LocalId id)
    {
        m_graphEdges.erase(m_graphEdges.begin()+id);
    }

private:
    std::vector<GraphEdgesForElement> m_graphEdges;
    size_t m_numEdges = 0;
};

class ParallelInfoForGraphEdges
{
public:
    ParallelInfoForGraphEdges(int procRank) : m_procRank(procRank) {}
    impl::ParallelInfo& get_parallel_info_for_graph_edge(const GraphEdge& graphEdge);
    const impl::ParallelInfo& get_parallel_info_for_graph_edge(const GraphEdge& graphEdge) const;
    impl::ParallelGraphInfo &get_parallel_graph_info() { return m_parallel_graph_info; }
    const impl::ParallelGraphInfo &get_parallel_graph_info() const { return m_parallel_graph_info; }

    void insert_parallel_info_for_graph_edge(const GraphEdge& graphEdge, const impl::ParallelInfo& p_info);
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
