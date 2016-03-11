#ifndef GRAPHEDGEDATA_HPP_
#define GRAPHEDGEDATA_HPP_

#include "ElemElemGraphImpl.hpp"

namespace stk
{
namespace mesh
{

class GraphEdge;

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
    const GraphEdge & get_edge_at_index(impl::LocalId elem) const
    {
        return graphEdges[elem];
    }
    void push_back(const GraphEdge &graphEdge)
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

    void add_edge(const GraphEdge &graphEdge);
    void delete_edge_from_graph(impl::LocalId local_elem_id, int offset);
    void delete_edge(const GraphEdge &graphEdge);
    void delete_all_edges(impl::LocalId elem);
    void clear();
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
