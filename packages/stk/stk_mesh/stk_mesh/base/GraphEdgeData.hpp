#ifndef GRAPHEDGEDATA_HPP_
#define GRAPHEDGEDATA_HPP_

#include "ElemElemGraphImpl.hpp"

namespace stk
{
namespace mesh
{

class GraphEdge;

class Graph
{
public:
    void set_num_local_elements(size_t n);
    void add_new_element();
    size_t get_num_elements_in_graph() const;
    size_t get_num_edges() const;
    size_t get_num_edges_for_element(impl::LocalId elem) const;
    const GraphEdge & get_edge_for_element(impl::LocalId elem1, size_t index) const;
    std::vector<GraphEdge> get_edges_for_element_side(impl::LocalId elem, int side) const;

    void add_edge(const GraphEdge &graphEdge);
    void delete_edge_from_graph(impl::LocalId local_elem_id, int offset);
    void delete_edge(const GraphEdge &graphEdge);
    void delete_all_connections(impl::LocalId elem);
    void change_elem2_id_for_edge(const stk::mesh::GraphEdge& edgeToChange, impl::LocalId newElem2Id);
    void clear();
private:
    std::vector<std::vector<GraphEdge> > m_graphEdges;
    size_t m_numEdges = 0;
};

class ParallelInfoForGraphEdges
{
public:
    ParallelInfoForGraphEdges(int procRank) : m_procRank(procRank) {}
    impl::parallel_info& get_parallel_info_for_graph_edge(const GraphEdge& graphEdge);
    const impl::parallel_info& get_parallel_info_for_graph_edge(const GraphEdge& graphEdge) const;
    impl::ParallelGraphInfo &get_parallel_graph_info() { return m_parallel_graph_info; }

    void insert_parallel_info_for_graph_edge(const GraphEdge& graphEdge, const impl::parallel_info&);
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
