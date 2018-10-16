#include "GraphEdgeData.hpp"
#include "ElemElemGraphImpl.hpp"
#include <stk_util/util/ReportHandler.hpp>

namespace stk
{
namespace mesh
{

void Graph::set_num_local_elements(size_t n)
{
    m_graphEdges.resize(n);
}

void Graph::add_new_element()
{
    m_graphEdges.push_back(GraphEdgesForElement());
}

size_t Graph::get_num_elements_in_graph() const
{
    return m_graphEdges.size();
}

size_t Graph::get_num_edges() const
{
    return m_numEdges;
}

size_t Graph::get_num_edges_for_element(impl::LocalId elem) const
{
    return m_graphEdges[elem].size();
}

const GraphEdge & Graph::get_edge_for_element(impl::LocalId elem1, size_t index) const
{
    return m_graphEdges[elem1].get_edge_at_index(index);
}

void fill_graph_edges_for_elem_side(const GraphEdgesForElement &graphEdgesForElement, int side, std::vector<GraphEdge>& edges)
{
    for(size_t i = 0; i < graphEdgesForElement.size(); ++i)
    {
        if(graphEdgesForElement.get_edge_at_index(i).side1() == side)
        {
            edges.push_back(graphEdgesForElement.get_edge_at_index(i));
        }
    }
}

std::vector<GraphEdge> Graph::get_edges_for_element_side(impl::LocalId elem, int side) const
{
    std::vector<GraphEdge> edges;
    fill_graph_edges_for_elem_side(m_graphEdges[elem], side, edges);
    return edges;
}

const GraphEdgesForElement& Graph::get_edges_for_element(impl::LocalId elem) const
{
    return m_graphEdges[elem];
}

void Graph::add_edge(const GraphEdge &graphEdge)
{
    m_graphEdges[graphEdge.elem1()].push_back(graphEdge);
    ++m_numEdges;
}

void Graph::delete_edge_from_graph(impl::LocalId elem, int offset)
{
    m_graphEdges[elem].erase_at_index(offset);
    --m_numEdges;
}

void Graph::delete_edge(const GraphEdge &graphEdge)
{
    const size_t numConnected = m_graphEdges[graphEdge.elem1()].size();
    for(size_t i=0; i<numConnected; ++i)
        if(m_graphEdges[graphEdge.elem1()].get_edge_at_index(i) == graphEdge)
            delete_edge_from_graph(graphEdge.elem1(), i);
}

void Graph::delete_all_edges(impl::LocalId elem)
{
    m_numEdges -= m_graphEdges[elem].size();
    m_graphEdges[elem].clear();
}

void Graph::clear()
{
    m_numEdges = 0;
    m_graphEdges.clear();
}

impl::ParallelInfo& ParallelInfoForGraphEdges::get_parallel_info_for_graph_edge(const GraphEdge& graphEdge)
{
    return const_cast<impl::ParallelInfo&>(get_parallel_info_iterator_for_graph_edge(graphEdge)->second);
}

const impl::ParallelInfo& ParallelInfoForGraphEdges::get_parallel_info_for_graph_edge(const GraphEdge& graphEdge) const
{
    return get_parallel_info_iterator_for_graph_edge(graphEdge)->second;
}

void ParallelInfoForGraphEdges::erase_parallel_info_for_graph_edge(const GraphEdge& graphEdge)
{
    m_parallel_graph_info.erase(graphEdge);
}

impl::ParallelGraphInfo::const_iterator ParallelInfoForGraphEdges::get_parallel_info_iterator_for_graph_edge(const GraphEdge& graphEdge) const
{
    impl::ParallelGraphInfo::const_iterator iter = m_parallel_graph_info.find(graphEdge);
    ThrowRequireMsg( iter != m_parallel_graph_info.end(), "ERROR: Proc " << m_procRank << " failed to find parallel graph info for edge "
                     << graphEdge << ".");
    return iter;
}

std::string get_par_info_description(const impl::ParallelInfo &parInfo)
{
    std::ostringstream s;
    s << "    other proc: " << parInfo.get_proc_rank_of_neighbor() << std::endl;
    s << "    permutation: " << parInfo.m_permutation << std::endl;
    s << "    remote topology: " << parInfo.m_remote_element_toplogy << std::endl;
    return s.str();
}

void ParallelInfoForGraphEdges::insert_parallel_info_for_graph_edge(const GraphEdge& graphEdge, const impl::ParallelInfo &parInfo)
{
    std::pair<impl::ParallelGraphInfo::iterator, bool> inserted = m_parallel_graph_info.emplace(graphEdge, parInfo);
    if (!inserted.second)
    {
        const impl::ParallelInfo &existingParInfo = inserted.first->second;
        if (existingParInfo != parInfo) {
            ThrowErrorMsg("Program error. local elem/remote elem pair"
                            << " (" << graphEdge.elem1() << "," << graphEdge.side1() << "/" << convert_negative_local_id_to_remote_global_id(graphEdge.elem2()) << "," << graphEdge.side2() << ")"
                            << " on procs (" << m_procRank << "," << parInfo.get_proc_rank_of_neighbor() << ")"
                            << " already exists in map. Please contact sierra-help@sandia.gov for support." << std::endl
                            << "existing par info " << std::endl
                            << get_par_info_description(existingParInfo)
                            << "new par info " << std::endl
                            << get_par_info_description(parInfo));
        }
    }
}

impl::LocalId ParallelInfoForGraphEdges::convert_remote_global_id_to_negative_local_id(stk::mesh::EntityId remoteElementId) const
{
    return -remoteElementId;
}

stk::mesh::EntityId ParallelInfoForGraphEdges::convert_negative_local_id_to_remote_global_id(impl::LocalId negativeLocalId) const
{
    return -negativeLocalId;
}

void ParallelInfoForGraphEdges::clear()
{
    m_parallel_graph_info.clear();
}

GraphEdge create_symmetric_edge(const GraphEdge& graphEdge)
{
    return GraphEdge(graphEdge.elem2(), graphEdge.side2(), graphEdge.elem1(), graphEdge.side1());
}

}
}
