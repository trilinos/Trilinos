#include "GraphEdgeData.hpp"
#include "ElemElemGraphImpl.hpp"
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_util/util/ReportHandler.hpp>

namespace stk
{
namespace mesh
{

void Graph::set_num_local_elements(size_t n)
{
    m_elemOffsets.resize(n+1);
}

void Graph::add_new_element()
{
    if (m_elemOffsets.empty()) {
      m_elemOffsets.assign(1, 0);
    }
    m_elemOffsets.push_back(m_graphEdges.size());
}

size_t Graph::get_num_elements_in_graph() const
{
    return m_elemOffsets.size() - 1;
}

size_t Graph::get_num_edges() const
{
    return m_graphEdges.size();
}

size_t Graph::get_num_edges_for_element(impl::LocalId elem) const
{
    return m_elemOffsets[elem+1] - m_elemOffsets[elem];
}

const GraphEdge & Graph::get_edge_for_element(impl::LocalId elem1, size_t index) const
{
    return m_graphEdges[m_elemOffsets[elem1]+index];
}

void fill_graph_edges_for_elem_side(const GraphEdgesForElement &graphEdgesForElement, int side, std::vector<GraphEdge>& edges)
{
    for(size_t i = 0; i < graphEdgesForElement.size(); ++i)
    {
        if(graphEdgesForElement[i].side1() == side)
        {
            edges.push_back(graphEdgesForElement[i]);
        }
    }
}

std::vector<GraphEdge> Graph::get_edges_for_element_side(impl::LocalId elem, int side) const
{
    std::vector<GraphEdge> edges;
    fill_graph_edges_for_elem_side(get_edges_for_element(elem), side, edges);
    return edges;
}

GraphEdgesForElement Graph::get_edges_for_element(impl::LocalId elem) const
{
    const unsigned begin = m_elemOffsets[elem];
    const unsigned end = m_elemOffsets[elem+1];
    return GraphEdgesForElement(&m_graphEdges[begin], &m_graphEdges[end]);
}

void Graph::set_offsets()
{
  const unsigned numOffsets = m_elemOffsets.size();
  m_elemOffsets.assign(std::max(1u, numOffsets), 0);

  impl::LocalId prevElem = impl::INVALID_LOCAL_ID;
  unsigned edgeCounter = 0;
  for(const GraphEdge& edge : m_graphEdges) {
    impl::LocalId elem1 = edge.elem1();
    if (elem1 != prevElem) {
      if (prevElem != impl::INVALID_LOCAL_ID) {
        m_elemOffsets[prevElem] = edgeCounter;
      }
      edgeCounter = 0;
      prevElem = elem1;
    }
    ++edgeCounter;
  }

  if (prevElem != impl::INVALID_LOCAL_ID) {
    m_elemOffsets[prevElem] = edgeCounter;
  }

  unsigned edgeOffset = 0;
  size_t numElems = m_elemOffsets.size()-1;
  for(size_t i=0; i<numElems; ++i) {
    unsigned count = m_elemOffsets[i];
    m_elemOffsets[i] = edgeOffset;
    edgeOffset += count;
  }
  m_elemOffsets.back() = edgeOffset;
}

using IterType = std::vector<GraphEdge>::iterator;

void Graph::add_sorted_edges(const std::vector<GraphEdge>& graphEdges)
{
  ThrowAssertMsg(stk::util::is_sorted_and_unique(graphEdges, GraphEdgeLessByElem1()),"Input vector 'graphEdges' is expected to be sorted-and-unique");
  if (!graphEdges.empty()) {
    stk::util::insert_keep_sorted(graphEdges, m_graphEdges, GraphEdgeLessByElem1());
    set_offsets();
  }
}

void Graph::replace_sorted_edges(std::vector<GraphEdge>& graphEdges)
{
  m_graphEdges.swap(graphEdges);
  set_offsets();
}

void Graph::delete_sorted_edges(const std::vector<GraphEdge>& edgesToDelete)
{
  for(const GraphEdge& edgeToDelete : edgesToDelete) {
    impl::LocalId elem1 = edgeToDelete.elem1();
    for(unsigned offset = m_elemOffsets[elem1]; offset < m_elemOffsets[elem1+1]; ++offset) {
      GraphEdge& thisEdge = m_graphEdges[offset];
      if (thisEdge == edgeToDelete) {
        thisEdge.vertex1 = impl::INVALID_LOCAL_ID;
      }
    }
  }

  if (!edgesToDelete.empty()) {
    const unsigned offset = m_elemOffsets[edgesToDelete[0].elem1()];
    m_graphEdges.erase(std::remove_if(m_graphEdges.begin()+offset, m_graphEdges.end(),
                                      [](const GraphEdge& edge)
                                      { return edge.vertex1 == impl::INVALID_LOCAL_ID; }),
                       m_graphEdges.end());
    set_offsets();
  }
}

void Graph::clear()
{
    m_graphEdges.clear();
    m_elemOffsets.clear();
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
    s << "    remote topology: " << parInfo.m_remote_element_topology << std::endl;
    return s.str();
}

bool ParallelInfoForGraphEdges::insert_parallel_info_for_graph_edge(const GraphEdge& graphEdge, const impl::ParallelInfo &parInfo)
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

    return inserted.second;
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
