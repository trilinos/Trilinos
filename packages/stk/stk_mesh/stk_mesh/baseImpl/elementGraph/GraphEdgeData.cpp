#include "GraphEdgeData.hpp"
#include "ElemElemGraphImpl.hpp"
#include "stk_mesh/baseImpl/elementGraph/GraphTypes.hpp"
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_util/util/ReportHandler.hpp>

namespace stk
{
namespace mesh
{

void Graph::set_num_local_elements(size_t n)
{
    m_elemOffsets.resize(n, IndexRange(m_graphEdges.size()+1, m_graphEdges.size()+1));
}

void Graph::add_new_element()
{
    m_elemOffsets.push_back({m_graphEdges.size()+1, m_graphEdges.size()+1});
}


size_t Graph::get_num_elements_in_graph() const
{
    return m_elemOffsets.size();
}

size_t Graph::get_num_edges() const
{
    return m_graphEdges.size() - m_numUnusedEntries;    
}

size_t Graph::get_num_edges_for_element(impl::LocalId elem) const
{
    auto& indices = m_elemOffsets[elem];
    return indices.second - indices.first;
}

const GraphEdge & Graph::get_edge_for_element(impl::LocalId elem1, size_t index) const
{
    STK_ThrowAssertMsg(get_num_edges_for_element(elem1) != 0, "Cannot retrieve graph edge for element that has no faces");
    STK_ThrowAssertMsg(get_num_edges_for_element(elem1) > index, "index out of range");

    return m_graphEdges[m_elemOffsets[elem1].first+index];
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
    if (m_graphEdges.data() != nullptr) {
      const unsigned beginOffset = m_elemOffsets[elem].first;
      const unsigned endOffset   = m_elemOffsets[elem].second;

      const GraphEdge* beginEdge = m_graphEdges.data() + beginOffset;
      const GraphEdge* endEdge   = m_graphEdges.data() + endOffset;
      return GraphEdgesForElement(beginEdge, endEdge);
    }

    return GraphEdgesForElement(nullptr, nullptr);
}


void Graph::set_offsets()
{
  if (m_graphEdges.size() == 0)
  {
    return;
  }

  impl::LocalId currElem = m_graphEdges[0].elem1();
  unsigned startIdx = 0;
  for (unsigned i=0; i < m_graphEdges.size(); ++i)
  {
    impl::LocalId nextElem = m_graphEdges[i].elem1();
    if (nextElem != currElem)
    {
      STK_ThrowAssertMsg(currElem >= 0 && size_t(currElem) <= m_elemOffsets.size(), "element out of range");
      m_elemOffsets[currElem] = IndexRange(startIdx, i);
      for (impl::LocalId elem=currElem+1; elem < nextElem; elem++)
      {
        m_elemOffsets[elem] = IndexRange(0, 0);
      }

      currElem = nextElem;
      startIdx = i;
    }
  }

  m_elemOffsets[currElem] = IndexRange(startIdx, m_graphEdges.size());
}


using IterType = std::vector<GraphEdge>::iterator;

void Graph::add_sorted_edges(const std::vector<GraphEdge>& graphEdges)
{
  STK_ThrowAssertMsg(stk::util::is_sorted_and_unique(graphEdges, GraphEdgeLessByElem1()),"Input vector 'graphEdges' is expected to be sorted-and-unique");

  for (auto& edge : graphEdges)
  {
    insert_edge(edge);
  }
}


void Graph::insert_edge(const GraphEdge& graphEdge)
{
  auto elem1 = graphEdge.elem1();
  auto& indices = m_elemOffsets[elem1];

  if (check_for_edge(graphEdge))
  {
    return;
  }

  if (m_graphEdges.size() > 0 && double(m_numUnusedEntries) / m_graphEdges.size() > m_compressionThreshold)
  {
    compress_graph();
  }

  if (get_num_edges_for_element(elem1) == 0)
  {
      m_graphEdges.push_back(graphEdge);
      indices.first  = m_graphEdges.size()-1;
      indices.second = m_graphEdges.size();
  } else if (indices.second >= m_graphEdges.size())
  {
    m_graphEdges.emplace_back();
    insert_edge_into_sorted_range_or_next_entry(indices, graphEdge);
  } else if (is_valid(m_graphEdges[indices.second]))
  {
    move_edges_to_end(elem1);

    m_graphEdges.emplace_back();
    insert_edge_into_sorted_range_or_next_entry(indices, graphEdge);
  } else if (!is_valid(m_graphEdges[indices.second]))
  {
    insert_edge_into_sorted_range_or_next_entry(indices, graphEdge);
    m_numUnusedEntries--;
  } else
  {
    throw std::runtime_error("unreachable case");
  }
}

void Graph::insert_edge_into_sorted_range_or_next_entry(IndexRange& indices, const GraphEdge& graphEdge)
{
    unsigned idxToInsert = find_sorted_insertion_index(indices, graphEdge);

    for (unsigned i=indices.second; i > idxToInsert; i--)
    {
      m_graphEdges[i] = m_graphEdges[i-1];
    }

    m_graphEdges[idxToInsert] = graphEdge;
    indices.second++;
}


unsigned Graph::find_sorted_insertion_index(IndexRange indices, const GraphEdge& graphEdge)
{
    GraphEdgeLessByElem2Only isLess;
    for (unsigned i=indices.first; i < indices.second; ++i)
    {
      if (isLess(graphEdge, m_graphEdges[i]))
      {
          return i;
      }
    }

    return indices.second;
}

void Graph::replace_sorted_edges(std::vector<GraphEdge>& graphEdges)
{
  STK_ThrowAssertMsg(stk::util::is_sorted_and_unique(graphEdges, GraphEdgeLessByElem1()),"Input vector 'graphEdges' is expected to be sorted-and-unique");

  m_graphEdges.swap(graphEdges);
  set_offsets();
  m_numUnusedEntries = 0;
}


void Graph::delete_sorted_edges(const std::vector<GraphEdge>& edgesToDelete)
{
  STK_ThrowAssertMsg(std::is_sorted(edgesToDelete.begin(), edgesToDelete.end(), GraphEdgeLessByElem1()),
                "Input vector is expected to be sorted");

  int startIdx = 0;
  while (size_t(startIdx) != edgesToDelete.size())
  {
    int endIdx = get_end_of_element_range_for_sorted_edges(edgesToDelete, startIdx);
    for (int idx=endIdx; idx >= startIdx; idx--)
    {
      delete_edge(edgesToDelete[idx]);
    }

    startIdx = endIdx + 1;
  }
}

unsigned Graph::get_end_of_element_range_for_sorted_edges(const std::vector<GraphEdge>& edges, unsigned startIdx)
{
    unsigned currElement = edges[startIdx].elem1();
    unsigned endIdx = startIdx;
    while (endIdx < edges.size() && edges[endIdx].elem1() == currElement)
    {
      endIdx++;
    }
    endIdx--;

    return endIdx;
}

void Graph::delete_edge(const GraphEdge& edgeToDelete)
{
  impl::LocalId elem1 = edgeToDelete.elem1();
  auto& indices = m_elemOffsets[elem1];
  for(unsigned offset = indices.first; offset < indices.second; ++offset) {
    if (m_graphEdges[offset] == edgeToDelete) 
    {
      for (unsigned i=offset; i < indices.second-1; ++i)
      {
        m_graphEdges[i] = m_graphEdges[i+1];
      }
      indices.second--;
      m_graphEdges[indices.second] = GraphEdge();
      m_numUnusedEntries++;
      break;
    }
  }
}

void Graph::clear()
{
    m_graphEdges.clear();
    m_elemOffsets.clear();
    m_numUnusedEntries = 0;
}

void Graph::print(std::ostream& os)
{
  for(impl::LocalId i=0; i<(impl::LocalId)get_num_elements_in_graph(); ++i) {
    GraphEdgesForElement graphEdges = get_edges_for_element(i);
    os << "elem "<<i<<":\n   ";
    for(const GraphEdge& graphEdge : graphEdges) {
      os << graphEdge << "; "<<std::endl;
    }
  }
  os << std::endl;
}

void Graph::move_edges_to_end(impl::LocalId elem)
{
  auto& indices = m_elemOffsets[elem];
  size_t newStartIdx = m_graphEdges.size();
  for (unsigned i=indices.first; i < indices.second; ++i)
  {
    m_graphEdges.push_back(m_graphEdges[i]);
    m_graphEdges[i] = GraphEdge();
    m_numUnusedEntries++;
  }

  m_elemOffsets[elem] = IndexRange(newStartIdx, m_graphEdges.size());
}

void Graph::compress_graph()
{
  if (m_graphEdges.size() == 0 || m_graphEdges.size() == m_numUnusedEntries)
    return;

  impl::LocalId prevElement = 0;
  unsigned offset = 0;
  for (unsigned i=0; i < m_graphEdges.size(); ++i)
  {
    if (is_valid(m_graphEdges[i]))
    {
      prevElement = m_graphEdges[i].elem1();
      break;
    } else
    {
      offset++;
    }
  }

  {
    auto& indices = m_elemOffsets[prevElement];
    indices.first  -= offset;
    indices.second -= offset;
  }

  for (unsigned idx=offset; idx < m_graphEdges.size(); ++idx)
  {
    if (is_valid(m_graphEdges[idx]))
    {
      m_graphEdges[idx - offset] = m_graphEdges[idx];
      
      impl::LocalId currElement = m_graphEdges[idx].elem1();
      if (currElement != prevElement)
      {
        auto& indices = m_elemOffsets[currElement];
        if (indices.first != indices.second)
        {
          indices.first  -= offset;
          indices.second -= offset;
        } 
        prevElement = currElement;
      }

    } else
    {
      offset++;
    }
  }

  STK_ThrowRequireMsg(is_valid(m_graphEdges[m_graphEdges.size() - offset - 1]), "The count of unused edges is incorrect");
  m_graphEdges.resize(m_graphEdges.size() - offset);
  m_numUnusedEntries = 0;
}


bool Graph::check_for_edge(const GraphEdge& edge)
{
  auto& indices = m_elemOffsets[edge.elem1()];
  for (unsigned i=indices.first; i < indices.second; ++i)
    if (m_graphEdges[i] == edge)
    {
      return true;
    }

  return false;
}

impl::ParallelInfo& ParallelInfoForGraphEdges::get_parallel_info_for_graph_edge(const GraphEdge& graphEdge)
{
  return const_cast<impl::ParallelInfo&>(get_parallel_info_iterator_for_graph_edge(graphEdge)->second);
}

const impl::ParallelInfo& ParallelInfoForGraphEdges::get_parallel_info_for_graph_edge(const GraphEdge& graphEdge) const
{
  return get_parallel_info_iterator_for_graph_edge(graphEdge)->second;
}

void ParallelInfoForGraphEdges::erase_edges(const std::vector<GraphEdge>& edges)
{
  for(const GraphEdge& edge : edges) {
    auto iter = get_parallel_info_iterator_for_graph_edge(edge);
    if (iter != m_parallel_graph_info.end()) {
      iter->second.set_proc_rank(-1);
    }
  }
  m_parallel_graph_info.erase(std::remove_if(m_parallel_graph_info.begin(), m_parallel_graph_info.end(),
                                             [&](const std::pair<GraphEdge,impl::ParallelInfo>& info){ return info.second.get_proc_rank_of_neighbor() == -1; }),
                              m_parallel_graph_info.end());
}

void ParallelInfoForGraphEdges::erase_parallel_info_for_graph_edge(const GraphEdge& graphEdge)
{
  auto iter = get_parallel_info_iterator_for_graph_edge(graphEdge);
  if (iter != m_parallel_graph_info.end()) {
    m_parallel_graph_info.erase(iter);
  }
}

impl::ParallelGraphInfo::const_iterator ParallelInfoForGraphEdges::get_parallel_info_iterator_for_graph_edge(const GraphEdge& graphEdge, bool throwIfNotFound) const
{
    impl::ParallelGraphInfo::const_iterator iter = std::lower_bound(m_parallel_graph_info.begin(), m_parallel_graph_info.end(), graphEdge, GraphEdgeLessByElem2());
    if (throwIfNotFound) {
      STK_ThrowRequireMsg( iter != m_parallel_graph_info.end() && iter->first == graphEdge,
        "ERROR: Proc " << m_procRank << " failed to find parallel graph info for edge "
                     << graphEdge << ".");
    }
    return iter;
}

impl::ParallelGraphInfo::iterator ParallelInfoForGraphEdges::get_parallel_info_iterator_for_graph_edge(const GraphEdge& graphEdge, bool throwIfNotFound)
{
    impl::ParallelGraphInfo::iterator iter = std::lower_bound(m_parallel_graph_info.begin(), m_parallel_graph_info.end(), graphEdge, GraphEdgeLessByElem2());
    if (throwIfNotFound) {
      STK_ThrowRequireMsg( iter != m_parallel_graph_info.end() && iter->first == graphEdge, "ERROR: Proc " << m_procRank << " failed to find parallel graph info for edge "
                     << graphEdge << ".");
    }
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

void ParallelInfoForGraphEdges::insert_sorted_edges(const impl::ParallelGraphInfo& newParallelEdges)
{
  m_parallel_graph_info.reserve(m_parallel_graph_info.size() + newParallelEdges.size());
  stk::util::insert_keep_sorted(newParallelEdges, m_parallel_graph_info, GraphEdgeLessByElem2());
}

bool ParallelInfoForGraphEdges::find_parallel_info_for_graph_edge(const GraphEdge& graphEdge) const
{
    impl::ParallelGraphInfo::const_iterator iter = std::lower_bound(m_parallel_graph_info.begin(), m_parallel_graph_info.end(), graphEdge, GraphEdgeLessByElem2());
    return iter != m_parallel_graph_info.end() && iter->first == graphEdge;
}

bool ParallelInfoForGraphEdges::insert_parallel_info_for_graph_edge(const GraphEdge& graphEdge, const impl::ParallelInfo &parInfo)
{
    impl::ParallelGraphInfo::iterator iter = std::lower_bound(m_parallel_graph_info.begin(), m_parallel_graph_info.end(), graphEdge, GraphEdgeLessByElem2());
    if (iter == m_parallel_graph_info.end() || iter->first != graphEdge)
    {
        m_parallel_graph_info.insert(iter, std::make_pair(graphEdge, parInfo));
    }
    else {
        if (iter->second != parInfo) {
            STK_ThrowErrorMsg("Program error. local elem/remote elem pair"
                            << " (" << graphEdge.elem1() << "," << graphEdge.side1() << "/" << convert_negative_local_id_to_remote_global_id(graphEdge.elem2()) << "," << graphEdge.side2() << ")"
                            << " on procs (" << m_procRank << "," << parInfo.get_proc_rank_of_neighbor() << ")"
                            << " already exists in map. Please contact sierra-help@sandia.gov for support." << std::endl
                            << "existing par info " << std::endl
                            << get_par_info_description(iter->second)
                            << "new par info " << std::endl
                            << get_par_info_description(parInfo));
        }
    }

    return true;
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
