#ifndef STK_ELEM_SPARSE_GRAPH_HPP
#define STK_ELEM_SPARSE_GRAPH_HPP

#include <vector>
#include <map>
#include "ElemElemGraphImpl.hpp"

namespace stk { namespace mesh { class ParallelInfoForGraphEdges; } }
namespace stk { namespace mesh { namespace impl { class IdMapper; } } }

namespace stk
{
namespace mesh
{
namespace impl
{

class SparseGraph
{
public:
    typedef std::map<LocalId, std::vector<GraphEdge>> ElemToCoincidentEdgesMap;
    typedef ElemToCoincidentEdgesMap::value_type value_type;

    size_t get_num_elements_in_graph() const { return m_graphEdges.size(); }
    const std::vector<GraphEdge> & get_edges_for_element(impl::LocalId elem) const;
    void add_edge(const GraphEdge &graphEdge);
    void delete_edge(const GraphEdge &graphEdge);
    void delete_all_edges(impl::LocalId elem);
    void clear()
    {
        m_graphEdges.clear();
    }

    ElemToCoincidentEdgesMap::const_iterator begin() const { return m_graphEdges.begin(); }
    ElemToCoincidentEdgesMap::const_iterator end() const { return m_graphEdges.end(); }

private:
    ElemToCoincidentEdgesMap m_graphEdges;
    static const std::vector<GraphEdge> emptyVector;

    void remove_coincident_edge(const GraphEdge& graphEdge, std::vector<GraphEdge>& graphEdgesForElem, size_t i);
};


}}} // end namespaces stk mesh

#endif
