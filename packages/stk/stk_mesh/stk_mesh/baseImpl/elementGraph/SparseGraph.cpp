#include "ElemElemGraphImpl.hpp"
#include "ElemElemGraph.hpp"
#include "SparseGraph.hpp"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>


namespace stk
{
namespace mesh
{
namespace impl
{

const std::vector<GraphEdge> SparseGraph::emptyVector;

const std::vector<GraphEdge> & SparseGraph::get_edges_for_element(LocalId elemId) const
{
    auto iter = m_graphEdges.find(elemId);
    if(iter != m_graphEdges.end())
        return iter->second;

    return emptyVector;
}

void SparseGraph::remove_coincident_edge(const GraphEdge& graphEdge, std::vector<GraphEdge>& graphEdgesForElem, size_t edgeIndex)
{
    graphEdgesForElem.erase(graphEdgesForElem.begin() + edgeIndex);
    if(graphEdgesForElem.empty())
        m_graphEdges.erase(graphEdge.elem1());
}

void SparseGraph::delete_edge(const GraphEdge &graphEdge)
{
    std::vector<GraphEdge> &graphEdgesForElem = m_graphEdges[graphEdge.elem1()];
    for(size_t edgeIndex=0; edgeIndex<graphEdgesForElem.size(); ++edgeIndex)
    {
        if(graphEdgesForElem[edgeIndex] == graphEdge)
        {
            remove_coincident_edge(graphEdge, graphEdgesForElem, edgeIndex);
            break;
        }
    }
}

void SparseGraph::add_edge(const GraphEdge &graphEdge)
{
    m_graphEdges[graphEdge.elem1()].push_back(graphEdge);
}

void SparseGraph::delete_all_edges(LocalId elemId)
{
    m_graphEdges.erase(elemId);
}

} } }
