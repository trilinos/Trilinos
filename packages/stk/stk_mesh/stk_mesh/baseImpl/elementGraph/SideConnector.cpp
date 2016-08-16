#include "SideConnector.hpp"
#include "ElemElemGraphImpl.hpp"
#include "ElemGraphCoincidentElems.hpp"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/baseImpl/elementGraph/BulkDataIdMapper.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>

namespace stk
{
namespace mesh
{

stk::mesh::Entity get_entity_for_local_id(const stk::mesh::BulkData &bulk, const stk::mesh::impl::ElementLocalIdMapper &localMapper, stk::mesh::impl::LocalId localId)
{
    if(impl::is_local_element(localId))
        return localMapper.local_to_entity(localId);
    else
        return bulk.get_entity(stk::topology::ELEM_RANK, -localId);
}


template <typename GraphEdgeCollection>
GraphEdge get_graph_edge_of_min_id(const GraphEdgeCollection& graphEdges,
                                                 int elemSide,
                                                 const stk::mesh::impl::IdMapper& idMapper)
{
    stk::mesh::impl::LocalId maxId = std::numeric_limits<stk::mesh::impl::LocalId>::max();
    int maxInt = std::numeric_limits<int>::max();
    GraphEdge edgeWithMinId(maxId, maxInt, -maxId, maxInt);
    stk::mesh::EntityId currentMinId = maxId;
    for(const GraphEdge& graphEdge : graphEdges)
    {
        if(graphEdge.side1() == elemSide)
        {
            stk::mesh::EntityId elem2Id = idMapper.localToGlobal(graphEdge.elem2());
            if(currentMinId > elem2Id)
            {
                currentMinId = elem2Id;
                edgeWithMinId = graphEdge;
            }
        }
    }
    return edgeWithMinId;
}

GraphEdge SideCreationElementChooser::get_chosen_graph_edge(stk::mesh::Entity elemEntity, int elemSide) const
{
    stk::mesh::impl::BulkDataIdMapper idMapper(bulk, localIdMapper);
    stk::mesh::impl::LocalId elemLocalId = localIdMapper.entity_to_local(elemEntity);

    GraphEdge localEdge(elemLocalId, elemSide, elemLocalId, elemSide);
    GraphEdge adjEdgeWithMinId = get_graph_edge_of_min_id(graph.get_edges_for_element(elemLocalId), elemSide, idMapper);
    GraphEdge coinEdgeWithMinId = get_graph_edge_of_min_id(coincidentGraph.get_edges_for_element(elemLocalId), elemSide, idMapper);

    stk::mesh::EntityId minLocId = idMapper.localToGlobal(localEdge.elem2());
    stk::mesh::EntityId minAdjId = idMapper.localToGlobal(adjEdgeWithMinId.elem2());
    stk::mesh::EntityId minCoinId = idMapper.localToGlobal(coinEdgeWithMinId.elem2());

    GraphEdge chosenEdge;
    if(minLocId < minAdjId && minLocId < minCoinId)
        chosenEdge = localEdge;
    else if(minAdjId < minCoinId)
        chosenEdge = adjEdgeWithMinId;
    else
        chosenEdge = coinEdgeWithMinId;
    return chosenEdge;
}

stk::mesh::EntityId SideIdChooser::get_chosen_side_id(stk::mesh::Entity elem, int elemSide) const
{
    stk::mesh::GraphEdge graphEdgeForElementToCreateSideOn = elemChooser.get_chosen_graph_edge(elem, elemSide);
    return impl::side_id_formula(idMapper.localToGlobal(graphEdgeForElementToCreateSideOn.elem2()), graphEdgeForElementToCreateSideOn.side2());
}

void SideNodeConnector::connect_side_to_nodes(stk::mesh::Entity sideEntity, stk::mesh::Entity elemEntity, int elemSide)
{
    SideCreationElementChooser elementChooser(bulk, localMapper, graph, coincidentGraph);
    GraphEdge graphEdgeForElementToCreateSideOn = elementChooser.get_chosen_graph_edge(elemEntity, elemSide);
    connect_side_to_other_elements_nodes(graphEdgeForElementToCreateSideOn , sideEntity, elemEntity, elemSide);
}

void SideNodeConnector::connect_side_to_elements_nodes(stk::mesh::Entity sideEntity, stk::mesh::Entity elemEntity, int elemSide)
{
    stk::mesh::EntityVector sideNodes;
    stk::mesh::impl::fill_element_side_nodes_from_topology(bulk, elemEntity, elemSide, sideNodes);
    for(size_t i = 0; i < sideNodes.size(); i++)
        bulk.declare_relation(sideEntity, sideNodes[i], i);
}

void SideNodeConnector::connect_side_to_other_elements_nodes(const GraphEdge &edgeWithMinId, stk::mesh::Entity sideEntity, stk::mesh::Entity elemEntity, int elemSide)
{
    if(impl::is_local_element(edgeWithMinId.elem2()))
    {
        stk::mesh::Entity elemWithMinId = get_entity_for_local_id(bulk, localMapper, edgeWithMinId.elem2());
        connect_side_to_elements_nodes(sideEntity, elemWithMinId, edgeWithMinId.side2());
    }
    else
    {
        stk::mesh::EntityVector sideNodes;
        stk::topology sideTop = bulk.bucket(elemEntity).topology().side_topology(elemSide);
        const stk::mesh::impl::ParallelInfo &parInfo = parallelGraph.get_parallel_info_for_graph_edge(edgeWithMinId);
        stk::mesh::impl::fill_element_side_nodes_from_topology(bulk, elemEntity, elemSide, sideNodes);
        stk::mesh::EntityVector permutedSideNodes(sideTop.num_nodes());
        sideTop.permutation_nodes(sideNodes, parInfo.m_permutation, permutedSideNodes.begin());
        for(size_t i = 0; i < permutedSideNodes.size(); ++i)
            bulk.declare_relation(sideEntity, permutedSideNodes[i], i);
    }
}

void SideConnector::connect_side_to_all_elements(stk::mesh::Entity sideEntity, stk::mesh::Entity elemEntity, int elemSide)
{
    connect_side_to_elem(sideEntity, elemEntity, elemSide);
    stk::mesh::impl::LocalId elemLocalId = m_localMapper.entity_to_local(elemEntity);
    connect_side_to_coincident_elements(sideEntity, elemLocalId, elemSide);
    connect_side_to_adjacent_elements(sideEntity, elemLocalId, elemSide);
}

void SideConnector::connect_side_to_elem(stk::mesh::Entity sideEntity,
                                         stk::mesh::Entity element,
                                         int sideOrd)
{
    stk::mesh::Permutation perm = get_permutation_for_side(sideEntity, element, sideOrd);
    m_bulk_data.declare_relation(element, sideEntity, sideOrd, perm);
}

void SideConnector::connect_side_to_adjacent_elements(stk::mesh::Entity sideEntity,
                                                      stk::mesh::impl::LocalId elemLocalId,
                                                      int elemSide)
{
    for(const GraphEdge& graphEdge : m_graph.get_edges_for_element(elemLocalId))
    {
        if(graphEdge.side1() == elemSide)
        {
            connect_side_entity_to_other_element(sideEntity, graphEdge);
            connect_side_to_coincident_elements(sideEntity, graphEdge.elem2(), graphEdge.side2());
        }
    }
}

void SideConnector::connect_side_entity_to_other_element(stk::mesh::Entity sideEntity,
                                                         const stk::mesh::GraphEdge &graphEdge)
{
    stk::mesh::Entity otherElement = get_entity_for_local_id(m_bulk_data, m_localMapper, graphEdge.elem2());
    if(m_bulk_data.is_valid(otherElement))
        connect_side_to_elem(sideEntity, otherElement, graphEdge.side2());
}

void SideConnector::connect_side_to_coincident_elements(stk::mesh::Entity sideEntity,
                                                        stk::mesh::impl::LocalId elemLocalId,
                                                        int elemSide)
{
    for(const stk::mesh::GraphEdge &graphEdge : m_coincidentGraph.get_edges_for_element(elemLocalId))
        if(graphEdge.side1() == elemSide)
            connect_side_entity_to_other_element(sideEntity, graphEdge);
}


stk::mesh::EntityVector SideConnector::get_nodes_of_elem_side(stk::mesh::Entity element, int sideOrd)
{
    stk::topology elemTopology = m_bulk_data.bucket(element).topology();
    stk::mesh::EntityVector nodesOfSide(elemTopology.side_topology(sideOrd).num_nodes());
    elemTopology.side_nodes(m_bulk_data.begin_nodes(element), sideOrd, nodesOfSide.data());
    return nodesOfSide;
}

stk::mesh::Permutation SideConnector::get_permutation_for_side(stk::mesh::Entity sideEntity,
                                                               stk::mesh::Entity element,
                                                               int sideOrd)
{
    std::pair<bool,unsigned> result = stk::mesh::side_equivalent(m_bulk_data, element, sideOrd, m_bulk_data.begin_nodes(sideEntity));
    return static_cast<stk::mesh::Permutation>(result.second);
}

}
}
