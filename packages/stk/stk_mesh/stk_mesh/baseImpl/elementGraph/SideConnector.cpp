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

bool is_global_id_less(const stk::mesh::impl::BulkDataIdMapper &idMapper, stk::mesh::impl::LocalId a, stk::mesh::impl::LocalId b)
{
    return (idMapper.localToGlobal(a) < idMapper.localToGlobal(b));
}

template <typename GraphEdgeCollection>
void update_lowest_elem_id_among_edges(int elemSide, const stk::mesh::impl::BulkDataIdMapper &idMapper, const GraphEdgeCollection& graphEdges, stk::mesh::GraphEdge &currentMinEdge)
{
    for(const GraphEdge &graphEdge : graphEdges)
        if(graphEdge.side1() == elemSide)
            if(is_global_id_less(idMapper, graphEdge.elem2(), currentMinEdge.elem2()))
                currentMinEdge = graphEdge;
}

GraphEdge SideCreationElementChooser::get_chosen_graph_edge(stk::mesh::Entity elemEntity, int elemSide) const
{
    stk::mesh::impl::BulkDataIdMapper idMapper(bulk, localIdMapper);
    stk::mesh::impl::LocalId elemLocalId = localIdMapper.entity_to_local(elemEntity);

    GraphEdge currentMinEdge(elemLocalId, elemSide, elemLocalId, elemSide);
    update_lowest_elem_id_among_edges(elemSide, idMapper, graph.get_edges_for_element(elemLocalId), currentMinEdge);
    update_lowest_elem_id_among_edges(elemSide, idMapper, coincidentGraph.get_edges_for_element(elemLocalId), currentMinEdge);

    return currentMinEdge;
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

void SideNodeConnector::declare_relations_to_nodes(stk::mesh::Entity sideEntity, const stk::mesh::EntityVector &sideNodes)
{
    stk::mesh::Permutation perm = stk::mesh::Permutation::INVALID_PERMUTATION;
    for(size_t i = 0; i < sideNodes.size(); i++) {
        bulk.declare_relation(sideEntity, sideNodes[i], i, perm,
                              m_scratchOrdinals1, m_scratchOrdinals2, m_scratchOrdinals3);
    }
}

void SideNodeConnector::connect_side_to_elements_nodes(stk::mesh::Entity sideEntity, stk::mesh::Entity elemEntity, int elemSide)
{
    stk::mesh::EntityVector sideNodes;
    stk::mesh::impl::fill_element_side_nodes_from_topology(bulk.bucket(elemEntity).topology(), bulk.begin_nodes(elemEntity), elemSide, sideNodes);
    declare_relations_to_nodes(sideEntity, sideNodes);
}

stk::mesh::EntityVector SideNodeConnector::get_permuted_side_nodes(stk::mesh::Entity elemEntity, int elemSide, const stk::mesh::EntityVector &sideNodes, int permutation)
{
    stk::topology sideTop = bulk.bucket(elemEntity).topology().side_topology(elemSide);
    stk::mesh::EntityVector permutedSideNodes(sideTop.num_nodes());
    sideTop.permutation_nodes(sideNodes.data(), permutation, permutedSideNodes.data());
    return permutedSideNodes;
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
        stk::mesh::impl::fill_element_side_nodes_from_topology(bulk.bucket(elemEntity).topology(), bulk.begin_nodes(elemEntity), elemSide, sideNodes);

        const stk::mesh::impl::ParallelInfo &parInfo = parallelGraph.get_parallel_info_for_graph_edge(edgeWithMinId);
        stk::mesh::EntityVector permutedSideNodes = get_permuted_side_nodes(elemEntity, elemSide, sideNodes, parInfo.m_permutation);
        declare_relations_to_nodes(sideEntity, permutedSideNodes);
    }
}

void check_entity_has_local_id(const stk::mesh::BulkData &bulk, const stk::mesh::impl::ElementLocalIdMapper &localMapper, stk::mesh::Entity elemEntity)
{
  STK_ThrowRequireMsg(localMapper.does_entity_have_local_id(elemEntity),
      "no local id for " << (bulk.is_valid(elemEntity) ? "valid" : "invalid")
                << " elem " << bulk.identifier(elemEntity)
                << ", local_offset=" << elemEntity.local_offset());
}

void SideConnector::connect_side_to_all_elements(stk::mesh::Entity sideEntity, stk::mesh::Entity elemEntity, int elemSide)
{
    check_entity_has_local_id(m_bulk_data, m_localMapper, elemEntity);
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

    unsigned rankedOrdinal;
    stk::topology::rank_t sideRank;
    auto elemTopo = m_bulk_data.bucket(element).topology();
    elemTopo.ranked_side_ordinal(sideOrd, rankedOrdinal, sideRank);

    m_bulk_data.declare_relation(element, sideEntity, rankedOrdinal, perm, m_scratchOrdinals1, m_scratchOrdinals2, m_scratchOrdinals3);
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
    stk::EquivalentPermutation result = stk::mesh::side_equivalent(m_bulk_data, element, sideOrd, m_bulk_data.begin_nodes(sideEntity));
    return static_cast<stk::mesh::Permutation>(result.permutation_number);
}

}
}
