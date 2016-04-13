#include "SideConnector.hpp"
#include "ElemElemGraphImpl.hpp"
#include "ElemGraphCoincidentElems.hpp"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/baseImpl/elementGraph/BulkDataIdMapper.hpp>

namespace stk
{
namespace mesh
{

void SideConnector::connect_side_to_all_elements(stk::mesh::Entity sideEntity, stk::mesh::Entity elemEntity, int elemSide)
{
    //connect_side_to_elem(sideEntity, elemEntity, elemSide);
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
    stk::mesh::Entity otherElement = get_entity_for_local_id(graphEdge.elem2());
    if(m_bulk_data.is_valid(otherElement))
        connect_side_to_elem(sideEntity, otherElement, graphEdge.side2());
}

stk::mesh::Entity SideConnector::get_entity_for_local_id(stk::mesh::impl::LocalId localId) const
{
    if(impl::is_local_element(localId))
        return m_localMapper.local_to_entity(localId);
    else
        return m_bulk_data.get_entity(stk::topology::ELEM_RANK, -localId);
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

stk::mesh::Permutation get_permutation_of_side_nodes(stk::topology sideTopology,
                                                     const stk::mesh::Entity *sideNodes,
                                                     const stk::mesh::Entity *elemSideNodes)
{
    std::pair<bool, unsigned> result = sideTopology.equivalent(sideNodes, elemSideNodes);
    return static_cast<stk::mesh::Permutation>(result.second);
}

stk::mesh::Permutation SideConnector::get_permutation_for_side(stk::mesh::Entity sideEntity,
                                                               stk::mesh::Entity element,
                                                               int sideOrd)
{
    stk::mesh::EntityVector elemSideNodes = get_nodes_of_elem_side(element, sideOrd);
    return get_permutation_of_side_nodes(m_bulk_data.bucket(sideEntity).topology(),
                                         m_bulk_data.begin_nodes(sideEntity),
                                         elemSideNodes.data());
}

}
}
