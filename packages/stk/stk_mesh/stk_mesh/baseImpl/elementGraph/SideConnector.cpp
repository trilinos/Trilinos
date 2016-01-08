#include "SideConnector.hpp"
#include "ElemElemGraphImpl.hpp"
#include "ElemGraphCoincidentElems.hpp"
#include <stk_mesh/base/BulkData.hpp>

namespace stk
{
namespace mesh
{

void SideConnector::connect_side_to_all_elements(stk::mesh::Entity sideEntity,
                                                 impl::ElementSidePair skinnedElemSidePair)
{
    for(const GraphEdge & graphEdge : m_graph.get_edges_for_element(skinnedElemSidePair.first))
    {
        int skinnedElementSide = skinnedElemSidePair.second;
        if(graphEdge.side1 == skinnedElementSide)
        {
            connect_side_entity_to_other_element(sideEntity, graphEdge);
            connect_side_to_coincident_elements(sideEntity, {graphEdge.elem2, graphEdge.side2});
        }
    }
    connect_side_to_coincident_elements(sideEntity, skinnedElemSidePair);
}

void SideConnector::connect_side_entity_to_other_element(stk::mesh::Entity sideEntity,
                                                         const stk::mesh::GraphEdge &graphEdge)
{
    stk::mesh::Entity otherElement = get_entity_for_local_id(graphEdge.elem2);
    if(m_bulk_data.is_valid(otherElement))
    {
        stk::mesh::Permutation perm = static_cast<stk::mesh::Permutation>(m_bulk_data.bucket(otherElement).topology().num_positive_permutations());
        m_bulk_data.declare_relation(otherElement, sideEntity, graphEdge.side2, perm);
    }
}

stk::mesh::Entity SideConnector::get_entity_for_local_id(stk::mesh::impl::LocalId localId) const
{
    if(impl::is_local_element(localId))
        return m_local_id_to_element_entity[localId];
    else
        return m_bulk_data.get_entity(stk::topology::ELEM_RANK, -localId);
}

void SideConnector::connect_side_to_coincident_elements(stk::mesh::Entity sideEntity,
                                                        impl::ElementSidePair skinnedElemSidePair)
{
    auto iter = m_coincidentGraph.find(skinnedElemSidePair.first);
    if(iter != m_coincidentGraph.end())
        for(const stk::mesh::GraphEdge &graphEdge : iter->second)
        {
            int skinnedElementSide = skinnedElemSidePair.second;
            if(graphEdge.side1 == skinnedElementSide)
                connect_side_entity_to_other_element(sideEntity, graphEdge);
        }
}

}
}
