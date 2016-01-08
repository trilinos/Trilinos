#include "SideConnector.hpp"
#include "ElemElemGraphImpl.hpp"
#include "ElemGraphCoincidentElems.hpp"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace stk
{
namespace mesh
{

void SideConnector::connect_side_to_all_elements(stk::mesh::Entity sideEntity,
                                                 impl::ElementSidePair skinnedElemSidePair)
{
    connect_side_to_elem(sideEntity, m_local_id_to_element_entity[skinnedElemSidePair.first], skinnedElemSidePair.second);

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

stk::mesh::Permutation SideConnector::get_permutation_for_side(stk::mesh::Entity sideEntity,
                                                               stk::mesh::Entity element,
                                                               int sideOrd)
{
    stk::mesh::EntityRank sideRank = m_bulk_data.mesh_meta_data().side_rank();
    stk::topology elemTopology = m_bulk_data.bucket(element).topology();
    stk::topology sideTopology = elemTopology.sub_topology(sideRank, sideOrd);
    stk::mesh::EntityVector nodesOfSubTopo(sideTopology.num_nodes());
    elemTopology.sub_topology_nodes(m_bulk_data.begin_nodes(element), sideRank, sideOrd, nodesOfSubTopo.data());
    std::pair<bool, unsigned> result = sideTopology.equivalent(m_bulk_data.begin_nodes(sideEntity), nodesOfSubTopo);
    stk::mesh::Permutation perm = static_cast<stk::mesh::Permutation>(result.second);
    return perm;
}

void SideConnector::connect_side_to_elem(stk::mesh::Entity sideEntity,
                                         stk::mesh::Entity element,
                                         int sideOrd)
{
    stk::mesh::Permutation perm = get_permutation_for_side(sideEntity, element, sideOrd);
    m_bulk_data.declare_relation(element, sideEntity, sideOrd, perm);
}

void SideConnector::connect_side_entity_to_other_element(stk::mesh::Entity sideEntity,
                                                         const stk::mesh::GraphEdge &graphEdge)
{
    stk::mesh::Entity otherElement = get_entity_for_local_id(graphEdge.elem2);
    if(m_bulk_data.is_valid(otherElement))
        connect_side_to_elem(sideEntity, otherElement, graphEdge.side2);
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
