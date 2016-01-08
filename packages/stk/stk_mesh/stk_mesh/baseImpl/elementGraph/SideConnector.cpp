#include "SideConnector.hpp"
#include "ElemElemGraphImpl.hpp"
#include "ElemGraphCoincidentElems.hpp"
#include <stk_mesh/base/BulkData.hpp>

namespace stk
{
namespace mesh
{

bool do_sides_match_and_elem2_is_local(const stk::mesh::GraphEdge &graphEdge, int side)
{
    return graphEdge.side1 == side && impl::is_local_element(graphEdge.elem2);
}

void SideConnector::connect_side_to_all_elements(stk::mesh::Entity sideEntity,
                                                 impl::ElementSidePair skinnedElemSidePair,
                                                 stk::mesh::EntityVector &skinned_elements)
{
    for(const GraphEdge & graphEdge : m_graph.get_edges_for_element(skinnedElemSidePair.first))
    {
        if(do_sides_match_and_elem2_is_local(graphEdge, skinnedElemSidePair.second))
        {
            connect_side_entity_to_other_element(sideEntity, graphEdge, skinned_elements);
            connect_side_to_coincident_elements(sideEntity, {graphEdge.elem2, graphEdge.side2}, skinned_elements);
        }
    }
    connect_side_to_coincident_elements(sideEntity, skinnedElemSidePair, skinned_elements);
}

void SideConnector::connect_side_entity_to_other_element(stk::mesh::Entity sideEntity,
                                                         const stk::mesh::GraphEdge &graphEdge,
                                                         stk::mesh::EntityVector skinned_elements)
{
    stk::mesh::Entity other_element = m_local_id_to_element_entity[graphEdge.elem2];
    stk::mesh::Permutation perm =
            static_cast<stk::mesh::Permutation>(m_bulk_data.bucket(other_element).topology().num_positive_permutations());
    m_bulk_data.declare_relation(other_element, sideEntity, graphEdge.side2, perm);
    skinned_elements.push_back(other_element);
}

void SideConnector::connect_side_to_coincident_elements(stk::mesh::Entity sideEntity,
                                                        impl::ElementSidePair skinnedElemSidePair,
                                                        stk::mesh::EntityVector &skinned_elements)
{
    auto iter = m_coincidentGraph.find(skinnedElemSidePair.first);
    if(iter != m_coincidentGraph.end())
        for(const stk::mesh::GraphEdge &graphEdge : iter->second)
            if(do_sides_match_and_elem2_is_local(graphEdge, skinnedElemSidePair.second))
                connect_side_entity_to_other_element(sideEntity, graphEdge, skinned_elements);
}

impl::LocalId SideConnector::get_local_element_id(stk::mesh::Entity local_element, bool require_valid_id) const
{
    ThrowRequireMsg(m_entity_to_local_id.size() > local_element.local_offset(),"Program error. Contact sierra-help@sandia.gov for support.");
    impl::LocalId local_id = m_entity_to_local_id[local_element.local_offset()];
    if (require_valid_id)
    {
        ThrowRequireMsg(local_id != impl::INVALID_LOCAL_ID, "Program error. Contact sierra-help@sandia.gov for support.");
    }
    return local_id;
}

bool SideConnector::has_remote_graph_edge(stk::mesh::Entity localEntity,
                                          int side,
                                          stk::mesh::Entity remoteEntity)
{
    impl::LocalId localIdForLocalEntity = get_local_element_id(localEntity, true);
    stk::mesh::EntityId remoteEntityId = m_bulk_data.identifier(remoteEntity);
    std::vector<GraphEdge> graphEdges = m_graph.get_edges_for_element_side(localIdForLocalEntity, side);

    for(const GraphEdge &graphEdge : graphEdges)
    {
        impl::LocalId localIdForRemoteEntity = -remoteEntityId;
        if(!impl::is_local_element(graphEdge.elem2) && (graphEdge.elem2 == localIdForRemoteEntity))
            return true;
    }

    return false;
}

}
}
