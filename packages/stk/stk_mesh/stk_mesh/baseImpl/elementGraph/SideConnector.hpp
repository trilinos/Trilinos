#ifndef STK_SIDE_CONNECTOR_HPP
#define STK_SIDE_CONNECTOR_HPP

#include <vector>
#include <stk_mesh/base/Types.hpp>
#include "ElemElemGraphImpl.hpp"
#include "ElemGraphCoincidentElems.hpp"
#include "GraphEdgeData.hpp"

namespace stk { namespace mesh { class BulkData; } }

namespace stk
{
namespace mesh
{

class SideConnector
{
public:
    SideConnector(stk::mesh::BulkData &b,
                  const stk::mesh::Graph &g,
                  const stk::mesh::impl::SparseGraph &cg,
                  const stk::mesh::EntityVector &localToElement,
                  const std::vector<impl::LocalId> &elemToLocal) :
            m_bulk_data(b),
            m_graph(g),
            m_coincidentGraph(cg),
            m_local_id_to_element_entity(localToElement),
            m_entity_to_local_id(elemToLocal)
    {
    }
    void connect_side_to_all_elements(stk::mesh::Entity sideEntity,
                                      impl::ElementSidePair skinnedElemSidePair,
                                      stk::mesh::EntityVector &skinned_elements);
    bool has_remote_graph_edge(stk::mesh::Entity localEntity,
                               int side,
                               stk::mesh::Entity remoteEntity);


private:
    impl::LocalId get_local_element_id(stk::mesh::Entity local_element, bool require_valid_id) const;

    void connect_side_entity_to_other_element(stk::mesh::Entity sideEntity,
                                              const stk::mesh::GraphEdge &graphEdge,
                                              stk::mesh::EntityVector skinned_elements);
    void connect_side_to_coincident_elements(stk::mesh::Entity sideEntity,
                                             impl::ElementSidePair skinnedElemSidePair,
                                             stk::mesh::EntityVector &skinned_elements);
    stk::mesh::BulkData &m_bulk_data;
    const stk::mesh::Graph &m_graph;
    const stk::mesh::impl::SparseGraph &m_coincidentGraph;
    const stk::mesh::EntityVector &m_local_id_to_element_entity;
    const std::vector<impl::LocalId> &m_entity_to_local_id;
};

}
}

#endif
