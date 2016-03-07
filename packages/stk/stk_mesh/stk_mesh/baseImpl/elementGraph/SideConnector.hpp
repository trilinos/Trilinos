#ifndef STK_SIDE_CONNECTOR_HPP
#define STK_SIDE_CONNECTOR_HPP

#include <vector>
#include <stk_mesh/base/Types.hpp>
#include "BulkDataIdMapper.hpp"
#include "ElemElemGraphImpl.hpp"
#include "ElemGraphCoincidentElems.hpp"
#include "GraphEdgeData.hpp"

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { namespace impl { class ElementLocalIdMapper; } } }

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
                  const stk::mesh::impl::ElementLocalIdMapper & localMapper) :
            m_bulk_data(b),
            m_graph(g),
            m_coincidentGraph(cg),
            m_localMapper(localMapper)
    {
    }
    void connect_side_to_all_elements(stk::mesh::Entity sideEntity,
                                      stk::mesh::Entity elemEntity,
                                      int elemSide);
private:
    SideConnector();

    stk::mesh::Entity get_entity_for_local_id(stk::mesh::impl::LocalId localId) const;

    void connect_side_entity_to_other_element(stk::mesh::Entity sideEntity,
                                              const stk::mesh::GraphEdge &graphEdge);
    void connect_side_to_coincident_elements(stk::mesh::Entity sideEntity, stk::mesh::impl::LocalId elemLocalId, int elemSide);
    stk::mesh::Permutation get_permutation_for_side(stk::mesh::Entity sideEntity,
                                                    stk::mesh::Entity element,
                                                    int sideOrd);
    stk::mesh::EntityVector get_nodes_of_elem_side(stk::mesh::Entity element, int sideOrd);
    void connect_side_to_elem(stk::mesh::Entity sideEntity,
                              stk::mesh::Entity otherElement,
                              int sideOrd);
    void connect_side_to_adjacent_elements(stk::mesh::Entity sideEntity, stk::mesh::impl::LocalId elemLocalId, int elemSide);

    stk::mesh::BulkData &m_bulk_data;
    const stk::mesh::Graph &m_graph;
    const stk::mesh::impl::SparseGraph &m_coincidentGraph;
    const stk::mesh::impl::ElementLocalIdMapper & m_localMapper;
};

}
}

#endif
