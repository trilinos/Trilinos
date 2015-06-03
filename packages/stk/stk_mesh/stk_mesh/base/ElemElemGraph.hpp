#ifndef STK_ELEM_ELEM_GRAPH_HPP
#define STK_ELEM_ELEM_GRAPH_HPP

#include <vector>
#include <map>
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Types.hpp>

#include "ElemElemGraphImpl.hpp"

namespace stk { namespace mesh { class BulkData; } }

namespace stk { namespace mesh {


class ElemElemGraph
{
public:
    ElemElemGraph(stk::mesh::BulkData& bulkData);

    ~ElemElemGraph();

    size_t get_num_connected_elems(stk::mesh::Entity local_element) const;

    bool is_connected_elem_locally_owned(stk::mesh::Entity local_element, size_t index_conn_elem) const;

    stk::mesh::Entity get_connected_element(stk::mesh::Entity local_element, size_t index_conn_elem) const;

    stk::mesh::EntityId get_entity_id_of_remote_element(stk::mesh::Entity local_element, size_t index_conn_elem) const;

    int get_owning_proc_id_of_remote_element(stk::mesh::Entity local_element, stk::mesh::EntityId other_element_id) const;

    int get_side_id_to_connected_element(stk::mesh::Entity local_element, size_t index_conn_elem) const;

    int get_side_from_element1_to_remote_element2(stk::mesh::Entity local_element, stk::mesh::EntityId other_element_id) const;

    int get_side_from_element1_to_locally_owned_element2(stk::mesh::Entity local_element, stk::mesh::Entity other_element) const;

    impl::parallel_info& get_parallel_edge_info(stk::mesh::Entity element, stk::mesh::EntityId remote_id);

    size_t num_edges() const;

    size_t num_parallel_edges() const { return m_num_parallel_edges; }

private:

    impl::LocalId get_local_element_id(stk::mesh::Entity local_element) const;
    void size_data_members();

    stk::mesh::BulkData &m_bulk_data;
    impl::ElementGraph elem_graph;
    impl::SidesForElementGraph via_sides;
    impl::ParallelGraphInfo parallel_graph_info;
    stk::mesh::EntityVector local_id_to_element_entity;
    std::vector<stk::topology> element_topologies;
    size_t m_num_edges;
    size_t m_num_parallel_edges;
};

void perform_element_death(stk::mesh::BulkData& bulkData, ElemElemGraph& elementGraph, const stk::mesh::EntityVector& killedElements, stk::mesh::Part& active,
        const stk::mesh::PartVector& boundary_mesh_parts);

}} // end stk mesh namespaces

#endif
