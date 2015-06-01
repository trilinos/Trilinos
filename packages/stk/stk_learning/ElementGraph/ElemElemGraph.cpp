#include "ElemElemGraph.hpp"
#include "ElemElemGraphImpl.hpp"

#include <vector>
#include <algorithm>

#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>

#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/environment/ReportHandler.hpp>


ElemElemGraph::ElemElemGraph(stk::mesh::BulkData& bulkData) : m_bulk_data(bulkData)
{
    size_data_members();

    impl::set_local_ids_and_fill_element_entities_and_topologies(m_bulk_data, local_id_to_element_entity, element_topologies);
    impl::fill_graph(m_bulk_data, elem_graph, via_sides);
    impl::fill_parallel_graph(m_bulk_data, elem_graph, via_sides, parallel_graph_info);
}

ElemElemGraph::~ElemElemGraph() {}

size_t ElemElemGraph::get_num_connected_elems(stk::mesh::Entity local_element) const
{
    LocalId local_id = get_local_element_id(local_element);
    return elem_graph[local_id].size();
}

bool ElemElemGraph::is_connected_elem_locally_owned(stk::mesh::Entity local_element, size_t index_conn_elem) const
{
    LocalId local_id = get_local_element_id(local_element);
    return elem_graph[local_id][index_conn_elem] >= 0;
}

int ElemElemGraph::get_side_id_to_connected_element(stk::mesh::Entity local_element, size_t index_conn_elem) const
{
    LocalId local_id = get_local_element_id(local_element);
    return via_sides[local_id][index_conn_elem];
}

stk::mesh::Entity ElemElemGraph::get_connected_element(stk::mesh::Entity local_element, size_t index_conn_elem) const
{
    LocalId local_id = get_local_element_id(local_element);
    LocalId other_element_id = elem_graph[local_id][index_conn_elem];
    return local_id_to_element_entity[other_element_id];
}

stk::mesh::EntityId ElemElemGraph::get_entity_id_of_remote_element(stk::mesh::Entity local_element, size_t index_conn_elem) const
{
    ThrowRequireMsg(!is_connected_elem_locally_owned(local_element, index_conn_elem) , "Program error. Contact sierra-help@sandia.gov for support.");
    LocalId local_id = get_local_element_id(local_element);
    stk::mesh::EntityId id = -elem_graph[local_id][index_conn_elem];
    return id;
}

int ElemElemGraph::get_owning_proc_id_of_remote_element(stk::mesh::Entity local_element, stk::mesh::EntityId other_element_id) const
{
    LocalId local_id = get_local_element_id(local_element);
    ParallelGraphInfo::const_iterator iter = parallel_graph_info.find(std::make_pair(local_id, other_element_id));
    ThrowRequireMsg( iter != parallel_graph_info.end(), "Program error. Contact sierra-help@sandia.gov for support.");
    int other_proc = iter->second.m_other_proc;
    return other_proc;
}

int ElemElemGraph::get_side_from_element1_to_remote_element2(stk::mesh::Entity local_element, stk::mesh::EntityId other_element_id) const
{
    LocalId remote_element_local_id = -other_element_id;
    LocalId element1_local_id = get_local_element_id(local_element);

    int side = -1;
    const std::vector<LocalId>& conn_elements = elem_graph[element1_local_id];

    std::vector<LocalId>::const_iterator iter = std::find(conn_elements.begin(), conn_elements.end(), remote_element_local_id);
    if ( iter != conn_elements.end() )
    {
        int64_t index = iter - conn_elements.begin();
        side = via_sides[element1_local_id][index];
    }
    return side;
}

int ElemElemGraph::get_side_from_element1_to_locally_owned_element2(stk::mesh::Entity local_element, stk::mesh::Entity other_element) const
{
    LocalId other_element_id = get_local_element_id(other_element);
    LocalId element1_local_id = get_local_element_id(local_element);

    int side = -1;
    const std::vector<LocalId>& conn_elements = elem_graph[element1_local_id];

    std::vector<LocalId>::const_iterator iter = std::find(conn_elements.begin(), conn_elements.end(), other_element_id);
    if ( iter != conn_elements.end() )
    {
        int64_t index = iter - conn_elements.begin();
        side = via_sides[element1_local_id][index];
    }
    return side;
}

parallel_info& ElemElemGraph::get_parallel_edge_info(stk::mesh::Entity element, stk::mesh::EntityId remote_id)
{
    LocalId this_elem_local_id = get_local_element_id(element);

    ParallelGraphInfo::iterator iter = parallel_graph_info.find(std::make_pair(this_elem_local_id, remote_id));
    ThrowRequireMsg( iter != parallel_graph_info.end(), "Program error. Contact sierra-help@sandia.gov for support.");
    return iter->second;
}

LocalId ElemElemGraph::get_local_element_id(stk::mesh::Entity local_element) const
{
    ThrowRequireMsg(m_bulk_data.is_valid(local_element), "Program error. Contact sierra-help@sandia.gov for support.");
    LocalId local_id = m_bulk_data.local_id(local_element);
    return local_id;
}

void ElemElemGraph::size_data_members()
{
    std::vector<unsigned> counts;
    stk::mesh::count_entities(m_bulk_data.mesh_meta_data().locally_owned_part(), m_bulk_data, counts);
    int numElems = counts[stk::topology::ELEM_RANK];

    elem_graph.resize(numElems);
    via_sides.resize(numElems);
    local_id_to_element_entity.resize(numElems, 0);
    element_topologies.resize(numElems);
}
