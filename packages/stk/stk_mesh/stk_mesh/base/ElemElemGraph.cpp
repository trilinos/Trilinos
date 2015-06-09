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

namespace stk { namespace mesh {


ElemElemGraph::ElemElemGraph(stk::mesh::BulkData& bulkData) : m_bulk_data(bulkData)
{
    size_data_members();

    impl::set_local_ids_and_fill_element_entities_and_topologies(m_bulk_data, local_id_to_element_entity, element_topologies);
    impl::fill_graph(m_bulk_data, elem_graph, via_sides);

    impl::ElemSideToProcAndFaceId elem_side_comm = impl::get_element_side_ids_to_communicate(bulkData);
    size_t num_face_ids_needed = elem_side_comm.size();
    for(size_t i=0;i<via_sides.size();++i)
    {
        m_num_edges += via_sides[i].size();
    }
    num_face_ids_needed += m_num_edges;

    bulkData.generate_new_ids(stk::topology::FACE_RANK, num_face_ids_needed, suggested_face_ids);

    impl::fill_parallel_graph(m_bulk_data, elem_graph, via_sides, parallel_graph_info, elem_side_comm, suggested_face_ids);

    m_num_parallel_edges = parallel_graph_info.size();
    set_num_face_ids_used(m_num_parallel_edges);
    m_num_edges += m_num_parallel_edges;
}

const std::vector<stk::mesh::EntityId>& ElemElemGraph::get_suggested_face_ids() const
{
    return suggested_face_ids;
}

void ElemElemGraph::set_num_face_ids_used(size_t num_used)
{
    suggested_face_ids.erase(suggested_face_ids.begin(), suggested_face_ids.begin()+num_used);
}

ElemElemGraph::~ElemElemGraph() {}

size_t ElemElemGraph::get_num_connected_elems(stk::mesh::Entity local_element) const
{
    impl::LocalId local_id = get_local_element_id(local_element);
    return elem_graph[local_id].size();
}

bool ElemElemGraph::is_connected_elem_locally_owned(stk::mesh::Entity local_element, size_t index_conn_elem) const
{
    impl::LocalId local_id = get_local_element_id(local_element);
    return elem_graph[local_id][index_conn_elem] >= 0;
}

int ElemElemGraph::get_side_id_to_connected_element(stk::mesh::Entity local_element, size_t index_conn_elem) const
{
    impl::LocalId local_id = get_local_element_id(local_element);
    return via_sides[local_id][index_conn_elem];
}

stk::mesh::Entity ElemElemGraph::get_connected_element(stk::mesh::Entity local_element, size_t index_conn_elem) const
{
    impl::LocalId local_id = get_local_element_id(local_element);
    impl::LocalId other_element_id = elem_graph[local_id][index_conn_elem];
    return local_id_to_element_entity[other_element_id];
}

stk::mesh::EntityId ElemElemGraph::get_entity_id_of_remote_element(stk::mesh::Entity local_element, size_t index_conn_elem) const
{
    ThrowRequireMsg(!is_connected_elem_locally_owned(local_element, index_conn_elem) , "Program error. Contact sierra-help@sandia.gov for support.");
    impl::LocalId local_id = get_local_element_id(local_element);
    stk::mesh::EntityId id = -elem_graph[local_id][index_conn_elem];
    return id;
}

int ElemElemGraph::get_owning_proc_id_of_remote_element(stk::mesh::Entity local_element, stk::mesh::EntityId other_element_id) const
{
    impl::LocalId local_id = get_local_element_id(local_element);
    impl::ParallelGraphInfo::const_iterator iter = parallel_graph_info.find(std::make_pair(local_id, other_element_id));
    ThrowRequireMsg( iter != parallel_graph_info.end(), "Program error. Contact sierra-help@sandia.gov for support.");
    int other_proc = iter->second.m_other_proc;
    return other_proc;
}

int ElemElemGraph::get_side_from_element1_to_remote_element2(stk::mesh::Entity local_element, stk::mesh::EntityId other_element_id) const
{
    impl::LocalId remote_element_local_id = -other_element_id;
    impl::LocalId element1_local_id = get_local_element_id(local_element);

    int side = -1;
    const std::vector<impl::LocalId>& conn_elements = elem_graph[element1_local_id];

    std::vector<impl::LocalId>::const_iterator iter = std::find(conn_elements.begin(), conn_elements.end(), remote_element_local_id);
    if ( iter != conn_elements.end() )
    {
        int64_t index = iter - conn_elements.begin();
        side = via_sides[element1_local_id][index];
    }
    return side;
}

int ElemElemGraph::get_side_from_element1_to_locally_owned_element2(stk::mesh::Entity local_element, stk::mesh::Entity other_element) const
{
    impl::LocalId other_element_id = get_local_element_id(other_element);
    impl::LocalId element1_local_id = get_local_element_id(local_element);

    int side = -1;
    const std::vector<impl::LocalId>& conn_elements = elem_graph[element1_local_id];

    std::vector<impl::LocalId>::const_iterator iter = std::find(conn_elements.begin(), conn_elements.end(), other_element_id);
    if ( iter != conn_elements.end() )
    {
        int64_t index = iter - conn_elements.begin();
        side = via_sides[element1_local_id][index];
    }
    return side;
}

size_t ElemElemGraph::num_edges() const
{
    return m_num_edges;
}

impl::parallel_info& ElemElemGraph::get_parallel_edge_info(stk::mesh::Entity element, stk::mesh::EntityId remote_id)
{
    impl::LocalId this_elem_local_id = get_local_element_id(element);

    impl::ParallelGraphInfo::iterator iter = parallel_graph_info.find(std::make_pair(this_elem_local_id, remote_id));
    ThrowRequireMsg( iter != parallel_graph_info.end(), "Program error. Contact sierra-help@sandia.gov for support.");
    return iter->second;
}

impl::LocalId ElemElemGraph::get_local_element_id(stk::mesh::Entity local_element) const
{
    ThrowRequireMsg(m_bulk_data.is_valid(local_element), "Program error. Contact sierra-help@sandia.gov for support.");
    impl::LocalId local_id = m_bulk_data.local_id(local_element);
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
    m_num_edges = 0;
    m_num_parallel_edges = 0;
}

void perform_element_death(stk::mesh::BulkData& bulkData, ElemElemGraph& elementGraph, const stk::mesh::EntityVector& killedElements, stk::mesh::Part& active,
        const stk::mesh::PartVector& boundary_mesh_parts)
{
    const std::vector<stk::mesh::EntityId>& requestedIds = elementGraph.get_suggested_face_ids();
    size_t id_counter = 0;

    std::vector<stk::mesh::sharing_info> shared_modified;
    stk::mesh::EntityVector deletedEntities;

    std::vector<impl::graphEdgeProc> elements_to_comm = impl::get_elements_to_communicate(bulkData, killedElements, elementGraph);
    std::vector<std::pair<stk::mesh::EntityId, stk::mesh::EntityId> > remote_edges;

    impl::communicate_killed_entities(bulkData, elements_to_comm, remote_edges);

    stk::mesh::PartVector add_parts = boundary_mesh_parts;
    add_parts.push_back(&active);

    bulkData.modification_begin();

    for(size_t re = 0; re < remote_edges.size(); ++re)
    {
        stk::mesh::EntityId local_id = remote_edges[re].first;
        stk::mesh::EntityId remote_id = remote_edges[re].second;

        stk::mesh::Entity element = bulkData.get_entity(stk::topology::ELEM_RANK, local_id);
        bool create_face = true;
        if(!bulkData.bucket(element).member(active))
        {
            create_face = false;
        }

        impl::parallel_info &parallel_edge_info = elementGraph.get_parallel_edge_info(element, remote_id);
        impl::create_or_delete_shared_face(bulkData, parallel_edge_info, elementGraph, element, remote_id, create_face, add_parts,
                shared_modified, deletedEntities, id_counter, requestedIds[id_counter]);

        parallel_edge_info.m_in_part = false;
    }

    std::vector<impl::ElementSidePair> element_side_pairs;
    element_side_pairs.reserve(impl::get_element_face_multiplier() * killedElements.size());

    for(size_t k = 0; k < killedElements.size(); ++k)
    {
        stk::mesh::Entity this_elem_entity = killedElements[k];

        for(size_t j = 0; j < elementGraph.get_num_connected_elems(this_elem_entity); ++j)
        {
            if(elementGraph.is_connected_elem_locally_owned(this_elem_entity, j))
            {
                stk::mesh::Entity other_element = elementGraph.get_connected_element(this_elem_entity, j);
                int side_id = elementGraph.get_side_id_to_connected_element(this_elem_entity, j);
                ThrowRequireMsg(side_id != -1, "Program error. Please contact sierra-help@sandia.gov for support.");

                bool is_other_element_alive = bulkData.bucket(other_element).member(active);
                if(is_other_element_alive)
                {
                    // create or delete a face with a particular id

                    stk::topology side_top = bulkData.bucket(this_elem_entity).topology().side_topology(side_id);
                    stk::mesh::PartVector parts = add_parts;
                    parts.push_back(&bulkData.mesh_meta_data().get_topology_root_part(side_top));

                    std::string msg = "Program error. Please contact sierra-help@sandia.gov for support.";

                    ThrowRequireMsg(!impl::does_side_exist_with_different_permutation(bulkData, this_elem_entity,
                            static_cast<stk::mesh::ConnectivityOrdinal>(side_id), static_cast<stk::mesh::Permutation>(0)), msg);

                    ThrowRequireMsg(!impl::does_element_side_exist(bulkData, this_elem_entity, static_cast<stk::mesh::ConnectivityOrdinal>(side_id)), msg);

                    stk::mesh::Entity face = stk::mesh::impl::get_face_for_element_side(bulkData, this_elem_entity, side_id);
                    if(!bulkData.is_valid(face))
                    {
                        stk::mesh::EntityId face_global_id = requestedIds[id_counter];
                        ++id_counter;
                        ThrowRequireMsg(!impl::is_id_already_in_use_locally(bulkData, bulkData.mesh_meta_data().side_rank(), face_global_id), msg);

                        face = stk::mesh::declare_element_side(bulkData, face_global_id, this_elem_entity, side_id, parts);
                    }
                    else
                    {
                        bulkData.change_entity_parts(face, parts, stk::mesh::PartVector());
                    }

                    const stk::mesh::Entity* side_nodes = bulkData.begin_nodes(face);
                    unsigned num_side_nodes = bulkData.num_nodes(face);
                    stk::mesh::EntityVector side_nodes_vec(side_nodes, side_nodes + num_side_nodes);

                    std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> ord_and_perm =
                            stk::mesh::get_ordinal_and_permutation(bulkData, other_element, stk::topology::FACE_RANK, side_nodes_vec);

                    if(ord_and_perm.first == stk::mesh::INVALID_CONNECTIVITY_ORDINAL)
                    {
                        std::ostringstream yo;
                        yo << "Proc: " << bulkData.parallel_rank() << std::endl;
                        yo << "this element: " << bulkData.identifier(this_elem_entity) << std::endl;
                        yo << "other element: " << bulkData.identifier(other_element) << std::endl;
                        yo << "Nodes: ";

                        for(stk::mesh::Entity side_node : side_nodes_vec)
                        {
                            yo << bulkData.identifier(side_node) << " "; // nodes of elem 960 (other_elem: 1, this_elem: 401)
                        }

                        yo << std::endl;
                        std::cerr << yo.str();
                    }

                    ThrowRequireMsg(ord_and_perm.first != stk::mesh::INVALID_CONNECTIVITY_ORDINAL, "yikes!");
                    ThrowRequireMsg(ord_and_perm.second != stk::mesh::INVALID_PERMUTATION, "yikes!");

                    bulkData.declare_relation(other_element, face, ord_and_perm.first, ord_and_perm.second);
                }
                else
                {
                    stk::mesh::Entity face = stk::mesh::impl::get_face_for_element_side(bulkData, this_elem_entity, side_id);
                    if(bulkData.is_valid(face))
                    {
                        deletedEntities.push_back(face);
                    }
                }
            }
            else
            {
                // create or delete a face with a particular id

                // stk::mesh::EntityId local_id = bulkData.identifier(this_elem_entity);
                stk::mesh::EntityId remote_id = elementGraph.get_entity_id_of_remote_element(this_elem_entity, j);

                impl::parallel_info &parallel_edge_info = elementGraph.get_parallel_edge_info(this_elem_entity, remote_id);
                bool other_element_active = parallel_edge_info.m_in_part;
                bool create_face = false;
                if(other_element_active)
                {
                    create_face = true;
                }

                impl::create_or_delete_shared_face(bulkData, parallel_edge_info, elementGraph, this_elem_entity, remote_id, create_face, add_parts,
                        shared_modified, deletedEntities, id_counter, requestedIds[id_counter]);
            }
        }
    }

    ThrowRequireMsg(id_counter<requestedIds.size(), "Program error. Please contact sierra-help@sandia.gov for support.");
    elementGraph.set_num_face_ids_used(id_counter);
    stk::mesh::impl::delete_entities_and_upward_relations(bulkData, deletedEntities);
    bulkData.modification_end_for_face_creation_and_deletion(shared_modified, deletedEntities);
}

}} // end namespaces stk mesh

