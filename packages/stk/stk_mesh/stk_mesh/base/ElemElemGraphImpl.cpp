#include "ElemElemGraphImpl.hpp"
#include "ElemElemGraph.hpp"

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
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/environment/ReportHandler.hpp>
#include <stk_util/util/SortAndUnique.hpp>

namespace stk { namespace mesh {

namespace impl
{

void set_local_ids_and_fill_element_entities_and_topologies(stk::mesh::BulkData& bulkData, stk::mesh::EntityVector& local_id_to_element_entity, std::vector<stk::topology>& element_topologies)
{
    const stk::mesh::BucketVector& elemBuckets = bulkData.get_buckets(stk::topology::ELEM_RANK, bulkData.mesh_meta_data().locally_owned_part());
    size_t local_id = 0;
    for(size_t i=0; i<elemBuckets.size(); ++i)
    {
        const stk::mesh::Bucket& bucket = *elemBuckets[i];
        for(size_t j=0; j<bucket.size(); ++j)
        {
            local_id_to_element_entity[local_id] = bucket[j];
            element_topologies[local_id] = bucket.topology();
            bulkData.set_local_id(bucket[j], local_id);
            local_id++;
        }
    }
}

void fill_local_ids_and_fill_element_entities_and_topologies(stk::mesh::BulkData& bulkData, stk::mesh::EntityVector& local_id_to_element_entity,
                                                             std::vector<impl::LocalId>& entity_to_local_id, std::vector<stk::topology>& element_topologies)
{
    const stk::mesh::BucketVector & elemBuckets = bulkData.get_buckets(stk::topology::ELEM_RANK, bulkData.mesh_meta_data().locally_owned_part());
    size_t local_id = 0;
    for(size_t i=0; i<elemBuckets.size(); ++i)
    {
        const stk::mesh::Bucket& bucket = *elemBuckets[i];
        for(size_t j=0; j<bucket.size(); ++j)
        {
            local_id_to_element_entity[local_id] = bucket[j];
            element_topologies[local_id] = bucket.topology();
            entity_to_local_id[bucket[j].local_offset()] = local_id;
            local_id++;
        }
    }
}



ElemSideToProcAndFaceId get_element_side_ids_to_communicate(const stk::mesh::BulkData& bulkData, const stk::mesh::EntityVector &element_list)
{
    stk::mesh::EntityVector elements_to_communicate;

    for(const stk::mesh::Entity &element : element_list )
    {
        const stk::mesh::Entity* nodes = bulkData.begin_nodes(element);
        int numNodes = bulkData.num_nodes(element);

        for(int i=0; i<numNodes; ++i)
        {
            stk::mesh::Entity node = nodes[i];

            if(bulkData.bucket(node).shared())
            {
                elements_to_communicate.push_back(element);
                break;
            }
        }
    }

    return build_element_side_ids_to_proc_map(bulkData, elements_to_communicate);
}


ElemSideToProcAndFaceId build_element_side_ids_to_proc_map(const stk::mesh::BulkData& bulkData, const stk::mesh::EntityVector &elements_to_communicate)
{
    ElemSideToProcAndFaceId elem_side_comm;

    for(size_t i=0;i<elements_to_communicate.size();++i)
    {
        stk::mesh::Entity elem = elements_to_communicate[i];
        const stk::mesh::Entity* elem_nodes = bulkData.begin_nodes(elem);
        stk::topology elem_top = bulkData.bucket(elem).topology();
        unsigned num_sides = elem_top.num_sides();
        for(unsigned side=0;side<num_sides;++side)
        {
            stk::topology side_top = elem_top.side_topology(side);
            unsigned num_nodes_this_side = side_top.num_nodes();
            stk::mesh::EntityVector side_nodes(num_nodes_this_side);
            elem_top.side_nodes(elem_nodes, side, side_nodes.begin());
            std::vector<stk::mesh::EntityKey> keys(num_nodes_this_side);
            for(size_t j=0;j<keys.size();++j)
            {
                keys[j] = bulkData.entity_key(side_nodes[j]);
            }
            std::vector<int> sharing_procs;
            bulkData.shared_procs_intersection(keys, sharing_procs);
            for (int proc: sharing_procs) {
                elem_side_comm.insert(std::pair<EntitySidePair, ProcFaceIdPair>(EntitySidePair(elem, side), ProcFaceIdPair(proc,0)));
            }
        }
    }
    return elem_side_comm;
}

/*
size_t pack_shared_side_nodes_of_elements(stk::CommSparse& comm,
                                        const stk::mesh::BulkData& bulkData,
                                        ElemSideToProcAndFaceId &elements_to_communicate,
                                        const std::vector<stk::mesh::EntityId>& suggested_side_ids,
                                        const stk::mesh::Selector &sel, const stk::mesh::Selector* air)
{
    ElemSideToProcAndFaceId::iterator iter = elements_to_communicate.begin();
    ElemSideToProcAndFaceId::const_iterator end = elements_to_communicate.end();
    size_t counter = 0;

    for(; iter!= end; ++iter)
    {
        stk::mesh::Entity elem = iter->first.entity;
        unsigned side_index    = iter->first.side_id;
        int sharing_proc       = iter->second.proc;
        stk::mesh::EntityId element_id     = bulkData.identifier(elem);
        stk::mesh::EntityId suggested_side_id = suggested_side_ids[counter];
        ++counter;
        iter->second.side_id = suggested_side_id;

        stk::topology topology = bulkData.bucket(elem).topology();
        const bool isSelected = sel(bulkData.bucket(elem));
        const stk::mesh::Entity* elem_nodes = bulkData.begin_nodes(elem);

        unsigned num_nodes_this_side = topology.side_topology(side_index).num_nodes();
        stk::mesh::EntityVector side_nodes(num_nodes_this_side);
        topology.side_nodes(elem_nodes, side_index, side_nodes.begin());

        std::vector<stk::mesh::EntityKey> side_node_entity_keys(num_nodes_this_side);
        for(size_t i=0; i<num_nodes_this_side; ++i)
        {
            side_node_entity_keys[i] = bulkData.entity_key(side_nodes[i]);
        }

        comm.send_buffer(sharing_proc).pack<stk::mesh::EntityId>(element_id);
        comm.send_buffer(sharing_proc).pack<stk::topology>(topology);
        comm.send_buffer(sharing_proc).pack<unsigned>(side_index);
        comm.send_buffer(sharing_proc).pack<stk::mesh::EntityId>(suggested_side_id);
        comm.send_buffer(sharing_proc).pack<bool>(isSelected);
        if(air!=nullptr)
        {
            bool is_in_air = (*air)(bulkData.bucket(elem));
            comm.send_buffer(sharing_proc).pack<bool>(is_in_air);
        }

        comm.send_buffer(sharing_proc).pack<unsigned>(num_nodes_this_side);
        for(size_t i=0; i<num_nodes_this_side; ++i)
        {
            comm.send_buffer(sharing_proc).pack<stk::mesh::EntityKey>(side_node_entity_keys[i]);
        }
    }
    return counter;
}
*/

bool does_element_have_side(const stk::mesh::BulkData& bulkData, stk::mesh::Entity element)
{
    unsigned dimension_of_element = bulkData.bucket(element).topology().dimension();
    unsigned dimension_of_mesh = bulkData.mesh_meta_data().spatial_dimension();
    return dimension_of_element == dimension_of_mesh;
}

void pack_elements_to_comm(stk::CommSparse &comm, const std::vector<graphEdgeProc>& elements_to_comm)
{
    for(size_t i=0;i<elements_to_comm.size();++i)
    {
        int remote_proc = elements_to_comm[i].m_proc_id;
        comm.send_buffer(remote_proc).pack<stk::mesh::EntityId>(elements_to_comm[i].m_localElementId);
        comm.send_buffer(remote_proc).pack<int>(elements_to_comm[i].m_localSide);
        comm.send_buffer(remote_proc).pack<stk::mesh::EntityId>(elements_to_comm[i].m_remoteElementId);
        comm.send_buffer(remote_proc).pack<int>(elements_to_comm[i].m_remoteSide);
    }
}

std::vector<graphEdgeProc> communicate_killed_entities(stk::ParallelMachine communicator, const std::vector<graphEdgeProc>& elements_to_comm)
{
    std::vector<graphEdgeProc> remote_edges;
    stk::CommSparse comm(communicator);
    pack_elements_to_comm(comm, elements_to_comm);
    comm.allocate_buffers();
    pack_elements_to_comm(comm, elements_to_comm);
    comm.communicate();

    int num_procs = stk::parallel_machine_size(communicator);
    for(int i=0;i<num_procs;++i)
    {
        while(comm.recv_buffer(i).remaining())
        {
            stk::mesh::EntityId remoteId;
            int remoteSide;
            stk::mesh::EntityId localId;
            int localSide;
            comm.recv_buffer(i).unpack<stk::mesh::EntityId>(remoteId);
            comm.recv_buffer(i).unpack<int>(remoteSide);
            comm.recv_buffer(i).unpack<stk::mesh::EntityId>(localId);
            comm.recv_buffer(i).unpack<int>(localSide);
            remote_edges.push_back(graphEdgeProc(localId, localSide, remoteId, remoteSide, i));
        }
    }
    return remote_edges;
}

int get_element_side_multiplier()
{
    return 1000;
}

stk::mesh::EntityId get_side_id_for_remotely_connected_element(stk::mesh::EntityId local_element_id,
        stk::mesh::EntityId local_side_id, stk::mesh::EntityId remote_element_id,
        stk::mesh::EntityId remote_side_id, size_t& id_counter)
{
    stk::mesh::EntityId side_global_id = 0;
    if(local_element_id < remote_element_id)
    {
        side_global_id = local_side_id;
        id_counter++;
    }
    else
    {
        side_global_id = remote_side_id;
    }

    return side_global_id;
}

stk::mesh::Permutation get_permutation_for_new_side(const parallel_info& parallel_edge_info, stk::mesh::EntityId local_element_id, stk::mesh::EntityId remote_element_id)
{
    stk::mesh::Permutation perm;
    if(local_element_id < remote_element_id)
    {
        perm = static_cast<stk::mesh::Permutation>(0);
    }
    else
    {
        perm = static_cast<stk::mesh::Permutation>(parallel_edge_info.m_permutation);
    }
    return perm;
}

bool is_id_already_in_use_locally(stk::mesh::BulkData& bulkData, stk::mesh::EntityRank rank, stk::mesh::EntityId id)
{
    stk::mesh::Entity entity = bulkData.get_entity(rank, id);
    return bulkData.is_valid(entity);
}

bool does_side_exist_with_different_permutation(stk::mesh::BulkData& bulkData, stk::mesh::Entity element,
        stk::mesh::ConnectivityOrdinal side_ordinal, stk::mesh::Permutation perm)
{
    stk::mesh::EntityRank side_rank = bulkData.mesh_meta_data().side_rank();

    unsigned elem_num_sides = bulkData.num_connectivity(element, side_rank);
    const stk::mesh::ConnectivityOrdinal * elem_ord_it = bulkData.begin_ordinals(element, side_rank);
    const stk::mesh::Permutation * elem_perm_it = bulkData.begin_permutations(element, side_rank);

    for (unsigned i=0 ; i<elem_num_sides ; ++i)
    {
        if (elem_ord_it[i] == static_cast<unsigned>(side_ordinal))
        {
            if (perm != elem_perm_it[i])
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    }
    return false;
}

bool does_element_side_exist(stk::mesh::BulkData& bulkData, stk::mesh::Entity element, stk::mesh::ConnectivityOrdinal side_ordinal)
{
    stk::mesh::Entity side = stk::mesh::Entity();
    stk::mesh::EntityRank side_rank = bulkData.mesh_meta_data().side_rank();

    unsigned elem_num_sides = bulkData.num_connectivity(element, side_rank);
    const stk::mesh::Entity * elem_sides = bulkData.begin(element, side_rank);
    const stk::mesh::ConnectivityOrdinal * elem_ord_it = bulkData.begin_ordinals(element, side_rank);
    for (unsigned i=0 ; i<elem_num_sides ; ++i)
    {
        if (elem_ord_it[i] == static_cast<unsigned>(side_ordinal))
        {
            side = elem_sides[i];
            break;
        }
    }

    return bulkData.is_valid(side);
}

stk::mesh::Entity get_side_for_element(const stk::mesh::BulkData& bulkData, stk::mesh::Entity element, int side_id)
{
    stk::mesh::EntityRank side_rank = bulkData.mesh_meta_data().side_rank();

    unsigned num_sides = bulkData.num_connectivity(element, side_rank);
    const stk::mesh::Entity *sides = bulkData.begin(element, side_rank);
    const stk::mesh::ConnectivityOrdinal *ordinals = bulkData.begin_ordinals(element, side_rank);

    stk::mesh::Entity side;

    for(unsigned ii = 0; ii < num_sides; ++ii)
    {
        if(ordinals[ii] == static_cast<stk::mesh::ConnectivityOrdinal>(side_id))
        {
            side = sides[ii];
            break;
        }
    }
    return side;
}

stk::mesh::PartVector get_stk_parts_for_moving_parts_into_death_boundary(const stk::mesh::PartVector *bc_mesh_parts)
{
    stk::mesh::PartVector sideParts;
    if(bc_mesh_parts != nullptr)
    {
        const stk::mesh::PartVector * meshparts_to_apply = bc_mesh_parts;
        unsigned int number_of_meshparts = meshparts_to_apply->size();

        for(unsigned int index = 0; index < number_of_meshparts; ++index)
        {
            stk::mesh::Part * mp = (*meshparts_to_apply)[index];

            sideParts.push_back(mp);

            stk::mesh::PartVector::const_iterator isup = mp->supersets().begin();
            for(; isup != mp->supersets().end(); ++isup)
            {
                if(!stk::mesh::is_auto_declared_part(**isup))
                {
                    sideParts.push_back(*isup);
                }
            }
        }
    }
    return sideParts;
}

stk::mesh::Part* get_sides_created_during_death_part(const stk::mesh::MetaData &metaData)
{
    return metaData.get_part("sides_created_during_death");
}

void create_sides_created_during_death_part(stk::mesh::MetaData &metaData)
{
    stk::mesh::EntityRank side_rank = metaData.side_rank();
    const bool forceNoInduce = true;
    metaData.declare_part("sides_created_during_death", side_rank, forceNoInduce);
}

void add_parts_from_element(stk::mesh::BulkData& bulkData, stk::mesh::Entity element, stk::mesh::PartVector& side_parts)
{
    const stk::mesh::PartVector & supersets = bulkData.bucket(element).supersets();
    for (size_t part_i=0 ; part_i<supersets.size() ; ++part_i)
    {
        if(!stk::mesh::is_auto_declared_part(*supersets[part_i]))
        {
            side_parts.push_back(supersets[part_i]);
        }
    }
}

stk::mesh::PartVector get_parts_for_creating_side(stk::mesh::BulkData& bulkData, const stk::mesh::PartVector& parts_for_creating_side, stk::mesh::Entity element, int side_ord)
{
    stk::mesh::PartVector side_parts = parts_for_creating_side;
    add_parts_from_element(bulkData, element, side_parts);
    stk::topology side_top = bulkData.bucket(element).topology().side_topology(side_ord);
    side_parts.push_back(&bulkData.mesh_meta_data().get_topology_root_part(side_top));
    side_parts.push_back(get_sides_created_during_death_part(bulkData.mesh_meta_data()));

    return side_parts;
}

void add_side_into_exposed_boundary(stk::mesh::BulkData& bulkData, const parallel_info& parallel_edge_info,
        stk::mesh::Entity local_element, int side_id, stk::mesh::EntityId remote_id, const stk::mesh::PartVector& parts_for_creating_side,
        std::vector<stk::mesh::sharing_info> &shared_modified, const stk::mesh::PartVector *boundary_mesh_parts)
{
    stk::mesh::EntityId side_global_id = parallel_edge_info.m_chosen_side_id;
    stk::mesh::ConnectivityOrdinal side_ord = static_cast<stk::mesh::ConnectivityOrdinal>(side_id);
    std::string msg = "Program error. Contact sierra-help@sandia.gov for support.";

    // determine which element is active
    stk::mesh::Permutation perm = stk::mesh::DEFAULT_PERMUTATION;
    int other_proc = parallel_edge_info.m_other_proc;
    int owning_proc = std::min(other_proc, bulkData.parallel_rank());

    if(parallel_edge_info.m_in_body_to_be_skinned)
    {
        perm = static_cast<stk::mesh::Permutation>(parallel_edge_info.m_permutation);
    }

    stk::mesh::Entity side = stk::mesh::impl::get_side_for_element(bulkData, local_element, side_id);

    if(!bulkData.is_valid(side))
    {
        stk::mesh::PartVector side_parts = get_parts_for_creating_side(bulkData, parts_for_creating_side, local_element, side_id);
        ThrowRequireMsg(!impl::is_id_already_in_use_locally(bulkData, bulkData.mesh_meta_data().side_rank(), side_global_id), msg);
        side = connect_side_to_element(bulkData, local_element, side_global_id, side_ord, perm, side_parts);
        shared_modified.push_back(stk::mesh::sharing_info(side, other_proc, owning_proc));
    }
    else
    {
        if(bulkData.bucket(side).owned())
        {
            stk::mesh::PartVector parts = get_stk_parts_for_moving_parts_into_death_boundary(boundary_mesh_parts);
            bulkData.change_entity_parts(side, parts, stk::mesh::PartVector());
            shared_modified.push_back(stk::mesh::sharing_info(side, other_proc, bulkData.parallel_owner_rank(side)));
        }
    }
}

bool side_created_during_death(stk::mesh::BulkData& bulkData, stk::mesh::Entity side)
{
    stk::mesh::Part& sides_created_during_death = *get_sides_created_during_death_part(bulkData.mesh_meta_data());
    return bulkData.is_valid(side) && bulkData.bucket(side).member(sides_created_during_death);
}

void remove_side_from_death_boundary(stk::mesh::BulkData& bulkData, stk::mesh::Entity local_element,
        stk::mesh::Part &activePart, stk::mesh::EntityVector &deletedEntities, int side_id)
{
    stk::mesh::Entity side = stk::mesh::impl::get_side_for_element(bulkData, local_element, side_id);
    if(side_created_during_death(bulkData, side))
    {
        deletedEntities.push_back(side);
    }
    else if(bulkData.is_valid(side) && bulkData.bucket(side).owned())
    {
        bulkData.change_entity_parts(side, {}, {&activePart});
    }
}


stk::mesh::Entity connect_side_to_element(stk::mesh::BulkData& bulkData, stk::mesh::Entity element,
        stk::mesh::EntityId side_global_id, stk::mesh::ConnectivityOrdinal side_ordinal,
        stk::mesh::Permutation side_permutation, const stk::mesh::PartVector& parts)
{
    stk::mesh::EntityRank side_rank = bulkData.mesh_meta_data().side_rank();

    stk::mesh::Entity side = bulkData.declare_entity(side_rank, side_global_id, parts);

    // connect element to side
    bulkData.declare_relation(element, side, side_ordinal, side_permutation);

    // connect side to nodes
    const stk::mesh::Entity* nodes = bulkData.begin_nodes(element);
    stk::topology side_top = bulkData.bucket(element).topology().side_topology(side_ordinal);
    stk::mesh::EntityVector side_nodes(side_top.num_nodes());
    bulkData.bucket(element).topology().side_nodes(nodes, side_ordinal, side_nodes.begin());
    stk::mesh::EntityVector permuted_side_nodes(side_top.num_nodes());
    side_top.permutation_nodes(side_nodes, side_permutation, permuted_side_nodes.begin());
    for(size_t i=0;i<permuted_side_nodes.size();++i)
    {
        bulkData.declare_relation(side, permuted_side_nodes[i], i);
    }

    return side;
}

void add_shell_element_if_coincident(const stk::mesh::BulkData& mesh, const unsigned sideOrdinal, const stk::mesh::Entity localElement, const ConnectedElementData& connectedElem, ConnectedElementDataVector& filteredConnectedElements)
{
    const stk::mesh::EntityVector &sideNodesOfReceivedElement = connectedElem.m_sideNodes;
    stk::mesh::OrdinalAndPermutation localElemOrdAndPerm =
            stk::mesh::get_ordinal_and_permutation(mesh, localElement, mesh.mesh_meta_data().side_rank(), sideNodesOfReceivedElement);
    // for shell element, want the nodes of the solid to be in opposite order. So getting non-matching side ordinals
    // means the normals oppose
    bool does_local_shell_side_normal_oppose_other_element_side_normal = (localElemOrdAndPerm.first == sideOrdinal);

    if (does_local_shell_side_normal_oppose_other_element_side_normal)
    {
        filteredConnectedElements.push_back(connectedElem);
    }
}

void add_solid_element_if_normals_oppose_to_shell(const stk::mesh::BulkData& mesh, const unsigned sideOrdinal, const stk::mesh::Entity localElement, const ConnectedElementData& connectedElem, ConnectedElementDataVector& filteredConnectedElements)
{
    const stk::mesh::EntityVector &sideNodesOfReceivedElement = connectedElem.m_sideNodes;
    stk::mesh::OrdinalAndPermutation localElemOrdAndPerm =
            stk::mesh::get_ordinal_and_permutation(mesh, localElement, mesh.mesh_meta_data().side_rank(), sideNodesOfReceivedElement);
    // for shell element, want the nodes of the solid to be in opposite order. So getting non-matching side ordinals
    // means the normals oppose
    bool does_local_shell_side_normal_oppose_other_element_side_normal = (localElemOrdAndPerm.first != sideOrdinal);

    if (does_local_shell_side_normal_oppose_other_element_side_normal)
    {
        filteredConnectedElements.push_back(connectedElem);
    }
}

void add_shell_connections_to_this_solid_if_normals_oppose(const stk::mesh::BulkData& mesh, const stk::mesh::Entity localElement, const unsigned sideOrdinal, ConnectedElementDataVector& connectedElementData)
{
    ConnectedElementDataVector filteredConnectedElements;

    stk::topology localElemTopology = mesh.bucket(localElement).topology();
    stk::topology localSideTopology = localElemTopology.side_topology(sideOrdinal);
    bool foundAnySingleElementThatIsEquivalentToLocalElement = false;
    const stk::mesh::Entity* localElemNodes = mesh.begin_nodes(localElement);
    stk::mesh::EntityVector localElemSideNodes;
    localElemSideNodes.resize(localSideTopology.num_nodes());
    localElemTopology.side_nodes(localElemNodes, sideOrdinal, localElemSideNodes.begin());
    
    for (const ConnectedElementData & connectedElem: connectedElementData)
    {
        std::pair<bool,unsigned> result = localSideTopology.equivalent(localElemSideNodes, connectedElem.m_sideNodes);
        const bool isEquivalentNodes = result.first;
        foundAnySingleElementThatIsEquivalentToLocalElement = foundAnySingleElementThatIsEquivalentToLocalElement || isEquivalentNodes;
    
        if (connectedElem.m_elementTopology.is_shell() && isEquivalentNodes)
        {
            stk::mesh::OrdinalAndPermutation localElemOrdAndPerm = stk::mesh::get_ordinal_and_permutation(mesh, localElement, mesh.mesh_meta_data().side_rank(), connectedElem.m_sideNodes);
            bool localNegativeRelativeFacePolarity = !localSideTopology.is_positive_polarity(localElemOrdAndPerm.second);
    
            if (localNegativeRelativeFacePolarity)
            {
                filteredConnectedElements.push_back(connectedElem);
            }
        }
    }
    
    if (!filteredConnectedElements.empty()) 
    {
        connectedElementData.swap(filteredConnectedElements);
    }

    if (!foundAnySingleElementThatIsEquivalentToLocalElement)
    {
        connectedElementData.clear();
    }
}

void filter_out_invalid_solid_shell_connections(const stk::mesh::BulkData & mesh,
                                              const stk::mesh::Entity localElement,
                                              const unsigned sideOrdinal,
                                              ConnectedElementDataVector & connectedElementData)
{
    stk::topology localElemTopology = mesh.bucket(localElement).topology();

    if (localElemTopology.is_shell())
    {
        ConnectedElementDataVector filteredConnectedElements;
        for (const ConnectedElementData & connectedElem: connectedElementData)
        {
            if(mesh.identifier(localElement) != connectedElem.m_elementIdentifier)
            {
                if(connectedElem.m_elementTopology.is_shell())
                    add_shell_element_if_coincident(mesh, sideOrdinal, localElement, connectedElem, filteredConnectedElements);
                else
                    add_solid_element_if_normals_oppose_to_shell(mesh, sideOrdinal, localElement, connectedElem, filteredConnectedElements);
            }
        }
        connectedElementData.swap(filteredConnectedElements);
    }
    else
    {
        add_shell_connections_to_this_solid_if_normals_oppose(mesh, localElement, sideOrdinal, connectedElementData);
    }
}

void pack_newly_shared_remote_edges(stk::CommSparse &comm, const stk::mesh::BulkData &bulkData, const std::vector<SharedEdgeInfo> &newlySharedEdges)
{
    std::vector<SharedEdgeInfo>::const_iterator iter = newlySharedEdges.begin();
    std::vector<SharedEdgeInfo>::const_iterator endIter = newlySharedEdges.end();

    for(; iter!= endIter; ++iter)
    {
        stk::mesh::EntityId localId = iter->m_locaElementlId;
        stk::mesh::Entity localEntity = bulkData.get_entity(stk::topology::ELEM_RANK, localId);
        stk::mesh::EntityId remoteId = iter->m_remoteElementId;
        unsigned side_index    = iter->m_sideIndex;
        int sharing_proc       = iter->m_procId;
        stk::mesh::EntityId chosenId = iter->m_chosenSideId;
        const bool isInPart = iter->m_isInPart;
        const bool isAir    = iter->m_isInAir;

        size_t numNodes= iter->m_sharedNodes.size();
        std::vector<stk::mesh::EntityKey> side_node_entity_keys(numNodes);
        for(size_t i=0; i<numNodes; ++i)
        {
            side_node_entity_keys[i] = bulkData.entity_key(iter->m_sharedNodes[i]);
        }

        comm.send_buffer(sharing_proc).pack<stk::mesh::EntityId>(localId);
        comm.send_buffer(sharing_proc).pack<stk::mesh::EntityId>(remoteId);
        comm.send_buffer(sharing_proc).pack<unsigned>(side_index);
        comm.send_buffer(sharing_proc).pack<stk::mesh::EntityId>(chosenId);
        comm.send_buffer(sharing_proc).pack<bool>(isInPart);
        comm.send_buffer(sharing_proc).pack<bool>(isAir);
        comm.send_buffer(sharing_proc).pack<stk::topology>(bulkData.bucket(localEntity).topology());
        comm.send_buffer(sharing_proc).pack<unsigned>(numNodes);
        for(size_t i=0; i<numNodes; ++i)
        {
            comm.send_buffer(sharing_proc).pack<stk::mesh::EntityKey>(side_node_entity_keys[i]);
        }
    }
}

int count_shared_sides(const stk::mesh::Graph &graph, stk::mesh::impl::LocalId elem1,  stk::mesh::impl::LocalId elem2)
{
    int numSharedSides = 0;
    for(size_t i=0; i < graph.get_num_edges_for_element(elem1); i++)
    {
        const stk::mesh::GraphEdge &graphEdge = graph.get_edge_for_element(elem1, i);
        if(graphEdge.elem2 == elem2)
            numSharedSides++;
    }
    return numSharedSides;
}

bool are_elements_coincident(const stk::mesh::Graph &graph,
                             const CoincidentElementDescription &elemDesc)
{
    int numSharedSides = count_shared_sides(graph, elemDesc.elem1, elemDesc.elem2);
    return (numSharedSides == elemDesc.numSides);
}

void add_exposed_sides(LocalId elementId, size_t numElemSides,
                      const stk::mesh::Graph &graph, std::vector<int> &element_side_pairs)
{
    if (graph.get_num_edges_for_element(elementId) < numElemSides)
    {
        std::vector<int> elemSides(numElemSides, -1);
        for(size_t j = 0; j < graph.get_num_edges_for_element(elementId); ++j)
        {
            const stk::mesh::GraphEdge & graphEdge = graph.get_edge_for_element(elementId, j);
            int sideId = graphEdge.side1;
            elemSides[sideId] = sideId;
        }

        for(size_t sideId = 0; sideId < numElemSides; ++sideId)
            if (elemSides[sideId] == -1)
                element_side_pairs.push_back(sideId);
    }
}

void fill_coincident_edges(stk::mesh::Graph &graph, const std::vector<stk::topology> &topologies, size_t elemId, std::vector<stk::mesh::GraphEdge> &edgesToDelete)
{
    for(size_t edgeIndex = 0; edgeIndex < graph.get_num_edges_for_element(elemId); edgeIndex++)
    {
        const stk::mesh::GraphEdge &graphEdge = graph.get_edge_for_element(elemId, edgeIndex);
        if(stk::mesh::impl::are_elements_coincident(graph, {static_cast<int>(topologies[elemId].num_sides()), graphEdge.elem1, graphEdge.elem2}))
            edgesToDelete.push_back(graphEdge);
    }
}

//void delete_symmetric_connections(stk::mesh::Graph& graph, const stk::mesh::impl::LocalId elemId)
//{
//    for(size_t i = 0; i < graph.get_num_edges_for_element(elemId); i++)
//    {
//        const stk::mesh::GraphEdge& otherEdge = graph.get_edge_for_element(elemId, i);
//        if(is_local_element(otherEdge.elem2))
//            graph.delete_edge(create_symmetric_edge(otherEdge));
//    }
//}
//
//void delete_all_edges_involving_element(stk::mesh::Graph& graph, const stk::mesh::impl::LocalId elemId)
//{
//    delete_symmetric_connections(graph, elemId);
//    graph.delete_all_connections(elemId);
//    graph.invalidate(elemId);
//}

void delete_edges(stk::mesh::Graph& graph, const std::vector<stk::mesh::GraphEdge>& edgesToDelete)
{
    for(const stk::mesh::GraphEdge& edgeToDelete : edgesToDelete)
        graph.delete_edge(edgeToDelete);

//    for(const stk::mesh::GraphEdge& edgeToDelete : edgesToDelete)
//        if(is_local_element(edgeToDelete.elem2))
//            delete_all_edges_involving_element(graph, edgeToDelete.elem2);
}

void add_edges(const std::vector<stk::mesh::GraphEdge>& edgesToDelete, SparseGraph& extractedCoincidentSides)
{
    for(const stk::mesh::GraphEdge& edgeToDelete : edgesToDelete)
        extractedCoincidentSides[edgeToDelete.elem1].push_back(edgeToDelete);
}

void extract_coincident_sides_for_element(stk::mesh::Graph& graph,
                                          const std::vector<stk::topology>& topologies,
                                          size_t elemId,
                                          SparseGraph& extractedCoincidentSides)
{
    std::vector<stk::mesh::GraphEdge> coincidentEdges;
    fill_coincident_edges(graph, topologies, elemId, coincidentEdges);
    add_edges(coincidentEdges, extractedCoincidentSides);
    delete_edges(graph, coincidentEdges);
}

SparseGraph extract_coincident_sides(stk::mesh::Graph &graph, const std::vector<stk::topology> &topologies)
{
    SparseGraph extractedCoincidentSides;
    for(size_t elemId = 0; elemId < graph.get_num_elements_in_graph(); elemId++)
        extract_coincident_sides_for_element(graph, topologies, elemId, extractedCoincidentSides);
    return extractedCoincidentSides;
}

bool is_local_element(stk::mesh::impl::LocalId elemId)
{
    return (elemId >= 0);
}

std::map<int, std::vector<LocalId>> get_extracted_coincident_local_ids(const stk::mesh::Graph &graph,
                                                                       const stk::mesh::ParallelInfoForGraphEdges &parallelInfoForGraphEdges,
                                                                       const stk::mesh::impl::SparseGraph &extractedCoincidentElements)
{
    std::map<int, std::vector<LocalId>> extractedIdsAndConnectedProcs;
    for(const stk::mesh::impl::SparseGraph::value_type &extractedElementGraphEdges : extractedCoincidentElements)
    {
        std::vector<LocalId> extractedIds;
        for(const stk::mesh::GraphEdge &graphEdge : extractedElementGraphEdges.second)
        {
            if(is_local_element(graphEdge.elem2) && !graph.is_valid(graphEdge.elem2))
            {
                extractedIds.push_back(graphEdge.elem2);
            }
        }
        stk::util::sort_and_unique(extractedIds);
        if(!extractedIds.empty())
        {
            for(const stk::mesh::GraphEdge &graphEdge : graph.get_edges_for_element(extractedElementGraphEdges.first))
            {
                if(!is_local_element(graphEdge.elem2))
                {
                    int otherProc = parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(graphEdge).m_other_proc;
                    std::vector<LocalId> &ids = extractedIdsAndConnectedProcs[otherProc];
                    ids.insert(ids.end(), extractedIds.begin(), extractedIds.end());
                }
            }
        }
    }
    return extractedIdsAndConnectedProcs;
}

void remove_edges_to_extracted_coincident_elements_on_other_procs(const std::map<int, stk::mesh::EntityIdVector> &extractedEntityIdsByProc,
                                                                  stk::mesh::Graph &graph,
                                                                  MPI_Comm comm)
{
    stk::CommSparse commSparse(comm);
    stk::pack_and_communicate(commSparse,
        [&commSparse, &extractedEntityIdsByProc]()
        {
            stk::mesh::impl::pack_extracted_coincident_element_ids(commSparse, extractedEntityIdsByProc);
        });
    stk::unpack_communications(commSparse,
        [&commSparse,&graph](int procId)
        {
            stk::mesh::impl::unpack_extracted_coincident_element_ids_and_delete_connecting_graph_edges(commSparse, procId, graph);
        });
}

void pack_extracted_coincident_element_ids(stk::CommSparse &commSparse, const std::map<int, stk::mesh::EntityIdVector> &extractedEntityIdsByProc)
{
    for(const auto &procAndEntityIdVec : extractedEntityIdsByProc)
    {
        int proc = procAndEntityIdVec.first;
        const stk::mesh::EntityIdVector &idVec = procAndEntityIdVec.second;
        commSparse.send_buffer(proc).pack<size_t>(idVec.size());
        for(const stk::mesh::EntityId extractedId : idVec)
        {
            commSparse.send_buffer(proc).pack<stk::mesh::EntityId>(extractedId);
        }
    }
}

void delete_graph_edges_connected_to_coincident_element(stk::mesh::Graph &graph, std::vector<stk::mesh::GraphEdge> graphEdgesToDelete)
{
    for(const stk::mesh::GraphEdge &graphEdge : graphEdgesToDelete)
        graph.delete_edge(graphEdge);
}

size_t receive_number_of_coincident_elements(stk::CommSparse &commSparse, int procId)
{
    size_t numCoincidentElems;
    commSparse.recv_buffer(procId).unpack<size_t>(numCoincidentElems);
    return numCoincidentElems;
}

stk::mesh::EntityId unpack_next_coincident_element_id(stk::CommSparse &commSparse, int procId)
{
    stk::mesh::EntityId coincidentElem;
    commSparse.recv_buffer(procId).unpack<stk::mesh::EntityId>(coincidentElem);
    return coincidentElem;
}

bool is_element_remote(stk::mesh::EntityId connectedElementId)
{
    return (!stk::mesh::impl::is_local_element(connectedElementId));
}
bool is_element_coincident(stk::mesh::EntityId coincidentElemId, stk::mesh::EntityId connectedElementId)
{
    return (static_cast<stk::mesh::EntityId>(-connectedElementId) == coincidentElemId);
}
bool is_connected_element_coincident_and_remote(stk::mesh::EntityId coincidentElemId, stk::mesh::EntityId connectedElementId)
{
    return (is_element_coincident(coincidentElemId, connectedElementId) && is_element_remote(connectedElementId));
}

void collect_remote_coincident_graph_edges(stk::mesh::EntityId coincidentElemId, const GraphEdgesForElement &graphEdgesOfConcern, std::vector<stk::mesh::GraphEdge> &graphEdgesToDelete)
{
    for(const stk::mesh::GraphEdge &graphEdge : graphEdgesOfConcern)
        if(is_connected_element_coincident_and_remote(coincidentElemId, graphEdge.elem2))
            graphEdgesToDelete.push_back(graphEdge);
}

void delete_graph_edges_connected_to_coincident_element(stk::mesh::Graph &graph, stk::mesh::EntityId coincidentElemId)
{
    std::vector<stk::mesh::GraphEdge> graphEdgesToDelete;
    for(size_t i=0; i<graph.get_num_elements_in_graph(); i++)
        collect_remote_coincident_graph_edges(coincidentElemId, graph.get_edges_for_element(i), graphEdgesToDelete);
    delete_graph_edges_connected_to_coincident_element(graph, graphEdgesToDelete);
}

void unpack_extracted_coincident_element_ids_and_delete_connecting_graph_edges(stk::CommSparse &commSparse, int procId, stk::mesh::Graph &graph)
{
    size_t numCoincidentElems = receive_number_of_coincident_elements(commSparse, procId);
    for(size_t i=0; i<numCoincidentElems; i++)
    {
        stk::mesh::EntityId coincidentElemId = unpack_next_coincident_element_id(commSparse, procId);
        delete_graph_edges_connected_to_coincident_element(graph, coincidentElemId);
    }
}


typedef std::map<stk::mesh::GraphEdge, stk::mesh::impl::parallel_info *, GraphEdgeLessByElem2> GraphEdgeToParInfoMap;

void pack_graph_edge_and_chosen_side_id(stk::CommSparse& commSparse,
                                        const stk::mesh::GraphEdge& graphEdge,
                                        const parallel_info* parInfo,
                                        const IdMapper& idMapper)
{
    stk::mesh::EntityId elem1GlobalId = idMapper.localToGlobal(graphEdge.elem1);
    commSparse.send_buffer(parInfo->m_other_proc).pack<stk::mesh::EntityId>(-graphEdge.elem2);
    commSparse.send_buffer(parInfo->m_other_proc).pack<int>(graphEdge.side2);
    commSparse.send_buffer(parInfo->m_other_proc).pack<stk::mesh::EntityId>(elem1GlobalId);
    commSparse.send_buffer(parInfo->m_other_proc).pack<int>(graphEdge.side1);
    commSparse.send_buffer(parInfo->m_other_proc).pack<stk::mesh::EntityId>(parInfo->m_chosen_side_id);
}

void pack_changed_edges_and_chosen_side_ids(stk::CommSparse& commSparse,
                                            const GraphEdgeToParInfoMap& changedEdgesAndParInfos,
                                            const IdMapper& idMapper)
{
    for(const GraphEdgeToParInfoMap::value_type& changedEdgeAndParInfo : changedEdgesAndParInfos)
    {
        const stk::mesh::GraphEdge& graphEdge = changedEdgeAndParInfo.first;
        const parallel_info* parInfo = changedEdgeAndParInfo.second;
        pack_graph_edge_and_chosen_side_id(commSparse, graphEdge, parInfo, idMapper);
    }
}

void unpack_local_elem_side(stk::CommSparse& commSparse,
                            int procId,
                            const IdMapper& idMapper,
                            stk::mesh::GraphEdge& recvGraphEdge)
{
    stk::mesh::EntityId globalId;
    commSparse.recv_buffer(procId).unpack<stk::mesh::EntityId>(globalId);
    recvGraphEdge.elem1 = idMapper.globalToLocal(globalId);
    commSparse.recv_buffer(procId).unpack<int>(recvGraphEdge.side1);
}

void unpack_remote_elem_side(stk::CommSparse& commSparse, int procId, stk::mesh::GraphEdge& recvGraphEdge)
{
    stk::mesh::EntityId globalId;
    commSparse.recv_buffer(procId).unpack<stk::mesh::EntityId>(globalId);
    recvGraphEdge.elem2 = -globalId;
    commSparse.recv_buffer(procId).unpack<int>(recvGraphEdge.side2);
}

stk::mesh::GraphEdge unpack_graph_edge(stk::CommSparse& commSparse, int procId, const IdMapper& idMapper)
{
    stk::mesh::GraphEdge recvGraphEdge;
    unpack_local_elem_side(commSparse, procId, idMapper, recvGraphEdge);
    unpack_remote_elem_side(commSparse, procId, recvGraphEdge);
    return recvGraphEdge;
}

void update_chosen_side_id_for_graph_edge(stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges,
                                          const stk::mesh::GraphEdge& recvGraphEdge,
                                          stk::mesh::EntityId chosenSideId)
{
    stk::mesh::impl::parallel_info& parInfoToChange = parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(recvGraphEdge);
    parInfoToChange.m_chosen_side_id = chosenSideId;
}

stk::mesh::EntityId unpack_chosen_side_id(stk::CommSparse &commSparse, int procId)
{
    stk::mesh::EntityId chosenSideId;
    commSparse.recv_buffer(procId).unpack<stk::mesh::EntityId>(chosenSideId);
    return chosenSideId;
}

void unpack_changed_edges_and_chosen_side_ids(stk::CommSparse &commSparse,
                                              int procId,
                                              const IdMapper& idMapper,
                                              stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges)
{
    stk::mesh::GraphEdge recvGraphEdge = unpack_graph_edge(commSparse, procId, idMapper);
    stk::mesh::EntityId chosenSideId = unpack_chosen_side_id(commSparse, procId);
    update_chosen_side_id_for_graph_edge(parallelInfoForGraphEdges, recvGraphEdge, chosenSideId);
}

void fix_chosen_side_ids_on_other_procs(const GraphEdgeToParInfoMap &changedEdgesAndParInfos,
                                        stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges,
                                        const IdMapper& idMapper,
                                        MPI_Comm comm)
{
    stk::CommSparse commSparse(comm);
    stk::pack_and_communicate(commSparse,
        [&commSparse, &changedEdgesAndParInfos, &idMapper]()
        {
            pack_changed_edges_and_chosen_side_ids(commSparse, changedEdgesAndParInfos, idMapper);
        });
    stk::unpack_communications(commSparse,
       [&commSparse, &idMapper, &parallelInfoForGraphEdges](int procId)
       {
            unpack_changed_edges_and_chosen_side_ids(commSparse, procId, idMapper, parallelInfoForGraphEdges);
       });
}

std::set<stk::mesh::impl::LocalId> get_coincident_elems(const std::vector<stk::mesh::GraphEdge>& coincidentEdgesForElem)
{
    std::set<stk::mesh::impl::LocalId> coincidentElems;
    for(const stk::mesh::GraphEdge& graphEdge : coincidentEdgesForElem)
        coincidentElems.insert(graphEdge.elem2);
    return coincidentElems;
}

void add_graph_edge_and_par_info(const stk::mesh::GraphEdge& parallelGraphEdge,
                                 stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges,
                                 GraphEdgeToParInfoMap& parInfos)
{
    stk::mesh::impl::parallel_info& parInfoThisElem =
            parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(parallelGraphEdge);
    parInfos.insert(std::make_pair(parallelGraphEdge, &parInfoThisElem));
}

bool are_graph_edges_for_same_side(const stk::mesh::GraphEdge& parallelGraphEdge, const stk::mesh::GraphEdge& coincidentGraphEdge)
{
    return parallelGraphEdge.elem1 == coincidentGraphEdge.elem1 && parallelGraphEdge.side1 == coincidentGraphEdge.side1;
}

void add_coincident_graph_edges_and_par_infos(const stk::mesh::GraphEdge& parallelGraphEdge,
                                              const std::vector<stk::mesh::GraphEdge>& coincidentEdgesForElem,
                                              stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges,
                                              GraphEdgeToParInfoMap& parInfos)
{
    for(const stk::mesh::GraphEdge& coincidentGraphEdge : coincidentEdgesForElem)
    {
        if(are_graph_edges_for_same_side(parallelGraphEdge, coincidentGraphEdge) &&
                impl::is_local_element(coincidentGraphEdge.elem2))
        {
            stk::mesh::GraphEdge otherElemGraphEdge(coincidentGraphEdge.elem2,
                                                    coincidentGraphEdge.side2,
                                                    parallelGraphEdge.elem2,
                                                    parallelGraphEdge.side2);
            add_graph_edge_and_par_info(otherElemGraphEdge, parallelInfoForGraphEdges, parInfos);
        }
    }
}

GraphEdgeToParInfoMap get_coincident_parallel_infos(const stk::mesh::GraphEdge& parallelGraphEdge,
                                                    const std::vector<stk::mesh::GraphEdge> &coincidentEdgesForElem,
                                                    stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges)
{
    GraphEdgeToParInfoMap parInfos;
    add_graph_edge_and_par_info(parallelGraphEdge, parallelInfoForGraphEdges, parInfos);
    add_coincident_graph_edges_and_par_infos(parallelGraphEdge, coincidentEdgesForElem, parallelInfoForGraphEdges, parInfos);
    return parInfos;
}

stk::mesh::EntityId get_min_chosen_side_id(const GraphEdgeToParInfoMap &coincidentParInfos)
{

    stk::mesh::EntityId chosenSideId = std::numeric_limits<stk::mesh::EntityId>::max();
    for(const GraphEdgeToParInfoMap::value_type &graphEdgeAndParInfo : coincidentParInfos)
        chosenSideId = std::min(chosenSideId, graphEdgeAndParInfo.second->m_chosen_side_id);
    return chosenSideId;
}

void update_chosen_side_id_for_coincident_graph_edges(const GraphEdgeToParInfoMap &coincidentParInfos,
                                                      stk::mesh::EntityId chosenSideId,
                                                      GraphEdgeToParInfoMap& changedEdgesAndProcs)
{
    for(const GraphEdgeToParInfoMap::value_type &graphEdgeAndParInfo : coincidentParInfos)
    {
        if(graphEdgeAndParInfo.second->m_chosen_side_id != chosenSideId)
        {
            changedEdgesAndProcs.insert(graphEdgeAndParInfo);
            graphEdgeAndParInfo.second->m_chosen_side_id = chosenSideId;
        }
    }
}

void update_chosen_side_ids_for_remote_graph_edge(const stk::mesh::GraphEdge& graphEdge,
                                                  stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges,
                                                  const std::vector<stk::mesh::GraphEdge>& coincidentEdgesForElem,
                                                  GraphEdgeToParInfoMap& changedEdgesAndParInfos)
{
    GraphEdgeToParInfoMap coincidentParInfos = get_coincident_parallel_infos(graphEdge,
                                                                             coincidentEdgesForElem,
                                                                             parallelInfoForGraphEdges);
    stk::mesh::EntityId chosenSideId = get_min_chosen_side_id(coincidentParInfos);
    update_chosen_side_id_for_coincident_graph_edges(coincidentParInfos, chosenSideId, changedEdgesAndParInfos);
}

void update_chosen_side_ids(const stk::mesh::Graph& graph,
                            stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges,
                            const std::vector<stk::mesh::GraphEdge>& coincidentEdgesForElem,
                            const stk::mesh::impl::LocalId elemId,
                            GraphEdgeToParInfoMap& changedEdgesAndParInfos)
{
    for(const stk::mesh::GraphEdge& graphEdge : graph.get_edges_for_element(elemId))
    {
        if(!stk::mesh::impl::is_local_element(graphEdge.elem2))
        {
            update_chosen_side_ids_for_remote_graph_edge(graphEdge,
                                                         parallelInfoForGraphEdges,
                                                         coincidentEdgesForElem,
                                                         changedEdgesAndParInfos);
        }
    }
}

void update_chosen_side_ids_for_coincident_elems(const stk::mesh::Graph& graph,
                                                 stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges,
                                                 const stk::mesh::impl::SparseGraph& extractedCoincidentElements,
                                                 GraphEdgeToParInfoMap& changedEdgesAndParInfos)
{
    for(const stk::mesh::impl::SparseGraph::value_type& extractedEdgesForElem : extractedCoincidentElements)
    {
        const stk::mesh::impl::LocalId elemId = extractedEdgesForElem.first;
        const std::vector<stk::mesh::GraphEdge>& coincidentEdgesForElem = extractedEdgesForElem.second;
        update_chosen_side_ids(graph, parallelInfoForGraphEdges, coincidentEdgesForElem, elemId, changedEdgesAndParInfos);
    }
}

void choose_face_id_for_coincident_elements(const stk::mesh::Graph &graph,
                                            stk::mesh::ParallelInfoForGraphEdges &parInfosForEdges,
                                            const stk::mesh::impl::SparseGraph &coincidentEdges,
                                            const IdMapper &idMapper,
                                            MPI_Comm comm)
{
    GraphEdgeToParInfoMap changedEdgesAndParInfos;
    update_chosen_side_ids_for_coincident_elems(graph, parInfosForEdges, coincidentEdges, changedEdgesAndParInfos);
    fix_chosen_side_ids_on_other_procs(changedEdgesAndParInfos, parInfosForEdges, idMapper, comm);
}

}}} // end namespaces stk mesh impl

