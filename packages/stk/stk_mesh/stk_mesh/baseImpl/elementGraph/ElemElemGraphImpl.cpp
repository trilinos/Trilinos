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
#include <stk_mesh/baseImpl/EquivalentEntityBlocks.hpp>
#include <stk_mesh/baseImpl/elementGraph/BulkDataIdMapper.hpp>

#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/util/SortAndUnique.hpp>

namespace stk { namespace mesh {

namespace impl
{

unsigned get_num_local_elems(const stk::mesh::BulkData& bulkData)
{
    if(bulkData.mesh_meta_data().entity_rank_count() < stk::topology::ELEM_RANK)
        return 0;
    else
        return count_selected_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData.buckets(stk::topology::ELEM_RANK));
}

void fill_topologies(const stk::mesh::BulkData& bulkData,
                     const stk::mesh::impl::ElementLocalIdMapper & localMapper,
                     std::vector<stk::topology>& element_topologies)
{
    const stk::mesh::BucketVector & elemBuckets = bulkData.get_buckets(stk::topology::ELEM_RANK, bulkData.mesh_meta_data().locally_owned_part());
    for(const stk::mesh::Bucket* bucket : elemBuckets) {
        for(stk::mesh::Entity element : *bucket) {
            element_topologies[localMapper.entity_to_local(element)] = bucket->topology();
        }
    }
}

ElemSideToProcAndFaceId build_element_side_ids_to_proc_map(const stk::mesh::BulkData& bulkData,
                                                           const stk::mesh::EntityVector &elements_to_communicate)
{
    ElemSideToProcAndFaceId elem_side_comm;
    stk::mesh::EntityVector side_nodes;
    std::vector<stk::mesh::EntityKey> keys;
    std::vector<int> sharing_procs;
    for(size_t i=0;i<elements_to_communicate.size();++i)
    {
        stk::mesh::Entity elem = elements_to_communicate[i];
        stk::topology elem_top = bulkData.bucket(elem).topology();
        unsigned num_sides = elem_top.num_sides();
        for(unsigned side=0;side<num_sides;++side)
        {
            fill_element_side_nodes_from_topology(bulkData, elem, side, side_nodes);
            keys.resize(side_nodes.size());
            for(size_t j=0;j<keys.size();++j)
            {
                keys[j] = bulkData.entity_key(side_nodes[j]);
            }
            bulkData.shared_procs_intersection(keys, sharing_procs);
            for (int proc: sharing_procs) {
                elem_side_comm.insert(std::pair<EntitySidePair, ProcFaceIdPair>(EntitySidePair(elem, side), ProcFaceIdPair(proc,0)));
            }
        }
    }
    return elem_side_comm;
}

void fill_element_side_nodes_from_topology(const stk::mesh::BulkData& bulkData, stk::mesh::Entity element, unsigned side_index, stk::mesh::EntityVector& localElemSideNodes)
{
    stk::topology localElemTopology = bulkData.bucket(element).topology();
    const stk::mesh::Entity* localElemNodes = bulkData.begin_nodes(element);

    unsigned num_nodes_this_side = localElemTopology.side_topology(side_index).num_nodes();
    if(num_nodes_this_side == 0) return;

    localElemSideNodes.resize(num_nodes_this_side);

    localElemTopology.side_nodes(localElemNodes, side_index, localElemSideNodes.begin());
}

bool does_element_have_side(const stk::mesh::BulkData& bulkData, stk::mesh::Entity element)
{
    unsigned dimension_of_element = bulkData.bucket(element).topology().dimension();
    unsigned dimension_of_mesh = bulkData.mesh_meta_data().spatial_dimension();
    return dimension_of_element == dimension_of_mesh;
}

void pack_elements_to_comm(stk::CommSparse &comm, const std::vector<GraphEdgeProc>& elements_to_comm)
{
    for(size_t i=0;i<elements_to_comm.size();++i)
    {
        int remote_proc = elements_to_comm[i].get_remote_processor_rank();
        comm.send_buffer(remote_proc).pack<stk::mesh::EntityId>(elements_to_comm[i].get_local_element_global_id());
        comm.send_buffer(remote_proc).pack<int>(elements_to_comm[i].get_local_element_side_index());
        comm.send_buffer(remote_proc).pack<stk::mesh::EntityId>(elements_to_comm[i].get_remote_element_global_id());
        comm.send_buffer(remote_proc).pack<int>(elements_to_comm[i].get_remote_element_side_index());
    }
}

std::vector<GraphEdgeProc> communicate_killed_entities(stk::ParallelMachine communicator, const std::vector<GraphEdgeProc>& elements_to_comm)
{
    std::vector<GraphEdgeProc> remote_edges;
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
            remote_edges.push_back(GraphEdgeProc(localId, localSide, remoteId, remoteSide, i));
        }
    }
    return remote_edges;
}

int get_element_side_multiplier()
{
    return 1000;
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

stk::mesh::ConstPartVector get_stk_parts_for_moving_parts_into_death_boundary(const stk::mesh::PartVector *bc_mesh_parts)
{
    stk::mesh::ConstPartVector sideParts;
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

void add_side_into_exposed_boundary(stk::mesh::BulkData& bulkData, const ParallelInfo& parallel_edge_info,
        stk::mesh::Entity local_element, int side_id, stk::mesh::EntityId remote_id, const stk::mesh::PartVector& parts_for_creating_side,
        std::vector<stk::mesh::sharing_info> &shared_modified, stk::mesh::impl::ParallelSelectedInfo &remoteActiveSelector, const stk::mesh::PartVector *boundary_mesh_parts)
{
//    stk::mesh::EntityId side_global_id = parallel_edge_info.m_chosen_side_id;
    stk::mesh::ConnectivityOrdinal side_ord = static_cast<stk::mesh::ConnectivityOrdinal>(side_id);

    // determine which element is active
//    stk::mesh::Permutation perm = stk::mesh::DEFAULT_PERMUTATION;
    int other_proc = parallel_edge_info.get_proc_rank_of_neighbor();
    int owning_proc = std::min(other_proc, bulkData.parallel_rank());

//    if(remoteActiveSelector[-remote_id])
//    {
//        perm = static_cast<stk::mesh::Permutation>(parallel_edge_info.m_permutation);
//    }

    stk::mesh::Entity side = stk::mesh::get_side_entity_for_elem_side_pair(bulkData, local_element, side_id);

    if(!bulkData.is_valid(side))
    {
        stk::mesh::PartVector side_parts = get_parts_for_creating_side(bulkData, parts_for_creating_side, local_element, side_id);
//        ThrowRequireWithSierraHelpMsg(!impl::is_id_already_in_use_locally(bulkData, bulkData.mesh_meta_data().side_rank(), side_global_id));
//        side = connect_side_to_element(bulkData, local_element, side_global_id, side_ord, perm, side_parts);
        side = bulkData.declare_element_side(local_element, side_ord, side_parts);
        shared_modified.push_back(stk::mesh::sharing_info(side, other_proc, owning_proc));
    }
    else
    {
        if(bulkData.bucket(side).owned())
        {
            stk::mesh::ConstPartVector parts = get_stk_parts_for_moving_parts_into_death_boundary(boundary_mesh_parts);
            bulkData.change_entity_parts(side, parts);
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
    stk::mesh::Entity side = stk::mesh::get_side_entity_for_elem_side_pair(bulkData, local_element, side_id);
    if(side_created_during_death(bulkData, side))
    {
        deletedEntities.push_back(side);
    }
    else if(bulkData.is_valid(side) && bulkData.bucket(side).owned())
    {
        bulkData.change_entity_parts(side, stk::mesh::ConstPartVector{}, stk::mesh::ConstPartVector{&activePart});
    }
}


stk::mesh::Entity connect_side_to_element(stk::mesh::BulkData& bulkData, stk::mesh::Entity element,
        stk::mesh::EntityId side_global_id, stk::mesh::ConnectivityOrdinal side_ordinal,
        stk::mesh::Permutation side_permutation, const stk::mesh::PartVector& parts)
{
    stk::mesh::EntityRank side_rank = bulkData.mesh_meta_data().side_rank();

    stk::mesh::Entity side = bulkData.internal_declare_entity(side_rank, side_global_id, parts);

    // connect element to side
    bulkData.declare_relation(element, side, side_ordinal, side_permutation);

    // connect side to nodes
    stk::topology side_top = bulkData.bucket(element).topology().side_topology(side_ordinal);
    stk::mesh::EntityVector side_nodes;
    fill_element_side_nodes_from_topology(bulkData, element, side_ordinal, side_nodes);
    stk::mesh::EntityVector permuted_side_nodes(side_top.num_nodes());
    side_top.permutation_nodes(side_nodes, side_permutation, permuted_side_nodes.begin());
    for(size_t i=0;i<permuted_side_nodes.size();++i)
    {
        bulkData.declare_relation(side, permuted_side_nodes[i], i);
    }

    return side;
}

void pack_newly_shared_remote_edges(stk::CommSparse &comm, const stk::mesh::BulkData &bulkData, const std::vector<SharedEdgeInfo> &newlySharedEdges)
{
    std::vector<SharedEdgeInfo>::const_iterator iter = newlySharedEdges.begin();
    std::vector<SharedEdgeInfo>::const_iterator endIter = newlySharedEdges.end();

    for(; iter!= endIter; ++iter)
    {
        stk::mesh::EntityId localId = iter->get_local_element_global_id();
        stk::mesh::Entity localEntity = bulkData.get_entity(stk::topology::ELEM_RANK, localId);
        stk::mesh::EntityId remoteId = iter->get_remote_element_global_id();
        int local_side_index    = iter->get_local_element_side_index();
        int remote_side_index    = iter->get_remote_element_side_index();
        int sharing_proc       = iter->get_remote_processor_rank();

        size_t numNodes= iter->m_sharedNodes.size();
        std::vector<stk::mesh::EntityKey> side_node_entity_keys(numNodes);
        for(size_t i=0; i<numNodes; ++i)
        {
            side_node_entity_keys[i] = bulkData.entity_key(iter->m_sharedNodes[i]);
        }

        comm.send_buffer(sharing_proc).pack<stk::mesh::EntityId>(localId);
        comm.send_buffer(sharing_proc).pack<stk::mesh::EntityId>(remoteId);
        comm.send_buffer(sharing_proc).pack<int>(local_side_index);
        comm.send_buffer(sharing_proc).pack<int>(remote_side_index);
        comm.send_buffer(sharing_proc).pack<stk::topology>(bulkData.bucket(localEntity).topology());
        comm.send_buffer(sharing_proc).pack<unsigned>(numNodes);
        for(size_t i=0; i<numNodes; ++i)
        {
            comm.send_buffer(sharing_proc).pack<stk::mesh::EntityKey>(side_node_entity_keys[i]);
        }
    }
}

bool is_local_element(stk::mesh::impl::LocalId elemId)
{
    return (elemId >= 0);
}

void add_exposed_sides(LocalId elementId, size_t maxSidesThisElement,
                      const stk::mesh::Graph &graph, std::vector<int> &element_side_pairs)
{
    std::vector<int> elemSides(maxSidesThisElement, -1);
    for(size_t j = 0; j < graph.get_num_edges_for_element(elementId); ++j)
    {
        const stk::mesh::GraphEdge & graphEdge = graph.get_edge_for_element(elementId, j);
        int sideId = graphEdge.side1();
        elemSides[sideId] = sideId;
    }

    for(size_t sideId = 0; sideId < maxSidesThisElement; ++sideId)
        if (elemSides[sideId] == -1)
            element_side_pairs.push_back(sideId);
}


}}} // end namespaces stk mesh impl

