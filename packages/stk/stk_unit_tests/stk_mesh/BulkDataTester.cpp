#include "BulkDataTester.hpp"
#include "stk_mesh/base/Entity.hpp"     // for Entity, operator<<, etc
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey, etc
#include <stk_util/parallel/CommSparse.hpp>  // for CommSparse
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>

namespace stk { namespace mesh { namespace unit_test {

bool BulkDataTester::is_entity_in_ghosting_comm_map(stk::mesh::Entity entity)
{
    stk::mesh::EntityKey entityKey = this->entity_key(entity);
    bool is_entity_in_aura_comm_map = !this->internal_entity_comm_map(entityKey, this->aura_ghosting()).empty();
    return is_entity_in_aura_comm_map;
}

void BulkDataTester::check_sharing_comm_maps()
{
    stk::CommSparse comm(parallel());

    for(int phase = 0; phase < 2; ++phase)
    {
        for(stk::mesh::EntityCommListInfoVector::const_iterator i = this->my_internal_comm_list().begin(); i != this->my_internal_comm_list().end(); ++i)
        {
            for(stk::mesh::PairIterEntityComm ec = this->internal_entity_comm_map(i->key); !ec.empty(); ++ec)
            {
                int type = ec->ghost_id;
                if ( type == 0 )
                {
                    std::vector<int> sharingProcs;
                    this->comm_shared_procs(i->key, sharingProcs);
                    // pack shared info
                    int owner = -1;
                    if(bucket_ptr(i->entity) != 0)
                    {
                        owner = parallel_owner_rank(i->entity);
                    }

                    comm.send_buffer(ec->proc).pack<stk::mesh::EntityKey>(i->key).pack<int>(type).pack<int>(owner);
                    comm.send_buffer(ec->proc).pack<size_t>(sharingProcs.size());
                    for (size_t proc=0;proc<sharingProcs.size();++proc)
                    {
                        comm.send_buffer(ec->proc).pack<int>(sharingProcs[proc]);
                    }
                }
            }
        }

        if(phase == 0)
        {
            comm.allocate_buffers();
        }
        else
        {
            comm.communicate();
        }
    }

    // unpack

    std::ostringstream os;
    bool anyErrors = false;

    for(int i = 0; i < parallel_size(); ++i)
    {
        if ( i != parallel_rank() )
        {
            stk::mesh::EntityKey key;
            int type = -1;
            int from = -1;
            int owner = -1;
            while(comm.recv_buffer(i).remaining())
            {
                from = i;
                comm.recv_buffer(from).unpack<stk::mesh::EntityKey>(key).unpack<int>(type).unpack<int>(owner);

                size_t numSharingProcs = 0;
                comm.recv_buffer(from).unpack<size_t>(numSharingProcs);
                std::vector<int> sharingProcs(numSharingProcs);
                for (size_t proc=0;proc<numSharingProcs;++proc)
                {
                    comm.recv_buffer(from).unpack<int>(sharingProcs[proc]);
                }

                std::vector<int> localSharingProcs;
                this->comm_shared_procs(key, localSharingProcs);

                std::sort(localSharingProcs.begin(), localSharingProcs.end());
                std::sort(sharingProcs.begin(), sharingProcs.end());
                size_t maxNum = localSharingProcs.size() + sharingProcs.size();
                std::vector<int> unsharedProcs(maxNum);
                std::vector<int>::iterator iter = std::set_symmetric_difference( localSharingProcs.begin(), localSharingProcs.end(),
                        sharingProcs.begin(), sharingProcs.end(), unsharedProcs.begin());

                size_t numUnshared = iter - unsharedProcs.begin();
                unsharedProcs.resize(numUnshared);

                int counter = 0;
                {
                    for (size_t j=0;j<numUnshared;++j)
                    {
                        std::vector<int>::iterator localIter = std::find(sharingProcs.begin(), sharingProcs.end(), unsharedProcs[j]);
                        if ( localIter != sharingProcs.end() && *localIter != parallel_rank() )
                        {
                            if ( counter == 0 )
                            {
                                os << "Error in sharing between procs for entity " << key.id() << " with rank " << key.rank()  << "  between procs: " << this->parallel_rank() << " and " << from << std::endl;
                                counter++;
                            }
                            os << "\tEntity " << key << " is shared with proc " << unsharedProcs[j] << " from other proc: "
                                    << from << " but not from this proc: " << parallel_rank() << std::endl;
                            anyErrors = true;
                        }

                        localIter = std::find(localSharingProcs.begin(), localSharingProcs.end(), unsharedProcs[j]);
                        if ( localIter != localSharingProcs.end() && *localIter != from )
                        {
                            if ( counter == 0 )
                            {
                                os << "Error in sharing between procs for entity " << key.id() << " with rank " << key.rank()  << "  between procs: " << this->parallel_rank() << " and " << from << std::endl;
                                counter++;
                            }
                            os << "\tEntity " << key << " is shared with proc " << unsharedProcs[j] << " from this proc: "
                                    << parallel_rank() << " but not from other proc: " << from << std::endl;
                            anyErrors = true;
                        }
                    }
                }
            }
        }
    }

    ThrowRequireMsg(!anyErrors, os.str());
}

///////////////////////////

unsigned BulkDataFaceSharingTester::get_index_of_side_in_element_bucket(stk::mesh::Entity element, stk::mesh::Entity side)
{
    unsigned elements_edge_offset = stk::mesh::INVALID_CONNECTIVITY_ORDINAL;
    unsigned num_edges_or_faces = this->num_connectivity(element, this->entity_rank(side));
    const stk::mesh::Entity* entities = this->begin(element, this->entity_rank(side));
    for(unsigned j=0;j<num_edges_or_faces;++j)
    {
        if (entities[j]==side)
        {
            elements_edge_offset = static_cast<stk::mesh::ConnectivityOrdinal>(j);
            break;
        }
    }
    return elements_edge_offset;
}

stk::mesh::Permutation BulkDataFaceSharingTester::get_permutation(stk::mesh::Entity element, const stk::mesh::EntityVector& nodes)
{
    std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> ordinalAndPermutation =
                  stk::mesh::get_ordinal_and_permutation(*this, element, mesh_meta_data().side_rank(), nodes);
    return ordinalAndPermutation.second;
}

stk::mesh::EntityVector BulkDataFaceSharingTester::convert_keys_to_entities(const std::vector<stk::mesh::EntityKey>& node_keys)
{
    stk::mesh::EntityVector nodes(node_keys.size());
    for (size_t i=0;i<nodes.size();++i)
        nodes[i] = this->get_entity(node_keys[i]);
    return nodes;
}

void BulkDataFaceSharingTester::change_connectivity_for_edge_or_face(stk::mesh::Entity side, const std::vector<stk::mesh::EntityKey>& node_keys)
{
    stk::mesh::EntityVector nodes = convert_keys_to_entities(node_keys);

    unsigned num_elems = this->num_elements(side);
    const stk::mesh::Entity *elements = this->begin_elements(side);
    for (unsigned i=0;i<num_elems;++i)
    {
        if(bucket(elements[i]).owned())
        {
            stk::mesh::unit_test::BucketTester& bucket_edge = static_cast<stk::mesh::unit_test::BucketTester&>(this->bucket(side));
            bucket_edge.my_change_exisiting_connectivity(this->bucket_ordinal(side), &nodes[0]);

            stk::mesh::Permutation new_permutation = get_permutation(elements[i], nodes);
            ThrowRequireMsg(new_permutation!=stk::mesh::INVALID_PERMUTATION, "Program error. Please contact sierra-help@sandia.gov for support.");


            unsigned edges_element_offset = static_cast<stk::mesh::ConnectivityOrdinal>(i);
            bucket_edge.my_change_exisiting_permutation_for_connected_element(this->bucket_ordinal(side), edges_element_offset, new_permutation);

            unsigned elements_edge_offset = get_index_of_side_in_element_bucket(elements[i], side);
            stk::mesh::unit_test::BucketTester& bucket_elem = static_cast<stk::mesh::unit_test::BucketTester&>(this->bucket(elements[i]));
            bucket_elem.my_change_exisiting_permutation_for_connected_edge(this->bucket_ordinal(elements[i]), elements_edge_offset, new_permutation);
        }
    }
}

void BulkDataFaceSharingTester::change_entity_key_and_nodes(const std::vector<shared_entity_type> & potentially_shared_sides)
{
    for(size_t i = 0, e = potentially_shared_sides.size(); i < e; ++i)
    {
        Entity entity = potentially_shared_sides[i].entity;
        if(potentially_shared_sides[i].need_update_nodes)
        {
            if(potentially_shared_sides[i].global_key != potentially_shared_sides[i].local_key)
            {
                my_internal_change_entity_key(potentially_shared_sides[i].local_key, potentially_shared_sides[i].global_key, entity);
            }
            change_connectivity_for_edge_or_face(entity, potentially_shared_sides[i].nodes);
        }
    }
}

void BulkDataFaceSharingTester::change_entity_key_and_update_sharing_info(std::vector<shared_entity_type> & potentially_shared_sides)
{
   change_entity_key_and_nodes(potentially_shared_sides);
   insert_sharing_info_into_comm_map(potentially_shared_sides);
}

void BulkDataFaceSharingTester::update_shared_entity_this_proc2(EntityKey global_key_other_proc, shared_entity_type& shared_entity, int proc_id, const std::vector<EntityKey>& nodes)
{
    Entity entity = shared_entity.entity;
    shared_entity.sharing_procs.push_back(proc_id);
    if(proc_id < this->parallel_rank())
    {
        shared_entity.global_key = global_key_other_proc;
        shared_entity.nodes = nodes;
        shared_entity.need_update_nodes = true;
    }
    this->internal_mark_entity(entity, BulkData::IS_SHARED);
}

stk::mesh::EntityVector get_elements_connected_to_all_nodes(const BulkDataFaceSharingTester& bulkData, const stk::mesh::EntityVector& nodes)
{
    stk::mesh::EntityVector elements;
    stk::mesh::impl::find_locally_owned_elements_these_nodes_have_in_common(bulkData, nodes.size(), nodes.data(), elements);
    return elements;
}

void add_side_to_shared_entities(stk::mesh::Entity side, std::vector<stk::mesh::shared_entity_type>& shared_entities_this_proc, const stk::mesh::shared_entity_type &shared_entity_other_proc, int other_proc_id)
{
    stk::mesh::shared_entity_type sentity(shared_entity_other_proc.global_key, side, shared_entity_other_proc.topology);
    sentity.nodes = shared_entity_other_proc.nodes;
    sentity.sharing_procs.push_back(other_proc_id);
    shared_entities_this_proc.push_back(sentity);
}

stk::mesh::Entity create_side_and_add_to_shared_entity_list(stk::mesh::Entity element, stk::mesh::EntityRank side_rank, const stk::mesh::EntityVector& nodes, const stk::mesh::shared_entity_type &shared_entity_other_proc,
        stk::mesh::BulkData& bulkData, std::vector<stk::mesh::shared_entity_type>& shared_entities_this_proc, int other_proc_id, stk::mesh::Part& root_topo_part)
{
    stk::mesh::Entity side = stk::mesh::declare_element_to_sub_topology_with_nodes(bulkData, element, nodes, shared_entity_other_proc.global_key.id(), side_rank, root_topo_part);
    ThrowRequireMsg(bulkData.is_valid(side), "Program error. Contact sierra-help@sandia.gov for support.");
    add_side_to_shared_entities(side, shared_entities_this_proc, shared_entity_other_proc, other_proc_id);
    return side;
}

void BulkDataFaceSharingTester::connect_side_from_other_proc_to_local_elements(const stk::mesh::EntityVector& elements, const stk::mesh::EntityVector& nodes, const stk::mesh::shared_entity_type &shared_entity_other_proc,
        stk::mesh::BulkData& bulkData, std::vector<stk::mesh::shared_entity_type>& shared_entities_this_proc, int other_proc_id)
{
    for(stk::mesh::Entity element : elements)
    {
        stk::mesh::Entity side = create_side_and_add_to_shared_entity_list(element, side_rank(), nodes, shared_entity_other_proc, *this,
                shared_entities_this_proc, other_proc_id, get_topology_root_part(shared_entity_other_proc.topology));
        this->internal_mark_entity(side, BulkData::IS_SHARED);
    }
}

void BulkDataFaceSharingTester::create_and_connect_shared_face_on_this_proc(const stk::mesh::shared_entity_type &shared_entity_other_proc, std::vector<stk::mesh::shared_entity_type>& shared_entities_this_proc, int other_proc_id)
{
    stk::mesh::EntityVector nodes = convert_keys_to_entities(shared_entity_other_proc.nodes);
    stk::mesh::EntityVector elements = get_elements_connected_to_all_nodes(*this, nodes);
    connect_side_from_other_proc_to_local_elements(elements, nodes, shared_entity_other_proc, *this, shared_entities_this_proc, other_proc_id);
}

void BulkDataFaceSharingTester::check_if_entity_from_other_proc_exists_on_this_proc_and_update_info_if_shared(std::vector<stk::mesh::shared_entity_type>& shared_entities_this_proc,
        int proc_id, const stk::mesh::shared_entity_type &shared_entity_other_proc)
{
    int matching_index = does_entity_exist_in_list(shared_entities_this_proc, shared_entity_other_proc);
    bool entitiesAreTheSame = matching_index >= 0;
    if( entitiesAreTheSame )
        update_shared_entity_this_proc2(shared_entity_other_proc.global_key, shared_entities_this_proc[matching_index], proc_id, shared_entity_other_proc.nodes);
    else
        create_and_connect_shared_face_on_this_proc(shared_entity_other_proc, shared_entities_this_proc, proc_id);
}

///////////////////////////

stk::mesh::EntityVector BulkDataElemGraphFaceSharingTester::get_local_sides() const
{
    stk::mesh::EntityVector sides;
    stk::mesh::get_selected_entities(locally_owned_part(), buckets(side_rank()), sides);
    return sides;
}

stk::mesh::EntityVector get_nodes_of_entity(const BulkDataElemGraphFaceSharingTester& bulkData, stk::mesh::Entity entity)
{
    ThrowRequireMsg(bulkData.num_nodes(entity)>0, "Program error. Contact sierra-help@sandia.gov for support.");
    stk::mesh::EntityVector nodes(bulkData.begin_nodes(entity), bulkData.end_nodes(entity));
    return nodes;
}

std::vector<stk::mesh::EntityKey> convert_entities_to_keys(const BulkDataElemGraphFaceSharingTester& bulkData, const stk::mesh::EntityVector& nodes)
{
    std::vector<stk::mesh::EntityKey> node_keys(nodes.size());
    for(size_t i=0;i<node_keys.size();++i)
        node_keys[i] = bulkData.entity_key(nodes[i]);
    return node_keys;
}

std::vector<stk::mesh::EntityKey> get_node_keys_of_entity(const BulkDataElemGraphFaceSharingTester& bulkData, stk::mesh::Entity entity)
{
    stk::mesh::EntityVector nodes = get_nodes_of_entity(bulkData, entity);
    return convert_entities_to_keys(bulkData, nodes);
}

void add_shared_entity(const BulkDataElemGraphFaceSharingTester &bulkData, stk::mesh::Entity side, std::vector<shared_entity_type>& shared_entities, int proc_id)
{
    shared_entity_type sentity(bulkData.entity_key(side), side, bulkData.get_entity_topology(side));
    sentity.nodes = get_node_keys_of_entity(bulkData, side);
    sentity.sharing_procs.push_back(proc_id);
    shared_entities.push_back(sentity);
}

void BulkDataElemGraphFaceSharingTester::add_side_if_remote_element_connection_exists(int num_connected_elems, stk::mesh::Entity element, stk::mesh::Entity side, const stk::mesh::ElemElemGraph& egraph,
        std::vector<shared_entity_type>& shared_entities)
{
    for(int k=0;k<num_connected_elems;++k)
    {
        if(!egraph.is_connected_elem_locally_owned(element, k))
        {
            add_shared_entity(*this, side, shared_entities, egraph.get_owning_proc_id_of_remote_element(element, k));
            break;
        }
    }
}

void BulkDataElemGraphFaceSharingTester::determine_if_side_is_shared_across_proc_boundary_and_add(const stk::mesh::EntityVector& elements, stk::mesh::Entity side,
        const stk::mesh::ElemElemGraph& egraph, std::vector<shared_entity_type>& shared_entities)
{
    for(stk::mesh::Entity element : elements)
    {
        if(is_entity_owned(element))
        {
            int num_connected_elems = egraph.get_num_connected_elems(element);
            add_side_if_remote_element_connection_exists(num_connected_elems, element, side, egraph, shared_entities);
        }
    }
}

void BulkDataElemGraphFaceSharingTester::fill_shared_entities_that_need_fixing(const stk::mesh::EntityVector& sides, const stk::mesh::ElemElemGraph& egraph, std::vector<shared_entity_type>& shared_entities)
{
    for(stk::mesh::Entity side : sides)
    {
        if(entityWasCreatedThisModCycle(side))
        {
            std::vector<stk::mesh::Entity> elements(begin_elements(side),end_elements(side));
            determine_if_side_is_shared_across_proc_boundary_and_add(elements, side, egraph, shared_entities);
        }
    }
}

void BulkDataElemGraphFaceSharingTester::mark_entities_as_shared(const std::vector<shared_entity_type>& sentities)
{
    for(shared_entity_type sentity : sentities)
        internal_mark_entity(sentity.entity, BulkData::IS_SHARED);
}

void BulkDataElemGraphFaceSharingTester::use_elem_elem_graph_to_determine_shared_entities(std::vector<shared_entity_type>& shared_entities)
{
    stk::mesh::ElemElemGraph egraph(*this, mesh_meta_data().universal_part());
    stk::mesh::EntityVector sides = get_local_sides();
    fill_shared_entities_that_need_fixing(sides, egraph, shared_entities);
    mark_entities_as_shared(shared_entities);
}


void BulkDataElemGraphFaceSharingTester::markEntitiesForResolvingSharingInfoUsingNodes(stk::mesh::EntityRank entityRank, std::vector<shared_entity_type>& shared_entities)
{
    if(entityRank == mesh_meta_data().side_rank())
        use_elem_elem_graph_to_determine_shared_entities(shared_entities);
    else
        BulkData::markEntitiesForResolvingSharingInfoUsingNodes(entityRank, shared_entities);
}

} } } // namespace stk mesh unit_test
