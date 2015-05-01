// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
#ifndef _BulkDataTester_hpp_
#define _BulkDataTester_hpp_

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/FEMHelpers.hpp>   // for BulkData, etc
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Types.hpp>      // for MeshIndex, EntityRank, etc
#include <stk_mesh/baseImpl/BucketRepository.hpp>  // for BucketRepository
#include <stk_mesh/base/EntityLess.hpp>
#include "BucketTester.hpp"

namespace stk { namespace mesh { namespace unit_test {

inline int does_entity_exist_in_list(std::vector<stk::mesh::shared_entity_type>& shared_entity_map, stk::mesh::shared_entity_type &sentity)
{
    int matching_index = -1;
    for (size_t i=0;i<shared_entity_map.size();++i)
    {
        stk::topology topo1 = shared_entity_map[i].topology;
        stk::topology topo2 = sentity.topology;
        size_t num_nodes1 = shared_entity_map[i].nodes.size();
        size_t num_nodes2 = sentity.nodes.size();
        if (topo1 == topo2 && num_nodes1 == num_nodes2)
        {
            bool sameType = topo1.equivalent(shared_entity_map[i].nodes, sentity.nodes).first;
            if (sameType)
            {
                matching_index = i;
                break;
            }
        }
    }
    return matching_index;
}

class BulkDataTester : public stk::mesh::BulkData
{
public:

    BulkDataTester(stk::mesh::MetaData &mesh_meta_data, MPI_Comm comm) :
            stk::mesh::BulkData(mesh_meta_data, comm)
    {
    }

    BulkDataTester(stk::mesh::MetaData &mesh_meta_data, MPI_Comm comm, enum stk::mesh::BulkData::AutomaticAuraOption auto_aura_option) :
            stk::mesh::BulkData(mesh_meta_data, comm, auto_aura_option)
    {
    }

    BulkDataTester(stk::mesh::MetaData &mesh_meta_data, MPI_Comm comm, stk::mesh::ConnectivityMap const &conn_map) :
            stk::mesh::BulkData(mesh_meta_data, comm, stk::mesh::BulkData::AUTO_AURA, false, &conn_map)
    {
    }

    BulkDataTester(stk::mesh::MetaData &mesh_meta_data,
                   MPI_Comm comm,
                   enum stk::mesh::BulkData::AutomaticAuraOption auto_aura_option,
                   bool add_fmwk_data,
                   ConnectivityMap const* arg_connectivity_map,
                   FieldDataManager *field_data_manager,
                   unsigned bucket_capacity) :
            stk::mesh::BulkData(mesh_meta_data, comm, auto_aura_option, add_fmwk_data, arg_connectivity_map, field_data_manager, bucket_capacity)
    {
    }

    virtual ~BulkDataTester()
    {
    }

    bool is_entity_in_ghosting_comm_map(stk::mesh::Entity entity);

    void check_sharing_comm_maps();

    uint16_t closure_count(stk::mesh::Entity entity)
    {
        return m_closure_count[entity.local_offset()];
    }

    void reset_closure_count(stk::mesh::Entity entity)
    {
        m_closure_count[entity.local_offset()] = 0;
    }

    uint16_t my_orphaned_node_marking()
    {
        return orphaned_node_marking;
    }

    const stk::mesh::EntityCommDatabase my_entity_comm_map() const
    {
        return m_entity_comm_map;
    }

    bool is_ghosted_somewhere(stk::mesh::EntityKey key) const
    {
        return !internal_entity_comm_map(key, aura_ghosting()).empty();
    }

    void my_internal_change_entity_owner( const std::vector<stk::mesh::EntityProc> & arg_change, bool regenerate_aura = true, modification_optimization mod_optimization = MOD_END_SORT )
    {
        this->internal_change_entity_owner(arg_change,mod_optimization);
    }

    void my_resolve_ownership_of_modified_entities(const std::vector<stk::mesh::Entity> &shared_new)
    {
        this->resolve_ownership_of_modified_entities(shared_new);
    }

    bool my_entity_comm_map_insert(stk::mesh::Entity entity, const stk::mesh::EntityCommInfo & val)
    {
        return BulkData::entity_comm_map_insert(entity, val);
    }

    bool my_entity_comm_map_erase(const stk::mesh::EntityKey& key, const stk::mesh::EntityCommInfo& commInfo)
    {
        return BulkData::entity_comm_map_erase(key, commInfo);
    }

    bool my_entity_comm_map_erase(const stk::mesh::EntityKey& key, const stk::mesh::Ghosting& ghost)
    {
        return BulkData::entity_comm_map_erase(key, ghost);
    }

    void my_entity_comm_map_clear(const stk::mesh::EntityKey& key)
    {
        BulkData::entity_comm_map_clear(key);
    }

    void my_entity_comm_map_clear_ghosting(const stk::mesh::EntityKey& key)
    {
        BulkData::entity_comm_map_clear_ghosting(key);
    }

    bool my_internal_modification_end_for_change_entity_owner(modification_optimization opt )
    {
        return this->internal_modification_end_for_change_entity_owner(opt);
    }

    bool my_is_entity_in_sharing_comm_map(stk::mesh::Entity entity)
    {
        return this->is_entity_in_sharing_comm_map(entity);
    }

    void my_update_sharing_after_change_entity_owner()
    {
        this->update_sharing_after_change_entity_owner();
    }

    stk::mesh::impl::EntityRepository &my_get_entity_repository()
    {
        return get_entity_repository();
    }

    inline bool my_set_parallel_owner_rank_but_not_comm_lists(stk::mesh::Entity entity, int in_owner_rank)
    {
        return this->internal_set_parallel_owner_rank_but_not_comm_lists(entity, in_owner_rank);
    }

    bool my_internal_set_parallel_owner_rank_but_not_comm_lists(stk::mesh::Entity entity, int in_owner_rank)
    {
        return this->internal_set_parallel_owner_rank_but_not_comm_lists(entity, in_owner_rank);
    }

    void my_fix_up_ownership(stk::mesh::Entity entity, int new_owner)
    {
        this->fix_up_ownership(entity, new_owner);
    }

    stk::mesh::PairIterEntityComm my_internal_entity_comm_map_shared(const stk::mesh::EntityKey & key) const
    {
        return internal_entity_comm_map_shared(key);
    }

    int my_internal_entity_comm_map_owner(const EntityKey & key) const
    {
        return internal_entity_comm_map_owner(key);
    }

    const EntityCommListInfoVector & my_internal_comm_list() const
    {
        return internal_comm_list();
    }

    PairIterEntityComm my_internal_entity_comm_map(const EntityKey & key) const
    {
        return internal_entity_comm_map(key);
    }

    PairIterEntityComm my_internal_entity_comm_map(const EntityKey & key, const Ghosting & sub) const
    {
        return internal_entity_comm_map(key, sub);
    }

    bool my_comm_mesh_verify_parallel_consistency(std::ostream & error_log)
    {
        return comm_mesh_verify_parallel_consistency(error_log);
    }

    void my_internal_resolve_shared_modify_delete()
    {
        this->internal_resolve_shared_modify_delete();
    }

    void my_internal_resolve_ghosted_modify_delete()
    {
        this->internal_resolve_ghosted_modify_delete();
    }

    void my_internal_resolve_parallel_create()
    {
        this->internal_resolve_parallel_create();
    }

    void my_update_comm_list_based_on_changes_in_comm_map()
    {
        this->update_comm_list_based_on_changes_in_comm_map();
    }

    void my_internal_update_distributed_index(std::vector<stk::mesh::Entity> & shared_new )
    {
        this->internal_update_sharing_comm_map_and_fill_list_modified_shared_entities( shared_new );
    }

    void my_internal_update_distributed_index(stk::mesh::EntityRank entityRank, std::vector<stk::mesh::Entity> & shared_new )
    {
        this->internal_update_sharing_comm_map_and_fill_list_modified_shared_entities_of_rank(entityRank, shared_new);
    }

    void my_move_entities_to_proper_part_ownership( std::vector<stk::mesh::Entity> &shared_modified )
    {
        this->move_entities_to_proper_part_ownership( shared_modified );
    }

    void my_add_comm_list_entries_for_entities( std::vector<stk::mesh::Entity> &shared_modified )
    {
        this->add_comm_list_entries_for_entities( shared_modified );
    }

    void my_internal_resolve_shared_membership()
    {
        this->internal_resolve_shared_membership();
    }

    void my_internal_regenerate_aura()
    {
        this->internal_regenerate_aura();
    }

    void my_set_state(stk::mesh::Entity entity, stk::mesh::EntityState entity_state)
    {
        set_state(entity,entity_state);
    }

    void my_delete_shared_entities_which_are_no_longer_in_owned_closure()
    {
        delete_shared_entities_which_are_no_longer_in_owned_closure();
    }

    void my_ghost_entities_and_fields(Ghosting & ghosting, const std::set<EntityProc , EntityLess>& new_send)
    {
        ghost_entities_and_fields(ghosting, new_send);
    }

    void my_add_closure_entities(const stk::mesh::Ghosting& ghosting, const stk::mesh::EntityProcVec& entities, std::set <stk::mesh::EntityProc , stk::mesh::EntityLess > &entitiesWithClosure)
    {
        add_closure_entities(ghosting, entities, entitiesWithClosure);
    }

    void my_internal_modification_end_for_change_ghosting()
    {
        internal_modification_end_for_change_ghosting();
    }

    bool my_in_send_ghost(const stk::mesh::Ghosting& ghosting, stk::mesh::EntityKey key, int proc)
    {
        return in_send_ghost(ghosting, key, proc);
    }

    void my_markEntitiesForResolvingSharingInfoUsingNodes(stk::mesh::EntityRank entityRank, std::vector<stk::mesh::shared_entity_type>& shared_entities)
    {
        markEntitiesForResolvingSharingInfoUsingNodes(entityRank, shared_entities);
    }

    void my_fillSharedEntities(stk::mesh::Ghosting& ghost_id,
                            stk::mesh::BulkData &mesh,
                            std::vector<shared_entity_type> & shared_entity_map,
                            std::vector<std::vector<shared_entity_type> > &shared_entities)
    {
        fillSharedEntities(ghost_id, mesh, shared_entity_map, shared_entities);
    }

    void my_unpackEntityInfromFromOtherProcsAndMarkEntitiesAsSharedAndTrackProcessorsThatNeedAlsoHaveEntity(stk::CommSparse &comm,
            std::vector<stk::mesh::shared_entity_type> & shared_entity_map)
    {
        unpackEntityInfromFromOtherProcsAndMarkEntitiesAsSharedAndTrackProcessorsThatNeedAlsoHaveEntity(comm, shared_entity_map);
    }

    void my_internal_change_entity_key(EntityKey old_key, EntityKey new_key, Entity entity)
    {
        internal_change_entity_key(old_key, new_key, entity);
    }

    stk::mesh::BulkData::entitySharing my_internal_is_entity_marked(stk::mesh::Entity entity) const
    {
        return internal_is_entity_marked(entity);
    }
};

class BulkDataFaceSharingTester : public BulkDataTester
{
public:
    BulkDataFaceSharingTester(stk::mesh::MetaData &mesh_meta_data, MPI_Comm comm) :
            BulkDataTester(mesh_meta_data, comm)
    {
    }

    ~BulkDataFaceSharingTester(){}

    void change_connectivity_for_edge_or_face(stk::mesh::Entity edgeOrFace, std::vector<stk::mesh::EntityKey>& node_keys)
    {
        stk::mesh::EntityVector nodes(node_keys.size());
        for (size_t i=0;i<nodes.size();++i)
        {
            nodes[i] = this->get_entity(node_keys[i]);
        }

        unsigned edges_element_offset = stk::mesh::INVALID_CONNECTIVITY_ORDINAL;
        unsigned elements_edge_offset = stk::mesh::INVALID_CONNECTIVITY_ORDINAL;
        unsigned num_elems = this->num_elements(edgeOrFace);
        const stk::mesh::Entity *elements = this->begin_elements(edgeOrFace);
        for (unsigned i=0;i<num_elems;++i)
        {
            edges_element_offset = static_cast<stk::mesh::ConnectivityOrdinal>(i);
            std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> ordinalAndPermutation =
                          stk::mesh::get_ordinal_and_permutation(*this, elements[i], stk::topology::EDGE_RANK, nodes);
            stk::mesh::Permutation new_permutation = ordinalAndPermutation.second;

            stk::mesh::unit_test::BucketTester& bucket_edge = static_cast<stk::mesh::unit_test::BucketTester&>(this->bucket(edgeOrFace));

            bucket_edge.my_change_exisiting_connectivity(this->bucket_ordinal(edgeOrFace), &nodes[0]);
            bucket_edge.my_change_exisiting_permutation_for_connected_element(this->bucket_ordinal(edgeOrFace), edges_element_offset, new_permutation);

            unsigned num_edges_or_faces = this->num_connectivity(elements[i], this->entity_rank(edgeOrFace));
            const stk::mesh::Entity* entities = this->begin(elements[i], this->entity_rank(edgeOrFace));
            for(unsigned j=0;j<num_edges_or_faces;++j)
            {
                if (entities[j]==edgeOrFace)
                {
                    elements_edge_offset = static_cast<stk::mesh::ConnectivityOrdinal>(j);
                    break;
                }
            }

            stk::mesh::unit_test::BucketTester& bucket_elem = static_cast<stk::mesh::unit_test::BucketTester&>(this->bucket(elements[i]));
            bucket_elem.my_change_exisiting_permutation_for_connected_edge(this->bucket_ordinal(elements[i]), elements_edge_offset, new_permutation);
        }
    }

    virtual void resolveUniqueIdForSharedEntityAndCreateCommMapInfoForSharingProcs(std::vector<shared_entity_type> & shared_entity_map)
    {
       for(size_t i = 0, e = shared_entity_map.size(); i < e; ++i)
       {
           Entity entity = get_entity(shared_entity_map[i].local_key);
           if(shared_entity_map[i].need_update_nodes)
           {
               if(shared_entity_map[i].global_key != shared_entity_map[i].local_key)
               {
                   my_internal_change_entity_key(shared_entity_map[i].local_key, shared_entity_map[i].global_key, entity);
               }
               change_connectivity_for_edge_or_face(entity, shared_entity_map[i].nodes);

           }
           for(size_t j = 0; j < shared_entity_map[i].sharing_procs.size(); j++)
           {
               entity_comm_map_insert(entity, EntityCommInfo(stk::mesh::BulkData::SHARED, shared_entity_map[i].sharing_procs[j]));
           }
       }
    }

    virtual void is_entity_shared(std::vector<stk::mesh::shared_entity_type>& shared_entity_map, int proc_id, stk::mesh::shared_entity_type &sentity)
    {
        int matching_index = does_entity_exist_in_list(shared_entity_map, sentity);
        bool entitiesAreTheSame = matching_index >= 0;

        if( entitiesAreTheSame )
        {
            Entity entity = this->get_entity(shared_entity_map[matching_index].local_key);
            shared_entity_map[matching_index].sharing_procs.push_back(proc_id);
            if(proc_id < this->parallel_rank())
            {
                shared_entity_map[matching_index].global_key = sentity.global_key;
                shared_entity_map[matching_index].nodes = sentity.nodes;
                shared_entity_map[matching_index].need_update_nodes = true;
            }
            this->internal_mark_entity(entity, BulkData::IS_SHARED);
        }
    }

    virtual void sortNodesIfNeeded(std::vector<stk::mesh::EntityKey>& nodes)
    {
        // do not do any sorting
    }
};

} } } // namespace stk mesh unit_test

#endif
