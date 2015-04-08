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
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Types.hpp>      // for MeshIndex, EntityRank, etc
#include <stk_mesh/baseImpl/BucketRepository.hpp>  // for BucketRepository
#include <stk_mesh/base/EntityLess.hpp>

namespace stk { namespace mesh { namespace unit_test {

class BulkDataTester : public stk::mesh::BulkData
{
public:

    BulkDataTester(stk::mesh::MetaData &mesh_meta_data, MPI_Comm comm) :
            stk::mesh::BulkData(mesh_meta_data, comm)
    {
    }

    BulkDataTester(stk::mesh::MetaData &mesh_meta_data, MPI_Comm comm, stk::mesh::ConnectivityMap const &conn_map) :
            stk::mesh::BulkData(mesh_meta_data, comm, false, &conn_map)
    {
    }

    BulkDataTester(stk::mesh::MetaData &mesh_meta_data,
                   MPI_Comm comm,
                   bool add_fmwk_data,
                   ConnectivityMap const* arg_connectivity_map,
                   FieldDataManager *field_data_manager,
                   unsigned bucket_capacity) :
            stk::mesh::BulkData(mesh_meta_data, comm, add_fmwk_data, arg_connectivity_map, field_data_manager, bucket_capacity)
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


    bool my_internal_modification_end(bool regenerate_aura = true, modification_optimization opt = MOD_END_COMPRESS_AND_SORT)
    {
      return this->internal_modification_end(regenerate_aura, opt);
    }
    void my_internal_change_entity_owner( const std::vector<stk::mesh::EntityProc> & arg_change, bool regenerate_aura = true, modification_optimization mod_optimization = MOD_END_SORT )
    {
        this->internal_change_entity_owner(arg_change,regenerate_aura,mod_optimization);
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

    bool my_internal_modification_end_for_change_entity_owner( bool regenerate_aura, modification_optimization opt )
    {
        return this->internal_modification_end_for_change_entity_owner(regenerate_aura, opt);
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
};

} } } // namespace stk mesh unit_test

#endif
