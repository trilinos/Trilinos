// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

#ifndef stk_mesh_impl_MeshModification_hpp
#define stk_mesh_impl_MeshModification_hpp

#include <stk_mesh/base/Types.hpp>      // for MeshIndex, EntityRank, etc
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityCommListInfo.hpp>
#include "stk_mesh/base/EntityKey.hpp"
#include "stk_mesh/base/EntityParallelState.hpp"
#include "stk_mesh/baseImpl/DeletedEntityCache.hpp"

namespace stk {
class CommSparse;
namespace mesh {

class BulkData;

namespace impl {

class MeshModification
{
public:
    enum BulkDataSyncState { MODIFIABLE = 1 , SYNCHRONIZED = 2 };
    enum modification_optimization {MOD_END_SORT, MOD_END_NO_SORT };

    MeshModification(stk::mesh::BulkData& bulkData) : m_bulkData(bulkData), m_entity_states(),
            m_deleted_entity_cache(bulkData), m_sync_state(MODIFIABLE), m_sync_count(0), m_did_any_shared_entity_change_parts(false)
    {
        m_entity_states.push_back(Deleted);
    }

    ~MeshModification() {}

    BulkDataSyncState synchronized_state() const { return m_sync_state ; }
    bool in_modifiable_state() const { return m_sync_state == MODIFIABLE; }
    bool in_synchronized_state() const { return m_sync_state == SYNCHRONIZED; }
    void set_sync_state_synchronized() { m_sync_state = SYNCHRONIZED; }
    void set_sync_state_modifiable() { m_sync_state = MODIFIABLE; }

    size_t synchronized_count() const { return m_sync_count ; }
    void increment_sync_count() { ++m_sync_count; }
    void set_sync_count(size_t syncCount) { m_sync_count = syncCount; }

    bool modification_begin(const std::string description);

    bool modification_end(modification_optimization opt=MOD_END_SORT);
    bool resolve_node_sharing();
    bool modification_end_after_node_sharing_resolution();

    bool change_entity_owner( const EntityProcVec & arg_change);

    void internal_resolve_shared_modify_delete(
         const std::vector<EntityParallelState>& remotely_modified_shared_entities,
         stk::mesh::EntityProcVec& entitiesToRemoveFromSharing,
         stk::mesh::EntityVector & entitiesNoLongerShared);
    void internal_resolve_ghosted_modify_delete(const std::vector<EntityParallelState >& remotely_modified_ghosted_entities);

    bool did_any_shared_entity_change_parts () const { return m_did_any_shared_entity_change_parts; }
    void set_shared_entity_changed_parts() { m_did_any_shared_entity_change_parts = true; }

    //TODO: these should be Entity::entity_value_type
    bool is_entity_deleted(size_t entity_index) const { return m_entity_states[entity_index] == Deleted; }
    bool is_entity_modified(size_t entity_index) const { return m_entity_states[entity_index] == Modified; }
    bool is_entity_created(size_t entity_index) const { return m_entity_states[entity_index] == Created; }
    bool is_entity_unchanged(size_t entity_index) const { return m_entity_states[entity_index] == Unchanged; }

    stk::mesh::EntityState get_entity_state(size_t entity_index) const { return static_cast<stk::mesh::EntityState>(m_entity_states[entity_index]); }
    void set_entity_state(size_t entity_index, stk::mesh::EntityState state);

    void mark_entity_as_deleted(Entity entity, bool is_ghost)
    {
        set_entity_state(entity.local_offset(), Deleted);
        m_deleted_entity_cache.mark_entity_as_deleted(entity, is_ghost);
    }

    void mark_entity_as_created(size_t entity_index) {  set_entity_state(entity_index, Created); }

    void add_created_entity_state() { m_entity_states.push_back(Created); }

    DeletedEntityCache& get_deleted_entity_cache() { return m_deleted_entity_cache; }

    const DeletedEntityCache& get_deleted_entity_cache() const { return m_deleted_entity_cache; }
    void delete_shared_entities_which_are_no_longer_in_owned_closure(EntityProcVec& entitiesToRemoveFromSharing);
    void internal_change_entity_owner( const std::vector<EntityProc> & local_change,
                                             modification_optimization mod_optimization );
private:
    bool remote_owner_destroyed(EntityKey key, const std::vector<EntityParallelState>& pllStates) const;
    void reset_shared_entity_changed_parts() { m_did_any_shared_entity_change_parts = false; }
    void process_changed_ownership_and_sharing(bool remoteOwnerDestroyed,
                                               Entity entity,
                                               bool shouldRemoveAuraPart = false);
    void internal_establish_new_owner(Entity entity);
    void internal_update_parts_for_shared_entity(Entity entity,
                                                 const bool is_entity_shared,
                                                 const bool did_i_just_become_owner,
                                                 const bool should_remove_aura_part = false);
    void destroy_dependent_ghosts(Entity entity,
                                  EntityProcVec& entitiesToRemoveFromSharing,
                                  EntityVector& auraEntitiesToDestroy);
    void remove_dependent_ghosts(Entity entity,
                                 int remoteProc,
                                 EntityProcVec& entitiesToRemoveFromSharing,
                                 EntityVector& auraEntitiesToDestroy);
    void add_entity_to_same_ghosting(Entity entity, Entity connectedGhost);
    void remove_entities_from_sharing(const EntityProcVec& entitiesToRemoveFromSharing,
                                      EntityVector & entitiesNoLongerShared);
    void internal_resolve_formerly_shared_entities(const EntityVector& entitiesNoLongerShared);
    void reset_undeleted_entity_states_to_unchanged();
    void ensure_meta_data_is_committed();

    stk::mesh::BulkData &m_bulkData;
    std::vector<stk::mesh::EntityState> m_entity_states;
    DeletedEntityCache m_deleted_entity_cache;

    BulkDataSyncState m_sync_state;
    size_t m_sync_count;
    bool m_did_any_shared_entity_change_parts;
};

}}} // end namepsace stk mesh impl


#endif
