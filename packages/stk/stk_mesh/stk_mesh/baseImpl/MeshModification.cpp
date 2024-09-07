#include "MeshModification.hpp"
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include <stk_mesh/baseImpl/EntityKeyMapping.hpp>
#include <stk_mesh/base/EntityLess.hpp>
#include <stk_mesh/base/SideSetHelper.hpp>
#include <stk_mesh/base/Relation.hpp>
#include <stk_mesh/baseImpl/CommEntityMods.hpp>
#include <stk_mesh/baseImpl/Visitors.hpp>

namespace stk {
namespace mesh {
namespace impl {

bool MeshModification::modification_begin(const std::string description)
{
    if (m_bulkData.m_runConsistencyCheck) {
      parallel_machine_barrier( m_bulkData.parallel() );
    }

    if (this->synchronized_count() == 0)
    {
        m_bulkData.mesh_meta_data().set_mesh_on_fields(&m_bulkData);
        m_bulkData.m_entityKeyMapping->update_num_ranks(m_bulkData.mesh_meta_data().entity_rank_count());
        const unsigned numRanks = m_bulkData.mesh_meta_data().entity_rank_count();
        if (numRanks > m_bulkData.m_selector_to_buckets_maps.size()) {
          m_bulkData.m_selector_to_buckets_maps.resize(numRanks);
        }
    }

    if ( this->in_modifiable_state() ) return false ;

    if (this->synchronized_count() == 0)
    {
        this->ensure_meta_data_is_committed();

        if (m_bulkData.parallel_size() > 1) {
            verify_parallel_consistency( m_bulkData.mesh_meta_data() , m_bulkData.parallel() );
        }
    }
    else
    {
        this->reset_undeleted_entity_states_to_unchanged();
         m_bulkData.m_removedGhosts.clear();
    }

    this->set_sync_state_modifiable();
    this->reset_shared_entity_changed_parts();

    const stk::mesh::FieldVector allFields = m_bulkData.mesh_meta_data().get_fields();
    for (FieldBase * stkField : allFields) {
      stkField->sync_to_host();
      if (stkField->has_ngp_field()) {
        impl::get_ngp_field(*stkField)->debug_modification_begin();
      }
    }

    this->increment_sync_count();
    return true;
}

void MeshModification::set_entity_state(size_t entity_index, stk::mesh::EntityState state)
{
  m_entity_states[entity_index] = state;
}

bool MeshModification::modification_end(modification_optimization opt)
{
  if(this->in_synchronized_state())
  {
      return false;
  }

  STK_ThrowAssertMsg(impl::check_for_connected_nodes(m_bulkData)==0, "BulkData::modification_end ERROR, all entities with rank higher than node are required to have connected nodes.");

  STK_ThrowAssertMsg(m_bulkData.add_fmwk_data() || impl::check_no_shared_elements_or_higher(m_bulkData)==0, "BulkData::modification_end ERROR, Sharing of entities with rank ELEMENT_RANK or higher is not allowed.");

  m_bulkData.m_entityKeyMapping->clear_all_cache();

  if(m_bulkData.parallel_size() > 1)
  {
      // Resolve modification or deletion of shared entities
      // which can cause deletion of ghost entities.
      stk::mesh::EntityProcVec entitiesToRemoveFromSharing;
      delete_shared_entities_which_are_no_longer_in_owned_closure(entitiesToRemoveFromSharing);

      CommEntityMods commEntityMods(m_bulkData, m_bulkData.internal_comm_db(), m_bulkData.internal_comm_list());
      commEntityMods.communicate(CommEntityMods::PACK_ALL);
      stk::mesh::EntityVector entitiesNoLongerShared;
      internal_resolve_shared_modify_delete(commEntityMods.get_shared_mods(), entitiesToRemoveFromSharing, entitiesNoLongerShared);

      // Resolve modification or deletion of ghost entities
      // by destroying ghost entities that have been touched.
      internal_resolve_ghosted_modify_delete(commEntityMods.get_ghosted_mods());

      m_bulkData.update_comm_list_based_on_changes_in_comm_map();

      // Resolve creation of entities: discover sharing and set unique ownership.
      m_bulkData.internal_resolve_parallel_create();

      // Resolve part membership for shared entities.
      // This occurs after resolving creation so created and shared
      // entities are resolved along with previously existing shared entities.
      m_bulkData.internal_resolve_shared_membership(entitiesNoLongerShared);

      // Regenerate the ghosting aura around all shared mesh entities.
      if(m_bulkData.is_automatic_aura_on())
      {
          m_bulkData.internal_regenerate_aura();
      }
      else if (m_bulkData.m_turningOffAutoAura) {
          m_bulkData.internal_remove_aura();
      }

      m_bulkData.internal_resolve_send_ghost_membership();

      m_bulkData.m_modSummary.write_summary(synchronized_count());
  }

  m_bulkData.internal_finish_modification_end(opt);

  if(m_bulkData.parallel_size() > 1)
  {
      m_bulkData.check_mesh_consistency();
  }

  return true;
}

bool MeshModification::resolve_node_sharing()
{
    if(this->in_synchronized_state())
    {
        return false;
    }

    STK_ThrowAssertMsg(impl::check_for_connected_nodes(m_bulkData)==0, "BulkData::modification_end ERROR, all entities with rank higher than node are required to have connected nodes.");

    STK_ThrowAssertMsg(m_bulkData.add_fmwk_data() || impl::check_no_shared_elements_or_higher(m_bulkData)==0, "BulkData::modification_end ERROR, Sharing of entities with rank ELEMENT_RANK or higher is not allowed.");

    if(m_bulkData.parallel_size() > 1)
    {
        m_bulkData.internal_resolve_parallel_create_nodes();
    }
    else
    {
        m_bulkData.m_modSummary.write_summary(synchronized_count());
        if(!m_bulkData.add_fmwk_data())
        {
            std::vector<Entity> shared_modified;
            m_bulkData.internal_update_sharing_comm_map_and_fill_list_modified_shared_entities(shared_modified);
        }
    }

    return true;
}

bool MeshModification::modification_end_after_node_sharing_resolution()
{
    if(this->in_synchronized_state())
    {
        return false;
    }

    STK_ThrowAssertMsg(impl::check_for_connected_nodes(m_bulkData)==0, "BulkData::modification_end ERROR, all entities with rank higher than node are required to have connected nodes.");

    STK_ThrowAssertMsg(m_bulkData.add_fmwk_data() || impl::check_no_shared_elements_or_higher(m_bulkData)==0, "BulkData::modification_end ERROR, Sharing of entities with rank ELEMENT_RANK or higher is not allowed.");

    if(m_bulkData.parallel_size() > 1)
    {
        m_bulkData.internal_resolve_parallel_create_edges_and_faces();
        stk::mesh::EntityVector entitiesNoLongerShared;
        stk::mesh::EntityProcVec entitiesToRemoveFromSharing;
        delete_shared_entities_which_are_no_longer_in_owned_closure(entitiesToRemoveFromSharing);

        CommEntityMods commEntityMods(m_bulkData, m_bulkData.internal_comm_db(), m_bulkData.internal_comm_list());
        commEntityMods.communicate(CommEntityMods::PACK_SHARED);
        internal_resolve_shared_modify_delete(commEntityMods.get_shared_mods(), entitiesToRemoveFromSharing, entitiesNoLongerShared);
        m_bulkData.internal_resolve_shared_membership(entitiesNoLongerShared);

        if(m_bulkData.is_automatic_aura_on())
        {
            m_bulkData.internal_regenerate_aura();
        }
        else if (m_bulkData.m_turningOffAutoAura) {
            m_bulkData.internal_remove_aura();
        }

        m_bulkData.internal_resolve_send_ghost_membership();

        m_bulkData.m_modSummary.write_summary(synchronized_count());
        m_bulkData.check_mesh_consistency();
    }
    else
    {
        m_bulkData.m_modSummary.write_summary(synchronized_count());
        if(!m_bulkData.add_fmwk_data())
        {
            std::vector<Entity> shared_modified;
            m_bulkData.internal_update_sharing_comm_map_and_fill_list_modified_shared_entities(shared_modified);
        }
    }

    m_bulkData.internal_finish_modification_end(MOD_END_SORT);

    return true;
}

bool MeshModification::change_entity_owner( const EntityProcVec & arg_change)
{
    STK_ThrowRequireMsg(in_synchronized_state(), "BulkData::change_entity_owner() must not be called from within a modification cycle.");
    std::vector<EntityProc> local_change( arg_change );
    const bool validChangesOnAnyProc = impl::internal_clean_and_verify_parallel_change(m_bulkData, local_change);
    if (!validChangesOnAnyProc) {
      return false;
    }

    m_bulkData.notifier.notify_elements_about_to_move_procs(local_change);

    modification_optimization mod_optimization = MOD_END_SORT;
    modification_begin("change_entity_owner");
    internal_change_entity_owner(local_change, mod_optimization);
    m_bulkData.update_sharing_after_change_entity_owner();
    m_bulkData.internal_modification_end_for_change_entity_owner(mod_optimization);

    m_bulkData.notifier.notify_elements_moved_procs(local_change);

    return true;
}

void MeshModification::internal_change_entity_owner( const std::vector<EntityProc> & local_change,
                                             modification_optimization mod_optimization )
{
  m_bulkData.require_ok_to_modify();
  m_bulkData.m_modSummary.track_change_entity_owner(local_change);

  const MetaData  & meta = m_bulkData.mesh_meta_data() ;
  const int       p_rank = m_bulkData.parallel_rank() ;
  const int       p_size = m_bulkData.parallel_size() ;
  ParallelMachine p_comm = m_bulkData.parallel() ;

  //------------------------------
  // internal_change_entity_owner can assume a clean local change list, it was
  // checked in MeshModification::change_entity_owner, which called this method.

  //----------------------------------------
  // Parallel synchronous determination of changing shared and ghosted.

  // The two vectors below will contain changes to ghosted and shared
  // entities on this process coming from change-entity-owner requests
  // on other processes.
  std::vector<EntityProc> ghosted_change ;
  std::vector<EntityProc> shared_change ;

  impl::internal_generate_parallel_change_lists( m_bulkData , local_change ,
                            shared_change , ghosted_change ); 

  //------------------------------
  // Have enough information to delete all effected ghosts.
  // If the closure of a ghost contains a changing entity
  // then that ghost must be deleted.
  // Request that all ghost entities in the closure of the ghost be deleted.

  std::set<EntityProc,EntityLess> send_closure(m_bulkData);
  impl::StoreInEntityProcSet store_entity_proc_in_set(m_bulkData, send_closure);

  // Compute the closure of all the locally changing entities
  for (const EntityProc& entityProc : local_change) {
      store_entity_proc_in_set.proc = entityProc.second; 
      impl::VisitClosureGeneral(m_bulkData, entityProc.first, m_bulkData.entity_rank(entityProc.first), store_entity_proc_in_set, store_entity_proc_in_set);
  }

  // Calculate all the ghosts that are impacted by the set of ownership
  // changes. We look at ghosted, shared, and local changes looking for ghosts
  // that are either in the closure of the changing entity, or have the
  // changing entity in their closure. All modified ghosts will be removed.
  {
    impl::OnlyVisitGhostsOnce only_visit_ghosts_once(m_bulkData);
    impl::StoreEntity store_entity(m_bulkData);

    std::vector<EntityProc>& allChanges = ghosted_change;
    allChanges.reserve(allChanges.size()+shared_change.size()+send_closure.size());
    allChanges.insert(allChanges.end(), shared_change.begin(), shared_change.end());
    allChanges.insert(allChanges.end(), local_change.begin(), local_change.end());
    impl::VisitAuraClosureGeneral(m_bulkData,allChanges.begin(),allChanges.end(),store_entity,only_visit_ghosts_once);

    std::vector<Entity> remove_modified_ghosts;
    store_entity.store_visited_entities_in_vec(remove_modified_ghosts);

    std::vector<EntityProc> empty_add ;
    std::vector<Entity> removesForThisGhosting;
    removesForThisGhosting.reserve(remove_modified_ghosts.size());
    const bool notAddingSendGhosts = true;

    // Skip 'm_ghosting[0]' which is the shared subset.
    for (unsigned i=1; i<m_bulkData.m_ghosting.size(); ++i) {
      removesForThisGhosting.clear();
      for(Entity entity : remove_modified_ghosts) {
        if (m_bulkData.in_receive_ghost(*m_bulkData.m_ghosting[i], entity)) {
          removesForThisGhosting.push_back(entity);
        }
      }

      m_bulkData.internal_change_ghosting(*m_bulkData.m_ghosting[i], empty_add, removesForThisGhosting, notAddingSendGhosts);
    }
  }

  //------------------------------
  // Consistently change the owner on all processes.
  // 1) The local_change list is giving away ownership.
  // 2) The shared_change may or may not be receiving ownership
  {
    ConstPartVector owned;
    owned.push_back(& meta.locally_owned_part());
    OrdinalVector scratchOrdinalVec, scratchSpace;

    for (const EntityProc& entityProc : local_change) {
      // Giving ownership, change the parts first and then
      // the owner rank to pass the ownership test.
      Entity entity = entityProc.first;

      m_bulkData.internal_verify_and_change_entity_parts( entity , ConstPartVector() , owned,
                        scratchOrdinalVec, scratchSpace );

      m_bulkData.internal_set_owner(entity, entityProc.second);
    }

    for (const EntityProc& entityProc : shared_change) {
      Entity entity = entityProc.first;
      m_bulkData.internal_set_owner(entity, entityProc.second);
      if ( p_rank == entityProc.second ) { // I received ownership
          m_bulkData.internal_verify_and_change_entity_parts( entity , owned , ConstPartVector(),
                        scratchOrdinalVec, scratchSpace );
      }
    }
  }

  //------------------------------
  // Send entities, along with their closure, to the new owner processes
  {
    std::ostringstream error_msg ;
    int error_count = 0 ;

    stk::CommSparse comm( p_comm );

    EntityVector unique_list_of_send_closure;
    unique_list_of_send_closure.reserve(send_closure.size());

    const bool onlyPackDownwardRelations = true;
    for ( std::set<EntityProc,EntityLess>::iterator
          i = send_closure.begin() ; i != send_closure.end() ; ++i ) {
      CommBuffer & buffer = comm.send_buffer( i->second );
      Entity entity = i->first;
      pack_entity_info(m_bulkData, buffer, entity, onlyPackDownwardRelations);
      if (!m_bulkData.is_communicated_with_proc(entity, i->second) ||
          std::binary_search(local_change.begin(), local_change.end(), *i, EntityLess(m_bulkData))) {
        buffer.pack<int>(1);
        pack_field_values(m_bulkData, buffer , entity );
      }
      else {
        buffer.pack<int>(0);
      }
      pack_sideset_info(m_bulkData, buffer , entity );

      if (unique_list_of_send_closure.empty() || m_bulkData.entity_key(unique_list_of_send_closure.back()) != m_bulkData.entity_key(entity)) {
        unique_list_of_send_closure.push_back(entity);
      }
    }

    comm.allocate_buffers();

    for ( std::set<EntityProc,EntityLess>::iterator
          i = send_closure.begin() ; i != send_closure.end() ; ++i ) {
      CommBuffer & buffer = comm.send_buffer( i->second );
      Entity entity = i->first;
      pack_entity_info(m_bulkData, buffer, entity, onlyPackDownwardRelations);
      if (!m_bulkData.is_communicated_with_proc(entity, i->second) ||
          std::binary_search(local_change.begin(), local_change.end(), *i, EntityLess(m_bulkData))) {
        buffer.pack<int>(1);
        pack_field_values(m_bulkData, buffer , entity );
      }
      else {
        buffer.pack<int>(0);
      }
      pack_sideset_info(m_bulkData, buffer , entity );
    }

    const bool deallocateSendBuffers = true;
    comm.communicate(deallocateSendBuffers);

    SideSetHelper helper(m_bulkData, m_bulkData.mesh_meta_data().universal_part());
    for ( std::set<EntityProc,EntityLess>::iterator
          i = send_closure.begin() ; i != send_closure.end() ; ++i ) {
      Entity entity = i->first;
      helper.remove_element_entries_from_sidesets(entity);
    }

    OrdinalVector partOrdinals;
    OrdinalVector scratchOrdinalVec, scratchSpace;
    PartVector parts ;
    std::vector<Relation> relations ;

    OrdinalVector removeCustomGhostParts;
    const std::vector<Ghosting*>& ghostingObjs = m_bulkData.ghostings();
    const unsigned firstCustomGhosting = 2;
    for(unsigned i=firstCustomGhosting; i<ghostingObjs.size(); ++i) {
      removeCustomGhostParts.push_back(m_bulkData.ghosting_part(*ghostingObjs[i]).mesh_meta_data_ordinal());
    }

    for ( int p = 0 ; p < p_size ; ++p ) {
      CommBuffer & buf = comm.recv_buffer(p);
      while ( buf.remaining() ) {
        EntityKey key ;
        int owner = ~0u ;

        parts.clear();
        relations.clear();
        unpack_entity_info( buf, m_bulkData, key, owner, parts, relations );

        // Received entity information will be correct,
        // modulo the owned and shared parts

        remove( parts , meta.globally_shared_part() );

        if ( owner == p_rank ) {
          // Must have the locally_owned_part
          insert( parts , meta.locally_owned_part() );
        }
        else {
          // Must not have the locally_owned_part
          remove( parts , meta.locally_owned_part() );
        }

        std::pair<Entity ,bool> result = m_bulkData.internal_create_entity( key );

        Entity entity = result.first;

        // The entity was copied and not created.
        partOrdinals.clear();
        for(const stk::mesh::Part* part : parts) {
            partOrdinals.push_back(part->mesh_meta_data_ordinal());
        }

        m_bulkData.internal_change_entity_parts( entity , partOrdinals , removeCustomGhostParts, scratchOrdinalVec, scratchSpace );
        for(unsigned i=firstCustomGhosting; i<ghostingObjs.size(); ++i) {
          m_bulkData.entity_comm_map_erase(key, EntityCommInfo(ghostingObjs[i]->ordinal(), p));
        }

        if (m_bulkData.state(entity) == Created) {
          set_entity_state(entity.local_offset(), Modified);
        }

        m_bulkData.internal_set_owner(entity, owner);

        m_bulkData.internal_declare_relation( entity , relations, scratchOrdinalVec );

        int shouldUnpackFieldValues = 0;
        buf.unpack<int>(shouldUnpackFieldValues);
        if ( shouldUnpackFieldValues==1 ) {
          if ( ! unpack_field_values(m_bulkData, buf , entity , error_msg ) ) {
            ++error_count ;
          }
        }

        unpack_sideset_info( buf, m_bulkData, entity);
      }
    }

#ifndef NDEBUG
    all_reduce( p_comm , ReduceSum<1>( & error_count ) );
#endif
    STK_ThrowRequireMsg(error_count==0, error_msg.str() );

    // Any entity that I sent and is not in an owned closure is deleted.
    // The owned closure will be effected by received entities, so can
    // only clean up after the newly owned entities have been received.
    // Destroy backwards so as not to invalidate closures in the process.
    {
        for ( EntityVector::reverse_iterator i = unique_list_of_send_closure.rbegin() ; i != unique_list_of_send_closure.rend() ; ++i) {
            stk::mesh::Entity entity = *i;
            if ( ! m_bulkData.owned_closure(entity) ) {
                for(unsigned ig=firstCustomGhosting; ig<ghostingObjs.size(); ++ig) {
                  m_bulkData.entity_comm_map_erase(m_bulkData.entity_key(entity), *ghostingObjs[ig]);
                }
                STK_ThrowRequireMsg( m_bulkData.internal_destroy_entity( entity ), "Failed to destroy entity " << m_bulkData.identifier(entity) );
            }
        }
    }
    send_closure.clear(); // Has been invalidated
  }

  m_bulkData.update_comm_list_based_on_changes_in_comm_map();
}

// Resolve modifications for shared entities:
// If not locally destroyed and remotely modified
// then set to locally modified.
// If remotely destroyed then determine the new owner.
//
// Post condition:
//  Shared entities are in-sync with respect to modification state.
//  Shared communication lists are updated to reflect all deletions.
//  Ownership has been re-assigned as necessary for deletion
//  of shared entities.

void MeshModification::internal_resolve_shared_modify_delete(
         const std::vector<EntityParallelState>& remotely_modified_shared_entities,
         stk::mesh::EntityProcVec& entitiesToRemoveFromSharing,
         stk::mesh::EntityVector & entitiesNoLongerShared)
{
    STK_ThrowRequireMsg(m_bulkData.parallel_size() > 1, "Do not call this in serial");
    EntityVector auraEntitiesToDestroy;

    // We iterate backwards to ensure that we hit the higher-ranking entities first.
    for(std::vector<EntityParallelState>::const_reverse_iterator
    i = remotely_modified_shared_entities.rbegin(); i != remotely_modified_shared_entities.rend();)
    {
        Entity entity = i->comm_info.entity;
        EntityKey key = i->comm_info.key;
        const int owner = m_bulkData.parallel_owner_rank(entity);
        const bool locally_destroyed = !m_bulkData.is_valid(entity);
        bool remote_owner_destroyed = false;

        // Iterate over all of this entity's remote changes
        for(; i != remotely_modified_shared_entities.rend() && i->comm_info.entity == entity; ++i)
        {

            const int remote_proc = i->from_proc;
            const bool remotely_destroyed = Deleted == i->state;

            // When a shared entity is remotely modified or destroyed
            // then the local copy is also modified.  This modification
            // status is applied to all related higher ranking entities.
            if(!locally_destroyed)
            {
                m_bulkData.mark_entity_and_upward_related_entities_as_modified(entity);
            }

            // A shared entity is being deleted on the remote process.
            // Remove it from the sharing communication list.
            // Ownership changes are processed later, but we'll need
            // to know if the remote owner destroyed the entity in order
            // to correctly resolve ownership (it is not sufficient to just
            // look at the comm list of the entity since there is no
            // guarantee that the comm list is correct or up-to-date).

            const bool not_shared_remotely = !(i->remote_owned_closure);
            if(remotely_destroyed || not_shared_remotely)
            {
                if (remotely_destroyed && remote_proc == owner) {
                  m_bulkData.entity_comm_map_clear_ghosting(key);
                }

                m_bulkData.entity_comm_map_erase(key, EntityCommInfo(stk::mesh::BulkData::SHARED, remote_proc));
                if (!locally_destroyed && not_shared_remotely) {
                    entitiesToRemoveFromSharing.push_back(EntityProc(entity, remote_proc));
                    remove_dependent_ghosts(entity, remote_proc, entitiesToRemoveFromSharing, auraEntitiesToDestroy);
                }
            }
            // check if owner is destroying
            if(remotely_destroyed && (owner == remote_proc))
            {
                remote_owner_destroyed = true;
            }
        }

        // Have now processed all remote changes knowledge for this entity.

        if(!locally_destroyed)
        {
            process_changed_ownership_and_sharing(remote_owner_destroyed, entity);
        }
    } // remote mod loop

    // Erase all sharing communication lists for Destroyed entities:
    for(EntityCommListInfoVector::const_reverse_iterator
    i = m_bulkData.internal_comm_list().rbegin(); i != m_bulkData.internal_comm_list().rend(); ++i)
    {
        if(!m_bulkData.is_valid(i->entity))
        {
            m_bulkData.entity_comm_map_erase(i->key, m_bulkData.shared_ghosting());
        }
    }

    stk::util::sort_and_unique(entitiesToRemoveFromSharing);
    m_bulkData.remove_entities_from_sharing(entitiesToRemoveFromSharing, entitiesNoLongerShared);
    internal_resolve_formerly_shared_entities(entitiesNoLongerShared);

    stk::util::sort_and_unique(auraEntitiesToDestroy, EntityLess(m_bulkData));

    for(EntityVector::const_reverse_iterator iter = auraEntitiesToDestroy.rbegin();
        iter != auraEntitiesToDestroy.rend(); ++iter) {
      m_bulkData.destroy_entity(*iter);
    }
}

void MeshModification::process_changed_ownership_and_sharing(bool remoteOwnerDestroyed,
                                                             Entity entity,
                                                             bool shouldRemoveAuraPart)
{
  const bool am_i_old_local_owner = m_bulkData.parallel_rank() == m_bulkData.parallel_owner_rank(entity);

  if(remoteOwnerDestroyed) {
      internal_establish_new_owner(entity);
  }

  const bool am_i_new_local_owner = m_bulkData.parallel_rank() == m_bulkData.parallel_owner_rank(entity);
  const bool did_i_just_become_owner = (!am_i_old_local_owner && am_i_new_local_owner );
  const bool is_entity_shared = !m_bulkData.internal_entity_comm_map_shared(m_bulkData.entity_key(entity)).empty();
  internal_update_parts_for_shared_entity(entity, is_entity_shared, did_i_just_become_owner, shouldRemoveAuraPart);
}

void MeshModification::internal_establish_new_owner(stk::mesh::Entity entity)
{
  const int new_owner = m_bulkData.determine_new_owner(entity);
  m_bulkData.internal_set_owner(entity, new_owner);
}

void MeshModification::internal_update_parts_for_shared_entity(Entity entity,
                                                               const bool is_entity_shared,
                                                               const bool did_i_just_become_owner,
                                                               const bool should_remove_aura_part)
{
  OrdinalVector parts_to_add_entity_to , parts_to_remove_entity_from, scratchOrdinalVec, scratchSpace;

  if (is_entity_shared) {
    const bool notPreviouslyShared = !m_bulkData.bucket(entity).shared();
    if (notPreviouslyShared) {
      parts_to_add_entity_to.push_back(m_bulkData.mesh_meta_data().globally_shared_part().mesh_meta_data_ordinal());
    }
  }
  else {
    parts_to_remove_entity_from.push_back(m_bulkData.mesh_meta_data().globally_shared_part().mesh_meta_data_ordinal());
  }

  if (should_remove_aura_part) {
    parts_to_remove_entity_from.push_back(m_bulkData.mesh_meta_data().aura_part().mesh_meta_data_ordinal());
  }

  if ( did_i_just_become_owner ) {
    parts_to_add_entity_to.push_back(m_bulkData.mesh_meta_data().locally_owned_part().mesh_meta_data_ordinal());
  }

  if ( ! parts_to_add_entity_to.empty() || ! parts_to_remove_entity_from.empty() ) {
    m_bulkData.internal_change_entity_parts( entity , parts_to_add_entity_to , parts_to_remove_entity_from, scratchOrdinalVec, scratchSpace );
  }
}

void MeshModification::destroy_dependent_ghosts(Entity entity,
                                                EntityProcVec& entitiesToRemoveFromSharing,
                                                EntityVector& auraEntitiesToDestroy)
{
  EntityRank entity_rank = m_bulkData.entity_rank(entity);
  const EntityRank end_rank = static_cast<EntityRank>(m_bulkData.mesh_meta_data().entity_rank_count());
  for (EntityRank irank = static_cast<EntityRank>(end_rank - 1); irank > entity_rank; --irank)
  {
    int num_rels = m_bulkData.num_connectivity(entity, irank);
    const Entity* rels     = m_bulkData.begin(entity, irank);

    for (int r = num_rels - 1; r >= 0; --r)
    {
      Entity e = rels[r];

      const bool upwardRelationOfEntityIsInClosure = m_bulkData.owned_closure(e);
      STK_ThrowRequireMsg( !upwardRelationOfEntityIsInClosure, m_bulkData.entity_rank(e) << " with id " << m_bulkData.identifier(e) << " should not be in closure." );

      // Recursion
      if (m_bulkData.is_valid(e) && m_bulkData.bucket(e).in_aura()) {
          destroy_dependent_ghosts(e, entitiesToRemoveFromSharing, auraEntitiesToDestroy);
      }
    }
  }

  for(EntityRank downwardRank=stk::topology::NODE_RANK; downwardRank < entity_rank; ++downwardRank) {
    const unsigned numConnected = m_bulkData.num_connectivity(entity, downwardRank);
    const Entity* connected = m_bulkData.begin(entity, downwardRank);
    for(unsigned ic=0; ic<numConnected; ++ic) {
      if (m_bulkData.is_valid(connected[ic]) && m_bulkData.bucket(connected[ic]).in_aura()) {
        auraEntitiesToDestroy.push_back(connected[ic]);
      }
    }
  }

  const bool successfully_destroyed_entity = m_bulkData.destroy_entity(entity);
  if (!successfully_destroyed_entity) {
    std::vector<int> sharing_procs;
    m_bulkData.comm_shared_procs(m_bulkData.entity_key(entity), sharing_procs);
    for(int p : sharing_procs) {
      entitiesToRemoveFromSharing.emplace_back(entity, p);
    }
  }
}

void MeshModification::remove_dependent_ghosts(Entity entity,
                                               int remoteProc,
                                               EntityProcVec& entitiesToRemoveFromSharing,
                                               EntityVector& auraEntitiesToDestroy)
{
  EntityRank entity_rank = m_bulkData.entity_rank(entity);
  const EntityRank end_rank = static_cast<EntityRank>(m_bulkData.mesh_meta_data().entity_rank_count());
  for (EntityRank irank = static_cast<EntityRank>(end_rank - 1); irank > entity_rank; --irank)
  {
    int num_rels = m_bulkData.num_connectivity(entity, irank);
    const Entity* rels     = m_bulkData.begin(entity, irank);

    for (int r = num_rels - 1; r >= 0; --r)
    {
      Entity e = rels[r];

      // Recursion
      if (m_bulkData.is_valid(e)) {
        remove_dependent_ghosts(e, remoteProc, entitiesToRemoveFromSharing, auraEntitiesToDestroy);
      }
    }
  }

  for(EntityRank downwardRank=stk::topology::NODE_RANK; downwardRank < entity_rank; ++downwardRank) {
    const unsigned numConnected = m_bulkData.num_connectivity(entity, downwardRank);
    const Entity* connected = m_bulkData.begin(entity, downwardRank);
    for(unsigned ic=0; ic<numConnected; ++ic) {
      if (m_bulkData.is_valid(connected[ic]) && m_bulkData.in_ghost(m_bulkData.aura_ghosting(), connected[ic], remoteProc)) {
        if (m_bulkData.bucket(connected[ic]).owned()) {
          m_bulkData.entity_comm_map_erase(m_bulkData.entity_key(connected[ic]), m_bulkData.aura_ghosting());
        }
        else {
          auraEntitiesToDestroy.push_back(connected[ic]);
        }
      }
    }
  }
}

// Entities with sharing information that are not in the owned closure
// have been modified such that they are no longer shared.
// These may no longer be needed or may become ghost entities.
// There is not enough information so assume they are to be deleted
// and let these entities be re-ghosted if they are needed.

// Open question: Should an owned and shared entity that does not
// have an upward relation to an owned entity be destroyed so that
// ownership transfers to another process?

void MeshModification::delete_shared_entities_which_are_no_longer_in_owned_closure(EntityProcVec& entitiesToRemoveFromSharing)
{
  EntityVector auraEntitiesToDestroy;
  auraEntitiesToDestroy.reserve(m_bulkData.internal_comm_list().size()/2);

  for ( EntityCommListInfoVector::const_reverse_iterator
        i =  m_bulkData.internal_comm_list().rbegin() ;
        i != m_bulkData.internal_comm_list().rend() ; ++i)
  {
    Entity entity = i->entity;
    if (m_bulkData.is_valid(entity) && !m_bulkData.owned_closure(entity)) {
      if (i->entity_comm != -1 && m_bulkData.in_shared(entity)) {
        destroy_dependent_ghosts(entity, entitiesToRemoveFromSharing, auraEntitiesToDestroy);
      }
    }
  }

  stk::util::sort_and_unique(auraEntitiesToDestroy, EntityLess(m_bulkData));

  for(EntityVector::const_reverse_iterator iter = auraEntitiesToDestroy.rbegin();
      iter != auraEntitiesToDestroy.rend(); ++iter) {
    m_bulkData.destroy_entity(*iter);
  }
}

bool MeshModification::remote_owner_destroyed(EntityKey key,
                                              const std::vector<EntityParallelState>& pllStates) const
{
  std::vector<EntityParallelState>::const_iterator iter = std::lower_bound(pllStates.begin(), pllStates.end(), key);
  if (iter != pllStates.end() && iter->comm_info.key == key) {
    return iter->from_proc == m_bulkData.parallel_owner_rank(iter->comm_info.entity)
         && iter->state == Deleted;
  }
  return false;
}

//----------------------------------------------------------------------
// Resolve modifications for ghosted entities:
// If a ghosted entity is modified or destroyed on the owning
// process then the ghosted entity must be destroyed.
//
// Post condition:
//  Ghosted entities of modified or deleted entities are destroyed.
//  Ghosted communication lists are cleared to reflect all deletions.

void MeshModification::internal_resolve_ghosted_modify_delete(const std::vector<EntityParallelState >& remotely_modified_ghosted_entities)
{
  STK_ThrowRequireMsg(m_bulkData.parallel_size() > 1, "Do not call this in serial");
  // Resolve modifications for ghosted entities:

  const size_t ghosting_count = m_bulkData.m_ghosting.size();
  const size_t ghosting_count_minus_shared = ghosting_count - 1;

  std::vector<Entity> promotingToShared;

  // We iterate backwards over remote_mod to ensure that we hit the
  // higher-ranking entities first. This is important because higher-ranking
  // entities like element must be deleted before the nodes they have are
  // deleted.
  for ( std::vector<EntityParallelState>::const_reverse_iterator
        i = remotely_modified_ghosted_entities.rbegin(); i != remotely_modified_ghosted_entities.rend() ; ++i )
  {
    Entity entity                 = i->comm_info.entity;
    const EntityKey key           = i->comm_info.key;
    const int      remote_proc    = i->from_proc;
    const bool     local_owner    = m_bulkData.parallel_owner_rank(entity) == m_bulkData.parallel_rank() ;
    const bool remotely_destroyed = Deleted == i->state ;
    const bool remote_proc_is_owner = remote_proc == m_bulkData.parallel_owner_rank(entity);
    const bool isAlreadyDestroyed  = !m_bulkData.is_valid(entity);

    if ( local_owner ) { // Sending to 'remote_proc' for ghosting

      if ( remotely_destroyed ) {

        // remove from ghost-send list

        for ( size_t j = ghosting_count_minus_shared ; j>=1 ; --j) {
          m_bulkData.entity_comm_map_erase( key, EntityCommInfo( j , remote_proc ) );
        }
      }
      else {
        const bool shouldPromoteToShared = !isAlreadyDestroyed && i->remote_owned_closure==1 && key.rank() < stk::topology::ELEM_RANK;
        if ((shouldPromoteToShared || !isAlreadyDestroyed) && m_bulkData.state(entity)==Unchanged) {
          m_bulkData.set_state(entity, Modified);
        }

        if (shouldPromoteToShared) {
          m_bulkData.entity_comm_map_insert(entity, EntityCommInfo(BulkData::SHARED, remote_proc));
          promotingToShared.push_back(entity);
        }
      }
    }
    else if (remote_proc_is_owner) { // Receiving from 'remote_proc' for ghosting

      const bool hasBeenPromotedToSharedOrOwned = m_bulkData.owned_closure(entity);
      bool isAuraGhost = false;
      bool isCustomGhost = false;
      PairIterEntityComm pairIterEntityComm = m_bulkData.internal_entity_comm_map(entity);
      for(unsigned j=0; j<pairIterEntityComm.size(); ++j) {
        if (pairIterEntityComm[j].ghost_id == BulkData::AURA) {
          isAuraGhost = true;
        }
        else if (pairIterEntityComm[j].ghost_id > BulkData::AURA) {
          isCustomGhost = true;
        }
      }

      if ( isAuraGhost ) {
        if (!isAlreadyDestroyed && hasBeenPromotedToSharedOrOwned) {
          if (!remotely_destroyed) {
            m_bulkData.entity_comm_map_insert(entity, EntityCommInfo(BulkData::SHARED, remote_proc));
          }
          promotingToShared.push_back(entity);
        }
        m_bulkData.entity_comm_map_erase(key, m_bulkData.aura_ghosting());
      }

      if(!isAlreadyDestroyed) {
        const bool wasDestroyedByOwner = remotely_destroyed;
        const bool shouldDestroyGhost = wasDestroyedByOwner || (isAuraGhost && !isCustomGhost && !hasBeenPromotedToSharedOrOwned);
        const bool shouldRemoveFromGhosting = remotely_destroyed && !isAuraGhost && hasBeenPromotedToSharedOrOwned;

        if (shouldRemoveFromGhosting) {
          for ( size_t j = ghosting_count_minus_shared ; j >=1 ; --j ) {
            m_bulkData.entity_comm_map_erase( key, *m_bulkData.m_ghosting[j] );
          }
        }

        if ( shouldDestroyGhost ) {
          const bool was_ghost = true;
          m_bulkData.internal_destroy_entity_with_notification(entity, was_ghost);
        }

        m_bulkData.entity_comm_list_insert(entity);
      }
    }
  } // end loop on remote mod

  // Erase all ghosting communication lists for:
  // 1) Destroyed entities.
  // 2) Owned and modified entities.

  for ( EntityCommListInfoVector::const_reverse_iterator
        i = m_bulkData.internal_comm_list().rbegin() ; i != m_bulkData.internal_comm_list().rend() ; ++i) {

    Entity entity = i->entity;

    const bool locally_destroyed = !m_bulkData.is_valid(entity);
    const bool locally_owned_and_modified = locally_destroyed ? false :
      (Modified == m_bulkData.state(entity) && (m_bulkData.parallel_rank() == m_bulkData.parallel_owner_rank(entity)));

    if ( locally_destroyed ) {
      for ( size_t j = ghosting_count_minus_shared ; j >=1 ; --j ) {
        m_bulkData.entity_comm_map_erase( i->key, *m_bulkData.m_ghosting[j] );
      }
    }
    else if ( locally_owned_and_modified ) {
      m_bulkData.entity_comm_map_erase( i->key, m_bulkData.aura_ghosting() );
    }
  }

  if (!promotingToShared.empty()) {
    const bool shouldRemoveAuraPart = true;
    for(Entity entity : promotingToShared) {
      const bool remoteOwnerDestroyed = remote_owner_destroyed(m_bulkData.entity_key(entity), remotely_modified_ghosted_entities);
      process_changed_ownership_and_sharing(remoteOwnerDestroyed, entity, shouldRemoveAuraPart);
    }
    m_bulkData.add_comm_list_entries_for_entities(promotingToShared);
  }
}

void MeshModification::add_entity_to_same_ghosting(Entity entity, Entity connectedGhost)
{
  std::vector<EntityCommInfo> to_insert;
  for(PairIterEntityComm ec(m_bulkData.internal_entity_comm_map(connectedGhost)); ! ec.empty(); ++ec) {
    if (ec->ghost_id > BulkData::AURA) {
      to_insert.emplace_back(ec->ghost_id, ec->proc);
    }
  }
  if(!to_insert.empty()) {
    m_bulkData.entity_comm_list_insert(entity);
    for(const auto & entry : to_insert) {
      m_bulkData.entity_comm_map_insert(entity, EntityCommInfo(entry.ghost_id, entry.proc));
    }
  }
}

void MeshModification::internal_resolve_formerly_shared_entities(const EntityVector& entitiesNoLongerShared)
{
  for(Entity entity : entitiesNoLongerShared) {
    EntityVector ghostRelations = m_bulkData.get_upward_send_ghost_relations(entity);

    for(Entity ghost : ghostRelations) {
      add_entity_to_same_ghosting(entity, ghost);
    }
  }
}

void MeshModification::ensure_meta_data_is_committed()
{
  if (!m_bulkData.mesh_meta_data().is_commit())
  {
      m_bulkData.mesh_meta_data().commit();
  }
}

void MeshModification::reset_undeleted_entity_states_to_unchanged()
{
    for(unsigned i=0, iend=m_entity_states.size(); i<iend; ++i)
    {
        if(m_entity_states[i] != Deleted)
        {
            m_entity_states[i] = Unchanged;
        }
    }
}

}
}
}
