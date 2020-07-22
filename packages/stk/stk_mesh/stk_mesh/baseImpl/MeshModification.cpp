#include "MeshModification.hpp"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include <stk_mesh/baseImpl/EntityRepository.hpp>

namespace stk {
namespace mesh {
namespace impl {

bool MeshModification::modification_begin(const std::string description)
{
    parallel_machine_barrier( m_bulkData.parallel() );

    if (this->synchronized_count() == 0)
    {
        m_bulkData.mesh_meta_data().set_mesh_on_fields(&m_bulkData);
        m_bulkData.m_entity_repo->update_num_ranks(m_bulkData.mesh_meta_data().entity_rank_count());
    }

    if ( this->in_modifiable_state() ) return false ;

    if (this->synchronized_count() == 0)
    {
        this->ensure_meta_data_is_committed();

        if (m_bulkData.parallel_size() > 1) {
            verify_parallel_consistency( m_bulkData.mesh_meta_data() , m_bulkData.parallel() );
        }

        this->increment_sync_count();
    }
    else
    {
        this->increment_sync_count();
        this->reset_undeleted_entity_states_to_unchanged();
    }

    this->set_sync_state_modifiable();
    this->reset_shared_entity_changed_parts();

    for (FieldBase * stkField : m_bulkData.mesh_meta_data().get_fields()) {
      stkField->sync_to_host();
    }

    return true;
}

bool MeshModification::modification_end(modification_optimization opt)
{
    return this->internal_modification_end( opt );
}

bool MeshModification::resolve_node_sharing()
{
    return this->internal_resolve_node_sharing( MOD_END_SORT );
}

bool MeshModification::modification_end_after_node_sharing_resolution()
{
    return this->internal_modification_end_after_node_sharing_resolution( MOD_END_SORT );
}

bool MeshModification::internal_modification_end(modification_optimization opt)
{
    if(this->in_synchronized_state())
    {
        return false;
    }

    ThrowAssertMsg(impl::check_for_connected_nodes(m_bulkData)==0, "BulkData::modification_end ERROR, all entities with rank higher than node are required to have connected nodes.");

    ThrowAssertMsg(m_bulkData.add_fmwk_data() || impl::check_no_shared_elements_or_higher(m_bulkData)==0, "BulkData::modification_end ERROR, Sharing of entities with rank ELEMENT_RANK or higher is not allowed.");

    m_bulkData.m_entity_repo->clear_all_cache();

    if(m_bulkData.parallel_size() > 1)
    {
        // Resolve modification or deletion of shared entities
        // which can cause deletion of ghost entities.
        stk::mesh::EntityVector entitiesNoLongerShared;
        internal_resolve_shared_modify_delete(entitiesNoLongerShared);

        // Resolve modification or deletion of ghost entities
        // by destroying ghost entities that have been touched.
        m_bulkData.internal_resolve_ghosted_modify_delete(entitiesNoLongerShared);
        m_bulkData.update_comm_list_based_on_changes_in_comm_map();

        // Resolve creation of entities: discover sharing and set unique ownership.
        m_bulkData.internal_resolve_parallel_create();

        // Manoj: consider adding check_sharing_comm_maps here which is currently
        // in BulkDataTester in UnitTestModificationEnd.cpp

        // Resolve part membership for shared entities.
        // This occurs after resolving creation so created and shared
        // entities are resolved along with previously existing shared entities.
        m_bulkData.internal_resolve_shared_membership(entitiesNoLongerShared);

        // Regenerate the ghosting aura around all shared mesh entities.
        if(m_bulkData.is_automatic_aura_on())
        {
            m_bulkData.internal_regenerate_aura();
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

    m_bulkData.internal_finish_modification_end(opt);

    return true;
}

bool MeshModification::internal_resolve_node_sharing(modification_optimization opt)
{
    if(this->in_synchronized_state())
    {
        return false;
    }

    ThrowAssertMsg(impl::check_for_connected_nodes(m_bulkData)==0, "BulkData::modification_end ERROR, all entities with rank higher than node are required to have connected nodes.");

    ThrowAssertMsg(m_bulkData.add_fmwk_data() || impl::check_no_shared_elements_or_higher(m_bulkData)==0, "BulkData::modification_end ERROR, Sharing of entities with rank ELEMENT_RANK or higher is not allowed.");

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

bool MeshModification::internal_modification_end_after_node_sharing_resolution(modification_optimization opt)
{
    if(this->in_synchronized_state())
    {
        return false;
    }

    ThrowAssertMsg(impl::check_for_connected_nodes(m_bulkData)==0, "BulkData::modification_end ERROR, all entities with rank higher than node are required to have connected nodes.");

    ThrowAssertMsg(m_bulkData.add_fmwk_data() || impl::check_no_shared_elements_or_higher(m_bulkData)==0, "BulkData::modification_end ERROR, Sharing of entities with rank ELEMENT_RANK or higher is not allowed.");

    if(m_bulkData.parallel_size() > 1)
    {
        m_bulkData.internal_resolve_parallel_create_edges_and_faces();
        stk::mesh::EntityVector entitiesNoLongerShared;
        internal_resolve_shared_modify_delete(entitiesNoLongerShared);
        m_bulkData.internal_resolve_shared_membership(entitiesNoLongerShared);

        if(m_bulkData.is_automatic_aura_on())
        {
            m_bulkData.internal_regenerate_aura();
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

    m_bulkData.internal_finish_modification_end(opt);

    return true;
}

void MeshModification::change_entity_owner( const EntityProcVec & arg_change)
{
    ThrowRequireMsg(in_synchronized_state(), "BulkData::change_entity_owner() must not be called from within a modification cycle.");
    modification_optimization mod_optimization = MOD_END_SORT;
    modification_begin("change_entity_owner");
    m_bulkData.internal_change_entity_owner(arg_change, mod_optimization);
    m_bulkData.update_sharing_after_change_entity_owner();
    m_bulkData.internal_modification_end_for_change_entity_owner(mod_optimization);
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

void MeshModification::internal_resolve_shared_modify_delete(stk::mesh::EntityVector & entitiesNoLongerShared)
{
    ThrowRequireMsg(m_bulkData.parallel_size() > 1, "Do not call this in serial");

    stk::mesh::EntityProcVec entitiesToRemoveFromSharing;
    m_bulkData.delete_shared_entities_which_are_no_longer_in_owned_closure(entitiesToRemoveFromSharing);

    std::vector<stk::mesh::BulkData::EntityParallelState> remotely_modified_shared_entities;

    // Communicate entity modification state for shared entities
    // the resulting vector is sorted by entity and process.
    const bool communicate_shared = true;
    m_bulkData.communicate_entity_modification(communicate_shared, remotely_modified_shared_entities);

    // We iterate backwards to ensure that we hit the higher-ranking entities first.
    for(std::vector<stk::mesh::BulkData::EntityParallelState>::reverse_iterator
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
                m_bulkData.entity_comm_map_erase(key, EntityCommInfo(stk::mesh::BulkData::SHARED, remote_proc));
                if (!locally_destroyed && not_shared_remotely) {
                    entitiesToRemoveFromSharing.push_back(EntityProc(entity, remote_proc));
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
            const bool am_i_old_local_owner = m_bulkData.parallel_rank() == owner;

            if(remote_owner_destroyed)
            {
                m_bulkData.internal_establish_new_owner(entity);
            }

            const bool am_i_new_local_owner = m_bulkData.parallel_rank() == m_bulkData.parallel_owner_rank(entity);
            const bool did_i_just_become_owner = (!am_i_old_local_owner && am_i_new_local_owner );

            const bool is_entity_shared = !m_bulkData.internal_entity_comm_map_shared(key).empty();
            m_bulkData.internal_update_parts_for_shared_entity(entity, is_entity_shared, did_i_just_become_owner);
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

    m_bulkData.remove_entities_from_sharing(entitiesToRemoveFromSharing, entitiesNoLongerShared);
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
