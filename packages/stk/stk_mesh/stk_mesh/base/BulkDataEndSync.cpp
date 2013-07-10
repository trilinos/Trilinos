/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

/**
 * @author H. Carter Edwards
 */

//----------------------------------------------------------------------

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <set>
#include <vector>
#include <algorithm>

#include <stk_util/environment/ReportHandler.hpp>

#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/memory_util.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityCommDatabase.hpp>
#include <stk_mesh/base/Trace.hpp>

#ifdef STK_PROFILE_MEMORY
  #include <stk_mesh/baseImpl/Partition.hpp>
#endif

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

bool comm_mesh_verify_parallel_consistency(
  BulkData & M , std::ostream & error_log );

//----------------------------------------------------------------------

int BulkData::determine_new_owner( Entity entity ) const
{
  // We will decide the new owner by looking at all the processes sharing
  // this entity. The new owner will be the sharing process with lowest rank.

  // The local process is a candidate only if the entity is not destroyed.
  int new_owner = is_valid(entity) ? m_parallel_rank : ~0u;

  for ( PairIterEntityComm
        share = entity_comm_sharing(entity_key(entity)); ! share.empty() ; ++share ) {
    if ( share->proc < m_parallel_size &&
         ( new_owner < share->proc || m_parallel_size <= new_owner ) ) {
      new_owner = share->proc ;
    }
  }

  return new_owner ;
}

//----------------------------------------------------------------------

namespace {

// A method for quickly finding an entity within a comm list
EntityCommListInfo find_entity(const BulkData& mesh,
                               const std::vector<EntityCommListInfo>& entities,
                               const EntityKey& key)
{
  std::vector<EntityCommListInfo>::const_iterator lb_itr = std::lower_bound(entities.begin(), entities.end(), key);
  ThrowAssertMsg(lb_itr != entities.end() && lb_itr->key == key,
                 "Cannot find id: " << key.id() << " in comm-list" );
  return *lb_itr;
}

struct EntityParallelState {
  int                 from_proc;
  EntityState         state;
  EntityCommListInfo  comm_info;
  const BulkData* mesh;

  bool operator<(const EntityParallelState& rhs) const
  { return EntityLess(*mesh)(comm_info.entity, rhs.comm_info.entity); }
};

bool pack_entity_modification( const BulkData & mesh ,
                               const bool pack_shared ,
                               CommAll & comm )
{
  bool flag = false;

  const std::vector<EntityCommListInfo> & entity_comm = mesh.comm_list();

  for ( std::vector<EntityCommListInfo>::const_iterator
        i = entity_comm.begin() ; i != entity_comm.end() ; ++i ) {

    Entity entity = i->entity;
    EntityState status = mesh.is_valid(entity) ? mesh.state(entity) : Deleted;

    if ( status == Modified || status == Deleted ) {

      for ( PairIterEntityComm ec = mesh.entity_comm(i->key); ! ec.empty() ; ++ec ) {
        const bool shared = 0 == ec->ghost_id ;
        if ( pack_shared == shared ) {
          comm.send_buffer( ec->proc )
              .pack<EntityKey>( i->key )
              .pack<EntityState>( status );

          flag = true ;
        }
      }
    }
  }

  return flag ;
}

void communicate_entity_modification( const BulkData & mesh ,
                                      const bool shared ,
                                      std::vector<EntityParallelState > & data )
{
  CommAll comm( mesh.parallel() );
  const int p_size = comm.parallel_size();

  // Sizing send buffers:
  const bool local_mod = pack_entity_modification( mesh , shared , comm );

  // Allocation of send and receive buffers:
  const bool global_mod =
    comm.allocate_buffers( comm.parallel_size() / 4 , false , local_mod );

  if ( global_mod ) {
    const std::vector<EntityCommListInfo> & entity_comm = mesh.comm_list();

    // Packing send buffers:
    pack_entity_modification( mesh , shared , comm );

    comm.communicate();

    for ( int p = 0 ; p < p_size ; ++p ) {
      CommBuffer & buf = comm.recv_buffer( p );
      EntityKey key;
      EntityState state;

      while ( buf.remaining() ) {

        buf.unpack<EntityKey>( key )
           .unpack<EntityState>( state );

        // search through entity_comm, should only receive info on entities
        // that are communicated.
        EntityCommListInfo info = find_entity(mesh, entity_comm, key);
        EntityParallelState parallel_state = {p, state, info, &mesh};
        data.push_back( parallel_state );
      }
    }
  }

  std::sort( data.begin() , data.end() );
}

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

// Postconditions:
//  * DistributedIndex is updated based on entity creation/deletions in the
//    last modification cycle.
//  * Comm lists for shared entities are up-to-date.
//  * shared_new contains all entities that were modified/created on a
//    different process
void BulkData::internal_update_distributed_index(
  std::vector<Entity> & shared_new )
{
  Trace_("stk::mesh::BulkData::internal_update_distributed_index");

  parallel::DistributedIndex::KeyTypeVector
    local_created_or_modified ; // only store locally owned/shared entities

  // Iterate over all entities known to this process, putting
  // modified shared/owned entities in local_created_or_modified.
  size_t num_created_or_modified = 0;
  for ( impl::EntityRepository::const_iterator
        i = m_entity_repo.begin() ; i != m_entity_repo.end() ; ++i ) {

    Entity entity = i->second ;

    if ( state(entity) != Unchanged &&
         in_owned_closure( *this, entity , m_parallel_rank ) ) {
      // Has been changed and is in owned closure, may be shared
      ++num_created_or_modified;
    }
  }

  local_created_or_modified.reserve(num_created_or_modified);

  for ( impl::EntityRepository::const_iterator
        i = m_entity_repo.begin() ; i != m_entity_repo.end() ; ++i ) {

    Entity entity = i->second ;

    if ( state(entity) != Unchanged &&
         in_owned_closure(*this, entity , m_parallel_rank ) ) {
      // Has been changed and is in owned closure, may be shared
      local_created_or_modified.push_back( entity_key(entity) );
    }
  }

  {
    // Update distributed index. Note that the DistributedIndex only
    // tracks ownership and sharing information.
    parallel::DistributedIndex::KeyTypeVector::const_iterator begin = local_created_or_modified.begin();
    parallel::DistributedIndex::KeyTypeVector::const_iterator end = local_created_or_modified.end();
    m_entities_index.update_keys( begin, end );
  }

  if (parallel_size() > 1) {
    // Retrieve data regarding which processes use the local_created_or_modified
    // including this process.
    parallel::DistributedIndex::KeyProcVector
      global_created_or_modified ;
    m_entities_index.query_to_usage( local_created_or_modified ,
                                     global_created_or_modified );

    //------------------------------
    // Take the usage data and update the sharing comm lists
    {
      Entity entity = Entity();

      // Iterate over all global modifications to this entity, this vector is
      // sorted, so we're guaranteed that all modifications to a particular
      // entities will be adjacent in this vector.
      for ( parallel::DistributedIndex::KeyProcVector::iterator
              i =  global_created_or_modified.begin() ;
            i != global_created_or_modified.end() ; ++i ) {

        EntityKey key( static_cast<EntityKey::entity_key_t>(i->first) );
        int modifying_proc = i->second;

        if ( m_parallel_rank != modifying_proc ) {
          // Another process also created or updated this entity.

          // Only want to look up entities at most once
          if ( !is_valid(entity) || entity_key(entity) != key ) {
            // Have not looked this entity up by key
            entity = get_entity( key );

            shared_new.push_back( entity );
          }

          // Add the other_process to the entity's sharing info.
          entity_comm_insert(entity, EntityCommInfo( 0 /*sharing*/, modifying_proc ));
        }
      }
    }
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

// Enforce that shared entities must be in the owned closure:

void destroy_dependent_ghosts( BulkData & mesh , Entity entity )
{
  EntityRank entity_rank = mesh.entity_rank(entity);

  const EntityRank end_rank = mesh.mesh_meta_data().entity_rank_count();
  for (EntityRank irank = end_rank; irank > entity_rank; )
  {
    --irank;

    Entity const *irels_ri = mesh.end(entity, irank);
    Entity const *irels_b  = mesh.begin(entity, irank);
    for ( ; irels_ri != irels_b ; --irels_ri)
    {
      Entity e = *(irels_ri - 1);

      ThrowRequireMsg( !in_owned_closure(mesh, e , mesh.parallel_rank()),
          "Entity " << mesh.identifier(e) << " should not be in closure." );

      destroy_dependent_ghosts( mesh , e );
    }
  }

  mesh.destroy_entity( entity );
}

// Entities with sharing information that are not in the owned closure
// have been modified such that they are no longer shared.
// These may no longer be needed or may become ghost entities.
// There is not enough information so assume they are to be deleted
// and let these entities be re-ghosted if they are needed.

// Open question: Should an owned and shared entity that does not
// have an upward relation to an owned entity be destroyed so that
// ownership transfers to another process?

void resolve_shared_removed_from_owned_closure( BulkData & mesh )
{
  for ( std::vector<EntityCommListInfo>::const_reverse_iterator
        i =  mesh.comm_list().rbegin() ;
        i != mesh.comm_list().rend() ; ++i) {

    Entity entity = i->entity;

    if ( mesh.is_valid(entity) &&
         ! mesh.entity_comm_sharing(i->key).empty() &&
         ! in_owned_closure(mesh, entity , mesh.parallel_rank() ) ) {

      destroy_dependent_ghosts( mesh , entity );
    }
  }
}

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

void BulkData::internal_resolve_shared_modify_delete()
{
  Trace_("stk::mesh::BulkData::internal_resolve_shared_modify_delete");

  ThrowRequireMsg(parallel_size() > 1, "Do not call this in serial");

  resolve_shared_removed_from_owned_closure( *this );

  std::vector< EntityParallelState > remote_mod ;

  // Communicate entity modification state for shared entities
  // the resulting vector is sorted by entity and process.
  const bool communicate_shared = true ;
  communicate_entity_modification( *this , communicate_shared , remote_mod );

  // We iterate backwards over remote_mod to ensure that we hit the
  // higher-ranking entities first.
  for ( std::vector<EntityParallelState>::reverse_iterator
        i = remote_mod.rbegin(); i != remote_mod.rend() ; ) {

    Entity entity                = i->comm_info.entity;
    EntityKey key                = i->comm_info.key;
    int owner                    = i->comm_info.owner;
    const bool locally_destroyed = !is_valid(entity);
    bool remote_owner_destroyed  = false;

    // Iterate over all of this entity's remote changes
    for ( ; i != remote_mod.rend() && i->comm_info.entity == entity ; ++i ) {

      const int remote_proc    = i->from_proc ;
      const bool remotely_destroyed = Deleted == i->state ;

      // When a shared entity is remotely modified or destroyed
      // then the local copy is also modified.  This modification
      // status is applied to all related higher ranking entities.

      if ( ! locally_destroyed ) {
        this->modified( entity );
      }

      // A shared entity is being deleted on the remote process.
      // Remove it from the sharing communication list.
      // Ownership changes are processed later, but we'll need
      // to know if the remote owner destroyed the entity in order
      // to correctly resolve ownership (it is not sufficient to just
      // look at the comm list of the entity since there is no
      // guarantee that the comm list is correct or up-to-date).

      if ( remotely_destroyed ) {
        entity_comm_erase( key, EntityCommInfo(0,remote_proc) );

        // check if owner is destroying
        if ( owner == remote_proc ) {
          remote_owner_destroyed = true ;
        }
      }
    }

    // Have now processed all remote changes knowledge for this entity.

    PairIterEntityComm new_sharing = entity_comm_sharing(key);
    const bool   exists_somewhere = ! ( remote_owner_destroyed &&
                                        locally_destroyed &&
                                        new_sharing.empty() );

    // If the entity has been deleted everywhere, nothing left to do
    if ( exists_somewhere && !locally_destroyed ) {

      const bool old_local_owner = m_parallel_rank == parallel_owner_rank(entity);

      // If we are giving away ownership or the remote owner destroyed
      // the entity, then we need to establish a new owner
      if ( remote_owner_destroyed ) {

        const int new_owner = determine_new_owner( entity );

        const bool changed = this->set_parallel_owner_rank( entity, new_owner );
        if (changed) {
          internal_change_owner_in_comm_data(entity_key(entity), new_owner);
        }
        set_synchronized_count( entity, m_sync_count );
      }

      PartVector add_part , remove_part ;

      if ( new_sharing.empty() ) {
        // Is no longer shared, remove the shared part.
        remove_part.push_back(& m_mesh_meta_data.globally_shared_part());
      }

      const bool new_local_owner = m_parallel_rank == parallel_owner_rank(entity);

      const bool local_claimed_ownership =
        ( ! old_local_owner && new_local_owner );

      if ( local_claimed_ownership ) {
        // Changing remotely owned to locally owned
        add_part.push_back( & m_mesh_meta_data.locally_owned_part() );
      }

      if ( ! add_part.empty() || ! remove_part.empty() ) {
        internal_change_entity_parts( entity , add_part , remove_part );
      }
    } // if ( exists_somewhere )
  } // remote mod loop

  // Erase all sharing communication lists for Destroyed entities:
  for ( std::vector<EntityCommListInfo>::const_reverse_iterator
        i = comm_list().rbegin() ; i != comm_list().rend() ; ++i) {
    if ( !is_valid(i->entity) ) {
      // m_ghosting[0] is the SHARED communication
      entity_comm_erase( i->key, *m_ghosting[0] );
    }
  }
}


//----------------------------------------------------------------------
// Resolve modifications for ghosted entities:
// If a ghosted entity is modified or destroyed on the owning
// process then the ghosted entity must be destroyed.
//
// Post condition:
//  Ghosted entities of modified or deleted entities are destroyed.
//  Ghosted communication lists are cleared to reflect all deletions.

void BulkData::internal_resolve_ghosted_modify_delete()
{
  Trace_("stk::mesh::BulkData::internal_resolve_ghosted_modify_delete");

  ThrowRequireMsg(parallel_size() > 1, "Do not call this in serial");
  // Resolve modifications for ghosted entities:

  std::vector<EntityParallelState > remote_mod ;

  // Communicate entity modification state for ghost entities
  const bool communicate_shared = false ;
  communicate_entity_modification( *this , communicate_shared , remote_mod );

  const size_t ghosting_count = m_ghosting.size();

  std::vector< int > ghosting_change_flags( ghosting_count , 0 );

  // We iterate backwards over remote_mod to ensure that we hit the
  // higher-ranking entities first. This is important because higher-ranking
  // entities like element must be deleted before the nodes they have are
  // deleted.
  for ( std::vector<EntityParallelState>::reverse_iterator
        i = remote_mod.rbegin(); i != remote_mod.rend() ; ++i ) {
    Entity entity                 = i->comm_info.entity;
    const EntityKey key           = i->comm_info.key;
    const int      remote_proc    = i->from_proc;
    const bool     local_owner    = i->comm_info.owner == m_parallel_rank ;
    const bool remotely_destroyed = Deleted == i->state ;
    const bool locally_destroyed  = !is_valid(entity);

    if ( local_owner ) { // Sending to 'remote_proc' for ghosting

      if ( remotely_destroyed ) {

        // remove from ghost-send list

        for ( size_t j = ghosting_count ; j-- ; ) {
          if ( entity_comm_erase( key, EntityCommInfo( j , remote_proc ) ) ) {
            ghosting_change_flags[ j ] = true ;
          }
        }
      }

      // Remotely modified ghosts are ignored

    }
    else { // Receiving from 'remote_proc' for ghosting

      // Owner modified or destroyed, must locally destroy.

      for ( PairIterEntityComm ec = entity_comm(key); !ec.empty() ; ++ec ) {
        ghosting_change_flags[ ec->ghost_id ] = true ;
      }

      // This is a receive ghost so the only communication information
      // is the ghosting information, can clear it all out.
      entity_comm_clear(key);

      if ( ! locally_destroyed ) {

        // If mesh modification causes a ghost entity to become
        // a member of an owned-closure then do not automatically
        // destroy it.  The new sharing status will be resolved
        // in 'internal_resolve_parallel_create'.

        if ( ! in_owned_closure(*this, entity , m_parallel_rank ) ) {

          const bool destroy_entity_successful = destroy_entity(entity);
          ThrowRequireMsg(destroy_entity_successful,
              "Could not destroy ghost entity " << identifier(entity));
        }
      }
    }
  } // end loop on remote mod

  // Erase all ghosting communication lists for:
  // 1) Destroyed entities.
  // 2) Owned and modified entities.

  for ( std::vector<EntityCommListInfo>::const_reverse_iterator
        i = comm_list().rbegin() ; i != comm_list().rend() ; ++i) {

    Entity entity = i->entity;

    const bool locally_destroyed = !is_valid(entity);
    const bool locally_owned_and_modified = locally_destroyed ? false :
      Modified == state(entity) &&
      m_parallel_rank   == i->owner ;

    if ( locally_destroyed || locally_owned_and_modified ) {

      // m_ghosting[0] is the SHARED communication

      for ( size_t j = ghosting_count ; j-- ; ) {
        if ( entity_comm_erase( i->key, *m_ghosting[j] ) ) {
          ghosting_change_flags[ j ] = true ;
        }
      }
    }
  }

  std::vector< int > ghosting_change_flags_global( ghosting_count , 0 );

  all_reduce_sum( m_parallel_machine ,
                  & ghosting_change_flags[0] ,
                  & ghosting_change_flags_global[0] ,
                  ghosting_change_flags.size() );

  for ( unsigned ic = 0 ; ic < ghosting_change_flags_global.size() ; ++ic ) {
    if ( ghosting_change_flags_global[ic] ) {
      m_ghosting[ic]->m_sync_count = m_sync_count ;
    }
  }
}

//----------------------------------------------------------------------

// Postconditions:
//  * All shared entities have parallel-consistent owner
//  * Part membership of shared entities is up-to-date
//  * m_entity_comm is up-to-date
void BulkData::internal_resolve_parallel_create()
{
  Trace_("stk::mesh::BulkData::internal_resolve_parallel_create");

  ThrowRequireMsg(parallel_size() > 1, "Do not call this in serial");

  std::vector<Entity> shared_modified ;

  // Update the parallel index and
  // output shared and modified entities.
  internal_update_distributed_index( shared_modified );

  // ------------------------------------------------------------
  // Claim ownership on all shared_modified entities that I own
  // and which were not created in this modification cycle. All
  // sharing procs will need to be informed of this claim.
  CommAll comm_all( m_parallel_machine );

  for ( int phase = 0; phase < 2; ++phase ) {
    for ( std::vector<Entity>::iterator
            i = shared_modified.begin() ; i != shared_modified.end() ; ++i ) {
      Entity entity = *i ;
      if ( parallel_owner_rank(entity) == m_parallel_rank &&
           state(entity)  != Created ) {

        for ( PairIterEntityComm
              jc = entity_comm_sharing(entity_key(entity)) ; ! jc.empty() ; ++jc ) {
          comm_all.send_buffer( jc->proc ) .pack<EntityKey>( entity_key(entity) );
        }
      }
    }

    if (phase == 0) { //allocation phase
      comm_all.allocate_buffers( m_parallel_size / 4 );
    }
    else { // communication phase
      comm_all.communicate();
    }
  }

  for ( int p = 0 ; p < m_parallel_size ; ++p ) {
    CommBuffer & buf = comm_all.recv_buffer( p );
    EntityKey key ;
    while ( buf.remaining() ) {
      buf.unpack<EntityKey>( key );

      Entity entity = get_entity( key );

      // Set owner, will correct part membership later
      const bool changed = this->set_parallel_owner_rank( entity, p);
      if (changed) {
        internal_change_owner_in_comm_data(key, p);
      }
    }
  }

  // ------------------------------------------------------------
  // Update shared created entities.
  // - Revise ownership to selected processor
  // - Update sharing.
  // - Work backward so the 'in_owned_closure' function
  //   can evaluate related higher ranking entities.

  std::ostringstream error_msg ;
  int error_flag = 0 ;

  PartVector shared_part , owned_part ;
  shared_part.push_back( & m_mesh_meta_data.globally_shared_part() );
  owned_part.push_back(  & m_mesh_meta_data.locally_owned_part() );

  std::vector<Entity>::const_reverse_iterator iend = shared_modified.rend();
  for ( std::vector<Entity>::const_reverse_iterator
        i = shared_modified.rbegin() ; i != iend ; ++i) {

    Entity entity = *i ;

    if ( parallel_owner_rank(entity) == m_parallel_rank &&
         state(entity) == Created ) {

      // Created and not claimed by an existing owner

      const int new_owner = determine_new_owner( entity );

      const bool changed = this->set_parallel_owner_rank( entity, new_owner);
      if (changed) {
        internal_change_owner_in_comm_data(entity_key(entity), new_owner);
      }
    }

    if ( parallel_owner_rank(entity) != m_parallel_rank ) {
      // Do not own it and still have it.
      // Remove the locally owned, add the globally_shared
      set_synchronized_count( entity, m_sync_count);
      internal_change_entity_parts( entity , shared_part /*add*/, owned_part /*remove*/);
    }
    else if ( ! entity_comm_sharing(entity_key(entity)).empty() ) {
      // Own it and has sharing information.
      // Add the globally_shared
      internal_change_entity_parts( entity , shared_part /*add*/, PartVector() /*remove*/ );
    }
    else {
      // Own it and does not have sharing information.
      // Remove the globally_shared
      internal_change_entity_parts( entity , PartVector() /*add*/, shared_part /*remove*/);
    }

    // Newly created shared entity had better be in the owned closure
    if ( ! in_owned_closure(*this, entity , m_parallel_rank ) ) {
      if ( 0 == error_flag ) {
        error_flag = 1 ;
        error_msg
          << "\nP" << m_parallel_rank << ": " << " FAILED\n"
          << "  The following entities were declared on multiple processors,\n"
          << "  cannot be parallel-shared, and were declared with"
          << "  parallel-ghosting information. {\n";
      }
      error_msg << "    " << print_entity_key(m_mesh_meta_data, entity_key(entity));
      error_msg << " also declared on" ;
      for ( PairIterEntityComm ec = entity_comm_sharing(entity_key(entity)); ! ec.empty() ; ++ec ) {
        error_msg << " P" << ec->proc ;
      }
      error_msg << "\n" ;
    }
  }

  // Parallel-consistent error checking of above loop
  if ( error_flag ) { error_msg << "}\n" ; }
  all_reduce( m_parallel_machine , ReduceMax<1>( & error_flag ) );
  ThrowErrorMsgIf( error_flag, error_msg.str() );

  // ------------------------------------------------------------
  // Update m_entity_comm based on shared_modified

  const size_t n_old = m_entity_comm_list.size();

  m_entity_comm_list.reserve(m_entity_comm_list.size() + shared_modified.size());
  for (size_t i = 0, e = shared_modified.size(); i < e; ++i) {
    Entity entity = shared_modified[i];
    EntityCommListInfo new_comm = {entity_key(entity), entity, parallel_owner_rank(entity)};
    m_entity_comm_list.push_back(new_comm);
  }

  std::inplace_merge( m_entity_comm_list.begin() ,
                      m_entity_comm_list.begin() + n_old ,
                      m_entity_comm_list.end() );

  {
    std::vector<EntityCommListInfo>::iterator i =
      std::unique( m_entity_comm_list.begin() , m_entity_comm_list.end() );

    m_entity_comm_list.erase( i , m_entity_comm_list.end() );

    internal_sync_comm_list_owners();
  }
}

//----------------------------------------------------------------------
namespace {



}
//----------------------------------------------------------------------

bool BulkData::modification_end()
{
  Trace_("stk::mesh::BulkData::modification_end");

  bool return_value = internal_modification_end( true );

#ifdef STK_PROFILE_MEMORY

  std::cout << "Modification cycle: " << synchronized_count() << std::endl;
  profile_memory_usage<parallel::DistributedIndex>("Distributed Index", parallel(),parallel_rank());
  profile_memory_usage<BucketRelationTag>("Fixed Relation", parallel(),parallel_rank());
  profile_memory_usage<DynamicBucketRelationTag>("Dynamic Relation", parallel(),parallel_rank());
  profile_memory_usage<FieldBase>("Fields", parallel(),parallel_rank());
  profile_memory_usage<impl::Partition>("Partitions", parallel(),parallel_rank());
  profile_memory_usage<Bucket>("Buckets", parallel(),parallel_rank());

#endif

  return return_value;
}

#if 0

namespace {

// Very, very handy for debugging parallel resolution...

void print_comm_list( const BulkData & mesh , bool doit )
{
  if ( doit ) {
    std::ostringstream msg ;

    msg << std::endl ;

    for ( std::vector<EntityCommListInfo>::const_iterator
          i =  mesh.comm_list().begin() ;
          i != mesh.comm_list().end() ; ++i ) {

      Entity entity = i->entity;
      msg << "P" << mesh.parallel_rank() << ": " ;

      print_entity_key( msg , MetaData::get(mesh) , i->key );

      msg << " owner(" << i->owner << ")" ;

      if ( !entity.is_valid() ) { msg << " del" ; }
      else if ( Modified == entity.state() ) { msg << " mod" ; }
      else { msg << "    " ; }

      for ( PairIterEntityComm ec = mesh.comm_list(i->key); ! ec.empty() ; ++ec ) {
        msg << " (" << ec->ghost_id << "," << ec->proc << ")" ;
      }
      msg << std::endl ;
    }

    std::cout << msg.str();
  }
}

}

#endif

bool BulkData::internal_modification_end( bool regenerate_aura )
{
  Trace_("stk::mesh::BulkData::internal_modification_end");

  if ( m_sync_state == SYNCHRONIZED ) { return false ; }

  verify_SOA_relations();

  if (parallel_size() > 1) {
    // Resolve modification or deletion of shared entities
    // which can cause deletion of ghost entities.
    internal_resolve_shared_modify_delete();

    // Resolve modification or deletion of ghost entities
    // by destroying ghost entities that have been touched.
    internal_resolve_ghosted_modify_delete();

    // Resolution of shared and ghost modifications can empty
    // the communication information for entities.
    // If there is no communication information then the
    // entity must be removed from the communication list.
    {
      std::vector<EntityCommListInfo>::iterator i = m_entity_comm_list.begin();
      bool changed = false ;
      for ( ; i != m_entity_comm_list.end() ; ++i ) {
        if ( entity_comm(i->key).empty() ) {
          i->key = EntityKey();
          changed = true;
        }
      }
      if ( changed ) {
        i = std::remove_if( m_entity_comm_list.begin() ,
                            m_entity_comm_list.end() , IsInvalid() );
        m_entity_comm_list.erase( i , m_entity_comm_list.end() );
      }
    }

    // Resolve creation of entities: discover sharing and set unique ownership.
    internal_resolve_parallel_create();

    // Resolve part membership for shared entities.
    // This occurs after resolving creation so created and shared
    // entities are resolved along with previously existing shared entities.
    internal_resolve_shared_membership();

    // Regenerate the ghosting aura around all shared mesh entities.
    if ( regenerate_aura ) { internal_regenerate_shared_aura(); }

    // ------------------------------
    // Verify parallel consistency of mesh entities.
    // Unique ownership, communication lists, sharing part membership,
    // application part membership consistency.
#ifndef NDEBUG
    std::ostringstream msg ;
    bool is_consistent = true;
    is_consistent = comm_mesh_verify_parallel_consistency( *this , msg );
    ThrowErrorMsgIf( !is_consistent, msg.str() );
#endif
  }
  else {
    std::vector<Entity> shared_modified ;
    internal_update_distributed_index( shared_modified );
  }

  // ------------------------------
  // Now sort the bucket entities.
  // This does not change the entities, relations, or field data.
  // However, it insures that the ordering of entities and buckets
  // is independent of the order in which a set of changes were
  // performed.
  //
  //optimize_buckets combines multiple buckets in a bucket-family into
  //a single larger bucket, and also does a sort.
  //If optimize_buckets has not been requested, still do the sort.

  m_bucket_repository.optimize_buckets();

  // ------------------------------

  m_bucket_repository.internal_modification_end();
  m_sync_state = SYNCHRONIZED ;

  verify_SOA_relations();

  update_deleted_entities_container();

  return true ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

template <typename T>
T const* get_begin_itr(const Bucket& bucket, Bucket::size_type bucket_ordinal, EntityRank rank);

template <typename T>
T const* get_end_itr(const Bucket& bucket, Bucket::size_type bucket_ordinal, EntityRank rank);


template <>
Entity const* get_begin_itr<Entity>(const Bucket& bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.begin(bucket_ordinal, rank); }

template <>
Entity const* get_end_itr<Entity>(const Bucket& bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.end(bucket_ordinal, rank); }


template <>
ConnectivityOrdinal const* get_begin_itr<ConnectivityOrdinal>(const Bucket& bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.begin_ordinals(bucket_ordinal, rank); }

template <>
ConnectivityOrdinal const* get_end_itr<ConnectivityOrdinal>(const Bucket& bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.end_ordinals(bucket_ordinal, rank); }


template <>
Permutation const* get_begin_itr<Permutation>(const Bucket& bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.begin_permutations(bucket_ordinal, rank); }

template <>
Permutation const* get_end_itr<Permutation>(const Bucket& bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.end_permutations(bucket_ordinal, rank); }

//
// Because the connectivity API in Bucket is not templated...
//

template <typename T>
T const *get_begin_relation_data(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank);

template <typename T>
T const *get_end_relation_data(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank);

template <>
Entity const *get_begin_relation_data<Entity>(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.begin(bucket_ordinal, rank); }

template <>
Entity const *get_end_relation_data<Entity>(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.end(bucket_ordinal, rank); }

template <>
ConnectivityOrdinal const *get_begin_relation_data<ConnectivityOrdinal>(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.begin_ordinals(bucket_ordinal, rank); }

template <>
ConnectivityOrdinal const *get_end_relation_data<ConnectivityOrdinal>(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.end_ordinals(bucket_ordinal, rank); }

template <>
Permutation const *get_begin_relation_data<Permutation>(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.begin_permutations(bucket_ordinal, rank); }

template <>
Permutation const *get_end_relation_data<Permutation>(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.end_permutations(bucket_ordinal, rank); }

template <typename T>
void verify_relation_data(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{
  T const * old_it = get_begin_relation_data<T>(bucket, bucket_ordinal, rank);
  T const * old_it_end = get_end_relation_data<T>(bucket, bucket_ordinal, rank);

  T const * new_data_begin = get_begin_itr<T>(bucket, bucket_ordinal, rank);
  T const * new_data_end   = get_end_itr<T>(bucket, bucket_ordinal, rank);

  ThrowRequire( std::distance(new_data_begin, new_data_end ) ==
                std::distance(old_it, old_it_end ) );

  ThrowRequire( std::distance(new_data_begin, new_data_end) ==
                bucket.num_connectivity(bucket_ordinal, rank) );

  T const * new_it = new_data_begin;
  for (; old_it != old_it_end ; ++old_it , ++new_it) {
    T old_data = *old_it;
    T new_data = *new_it;
    ThrowRequire(old_data == new_data);
  }
}

}

void BulkData::verify_relations(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank) const
{
  verify_relation_data<Entity>(bucket, bucket_ordinal, rank);
  verify_relation_data<ConnectivityOrdinal>(bucket, bucket_ordinal, rank);
  if (bucket.has_permutation(rank)) {
    verify_relation_data<Permutation>(bucket, bucket_ordinal, rank);
  }
}

void BulkData::verify_SOA_relations() const
{

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

enum { PART_ORD_UNIVERSAL = 0 };
enum { PART_ORD_OWNED     = 1 };
enum { PART_ORD_SHARED    = 2 };

namespace {

void pack_induced_memberships( BulkData& bulk_data,
                               CommAll & comm ,
                               const std::vector<EntityCommListInfo> & entity_comm )
{
  for ( std::vector<EntityCommListInfo>::const_iterator
        i = entity_comm.begin() ; i != entity_comm.end() ; ++i ) {

    if ( bulk_data.in_shared( i->key , i->owner ) ) {
      // Is shared with owner, send to owner.

      OrdinalVector empty , induced ;

      induced_part_membership(bulk_data, i->entity , empty , induced );

      CommBuffer & buf = comm.send_buffer( i->owner );

      unsigned tmp = induced.size();

      buf.pack<unsigned>( tmp );

      for ( OrdinalVector::iterator
            j = induced.begin() ; j != induced.end() ; ++j ) {
        buf.pack<unsigned>( *j );
      }
    }
  }
}

void generate_send_list( BulkData& bulk_data,
                         const size_t sync_count ,
                         const int p_rank ,
                         const std::vector<EntityCommListInfo> & entity_comm ,
                               std::vector<EntityProc> & send_list )
{
  for ( std::vector<EntityCommListInfo>::const_iterator
        i = entity_comm.begin() ; i != entity_comm.end() ; ++i ) {

    if ( i->owner == p_rank &&
         bulk_data.synchronized_count(i->entity) == sync_count ) {

      for ( PairIterEntityComm ec = bulk_data.entity_comm(i->key); ! ec.empty(); ++ec ) {
        EntityProc tmp( i->entity , ec->proc );
        send_list.push_back( tmp );
      }
    }
  }

  {
    std::sort( send_list.begin() , send_list.end() , EntityLess(bulk_data) );
    std::vector<EntityProc>::iterator i =
      std::unique( send_list.begin() , send_list.end() );
    send_list.erase( i , send_list.end() );
  }
}

void pack_part_memberships( BulkData& meshbulk, CommAll & comm ,
                            const std::vector<EntityProc> & send_list )
{
  for ( std::vector<EntityProc>::const_iterator
        i = send_list.begin() ; i != send_list.end() ; ++i ) {

    Entity entity = i->first;

    std::pair<const unsigned *, const unsigned *>
      part_ord = meshbulk.bucket(entity).superset_part_ordinals();

    // I am the owner; therefore, the first three members are
    // universal, uses, and owns.  Don't send them.

    // I am the owner.  The first two memberships are
    // universal_part and locally_owned_part.  The third
    // membership may be globally_shared_part ;

    const unsigned count_all  = part_ord.second - part_ord.first ;
    const unsigned count_skip =
      ( 2 < count_all && part_ord.first[2] == PART_ORD_SHARED ) ? 3 : 2 ;

    const unsigned count_send = count_all - count_skip ;

    const unsigned * const start_send = part_ord.first + count_skip ;

    comm.send_buffer( i->second ).pack<EntityKey>( meshbulk.entity_key(entity) )
                                 .pack<unsigned>( count_send )
                                 .pack<unsigned>( start_send , count_send );
  }
}

}

//  Mesh entity membership changes must be synchronized among
//  processes that share mesh entities and propagated to
//  processes that ghost copies of the mesh entities.
//
//  Precondition: correct shared and ghosting lists.
//
//  Part memberships may have been added or removed
//  either explicitly or indirectly via entity relationships
//  being added or removed.

void BulkData::internal_resolve_shared_membership()
{
  Trace_("stk::mesh::BulkData::internal_resolve_shared_membership");

  ThrowRequireMsg(parallel_size() > 1, "Do not call this in serial");

  const MetaData & meta        = m_mesh_meta_data ;
  ParallelMachine p_comm       = m_parallel_machine ;
  const int p_rank             = m_parallel_rank ;
  const int p_size             = m_parallel_size ;
  const PartVector & all_parts = meta.get_parts();

  const Part & part_universal = meta.universal_part();
  const Part & part_owned  = meta.locally_owned_part();
  const Part & part_shared = meta.globally_shared_part();

  // Quick verification of part ordinal assumptions

  ThrowRequireMsg(PART_ORD_UNIVERSAL == part_universal.mesh_meta_data_ordinal(),
                  "Universal part ordinal is wrong, expected "
                  << PART_ORD_UNIVERSAL << ", got: "
                  << part_universal.mesh_meta_data_ordinal());

  ThrowRequireMsg(PART_ORD_OWNED == part_owned.mesh_meta_data_ordinal(),
                  "Owned part ordinal is wrong, expected "
                  << PART_ORD_OWNED << ", got: "
                  << part_owned.mesh_meta_data_ordinal());

  ThrowRequireMsg(PART_ORD_SHARED == part_shared.mesh_meta_data_ordinal(),
                  "Shared part ordinal is wrong, expected "
                  << PART_ORD_SHARED << ", got: "
                  << part_shared.mesh_meta_data_ordinal());

  //  Shared entities may have been modified due to relationship changes.
  //  Send just the current induced memberships from the sharing to
  //  the owning processes.
  {
    CommAll comm( p_comm );

    pack_induced_memberships( *this, comm , m_entity_comm_list );

    comm.allocate_buffers( p_size / 4 );

    pack_induced_memberships( *this, comm , m_entity_comm_list );

    comm.communicate();

    for ( std::vector<EntityCommListInfo>::iterator
          i = m_entity_comm_list.begin() ; i != m_entity_comm_list.end() ; ++i ) {

      if ( i->owner == p_rank ) {
        // Receiving from all sharing processes

        OrdinalVector empty , induced_parts , current_parts , remove_parts ;

        induced_part_membership(*this, i->entity , empty , induced_parts );

        for ( PairIterEntityComm ec = entity_comm_sharing(i->key) ; ! ec.empty() ; ++ec ) {

          CommBuffer & buf = comm.recv_buffer( ec->proc );

          unsigned count = 0 ; buf.unpack<unsigned>( count );
          for ( unsigned j = 0 ; j < count ; ++j ) {
            unsigned part_ord = 0 ; buf.unpack<unsigned>( part_ord );
            insert_ordinal( induced_parts , part_ord );
          }
        }

        // Remove any part that is an induced part but is not
        // in the induced parts list.

        this->bucket(i->entity).supersets( current_parts );

        OrdinalVector::const_iterator induced_parts_begin = induced_parts.begin(),
                                      induced_parts_end   = induced_parts.end();

        for ( OrdinalVector::iterator
              p = current_parts.begin() ; p != current_parts.end() ; ++p ) {
          if ( membership_is_induced( *meta.get_parts()[*p] , i->key.rank() ) &&
               ! contains_ordinal( induced_parts_begin, induced_parts_end , *p ) ) {
            remove_parts.push_back( *p );
          }
        }

        internal_change_entity_parts( i->entity, induced_parts, remove_parts );
      }
    }
  }

  //------------------------------
  // The owners have complete knowledge of memberships.
  // Send membership information to sync the shared and ghosted copies.
  // Only need to do this for entities that have actually changed.

  {
    std::vector<EntityProc> send_list ;

    generate_send_list( *this, m_sync_count, p_rank, m_entity_comm_list, send_list);

    CommAll comm( p_comm );

    pack_part_memberships( *this, comm , send_list );

    comm.allocate_buffers( p_size / 4 );

    pack_part_memberships( *this, comm , send_list );

    comm.communicate();

    for ( int p = 0 ; p < p_size ; ++p ) {
      CommBuffer & buf = comm.recv_buffer( p );
      while ( buf.remaining() ) {

        PartVector owner_parts , current_parts , remove_parts ;

        EntityKey key ; buf.unpack<EntityKey>( key );
        unsigned count = 0 ; buf.unpack<unsigned>( count );
        for ( unsigned j = 0 ; j < count ; ++j ) {
          unsigned part_ord = 0 ; buf.unpack<unsigned>( part_ord );
          insert( owner_parts , * all_parts[ part_ord ] );
        }

        // Any current part that is not a member of owners_parts
        // must be removed.

        Entity const entity = find_entity(*this, m_entity_comm_list, key).entity;

        this->bucket(entity).supersets( current_parts );

        for ( PartVector::iterator
              ip = current_parts.begin() ; ip != current_parts.end() ; ++ip ) {
          Part * const part = *ip ;
          const unsigned part_ord = part->mesh_meta_data_ordinal();
          if ( PART_ORD_UNIVERSAL != part_ord &&
               PART_ORD_OWNED     != part_ord &&
               PART_ORD_SHARED    != part_ord &&
               ! contain( owner_parts , *part ) ) {
            remove_parts.push_back( part );
          }
        }

//        std::cout << "P" << p_rank
//                  << " irsm calling internal_change_entity_parts(..) on entity "
//                  << entity << " to obey owner" << std::endl;

        internal_change_entity_parts( entity , owner_parts , remove_parts );
      }
    }
  }
}

} // namespace mesh
} // namespace stk
