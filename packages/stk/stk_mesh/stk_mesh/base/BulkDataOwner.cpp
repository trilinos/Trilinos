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

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <set>
#include <vector>
#include <algorithm>

#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityCommDatabase.hpp>
#include <stk_mesh/base/Trace.hpp>

namespace stk {
namespace mesh {

//----------------------------------------------------------------------

namespace {

// Given an entity, if it's a ghost, insert the closure of the entity
// into work_list.
void insert_closure_ghost(const BulkData& mesh, Entity const entity ,
                           const int proc_local ,
                           std::set<EntityKey> & work_list )
{
  if ( ! in_owned_closure( mesh, entity , proc_local ) ) {
    // This entity is a ghost, put it on the work_list
    // along with all ghosts in its closure

    const bool was_inserted = work_list.insert(mesh.entity_key(entity)).second;

    if ( was_inserted ) {
      // This ghost entity is new to the list, traverse its closure.

      const unsigned erank = mesh.entity_rank(entity);

      // Recurse over downward relations
      for (EntityRank irank = stk::topology::BEGIN_RANK; irank < erank; ++irank)
      {
        Entity const *irels_j = mesh.begin(entity, irank);
        Entity const *irels_e = mesh.end(entity, irank);
        for (; irels_j != irels_e; ++irels_j)
        {
          insert_closure_ghost(mesh, *irels_j, proc_local, work_list);
        }
      }
    }
  }
}

// Given an entity, insert the closures of every entity that has this entity
// in its closure. Only ghosts will be inserted.
void insert_transitive_ghost(const BulkData& mesh, Entity const entity ,
                              const int proc_local ,
                              std::set<EntityKey> & work_list )
{
  insert_closure_ghost(mesh, entity , proc_local , work_list );

  // Transitive:
  // If this entity is a member of another entity's closure
  // then that other entity is part of the traversal.

  const unsigned erank = mesh.entity_rank(entity);

  // Recurse over upward relations
  const EntityRank end_rank = mesh.mesh_meta_data().entity_rank_count();
  EntityVector temp_entities;
  Entity const* rels;
  int num_rels;
  for (EntityRank irank = erank + 1; irank < end_rank; ++irank)
  {
    if (mesh.connectivity_map().valid(erank, irank)) {
      num_rels = mesh.num_connectivity(entity, irank);
      rels     = mesh.begin(entity, irank);
    }
    else {
      num_rels = get_connectivity( mesh, entity, irank, temp_entities);
      rels     = &*temp_entities.begin();
    }

    for (int r = 0; r < num_rels; ++r)
    {
      insert_transitive_ghost(mesh, rels[r] , proc_local , work_list );
    }
  }
}

//----------------------------------------------------------------------

// Add EntityProc pairs to send_list for every entity in the closure of the
// entity in send_entry. All these entities will be sent to the same proc as
// the original send_entry.
void insert_closure_send(
  const BulkData &mesh,
  const EntityProc                  send_entry ,
  std::set<EntityProc,EntityLess> & send_list )
{
  ThrowRequireMsg( mesh.is_valid(send_entry.first),
                   "Cannot send destroyed entity");

  std::pair< std::set<EntityProc,EntityLess>::iterator , bool >
    result = send_list.insert( send_entry );

  if ( result.second ) {
    // First time this entity was inserted into the send_list.

    const unsigned erank  = mesh.entity_rank(send_entry.first);
    const Bucket &ebucket = mesh.bucket(send_entry.first);
    const Ordinal ebordinal = mesh.bucket_ordinal(send_entry.first);

    // Recurse over downward relations
    for (EntityRank irank = stk::topology::BEGIN_RANK; irank < erank; ++irank)
    {
      Entity const *rels_itr = ebucket.begin(ebordinal, irank);
      Entity const *rels_end= ebucket.end(ebordinal, irank);
      for (; rels_itr != rels_end; ++rels_itr)
      {
        const EntityProc rel_send_entry( *rels_itr, send_entry.second );
        insert_closure_send(mesh, rel_send_entry , send_list );
      }
    }
  }
}

//----------------------------------------------------------------------

bool member_of_owned_closure(const BulkData& mesh, const Entity e , const int p_rank )
{
  if (p_rank == mesh.parallel_owner_rank(e)) {
    return true;
  }

  const unsigned erank = mesh.entity_rank(e);
  const EntityRank end_rank = mesh.mesh_meta_data().entity_rank_count();

  // Any higher ranking entities locally owned?
  EntityVector temp_entities;
  Entity const* rels;
  int num_rels;
  for (EntityRank irank = end_rank - 1; irank > erank; --irank)
  {
    if (mesh.connectivity_map().valid(erank, irank)) {
      num_rels = mesh.num_connectivity(e, irank);
      rels     = mesh.begin(e, irank);
    }
    else {
      num_rels = get_connectivity( mesh, e, irank, temp_entities);
      rels     = &*temp_entities.begin();
    }

    for (int r = 0; r < num_rels; ++r) {
      if (p_rank == mesh.parallel_owner_rank(rels[r]) ||  member_of_owned_closure(mesh, rels[r], p_rank) ) {
        return true;
      }
    }
  }

  return false;
}

//----------------------------------------------------------------------

// Given a vector of local ownership changes, remove duplicates and
// sanity check.
void clean_and_verify_parallel_change(
  const BulkData & mesh ,
  std::vector<EntityProc> & local_change )
{
  const int             p_rank = mesh.parallel_rank();
  const int             p_size = mesh.parallel_size();
  const ParallelMachine p_comm = mesh.parallel();

  size_t error_count = 0 ;

  std::ostringstream error_msg ;

  // Order and eliminate redundancies:
  {
    std::vector<EntityProc>::iterator i = local_change.begin() ,
                                      j = local_change.end() ;
    std::sort( i , j , EntityLess(mesh) );
    i = std::unique( i , j );
    local_change.erase( i , j );
  }

  for ( std::vector<EntityProc>::iterator
        i = local_change.begin() ; i != local_change.end() ; ++i ) {
    std::vector<EntityProc>::iterator next = i+1 ;
    Entity const entity    = i->first ;
    const int new_owner = i->second ;

    // Verification:
    // 1) Cannot change the ownership of an entity you do not own
    // 2) New owner must exist
    // 3) Cannot grant ownership to two different owners

    const bool bad_null = !mesh.is_valid(entity);

    // Cannot change the ownership of an entity you do not own
    const bool bad_process_not_entity_owner = ! bad_null && mesh.parallel_owner_rank(entity) != p_rank ;

    // New owner must exist
    const bool bad_new_owner_does_not_exist = p_size <= new_owner ;

    // Cannot grant ownership to two different owners
    const bool bad_inconsistent_change = ! bad_null && next != local_change.end() && entity == next->first ;

    if ( bad_null ||
         bad_process_not_entity_owner ||
         bad_new_owner_does_not_exist ||
         bad_inconsistent_change)
    {
      ++error_count ;

      error_msg << "  P" << p_rank << ": " ;
      if ( bad_null ) { error_msg << " NULL ENTITY" ; }
      else { error_msg << mesh.identifier(entity); }
      if ( bad_process_not_entity_owner ) { error_msg << " NOT_CURRENT_OWNER" ; }
      if ( bad_new_owner_does_not_exist ) {
        error_msg << " BAD_NEW_OWNER( " << new_owner << " )" ;
      }
      if ( bad_inconsistent_change ) {
        error_msg << " CONFLICTING_NEW_OWNER( " << new_owner ;
        error_msg << " != " << next->second << " )" ;
      }
      error_msg << std::endl ;
    }
    else if ( new_owner == p_rank ) {
      // Eliminate non-changes
      i->first = Entity();
      i->second = 0;
    }
  }

  all_reduce( p_comm , ReduceSum<1>( & error_count ) );

  if ( error_count ) {
    all_write_string( p_comm , std::cerr , error_msg.str() );

    ThrowErrorMsg("Bad change ownership directives\n");
  }

  // Filter out non-changes (entity will be NULL
  {
    std::vector<EntityProc>::iterator i = local_change.begin(),
                                      j = local_change.end();
    i = std::remove( i , j , EntityProc(Entity(), 0) );
    local_change.erase( i , j );
  }
}

//----------------------------------------------------------------------
// Generate a parallel consistent list of ownership changes:
// 1) Shared entities (not owned but in closure of an owned entity),
// 2) Ghosted entities (not owned and not in closure of an owned entity), and
// 3) Parallel index.

void generate_parallel_change( const BulkData & mesh ,
                               const std::vector<EntityProc> & local_change ,
                                     std::vector<EntityProc> & shared_change ,
                                     std::vector<EntityProc> & ghosted_change )
{
  const int p_size = mesh.parallel_size();

  CommAll comm( mesh.parallel() );

  std::vector<int> procs ;

  // pack and communicate change owner information to all
  // processes that know about the entity
  for ( int phase = 0; phase < 2; ++phase) {
    for ( std::vector<EntityProc>::const_iterator
          ip = local_change.begin() ; ip != local_change.end() ; ++ip ) {
      Entity entity      = ip->first ;
      int new_owner = ip->second;
      mesh.comm_procs( mesh.entity_key(entity) , procs );
      for ( std::vector<int>::iterator
            j = procs.begin() ; j != procs.end() ; ++j )
      {
        comm.send_buffer( *j )
          .pack<EntityKey>( mesh.entity_key(entity) )
          .pack<int>(  new_owner );
      }
    }
    if (phase == 0) { // allocation phase
      comm.allocate_buffers( p_size / 4 , 0 );
    }
    else { // communication phase
      comm.communicate();
    }
  }

  // unpack communicated owner information into the
  // ghosted and shared change vectors.
  for ( int ip = 0 ; ip < p_size ; ++ip ) {
    CommBuffer & buf = comm.recv_buffer( ip );
    while ( buf.remaining() ) {
      EntityProc entry ;
      EntityKey key ;
      buf.unpack<EntityKey>( key )
         .unpack<int>( entry.second );

      entry.first = mesh.get_entity( key );

      if ( mesh.in_receive_ghost( mesh.entity_key(entry.first) ) ) {
        ghosted_change.push_back( entry );
      }
      else {
        shared_change.push_back( entry );
      }
    }
  }

  std::sort( shared_change.begin() , shared_change.end() , EntityLess(mesh) );
  std::sort( ghosted_change.begin() , ghosted_change.end() , EntityLess(mesh) );
}

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void BulkData::change_entity_owner( const std::vector<EntityProc> & arg_change )
{
  Trace_("stk::mesh::BulkData::change_entity_owner");
  DiagIf(LOG_ENTITY, "arg_change: " << arg_change);

  const MetaData  & meta = m_mesh_meta_data ;
  const int       p_rank = m_parallel_rank ;
  const int       p_size = m_parallel_size ;
  ParallelMachine p_comm = m_parallel_machine ;

  //------------------------------
  // Verify the input changes, generate a clean local change list, and
  // generate the remote change list so that all processes know about
  // pending changes.

  std::vector<EntityProc> local_change( arg_change );

  // Parallel synchronous clean up and verify the requested changes:
  clean_and_verify_parallel_change( *this , local_change );

  //----------------------------------------
  // Parallel synchronous determination of changing shared and ghosted.

  // The two vectors below will contain changes to ghosted and shared
  // entities on this process coming from change-entity-owner requests
  // on other processes.
  std::vector<EntityProc> ghosted_change ;
  std::vector<EntityProc> shared_change ;

  generate_parallel_change( *this , local_change ,
                            shared_change , ghosted_change );

  //------------------------------
  // Have enough information to delete all effected ghosts.
  // If the closure of a ghost contains a changing entity
  // then that ghost must be deleted.
  // Request that all ghost entities in the closure of the ghost be deleted.

  typedef std::set<EntityProc,EntityLess> EntityProcSet;

  // Compute the closure of all the locally changing entities
  std::set<EntityProc,EntityLess> send_closure(EntityLess(*this)) ;
  for ( std::vector<EntityProc>::iterator
        i = local_change.begin() ; i != local_change.end() ; ++i ) {
    insert_closure_send(*this, *i , send_closure );
  }

  // Calculate all the ghosts that are impacted by the set of ownership
  // changes. We look at ghosted, shared, and local changes looking for ghosts
  // that are either in the closure of the changing entity, or have the
  // changing entity in their closure. All modified ghosts will be removed.
  {
    std::set<EntityKey> modified_ghosts;

    for ( std::vector<EntityProc>::const_iterator
          i = ghosted_change.begin() ; i != ghosted_change.end() ; ++i ) {
      insert_transitive_ghost(*this, i->first , m_parallel_rank , modified_ghosts );
    }

    for ( std::vector<EntityProc>::const_iterator
          i = shared_change.begin() ; i != shared_change.end() ; ++i ) {
      insert_transitive_ghost(*this, i->first , m_parallel_rank , modified_ghosts );
    }

    for ( EntityProcSet::iterator
          i = send_closure.begin() ; i != send_closure.end() ; ++i ) {
      insert_transitive_ghost(*this, i->first , m_parallel_rank , modified_ghosts );
    }

    // The ghosted change list will become invalid
    ghosted_change.clear();

    std::vector<EntityProc> empty_add ;
    std::vector<EntityKey>  remove_modified_ghosts( modified_ghosts.begin() ,
                                                    modified_ghosts.end() );

    // Skip 'm_ghosting[0]' which is the shared subset.
    for ( std::vector<Ghosting*>::iterator
          ig = m_ghosting.begin() + 1; ig != m_ghosting.end() ; ++ig ) {
      // parallel synchronous:
      internal_change_ghosting( **ig , empty_add , remove_modified_ghosts );
    }
  }

  //------------------------------
  // Consistently change the owner on all processes.
  // 1) The local_change list is giving away ownership.
  // 2) The shared_change may or may not be receiving ownership

  {
    PartVector owned;
    owned.push_back(& meta.locally_owned_part());

    for ( std::vector<EntityProc>::iterator
          i = local_change.begin() ; i != local_change.end() ; ++i ) {
      // Giving ownership, change the parts first and then
      // the owner rank to pass the ownership test.
      Entity entity = i->first;

      change_entity_parts( entity , PartVector() , owned );

      const bool changed = this->set_parallel_owner_rank( entity, i->second );
      if (changed) {
        internal_change_owner_in_comm_data(entity_key(entity), i->second);
      }
    }

    for ( std::vector<EntityProc>::iterator
          i = shared_change.begin() ; i != shared_change.end() ; ++i ) {
      Entity entity = i->first;
      const bool changed = this->set_parallel_owner_rank( entity, i->second );
      if (changed) {
        internal_change_owner_in_comm_data(entity_key(entity), i->second);
      }
      if ( p_rank == i->second ) { // I receive ownership
        change_entity_parts( entity , owned , PartVector() );
      }
    }
  }

  //------------------------------
  // Send entities, along with their closure, to the new owner processes
  {
    std::ostringstream error_msg ;
    int error_count = 0 ;

    CommAll comm( p_comm );

    EntityVector unique_list_of_send_closure;
    unique_list_of_send_closure.reserve(send_closure.size());

    for ( std::set<EntityProc,EntityLess>::iterator
          i = send_closure.begin() ; i != send_closure.end() ; ++i ) {
      CommBuffer & buffer = comm.send_buffer( i->second );
      Entity entity = i->first;
      pack_entity_info(*this, buffer , entity );
      pack_field_values(*this, buffer , entity );

      if (unique_list_of_send_closure.empty() || entity_key(unique_list_of_send_closure.back()) != entity_key(entity)) {
        unique_list_of_send_closure.push_back(entity);
      }
    }

    comm.allocate_buffers( p_size / 4 );

    for ( std::set<EntityProc,EntityLess>::iterator
          i = send_closure.begin() ; i != send_closure.end() ; ++i ) {
      CommBuffer & buffer = comm.send_buffer( i->second );
      Entity entity = i->first;
      pack_entity_info(*this, buffer , entity );
      pack_field_values(*this, buffer , entity );
    }

    comm.communicate();

    for ( int p = 0 ; p < p_size ; ++p ) {
      CommBuffer & buf = comm.recv_buffer(p);
      while ( buf.remaining() ) {
        PartVector parts ;
        std::vector<Relation> relations ;
        EntityKey key ;
        int owner = ~0u ;

        unpack_entity_info( buf, *this, key, owner, parts, relations );

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

        std::pair<Entity ,bool> result =
          m_entity_repo.internal_create_entity( key );

        Entity entity = result.first;

        // The entity was copied and not created.

        internal_change_entity_parts( entity , parts , PartVector() );

        log_created_parallel_copy( entity );

        const bool changed = this->set_parallel_owner_rank( entity, owner );
        if (changed) {
          internal_change_owner_in_comm_data(entity_key(entity), owner);
        }

        declare_relation( entity , relations );

        if ( ! unpack_field_values(*this, buf , entity , error_msg ) ) {
          ++error_count ;
        }
      }
    }

    all_reduce( p_comm , ReduceSum<1>( & error_count ) );
    ThrowErrorMsgIf( error_count, error_msg.str() );

    // Any entity that I sent and is not in an owned closure is deleted.
    // The owned closure will be effected by received entities, so can
    // only clean up after the newly owned entities have been received.
    // Destroy backwards so as not to invalidate closures in the process.

    {
      for ( EntityVector::reverse_iterator i = unique_list_of_send_closure.rbegin() ;
            i != unique_list_of_send_closure.rend() ;
            ++i) {
        if ( ! member_of_owned_closure(*this, *i , p_rank ) ) {
          ThrowRequireMsg( destroy_entity( *i ),
                           "Failed to destroy entity " << identifier(*i) );
        }
      }
    }

    send_closure.clear(); // Has been invalidated
  }
}

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

