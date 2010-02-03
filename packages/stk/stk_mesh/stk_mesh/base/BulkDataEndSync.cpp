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
#include <stk_mesh/base/EntityComm.hpp>

namespace stk {
namespace mesh {

//----------------------------------------------------------------------

typedef std::pair< EntityKey , unsigned > KeyProc ;

namespace {

struct KeyProcLess {
  bool operator()( const KeyProc & lhs , const KeyProc & rhs ) const
    { return lhs.first != rhs.first ? lhs.first < rhs.first
                                    : lhs.second < rhs.second ; }
};

//----------------------------------------------------------------------

struct EntityProcClean {
  enum { flag = ~0u };

  bool operator()( const KeyProc & value ) const
  { return value.second == flag ; }

  bool operator()( const EntityProc & value ) const
  { return value.first == NULL || value.second == flag ; }

  bool operator()( Entity * const value ) const
  { return value == NULL ; }

  static void mark( KeyProc & kp ) { kp.second = flag ; }

  static void mark( EntityProc & ep ) { ep.second = flag ; }

  static void mark( Entity * & ep ) { ep = NULL ; }

  static bool clean( std::vector<KeyProc> & kpv )
  {
    std::vector<KeyProc>::iterator i = kpv.begin() , j = kpv.end();
    i = std::remove_if( i , j , EntityProcClean() );
    const bool changed = i != j ;
    if ( changed ) { kpv.erase( i , j ); }
    return changed ;
  }

  static bool clean( std::vector<EntityProc> & epv )
  {
    std::vector<EntityProc>::iterator i = epv.begin() , j = epv.end();
    i = std::remove_if( i , j , EntityProcClean() );
    const bool changed = i != j ;
    if ( changed ) { epv.erase( i , j ); }
    return changed ;
  }

  static bool clean( std::vector<Entity*> & ev )
  {
    std::vector<Entity*>::iterator i = ev.begin() , j = ev.end();
    i = std::remove_if( i , j , EntityProcClean() );
    const bool changed = i != j ;
    if ( changed ) { ev.erase( i , j ); }
    return changed ;
  }
};

//----------------------------------------------------------------------

unsigned determine_new_owner(
  Entity & entity , const unsigned p_rank , const unsigned p_size )
{
  const bool keeping = entity.bucket().capacity();

  unsigned new_owner = keeping ? p_rank : ~0u ;

  for ( PairIterEntityProc
        share = entity.sharing(); ! share.empty() ; ++share ) {
    if ( share->second < p_size &&
         ( new_owner < share->second || p_size <= new_owner ) ) {
      new_owner = share->second ;
    }
  }

  return new_owner ;
}

//----------------------------------------------------------------------

void send_owned_shared_ghosted( const unsigned                  p_rank ,
                                const std::vector<Ghosting*>  & all_ghosts ,
                                const std::vector<Entity*>    & entities ,
                                      std::vector<EntityProc> & send )
{
  // 'entities' must be sorted and unique:

  for ( std::vector<Entity*>::const_iterator
        i = entities.begin(); i != entities.end() ; ++i ) {

    const size_t n = send.size();

    Entity * const entity = *i ;
    const unsigned p_owner = entity->owner_rank();

    // Shared info is sorted, unique, and not local:

    PairIterEntityProc share = entity->sharing();
    send.insert( send.end() , share.begin() , share.end() );

    if ( p_owner != p_rank ) {
      // If owner is not local then send to owner:

      EntityProc tmp( entity , p_owner );
      send.push_back( tmp );
    }
    else {

      // Owned entity may be ghosted elsewhere:

      for ( std::vector<Ghosting*>::const_iterator
            j = all_ghosts.begin() ; j != all_ghosts.end() ; ++j ) {
 
        const std::vector<EntityProc> & ghost_send_list = (*j)->send();
 
        std::vector<EntityProc>::const_iterator
          k_beg = ghost_send_list.begin(),
          k_end = ghost_send_list.end();

        k_beg = std::lower_bound( k_beg , k_end , *entity , EntityLess() );
 
        for ( k_end = k_beg ;
              k_end != ghost_send_list.end() &&
              k_end->first == entity ; ++k_end );
 
        send.insert( send.end() , k_beg , k_end );
      }
    }

    // Merge and unique the ghosted information:

    std::vector<EntityProc>::iterator ii = send.begin() + n ;
    std::vector<EntityProc>::iterator jj = send.end();
    std::sort( ii , jj , EntityLess() );
    ii = std::unique( ii , jj );
    send.erase( ii , jj );
  }
}

//----------------------------------------------------------------------

bool communicate_entity_vector( const BulkData                & mesh ,
                                const std::vector<EntityProc> & send ,
                                      std::vector<EntityProc> & recv )
{
  const unsigned p_size = mesh.parallel_size();
  CommAll comm( mesh.parallel() );

  for ( std::vector<EntityProc>::const_iterator
        i = send.begin(); i != send.end() ; ++i ) {
    comm.send_buffer( i->second ).pack<EntityKey>( i->first->key() );
  }

  const bool local  = ! send.empty();
  const bool global = comm.allocate_buffers( p_size/4, false, local );

  if ( ! global ) { return false ; }

  for ( std::vector<EntityProc>::const_iterator
        i = send.begin(); i != send.end() ; ++i ) {
    comm.send_buffer( i->second ).pack<EntityKey>( i->first->key() );
  }

  comm.communicate();

  for ( unsigned p = 0 ; p < p_size ; ++p ) {
    CommBuffer & buf = comm.recv_buffer( p );

    EntityProc ep ; ep.second = p ;

    while ( buf.remaining() ) {
      EntityKey key ; buf.unpack<EntityKey>( key );

      ep.first = mesh.get_entity( entity_type(key) , entity_id(key) );

      recv.push_back( ep );
    }
  }

  std::sort( recv.begin() , recv.end() , EntityLess() );

  return true ;
}

}

//----------------------------------------------------------------------

void BulkData::internal_resolve_destroy_ghosted(
  const std::vector<Entity*>    & del_local ,
  const std::vector<EntityProc> & del_remote )
{
  const EntityLess compare ;

  //  Remove destroyed owned and ghosts from
  //  the ghost's send and receive lists

  const unsigned count = m_ghosting.size();
  std::vector<int> local_flags(  count , 0 );
  std::vector<int> global_flags( count , 0 );
  std::vector<Entity*> ghosts_to_destroy ;

  for ( unsigned ic = 0 ; ic < count ; ++ic ) {

    std::vector<EntityProc> & send = m_ghosting[ic]->m_send ;
    std::vector<Entity*>    & recv = m_ghosting[ic]->m_recv ;

    std::vector<Entity*   >::const_iterator ilocal  = del_local .begin();
    std::vector<EntityProc>::const_iterator iremote = del_remote.begin();

    // Entities in the send list are locally owned

    for ( std::vector<EntityProc>::iterator
          i = send.begin() ; i != send.end() ; ) {

      Entity * const entity = i->first ;

      for ( ; ilocal != del_local.end() &&
              compare( *ilocal , *entity ) ; ++ilocal );

      if ( ilocal != del_local.end() && entity == *ilocal ) {
        // Local deletion of an owned entity,
        // mark all send entries of this entity for deletion.

        for ( ; i != send.end() && entity == i->first ; ++i ) {
          EntityProcClean::mark( *i );
        }
        ++ilocal ;
      }
      else {
        for ( ; iremote != del_remote.end() &&
                compare( *iremote , *entity ) ; ++iremote );

        if ( iremote != del_remote.end() && *iremote == *i ) {
          // Remote deletion of a sent ghost, don't send it anymore

          EntityProcClean::mark( *i );
          ++iremote ;
        }
        ++i ;
      }
    }

    ilocal  = del_local .begin();
    iremote = del_remote.begin();

    // Entities in the receive list are not locally owned

    for ( std::vector<Entity*>::iterator
          i = recv.begin() ; i != recv.end() ; ++i ) {

      Entity * const entity = *i ;

      for ( ; ilocal != del_local.end() &&
              compare( *ilocal , *entity ) ; ++ilocal );

      if ( ilocal != del_local.end() && entity == *ilocal ) {
        // Local deletion of a ghosted entity:
        EntityProcClean::mark( *i );
        ++ilocal ;
      }
      else {
        for ( ; iremote != del_remote.end() &&
                compare( *iremote , *entity ) ; ++iremote );

        if ( iremote != del_remote.end() &&
             iremote->first  == entity &&
             iremote->second == entity->owner_rank() ) {
          // Remote deletion by owner of a ghosted entity:
          ghosts_to_destroy.push_back( *i );
          EntityProcClean::mark( *i );
          ++iremote ;
        }
      }
    }

    const bool change_send = EntityProcClean::clean( send );
    const bool change_recv = EntityProcClean::clean( recv );
    local_flags[ic] = change_send || change_recv ;
  }

  // Destroy received ghosts of destroyed entities.
  // Destroy from highest to lowest ranking for proper
  // removal of relations.
  {
    std::sort( ghosts_to_destroy.begin() ,
               ghosts_to_destroy.end() , EntityLess() );
    std::vector<Entity*>::iterator
      i = std::unique( ghosts_to_destroy.begin() , ghosts_to_destroy.end() );
    ghosts_to_destroy.erase( i , ghosts_to_destroy.end() );
    while ( ! ghosts_to_destroy.empty() ) {
      internal_destroy_entity( ghosts_to_destroy.back() );
      ghosts_to_destroy.pop_back();
    }
  }

  all_reduce_sum( m_parallel_machine ,
                  & local_flags[0] , & global_flags[0] , count );

  for ( unsigned ic = 0 ; ic < count ; ++ic ) {
    if ( global_flags[ic] ) {
      m_ghosting[ic]->m_sync_count = m_sync_count ;
    }
  }
}

//----------------------------------------------------------------------

void BulkData::internal_resolve_destroy_shared(
  const std::vector<Entity*>    & del_local ,
  const std::vector<EntityProc> & del_remote )
{
  // Update all of the sharing information and set new owners as needed.

  for ( std::vector<EntityProc>::const_iterator
        i =  del_remote.begin(); i != del_remote.end(); ++i ) {

    Entity * const entity = i->first ;

    PairIterEntityProc shares = entity->sharing();

    if ( ! shares.empty() ) { // Is a shared entity, not a ghost entity

      // A shared entity that is being deleted on the remote process.

      {
        const std::vector<EntityProc>::const_iterator
          shares_all_begin = m_shares_all.begin();

        // Get the non-const iterator for this entry so that it can be updated.
        // Invalidate the shared entry for this process,
        // entry will be removed later as a batch operation.
        std::vector<EntityProc>::iterator ish =
          m_shares_all.begin() +
            std::distance( shares_all_begin , shares.begin() );

        for ( ; ish->second != i->second ; ++ish );

        EntityProcClean::mark( *ish );
      }

      if ( entity->owner_rank() == i->second ) {
        // Remote owner is deleting a shared entity.
        // Have to determine new owner from remaining sharing processes.

        const unsigned new_owner =
          determine_new_owner( *entity , m_parallel_rank , m_parallel_size );

        entity->m_owner_rank = new_owner ;
        entity->m_sync_count = m_sync_count ;

        if ( new_owner == m_parallel_rank ) {
          // Changing remotely owned to locally owned
          PartVector add( 1 );
          add[0] = & m_mesh_meta_data.locally_owned_part();
          internal_change_entity_parts( *entity , add , PartVector() );
        }
      }
    }
  }

  //------------------------------
  // If an owned and shared entity is being deleted
  // and a sharing process is going to keep that entity
  // then transfer ownership.  The update of the 
  // parallel index must be informed.

  for ( std::vector<Entity*>::const_iterator
        i = del_local.begin(); i != del_local.end() ; ++i ) {

    Entity * const entity = *i ;

    const unsigned new_owner =
      determine_new_owner( *entity , m_parallel_rank , m_parallel_size );

    entity->m_owner_rank = new_owner ;
    entity->m_sync_count = m_sync_count ;

    PairIterEntityProc shares = entity->sharing();

    if ( ! shares.empty() ) {
      const std::vector<EntityProc>::const_iterator
        shares_all_begin = m_shares_all.begin();

      // Is shared and have extracted all information from the sharing.
      // Now mark all of the sharing for removal.

      std::vector<EntityProc>::iterator ish =
        m_shares_all.begin() +
          std::distance( shares_all_begin , shares.begin() );

      while ( ish != m_shares_all.end() )
      {
        if ( entity != ish->first ) break;
          EntityProcClean::mark( *ish );
        ish++; 
      }
    }
  }

  for ( std::vector<EntityProc>::iterator
        ish =  m_shares_all.begin() ; 
        ish != m_shares_all.end() ; ++ish ) {
    PairIterEntityProc & pit = ish->first->m_sharing ;
    pit.second = pit.first = std::vector< EntityProc >::const_iterator();
  }

  EntityProcClean::clean( m_shares_all );

  {
    const std::vector<EntityProc>::iterator es = m_shares_all.end();
          std::vector<EntityProc>::iterator js = m_shares_all.begin();
 
    while ( js != es ) {
      const std::vector<EntityProc>::iterator is = js ;
      for ( ; js != es && is->first == js->first ; ++js );
      PairIterEntityProc & pit = is->first->m_sharing ;
      pit.first  = is ;
      pit.second = js ;
    }
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void BulkData::internal_update_parallel_index(
  std::vector<EntityProc> & shared_new )
{
  static const char method[] =
    "stk::mesh::BulkData::internal_update_parallel_index" ;

  const size_t new_count = m_new_entities.size();
  const size_t del_count = m_del_entities.size();

  std::vector< parallel::DistributedIndex::KeyType >
    new_entities_keys(new_count),
    del_entities_keys(del_count);

  std::vector< parallel::DistributedIndex::KeyProc >
    new_entities_keyprocs ;

  for ( size_t i = 0 ; i < new_count ; ++i ) {
    new_entities_keys[i] = m_new_entities[i]->key().raw_key();
  }

  for ( size_t i = 0 ; i < del_count ; ++i ) {
    del_entities_keys[i] = m_del_entities[i]->key().raw_key();
  }

  //------------------------------
  // Verify created did not already exist:
  {
    m_entities_index.query( new_entities_keys , new_entities_keyprocs );

    unsigned error_count = new_entities_keyprocs.size();

    all_reduce( m_parallel_machine , ReduceSum<1>( & error_count ) );

    if ( error_count ) {
      std::ostringstream msg ;
      msg << method ;
      msg << " FAILED: Declare entities already existed" ;

      for ( std::vector< parallel::DistributedIndex::KeyProc >::iterator
            i =  new_entities_keyprocs.begin() ;
            i != new_entities_keyprocs.end() ; ++i ) {
        EntityKey key( & i->first );
        msg << std::endl << "  ( " ;
        print_entity_key( msg , m_mesh_meta_data , key );
        msg << "  Already exists on P" << i->second << " )" ;
      }
 
      throw std::runtime_error( msg.str() );
    }
  }

  //------------------------------
  // Update distributed index.

  m_entities_index.update_keys( new_entities_keys , del_entities_keys );

  //------------------------------
  // Inform creating processes:

  m_entities_index.query( new_entities_keys , new_entities_keyprocs );

  {
    EntityProc entry( (Entity*) NULL , 0 );

    for ( std::vector< parallel::DistributedIndex::KeyProc >::iterator
          i =  new_entities_keyprocs.begin() ;
          i != new_entities_keyprocs.end() ; ++i ) {

      EntityKey key( & i->first );

      if ( (int) m_parallel_rank != i->second ) {
        // Another process also created this entity.

        if ( entry.first == NULL || entry.first->key() != key ) {
          // Have not looked this entity up by key 
          const EntitySet::iterator j = m_entities.find( key );
          entry.first = j->second ;
        }
        entry.second = i->second ;
        // Add the ( entity , other_process ) entry to the shared list.
        shared_new.push_back( entry );
      }
    }
  }
}

void BulkData::internal_resolve_created_shared(
  const std::vector<EntityProc> & new_shared )
{
  // Revise ownership to highest ranking processor:
  {
    Part * const owns_part = & m_mesh_meta_data.locally_owned_part();
     
    const std::vector<Part*> add_parts , remove_parts( 1 , owns_part );
     
    for ( std::vector<EntityProc>::const_iterator 
          ip = new_shared.begin() ; ip != new_shared.end() ; ) {
     
      EntityProc tmp = *ip ;
     
      for ( ; ip != new_shared.end() && tmp.first == ip->first ; ++ip ) {
        tmp.second = ip->second ;
      }
     
      if ( m_parallel_rank < tmp.second ) {
        // Shared with a higher ranking process
        tmp.first->m_owner_rank = tmp.second ; 
        tmp.first->m_sync_count = m_sync_count ;
        internal_change_entity_parts( *tmp.first , add_parts , remove_parts );
      }
    }
  }


  {  
    const size_t n_old = m_shares_all.size();
     
    m_shares_all.insert( m_shares_all.end() ,  
                         new_shared.begin() , new_shared.end() );
     
    std::inplace_merge( m_shares_all.begin() ,  
                        m_shares_all.begin() + n_old ,
                        m_shares_all.end() ,
                        EntityLess() );
  }
   
  {
    const std::vector<EntityProc>::iterator es = m_shares_all.end();
          std::vector<EntityProc>::iterator js = m_shares_all.begin();
 
    while ( js != es ) {
      const std::vector<EntityProc>::iterator is = js ;
      for ( ; js != es && is->first == js->first ; ++js );
      PairIterEntityProc & pit = is->first->m_sharing ;
      pit.first  = is ;
      pit.second = js ;
    }
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void BulkData::internal_resolve_parallel_create_delete()
{
  std::vector<EntityProc> send_entities ;
  std::vector<EntityProc> del_entities_remote ;
  std::vector<EntityProc> owner_change ;
 
  // 'm_del_entities' is guaranteed unique, but is not sorted.
  // Sort it to for call to 'send_shared_ghosted'.

  std::sort( m_del_entities.begin(), m_del_entities.end(), EntityLess() );

  send_owned_shared_ghosted( m_parallel_rank ,
                             m_ghosting , m_del_entities , send_entities );

  bool global_delete_flag =
    communicate_entity_vector( *this , send_entities , del_entities_remote );

  if ( global_delete_flag ) {
    // At least one process has deleted an entity that is
    // not owned, shared, or ghosted.

    send_entities.clear();

    //------------------------------
    // Clear the ghost communication lists of deleted send and receive entities

    internal_resolve_destroy_ghosted( m_del_entities , del_entities_remote );

    //------------------------------
    // Clear the shared communication list of deleted send and receive entities.
    // Generate list for parallel ownership change.

    internal_resolve_destroy_shared( m_del_entities , del_entities_remote );

    del_entities_remote.clear();
  }

  //------------------------------

  std::vector<EntityProc> shared_created ;

  // Update the parallel index and
  // output sharing of created entities.
  internal_update_parallel_index( shared_created );

  //------------------------------
  // Update shared created entities.
  // - Set owner to highest ranking process, removed owned part as necessary.
  // - Update sharing.
  internal_resolve_created_shared( shared_created );
}

bool BulkData::modification_end()
{
  if ( m_sync_state == SYNCHRONIZED ) { return false ; }

  int local_change_count[2] ;
  int global_change_count[2] ;
  local_change_count[0] = m_new_entities.size();
  local_change_count[1] = m_del_entities.size();

  all_reduce_sum( m_parallel_machine ,
                  local_change_count , global_change_count , 2 );

  if ( global_change_count[0] || global_change_count[1] ) {

    if ( 1 < m_parallel_size ) {
      internal_resolve_parallel_create_delete();
    }

    m_new_entities.clear();

    while ( ! m_del_entities.empty() ) {
      internal_destroy_entity( m_del_entities.back() );
      m_del_entities.pop_back();
    }
  }

  // Parallel distribution considerations:
  if ( 1 < m_parallel_size ) {

    // Regenerate the ghosting aura around all shared mesh entities.
    internal_regenerate_shared_aura();

    // Resolve changed part membership:
    internal_resolve_shared_membership();
  }
  // Parallel consistency is restored.
  // ------------------------------
  // The very last operation performed is to sort the bucket entities.
  // This does not change the entities, relations, or field data.
  // However, it insures that the ordering of entities and buckets
  // is independent of the order in which a set of changes were
  // performed.

  internal_sort_bucket_entities();

  // ------------------------------

  ++m_sync_count ;

  m_sync_state = SYNCHRONIZED ;

  return true ;
}


//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

