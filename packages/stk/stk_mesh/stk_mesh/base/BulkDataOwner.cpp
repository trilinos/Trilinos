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
// Set the shared entities' pair of sharing iterators.
// The m_sharing is guaranteed to be consistent with
// the locally owned and locally used part membership.

void BulkData::internal_set_shared_entities()
{
  if (parallel_size() == 1) return;

  // Clear the entities' sharing in preparation for setting it.
 
  {
    const EntitySet::iterator e = m_entities.end();
          EntitySet::iterator i = m_entities.begin();
    for ( ; e != i ; ++i ) {
      PairIterEntityProc & pit = i->second->m_sharing ;
      pit.second = pit.first = std::vector< EntityProc >::const_iterator();
    }
  }
 
  // Set the entities' sharing.
 
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

typedef std::pair< EntityKey , unsigned > KeyProc ;

namespace {

struct KeyProcLess {
  bool operator()( const KeyProc & lhs , const KeyProc & rhs ) const
    { return lhs.first != rhs.first ? lhs.first < rhs.first
                                    : lhs.second < rhs.second ; }
};

}

//----------------------------------------------------------------------

namespace {

void insert_closure_ghost( const MetaData & meta ,
                           Entity * const entity ,
                           std::set<Entity*,EntityLess> & remove_list )
{
  if ( ! entity->bucket().member( meta.locally_used_part() ) ) {
    // This entity is a ghost, put it on the remove_list
    // along with all ghosts in its closure

    std::pair< std::set<Entity*,EntityLess>::iterator , bool >
      result = remove_list.insert( entity );

    if ( result.second ) {
      // This ghost entity is new to the list, traverse its closure.

      const unsigned etype = entity->entity_type();

      for ( PairIterRelation
            irel = entity->relations() ; ! irel.empty() ; ++irel ) {

        if ( irel->entity_type() < etype ) {
          insert_closure_ghost( meta , irel->entity() , remove_list );
        }
      }
    }
  }
}

void insert_transitive_ghost( const MetaData & meta ,
                              Entity * const entity ,
                              std::set<Entity*,EntityLess> & remove_list )
{
  insert_closure_ghost( meta , entity , remove_list );

  // Transitive:
  // If this entity is a member of another entity's closure
  // then that other entity is part of the traversal.

  const unsigned etype = entity->entity_type();

  for ( PairIterRelation rel = entity->relations(); ! rel.empty() ; ++rel ) {
    if ( etype < rel->entity_type() ) {
      insert_transitive_ghost( meta , rel->entity() , remove_list );
    }
  }
}

//----------------------------------------------------------------------

void insert_closure_send(
  const EntityProc                  send_entry ,
  std::set<EntityProc,EntityLess> & send_list )
{
  std::pair< std::set<EntityProc,EntityLess>::iterator , bool >
    result = send_list.insert( send_entry );

  if ( result.second ) {
    // First time this entity was inserted into the send_list.

    const unsigned etype  = send_entry.first->entity_type();
    PairIterRelation irel = send_entry.first->relations();

    for ( ; ! irel.empty() ; ++irel ) {
      if ( irel->entity_type() < etype ) {
        const EntityProc rel_send_entry( irel->entity(), send_entry.second );

        insert_closure_send( rel_send_entry , send_list );
      }
    }
  }
}

//----------------------------------------------------------------------

bool member_of_owned_closure( const Entity & e , const unsigned p_rank )
{
  bool result = p_rank == e.owner_rank();

  const unsigned etype = e.entity_type();

  // Any higher ranking entities locally owned?
  for ( PairIterRelation
        irel = e.relations(); ! result && ! irel.empty() ; ++irel ) {
    result = etype  <  irel->entity_type() &&
             p_rank == irel->entity()->owner_rank();
  }

  // Any higher ranking entity member of an owned closure?
  for ( PairIterRelation
        irel = e.relations(); ! result && ! irel.empty() ; ++irel ) {
    result = etype < irel->entity_type() &&
             member_of_owned_closure( * irel->entity() , p_rank );
  }

  return result ;
}

//----------------------------------------------------------------------

void clean_and_verify_parallel_change(
  const char method[] ,
  const BulkData & mesh ,
  std::vector<EntityProc> & local_change )
{
  const MetaData      & meta   = mesh.mesh_meta_data() ;
  const unsigned        p_rank = mesh.parallel_rank();
  const unsigned        p_size = mesh.parallel_size();
  const ParallelMachine p_comm = mesh.parallel();

  size_t error_count = 0 ;

  std::ostringstream error_msg ;

  // Order and eliminate redundancies:
  {
    std::vector<EntityProc>::iterator i = local_change.begin() ,
                                      j = local_change.end() ;
    std::sort( i , j , EntityLess() );
    i = std::unique( i , j );
    local_change.erase( i , j );
  }

  for ( std::vector<EntityProc>::iterator 
        i = local_change.begin() ; i != local_change.end() ; ) {
    std::vector<EntityProc>::iterator j = i ; ++i ;
    Entity * const entity    = j->first ;
    const unsigned new_owner = j->second ;
  
    // Verification:
    // 1) If bucket has no capacity then is destined for deletion
    // 2) If not locally owned then not allowed grant ownership
    // 3) New owner must be legit
    // 4) Cannot grant to two different owners

    const bool bad_delete = 0 == entity->bucket().capacity();
    const bool bad_entity = entity->owner_rank() != p_rank ;
    const bool bad_owner  = p_size <= new_owner ;
    const bool bad_dup    = i != local_change.end() && entity == i->first ;

    if ( bad_entity || bad_owner || bad_dup || bad_delete ) {
      ++error_count ;

      error_msg << "  P" << p_rank << ": " ;
      print_entity_key( error_msg , meta , entity->key() );
      if ( bad_delete ) { error_msg << " HAS_BEEN_DELETED" ; }
      if ( bad_entity ) { error_msg << " NOT_CURRENT_OWNER" ; }
      if ( bad_owner ) {
        error_msg << " BAD_NEW_OWNER( " << new_owner << " )" ;
       }
      if ( bad_dup ) {
        error_msg << " CONFLICTING_NEW_OWNER( " << new_owner ;
        error_msg << " != " << i->second << " )" ;
      }
      error_msg << std::endl ;
    }
    else if ( new_owner == p_rank ) {
      // Eliminate non-changes
      j->first = NULL ;
      j->second = 0 ;
    }
  }

  all_reduce( p_comm , ReduceSum<1>( & error_count ) );

  if ( error_count ) {
    std::string msg_throw ;
    msg_throw.append( method );
    msg_throw.append( " FAILED: Bad change ownership directives" );

    if ( 0 == p_rank ) { std::cerr << msg_throw ; }

    all_write_string( p_comm , std::cerr , error_msg.str() );

    throw std::runtime_error( msg_throw );
  }

  {
    std::vector<EntityProc>::iterator i = local_change.begin(),
                                      j = local_change.end();
    i = std::remove( i , j , EntityProc(NULL,0) );
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
                               const std::vector<EntityProc> & send_change ,
                                     std::vector<EntityProc> & shared_change ,
                                     std::vector<EntityProc> & ghosted_change ,
                                     std::vector<KeyProc>    & index_change )
{
  const MetaData & meta = mesh.mesh_meta_data() ;
  const unsigned p_rank = mesh.parallel_rank();
  const unsigned p_size = mesh.parallel_size();

  // Local changes to parallel index:
  for ( std::vector<EntityProc>::const_iterator
        ip = local_change.begin() ; ip != local_change.end() ; ++ip ) {
    if ( p_rank == mesh.parallel_index_rank( ip->first->key() ) ) {
      const KeyProc entry( ip->first->key() , ip->second );
      index_change.push_back( entry );
    }
  }

  CommAll comm( mesh.parallel() );

  // Sizing:

  for ( std::vector<EntityProc>::const_iterator
        ip = send_change.begin() ; ip != send_change.end() ; ++ip ) {
    comm.send_buffer( ip->second ).skip<EntityKey>( 1 ).skip<unsigned>( 1 );
  }

  // Allocation:

  comm.allocate_buffers( p_size / 4 , 0 );

  // Packing new owner info:
  {
    std::vector<EntityProc>::const_iterator i = local_change.begin();
    std::vector<EntityProc>::const_iterator ip = send_change.begin() ;
    for ( ; i != local_change.end() ; ++i ) {
      for ( ; ip != send_change.end() && ip->first == i->first ; ++ip ) {
        comm.send_buffer( ip->second ).pack<EntityKey>( i->first->key() )
                                      .pack<unsigned>( i->second );
      }
    }
  }

  comm.communicate();

  for ( unsigned ip = 0 ; ip < p_size ; ++ip ) {
    CommBuffer & buf = comm.recv_buffer( ip );
    while ( buf.remaining() ) {
      EntityProc entry ;
      EntityKey key ;
      buf.unpack<EntityKey>( key ).unpack<unsigned>( entry.second );

      // Remote changes to parallel index:
      if ( p_rank == mesh.parallel_index_rank( key ) ) {
        index_change.push_back( KeyProc( key , entry.second ) );
      }

      entry.first = mesh.get_entity( entity_type(key), entity_id(key) );

      if ( NULL != entry.first ) {
        if ( entry.first->bucket().member( meta.locally_used_part() ) ) {
          shared_change.push_back( entry );
        }
        else {
          ghosted_change.push_back( entry );
        }
      }
    }
  }

  std::sort( shared_change.begin() , shared_change.end() , EntityLess() );
  std::sort( ghosted_change.begin() , ghosted_change.end() , EntityLess() );
  std::sort( index_change.begin() , index_change.end() );
}

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

unsigned BulkData::parallel_index_rank( const EntityKey & key ) const
{
  return ( entity_id( key ) >> 12 ) % m_parallel_size ;
}

void BulkData::owner_send_to_all( Entity * const entity ,
                                  std::vector<EntityProc> & send_change ) const
{
  const unsigned p_rank = m_parallel_rank ;

  const size_t n_old = send_change.size();

  // Shared processes:
  {
    PairIterEntityProc share = entity->sharing() ;
    send_change.insert( send_change.end() , share.begin() , share.end() );
  }

  // Ghosting processes:
  for ( std::vector<Ghosting*>::const_iterator
        j =  m_ghosting.begin() ;
        j != m_ghosting.end() ; ++j ) {

    const std::vector<EntityProc> & ghost_send_list = (*j)->send();

    std::vector<EntityProc>::const_iterator k_beg = ghost_send_list.begin(),
                                            k_end = ghost_send_list.end();
    k_beg = std::lower_bound( k_beg , k_end ,
                              EntityProc( entity , 0 ),
                              EntityLess() );

    for ( k_end = k_beg ;
          k_end != ghost_send_list.end() &&
          k_end->first == entity ; ++k_end );

    send_change.insert( send_change.end() , k_beg , k_end );
  }

  // Parallel index process:
  {
    EntityProc entry ;
    entry.first  = entity ;
    entry.second = parallel_index_rank( entity->key() );
    if ( p_rank != entry.second ) { send_change.push_back( entry ); }
  }

  // Sort to proper order:
  {
    std::vector<EntityProc>::iterator ii = send_change.begin() + n_old ,
                                      jj = send_change.end();
    std::sort( ii , jj , EntityLess() );
    ii = std::unique( ii , jj );
    send_change.erase( ii , jj );
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void BulkData::change_entity_owner( const std::vector<EntityProc> & arg_change )
{
  static const char method[] = "stk::mesh::BulkData::change_entity_owner" ;

  const MetaData  & meta = m_mesh_meta_data ;
  const unsigned  p_rank = m_parallel_rank ;
  const unsigned  p_size = m_parallel_size ;
  ParallelMachine p_comm = m_parallel_machine ;

  //------------------------------
  // Verify the input changes, generate a clean local change list, and
  // generate the remote change list so that all processes know about
  // pending changes.

  std::vector<EntityProc> local_change( arg_change );

  // Parallel synchronous clean up and verify the requested changes:
  clean_and_verify_parallel_change( method , *this , local_change );

  std::vector<EntityProc> send_change ;
  std::vector<EntityProc> ghosted_change ;
  std::vector<EntityProc> shared_change ;
  std::vector<KeyProc>    index_change ;

  //----------------------------------------
  // Generate complete send list from
  // owner to shared, ghosted, and parallel index.

  // Local change is properly ordered by entity key:

  for ( std::vector<EntityProc>::iterator
        i = local_change.begin() ; i != local_change.end() ; ++i ) {

    // Append shared, ghosted, and indexing processes:
    owner_send_to_all( i->first , send_change );
  }

  // Parallel synchronous determination of changing
  // shared, ghosted, and index:

  generate_parallel_change( *this , local_change , send_change ,
                            shared_change , ghosted_change , index_change );

  //------------------------------
  // Have enough information to delete all effected ghosts.
  // If the closure of a ghost contains a changing entity
  // then that ghost must be deleted.
  // Request that all ghost entities in the closure of the ghost be deleted.

  {
    std::set<Entity*,EntityLess> work ;

    for ( std::vector<EntityProc>::const_iterator
          i = ghosted_change.begin() ; i != ghosted_change.end() ; ++i ) {
      insert_transitive_ghost( meta , i->first , work );
    }

    for ( std::vector<EntityProc>::const_iterator
          i = shared_change.begin() ; i != shared_change.end() ; ++i ) {
      insert_transitive_ghost( meta , i->first , work );
    }

    for ( std::vector<EntityProc>::const_iterator
          i = local_change.begin() ; i != local_change.end() ; ++i ) {
      insert_transitive_ghost( meta , i->first , work );
    }

    // The ghosted change list will become invalid
    ghosted_change.clear();

    std::vector<EntityProc> empty ;
    std::vector<Entity*> effected_ghosts( work.begin() , work.end() );

    for ( std::vector<Ghosting*>::iterator
          ig = m_ghosting.begin() ; ig != m_ghosting.end() ; ++ig ) {
      // parallel synchronous:
      internal_change_ghosting( **ig , empty , effected_ghosts );
    }
  }

  //------------------------------
  // Consistently change the owner on all processes.
  // 1) The local_change list is giving away ownership.
  // 2) The shared_change may or may not be receiving ownership
  // 3) The index_change updates m_entities_owner_index

  {
    PartVector owned( 1 );
    owned[0] = & meta.locally_owned_part();

    for ( std::vector<EntityProc>::iterator
          i = local_change.begin() ; i != local_change.end() ; ++i ) {
      // Giving ownership, change the parts first and then
      // the owner rank to pass the ownership test.
      change_entity_parts( * i->first , PartVector() , owned );
      i->first->m_owner_rank = i->second ;
    }

    for ( std::vector<EntityProc>::iterator
          i = shared_change.begin() ; i != shared_change.end() ; ++i ) {
      i->first->m_owner_rank = i->second ;
      if ( p_rank == i->second ) { // I receive ownership
        change_entity_parts( * i->first , owned , PartVector() );
      }
    }
  }

  {
    for ( std::vector<KeyProc>::iterator
          i = index_change.begin() ; i != index_change.end() ; ++i ) {
      KeyProc tmp( i->first , 0 );
      std::vector<KeyProc>::iterator j =
        std::lower_bound( m_entities_owner_index.begin() ,
                          m_entities_owner_index.end() ,
                          tmp , KeyProcLess() );
      if ( j == m_entities_owner_index.end() || j->first != i->first ) {
        // Theoretically cannot happen
        throw std::logic_error( std::string( method ) );
      }
      j->second = i->second ;
    }
  }

  //------------------------------
  // Send entities, along with their closure, to the new owner processes
  {
    std::set< EntityProc , EntityLess > send_closure ;
    std::ostringstream error_msg ;
    int error_count = 0 ;

    for ( std::vector<EntityProc>::iterator
          i = local_change.begin() ; i != local_change.end() ; ++i ) {
      insert_closure_send( *i , send_closure );
    }

    CommAll comm( p_comm );

    for ( std::set<EntityProc,EntityLess>::iterator
          i = send_closure.begin() ; i != send_closure.end() ; ++i ) {
      CommBuffer & buffer = comm.send_buffer( i->second );
      Entity     & entity = * i->first ;
      pack_entity_info( buffer , entity );
      pack_field_values( buffer , entity );
    }

    comm.allocate_buffers( p_size / 4 );

    for ( std::set<EntityProc,EntityLess>::iterator
          i = send_closure.begin() ; i != send_closure.end() ; ++i ) {
      CommBuffer & buffer = comm.send_buffer( i->second );
      Entity     & entity = * i->first ;
      pack_entity_info( buffer , entity );
      pack_field_values( buffer , entity );
    }

    // Any entity that I sent and is not in an owned closure is deleted.
    // Destroy backwards so as not to invalidate closures in the process.

    for ( std::set<EntityProc,EntityLess>::iterator
          i = send_closure.end() ; i != send_closure.begin() ; ) {
      --i ;
      if ( ! member_of_owned_closure( * i->first , p_rank ) ) {
        internal_destroy_entity( i->first );
      }
    }

    send_closure.clear(); // Has been invalidated

    comm.communicate();

    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      CommBuffer & buf = comm.recv_buffer(p);
      while ( buf.remaining() ) {
        PartVector parts ;
        std::vector<Relation> relations ;
        EntityKey key ;
        unsigned  owner = ~0u ;

        unpack_entity_info( buf, *this, key, owner, parts, relations );

        // Received entity information will be correct, modulo the owned part
        if ( owner == p_rank ) {
          // Must have the locally_owned_part
          insert( parts , meta.locally_owned_part() );
        }
        else {
          // Must not have the locally_owned_part
          remove( parts , meta.locally_owned_part() );
        }

        std::pair<Entity*,bool> result = internal_create_entity( key , owner );
        internal_change_entity_parts( *result.first , parts , PartVector() );

        declare_relation( *result.first , relations );

        if ( ! unpack_field_values( buf , * result.first , error_msg ) ) {
          ++error_count ;
        }
      }
    }

    all_reduce( p_comm , ReduceSum<1>( & error_count ) );

    if ( error_count ) { throw std::runtime_error( error_msg.str() ); }
  }

  //------------------------------
  // Regnerate sharing information for existing entities.
  // Send to owner and owner sends to sharing processes
  // Sharing changes go beyond ownership change due to including
  // the closure of the changed owners.

  internal_regenerate_sharing();
}

//----------------------------------------------------------------------

namespace {

void pack_owned( const std::set<EntityProc,EntityLess> & work ,
                 const unsigned p_rank ,
                 CommAll & comm )
{
  for ( std::set<EntityProc,EntityLess>::const_iterator
        j = work.begin() ; j != work.end() ; ) {

    const std::set<EntityProc,EntityLess>::const_iterator i = j ;

    for ( ; j != work.end() && j->first == i->first ; ++j );

    if ( i->first->owner_rank() == p_rank ) {

      for ( std::set<EntityProc,EntityLess>::const_iterator
            k = i ; k != j ; ++k ) {
        CommBuffer & buf = comm.send_buffer( k->second );

        for ( std::set<EntityProc,EntityLess>::const_iterator
              kk = i ; kk != j ; ++kk ) {
          buf.pack<EntityKey>( kk->first->key() );
          buf.pack<unsigned>(  kk->second );
        }
      }
    }
  }
}

}

//----------------------------------------------------------------------

void BulkData::internal_regenerate_sharing()
{
  static const char method[] =
    "stk::mesh::BulkData::internal_regenerate_sharing" ;

  std::set<EntityProc,EntityLess> work ;

  bool local_error_flag = false ;
  bool global_error_flag = false ;

  {
    CommAll comm_to_owner( parallel() );

    for ( std::vector< std::vector<Bucket*> >::iterator
          ib = m_buckets.begin() ; ib != m_buckets.end() ; ++ib ) {
      for ( std::vector<Bucket*>::iterator
            jb = ib->begin() ; jb != ib->end() ; ++jb ) {

        Bucket & bucket = **jb ;

        const bool locally_used =
          bucket.member( m_mesh_meta_data.locally_used_part() );

        const bool locally_owned =
          bucket.member( m_mesh_meta_data.locally_owned_part() );

        if ( ! locally_owned ) { // Does not own these:

          for ( unsigned k = 0 ; k < bucket.size() ; ++k ) {
            if ( bucket[k].owner_rank() == m_parallel_rank ) {
              local_error_flag = true ;
            }
          }

          if ( locally_used ) { // These are shared:
            for ( unsigned k = 0 ; k < bucket.size() ; ++k ) {
              EntityProc entry ;
              entry.first = & bucket[k] ;
              entry.second = entry.first->owner_rank();
              work.insert( entry );
              comm_to_owner.send_buffer( entry.second ).skip<EntityKey>(1);
            }
          }
          else { // These are ghosted:
            ;
          }
        }
        else { // Owns these
          for ( unsigned k = 0 ; k < bucket.size() ; ++k ) {
            if ( bucket[k].owner_rank() != m_parallel_rank ) {
              local_error_flag = true ;
            }
          }
        }
      }
    }

    global_error_flag =
      comm_to_owner.allocate_buffers( parallel_size() / 4 ,
                                      false /* not symmetric */ ,
                                      local_error_flag );

    if ( global_error_flag ) {
      // Should not happen:
      throw std::logic_error( std::string( method ) );
    }
  
    for ( std::set<EntityProc,EntityLess>::iterator
          i = work.begin() ; i != work.end() ; ++i ) {
      comm_to_owner.send_buffer( i->second ).pack<EntityKey>( i->first->key() );
    }

    comm_to_owner.communicate();

    for ( unsigned p = 0 ; p < m_parallel_size ; ++p ) {
      CommBuffer & buf = comm_to_owner.recv_buffer( p );
      while ( buf.remaining() ) {
        EntityKey key ; buf.unpack<EntityKey>( key );
        EntityProc entry ;
        entry.first = get_entity( entity_type(key), entity_id(key) );
        entry.second = p ;

        // This process must own every entity received:
        if ( entry.first->owner_rank() != m_parallel_rank ) {
          local_error_flag = true ;
        }
        work.insert( entry );
      }
    }
  }

  // The owner of an entity now knows which other processes share that entity

  {
    CommAll comm_from_owner( parallel() );

    pack_owned( work , m_parallel_rank , comm_from_owner );

    global_error_flag =
      comm_from_owner.allocate_buffers( m_parallel_size / 4 ,
                                        false /* not symmetric */ ,
                                        local_error_flag );

    if ( global_error_flag ) {
      // Should not happen:
      throw std::logic_error( std::string( method ) );
    }

    pack_owned( work , m_parallel_rank , comm_from_owner );

    comm_from_owner.communicate();

    for ( unsigned p = 0 ; p < m_parallel_size ; ++p ) {
      CommBuffer & buf = comm_from_owner.recv_buffer( p );
      while ( buf.remaining() ) {
        EntityKey key ;    buf.unpack<EntityKey>( key );
        EntityProc entry ; buf.unpack<unsigned>( entry.second );
        // Don't insert myself
        if ( entry.second != m_parallel_rank ) {
          entry.first = get_entity( entity_type(key), entity_id(key) );
          work.insert( entry );
        }
      }
    }
  }

  m_shares_all.assign( work.begin() , work.end() );

  internal_set_shared_entities();
}

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

