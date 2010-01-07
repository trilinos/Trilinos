
/**
 * @author H. Carter Edwards
 */

#include <set>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <assert.h>
#include <string.h>

#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/Ghosting.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/EntityComm.hpp>
#include <stk_mesh/base/Comm.hpp>

namespace stk {
namespace mesh {

//----------------------------------------------------------------------

Ghosting & BulkData::create_ghosting( const std::string & name )
{
  static const char method[] = "stk::mesh::BulkData::create_ghosting" ;

  assert_ok_to_modify( method );

  // Verify name is the same on all processors,
  // if not then throw an exception on all processors.
  {
    CommBroadcast bc( parallel() , 0 );

    if ( bc.parallel_rank() == 0 ) {
      bc.send_buffer().skip<char>( name.size() + 1 );
    }

    bc.allocate_buffer();

    if ( bc.parallel_rank() == 0 ) {
      bc.send_buffer().pack<char>( name.c_str() , name.size() + 1 );
    }

    bc.communicate();

    const char * const bc_name = reinterpret_cast<const char *>( bc.recv_buffer().buffer() );

    int error = 0 != strcmp( bc_name , name.c_str() );

    all_reduce( parallel() , ReduceMax<1>( & error ) );

    if ( error ) {
      std::string msg("stk::mesh::BulkData::create_ghosting ERROR: Parallel name inconsistency");
      throw std::runtime_error( msg );
    }
  }

  Ghosting * const g = new Ghosting( *this , name , m_sync_count );
  m_ghosting.push_back( g );
  return *g ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void BulkData::destroy_all_ghosting()
{
  static const char method[] = "stk::mesh::BulkData::destroy_all_ghosting" ;

  assert_ok_to_modify( method );

  Part & uses_part = m_mesh_meta_data.locally_used_part();

  // Clear Ghosting data

  for ( std::vector<Ghosting*>::iterator
        ig = m_ghosting.begin() ; ig != m_ghosting.end() ; ++ig ) {
    Ghosting & gh = **ig ;
    gh.m_send.clear();
    gh.m_recv.clear();
    gh.m_sync_count = m_sync_count ;
  }

  // Destroy all entities that are not in the uses_part.
  // Iterate backwards so as not to invalidate a closure.

  for ( std::vector< std::vector<Bucket*> >::iterator
        iib = m_buckets.end() ; iib != m_buckets.begin() ; ) {
    std::vector<Bucket*> & buckets = *--iib ;

    for ( size_t ib = buckets.size() ; ib ; ) {
      Bucket * const b = buckets[--ib] ;

      if ( ! has_superset( *b , uses_part ) ) {
        for ( unsigned ie = b->size() ; ie ; ) {
          internal_destroy_entity( & (*b)[--ie] );
        }
        // The bucket reference is now invalid.
      }
    }
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

bool verify_change_ghosting_arguments(
  BulkData & mesh ,
  Ghosting & ghosts ,
  const std::vector<EntityProc> & add_send ,
  const std::vector<Entity*> & remove_receive ,
  std::ostream & error_log );

void insert_transitive_closure( std::set<EntityProc,EntityLess> & new_send ,
                                const EntityProc & entry );

void comm_recv_to_send(
  BulkData & mesh ,
  const std::set< Entity * , EntityLess > & new_recv ,
        std::set< EntityProc , EntityLess > & new_send );

void comm_sync_send_recv(
  BulkData & mesh ,
  std::set< EntityProc , EntityLess > & new_send ,
  std::set< Entity * , EntityLess > & new_recv );

inline
bool is_receive_member( const Ghosting & gh , Entity * const entity )
{
  const std::vector<Entity*>::const_iterator j = gh.receive().end(); 
        std::vector<Entity*>::const_iterator i = gh.receive().begin(); 

  i = std::lower_bound( i , j , entity , EntityLess() );

  return j != i && entity == *i ;
}

} // namespace <>

//----------------------------------------------------------------------

void BulkData::change_ghosting(
  Ghosting & ghosts ,
  const std::vector<EntityProc> & add_send ,
  const std::vector<Entity*> & remove_receive )
{
  static const char method[] = "stk::mesh::BulkData::change_ghosting" ;

  assert_ok_to_modify( method );

  bool status = true ;
  
  std::ostringstream error_log ;

  if ( & ghosts.mesh() != this ) {
    status = false ;
    error_log << method
              << " ERROR: Mesh does not own Ghosting(" 
              << ghosts.name() << ")" << std::endl ;
  }
  else if ( & ghosts == m_ghosting[0] ) {
    status = false ;
    error_log << method
              << " ERROR: Cannot change shared_aura" << std::endl ;
  }
  else {
    status = verify_change_ghosting_arguments(
               *this, ghosts, add_send, remove_receive, error_log );
  }

  if ( ! status ) {
    throw std::runtime_error( error_log.str() );
  }

  internal_change_ghosting( ghosts , add_send , remove_receive );
}

void BulkData::internal_change_ghosting(
  Ghosting & ghosts ,
  const std::vector<EntityProc> & add_send ,
  const std::vector<Entity*> & remove_receive )
{
  const MetaData & meta = m_mesh_meta_data ;
  const unsigned p_size = m_parallel_size ;

  //------------------------------------
  // Copy ghosting lists into more efficiently editted container.

  std::set< EntityProc , EntityLess > new_send ;
  std::set< Entity * ,   EntityLess > new_recv ;

  //------------------------------------
  // Insert the current ghost receives and then remove from that list.

  if ( & ghosts.receive() != & remove_receive ) {

    new_recv.insert( ghosts.receive().begin() , ghosts.receive().end() );

    for ( std::vector< Entity * >::const_iterator
          i = remove_receive.begin() ; i != remove_receive.end() ; ++i ) {
      new_recv.erase( *i );
    }

    //  Keep the closure of the remaining received ghosts.
    //  Working from highest-to-lowest key (rank entity type)
    //  results in insertion of the transitive closure.
    //  Insertion will not invalidate the associative container's iterator.

    for ( std::set< Entity * , EntityLess >::iterator
          i = new_recv.end() ; i != new_recv.begin() ; ) {
      --i ;

      const unsigned etype = (*i)->entity_type();

      for ( PairIterRelation
            irel = (*i)->relations(); ! irel.empty() ; ++irel ) {
        if ( irel->entity_type() < etype &&
             is_receive_member( ghosts , irel->entity() ) ) {
          new_recv.insert( irel->entity() );
        }
      }
    }
  }

  //  Initialize the new_send from the new_recv
  comm_recv_to_send( *this , new_recv , new_send );

  //------------------------------------
  // Add the specified entities and their closure to the send ghosting
  //
  for ( std::vector< EntityProc >::const_iterator
        i = add_send.begin() ; i != add_send.end() ; ++i ) {
    insert_transitive_closure( new_send , *i );
  }

  comm_sync_send_recv( *this , new_send , new_recv );

  // The new_send list is now complete and accurate
  // The new_recv has those ghost entities that are to be kept.
  //------------------------------------
  // Remove the ghost entities defined by: ghosts.m_recv - new_recv
  // If the last reference to the ghost entity then delete it.

  while ( ! ghosts.m_recv.empty() ) {
    Entity * const entity = ghosts.m_recv.back();
    ghosts.m_recv.pop_back();

    bool destroy_it = ! new_recv.count( entity );

    for ( std::vector<Ghosting*>::const_iterator
          i = m_ghosting.begin() ; destroy_it && i != m_ghosting.end() ; ++i ) {
      destroy_it = ! is_receive_member( **i , entity );
    }

    if ( destroy_it ) {
      internal_destroy_entity( entity );
    }
  }

  //------------------------------------
  // Push newly ghosted entities to the receivers and update the new_recv list

  {
    const EntityLess less ;

    std::vector< EntityProc > send ;
    std::vector< EntityProc > recv ;

    for ( std::set< EntityProc , EntityLess >::iterator
          j = new_send.begin(); j != new_send.end() ; ++j ) {

      std::vector< EntityProc >::const_iterator i = ghosts.send().begin();

      for ( ; i != ghosts.send().end() && less(*i,*j) ; ++i );

      if ( i == ghosts.send().end() || *i != *j ) {
        send.push_back( *j );
      }
    }

    CommAll comm( m_parallel_machine );

    for ( std::vector< EntityProc >::iterator
          k = send.begin() ; k != send.end() ; ++k ) {
      pack_entity_info( comm.send_buffer( k->second ) , * k->first );
    }

    comm.allocate_buffers( p_size / 4 );

    for ( std::vector< EntityProc >::iterator
          k = send.begin() ; k != send.end() ; ++k ) {
      pack_entity_info( comm.send_buffer( k->second ) , * k->first );
    }

    comm.communicate();

    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      CommBuffer & buf = comm.recv_buffer(p);
      while ( buf.remaining() ) {
        PartVector parts ;
        std::vector<Relation> relations ;
        EntityKey key ;
        unsigned  owner = ~0u ;
 
        unpack_entity_info( buf, *this, key, owner, parts, relations );
 
        // Must not have the locally_owned_part or locally_used_part

        remove( parts , meta.locally_owned_part() );
        remove( parts , meta.locally_used_part() );
 
        std::pair<Entity*,bool> result = internal_create_entity( key , owner );

        internal_change_entity_parts( * result.first , parts , PartVector() );
 
        declare_relation( * result.first , relations );

        new_recv.insert( result.first );
      }
    }
  }

  ghosts.m_send.assign( new_send.begin() , new_send.end() );
  ghosts.m_recv.assign( new_recv.begin() , new_recv.end() );
  ghosts.m_sync_count = m_sync_count ;
}

//----------------------------------------------------------------------

namespace {

//  Locally verify the input arguments,
//  globally synchronize the results, and
//  if an error is detected globally throw
//  the same exception - with the local text.

bool verify_change_ghosting_arguments(
  BulkData & mesh ,
  Ghosting & ghosts ,
  const std::vector<EntityProc> & add_send ,
  const std::vector<Entity*> & remove_receive ,
  std::ostream & error_log )
{
  int status = 0 ;

  // Verify all 'add' are locally owned.
  for ( std::vector<EntityProc>::const_iterator
        i = add_send.begin() ; ! status && i != add_send.end() ; ++i ) {
    if ( i->first->owner_rank() != mesh.parallel_rank() ) {
      print_entity_key( error_log , mesh.mesh_meta_data() , i->first->key() );
      error_log << " ERROR: add_send is not locally owned" ;
      status = -1 ;
    }
  }

  // Verify all 'remove' are members of the ghosting.
  if ( & ghosts.receive() != & remove_receive ) {
    for ( std::vector<Entity*>::const_iterator
          i = remove_receive.begin() ;
          ! status && i != remove_receive.end() ; ++i ) {

      if ( ! is_receive_member( ghosts , *i ) ) {
        print_entity_key( error_log , mesh.mesh_meta_data() , (*i)->key() );
        error_log << " ERROR: remove_receive is not in the ghosting" ;
        status = -1 ;
      }
    }
  }

  all_reduce( mesh.parallel() , ReduceMin<1>( & status ) );

  return 0 == status ;
}

void insert_transitive_closure( std::set<EntityProc,EntityLess> & new_send ,
                                const EntityProc & entry )
{
  // Do not insert if I can determine that this entity is already
  // owned or shared by the receiving processor.

  PairIterEntityProc sharing = entry.first->sharing();

  for ( ; ! sharing.empty() && entry != *sharing ; ++sharing );
 
  if ( entry.second != entry.first->owner_rank() && sharing.empty() ) {

    std::pair< std::set<EntityProc,EntityLess>::iterator , bool > 
      result = new_send.insert( entry );

    if ( result.second ) {
      // A new insertion, must also insert the closure

      const unsigned etype = entry.first->entity_type();
      PairIterRelation irel  = entry.first->relations();

      for ( ; ! irel.empty() ; ++irel ) {
        if ( irel->entity_type() < etype ) {
          EntityProc tmp( irel->entity() , entry.second );
          insert_transitive_closure( new_send , tmp );
        }
      }
    }
  }
}

// Fill a new send list from the receive list.

void comm_recv_to_send(
  BulkData & mesh ,
  const std::set< Entity * , EntityLess > & new_recv ,
        std::set< EntityProc , EntityLess > & new_send )
{
  static const char method[] = "stk::mesh::BulkData::change_ghosting" ;

  const unsigned parallel_size = mesh.parallel_size();

  CommAll all( mesh.parallel() );

  for ( std::set< Entity * , EntityLess >::const_iterator
        i = new_recv.begin() ; i != new_recv.end() ; ++i ) {
    const unsigned owner = (*i)->owner_rank();
    all.send_buffer( owner ).skip<EntityKey>(1);
  }

  all.allocate_buffers( parallel_size / 4 , false /* Not symmetric */ );

  for ( std::set< Entity * , EntityLess >::const_iterator
        i = new_recv.begin() ; i != new_recv.end() ; ++i ) {
    const unsigned owner = (*i)->owner_rank();
    const EntityKey key = (*i)->key();
    all.send_buffer( owner ).pack<EntityKey>( & key , 1 );
  }

  all.communicate();

  for ( unsigned p = 0 ; p < parallel_size ; ++p ) {
    CommBuffer & buf = all.recv_buffer(p);
    while ( buf.remaining() ) {
      EntityKey key ;
      buf.unpack<EntityKey>( & key , 1 );
      EntityProc tmp( mesh.get_entity( entity_type(key), entity_id(key) , method ) , p );
      new_send.insert( tmp );
    }
  }
}

// Synchronize the send list to the receive list.

void comm_sync_send_recv(
  BulkData & mesh ,
  std::set< EntityProc , EntityLess > & new_send ,
  std::set< Entity * , EntityLess > & new_recv )
{
  static const char method[] = "stk::mesh::BulkData::change_ghosting" ;
  const unsigned parallel_rank = mesh.parallel_rank();
  const unsigned parallel_size = mesh.parallel_size();

  CommAll all( mesh.parallel() );

  // Communication sizing:

  for ( std::set< EntityProc , EntityLess >::iterator
        i = new_send.begin() ; i != new_send.end() ; ++i ) {
    const unsigned owner = i->first->owner_rank();
    all.send_buffer( i->second ).skip<EntityKey>(2);
    if ( owner != parallel_rank ) {
      all.send_buffer( owner ).skip<EntityKey>(2);
    }
  }

  all.allocate_buffers( parallel_size / 4 , false /* Not symmetric */ );

  // Communication packing (with message content comments):
  for ( std::set< EntityProc , EntityLess >::iterator
        i = new_send.begin() ; i != new_send.end() ; ) {
    const unsigned owner = i->first->owner_rank();

    // Inform receiver of ghosting, the receiver does not own
    // and does not share this entity.
    // The ghost either already exists or is a to-be-done new ghost.
    // This status will be resolved on the final communication pass
    // when new ghosts are packed and sent.

    const EntityKey &entity_key = i->first->key();
    const uint64_t &proc = i->second;
    
    all.send_buffer( i->second ).pack(entity_key).pack(proc);

    if ( owner != parallel_rank ) {
      // I am not the owner of this entity.
      // Inform the owner of this ghosting need.
      all.send_buffer( owner ).pack(entity_key).pack(proc);

      // Erase it from my processor's ghosting responsibility:
      // The iterator passed to the erase method will be invalidated.
      std::set< EntityProc , EntityLess >::iterator jrem = i ; ++i ;
      new_send.erase( jrem );
    }
    else {
      ++i ;
    }
  }

  all.communicate();

  // Communication unpacking:
  for ( unsigned p = 0 ; p < parallel_size ; ++p ) {
    CommBuffer & buf = all.recv_buffer(p);
    while ( buf.remaining() ) {

      EntityKey entity_key;
      uint64_t proc(0);
    
      buf.unpack(entity_key).unpack(proc);

      if ( parallel_rank != proc ) {
        //  Receiving a ghosting need for an entity I own.
        //  Add it to my send list.
        EntityProc tmp( mesh.get_entity( entity_type(entity_key), entity_id(entity_key) , method ) , proc );
        new_send.insert( tmp );
      }
      else {
        //  I am the receiver for this ghost.
        //  If I already have it add it to the receive list,
        //  otherwise don't worry about it - I will receive
        //  it in the final new-ghosting communication.
        Entity * e = mesh.get_entity( entity_type(entity_key), entity_id(entity_key) );
        if ( e ) {
          new_recv.insert( e );
        }
      }
    }
  }
}

} // namespace <>

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void BulkData::internal_regenerate_shared_aura()
{
  static const char method[] =
    "stk::mesh::BulkData::internal_regenerate_shared_aura" ;

  assert_ok_to_modify( method );

  std::vector<EntityProc> send ;

  const std::vector<EntityProc>::const_iterator
    shares_beg = m_shares_all.begin(),
    shares_end = m_shares_all.end();

  for ( std::vector<EntityProc>::const_iterator
        i = shares_beg ; i != shares_end ; ++i ) {

    const unsigned etype = i->first->entity_type();

    for ( PairIterRelation
          rel = i->first->relations() ; ! rel.empty() ; ++rel ) {

      Entity * const entity = rel->entity();

      // Higher rank and I own it, ghost to the sharing processor
      if ( etype < entity->entity_type() &&
                   entity->owner_rank() == m_parallel_rank ) {

        // ... unless already shared by the destination processor

        PairIterEntityProc sh = entity->sharing();

        while ( ! sh.empty() && sh->second != i->second ) { ++sh ; }

        if ( sh.empty() ) {
          EntityProc entry( entity , i->second );
          send.push_back( entry );
        }
      }
    }
  }

  Ghosting & aura = * m_ghosting[0] ;

  // Add new aura, remove all of the old aura.
  // The change_ghosting figures out what to actually delete and add.

  internal_change_ghosting( aura , send , aura.receive() );
}

void BulkData::destroy_shared_aura()
{
  static const char method[] = "stk::mesh::BulkData::destroy_shared_aura" ;

  assert_ok_to_modify( method );

  Ghosting & aura = * m_ghosting[0] ;

  internal_change_ghosting( aura, std::vector<EntityProc>(), aura.receive() );
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

