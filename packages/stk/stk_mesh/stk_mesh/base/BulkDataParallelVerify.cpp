
/**
 * @author H. Carter Edwards
 */

#include <set>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <assert.h>

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

namespace {

bool verify_parallel_attributes( BulkData & M , std::ostream & error_log );

bool verify_shared_comm( ParallelMachine parallel ,
                         const std::vector< EntityProc > & v ,
                         std::ostream & error_log );

bool verify_ghosting_comm( ParallelMachine parallel ,
                           const Ghosting & gh ,
                           std::ostream & error_log );

bool verify_not_owned_entities( BulkData & M , std::ostream & error_log );

}

bool comm_mesh_verify_parallel_consistency(
  BulkData & M ,
  std::ostream & error_log )
{
  int result = 1 ;

  // Verify consistency of parallel attributes

  result = verify_parallel_attributes( M , error_log );

  all_reduce( M.parallel() , ReduceMin<1>( & result ) );

  // Verify consistency of sharing list

  if ( result ) {
    result = verify_shared_comm( M.parallel(), M.shared_entities(), error_log );
  }

  // Verify consistency of ghosting lists

  for ( std::vector<Ghosting*>::const_iterator
        i = M.ghostings().begin() ; 
        i != M.ghostings().end() && result ; ++i ) {
    result = verify_ghosting_comm( M.parallel() , **i , error_log );
  }

  if ( result ) {
    // Lists are all good, verify owner vs. not-owner attributes and relations
    result = verify_not_owned_entities( M , error_log );
  }

  return result ;
}

namespace {

// Query if entity is within a receive ghost list

bool in_recv_ghost( BulkData & M , Entity * e )
{
  const std::vector< Ghosting * > & all_ghosts = M.ghostings();

  bool result = false ;
  for ( std::vector< Ghosting * >::const_iterator
        ig = all_ghosts.begin() ; ! result && ig != all_ghosts.end() ; ++ig ) {
    Ghosting & gh = **ig ;

    std::vector< Entity * >::const_iterator ie =
      std::lower_bound( gh.receive().begin() ,
                        gh.receive().end() , e , EntityLess() );
    result = ie != gh.receive().end() && e == *ie ;
  }
  return result ;
}

bool in_send_ghost( BulkData & M , Entity * e )
{
  const std::vector< Ghosting * > & all_ghosts = M.ghostings();

  bool result = false ;
  for ( std::vector< Ghosting * >::const_iterator
        ig = all_ghosts.begin() ; ! result && ig != all_ghosts.end() ; ++ig ) {
    Ghosting & gh = **ig ;

    std::vector< EntityProc >::const_iterator ie =
      std::lower_bound( gh.send().begin() ,
                        gh.send().end() , e , EntityLess() );
    result = ie != gh.send().end() && e == ie->first ;
  }
  return result ;
}

bool verify_parallel_attributes( BulkData & M , std::ostream & error_log )
{
  bool result = true ;

  const MetaData & S = M.mesh_meta_data();
  Part & uses_part = S.locally_used_part();
  Part & owns_part = S.locally_owned_part();

  const unsigned p_rank = M.parallel_rank();

  const size_t EntityTypeEnd = M.mesh_meta_data().entity_type_count();

  for ( size_t itype = 0 ; itype < EntityTypeEnd ; ++itype ) {
    const std::vector< Bucket * > & all_buckets = M.buckets( itype );

    const std::vector<Bucket*>::const_iterator i_end = all_buckets.end();
          std::vector<Bucket*>::const_iterator i     = all_buckets.begin();

    while ( i != i_end ) {
      Bucket & bucket = **i ; ++i ;

      const bool uses = has_superset( bucket , uses_part );
      const bool owns = has_superset( bucket , owns_part );

      const Bucket::iterator j_end = bucket.end();
            Bucket::iterator j     = bucket.begin();

      while ( j != j_end ) {
        Entity & entity = *j ; ++j ;

        bool this_result = true ;

        const unsigned p_owner =   entity.owner_rank();
        const bool     shares  = ! entity.sharing().empty();
        const bool     recv_ghost = in_recv_ghost( M , & entity );
        const bool     send_ghost = in_send_ghost( M , & entity );

        // Owner consistency:

        if (   owns && p_owner != p_rank ) { this_result = false ; }
        if ( ! owns && p_owner == p_rank ) { this_result = false ; }

        // Owns is defined to be a subset of uses

        if ( owns && ! uses ) { this_result = false ; }

        // Shared is defined to be a subset of uses.

        if ( shares && ! uses ) { this_result = false ; }

        // Must be either uses or recv_ghost but not both.

        if (   uses &&   recv_ghost ) { this_result = false ; }
        if ( ! uses && ! recv_ghost ) { this_result = false ; }

        // If sending as a ghost then I must own it

        if ( ! owns && send_ghost ) { this_result = false ; }
 
        // If shared then I am owner or owner is in the shared list

        if ( shares && p_owner != p_rank ) {
          PairIterEntityProc ip = entity.sharing();
          for ( ; ! ip.empty() && p_owner != ip->second ; ++ip );
          if ( ip.empty() ) { result = false ; }
        }

        if ( ! this_result ) {
          result = false ;
          print_entity_key( error_log , M.mesh_meta_data(), entity.key() );
          error_log << " ERROR: owner(" << p_owner << ")"
                    << ") owns(" << owns
                    << ") uses(" << uses
                    << ") recv_ghost(" << recv_ghost
                    << ") send_ghost(" << send_ghost
                    << ") shares(" ;
          PairIterEntityProc ip = entity.sharing();
          for ( ; ! ip.empty() ; ++ip ) { error_log << " " << ip->second ; }
          error_log << " )" << std::endl ;
        }
      }
    }
  }

  for ( std::vector<EntityProc>::const_iterator
        i =  M.shared_entities().begin() ;
        i != M.shared_entities().end() ; ) {

    PairIterEntityProc ip = i->first->sharing();

    for ( ; ! ip.empty() && i == ip.first ; ++ip , ++i );

    if ( ! ip.empty() ) {
      error_log << " ERROR: Shared list subsetting" << std::endl ;
      result = false ;
    }
  }

  return result ;
}

// Unpacking all of my not-owned entities.

int unpack_entities_verify( CommAll & comm_all ,
                            const BulkData & mesh ,
                            const std::vector<Entity*> & recv_list ,
                            std::ostream & error_log )
{
  const MetaData & meta = mesh.mesh_meta_data();
  Part * const       owns_part  = & meta.locally_owned_part();
  Part * const       uses_part  = & meta.locally_used_part();
  const PartVector & mesh_parts = meta.get_parts();

  int result = 0 ;

  EntityKey             recv_entity_key ;
  unsigned              recv_owner_rank = 0 ;
  std::vector<Part*>    recv_parts ;
  std::vector<Relation> recv_relations ;

  for ( std::vector<Entity*>::const_iterator
        i = recv_list.begin() ; i != recv_list.end() ; ) {

    Entity * const entity = *i ;
    const Bucket & bucket = entity->bucket();
    std::pair<const unsigned *,const unsigned *>
      part_ordinals = bucket.superset_part_ordinals();

    CommBuffer & buf = comm_all.recv_buffer( entity->owner_rank() );

    unpack_entity_info( buf , mesh ,
                        recv_entity_key , recv_owner_rank ,
                        recv_parts , recv_relations );

    // Match key and owner

    const bool bad_key = entity->key()        != recv_entity_key ;
    const bool bad_own = entity->owner_rank() != recv_owner_rank ;
    bool bad_part = false ;
    bool bad_rel  = false ;

    // Compare everything but the owns part and uses part

    if ( ! bad_key && ! bad_own ) {

      const unsigned * k = part_ordinals.first ;

      std::vector<Part*>::iterator ip = recv_parts.begin();

      for ( ; ! bad_part && ip != recv_parts.end() ; ++ip ) {
        if ( owns_part != *ip && uses_part != *ip ) {
          bad_part = k == part_ordinals.second ||
                     (*ip)->mesh_meta_data_ordinal() != *k ;
          ++k ;
        }
      }
    }

    // Compare the closure relations:
    if ( ! bad_key && ! bad_own && ! bad_part ) {

      PairIterRelation ir = entity->relations();

      std::vector<Relation>::iterator jr = recv_relations.begin() ;

      for ( ; ! bad_rel && jr != recv_relations.end() &&
              jr->entity_type() < entity->entity_type() ; ++jr , ++ir ) {
        bad_rel = ir.empty() || *jr != *ir ;
      }
    }

    if ( bad_key || bad_own || bad_part || bad_rel ) {
      print_entity_key( error_log , meta, entity->key() );
      error_log << " owner(" << entity->owner_rank() << ")" ;

      if ( bad_key || bad_own ) {
        error_log << " != " ;
        print_entity_key( error_log , meta, recv_entity_key );
        error_log << " owner(" << recv_owner_rank
                  << ") received" << std::endl ;
      }
      else if ( bad_part ) {
        error_log << " Parts( " << std::endl ;

        for ( const unsigned * k = part_ordinals.first ;
                               k < part_ordinals.second ; ++k ) {
          error_log << " " << mesh_parts[ *k ]->name();
        }
        error_log << " ) != Parts( " ;

        for ( std::vector<Part*>::iterator
              ip =  recv_parts.begin();
              ip != recv_parts.end() ; ++ip ) {
          error_log << " " << (*ip)->name();
        }
        error_log << " ) received" << std::endl ;
      }
      else if ( bad_rel ) {
        error_log << " Relations(" ;
        PairIterRelation ir = entity->relations();
        for ( ; ! ir.empty() &&
                ir->entity_type() < entity->entity_type() ; ++ir ) {
          print_relation( error_log , meta, ir->attribute() , ir->entity()->key() );
        }
        error_log << " ) != Relations(" ;
        std::vector<Relation>::iterator jr = recv_relations.begin() ;
        for ( ; jr != recv_relations.end() &&
                jr->entity_type() < entity->entity_type() ; ++jr ) {
          print_relation( error_log , meta, jr->attribute() , jr->entity()->key() );
        }
        error_log << " ) received" << std::endl ;
      }
      result = -1 ;
    }
  }

  return result ;
}

// Verify not-owned entities against the owner's entity

bool verify_not_owned_entities( BulkData & M , std::ostream & error_log )
{
  const unsigned p_size = M.parallel_size();
  const unsigned p_rank = M.parallel_rank();

  std::vector< EntityProc > owner_send ;
  std::vector< Entity * > receive ;

  int result = 1 ;

  //----------------------------------------
  // Build not-owned send and receive lists:

  for ( std::vector< EntityProc >::const_iterator
        j =  M.shared_entities().begin() ;
        j != M.shared_entities().end() ; ++j ) {
    if ( p_rank == j->first->owner_rank() ) { owner_send.push_back( *j ); }
    else if ( p_rank == j->second )         { receive.push_back( j->first ); }
  }

  for ( std::vector< Ghosting * >::const_iterator
        j =  M.ghostings().begin() ;
        j != M.ghostings().end() ; ++j ) {
    Ghosting & gh = **j ;
    owner_send.insert( owner_send.end(), gh.send().begin(), gh.send().end() );
    receive.insert( receive.end(), gh.receive().begin(), gh.receive().end());
  }

  {
    std::vector<EntityProc>::iterator i = owner_send.begin();
    std::vector<EntityProc>::iterator e = owner_send.end();
    std::sort( i , e , EntityLess() );
    i = std::unique( i , e );
    owner_send.erase( i , e );
  }
  {
    std::vector<Entity*>::iterator i = receive.begin();
    std::vector<Entity*>::iterator e = receive.end();
    std::sort( i , e , EntityLess() );
    i = std::unique( i , e );
    receive.erase( i , e );
  }
  //----------------------------------------
  // Send owned to not-owned and verify consistency

  CommAll all_info( M.parallel() );

  for ( std::vector<EntityProc>::const_iterator
        i = owner_send.begin() ; i != owner_send.end() ; ++i ) {
    pack_entity_info( all_info.send_buffer( i->second ) , * i->first );
  }

  all_info.allocate_buffers( p_size / 4 , false /* Not symmetric */ );

  for ( std::vector<EntityProc>::const_iterator
        i = owner_send.begin() ; i != owner_send.end() ; ++i ) {
    pack_entity_info( all_info.send_buffer( i->second ) , * i->first );
  }

  all_info.communicate();

  result = unpack_entities_verify( all_info , M , receive , error_log );

  all_reduce( M.parallel() , ReduceMin<1>( & result ) );

  return result ;
}

// The local parallel attributes have already been verified on all processors

bool verify_shared_comm( ParallelMachine parallel ,
                         const std::vector< EntityProc > & v ,
                         std::ostream & error_log )
{
  bool result = true ;
  CommAll all( parallel );

  std::vector< EntityProc >::const_iterator i , j , k , kk ;

  for ( i = v.begin() ; i != v.end() ; i = j ) {
    for ( j = i ; j != v.end() && j->first == i->first ; ++j );
    if ( i->first->owner_rank() == all.parallel_rank() ) {
      for ( k = i ; k != j ; ++k ) {
        CommBuffer & buf = all.send_buffer( k->second );
        buf.skip<EntityKey>(1);
        buf.skip<unsigned>( 1 + std::distance( i , j ) );
      }
    }
  }

  all.allocate_buffers( all.parallel_size() / 4 , false /* Not symmetric */ );

  for ( i = v.begin() ; i != v.end() ; i = j ) {
    for ( j = i ; j != v.end() && j->first == i->first ; ++j );
    if ( i->first->owner_rank() == all.parallel_rank() ) {
      const unsigned count = static_cast<unsigned>(std::distance( i , j ));
      for ( k = i ; k != j ; ++k ) {
        CommBuffer & buf = all.send_buffer( k->second );
        buf.pack<EntityKey>( i->first->key() );
        buf.pack<unsigned>( count );
        // Pack what the receiver should have:
        // My list, remove receiver, add myself
        for ( kk = i ; kk != j && kk->second < all.parallel_rank() ; ++kk ) {
          if ( k != kk ) { buf.pack<unsigned>( kk->second ); }
        }
        buf.pack<unsigned>( all.parallel_rank() );
        for ( ++kk ; kk != j ; ++kk ) {
          if ( k != kk ) { buf.pack<unsigned>( kk->second ); }
        }
      }
    }
  }

  all.communicate();

  for ( i = v.begin() ; i != v.end() ; i = j ) {
    for ( j = i ; j != v.end() && j->first == i->first ; ++j );
    if ( i->first->owner_rank() != all.parallel_rank() ) {

      CommBuffer & buf = all.recv_buffer( i->first->owner_rank() );
      EntityKey key ;
      unsigned count = 0 ;
      unsigned rank = 0 ;

      bool ok = static_cast<int>( sizeof(EntityKey) +
                                  sizeof(unsigned) * (1+std::distance(i,j)) )
                <= buf.remaining();

      if ( ok ) {
        buf.unpack<EntityKey>( & key , 1 );
        buf.unpack<unsigned>( & count , 1 );

        ok = i->first->key() == key &&
             static_cast<int>(count) == std::distance( i , j );

        for ( kk = i ; ok && kk != j ; ++kk ) {
          buf.unpack<unsigned>( & rank , 1 );
          --count ;
          ok = rank == kk->second ;
        }
        buf.skip<unsigned>( count );
      }

      if ( ! ok ) {
        const MetaData & meta_data = i->first->bucket().mesh().mesh_meta_data();
        print_entity_key( error_log , meta_data, i->first->key() );
        error_log << " owner(" << i->first->owner_rank()
                  << ") share#(" << std::distance(i,j)
                  << ") != " ;
        print_entity_key( error_log , meta_data, key );
        error_log << " sharing processor info" << std::endl ;
        result = false ;
      }
    }
  }
  return result ;
}

bool verify_ghosting_comm( ParallelMachine parallel ,
                           const Ghosting & gh ,
                           std::ostream & error_log )
{
  bool result = true ;

  CommAll all( parallel );

  for ( std::vector<EntityProc>::const_iterator
        i = gh.send().begin() ; i != gh.send().end() ; ++i ) {
    all.send_buffer( i->second ).skip<EntityKey>(1);
  }

  all.allocate_buffers( all.parallel_rank() / 4 , false );

  for ( std::vector<EntityProc>::const_iterator
        i = gh.send().begin() ; i != gh.send().end() ; ++i ) {
    all.send_buffer( i->second ).pack<EntityKey>( i->first->key() );
  }

  all.communicate();

  for ( std::vector<Entity*>::const_iterator
        i = gh.receive().begin() ; i != gh.receive().end() ; ++i ) {

    CommBuffer & buf = all.recv_buffer( (*i)->owner_rank() );

    EntityKey key ;

    if ( static_cast<int>( sizeof(EntityKey) ) <= buf.remaining() ) {
      buf.unpack<EntityKey>( key );
    }

    if ( (*i)->key() != key ) {
      const MetaData & meta_data = (*i)->bucket().mesh().mesh_meta_data();
      print_entity_key( error_log , meta_data , (*i)->key() );
      error_log << " owner(" << (*i)->owner_rank() << ") != " ;
      print_entity_key( error_log , meta_data , key );
      error_log << " ERROR in Ghosting(" << gh.name() << ")" << std::endl ;
      print_entity_key( error_log , meta_data , key );
      result = false ;
    }
  }

  return result ;
}

} // namespace<>

} // namespace mesh
} // namespace stk

