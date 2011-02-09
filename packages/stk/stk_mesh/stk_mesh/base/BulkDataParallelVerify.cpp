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

#include <set>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

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

void pack_owned_verify( CommAll & all , const BulkData & mesh );

bool unpack_not_owned_verify( CommAll & comm_all ,
                              const BulkData & mesh ,
                              std::ostream & error_log );

}

bool comm_mesh_verify_parallel_consistency(
  BulkData & M , std::ostream & error_log )
{
  int result = 1 ;

  // Verify consistency of parallel attributes

  result = verify_parallel_attributes( M , error_log );

  all_reduce( M.parallel() , ReduceMin<1>( & result ) );

  // Verify entities against owner.

  if ( result ) {
    CommAll all( M.parallel() );

    pack_owned_verify( all , M );

    all.allocate_buffers( all.parallel_size() / 4 );

    pack_owned_verify( all , M );

    all.communicate();

    result = unpack_not_owned_verify( all , M , error_log );

    all_reduce( M.parallel() , ReduceMin<1>( & result ) );
  }

  return result == 1 ;
}

namespace {

bool ordered_comm( const Entity & entity )
{
  const PairIterEntityComm ec = entity.comm();
  const size_t n = ec.size();
  for ( size_t i = 1 ; i < n ; ++i ) {
    if ( ! ( ec[i-1] < ec[i] ) ) {
      return false ;
    }
  }
  return true ;
}

bool verify_parallel_attributes( BulkData & M , std::ostream & error_log )
{
  bool result = true ;

  const MetaData & S = MetaData::get(M);
  Part & owns_part = S.locally_owned_part();
  Part & shares_part = S.globally_shared_part();

  const unsigned p_rank = M.parallel_rank();

  const size_t EntityRankEnd = MetaData::get(M).entity_rank_count();

  size_t comm_count = 0 ;

  for ( size_t itype = 0 ; itype < EntityRankEnd ; ++itype ) {
    const std::vector< Bucket * > & all_buckets = M.buckets( itype );

    const std::vector<Bucket*>::const_iterator i_end = all_buckets.end();
          std::vector<Bucket*>::const_iterator i     = all_buckets.begin();

    while ( i != i_end ) {
      Bucket & bucket = **i ; ++i ;

      const bool has_owns_part   = has_superset( bucket , owns_part );
      const bool has_shares_part = has_superset( bucket , shares_part );

      const Bucket::iterator j_end = bucket.end();
            Bucket::iterator j     = bucket.begin();

      while ( j != j_end ) {
        Entity & entity = *j ; ++j ;

        bool this_result = true ;

        const unsigned p_owner    = entity.owner_rank();
        const bool     ordered    = ordered_comm( entity );
        const bool     shares     = in_shared( entity );
        const bool     recv_ghost = in_receive_ghost( entity );
        const bool     send_ghost = in_send_ghost( entity );
        const bool     owned_closure = in_owned_closure( entity , p_rank );

        if ( ! ordered ) {
          error_log << "Problem is unordered" << std::endl;
          this_result = false ;
        }

        // Owner consistency:

        if (   has_owns_part && p_owner != p_rank ) {
          error_log << "problem is owner-consistency check 1: "
                    << "has_owns_part: " << has_owns_part << ", "
                    << "p_owner: " << p_owner << ", "
                    << "p_rank: " << p_rank << std::endl;
          this_result = false ;
        }

        if ( ! has_owns_part && p_owner == p_rank ) {
          error_log << "problem is owner-consistency check 2: "
                    << "has_owns_part: " << has_owns_part << ", "
                    << "p_owner: " << p_owner << ", "
                    << "p_rank: " << p_rank << std::endl;
          this_result = false ;
        }

        if ( has_shares_part != shares ) {
          error_log << "problem is owner-consistency check 3: "
                    << "has_shares_part: " << has_shares_part << ", "
                    << "shares: " << shares << std::endl;
          this_result = false ;
        }

        // Definition of 'closure'

        if ( ( has_owns_part || has_shares_part ) != owned_closure ) {
          error_log << "problem is closure check 1: "
                    << "has_owns_part: " << has_owns_part << ", "
                    << "has_shares_part: " << has_shares_part << ", "
                    << "owned_closure: " << owned_closure << std::endl;
          this_result = false ;
        }

        // Must be either owned_closure or recv_ghost but not both.

        if (   owned_closure &&   recv_ghost ) {
          error_log << "problem is recv ghost check 1: "
                    << "owned_closure: " << owned_closure << ", "
                    << "recv_ghost: " << recv_ghost << std::endl;
          this_result = false ;
        }
        if ( ! owned_closure && ! recv_ghost ) {
          error_log << "problem is recv ghost check 2: "
                    << "owned_closure: " << owned_closure << ", "
                    << "recv_ghost: " << recv_ghost << std::endl;
          this_result = false ;
        }

        // If sending as a ghost then I must own it

        if ( ! has_owns_part && send_ghost ) {
          error_log << "problem is send ghost check 1: "
                    << "has_owns_part: " << has_owns_part << ", "
                    << "send_ghost: " << send_ghost << std::endl;
          this_result = false ;
        }

        // If shared then I am owner or owner is in the shared list

        if ( shares && p_owner != p_rank ) {
          PairIterEntityComm ip = entity.sharing();
          for ( ; ! ip.empty() && p_owner != ip->proc ; ++ip );
          if ( ip.empty() ) {
            error_log << "problem is shared check 1" << std::endl;
            this_result = false ;
          }
        }

        if ( shares || recv_ghost || send_ghost ) { ++comm_count ; }

        if ( ! this_result ) {
          result = false ;
          error_log << "P" << M.parallel_rank() << ": " ;
          print_entity_key( error_log , MetaData::get(M), entity.key() );
          error_log << " ERROR: owner(" << p_owner
                    << ") owns(" << has_owns_part
                    << ") shares(" << has_shares_part
                    << ") owned_closure(" << owned_closure
                    << ") recv_ghost(" << recv_ghost
                    << ") send_ghost(" << send_ghost
                    << ") comm(" ;
          PairIterEntityComm ip = entity.comm();
          for ( ; ! ip.empty() ; ++ip ) {
            error_log << " " << ip->ghost_id << ":" << ip->proc ;
          }
          error_log << " )" << std::endl ;
        }
      }
    }
  }

  for ( std::vector<Entity*>::const_iterator
        i =  M.entity_comm().begin() ;
        i != M.entity_comm().end() ; ++i ) {

    const PairIterEntityComm ec = (*i)->comm();

    if ( ec.empty() ) {
      print_entity_key( error_log , MetaData::get(M), (*i)->key() );
      error_log << " ERROR: in entity_comm but has no comm info" << std::endl ;
      result = false ;
    }
  }

  if ( M.entity_comm().size() != comm_count ) {
    error_log << " ERROR: entity_comm.size() = " << M.entity_comm().size();
    error_log << " != " << comm_count << " = entities with comm info" ;
    error_log << std::endl ;
    result = false ;
  }

  return result ;
}

//----------------------------------------------------------------------------
// Packing my owned entities.

void insert( std::vector<unsigned> & vec , unsigned val )
{
  std::vector<unsigned>::iterator j =
    std::lower_bound( vec.begin() , vec.end() , val );
  if ( j == vec.end() || *j != val ) {
    vec.insert( j , val );
  }
}

void pack_owned_verify( CommAll & all , const BulkData & mesh )
{
  const std::vector<Entity*> & entity_comm = mesh.entity_comm();
  const unsigned p_rank = all.parallel_rank();

  for ( std::vector<Entity*>::const_iterator
        i = entity_comm.begin() ; i != entity_comm.end() ; ++i ) {

    Entity & entity = **i ;

    if ( entity.owner_rank() == p_rank ) {

      std::vector<unsigned> share_proc ;
      std::vector<unsigned> ghost_proc ;

      const PairIterEntityComm comm = entity.comm();

      for ( size_t j = 0 ; j < comm.size() ; ++j ) {
        if ( comm[j].ghost_id == 0 ) {
          // Will be ordered by proc
          share_proc.push_back( comm[j].proc );
        }
        else {
          // No guarantee of ordering by proc
          insert( ghost_proc , comm[j].proc );
        }
      }

      const unsigned share_count = share_proc.size();

      for ( size_t j = 0 ; j < share_proc.size() ; ++j ) {

        // Sharing process, send sharing process list

        const unsigned p = share_proc[j] ;

        CommBuffer & buf = all.send_buffer( p );

        pack_entity_info( buf , entity );

        buf.pack<unsigned>( share_count );

        // Pack what the receiver should have:
        // My list, remove receiver, add myself
        size_t k = 0 ;
        for ( ; k < share_count && share_proc[k] < p_rank ; ++k ) {
          if ( k != j ) { buf.pack<unsigned>( share_proc[k] ); }
        }
        buf.pack<unsigned>( p_rank );
        for ( ; k < share_count ; ++k ) {
          if ( k != j ) { buf.pack<unsigned>( share_proc[k] ); }
        }
      }

      for ( size_t j = 0 ; j < ghost_proc.size() ; ++j ) {
        const unsigned p = ghost_proc[j] ;

        CommBuffer & buf = all.send_buffer( p );

        pack_entity_info( buf , entity );

        // What ghost subsets go to this process?
        unsigned count = 0 ;
        for ( size_t k = 0 ; k < comm.size() ; ++k ) {
          if ( comm[k].ghost_id != 0 && comm[k].proc == p ) {
            ++count ;
          }
        }
        buf.pack<unsigned>( count );
        for ( size_t k = 0 ; k < comm.size() ; ++k ) {
          if ( comm[k].ghost_id != 0 && comm[k].proc == p ) {
            buf.pack<unsigned>( comm[k].ghost_id );
          }
        }
      }
    }
  }
}

//----------------------------------------------------------------------------
// Unpacking all of my not-owned entities.

bool unpack_not_owned_verify( CommAll & comm_all ,
                              const BulkData & mesh ,
                              std::ostream & error_log )
{
  const MetaData & meta = MetaData::get(mesh);
  Part * const       owns_part   = & meta.locally_owned_part();
  Part * const       shares_part = & meta.globally_shared_part();
  const PartVector & mesh_parts  = meta.get_parts();
  const unsigned     p_rank = mesh.parallel_rank();
  const std::vector<Entity*> & entity_comm = mesh.entity_comm();

  bool result = true ;

  EntityKey             recv_entity_key ;
  unsigned              recv_owner_rank = 0 ;
  unsigned              recv_comm_count = 0 ;
  std::vector<Part*>    recv_parts ;
  std::vector<Relation> recv_relations ;
  std::vector<unsigned> recv_comm ;

  for ( std::vector<Entity*>::const_iterator
        i = entity_comm.begin() ; i != entity_comm.end() ; ++i ) {

    Entity & entity = **i ;

    if ( entity.owner_rank() != p_rank ) {

      const Bucket & bucket = entity.bucket();

      std::pair<const unsigned *,const unsigned *>
        part_ordinals = bucket.superset_part_ordinals();

      CommBuffer & buf = comm_all.recv_buffer( entity.owner_rank() );

      unpack_entity_info( buf , mesh ,
                          recv_entity_key , recv_owner_rank ,
                          recv_parts , recv_relations );

      recv_comm_count = 0 ;
      buf.unpack<unsigned>( recv_comm_count );
      recv_comm.resize( recv_comm_count );
      buf.unpack<unsigned>( & recv_comm[0] , recv_comm_count );

      // Match key and owner

      const bool bad_key = entity.key()        != recv_entity_key ;
      const bool bad_own = entity.owner_rank() != recv_owner_rank ;
      bool bad_part = false ;
      bool bad_rel  = false ;
      bool bad_comm = false ;

      // Compare communication information:

      if ( ! bad_key && ! bad_own ) {
        const PairIterEntityComm ec = entity.comm();
        const unsigned ec_size = ec.size();
        bad_comm = ec_size != recv_comm.size();
        if ( ! bad_comm ) {
          size_t j = 0 ;
          if ( in_shared( entity ) ) {
            for ( ; j < ec_size &&
                    ec[j].ghost_id == 0 &&
                    ec[j].proc   == recv_comm[j] ; ++j );
            bad_comm = j != ec_size ;
          }
          else {
            for ( ; j < ec_size &&
                    ec[j].ghost_id == recv_comm[j] &&
                    ec[j].proc   == entity.owner_rank() ; ++j );
            bad_comm = j != ec_size ;
          }
        }
      }

      // Compare everything but the owns part and uses part

      if ( ! bad_key && ! bad_own && ! bad_comm ) {

        const unsigned * k = part_ordinals.first ;

        std::vector<Part*>::iterator ip = recv_parts.begin();

        for ( ; ! bad_part && ip != recv_parts.end() ; ++ip ) {
          if ( owns_part != *ip ) {
            if ( shares_part != *ip ) {
              // All not-owned and not-shares parts must match:
              bad_part = k == part_ordinals.second ||
                         (*ip)->mesh_meta_data_ordinal() != *k ;
              ++k ;
            }
            else if ( k != part_ordinals.second &&
                     *k == shares_part->mesh_meta_data_ordinal() ) {
              // shares-part matches
              ++k ;
            }
          }
        }
      }

      // Compare the closure relations:
      if ( ! bad_key && ! bad_own && ! bad_comm && ! bad_part ) {

        PairIterRelation ir = entity.relations();

        std::vector<Relation>::iterator jr = recv_relations.begin() ;

        for ( ; ! bad_rel && jr != recv_relations.end() &&
                jr->entity_rank() < entity.entity_rank() ; ++jr , ++ir ) {
          bad_rel = ir.empty() || *jr != *ir ;
        }
      }

      if ( bad_key || bad_own || bad_comm || bad_part || bad_rel ) {
        error_log << "P" << p_rank << ": " ;
        print_entity_key( error_log , meta, entity.key() );
        error_log << " owner(" << entity.owner_rank() << ")" ;

        if ( bad_key || bad_own ) {
          error_log << " != received " ;
          print_entity_key( error_log , meta, recv_entity_key );
          error_log << " owner(" << recv_owner_rank
                    << ")" << std::endl ;
        }
        else if ( bad_comm ) {
          const PairIterEntityComm ec = entity.comm();
          if ( in_shared( entity ) ) {
            error_log << " sharing(" ;
            for ( size_t j = 0 ; j < ec.size() &&
                                 ec[j].ghost_id == 0 ; ++j ) {
              error_log << " " << ec[j].proc ;
            }
            error_log << " ) != received sharing(" ;
            for ( size_t j = 0 ; j < recv_comm.size() ; ++j ) {
              error_log << " " << recv_comm[j] ;
            }
            error_log << " )" << std::endl ;
          }
          else {
            error_log << " ghosting(" ;
            for ( size_t j = 0 ; j < ec.size() ; ++j ) {
              error_log << " (g" << ec[j].ghost_id ;
              error_log << ",p" << ec[j].proc ;
              error_log << ")" ;
            }
            error_log << " ) != received ghosting(" ;
            for ( size_t j = 0 ; j < recv_comm.size() ; ++j ) {
              error_log << " (g" << recv_comm[j] ;
              error_log << ",p" << entity.owner_rank();
              error_log << ")" ;
            }
            error_log << " )" << std::endl ;
          }
        }
        else if ( bad_part ) {
          error_log << " Parts( " ;

          for ( const unsigned * k = part_ordinals.first ;
                                 k < part_ordinals.second ; ++k ) {
            error_log << " \"" << mesh_parts[ *k ]->name() << "\"" ;
          }
          error_log << " ) != received Parts( " ;

          for ( std::vector<Part*>::iterator
                ip =  recv_parts.begin();
                ip != recv_parts.end() ; ++ip ) {
            error_log << " \"" << (*ip)->name() << "\"" ;
          }
          error_log << " )" << std::endl ;
        }
        else if ( bad_rel ) {
          error_log << " Relations(" ;
          PairIterRelation ir = entity.relations();
          for ( ; ! ir.empty() &&
                  ir->entity_rank() < entity.entity_rank() ; ++ir ) {
            error_log << " " << *ir ;
          }
          error_log << " ) != received Relations(" ;
          std::vector<Relation>::iterator jr = recv_relations.begin() ;
          for ( ; jr != recv_relations.end() &&
                  jr->entity_rank() < entity.entity_rank() ; ++jr ) {
            error_log << " " << *jr ;
          }
          error_log << " )" << std::endl ;
        }
        result = false ;
      }
    }
  }

  return result ;
}

} // namespace<>

} // namespace mesh
} // namespace stk

