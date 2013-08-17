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
#include <stk_mesh/base/EntityCommDatabase.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <fstream>

namespace stk {
namespace mesh {

namespace {

#define DEBUG_PRINT_COMM_LIST 0
#define DEBUG_PRINT_COMM_LIST_UNPACK 0

#if DEBUG_PRINT_COMM_LIST

// Very, very handy for debugging parallel resolution...
static int s_step = 0;

void par_verify_print_comm_list( const BulkData & mesh , bool doit, const std::string message )
{
  ++s_step;
  if ( doit ) {
    std::ostringstream file ;
    file << "comm-list." << s_step << "." << mesh.parallel_rank() << ".dat";
    std::ofstream fout(file.str().c_str());
    std::ostringstream msg ;
    msg << message;
    msg << std::endl ;

    for ( EntityCommListInfoVector::const_iterator
          i =  mesh.comm_list().begin() ;
          i != mesh.comm_list().end() ; ++i ) {

      Entity entity = i->entity;
      msg << "S< " << s_step << " > P" << mesh.parallel_rank() << ": " ;

      print_entity_key( msg , MetaData::get(mesh) , i->key );

      msg << " owner(" << i->owner << ")" ;

      if ( !mesh.is_valid(entity) ) { msg << " del" ; }
      else if ( Modified == mesh.state(entity) ) { msg << " mod" ; }
      else { msg << "    " ; }

      for ( PairIterEntityComm ec = mesh.entity_comm(i->key); ! ec.empty() ; ++ec ) {
        msg << " (" << ec->ghost_id << "," << ec->proc << ")" ;
      }
      msg << std::endl ;
    }

    fout << msg.str();
  }
}

#endif

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

  if (M.parallel_size() > 1) {
    all_reduce( M.parallel() , ReduceMin<1>( & result ) );
  }

  // Verify entities against owner.

  if ( result ) {
    CommAll all( M.parallel() );

    pack_owned_verify( all , M );

    all.allocate_buffers( all.parallel_size() / 4 );

    pack_owned_verify( all , M );

    all.communicate();

    result = unpack_not_owned_verify( all , M , error_log );

    if (M.parallel_size() > 1) {
      all_reduce( M.parallel() , ReduceMin<1>( & result ) );
    }
  }

  return result == 1 ;
}

namespace {

bool ordered_comm(const BulkData& bulk, const Entity entity )
{
  const PairIterEntityComm ec = bulk.entity_comm(bulk.entity_key(entity));
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

  const int p_rank = M.parallel_rank();

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
        Entity entity = *j ; ++j ;

        bool this_result = true ;

        const int      p_owner    = M.parallel_owner_rank(entity);
        const bool     ordered    = ordered_comm(M, entity );
        const bool     shares     = M.in_shared( M.entity_key(entity) );
        const bool     recv_ghost = M.in_receive_ghost( M.entity_key(entity) );
        const bool     send_ghost = M.in_send_ghost( M.entity_key(entity) );
        const bool     owned_closure = in_owned_closure( M, entity , p_rank );

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
          error_log << "problem is closure check: "
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
          error_log << "problem is send ghost check: "
                    << "has_owns_part: " << has_owns_part << ", "
                    << "send_ghost: " << send_ghost << std::endl;
          this_result = false ;
        }

        // If shared then I am owner or owner is in the shared list

        if ( shares && p_owner != p_rank ) {
          PairIterEntityComm ip = M.entity_comm_sharing(M.entity_key(entity));
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
          error_log << M.identifier(entity);
          error_log << " ERROR: owner(" << p_owner
                    << ") owns(" << has_owns_part
                    << ") shares(" << has_shares_part
                    << ") owned_closure(" << owned_closure
                    << ") recv_ghost(" << recv_ghost
                    << ") send_ghost(" << send_ghost
                    << ") comm(" ;
          PairIterEntityComm ip = M.entity_comm(M.entity_key(entity));
          for ( ; ! ip.empty() ; ++ip ) {
            error_log << " " << ip->ghost_id << ":" << ip->proc ;
          }
          error_log << " )" << std::endl ;
        }
      }
    }
  }

  for ( EntityCommListInfoVector::const_iterator
        i =  M.comm_list().begin() ;
        i != M.comm_list().end() ; ++i ) {

    const PairIterEntityComm ec = M.entity_comm(i->key);

    if ( ec.empty() ) {
      error_log << i->key.id();
      error_log << " ERROR: in entity_comm but has no comm info" << std::endl ;
      result = false ;
    }

    if (i->key != M.entity_key(i->entity)) {
      error_log << i->key.id();
      error_log << " ERROR: out of sync entity keys in comm list, real key is " << M.entity_key(i->entity).id() << std::endl ;
      result = false ;
    }

    if (i->owner != M.parallel_owner_rank(i->entity)) {
      error_log << i->key.id();
      error_log << " ERROR: out of sync owners, in comm-info " << i->owner << ", in entity " << M.parallel_owner_rank(i->entity) << std::endl ;
      result = false ;
    }
  }

  if ( M.comm_list().size() != comm_count ) {
    error_log << " ERROR: entity_comm.size() = " << M.comm_list().size();
    error_log << " != " << comm_count << " = entities with comm info" ;
    error_log << std::endl ;
    result = false ;
  }

  return result ;
}

//----------------------------------------------------------------------------
// Packing my owned entities.

void insert( std::vector<int> & vec , int val )
{
  std::vector<int>::iterator j =
    std::lower_bound( vec.begin() , vec.end() , val );
  if ( j == vec.end() || *j != val ) {
    vec.insert( j , val );
  }
}

// these are for debugging, they're used to mark where we are in the packing/unpacking process
#define USE_PACK_TAGS !defined(NDEBUG)
enum PackTags {
  PACK_TAG_INVALID = 12345600,
  PACK_TAG_SHARED_COUNT,
  PACK_TAG_GHOST_COUNT,
  PACK_TAG_GHOST_COUNT_AFTER_SHARED,
  PACK_TAG_ENTITY_SHARED,
  PACK_TAG_ENTITY_GHOST
};

static void check_tag(const BulkData& mesh, CommBuffer& buf, PackTags expected_tag, PackTags expected_tag2 = PACK_TAG_INVALID)
{
#if USE_PACK_TAGS
  int tag;
  buf.unpack<int>(tag);
  if (tag != expected_tag && tag != expected_tag2) {
    std::ostringstream msg;
    msg << "P[" << mesh.parallel_rank() << "] bad tag = " << tag << " expecting " << expected_tag << " or " << expected_tag2;
    ThrowRequireMsg(tag == expected_tag || tag == expected_tag2, msg.str());
  }
#endif
}

static void put_tag(CommBuffer& buf, PackTags tag)
{
#if USE_PACK_TAGS
  buf.pack<int>(tag);
#endif
}

void pack_owned_verify( CommAll & all , const BulkData & mesh )
{
  const EntityCommListInfoVector & entity_comm = mesh.comm_list();
  const int p_rank = all.parallel_rank();

  for ( EntityCommListInfoVector::const_iterator
        i = entity_comm.begin() ; i != entity_comm.end() ; ++i ) {

    if ( i->owner == p_rank ) {

      std::vector<int> share_procs ;
      std::vector<int> ghost_procs ;

      const PairIterEntityComm comm = mesh.entity_comm(i->key);

      for ( size_t j = 0 ; j < comm.size() ; ++j ) {
        if ( comm[j].ghost_id == 0 ) {
          // Will be ordered by proc
          share_procs.push_back( comm[j].proc );
        }
        else {
          // No guarantee of ordering by proc
          insert( ghost_procs , comm[j].proc );
        }
      }

      const unsigned share_count = share_procs.size();

      for ( size_t j = 0 ; j < share_procs.size() ; ++j ) {

        // Sharing process, send sharing process list

        const int share_proc = share_procs[j] ;

        CommBuffer & buf = all.send_buffer( share_proc );

        put_tag(buf,PACK_TAG_ENTITY_SHARED);
        pack_entity_info(mesh, buf , i->entity );

        put_tag(buf,PACK_TAG_SHARED_COUNT);
        buf.pack<unsigned>( share_count );

        // Pack what the receiver should have:
        // My list, remove receiver, add myself
        size_t k = 0 ;
        for ( ; k < share_count && share_procs[k] < p_rank ; ++k ) {
          if ( k != j ) { buf.pack<int>( share_procs[k] ); }
        }
        buf.pack<int>( p_rank );
        for ( ; k < share_count ; ++k ) {
          if ( k != j ) { buf.pack<int>( share_procs[k] ); }
        }

        // see if we also have ghosts
        unsigned ghost_count = 0 ;
        for ( size_t kk = 0 ; kk < comm.size() ; ++kk ) {
          ThrowRequireMsg( !(comm[kk].ghost_id == 1 && comm[kk].proc == share_proc ) ,
                           "error - shouldn't have shared and aura, only shared and custom ghost");
          if ( comm[kk].ghost_id > 1 && comm[kk].proc == share_proc ) {
            ++ghost_count ;
          }
        }
        put_tag(buf,PACK_TAG_GHOST_COUNT_AFTER_SHARED);
        buf.pack<unsigned>(ghost_count);
      }

      for ( size_t j = 0 ; j < ghost_procs.size() ; ++j ) {
        const int ghost_proc = ghost_procs[j] ;

        CommBuffer & buf = all.send_buffer( ghost_proc );

        put_tag(buf,PACK_TAG_ENTITY_GHOST);
        pack_entity_info(mesh, buf , i->entity );

        // What ghost subsets go to this process?
        unsigned count = 0 ;
        for ( size_t k = 0 ; k < comm.size() ; ++k ) {
          if ( comm[k].ghost_id != 0 && comm[k].proc == ghost_proc ) {
            ++count ;
          }
        }
        put_tag(buf,PACK_TAG_GHOST_COUNT);
        buf.pack<unsigned>( count );
        for ( size_t k = 0 ; k < comm.size() ; ++k ) {
          if ( comm[k].ghost_id != 0 && comm[k].proc == ghost_proc ) {
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
  const int               p_rank = mesh.parallel_rank();
  const EntityCommListInfoVector & entity_comm = mesh.comm_list();
  const EntityRank      end_rank = meta.entity_rank_count();

#if DEBUG_PRINT_COMM_LIST && DEBUG_PRINT_COMM_LIST_UNPACK
  par_verify_print_comm_list(mesh, true, "unpack_not_owned_verify");
#endif

  bool result = true ;

  EntityKey             recv_entity_key ;
  int                   recv_owner_rank = 0 ;
  unsigned              recv_comm_count = 0 ;
  std::vector<Part*>    recv_parts ;
  std::vector<Relation> recv_relations ;
  std::vector<int>      recv_comm ;

  for ( EntityCommListInfoVector::const_iterator
        i = entity_comm.begin() ; i != entity_comm.end() ; ++i ) {

    EntityKey key = i->key;
    Entity entity = i->entity;
    EntityRank erank = mesh.entity_rank(entity);

    if ( i->owner != p_rank ) {

      const Bucket & bucket = mesh.bucket(entity);
      const Ordinal bucket_ordinal = mesh.bucket_ordinal(entity);

      std::pair<const unsigned *,const unsigned *>
        part_ordinals = bucket.superset_part_ordinals();

      CommBuffer & buf = comm_all.recv_buffer( i->owner );

      check_tag(mesh, buf, PACK_TAG_ENTITY_SHARED, PACK_TAG_ENTITY_GHOST);
      unpack_entity_info( buf , mesh ,
                          recv_entity_key , recv_owner_rank ,
                          recv_parts , recv_relations );

      if (mesh.in_shared(key))
        check_tag(mesh, buf, PACK_TAG_SHARED_COUNT);
      else
        check_tag(mesh, buf, PACK_TAG_GHOST_COUNT);
      recv_comm_count = 0 ;
      buf.unpack<unsigned>( recv_comm_count );
      recv_comm.resize( recv_comm_count );
      buf.unpack<int>( & recv_comm[0] , recv_comm_count );

      // Match key and owner

      const bool bad_key = key                 != recv_entity_key ;
      const bool bad_own = mesh.parallel_owner_rank(entity) != recv_owner_rank ;
      bool bad_part = false ;
      bool bad_rel  = false ;
      bool bad_comm = false ;

      // Compare communication information:

      if ( ! bad_key && ! bad_own ) {
        const PairIterEntityComm ec = mesh.entity_comm(key);
        const unsigned ec_size = ec.size();
        std::vector<unsigned> ec_idx_shared;
        std::vector<unsigned> ec_idx_not_shared;
        for (unsigned iec=0; iec < ec_size; iec++)
          {
            if (0 == ec[iec].ghost_id)
              ec_idx_shared.push_back(iec);
            else
              ec_idx_not_shared.push_back(iec);
          }
        //bad_comm = ec_size != recv_comm.size();
        unsigned ghost_after_shared_count=0;
        if ( mesh.in_shared( key ) ) {
          // only packed shared size, so only compare with shared here
          bad_comm = ec_idx_shared.size() != recv_comm.size();
          if ( ! bad_comm ) {
            size_t j = 0 ;
            for ( ; j < ec_idx_shared.size() &&
                    ec[ec_idx_shared[j]].ghost_id == 0 &&
                    ec[ec_idx_shared[j]].proc   == recv_comm[j] ; ++j );
            bad_comm = j != ec_idx_shared.size() ;

            // unpack count of additional ghosts
            check_tag(mesh, buf, PACK_TAG_GHOST_COUNT_AFTER_SHARED);
            buf.unpack<unsigned>( ghost_after_shared_count);
          }
        }

        if ( ! bad_comm ) {

          if (ghost_after_shared_count)
            {
              check_tag(mesh, buf, PACK_TAG_ENTITY_GHOST);
              unpack_entity_info( buf , mesh ,
                                  recv_entity_key , recv_owner_rank ,
                                  recv_parts , recv_relations );

              check_tag(mesh, buf, PACK_TAG_GHOST_COUNT);
              buf.unpack<unsigned>(recv_comm_count);
              recv_comm.resize( recv_comm_count);
              buf.unpack<int>( & recv_comm[0] , recv_comm_count);
            }

          if ( !mesh.in_shared( key ) || ghost_after_shared_count) {

            size_t j = 0;
            // recv_comm contains ghost_ids for ghosted entities
            for ( ; j < ec_idx_not_shared.size() &&
                    static_cast<int>(ec[ec_idx_not_shared[j]].ghost_id) == recv_comm[j] &&
                      ec[ec_idx_not_shared[j]].proc   == mesh.parallel_owner_rank(entity) ; ++j );
            bad_comm = j != ec_idx_not_shared.size() ;
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
      if ( ! bad_key && ! bad_own && ! bad_comm && ! bad_part )
      {
        EntityRank irank = stk::topology::BEGIN_RANK;

        Entity const *rels_itr = bucket.begin(bucket_ordinal, irank);
        Entity const *rels_end = bucket.end(bucket_ordinal, irank);
        ConnectivityOrdinal const *ords_itr = bucket.begin_ordinals(bucket_ordinal, irank);
        Permutation const *perms_itr = bucket.begin_permutations(bucket_ordinal, irank);

        std::vector<Relation>::iterator jr = recv_relations.begin() ;

        for ( ; ! bad_rel && jr != recv_relations.end() &&
                jr->entity_rank() < erank; ++jr , ++rels_itr, ++ords_itr )
        {
          while ((rels_itr == rels_end) && (irank < end_rank))
          {
            // There are no more relations of the current, so try the next
            // higher rank if there is one.
            ++irank;
            rels_itr = bucket.begin(bucket_ordinal, irank);
            rels_end = bucket.end(bucket_ordinal, irank);
            ords_itr = bucket.begin_ordinals(bucket_ordinal, irank);
            perms_itr = bucket.begin_permutations(bucket_ordinal, irank);
          }
          bad_rel = (rels_itr == rels_end) || (jr->entity() != *rels_itr)
                    || (static_cast<ConnectivityOrdinal>(jr->getOrdinal()) != *ords_itr);

          if (perms_itr)
          {
            bad_rel = (bad_rel || (static_cast<Permutation>(jr->permutation()) != *perms_itr));
            ++perms_itr;
          }
        }
      }

      // The rest of this code is just error handling
      if ( bad_key || bad_own || bad_comm || bad_part || bad_rel ) {
        error_log << "P" << p_rank << ": " ;
        error_log << key.id();
        error_log << " owner(" << mesh.parallel_owner_rank(entity) << ")" ;

        if ( bad_key || bad_own ) {
          error_log << " != received " ;
          error_log << recv_entity_key.id();
          error_log << " owner(" << recv_owner_rank
                    << ")" << std::endl ;
        }
        else if ( bad_comm ) {
          const PairIterEntityComm ec = mesh.entity_comm(key);
          if ( mesh.in_shared( key ) ) {
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
              error_log << ",p" << mesh.parallel_owner_rank(entity);
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
          for (EntityRank irank = stk::topology::BEGIN_RANK;
                irank < erank; ++irank)
          {
            Entity const *ir_itr = bucket.begin(bucket_ordinal, irank);
            Entity const *ir_end = bucket.end(bucket_ordinal, irank);
            for ( ; ir_itr != ir_end; ++ir_itr ) {
              error_log << " " << *ir_itr ;
            }
          }
          error_log << " ) != received Relations(" ;
          std::vector<Relation>::iterator jr = recv_relations.begin() ;
          for ( ; jr != recv_relations.end() &&
                  jr->entity_rank() < erank ; ++jr ) {
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

