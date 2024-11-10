// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Relation.hpp>
#include <stk_topology/topology.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/util/SameType.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_util/util/StaticAssert.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/baseImpl/MeshCommImplUtils.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include <stk_mesh/baseImpl/check_comm_list.hpp>
#include <stk_mesh/base/DumpMeshInfo.hpp>

#include <vector>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {
namespace impl {

void unpack_not_owned_verify_compare_closure_relations( const BulkData &    mesh,
                                               Entity                       entity,
                                               std::vector<Relation> const& recv_relations,
                                               bool&                        bad_rel)
{
    const Bucket & bucket = mesh.bucket(entity);
    const Ordinal bucket_ordinal = mesh.bucket_ordinal(entity);
    const EntityRank end_rank = mesh.entity_rank(entity);

    if (!mesh.in_shared(entity) && !mesh.in_ghost(mesh.aura_ghosting(), entity)) {
      //don't do consistency-check for custom-ghosts
      bad_rel = false;
      return;
    }      

    std::vector<Relation>::const_iterator jr = recv_relations.begin();

    for(EntityRank irank=stk::topology::BEGIN_RANK; !bad_rel && irank<end_rank && jr != recv_relations.end();++irank)
    {
        Entity const *rels_itr = bucket.begin(bucket_ordinal, irank);
        Entity const *rels_end = bucket.end(bucket_ordinal, irank);
        ConnectivityOrdinal const *ords_itr = bucket.begin_ordinals(bucket_ordinal, irank);

        STK_ThrowAssertMsg((rels_itr != rels_end && ords_itr == nullptr) == false, "Relations found without ordinals");      


        for(;rels_itr!=rels_end;++rels_itr,++ords_itr)
        {
          bool is_this_relation_the_same = jr->entity() == *rels_itr;
          bool is_this_ordinal_the_same  = static_cast<ConnectivityOrdinal>(jr->getOrdinal()) == *ords_itr;
          bad_rel = !is_this_relation_the_same || !is_this_ordinal_the_same;
          ++jr;
          if (bad_rel) { 
            break;
          }
        }
        bool recv_relation_still_has_entity_of_irank = jr != recv_relations.end() && jr->entity_rank() == irank;
        bad_rel = bad_rel || recv_relation_still_has_entity_of_irank;
    }
}

void unpack_not_owned_verify_compare_parts(const BulkData &  mesh,
                                           Entity            entity,
                                           PartVector const& recv_parts,
                                           bool&             bad_part)
{
  if (!mesh.in_shared(entity) && !mesh.in_ghost(mesh.aura_ghosting(), entity)) {
    //don't do consistency-check for custom-ghosts
    bad_part = false;
    return;
  }

  std::set<std::string> thisProcExtraParts;
  std::set<std::string> otherProcExtraParts;
  impl::fill_part_list_differences(mesh, entity, recv_parts, thisProcExtraParts, otherProcExtraParts);
  if(!thisProcExtraParts.empty() || !otherProcExtraParts.empty())
  {
      bad_part = true;
  }
}

bool verify_parallel_attributes_for_bucket(const Bucket& bucket,
                      const std::function<PairIterEntityComm(Entity)>& getEntityComm,
                                           std::ostream& error_log)
{
  const BulkData& mesh = bucket.mesh();
  const MetaData& meta = mesh.mesh_meta_data();
  const int p_rank = mesh.parallel_rank();

  bool result = true;

  Part & owns_part = meta.locally_owned_part();
  Part & shares_part = meta.globally_shared_part();

  const bool has_owns_part   = has_superset( bucket , owns_part );
  const bool has_shares_part = has_superset( bucket , shares_part );

  const Bucket::iterator j_end = bucket.end();
  Bucket::iterator j           = bucket.begin();

  while ( j != j_end ) {
    Entity entity = *j ; ++j ;

    bool this_result = true;

    const int      p_owner    = mesh.parallel_owner_rank(entity);
    const bool     ordered    = impl::is_comm_ordered(getEntityComm(entity));
    const bool     shares     = mesh.in_shared( entity );
    const bool     recv_aura = mesh.in_receive_ghost( mesh.aura_ghosting(), entity );
    const bool     recv_any_ghost = mesh.in_receive_ghost( entity );
    const bool     custom_ghost = recv_any_ghost && !recv_aura;
    const bool     send_ghost = mesh.in_send_ghost( entity );
    const bool     ownedClosure = mesh.owned_closure(entity);

    if ( ! ordered ) {
      error_log << __FILE__ << ":" << __LINE__ << ": ";
      error_log << "Problem is unordered" << std::endl;
      this_result = false ;
    }

    // Owner consistency:

    if ( has_owns_part != (p_owner == p_rank) ) {
      error_log << __FILE__ << ":" << __LINE__ << ": ";
      error_log << "problem is owner-consistency (entity in locally-owned part iff owned by current proc): "
                << "has_owns_part: " << (has_owns_part?"true":"false") << ", "
                << "p_owner: " << p_owner << ", "
                << "p_rank: " << p_rank << std::endl;
      this_result = false ;
    }

    if ( has_shares_part != shares ) {
      error_log << __FILE__ << ":" << __LINE__ << ": ";
      error_log << "problem is sharing-consistency (entity in shared part iff it is in comm-list): "
                << "has_shares_part: " << (has_shares_part?"true":"false") << ", "
                << "in comm-list: " << (shares?"true":"false") << ", entity key " << mesh.entity_key(entity) <<" "<<mesh.bucket(entity).topology() << std::endl;
      this_result = false ;
    }

    // Definition of 'closure'

    if (!custom_ghost && (( has_owns_part || has_shares_part ) != ownedClosure)) {
      error_log << __FILE__ << ":" << __LINE__ << ": ";
      error_log << "problem is closure check: "
                << "has_owns_part: " << (has_owns_part?"true":"false") << ", "
                << "has_shares_part: " << (has_shares_part?"true":"false") << ", "
                << "owned_closure: " << (ownedClosure?"true":"false") << std::endl;
      this_result = false ;
    }

    // Must be either owned_closure or recv_aura but not both.


    if (   ownedClosure &&   recv_aura ) {
      error_log << __FILE__ << ":" << __LINE__ << ": ";
      error_log << "problem: entity "<<mesh.entity_key(entity)<<" (with topology "<<bucket.topology()<<") is both recv aura ghost and in owned_closure;"<<std::endl;
      this_result = false ;
    }
    if ( ! ownedClosure && ! recv_any_ghost ) {
      error_log << __FILE__ << ":" << __LINE__ << ": ";
      error_log << "problem: entity is neither a recv ghost nor in owned_closure;"<<std::endl;
      this_result = false ;
    }

    // If sending as a ghost then I must own it

    if ( ! has_owns_part && send_ghost ) {
      error_log << __FILE__ << ":" << __LINE__ << ": ";
      error_log << "problem is send ghost check: "
                << "has_owns_part: " << has_owns_part << ", "
                << "send_ghost: " << send_ghost << std::endl;
      this_result = false ;
    }

    // If shared then I am owner or owner is in the shared list

    if ( shares && p_owner != p_rank ) {
      std::vector<int> shared_procs;
      mesh.comm_shared_procs(mesh.entity_key(entity),shared_procs);
      std::vector<int>::const_iterator it = std::find(shared_procs.begin(),shared_procs.end(),p_owner);
      if (it == shared_procs.end()) {
        error_log << __FILE__ << ":" << __LINE__ << ": ";
        error_log << "problem: entity shared-not-owned, but comm_shared_procs does not contain owner;" << std::endl;
        this_result = false ;
      }
    }

    if ( ! this_result ) {
      result = false ;
      error_log << __FILE__ << ":" << __LINE__ << ": ";
      error_log << "P" << mesh.parallel_rank() << " sync-count="<<mesh.synchronized_count();
      error_log << ": entity " << mesh.entity_key(entity) << " " << mesh.bucket(entity).topology();
      error_log << " details: owner(" << p_owner<<"), shared=" << (bucket.shared() ? "true" : "false");
      error_log << ", aura=" << (bucket.in_aura() ? "true" : "false");
      error_log <<", custom-recv-ghost="<<mesh.in_receive_custom_ghost(mesh.entity_key(entity));
      error_log << std::endl;

      stk::mesh::EntityRank numRanks = static_cast<stk::mesh::EntityRank>(meta.entity_rank_count());
      for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<numRanks; rank++)
      {
          impl::print_connectivity_of_rank(mesh, entity, rank, error_log);
      }

      error_log<<"comm(";
      PairIterEntityComm ip = getEntityComm(entity);
      for ( ; ! ip.empty() ; ++ip ) {
        error_log << " ghost_id=" << ip->ghost_id << ":proc=" << ip->proc ;
      }
      error_log << " )" << std::endl ;
    }
  }

  return result;
}

// these are for debugging, they're used to mark where we are in the packing/unpacking process
//#define USE_PACK_TAGS !defined(NDEBUG)

enum PackTags {
  PACK_TAG_INVALID = 12345600,
  PACK_TAG_SHARED_COUNT,
  PACK_TAG_GHOST_COUNT,
  PACK_TAG_GHOST_COUNT_AFTER_SHARED,
  PACK_TAG_ENTITY_SHARED,
  PACK_TAG_ENTITY_GHOST
};

static bool check_tag(const BulkData& mesh, CommBuffer& buf, PackTags expected_tag, PackTags expected_tag2 = PACK_TAG_INVALID)
{
  bool badTag = false;
#if USE_PACK_TAGS
  int tag = -1;
  try
  {
      buf.unpack<int>(tag);
  }
  catch(...)
  {
      badTag = true;
  }
  badTag = badTag || ( tag != expected_tag && tag != expected_tag2);
  if(badTag)
  {
    std::ostringstream msg;
    msg << "P[" << mesh.parallel_rank() << "] bad tag = " << tag << " expecting " << expected_tag << " or " << expected_tag2 << std::endl;
    std::cerr << msg.str();
  }
#endif
  return badTag;
}

static void put_tag(CommBuffer& buf, PackTags tag)
{
#if USE_PACK_TAGS
  buf.pack<int>(tag);
#endif
}

void pack_owned_verify(const BulkData& mesh,
                       const EntityCommDatabase& commDB,
                       const EntityCommListInfoVector& commList,
                       CommSparse& commSparse)
{
  const int p_rank = commSparse.parallel_rank();

  for (const EntityCommListInfo& info : commList) { 

    if ( mesh.parallel_owner_rank(info.entity) == p_rank ) {

      std::vector<int> share_procs ;
      std::vector<int> ghost_procs ;

      STK_ThrowAssert(info.entity_comm != -1);
      STK_ThrowAssert(info.entity_comm == commDB.entity_comm(mesh.entity_key(info.entity)));
      const PairIterEntityComm comm = commDB.comm(info.entity_comm);

      for ( size_t j = 0 ; j < comm.size() ; ++j ) {
        if ( comm[j].ghost_id == BulkData::SHARED ) {
          // Will be ordered by proc
          share_procs.push_back( comm[j].proc );
        }
        else {
          // No guarantee of ordering by proc
          stk::util::insert_keep_sorted_and_unique(comm[j].proc, ghost_procs);
        }
      }

      const unsigned share_count = share_procs.size();

      for ( size_t j = 0 ; j < share_procs.size() ; ++j ) {

        // Sharing process, send sharing process list

        const int share_proc = share_procs[j] ;

        CommBuffer & buf = commSparse.send_buffer( share_proc );

        put_tag(buf,PACK_TAG_ENTITY_SHARED);

        pack_entity_info(mesh, buf , info.entity );

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
          if ( comm[kk].ghost_id > BulkData::AURA && comm[kk].proc == share_proc ) {
            ++ghost_count ;
          }
        }
        put_tag(buf,PACK_TAG_GHOST_COUNT_AFTER_SHARED);
        buf.pack<unsigned>(ghost_count);
      }

      for ( size_t j = 0 ; j < ghost_procs.size() ; ++j ) {
        const int ghost_proc = ghost_procs[j] ;

        CommBuffer & buf = commSparse.send_buffer( ghost_proc );

        put_tag(buf,PACK_TAG_ENTITY_GHOST);
        pack_entity_info(mesh, buf , info.entity );

        // What ghost subsets go to this process?
        unsigned count = 0 ;
        for ( size_t k = 0 ; k < comm.size() ; ++k ) {
          if ( comm[k].ghost_id != BulkData::SHARED && comm[k].proc == ghost_proc ) {
            ++count ;
          }
        }
        put_tag(buf,PACK_TAG_GHOST_COUNT);
        buf.pack<unsigned>( count );
        for ( size_t k = 0 ; k < comm.size() ; ++k ) {
          if ( comm[k].ghost_id != BulkData::SHARED && comm[k].proc == ghost_proc ) {
            buf.pack<unsigned>( comm[k].ghost_id );
          }
        }
      }
    }
  }
}

void unpack_not_owned_verify_compare_comm_info( const BulkData& mesh,
                      const std::function<PairIterEntityComm(Entity)>& getEntityComm,
                                                CommBuffer&            buf,
                                                Entity                 entity,
                                                EntityKey &            recv_entity_key,
                                                int       &            recv_owner_rank,
                                                unsigned  &            recv_comm_count,
                                                PartVector&    recv_parts,
                                                std::vector<Relation>& recv_relations,
                                                std::vector<int>    &  recv_comm,
                                                bool&                  bad_comm)
{
  if (!mesh.in_shared(entity) && !mesh.in_ghost(mesh.aura_ghosting(), entity)) {
    //don't do consistency-check for custom-ghosts
    bad_comm = false;
    return;
  }

  const PairIterEntityComm ec = getEntityComm(entity);
  const unsigned ec_size = ec.size();
  std::vector<unsigned> ec_idx_shared;
  std::vector<unsigned> ec_idx_not_shared;
  for (unsigned iec=0; iec < ec_size; iec++) {
    if (BulkData::SHARED == ec[iec].ghost_id) {
      ec_idx_shared.push_back(iec);
    }
    else {
      ec_idx_not_shared.push_back(iec);
    }
  }

  //bad_comm = ec_size != recv_comm.size();
  unsigned ghost_after_shared_count=0;
  if ( mesh.in_shared( entity ) ) {
    // only packed shared size, so only compare with shared here
    bad_comm = ec_idx_shared.size() != recv_comm.size();
    if ( ! bad_comm ) {
      size_t j = 0 ;
      for ( ; j < ec_idx_shared.size() &&
              ec[ec_idx_shared[j]].ghost_id == BulkData::SHARED &&
              ec[ec_idx_shared[j]].proc   == recv_comm[j] ; ++j );
      bad_comm = j != ec_idx_shared.size() ;

      // unpack count of additional ghosts
      bad_comm = bad_comm || check_tag(mesh, buf, PACK_TAG_GHOST_COUNT_AFTER_SHARED);
      buf.unpack<unsigned>( ghost_after_shared_count);
    }
  }

  if ( ! bad_comm ) {

    if (ghost_after_shared_count) {
      bad_comm = bad_comm || check_tag(mesh, buf, PACK_TAG_ENTITY_GHOST);
      unpack_entity_info( buf , mesh ,
                          recv_entity_key , recv_owner_rank ,
                          recv_parts , recv_relations );

      bad_comm = bad_comm || check_tag(mesh, buf, PACK_TAG_GHOST_COUNT);
      buf.unpack<unsigned>(recv_comm_count);
      recv_comm.resize( recv_comm_count);
      buf.unpack<int>( recv_comm.data() , recv_comm_count);
    }

    if ( !mesh.in_shared( entity ) || ghost_after_shared_count) {
        // recv_comm contains ghost_ids for ghosted entities
        bad_comm = !impl::all_ghost_ids_are_found_in_comm_data(ec, mesh.parallel_owner_rank(entity), recv_comm);
    }
  }
}

std::string bool_str(bool flag)
{
  return flag ? "true" : "false";
}

void unpack_not_owned_verify_report_errors(const BulkData& mesh,
                      const std::function<PairIterEntityComm(Entity)>& getEntityComm,
                                           Entity entity,
                                           bool bad_key,
                                           bool bad_own,
                                           bool bad_part,
                                           bool bad_rel,
                                           bool bad_comm,
                                           EntityKey            recv_entity_key,
                                           int                  recv_owner_rank,
                                           PartVector const&    recv_parts,
                                           std::vector<Relation> const& recv_relations,
                                           std::vector<int>    const&  recv_comm,
                                           std::ostream & error_log)
{
  const int p_rank = mesh.parallel_rank();

  const Ordinal bucketOrdinal = mesh.bucket_ordinal(entity);
  const EntityRank erank = mesh.entity_rank(entity);
  const EntityKey key = mesh.entity_key(entity);

  error_log << __FILE__ << ":" << __LINE__ << ": ";
  error_log << "P" << p_rank << ": " ;
  error_log << key;
  error_log << " owner(P" << mesh.parallel_owner_rank(entity) << ") shared: " << bool_str(mesh.bucket(entity).shared()) << " in aura: " << bool_str(mesh.bucket(entity).in_aura()) << " ";

  if ( bad_key || bad_own ) {
    error_log << " != received " ;
    error_log << recv_entity_key;
    error_log << " owner(" << recv_owner_rank
              << ")" << std::endl ;
  }
  else if ( bad_comm ) {
    const PairIterEntityComm ec = getEntityComm(entity);
    if ( mesh.in_shared( entity ) ) {
      error_log << " sharing(" ;
      for ( size_t j = 0 ; j < ec.size() &&
              ec[j].ghost_id == BulkData::SHARED ; ++j ) {
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
    error_log << " Comparing parts from this processor(" << mesh.parallel_rank() << ") against processor (" << recv_owner_rank << ")" << std::endl;

    std::set<std::string> thisProcExtraParts;
    std::set<std::string> otherProcExtraParts;
    impl::fill_part_list_differences(mesh, entity, recv_parts, thisProcExtraParts, otherProcExtraParts);

    if ( !thisProcExtraParts.empty() )
    {
        error_log << "\tParts on this proc, not on other proc:" << std::endl;
        std::set<std::string>::iterator iter = thisProcExtraParts.begin();
        for (;iter!=thisProcExtraParts.end();++iter)
        {
            const Part* part = mesh.mesh_meta_data().get_part(*iter);
            error_log << "\t\t" << *iter << ", rank: "<<part->primary_entity_rank()<<std::endl;
        }
    }

    if ( !otherProcExtraParts.empty() )
    {
        error_log << "\tParts on other proc, not on this proc:" << std::endl;
        std::set<std::string>::iterator iter = otherProcExtraParts.begin();
        for (;iter!=otherProcExtraParts.end();++iter)
        {
            error_log << "\t\t" << *iter << std::endl;
        }
    }
  }
  else if ( bad_rel ) {
    error_log << " Relations(" ;
    const Bucket & entityBucket = mesh.bucket(entity);
    for (EntityRank irank = stk::topology::BEGIN_RANK;
         irank < erank; ++irank)
    {
      error_log << " " << irank << ": ";
      Entity const *ir_itr = entityBucket.begin(bucketOrdinal, irank);
      Entity const *ir_end = entityBucket.end(bucketOrdinal, irank);
      for ( ; ir_itr != ir_end; ++ir_itr ) {
        error_log << mesh.identifier(*ir_itr)<<" " ;
        if (irank != stk::topology::NODE_RANK) {
          Entity const * nodes_begin = mesh.begin_nodes(*ir_itr);
          Entity const * nodes_end   = mesh.end_nodes(*ir_itr);
          error_log << "nodes (";
          for (Entity const* nodeId = nodes_begin; nodeId != nodes_end; ++nodeId)
          {
              error_log << mesh.identifier(*nodeId) << ", ";
          }
          error_log << ") ";
        }
      }
    }
    error_log << " ) != received Relations(" ;
    std::vector<Relation>::const_iterator jr = recv_relations.begin() ;
    EntityRank curRank = stk::topology::INVALID_RANK;
    for ( ; jr != recv_relations.end() &&
            jr->entity_rank() < erank ; ++jr ) {
      if (jr->entity_rank() != curRank) {
        error_log << " " << jr->entity_rank()<<": ";
        curRank = jr->entity_rank();
      }
      error_log << mesh.identifier(jr->entity()) << " ";
      if (jr->entity_rank() != stk::topology::NODE_RANK) {
        Entity const * nodes_begin = mesh.begin_nodes(jr->entity());
        Entity const * nodes_end   = mesh.end_nodes(jr->entity());
        error_log << " nodes (";
        for (Entity const* nodeId = nodes_begin; nodeId != nodes_end; ++nodeId)
        {
            error_log << mesh.identifier(*nodeId) << ", ";
        }
        error_log << ") ";
      }
    }
    error_log << " )" << std::endl ;
  }
}

bool unpack_not_owned_verify(const BulkData& mesh,
                             const EntityCommListInfoVector& commList,
                      const std::function<PairIterEntityComm(Entity)>& getEntityComm,
                             CommSparse& commSparse,
                             std::ostream& error_log)
{
  const int p_rank = mesh.parallel_rank();

  bool result = true ;

  EntityKey             recv_entity_key ;
  int                   recv_owner_rank = 0 ;
  unsigned              recv_comm_count = 0 ;
  PartVector    recv_parts ;
  std::vector<Relation> recv_relations ;
  std::vector<int>      recv_comm ;

  for ( EntityCommListInfoVector::const_iterator
        i = commList.begin() ; i != commList.end() ; ++i ) {

    EntityKey key = i->key;
    Entity entity = i->entity;
    STK_ThrowRequire( mesh.entity_key(entity) == key );


    if ( mesh.parallel_owner_rank(i->entity) != p_rank ) {

      bool broken_tag = false;
      CommBuffer & buf = commSparse.recv_buffer( mesh.parallel_owner_rank(i->entity));

      broken_tag = broken_tag || check_tag(mesh, buf, PACK_TAG_ENTITY_SHARED, PACK_TAG_ENTITY_GHOST);
      if(!broken_tag)
      {    
          unpack_entity_info( buf , mesh,
                              recv_entity_key , recv_owner_rank ,
                              recv_parts , recv_relations );

          if (mesh.in_shared(key)) {
              broken_tag = broken_tag || check_tag(mesh, buf, PACK_TAG_SHARED_COUNT);
          }
          else {
              broken_tag = broken_tag || check_tag(mesh, buf, PACK_TAG_GHOST_COUNT);
          }

          if(!broken_tag)
          {
              recv_comm_count = 0 ;
              buf.unpack<unsigned>( recv_comm_count );
              recv_comm.resize( recv_comm_count );
              buf.unpack<int>( recv_comm.data() , recv_comm_count );
          }
      }    

      // Match key and owner

      const bool bad_key = key != recv_entity_key ;
      const bool bad_own = mesh.parallel_owner_rank(entity) != recv_owner_rank ;
      bool bad_part = false ;
      bool bad_rel  = false ;
      bool bad_comm = false ;

      bool broken = broken_tag || bad_key || bad_own;

      // Compare communication information:

      if ( ! broken ) {
        impl::unpack_not_owned_verify_compare_comm_info(mesh,
                                                   getEntityComm,
                                                   buf,
                                                   entity,
                                                   recv_entity_key,
                                                   recv_owner_rank,
                                                   recv_comm_count,
                                                   recv_parts,
                                                   recv_relations,
                                                   recv_comm,
                                                   bad_comm);
        broken = bad_comm;
      }    

      // Compare everything but the owns part and uses part

      if ( ! broken ) {
        impl::unpack_not_owned_verify_compare_parts(mesh, entity,
                                              recv_parts, bad_part);
        broken = bad_part;
      }

      // Compare the closure relations:
      if ( ! broken )
      {
        impl::unpack_not_owned_verify_compare_closure_relations(mesh, entity,
                                                           recv_relations, bad_rel );
        broken = bad_rel;

      }

      // The rest of this code is just error handling
      if ( broken ) {
        impl::unpack_not_owned_verify_report_errors(mesh,
                                                   getEntityComm,
                                              entity,
                                              bad_key,
                                              bad_own,
                                              bad_part,
                                              bad_rel,
                                              bad_comm,
                                              recv_entity_key,
                                              recv_owner_rank,
                                              recv_parts,
                                              recv_relations,
                                              recv_comm,
                                              error_log);
        result = false ;
      }
    }
  }

  return result ;
}

bool verify_parallel_attributes(const BulkData& mesh,
                                const EntityCommDatabase& commDB,
                                const EntityCommListInfoVector& commList,
                                const std::function<PairIterEntityComm(Entity)>& getEntityComm,
                                std::ostream & error_log )
{
  bool result = true ;

  const EntityRank entityRankEnd = static_cast<EntityRank>(mesh.mesh_meta_data().entity_rank_count());

  for ( EntityRank rank = stk::topology::NODE_RANK ; rank < entityRankEnd ; ++rank ) {
    const BucketVector & all_buckets = mesh.buckets(rank);

    for(const Bucket* bucketptr : all_buckets)
    {    
      result = result && impl::verify_parallel_attributes_for_bucket(*bucketptr,
                               getEntityComm, error_log);
    }
  }

  bool isGloballyConsistentCommList = impl::is_comm_list_globally_consistent(mesh, commDB, commList, error_log);
  result = result && isGloballyConsistentCommList;

  return result ;
}

bool comm_mesh_verify_parallel_consistency(const BulkData& mesh,
                                           const EntityCommDatabase& commDB,
                                           const EntityCommListInfoVector& commList,
                                           const std::function<PairIterEntityComm(Entity)>& getEntityComm,
                                           std::ostream & error_log )
{
  int verified_ok = 1 ;

  // Verify consistency of parallel attributes

  verified_ok = verify_parallel_attributes(mesh, commDB, commList, getEntityComm, error_log );
  if (mesh.parallel_size() > 1) {
    all_reduce( mesh.parallel() , ReduceMin<1>( & verified_ok ) );
  }

  // Verify entities against owner.

  if ( verified_ok ) {
    CommSparse comm( mesh.parallel() );

    impl::pack_owned_verify(mesh, commDB, commList, comm);

    comm.allocate_buffers();

    impl::pack_owned_verify(mesh, commDB, commList, comm);

    comm.communicate();

    verified_ok = impl::unpack_not_owned_verify(mesh, commList, getEntityComm, comm , error_log );

    if (mesh.parallel_size() > 1) {
      all_reduce( mesh.parallel() , ReduceMin<1>( & verified_ok ) );
    }      
  }

  return verified_ok == 1 ;
}

void check_matching_parts_count(unsigned partsCount, int rank, int commSize, MPI_Comm comm)
{
  std::vector<unsigned> partsCounts(commSize);

  int status = MPI_Allgather(&partsCount, 1, MPI_UNSIGNED,
                             partsCounts.data(), 1, MPI_UNSIGNED,
                             comm);

  if(status != MPI_SUCCESS) {
    std::string error = "MPI_Allgather FAILED: MPI_Allgather = " + std::to_string(status);
    throw std::runtime_error(error);
  }

  for(auto count : partsCounts) {
    bool hasEqualPartCounts = count == partsCount;
    if(rank == 0) {
      STK_ThrowRequireMsg(hasEqualPartCounts,
                      "Rank: " + std::to_string(rank) +
                      " found unmatching part ordinal counts across procs in calls to change_entity_parts with selector");
    } else {
      STK_ThrowRequire(hasEqualPartCounts);
    }
  }
}

void check_matching_parts(const PartVector& parts, unsigned partsCount, int rank, int commSize, MPI_Comm comm)
{
  std::vector<unsigned> partOrds;
  for(Part* part : parts) {
    stk::util::insert_keep_sorted_and_unique(part->mesh_meta_data_ordinal(), partOrds);
  }

  std::vector<unsigned> recvPartOrds(partsCount * commSize);

  int status = MPI_Allgather(partOrds.data(), partsCount, MPI_UNSIGNED,
                             recvPartOrds.data(), partsCount, MPI_UNSIGNED,
                             comm);

  if(status != MPI_SUCCESS) {
    std::string error = "MPI_Allgather FAILED: MPI_Allgather = " + std::to_string(status);
    throw std::runtime_error(error);
  }

  for(unsigned i = 0; i < recvPartOrds.size(); i++) {
    bool hasEqualPartOrds = recvPartOrds[i] == partOrds[i % partsCount];
    if(rank == 0) {
      STK_ThrowRequireMsg(hasEqualPartOrds,
                      "Rank: " + std::to_string(rank) +
                      " found unmatching part ordinals in change_entity_parts with selector");
    } else {
      STK_ThrowRequire(hasEqualPartOrds);
    }
  }
}

void check_matching_parts_across_procs(const PartVector& parts, MPI_Comm comm)
{
  int rank, commSize;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &commSize);
  unsigned partsCount = parts.size();
  check_matching_parts_count(partsCount, rank, commSize, comm);
  check_matching_parts(parts, partsCount, rank, commSize, comm);
}

void check_matching_selectors_and_parts_across_procs(const Selector& selector,
                                                     const PartVector& add_parts,
                                                     const PartVector& remove_parts,
                                                     MPI_Comm comm)
{
  PartVector selectorParts;
  selector.get_parts(selectorParts);
  check_matching_parts_across_procs(selectorParts, comm);
  check_matching_parts_across_procs(add_parts, comm);
}

} // namespace impl
} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

