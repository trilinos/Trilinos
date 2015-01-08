// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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
#include <stk_topology/topology.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/SameType.hpp>
#include <stk_util/util/StaticAssert.hpp>
#include "stk_util/parallel/DistributedIndex.hpp"  // for DistributedIndex, etc

#include <vector>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {
namespace impl {

void find_elements_these_nodes_have_in_common(BulkData& mesh, unsigned numNodes, const Entity* nodes, std::vector<Entity>& elems)
{
  elems.clear();
  std::vector<Entity> tmp;
  for(unsigned i=0; i<numNodes; ++i) {
    const Entity* elements = mesh.begin_elements(nodes[i]);
    unsigned numElements = mesh.num_elements(nodes[i]);
    tmp.assign(elements, elements+numElements);
    std::sort(tmp.begin(), tmp.end());
    if (i==0) {
      elems.assign(tmp.begin(), tmp.end());
    }
    else {
       std::vector<Entity> intersect;
       std::back_insert_iterator<std::vector<Entity> > intersect_itr(intersect);
       std::set_intersection(elems.begin(), elems.end(),
                             tmp.begin(), tmp.end(),
                             intersect_itr);
       elems.swap(intersect);
    }
  }
}

void find_locally_owned_elements_these_nodes_have_in_common(BulkData& mesh, unsigned numNodes, const Entity* nodes, std::vector<Entity>& elems)
{
  find_elements_these_nodes_have_in_common(mesh, numNodes, nodes, elems);

  for(int i=elems.size()-1; i>=0; --i) {
    if (!mesh.bucket(elems[i]).owned()) {
      elems.erase(elems.begin()+i);
    }
  }
}

bool find_element_edge_ordinal_and_equivalent_nodes(BulkData& mesh, Entity element, unsigned numEdgeNodes, const Entity* edgeNodes, unsigned& elemEdgeOrdinal, Entity* elemEdgeNodes)
{
  stk::topology elemTopology = mesh.bucket(element).topology();
  stk::topology edgeTopology = elemTopology.edge_topology();
  const Entity* elemNodes = mesh.begin_nodes(element);
  ThrowAssertMsg(mesh.num_nodes(element) == elemTopology.num_nodes(), "findElementEdgeOrdinalAndNodes ERROR, element (id="<<mesh.identifier(element)<<") has wrong number of connected nodes ("<<mesh.num_nodes(element)<<"), expected elemTopology.num_nodes()="<<elemTopology.num_nodes());

  unsigned numEdgesPerElem = elemTopology.num_edges();
  for(elemEdgeOrdinal=0; elemEdgeOrdinal<numEdgesPerElem; ++elemEdgeOrdinal) {
    elemTopology.edge_nodes(elemNodes, elemEdgeOrdinal, elemEdgeNodes);
    if (edgeTopology.equivalent(edgeNodes, elemEdgeNodes).first) {
      //found element edge equivalent to edgeNodes.
      //output arguments elemEdgeOrdinal and elemEdgeNodes are set, let's get out of here.
      return true;
    }
  }

  return false;//didn't find element edge equivalent to input edgeNodes
}

int check_for_connected_nodes(const BulkData& mesh)
{
  //This function returns an old-fashioned int return-code which is 0 if all is well
  //and -1 if an error is found.
  //
  //All EDGE_RANK and FACE_RANK entities must have at least 1 connected node.
  for(stk::mesh::EntityRank rank=stk::topology::EDGE_RANK; rank<=stk::topology::ELEMENT_RANK; ++rank) {
    const stk::mesh::BucketVector& buckets = mesh.buckets(rank);
    for(size_t i=0; i<buckets.size(); ++i) {
      const stk::mesh::Bucket& bucket = *buckets[i];
      if (bucket.topology() == stk::topology::INVALID_TOPOLOGY)
      {
        std::cerr << "Entities on rank " << rank << " bucket " << i << " have no topology defined" << std::endl;
        return -1;
      }
      for(size_t j=0; j<bucket.size(); ++j) {
        if (bucket.num_nodes(j) < 1) {
          std::cerr << "Entity with rank="<<rank<<", identifier="<<mesh.identifier(bucket[j])<<" has no connected nodes."<<std::endl;
          return -1;
        }
        // NEED TO CHECK FOR EACH BUCKET INHABITANT THAT ALL ITS NODES ARE VALID.
        unsigned num_nodes = bucket.num_nodes(j);
        Entity const* nodes = bucket.begin_nodes(j);
        for (unsigned k = 0; k < num_nodes; ++k) {
          if (!mesh.is_valid(nodes[k])) {
            std::cerr << "Entity with rank="<<rank<<", identifier="<<mesh.identifier(bucket[j])<<" is connected to an invalid node."
                      << " via node relation " << k << std::endl;
            return -1;
          }
        }
      }
    }
  }
  return 0;
}

int check_no_shared_elements_or_higher(const BulkData& mesh)
{
  for(stk::mesh::EntityRank rank=stk::topology::ELEMENT_RANK; rank < mesh.mesh_meta_data().entity_rank_count(); ++rank) {
    const stk::mesh::BucketVector& buckets = mesh.buckets(rank);
    for(size_t j=0; j<buckets.size(); ++j) {
      if (buckets[j]->size() > 0 && buckets[j]->shared()) {
        stk::mesh::Entity entity = (*buckets[j])[0];
        std::cerr << "Entities with rank ELEMENT_RANK or higher must not be shared. Entity with rank="<<rank<<", identifier="<<mesh.identifier(entity)<<" is shared."<<std::endl;
        return -1;
      }
    }
  }
  return 0;
}



void connectEntityToEdge(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Entity entity,
        stk::mesh::Entity edge, const stk::mesh::Entity* nodes, size_t numNodes)
{
    // get node entity ids
    ThrowRequireMsg(numNodes==2,"connectEntityToEdge ERROR, numNodes must be 2 currently.");
    std::vector<stk::mesh::EntityId> nodeIdsForEdge(numNodes);
    nodeIdsForEdge[0] = stkMeshBulkData.identifier(nodes[0]);
    nodeIdsForEdge[1] = stkMeshBulkData.identifier(nodes[1]);

    // scratch space
    stk::mesh::OrdinalVector ordinal_scratch;
    ordinal_scratch.reserve(64);
    stk::mesh::PartVector part_scratch;
    part_scratch.reserve(64);
    stk::mesh::Permutation perm = static_cast<stk::mesh::Permutation>(0);

    // now what
    stk::topology entity_topology = stkMeshBulkData.bucket(entity).topology();
    std::vector<stk::mesh::EntityId> nodeIds(2);
    std::vector<stk::mesh::EntityId> entityNodes(entity_topology.num_nodes());
    unsigned edge_ordinal = 100000;
    stk::mesh::Entity const * elem_nodes = stkMeshBulkData.begin_nodes(entity);
    for (size_t k=0;k<entity_topology.num_nodes();k++)
    {
        entityNodes[k] =stkMeshBulkData.identifier(elem_nodes[k]);
    }

    for (size_t k=0;k<entity_topology.num_edges();k++)
    {
        entity_topology.edge_nodes(entityNodes, k, nodeIds.begin());
        stk::topology edge_top = entity_topology.edge_topology();
        if ( edge_top.equivalent(nodeIds, nodeIdsForEdge).first )
        {
            edge_ordinal = k;
            break;
        }
    }
    ThrowRequireMsg(edge_ordinal !=100000, "Program error. Contact sierra-help for support.");
    stkMeshBulkData.declare_relation(entity, edge, edge_ordinal, perm, ordinal_scratch, part_scratch);
}

//----------------------------------------------------------------------

// Given a vector of local ownership changes, remove duplicates and
// sanity check.
void internal_clean_and_verify_parallel_change(
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

void internal_generate_parallel_change_lists( const BulkData & mesh ,
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


void get_ghost_data( const BulkData& bulkData, Entity entity, std::vector<EntityGhostData> & dataVector )
{
    // Check locally owned:
    if (bulkData.bucket(entity).member(bulkData.mesh_meta_data().locally_owned_part())) {
        EntityGhostData data;
        data.bulkData = &bulkData;
        data.entity = entity;
        data.direction = EntityGhostData::NONE;
        data.ghostingLevel = EntityGhostData::LOCALLY_OWNED;
        ThrowAssert( bulkData.parallel_rank() == bulkData.parallel_owner_rank(entity) );
        data.processor = bulkData.parallel_rank();
        dataVector.push_back(data);
    }
    // Look through all the ghosting levels:
    std::vector<Ghosting *> ghosting = bulkData.ghostings();
    for (size_t z=0 ; z<ghosting.size() ; ++z) {
        std::vector<EntityProc> send_list;
        ghosting[z]->send_list(send_list);
        for (size_t si=0 ; si<send_list.size() ; ++si) {
            if (send_list[si].first == entity) {
                EntityGhostData data;
                data.bulkData = &bulkData;
                data.entity = entity;
                data.ghostingLevel = z;
                data.direction = EntityGhostData::SEND;
                data.processor = send_list[si].second;
                dataVector.push_back(data);
            }
        }
        std::vector<EntityKey> recv_list;
        ghosting[z]->receive_list(recv_list);
        for (size_t ri=0 ; ri<recv_list.size() ; ++ri) {
            if (recv_list[ri] == bulkData.entity_key(entity)) {
                EntityGhostData data;
                data.bulkData = &bulkData;
                data.entity = entity;
                data.ghostingLevel = z;
                data.direction = EntityGhostData::RECEIVE;
                data.processor = bulkData.parallel_owner_rank(entity);
                dataVector.push_back(data);
            }
        }
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

void unpack_not_owned_verify_compare_comm_info( CommBuffer&            buf,
                                                const BulkData &       mesh,
                                                Entity                 entity,
                                                EntityKey &            recv_entity_key,
                                                int       &            recv_owner_rank,
                                                unsigned  &            recv_comm_count,
                                                std::vector<Part*>&    recv_parts,
                                                std::vector<Relation>& recv_relations,
                                                std::vector<int>    &  recv_comm,
                                                bool&                  bad_comm)
{
  const EntityKey key = mesh.entity_key(entity);
  const PairIterEntityComm ec = mesh.entity_comm_map(key);
  const unsigned ec_size = ec.size();
  std::vector<unsigned> ec_idx_shared;
  std::vector<unsigned> ec_idx_not_shared;
  for (unsigned iec=0; iec < ec_size; iec++) {
    if (0 == ec[iec].ghost_id) {
      ec_idx_shared.push_back(iec);
    }
    else {
      ec_idx_not_shared.push_back(iec);
    }
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

    if (ghost_after_shared_count) {
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

void unpack_not_owned_verify_compare_closure_relations( const BulkData &             mesh,
                                                        Entity                       entity,
                                                        std::vector<Relation> const& recv_relations,
                                                        bool&                        bad_rel)
{
    const Bucket & bucket = mesh.bucket(entity);
    const Ordinal bucket_ordinal = mesh.bucket_ordinal(entity);
    EntityRank erank = mesh.entity_rank(entity);
    const EntityRank end_rank = static_cast<EntityRank>(MetaData::get(mesh).entity_rank_count());

    EntityRank irank = stk::topology::BEGIN_RANK;

    Entity const *rels_itr = bucket.begin(bucket_ordinal, irank);
    Entity const *rels_end = bucket.end(bucket_ordinal, irank);
    ConnectivityOrdinal const *ords_itr = bucket.begin_ordinals(bucket_ordinal, irank);
    Permutation const *perms_itr = bucket.begin_permutations(bucket_ordinal, irank);

    std::vector<Relation>::const_iterator jr = recv_relations.begin();

    for(; !bad_rel && jr != recv_relations.end() &&
            jr->entity_rank() < erank; ++jr, ++rels_itr, ++ords_itr)
    {
        while((rels_itr == rels_end) && (irank < end_rank))
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

        if(perms_itr)
        {
            bad_rel = (bad_rel || (static_cast<Permutation>(jr->permutation()) != *perms_itr));
            ++perms_itr;
        }
    }
}

void unpack_not_owned_verify_compare_parts(const BulkData &  mesh,
                                           Entity            entity,
                                           PartVector const& recv_parts,
                                           bool&             bad_part)
{
  const MetaData & meta = MetaData::get(mesh);
  Part * const       owns_part   = & meta.locally_owned_part();
  Part * const       shares_part = & meta.globally_shared_part();

  const Bucket & bucket = mesh.bucket(entity);
  std::pair<const unsigned *,const unsigned *>
    part_ordinals = bucket.superset_part_ordinals();

  const unsigned * k = part_ordinals.first ;

  std::vector<Part*>::const_iterator ip = recv_parts.begin();

  for ( ; ! bad_part && ip != recv_parts.end() ; ++ip ) {
    while ((k != part_ordinals.second) &&
            (!meta.get_part(*k).entity_membership_is_parallel_consistent())) {
      ++k;
    }
    if ( owns_part != *ip ) {
      if ( shares_part != *ip && (*ip)->entity_membership_is_parallel_consistent() ) {
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

void unpack_not_owned_verify_report_errors(const BulkData& mesh,
                                           Entity entity,
                                           bool bad_key,
                                           bool bad_own,
                                           bool bad_part,
                                           bool bad_rel,
                                           bool bad_comm,
                                           EntityKey            recv_entity_key,
                                           int                  recv_owner_rank,
                                           std::vector<Part*> const&    recv_parts,
                                           std::vector<Relation> const& recv_relations,
                                           std::vector<int>    const&  recv_comm,
                                           std::ostream & error_log)
{
  const int p_rank = mesh.parallel_rank();

  const MetaData & meta = MetaData::get(mesh);
  const PartVector & mesh_parts  = meta.get_parts();
  const Bucket & bucket = mesh.bucket(entity);
  const Ordinal bucket_ordinal = mesh.bucket_ordinal(entity);
  const EntityRank erank = mesh.entity_rank(entity);
  const EntityKey key = mesh.entity_key(entity);
  std::pair<const unsigned *,const unsigned *>
    part_ordinals = bucket.superset_part_ordinals();

  error_log << __FILE__ << ":" << __LINE__ << ": ";
  error_log << "P" << p_rank << ": " ;
  error_log << key;
  error_log << " owner(" << mesh.parallel_owner_rank(entity) << ")" ;

  if ( bad_key || bad_own ) {
    error_log << " != received " ;
    error_log << recv_entity_key;
    error_log << " owner(" << recv_owner_rank
              << ")" << std::endl ;
  }
  else if ( bad_comm ) {
    const PairIterEntityComm ec = mesh.entity_comm_map(key);
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

    for ( std::vector<Part*>::const_iterator
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
        error_log << " " << irank<<":"<<mesh.identifier(*ir_itr) ;
      }
    }
    error_log << " ) != received Relations(" ;
    std::vector<Relation>::const_iterator jr = recv_relations.begin() ;
    for ( ; jr != recv_relations.end() &&
            jr->entity_rank() < erank ; ++jr ) {
      error_log << " " << jr->entity_rank()<<":"<<mesh.identifier(jr->entity()) ;
      Entity const * nodes_begin = mesh.begin_nodes(jr->entity());
      Entity const * nodes_end   = mesh.end_nodes(jr->entity());
      error_log << " node-connectivity (";
      for (Entity const* nodeId = nodes_begin; nodeId != nodes_end; ++nodeId)
      {
          error_log << mesh.identifier(*nodeId) << ", ";
      }
      error_log << ") ";
    }
    error_log << " )" << std::endl ;
  }
}


//----------------------------------------------------------------------------
// Unpacking all of my not-owned entities.

bool unpack_not_owned_verify( CommAll & comm_all ,
                              const BulkData & mesh ,
                              std::ostream & error_log )
{
  const int               p_rank = mesh.parallel_rank();
  const EntityCommListInfoVector & entity_comm = mesh.comm_list();

#if (defined(DEBUG_PRINT_COMM_LIST)  && defined(DEBUG_PRINT_COMM_LIST_UNPACK))
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
    ThrowRequire( mesh.entity_key(entity) == key );


    if ( i->owner != p_rank ) {

      CommBuffer & buf = comm_all.recv_buffer( i->owner );

      check_tag(mesh, buf, PACK_TAG_ENTITY_SHARED, PACK_TAG_ENTITY_GHOST);
      unpack_entity_info( buf , mesh ,
                          recv_entity_key , recv_owner_rank ,
                          recv_parts , recv_relations );

      if (mesh.in_shared(key)) {
        check_tag(mesh, buf, PACK_TAG_SHARED_COUNT);
      }
      else {
        check_tag(mesh, buf, PACK_TAG_GHOST_COUNT);
      }
      recv_comm_count = 0 ;
      buf.unpack<unsigned>( recv_comm_count );
      recv_comm.resize( recv_comm_count );
      buf.unpack<int>( & recv_comm[0] , recv_comm_count );

      // Match key and owner

      const bool bad_key = key                              != recv_entity_key ;
      const bool bad_own = mesh.parallel_owner_rank(entity) != recv_owner_rank ;
      bool bad_part = false ;
      bool bad_rel  = false ;
      bool bad_comm = false ;

      bool broken = bad_key || bad_own;

      // Compare communication information:

      if ( ! broken ) {
        unpack_not_owned_verify_compare_comm_info( buf,
                                                   mesh,
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
        unpack_not_owned_verify_compare_parts(mesh,
                                              entity,
                                              recv_parts,
                                              bad_part);
        broken = bad_part;
      }

      // Compare the closure relations:
      if ( ! broken )
      {
        unpack_not_owned_verify_compare_closure_relations( mesh,
                                                           entity,
                                                           recv_relations,
                                                           bad_rel );
        broken = bad_rel;

      }

      // The rest of this code is just error handling
      if ( broken ) {
        unpack_not_owned_verify_report_errors(mesh,
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

void pack_owned_verify( CommAll & all , const BulkData & mesh )
{
  const EntityCommListInfoVector & entity_comm = mesh.comm_list();
  const int p_rank = all.parallel_rank();

  for ( EntityCommListInfoVector::const_iterator
        i = entity_comm.begin() ; i != entity_comm.end() ; ++i ) {

    if ( i->owner == p_rank ) {

      std::vector<int> share_procs ;
      std::vector<int> ghost_procs ;

      const PairIterEntityComm comm = mesh.entity_comm_map(i->key);

      for ( size_t j = 0 ; j < comm.size() ; ++j ) {
        if ( comm[j].ghost_id == stk::mesh::BulkData::SHARED ) {
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

namespace {
bool ordered_comm(const BulkData& bulk, const Entity entity )
{
  const PairIterEntityComm ec = bulk.entity_comm_map(bulk.entity_key(entity));
  const size_t n = ec.size();
  for ( size_t i = 1 ; i < n ; ++i ) {
    if ( ! ( ec[i-1] < ec[i] ) ) {
      return false ;
    }
  }
  return true ;
}
}

bool verify_parallel_attributes_for_bucket( BulkData& M, Bucket const& bucket, std::ostream & error_log, size_t& comm_count )
{
  const int p_rank = M.parallel_rank();

  bool result = true;

  const MetaData & S = MetaData::get(M);
  Part & owns_part = S.locally_owned_part();
  Part & shares_part = S.globally_shared_part();

  const bool has_owns_part   = has_superset( bucket , owns_part );
  const bool has_shares_part = has_superset( bucket , shares_part );

  const Bucket::iterator j_end = bucket.end();
  Bucket::iterator j           = bucket.begin();

  while ( j != j_end ) {
    Entity entity = *j ; ++j ;

    bool this_result = true;

    const int      p_owner    = M.parallel_owner_rank(entity);
    const bool     ordered    = ordered_comm(M, entity );
    const bool     shares     = M.in_shared( M.entity_key(entity) );
    const bool     recv_ghost = M.in_receive_ghost( M.entity_key(entity) );
    const bool     send_ghost = M.in_send_ghost( M.entity_key(entity) );
    const bool     owned_closure = M.owned_closure(entity);

    if ( ! ordered ) {
      error_log << __FILE__ << ":" << __LINE__ << ": ";
      error_log << "Problem is unordered" << std::endl;
      this_result = false ;
    }

    // Owner consistency:

    if (   has_owns_part && p_owner != p_rank ) {
      error_log << __FILE__ << ":" << __LINE__ << ": ";
      error_log << "problem is owner-consistency check 1: "
                << "has_owns_part: " << has_owns_part << ", "
                << "p_owner: " << p_owner << ", "
                << "p_rank: " << p_rank << std::endl;
      this_result = false ;
    }

    if ( ! has_owns_part && p_owner == p_rank ) {
      error_log << __FILE__ << ":" << __LINE__ << ": ";
      error_log << "problem is owner-consistency check 2: "
                << "has_owns_part: " << has_owns_part << ", "
                << "p_owner: " << p_owner << ", "
                << "p_rank: " << p_rank << std::endl;
      this_result = false ;
    }

    if ( has_shares_part != shares ) {
      error_log << __FILE__ << ":" << __LINE__ << ": ";
      error_log << "problem is owner-consistency check 3: "
                << "has_shares_part: " << has_shares_part << ", "
                << "shares: " << shares << " has entity key " << M.entity_key(entity).m_value << std::endl;
      this_result = false ;
    }

    // Definition of 'closure'

    if ( ( has_owns_part || has_shares_part ) != owned_closure ) {
      error_log << __FILE__ << ":" << __LINE__ << ": ";
      error_log << "problem is closure check: "
                << "has_owns_part: " << has_owns_part << ", "
                << "has_shares_part: " << has_shares_part << ", "
                << "owned_closure: " << owned_closure << std::endl;
      this_result = false ;
    }

    // Must be either owned_closure or recv_ghost but not both.

    if (   owned_closure &&   recv_ghost ) {
      error_log << __FILE__ << ":" << __LINE__ << ": ";
      error_log << "problem: entity is both recv ghost and in owned_closure;"<<std::endl;
      this_result = false ;
    }
    if ( ! owned_closure && ! recv_ghost ) {
      error_log << __FILE__ << ":" << __LINE__ << ": ";
      error_log << "problem: entity is neither recv_ghost nor in owned_closure;"<<std::endl;
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
      M.comm_shared_procs(M.entity_key(entity),shared_procs);
      std::vector<int>::const_iterator it = std::find(shared_procs.begin(),shared_procs.end(),p_owner);
      if (it == shared_procs.end()) {
        error_log << __FILE__ << ":" << __LINE__ << ": ";
        error_log << "problem: entity shared-not-owned, but comm_shared_procs does not contain owner;" << std::endl;
        this_result = false ;
      }
    }

    if ( shares || recv_ghost || send_ghost ) { ++comm_count ; }

    if ( ! this_result ) {
      result = false ;
      error_log << __FILE__ << ":" << __LINE__ << ": ";
      error_log << "P" << M.parallel_rank() << ": " << " entity " << M.entity_rank(entity)<< ",id="<<M.identifier(entity);
      error_log << " details: owner(" << p_owner<<"), node-connectivity(";
      const Entity* nodes = M.begin_nodes(entity);
      unsigned num_nodes = M.num_nodes(entity);
      for(unsigned i=0; i<num_nodes; ++i) {
        error_log<<M.identifier(nodes[i])<<" ";
      }
      error_log<<"), comm(";
      PairIterEntityComm ip = M.entity_comm_map(M.entity_key(entity));
      for ( ; ! ip.empty() ; ++ip ) {
        error_log << " ghost_id=" << ip->ghost_id << ":proc=" << ip->proc ;
      }
      error_log << " )" << std::endl ;
    }
  }

  return result;
}

bool verify_parallel_attributes_comm_list_info( BulkData & M , size_t comm_count, std::ostream & error_log )
{
  bool result = true;

  for ( EntityCommListInfoVector::const_iterator
        i =  M.comm_list().begin() ;
        i != M.comm_list().end() ; ++i ) {

    const PairIterEntityComm ec = M.entity_comm_map(i->key);

    if ( ec.empty() ) {
      error_log << __FILE__ << ":" << __LINE__ << ": ";
      error_log << i->key.id();
      error_log << " ERROR: in entity_comm but has no comm info" << std::endl ;
      result = false ;
    }

    if (i->key != M.entity_key(i->entity)) {
      error_log << __FILE__ << ":" << __LINE__ << ": ";
      error_log << i->key.id();
      error_log << " ERROR: out of sync entity keys in comm list, real key is " << M.entity_key(i->entity).id() << std::endl ;
      result = false ;
    }

    if (i->owner != M.parallel_owner_rank(i->entity)) {
      error_log << __FILE__ << ":" << __LINE__ << ": ";
      error_log << i->key.id();
      error_log << " ERROR: out of sync owners, in comm-info " << i->owner << ", in entity " << M.parallel_owner_rank(i->entity) << std::endl ;
      result = false ;
    }
  }

  if ( M.comm_list().size() != comm_count ) {
    error_log << __FILE__ << ":" << __LINE__ << ": ";
    error_log << " ERROR: entity_comm.size() = " << M.comm_list().size();
    error_log << " != " << comm_count << " = entities with comm info" ;
    error_log << std::endl ;
    result = false ;
  }

  return result;
}

bool verify_parallel_attributes( BulkData & M , std::ostream & error_log )
{
  bool result = true ;

  const size_t EntityRankEnd = MetaData::get(M).entity_rank_count();

  size_t comm_count = 0 ;

  for ( size_t itype = 0 ; itype < EntityRankEnd ; ++itype ) {
    const BucketVector & all_buckets = M.buckets( static_cast<EntityRank>(itype) );

    const BucketVector::const_iterator i_end = all_buckets.end();
          BucketVector::const_iterator i     = all_buckets.begin();

    while ( i != i_end ) {
      Bucket & bucket = **i ; ++i ;

      result = result && verify_parallel_attributes_for_bucket(M, bucket, error_log, comm_count);
    }
  }

  result = result && verify_parallel_attributes_comm_list_info(M, comm_count, error_log);

  return result ;
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


// Enforce that shared entities must be in the owned closure:

void destroy_dependent_ghosts( BulkData & mesh , Entity entity )
{
  EntityRank entity_rank = mesh.entity_rank(entity);

  const EntityRank end_rank = static_cast<EntityRank>(mesh.mesh_meta_data().entity_rank_count());
  EntityVector temp_entities;
  Entity const* rels = NULL;
  int num_rels = 0;

  for (EntityRank irank = static_cast<EntityRank>(end_rank - 1); irank > entity_rank; --irank)
  {
    bool canOneHaveConnectivityFromEntityRankToIrank = mesh.connectivity_map().valid(entity_rank, irank);
    if (canOneHaveConnectivityFromEntityRankToIrank)
    {
      num_rels = mesh.num_connectivity(entity, irank);
      rels     = mesh.begin(entity, irank);
    }
    else
    {
      num_rels = get_connectivity(mesh, entity, irank, temp_entities);
      rels     = &*temp_entities.begin();
    }

    for (int r = num_rels - 1; r >= 0; --r)
    {
      Entity e = rels[r];

      bool upwardRelationOfEntityIsInClosure = mesh.owned_closure(e);
      ThrowRequireMsg( !upwardRelationOfEntityIsInClosure, mesh.entity_rank(e) << " with id " << mesh.identifier(e) << " should not be in closure." );

      // Recursion
      destroy_dependent_ghosts( mesh , e );
    }
  }

  mesh.destroy_entity( entity );
}
}

// Entities with sharing information that are not in the owned closure
// have been modified such that they are no longer shared.
// These may no longer be needed or may become ghost entities.
// There is not enough information so assume they are to be deleted
// and let these entities be re-ghosted if they are needed.

// Open question: Should an owned and shared entity that does not
// have an upward relation to an owned entity be destroyed so that
// ownership transfers to another process?

void delete_shared_entities_which_are_no_longer_in_owned_closure( BulkData & mesh )
{
  for ( EntityCommListInfoVector::const_reverse_iterator
        i =  mesh.comm_list().rbegin() ;
        i != mesh.comm_list().rend() ; ++i) {

    Entity entity = i->entity;

    bool entityisValid = mesh.is_valid(entity);
    std::vector<int> shared_procs;
    mesh.comm_shared_procs(i->key,shared_procs);
    bool isSharedEntity = !shared_procs.empty();
    bool isNotInOwnedClosure = !mesh.owned_closure(entity);
    bool entityIsSharedButNotInClosure =  entityisValid && isSharedEntity && isNotInOwnedClosure;

    if ( entityIsSharedButNotInClosure )
    {
      destroy_dependent_ghosts( mesh , entity );
    }
  }
}

bool shared_entities_modified_on_any_proc(const BulkData& mesh, stk::ParallelMachine comm)
{
    Selector shared = mesh.mesh_meta_data().globally_shared_part();
    bool local_any_shared_entities_modified = false;
    for(stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank <= stk::topology::ELEM_RANK; ++rank ) {
        const stk::mesh::BucketVector& buckets = rank==stk::topology::ELEM_RANK ? mesh.buckets(rank) : mesh.get_buckets(rank, shared);
        for(size_t b=0; b<buckets.size(); ++b) {
            const stk::mesh::Bucket& bkt = *buckets[b];
            for(size_t i=0; i<bkt.size(); ++i) {
                if (mesh.state(bkt[i]) == Modified) {
                    if (rank == stk::topology::ELEM_RANK) {
                        unsigned num_nodes = mesh.num_nodes(bkt[i]);
                        const stk::mesh::Entity* nodes = mesh.begin_nodes(bkt[i]);
                        for(unsigned j=0; j<num_nodes; ++j) {
                            if (mesh.bucket(nodes[j]).shared()) {
                                local_any_shared_entities_modified = true;
                                break;
                            }
                        }
                    }
                    else { 
                        local_any_shared_entities_modified = true;
                        break;
                    }
                }
            }
            if (local_any_shared_entities_modified) {
                break;
            }
        }
        if (local_any_shared_entities_modified) {
            break;
        }
    }

    int local_shared_modified = local_any_shared_entities_modified ? 1 : 0;
    int global_shared_modified = 0;
    stk::all_reduce_max(comm, &local_shared_modified, &global_shared_modified, 1);
    return global_shared_modified > 0;
}

stk::parallel::DistributedIndex::KeySpanVector
convert_entity_keys_to_spans( const MetaData & meta )
{
  // Make sure the distributed index can handle the EntityKey

  enum { OK = StaticAssert< SameType< uint64_t, parallel::DistributedIndex::KeyType >::value >::OK };

  // Default constructed EntityKey has all bits set.

  const EntityKey invalid_key ;
  const EntityId  min_id = 1 ;
  const EntityId  max_id = invalid_key.id();

  const EntityRank rank_count = static_cast<EntityRank>(meta.entity_rank_count());

  parallel::DistributedIndex::KeySpanVector spans( rank_count );

  for ( EntityRank rank = stk::topology::NODE_RANK ; rank < rank_count ; ++rank ) {
    EntityKey key_min( rank , min_id );
    EntityKey key_max( rank , max_id );
    spans[rank].first  = key_min;
    spans[rank].second = key_max;
  }

  return spans ;
}

void internal_fix_node_sharing_delete_on_2015_03_06(stk::mesh::BulkData& bulk_data)
{
    stk::CommAll comm(bulk_data.parallel());

    for (int phase=0;phase<2;++phase)
    {
        for (int i=0;i<bulk_data.parallel_size();++i)
        {
            if ( i != bulk_data.parallel_rank() )
            {
                const stk::mesh::BucketVector& buckets = bulk_data.buckets(stk::topology::NODE_RANK);
                for (size_t j=0;j<buckets.size();++j)
                {
                    const stk::mesh::Bucket& bucket = *buckets[j];
                    if ( bucket.owned() )
                    {
                        for (size_t k=0;k<bucket.size();++k)
                        {
                            stk::mesh::EntityKey key = bulk_data.entity_key(bucket[k]);
                            comm.send_buffer(i).pack<stk::mesh::EntityKey>(key);
                        }
                    }
                }
            }
        }

        if (phase == 0 )
        {
            comm.allocate_buffers( bulk_data.parallel_size()/4 );
        }
        else
        {
            comm.communicate();
        }
    }

    for (int i=0;i<bulk_data.parallel_size();++i)
    {
        if ( i != bulk_data.parallel_rank() )
        {
            while(comm.recv_buffer(i).remaining())
            {
                stk::mesh::EntityKey key;
                comm.recv_buffer(i).unpack<stk::mesh::EntityKey>(key);
                stk::mesh::Entity node = bulk_data.get_entity(key);
                if ( bulk_data.is_valid(node) )
                {
                    bulk_data.add_node_sharing(node, i);
                }
            }
        }
    }
}


} // namespace impl
} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

