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
#include <stk_mesh/base/FindPermutation.hpp>
#include <stk_topology/topology.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/SameType.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_util/util/StaticAssert.hpp>
#include "stk_util/parallel/DistributedIndex.hpp"  // for DistributedIndex, etc
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include <stk_mesh/baseImpl/EntityGhostData.hpp>
#include <stk_mesh/baseImpl/Visitors.hpp>

#include <vector>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {
namespace impl {

void find_entities_with_larger_ids_these_nodes_have_in_common_and_locally_owned(stk::mesh::EntityId id, const BulkData& mesh, stk::mesh::EntityRank rank, unsigned numNodes, const Entity* nodes, std::vector<Entity>& elementsInCommon)
{
  find_entities_these_nodes_have_in_common_and(mesh, rank, numNodes, nodes, elementsInCommon,
  [&](const Entity& entity) {
    return mesh.identifier(entity) > id && mesh.bucket(entity).owned();
  });
}

const EntityCommListInfo& find_entity(const BulkData& mesh,
                               const EntityCommListInfoVector& entities,
                               const EntityKey& key)
{
  EntityCommListInfoVector::const_iterator lb_itr = std::lower_bound(entities.begin(), entities.end(), key);
  STK_ThrowAssertMsg(lb_itr != entities.end() && lb_itr->key == key,
                 "proc " << mesh.parallel_rank() << " Cannot find entity-key " << key << " in comm-list" );
  return *lb_itr;
}

bool do_these_nodes_have_any_shell_elements_in_common(BulkData& mesh, unsigned numNodes, const Entity* nodes)
{
  std::vector<Entity> elems;
  find_entities_these_nodes_have_in_common_and(mesh, stk::topology::ELEMENT_RANK, numNodes, nodes, elems,
  [&](const Entity& entity) {
    return mesh.bucket(entity).topology().is_shell();
  });
  return !elems.empty();
}

void find_locally_owned_elements_these_nodes_have_in_common(const BulkData& mesh, unsigned numNodes, const Entity* nodes, std::vector<Entity>& elems)
{
  find_entities_these_nodes_have_in_common_and(mesh, stk::topology::ELEMENT_RANK, numNodes, nodes, elems,
  [&](const Entity& entity) {
    return mesh.bucket(entity).owned();
  });
}

bool find_element_edge_ordinal_and_equivalent_nodes(BulkData& mesh, Entity element, unsigned numEdgeNodes, const Entity* edgeNodes, unsigned& elemEdgeOrdinal, Entity* elemEdgeNodes)
{
  stk::topology elemTopology = mesh.bucket(element).topology();
  const Entity* elemNodes = mesh.begin_nodes(element);
  STK_ThrowAssertMsg(mesh.num_nodes(element) == elemTopology.num_nodes(), "findElementEdgeOrdinalAndNodes ERROR, element (id="<<mesh.identifier(element)<<") has wrong number of connected nodes ("<<mesh.num_nodes(element)<<"), expected elemTopology.num_nodes()="<<elemTopology.num_nodes());

  unsigned numEdgesPerElem = elemTopology.num_edges();
  for(elemEdgeOrdinal=0; elemEdgeOrdinal<numEdgesPerElem; ++elemEdgeOrdinal) {
    elemTopology.edge_nodes(elemNodes, elemEdgeOrdinal, elemEdgeNodes);
    if (stk::mesh::is_edge_equivalent(mesh, element, elemEdgeOrdinal, edgeNodes)) {
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
      if (bucket.topology() == stk::topology::INVALID_TOPOLOGY && bucket.size() > 0)
      {
        std::cerr << "Entities on rank " << rank << " bucket " << i << " have no topology defined" << std::endl;
        return -1;
      }
      for(size_t j=0; j<bucket.size(); ++j) {
        if (mesh.is_valid(bucket[j])) {
          if(bucket.num_nodes(j) < 1) {
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

void delete_upward_relations(stk::mesh::BulkData& bulkData,
                             const stk::mesh::Entity& entity)
{
  stk::mesh::EntityRank entity_rank = bulkData.entity_rank(entity);
  const stk::mesh::EntityRank end_rank = static_cast<stk::mesh::EntityRank>(bulkData.mesh_meta_data().entity_rank_count() - 1);
  const stk::mesh::EntityRank begin_rank = static_cast<stk::mesh::EntityRank>(entity_rank);

  for(stk::mesh::EntityRank irank = end_rank; irank != begin_rank; --irank)
  {
    int num_conn = bulkData.num_connectivity(entity, irank);
    for(int j = num_conn - 1; j >= 0; --j)
    {
      const stk::mesh::Entity* rel_entities = bulkData.begin(entity, irank);
      const stk::mesh::ConnectivityOrdinal* rel_ordinals = bulkData.begin_ordinals(entity, irank);
      if(bulkData.is_valid(rel_entities[j]) && bulkData.state(rel_entities[j]) != Deleted)
      {
        bool relationDestoryed = bulkData.destroy_relation(rel_entities[j], entity, rel_ordinals[j]);
        STK_ThrowRequireWithSierraHelpMsg(relationDestoryed);
      }
    }
  }
}

void delete_entities_and_upward_relations(stk::mesh::BulkData &bulkData, const stk::mesh::EntityVector &entities)
{
    for (size_t i = 0; i < entities.size(); ++i)
    {
        const stk::mesh::Entity entity = entities[i];

        if(!bulkData.is_valid(entity)) {
          continue;
        }

        delete_upward_relations(bulkData, entity);

        bool successfully_destroyed = bulkData.destroy_entity(entity);
        STK_ThrowRequireWithSierraHelpMsg(successfully_destroyed);
    }
}

void connectUpwardEntityToEntity(stk::mesh::BulkData& mesh, stk::mesh::Entity upward_entity,
        stk::mesh::Entity entity, const stk::mesh::Entity* nodes)
{
    uint num_nodes = mesh.num_nodes(entity);
    EntityRank entity_rank = mesh.entity_rank(entity);
    stk::topology baseEntityTopology = mesh.bucket(entity).topology();

    // scratch space
    stk::mesh::OrdinalVector scratch1, scratch2, scratch3;
    stk::mesh::Permutation perm = stk::mesh::Permutation::INVALID_PERMUTATION;

    stk::topology upward_entity_topology = mesh.bucket(upward_entity).topology();
    std::vector<stk::mesh::Entity> nodes_of_this_side(num_nodes);
    unsigned entity_ordinal = 100000;
    stk::mesh::Entity const * upward_entity_nodes = mesh.begin_nodes(upward_entity);

    stk::topology entity_top;
    for (size_t k=0;k<upward_entity_topology.num_sub_topology(entity_rank);k++)
    {
        if(entity_rank == stk::topology::EDGE_RANK)
        {
          entity_top = upward_entity_topology.edge_topology(k);
          nodes_of_this_side.resize(entity_top.num_nodes());
          upward_entity_topology.edge_nodes(upward_entity_nodes, k, nodes_of_this_side.data());
        }
        else
        {
          entity_top = upward_entity_topology.face_topology(k);
          nodes_of_this_side.resize(entity_top.num_nodes());
          upward_entity_topology.face_nodes(upward_entity_nodes, k, nodes_of_this_side.data());
        }

        if (baseEntityTopology != entity_top) continue;

        if ( entity_top.is_equivalent(nodes, nodes_of_this_side.data()).is_equivalent )
        {
            entity_ordinal = k;
            break;
        }
    }
    STK_ThrowRequireMsg(entity_ordinal !=100000, "Program error. Contact sierra-help for support.");
    if ((entity_rank > stk::topology::NODE_RANK) && (mesh.entity_rank(upward_entity) > entity_rank))
    {
        perm = stk::mesh::find_permutation(mesh, upward_entity_topology, upward_entity_nodes, entity_top, nodes, entity_ordinal);
        STK_ThrowRequireMsg(perm != INVALID_PERMUTATION, "find_permutation could not find permutation that produces a match");
    }
    mesh.declare_relation(upward_entity, entity, entity_ordinal, perm, scratch1, scratch2, scratch3);
}

//----------------------------------------------------------------------

// Given a vector of local ownership changes, remove duplicates and
// sanity check.
bool internal_clean_and_verify_parallel_change(
  const BulkData & mesh ,
  std::vector<EntityProc> & local_change )
{
  const int             p_rank = mesh.parallel_rank();
  const int             p_size = mesh.parallel_size();
  const ParallelMachine p_comm = mesh.parallel();

  size_t error_and_count[2] = {0, 0} ;

  std::ostringstream error_msg ;

  // Order and eliminate redundancies:
  stk::util::sort_and_unique(local_change, EntityLess(mesh));

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
      ++error_and_count[0];

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

  // Filter out non-changes (entity will be NULL
  {
    std::vector<EntityProc>::iterator i = local_change.begin(),
                                      j = local_change.end();
    i = std::remove( i , j , EntityProc(Entity(), 0) );
    local_change.erase( i , j );
  }

  error_and_count[1] = local_change.size();

  all_reduce( p_comm , ReduceSum<2>( error_and_count ) );

  if ( error_and_count[0] ) {
    all_write_string( p_comm , std::cerr , error_msg.str() );

    STK_ThrowErrorMsg("Bad change ownership directives\n");
  }

  return error_and_count[1] > 0;
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

  CommSparse comm( mesh.parallel() );

  std::vector<int> procs ;

  // pack and communicate change owner information to all
  // processes that know about the entity
  for ( int phase = 0; phase < 2; ++phase) {
    for ( std::vector<EntityProc>::const_iterator
          ip = local_change.begin() ; ip != local_change.end() ; ++ip ) {
      Entity entity      = ip->first ;
      int new_owner = ip->second;
      mesh.comm_procs( entity , procs );
      for ( std::vector<int>::iterator
            j = procs.begin() ; j != procs.end() ; ++j )
      {
        comm.send_buffer( *j )
          .pack<EntityKey>( mesh.entity_key(entity) )
          .pack<int>(  new_owner );
      }
    }
    if (phase == 0) { // allocation phase
      comm.allocate_buffers();
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

      if ( mesh.in_receive_ghost( key ) ) {
        ghosted_change.push_back( entry );
      }

      if ( mesh.in_shared(key) )
      {
        shared_change.push_back( entry );
      }
    }
  }

  std::sort( shared_change.begin() , shared_change.end() , EntityLess(mesh) );
  std::sort( ghosted_change.begin() , ghosted_change.end() , EntityLess(mesh) );
}

stk::mesh::EntityVector convert_keys_to_entities(stk::mesh::BulkData &bulk, const std::vector<stk::mesh::EntityKey>& node_keys)
{
    stk::mesh::EntityVector nodes(node_keys.size());
    for (size_t i=0;i<nodes.size();++i)
        nodes[i] = bulk.get_entity(node_keys[i]);
    return nodes;
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
        STK_ThrowAssert( bulkData.parallel_rank() == bulkData.parallel_owner_rank(entity) );
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

void find_side_nodes(BulkData& mesh, Entity element, int side_ordinal, EntityVector & permuted_side_nodes)
{
    stk::topology elemTopology = mesh.bucket(element).topology();
    stk::topology sideTopology = elemTopology.side_topology(side_ordinal);
    const size_t num_elem_nodes = elemTopology.num_nodes();

    std::vector<EntityId> elem_node_ids(num_elem_nodes);
    Entity const *elem_nodes = mesh.begin_nodes(element);
    STK_ThrowRequire(mesh.num_nodes(element) == num_elem_nodes);
    for (size_t n=0; n<num_elem_nodes; ++n) {
        elem_node_ids[n] = mesh.identifier(elem_nodes[n]);
    }

    // Use node identifier instead of node local_offset for cross-processor consistency.
    typedef std::vector<EntityId>  EntityIdVector;
    EntityIdVector side_node_ids(sideTopology.num_nodes());
    elemTopology.side_nodes(elem_node_ids.data(), side_ordinal, side_node_ids.data());
    unsigned smallest_permutation;
    permuted_side_nodes.resize(sideTopology.num_nodes());
    //if this is a shell OR these nodes are connected to a shell
    EntityVector side_nodes(sideTopology.num_nodes());
    for (unsigned count=0 ; count<sideTopology.num_nodes() ; ++count) {
        side_nodes[count] = mesh.get_entity(stk::topology::NODE_RANK,side_node_ids[count]);
    }
    bool is_connected_to_shell = do_these_nodes_have_any_shell_elements_in_common(mesh,sideTopology.num_nodes(),side_nodes.data());

    if (elemTopology.is_shell() || is_connected_to_shell) {

        EntityIdVector element_node_id_vector(sideTopology.num_nodes());
        EntityIdVector element_node_ordinal_vector(sideTopology.num_nodes());
        EntityVector element_node_vector(sideTopology.num_nodes());
        elemTopology.side_node_ordinals(side_ordinal, element_node_ordinal_vector.data());
        for (unsigned count = 0; count < sideTopology.num_nodes(); ++count) {
            element_node_vector[count] = mesh.begin_nodes(element)[element_node_ordinal_vector[count]];
            element_node_id_vector[count] = mesh.identifier(element_node_vector[count]);
        }
        smallest_permutation = sideTopology.lexicographical_smallest_permutation_preserve_polarity(side_node_ids.data(), element_node_id_vector.data());
        sideTopology.permutation_nodes(element_node_vector.data(), smallest_permutation, permuted_side_nodes.data());
    }
    else {
        smallest_permutation = sideTopology.lexicographical_smallest_permutation(side_node_ids.data());
        EntityVector non_shell_side_nodes(sideTopology.num_nodes());
        elemTopology.side_nodes(elem_nodes, side_ordinal, non_shell_side_nodes.data());
        sideTopology.permutation_nodes(non_shell_side_nodes.data(), smallest_permutation, permuted_side_nodes.data());
    }
}

void get_part_ordinals_to_induce_on_lower_ranks_except_for_omits(const BulkData       & mesh,
                             const Bucket&          bucket_from,
                             const OrdinalVector  & omit,
                                   EntityRank       entity_rank_to,
                                   OrdinalVector  & induced_parts)
{
  const bool dont_check_owner     = mesh.parallel_size() == 1; // critical for fmwk

  // Only induce parts for normal (not back) relations. Can only trust
  // 'entity_from' to be accurate if it is owned by the local process.
  if ( dont_check_owner || bucket_from.owned() ) {
    const EntityRank entity_rank_from = bucket_from.entity_rank();
    STK_ThrowAssert(entity_rank_from > entity_rank_to);

    const stk::mesh::PartVector &superset_parts = bucket_from.supersets();

    // Contributions of the 'from' entity:
    for ( const Part* part : superset_parts ) {
      unsigned partOrd = part->mesh_meta_data_ordinal();
      if ( part->should_induce(entity_rank_from) && ! contains_ordinal( omit.begin(), omit.end() , partOrd )) {
          stk::util::insert_keep_sorted_and_unique(partOrd, induced_parts);
      }
    }
  }
}

void get_part_ordinals_to_induce_on_lower_ranks(const BulkData& mesh,
                                                const Bucket& bucket_from,
                                                EntityRank       entity_rank_to,
                                                OrdinalVector  & induced_parts)
{
  const bool dont_check_owner     = mesh.parallel_size() == 1; // critical for fmwk

  // Only induce parts for normal (not back) relations. Can only trust
  // 'entity_from' to be accurate if it is owned by the local process.
  if ( dont_check_owner || bucket_from.owned() ) {
    const EntityRank entity_rank_from = bucket_from.entity_rank();
    STK_ThrowAssert(entity_rank_from > entity_rank_to);

    const stk::mesh::PartVector &superset_parts = bucket_from.supersets();

    // Contributions of the 'from' entity:
    for ( const Part* part : superset_parts ) {
      if ( part->should_induce(entity_rank_from)) {
          unsigned partOrd = part->mesh_meta_data_ordinal();
          stk::util::insert_keep_sorted_and_unique(partOrd, induced_parts);
      }
    }
  }
}

template<typename PARTVECTOR>
Entity connect_element_to_entity(BulkData & mesh, Entity elem, Entity entity,
        const unsigned relationOrdinal, const PARTVECTOR& parts, stk::topology entity_top)
{
    stk::topology elem_top = mesh.bucket(elem).topology();

    OrdinalVector entity_node_ordinals(entity_top.num_nodes());
    elem_top.sub_topology_node_ordinals(mesh.entity_rank(entity), relationOrdinal, entity_node_ordinals.data());

    const stk::mesh::Entity *elem_nodes = mesh.begin_nodes(elem);
    EntityVector entity_top_nodes(entity_top.num_nodes());
    elem_top.sub_topology_nodes(elem_nodes, mesh.entity_rank(entity), relationOrdinal, entity_top_nodes.data());

    Permutation perm = stk::mesh::find_permutation(mesh, elem_top, elem_nodes, entity_top, entity_top_nodes.data(), relationOrdinal);

    OrdinalVector scratch1, scratch2, scratch3;

    PARTVECTOR initialParts;
    initialParts.reserve(parts.size() + 1);
    initialParts = parts;
    initialParts.push_back(&mesh.mesh_meta_data().get_topology_root_part(entity_top));
    mesh.change_entity_parts(entity, initialParts);

    const stk::mesh::ConnectivityOrdinal *side_ordinals = mesh.begin_ordinals(elem, mesh.entity_rank(entity));
    unsigned num_sides = mesh.count_valid_connectivity(elem, mesh.entity_rank(entity));

    bool elem_to_side_exists = false;
    for(unsigned i = 0; i < num_sides; ++i)
    {
        if(side_ordinals[i] == relationOrdinal)
        {
            elem_to_side_exists = true;
            break;
        }
    }

    if(!elem_to_side_exists)
    {
        mesh.declare_relation(elem, entity, relationOrdinal, perm, scratch1, scratch2, scratch3);
    }

    const unsigned num_side_nodes = mesh.count_valid_connectivity(entity, stk::topology::NODE_RANK);
    if(0 == num_side_nodes)
    {
        Permutation node_perm = stk::mesh::Permutation::INVALID_PERMUTATION;
        Entity const *elem_nodes_local = mesh.begin_nodes(elem);
        for(unsigned i = 0; i < entity_top.num_nodes(); ++i)
        {
            Entity node = elem_nodes_local[entity_node_ordinals[i]];
            mesh.declare_relation(entity, node, i, node_perm, scratch1, scratch2, scratch3);
        }
    }
    else
    {
        STK_ThrowAssertMsg(num_side_nodes == entity_top.num_nodes(),
                "declare_element_to_entity: " << mesh.entity_key(entity) << " already exists with different number of nodes.");
    }

    return entity;
}

template
Entity connect_element_to_entity<stk::mesh::PartVector>(BulkData & mesh, Entity elem, Entity entity,
        const unsigned relationOrdinal, const stk::mesh::PartVector& parts, stk::topology entity_top);
template
Entity connect_element_to_entity<stk::mesh::ConstPartVector>(BulkData & mesh, Entity elem, Entity entity,
        const unsigned relationOrdinal, const stk::mesh::ConstPartVector& parts, stk::topology entity_top);


stk::mesh::Entity get_or_create_face_at_element_side(stk::mesh::BulkData & bulk,
                                                     stk::mesh::Entity elem,
                                                     int side_ordinal,
                                                     stk::mesh::EntityId new_face_global_id,
                                                     const stk::mesh::PartVector & parts)
{
    stk::mesh::Entity new_face = stk::mesh::Entity();
    unsigned elem_num_faces = bulk.num_faces(elem);
    const stk::mesh::Entity * elem_faces = bulk.begin_faces(elem);
    const stk::mesh::ConnectivityOrdinal * elem_ord_it = bulk.begin_face_ordinals(elem);
    for (unsigned i=0 ; i<elem_num_faces ; ++i) {
        if (elem_ord_it[i] == static_cast<unsigned>(side_ordinal)) {
            new_face = elem_faces[i];
            break;
        }
    }
    if (!bulk.is_valid(new_face)) {
        new_face = bulk.declare_element_side(elem, side_ordinal, parts);
    } else {
        bulk.change_entity_parts(new_face, parts );
    }
    return new_face;
}

void create_shell_status(const std::vector<stk::topology> & elements_touching_surface, stk::topology original_element_topology, std::vector<ShellStatus> & element_shell_status) {
    element_shell_status.resize(elements_touching_surface.size(),NO_SHELLS);
    const bool original_element_is_shell = original_element_topology.is_shell();
    bool got_shells = original_element_is_shell;
    if (!original_element_is_shell) {
        for (size_t i=0 ; i<elements_touching_surface.size() ; ++i) {
            if (elements_touching_surface[i].is_shell()) {
                got_shells = true;
                break;
            }
        }
    }
    if (got_shells) {
        for (size_t i=0 ; i<elements_touching_surface.size() ; ++i) {
            if (original_element_is_shell == elements_touching_surface[i].is_shell() ) {
                element_shell_status[i] = YES_SHELLS_BOTH_SHELLS_OR_BOTH_SOLIDS;
            }
            else {
                element_shell_status[i] = YES_SHELLS_ONE_SHELL_ONE_SOLID;
            }
        }
    }
}



void connect_face_to_other_elements(stk::mesh::BulkData & bulk,
                                      stk::mesh::Entity face,
                                      stk::mesh::Entity elem_with_face,
                                      int elem_with_face_side_ordinal
)
{
    stk::topology elem_topology = bulk.bucket(elem_with_face).topology();
    stk::topology side_topology = elem_topology.face_topology(elem_with_face_side_ordinal);
    int num_side_nodes = side_topology.num_nodes();
    std::vector<stk::mesh::Entity> side_nodes(num_side_nodes);
    elem_topology.face_nodes(bulk.begin_nodes(elem_with_face),elem_with_face_side_ordinal,side_nodes.data());
    std::vector<stk::mesh::Entity> common_elements;
    stk::mesh::impl::find_entities_these_nodes_have_in_common(bulk,stk::topology::ELEMENT_RANK,num_side_nodes,side_nodes.data(),common_elements);

    std::vector<stk::topology> element_topology_touching_surface_vector(common_elements.size());
    for (size_t i=0 ; i<common_elements.size() ; ++i) {
        element_topology_touching_surface_vector[i] = bulk.bucket(common_elements[i]).topology();
    }
    std::vector<ShellStatus> element_shell_status;
    create_shell_status(element_topology_touching_surface_vector, elem_topology, element_shell_status);
    for (size_t count=0 ; count<common_elements.size() ; ++count) {
        if (common_elements[count] != elem_with_face) {
            stk::mesh::Entity other_elem = common_elements[count];
            stk::topology other_elem_topology = element_topology_touching_surface_vector[count];
            for (unsigned other_elem_side = 0; other_elem_side < other_elem_topology.num_faces() ; ++other_elem_side) {
                stk::topology other_elem_side_topology = other_elem_topology.face_topology(other_elem_side);
                std::vector<stk::mesh::Entity> other_elem_side_nodes(other_elem_side_topology.num_nodes());
                other_elem_topology.face_nodes(bulk.begin_nodes(other_elem),other_elem_side,other_elem_side_nodes.data());
                if (should_face_be_connected_to_element_side(side_nodes,other_elem_side_nodes,other_elem_side_topology,element_shell_status[count])) {
                    stk::mesh::connect_side_to_element_with_ordinal(bulk, other_elem, face, other_elem_side);
                    break;
                }
            }
        }
    }
}

void pack_downward_relations_if_valid_permutation_exists(stk::mesh::BulkData& mesh, stk::mesh::Entity some_entity, std::vector<stk::mesh::Relation>& recv_relations1)
{
    recv_relations1.clear();
    unsigned bucket_ordinal = mesh.bucket_ordinal(some_entity);
    const stk::mesh::Bucket& bucket = mesh.bucket(some_entity);
    for(EntityRank irank=stk::topology::EDGE_RANK; irank<mesh.entity_rank(some_entity);++irank)
    {
        if (bucket.has_permutation(irank))
        {
            Entity const *rels_itr = bucket.begin(bucket_ordinal, irank);
            Entity const *rels_end = bucket.end(bucket_ordinal, irank);
            stk::mesh::ConnectivityOrdinal const *ords_itr = bucket.begin_ordinals(bucket_ordinal, irank);
            stk::mesh::Permutation const *perms = bucket.begin_permutations(bucket_ordinal, irank);

            for(;rels_itr!=rels_end;++rels_itr,++ords_itr, ++perms)
            {
                if ( *perms != stk::mesh::Permutation::INVALID_PERMUTATION )
                {
                    recv_relations1.push_back(stk::mesh::Relation(*rels_itr, mesh.entity_rank(*rels_itr), *ords_itr));
                    uint32_t perm_attr = static_cast<uint32_t>(*perms);
                    recv_relations1.back().set_attribute(perm_attr);
                }
            }
        }
    }
}

bool check_permutations_on_all(stk::mesh::BulkData& mesh)
{
    std::ostringstream os;
    bool all_ok = true;
    std::vector<stk::mesh::Relation> entity_relations;
    stk::mesh::EntityRank end_rank = static_cast<stk::mesh::EntityRank>(mesh.mesh_meta_data().entity_rank_count());
    for (stk::mesh::EntityRank irank=stk::topology::FACE_RANK; irank<end_rank; ++irank)
    {
        const stk::mesh::BucketVector &buckets = mesh.buckets(irank);
        for (size_t i=0;i<buckets.size();++i)
        {
            const stk::mesh::Bucket &bucket = *buckets[i];
            for (size_t j=0;j<bucket.size();++j)
            {
                stk::mesh::Entity entity = bucket[j];
                pack_downward_relations_if_valid_permutation_exists(mesh, entity, entity_relations);
                for (size_t k=0;k<entity_relations.size();++k)
                {
                    if (mesh.entity_rank(entity) < stk::topology::ELEM_RANK)
                    {
                        bool valid_permutation = stk::mesh::check_permutation(mesh, entity,
                                                           entity_relations[k].entity(),
                                                           static_cast<unsigned>(entity_relations[k].getOrdinal()),
                                                           static_cast<stk::mesh::Permutation>(entity_relations[k].attribute())
                                                                   );
                        if (!valid_permutation)
                        {
                            all_ok = false;
                            os << "For processor " << mesh.parallel_rank() << " entity " << mesh.entity_key(entity) << " did not "
                                    << "get the right permutation for its relationship with entity " << mesh.entity_key(entity_relations[k].entity()) << " with invalid permutation " << entity_relations[k].attribute() <<  std::endl;
                        }
                    }
                }
            }
        }
    }

    if ( !all_ok )
    {
        std::cerr << os.str();
    }


    int verified_ok = all_ok ? 1 : 0;
    if (mesh.parallel_size() > 1)
    {
        all_reduce( mesh.parallel() , ReduceMin<1>( & verified_ok ) );
    }

    all_ok = verified_ok == 1;

    return all_ok;
}

void send_entity_keys_to_owners(
  BulkData & mesh ,
  const std::vector<Entity> & recvGhosts,
        std::set< EntityProc , EntityLess > & sendGhosts)
{
  const int parallel_size = mesh.parallel_size();

  stk::CommSparse sparse( mesh.parallel() );

  // For all entity keys in recvGhosts, send the entity key to the owning proc
  for ( int phase = 0; phase < 2; ++phase) {
    for (Entity recvGhost : recvGhosts) {
      const int owner = mesh.parallel_owner_rank(recvGhost);
      const EntityKey key = mesh.entity_key(recvGhost);
      sparse.send_buffer( owner ).pack<EntityKey>( key );
    }
    if (phase == 0) { //allocation phase
      sparse.allocate_buffers();
    }
    else { //communication phase
      sparse.communicate();
    }
  }

  // Insert into sendGhosts that entities need to be recvd on another proc
  for ( int proc_rank = 0 ; proc_rank < parallel_size ; ++proc_rank ) {
    stk::CommBuffer & buf = sparse.recv_buffer(proc_rank);
    while ( buf.remaining() ) {
      EntityKey key ;
      buf.unpack<EntityKey>( key );
      Entity entity = mesh.get_entity(key);
      EntityProc tmp( entity, proc_rank );
      sendGhosts.insert( tmp );
    }
  }
}

void comm_sync_send_recv(const BulkData & mesh ,
                         EntityProcVec& sendGhosts ,
                         std::set< EntityKey > & recvGhosts )
{
  const int parallel_rank = mesh.parallel_rank();
  const int parallel_size = mesh.parallel_size();

  stk::CommSparse all( mesh.parallel() );

  // Communication sizing:

  for (const EntityProc& ep : sendGhosts) {
    const int owner = mesh.parallel_owner_rank(ep.first);
    all.send_buffer(ep.second).skip<EntityKey>(1).skip<int>(1);
    if ( owner != parallel_rank ) {
      all.send_buffer(owner).skip<EntityKey>(1).skip<int>(1);
    }
  }

  all.allocate_buffers();

  // Loop thru all entities in sendGhosts, send the entity key to the sharing/ghosting proc
  // Also, if the owner of the entity is NOT me, also send the entity key to the owing proc

  // Communication packing (with message content comments):
  for (EntityProc& entProc : sendGhosts) {
    const int owner = mesh.parallel_owner_rank(entProc.first);

    // Inform receiver of ghosting, the receiver does not own
    // and does not share this entity.
    // The ghost either already exists or is a to-be-done new ghost.
    // This status will be resolved on the final communication pass
    // when new ghosts are packed and sent.

    const EntityKey entity_key = mesh.entity_key(entProc.first);
    const int proc = entProc.second;

    all.send_buffer(proc).pack(entity_key).pack(proc);

    if ( owner != parallel_rank ) {
      // I am not the owner of this entity.
      // Inform the owner of this ghosting need.
      all.send_buffer(owner).pack(entity_key).pack(proc);

      // Erase it from my processor's ghosting responsibility:
      // The iterator passed to the erase method will be invalidated.
      entProc.first = Entity();
    }
  }

  auto shouldRemove = [&](const EntityProc& ep){ return ep.first.local_offset() == 0; };
  sendGhosts.erase(std::remove_if(sendGhosts.begin(), sendGhosts.end(), shouldRemove),
                   sendGhosts.end());

  all.communicate();

  // Loop thru all the buffers, and insert ghosting request for entity e to other proc
  // if the proc sending me the data is me, then insert into new_recv.
  // Communication unpacking:
  for ( int p = 0 ; p < parallel_size ; ++p ) {
    CommBuffer & buf = all.recv_buffer(p);
    while ( buf.remaining() ) {

      EntityKey entity_key;
      int proc = 0;

      buf.unpack(entity_key).unpack(proc);

      Entity const e = mesh.get_entity( entity_key );

      if ( parallel_rank != proc ) {
        //  Receiving a ghosting need for an entity I own.
        //  Add it to my send list.
        STK_ThrowRequireMsg(mesh.is_valid(e),
            "Unknown entity key: " <<
            mesh.mesh_meta_data().entity_rank_name(entity_key.rank()) <<
            "[" << entity_key.id() << "]");
        sendGhosts.push_back(EntityProc(e, proc));
      }
      else if ( mesh.is_valid(e) ) {
        //  I am the receiver for this ghost.
        //  If I already have it add it to the receive list,
        //  otherwise don't worry about it - I will receive
        //  it in the final new-ghosting communication.
        recvGhosts.insert( mesh.entity_key(e) );
      }
    }
  }
}

struct OwnerEntityReceiver {
  OwnerEntityReceiver(int o, Entity e, int r) : owner(o), entity(e), receiver(r){}
  int owner;
  Entity entity;
  int receiver;
};

class GhostCommHelper
{
public:
  GhostCommHelper(const BulkData& mesh,
                  const std::vector<Entity>& removeRecvGhosts,
                  EntityProcVec& newSendGhosts)
  : m_mesh(mesh),
    m_commSparse(mesh.parallel()),
    m_rmRecvCounts(mesh.parallel_size(), 0),
    m_nonOwnedSendCounts(mesh.parallel_size(), 0),
    m_nonOwnedSends()
  {
    const int myRank = mesh.parallel_rank();

    for(Entity entity : removeRecvGhosts) {
      ++m_rmRecvCounts[mesh.parallel_owner_rank(entity)];
    }

    for(EntityProc& ep : newSendGhosts) {
      int owner = mesh.parallel_owner_rank(ep.first);
      if (owner != myRank) {
        ++m_nonOwnedSendCounts[owner];
        m_nonOwnedSends.emplace_back(owner,ep.first,ep.second);
        ep.first = Entity();
      }
    }
    auto shouldRemove = [&](const EntityProc& ep) { return ep.first.local_offset() == 0; };
    newSendGhosts.erase(
          std::remove_if(newSendGhosts.begin(), newSendGhosts.end(), shouldRemove),
          newSendGhosts.end());

    auto isLess = [](const OwnerEntityReceiver& lhs, const OwnerEntityReceiver& rhs)
                    { return lhs.owner < rhs.owner; };
    std::sort(m_nonOwnedSends.begin(), m_nonOwnedSends.end(), isLess);
  }

  void setup_send_buffers()
  {
    for(int p=0; p<m_commSparse.parallel_size(); ++p) {
      m_commSparse.send_buffer(p).skip<unsigned>(1);
      if (m_rmRecvCounts[p] > 0) {
        m_commSparse.send_buffer(p).skip<EntityKey>(m_rmRecvCounts[p]);
      }
      m_commSparse.send_buffer(p).skip<unsigned>(1);
      if (m_nonOwnedSendCounts[p] > 0) {
        for(unsigned i=0; i<m_nonOwnedSendCounts[p]; ++i) {
          m_commSparse.send_buffer(p).skip<EntityKey>(1);
          m_commSparse.send_buffer(p).skip<int>(1);
        }
      }
    }

    m_commSparse.allocate_buffers();
  }

  void pack_and_communicate_buffers(const std::vector<Entity>& removeRecvGhosts)
  {
    unsigned idx = 0;

    for(int p=0; p<m_commSparse.parallel_size(); ++p) {
      stk::CommBuffer& buf = m_commSparse.send_buffer(p);

      buf.pack<unsigned>(m_rmRecvCounts[p]);
      if (m_rmRecvCounts[p] > 0) {
        for(Entity entity : removeRecvGhosts) {
          int owner = m_mesh.parallel_owner_rank(entity);
          if (p == owner) {
            buf.pack<EntityKey>(m_mesh.entity_key(entity));
          }
        }
      }

      buf.pack<unsigned>(m_nonOwnedSendCounts[p]);
      if (m_nonOwnedSendCounts[p] > 0) {
        for(unsigned i=idx; i<idx+m_nonOwnedSendCounts[p]; ++i) {
          int owner = m_nonOwnedSends[i].owner;
          Entity ent = m_nonOwnedSends[i].entity;
          int proc = m_nonOwnedSends[i].receiver;
          STK_ThrowRequireMsg(owner == p, "Error, owner doesn't match");
          buf.pack<EntityKey>(m_mesh.entity_key(ent));
          buf.pack<int>(proc);
        }
        idx += m_nonOwnedSendCounts[p];
      }
    }

    m_commSparse.communicate();
  }

  void unpack_recv_buffers(EntityProcVec& newSendGhosts,
                           std::set<EntityKeyProc>& removeSendGhosts)
  {
    for ( int p = 0 ; p < m_commSparse.parallel_size() ; ++p ) {
      CommBuffer & buf = m_commSparse.recv_buffer(p);
      unsigned numRemoveSends = 0;
      buf.unpack<unsigned>(numRemoveSends);

      for(unsigned i=0; i<numRemoveSends; ++i) {
        EntityKey key;
        buf.unpack<EntityKey>(key);
        removeSendGhosts.insert(EntityKeyProc(key,p));
      }

      unsigned numReferredSends = 0;
      buf.unpack<unsigned>(numReferredSends);
      for(unsigned i=0; i<numReferredSends; ++i) {
        EntityKey key;
        buf.unpack<EntityKey>(key);
        int proc = 0;
        buf.unpack<int>(proc);
        Entity entity = m_mesh.get_entity(key);
        STK_ThrowRequireMsg(m_mesh.is_valid(entity), "Recvd invalid entity key");
        newSendGhosts.push_back(EntityProc(entity,proc));
      }
    }
  }

private:
  const BulkData& m_mesh;
  CommSparse m_commSparse;
  std::vector<unsigned> m_rmRecvCounts;
  std::vector<unsigned> m_nonOwnedSendCounts;
  std::vector<OwnerEntityReceiver> m_nonOwnedSends;
};

void comm_sync_send_recv(const BulkData & mesh ,
                         const std::vector<Entity>& removeRecvGhosts,
                         EntityProcVec& newSendGhosts,
                         std::set< EntityKeyProc> & removeSendGhosts)
{
  GhostCommHelper commHelper(mesh, removeRecvGhosts, newSendGhosts);

  commHelper.setup_send_buffers();

  commHelper.pack_and_communicate_buffers(removeRecvGhosts);

  commHelper.unpack_recv_buffers(newSendGhosts, removeSendGhosts);
}

void comm_sync_aura_send_recv(
  BulkData & mesh ,
  std::vector<EntityProc>& sendGhosts,
  EntityProcMapping& entityProcMapping,
  std::vector<bool>& ghostStatus )
{
  const int parallel_rank = mesh.parallel_rank();
  const int parallel_size = mesh.parallel_size();

  stk::CommSparse all( mesh.parallel() );

  // Communication sizing:

  for (const EntityProc& ep : sendGhosts) {
    const int owner = mesh.parallel_owner_rank(ep.first);
    all.send_buffer( ep.second ).skip<EntityKey>(1).skip<int>(1);
    if ( owner != parallel_rank ) {
      all.send_buffer( owner ).skip<EntityKey>(1).skip<int>(1);
    }
  }

  all.allocate_buffers();

  // Loop thru all entities in sendGhosts, send the entity key to the sharing/ghosting proc
  // Also, if the owner of the entity is NOT me, also send the entity key to the owing proc

  // Communication packing (with message content comments):
  for (const EntityProc& ep : sendGhosts) {
    const Entity entity = ep.first;
    const int owner = mesh.parallel_owner_rank(entity);

    // Inform receiver of ghosting, the receiver does not own
    // and does not share this entity.
    // The ghost either already exists or is a to-be-done new ghost.
    // This status will be resolved on the final communication pass
    // when new ghosts are packed and sent.

    const EntityKey entity_key = mesh.entity_key(entity);
    const int proc = ep.second;

    all.send_buffer( proc ).pack(entity_key).pack(proc);

    if (owner != parallel_rank) {
      all.send_buffer(owner).pack(entity_key).pack(proc);

      //Erase it from my processor's ghosting responsibility:
      entityProcMapping.eraseEntityProc(entity, proc);
    }
  }

  all.communicate();

  // Loop thru all the buffers, and insert ghosting request for entity e to other proc
  // if the proc sending me the data is me, then insert into new_recv.
  // Communication unpacking:
  for ( int p = 0 ; p < parallel_size ; ++p ) {
    CommBuffer & buf = all.recv_buffer(p);
    while ( buf.remaining() ) {

      EntityKey entity_key;
      int proc = 0;

      buf.unpack(entity_key).unpack(proc);

      Entity const e = mesh.get_entity( entity_key );

      if (parallel_rank != proc) {
        //Receiving a ghosting need for an entity I own, add it.
        entityProcMapping.addEntityProc(e, proc);
      }
      else if ( mesh.is_valid(e) ) {
        //  I am the receiver for this ghost.
        //  If I already have it add it to the receive list,
        //  otherwise don't worry about it - I will receive
        //  it in the final new-ghosting communication.
        ghostStatus[e.local_offset()] = true;
      }
    }
  }
}

void comm_sync_nonowned_sends(
  const BulkData & mesh ,
  std::vector<EntityProc> & nonOwnedSendGhosts,
  EntityProcMapping& entityProcMapping)
{
  const int parallel_size = mesh.parallel_size();
  const int parallel_rank = mesh.parallel_rank();
  stk::CommSparse commSparse( mesh.parallel() );

  for (const EntityProc& ep : nonOwnedSendGhosts) {
    const int owner = mesh.parallel_owner_rank(ep.first);
    if ( owner != mesh.parallel_rank() ) {
      commSparse.send_buffer( owner ).skip<EntityKey>(1).skip<int>(1);
    }
  }

  commSparse.allocate_buffers();

  for (const EntityProc& ep : nonOwnedSendGhosts) {
    const int owner = mesh.parallel_owner_rank(ep.first);
    if ( owner != parallel_rank ) {
      commSparse.send_buffer( owner ).pack<EntityKey>(mesh.entity_key(ep.first)).pack<int>(ep.second);
      entityProcMapping.eraseEntityProc(ep.first, ep.second);
    }
  }

  commSparse.communicate();

  for ( int p = 0 ; p < parallel_size ; ++p ) {
    CommBuffer & buf = commSparse.recv_buffer(p);
    while ( buf.remaining() ) {

      EntityKey entity_key;
      int proc = 0;

      buf.unpack(entity_key).unpack(proc);

      Entity const e = mesh.get_entity( entity_key );

      STK_ThrowAssert(parallel_rank != proc);
      STK_ThrowAssertMsg(mesh.is_valid(e), "comm_sync_nonowned_sends mod-cycle="<<mesh.synchronized_count()<<": P"<<parallel_rank
                                           <<" recvd "<<entity_key<<" from P"<<p<<" but valid entity not found.");

      //Receiving a ghosting need for an entity I own, add it.
      entityProcMapping.addEntityProc(e, proc);
    }
  }
}

std::pair<EntityRank,unsigned> get_highest_upward_connected_rank(const Bucket& bucket,
                                             unsigned bucketOrdinal,
                                             EntityRank entityRank,
                                             EntityRank maxRank)
{
  std::pair<EntityRank,unsigned> result(maxRank,0);
  while(result.first > entityRank &&
        (result.second = bucket.num_connectivity(bucketOrdinal, result.first)) == 0)
  {
    result.first = static_cast<EntityRank>(result.first-1);
  }
  return result;
}

void insert_upward_relations_for_owned(const BulkData& bulk_data,
                             const Entity entity,
                             const EntityRank entityRank,
                             const EntityRank maxRank,
                             const std::vector<int>& sharingProcs,
                             EntityProcMapping& send)
{
  // If related entity is higher rank, I own it, and it is not
  // already shared by proc, ghost it to the sharing processor.
  const MeshIndex& idx = bulk_data.mesh_index(entity);
  const Bucket& bucket = *idx.bucket;
  STK_ThrowAssert(bucket.owned());
  const unsigned bucketOrd = idx.bucket_ordinal;
  const std::pair<EntityRank,unsigned> rankAndNumConn = get_highest_upward_connected_rank(bucket, bucketOrd, entityRank, maxRank);
  const EntityRank upwardRank = rankAndNumConn.first;

  if (upwardRank > entityRank) {
    const int numRels = rankAndNumConn.second;
    const Entity* rels     = bucket.begin(bucketOrd, upwardRank);

    for (int r = 0; r < numRels; ++r) {
      Entity const upwardEntity = rels[r];
      if (bulk_data.is_valid(upwardEntity) && bulk_data.bucket(upwardEntity).owned()) {
        for(int sharingProc : sharingProcs) {
          if (upwardRank >= stk::topology::ELEM_RANK ||
              !bulk_data.in_shared(upwardEntity, sharingProc)) {
            send.addEntityProc(upwardEntity, sharingProc);
          }
        }
      }
    }
  }
  else {
    if (entityRank >= stk::topology::ELEM_RANK || !bulk_data.in_shared(entity)) {
      for(int sharingProc : sharingProcs) {
        send.addEntityProc(entity,sharingProc);
      }
    }
    else {
      for(int sharingProc : sharingProcs) {
        if (!bulk_data.in_shared(entity, sharingProc)) {
          send.addEntityProc(entity,sharingProc);
        }
      }
    }
  }
}

void move_unowned_entities_for_owner_to_ghost(BulkData & mesh, EntityProcVec& sendGhosts)
{
    const int myProcId = mesh.parallel_rank();
    const int parallel_size = mesh.parallel_size();

    stk::CommSparse all(mesh.parallel());

    for(const EntityProc& entProc : sendGhosts) {
        const int owner = mesh.parallel_owner_rank(entProc.first);
        if(owner != myProcId) {
            all.send_buffer(owner).skip < stk::mesh::EntityKey > (1).skip<int>(1);
        }
    }

    all.allocate_buffers();

    for(EntityProc& entProc : sendGhosts) {
        const int owner = mesh.parallel_owner_rank(entProc.first);

        const stk::mesh::EntityKey entity_key = mesh.entity_key(entProc.first);
        const int proc = entProc.second;

        if(owner != myProcId) {
            all.send_buffer(owner).pack(entity_key).pack(proc);
            entProc.first = Entity();
            entProc.second = -1;
        }
    }

    all.communicate();

    auto shouldRemove = [&](const EntityProc& ep) { return ep.first.local_offset() == 0; };

    sendGhosts.erase(std::remove_if(sendGhosts.begin(), sendGhosts.end(), shouldRemove),
                     sendGhosts.end());

    for(int p = 0; p < parallel_size; ++p)
    {
        stk::CommBuffer & buf = all.recv_buffer(p);
        while(buf.remaining())
        {
            EntityKey entity_key;
            int procToGhost = 0;

            buf.unpack(entity_key).unpack(procToGhost);

            Entity const e = mesh.get_entity(entity_key);

            STK_ThrowRequireWithSierraHelpMsg(mesh.is_valid(e) && mesh.parallel_owner_rank(e) == myProcId);

            if(myProcId != procToGhost) {
                sendGhosts.push_back(EntityProc(e, procToGhost));
            }
        }
    }
}

void convert_part_ordinals_to_parts(const stk::mesh::MetaData& meta,
                                    const OrdinalVector& input_ordinals,
                                    stk::mesh::PartVector& output_parts)
{
    output_parts.clear();
    output_parts.reserve(input_ordinals.size());
    for(unsigned ipart = 0; ipart < input_ordinals.size(); ++ipart)
    {
        output_parts.push_back(&meta.get_part(input_ordinals[ipart]));
    }
}

stk::mesh::ConnectivityOrdinal get_ordinal_from_side_entity(const std::vector<stk::mesh::Entity> &sides,
                                                            stk::mesh::ConnectivityOrdinal const * ordinals,
                                                            stk::mesh::Entity side)
{
    for(unsigned i = 0; i<sides.size(); ++i)
    {
        if(sides[i] == side)
            return ordinals[i];
    }

    return stk::mesh::INVALID_CONNECTIVITY_ORDINAL;
}

bool are_any_parts_ranked(const stk::mesh::MetaData& meta,
                          const OrdinalVector& partOrdinals)
{
  const PartVector& allParts = meta.get_parts();
  for(Ordinal ord : partOrdinals) {
    if (allParts[ord]->primary_entity_rank() != InvalidEntityRank) {
      return true;
    }
  }
  return false;
}

void filter_out( OrdinalVector & vec ,
                 const OrdinalVector & parts ,
                 OrdinalVector & removed ,
                 bool trackRemoved )
{
  OrdinalVector::iterator i , j ;
  i = j = vec.begin();

  OrdinalVector::const_iterator ip = parts.begin() ;

  while ( j != vec.end() && ip != parts.end() ) {
    if      ( *ip < *j ) { ++ip ; }
    else if ( *j < *ip ) { *i = *j ; ++i ; ++j ; }
    else {
      if (trackRemoved) {
        removed.push_back( *ip );
      }
      ++j ;
      ++ip ;
    }
  }

  if ( i != j ) { vec.erase( i , j ); }
}

void merge_in( OrdinalVector & vec , const OrdinalVector & parts )
{
  std::vector<unsigned>::iterator i = vec.begin();
  OrdinalVector::const_iterator ip = parts.begin() ;

  for ( ; i != vec.end() && ip != parts.end() ; ++i ) {

    const unsigned ord = *ip;

    if ( ord <= *i ) {
      if ( ord < *i ) { i = vec.insert( i , ord ); }
      ++ip ;
    }
  }

  for ( ; ip != parts.end() ; ++ip ) {
    const unsigned ord = *ip;
    vec.push_back( ord );
  }
}

stk::mesh::ConnectivityOrdinal get_ordinal_for_element_side_pair(const stk::mesh::BulkData &bulkData, stk::mesh::Entity element, stk::mesh::Entity side)
{
    const stk::mesh::Entity * sides = bulkData.begin(element, bulkData.mesh_meta_data().side_rank());
    stk::mesh::ConnectivityOrdinal const * ordinals = bulkData.begin_ordinals(element, bulkData.mesh_meta_data().side_rank());
    std::vector<stk::mesh::Entity> sideVector(sides, sides+bulkData.num_sides(element));
    return get_ordinal_from_side_entity(sideVector, ordinals, side);
}

void fill_inducible_parts_from_list(const MetaData& meta,
                                    const OrdinalVector & partList,
                                    EntityRank rank,
                                    OrdinalVector &induciblePartsFromList)
{
  const PartVector& allParts = meta.get_parts();
  for (size_t i = 0; i < partList.size(); i++) {
    if (allParts[partList[i]]->should_induce(rank)) {
      induciblePartsFromList.push_back(partList[i]);
    }
  }
}

void fill_part_list_differences(const BulkData &mesh,
                                Entity entity,
                                const PartVector &recv_parts,
                                std::set<std::string> &thisProcExtraParts,
                                std::set<std::string> &otherProcExtraParts)
{
  const Bucket & bucket = mesh.bucket(entity);
  std::pair<const unsigned *,const unsigned *> part_ordinals = bucket.superset_part_ordinals();
  const PartVector & mesh_parts = mesh.mesh_meta_data().get_parts();

  std::set<std::string> thisProcPartNames;
  for (const unsigned * k = part_ordinals.first; k < part_ordinals.second; ++k) {
    if (mesh_parts[*k]->entity_membership_is_parallel_consistent()) {
      if (mesh_parts[*k]->name() != "{OWNS}" && mesh_parts[*k]->name() != "{SHARES}" && mesh_parts[*k]->primary_entity_rank() != stk::mesh::InvalidEntityRank)
      {
        thisProcPartNames.insert(mesh_parts[*k]->name());
      }
    }
  }

  std::set<std::string> otherProcPartNames;
  for (const Part* part : recv_parts) {
    if (part->entity_membership_is_parallel_consistent()) {
      if (part->name() != "{OWNS}" && part->name() != "{SHARES}" && part->primary_entity_rank() != stk::mesh::InvalidEntityRank) {
        otherProcPartNames.insert(part->name());
      }
    }
  }

  std::set_difference(thisProcPartNames.begin(),
                      thisProcPartNames.end(),
                      otherProcPartNames.begin(),
                      otherProcPartNames.end(),
                      std::inserter(thisProcExtraParts, thisProcExtraParts.begin()));

  std::set_difference(otherProcPartNames.begin(),
                      otherProcPartNames.end(),
                      thisProcPartNames.begin(),
                      thisProcPartNames.end(),
                      std::inserter(otherProcExtraParts, otherProcExtraParts.begin()));
}

void check_size_of_types()
{
    const size_t sizeof_entity_in_bytes = sizeof(Entity);
#ifdef STK_32BIT_ENTITY
    STK_ThrowRequireMsg(4 == sizeof_entity_in_bytes, "sizeof(Entity) expected to be 4, is instead "<<sizeof_entity_in_bytes);
#else
    STK_ThrowRequireMsg(8 == sizeof_entity_in_bytes, "sizeof(Entity) expected to be 8, is instead "<<sizeof_entity_in_bytes);
#endif

    const size_t sizeof_connectivityordinal_in_bytes = sizeof(ConnectivityOrdinal);
#ifdef STK_16BIT_CONNECTIVITY_ORDINAL
    STK_ThrowRequireMsg(2 == sizeof_connectivityordinal_in_bytes, "sizeof(ConnectivityOrdinal) expected to be 2, is instead "<<sizeof_connectivityordinal_in_bytes);
#else
    STK_ThrowRequireMsg(4 == sizeof_connectivityordinal_in_bytes, "sizeof(ConnectivityOrdinal) expected to be 4, is instead "<<sizeof_connectivityordinal_in_bytes);
#endif
}

void require_valid_relation(const char action[],
                            const BulkData& mesh,
                            const Entity e_from,
                            const Entity e_to)
{
  const bool error_rank      = !(mesh.entity_rank(e_from) > mesh.entity_rank(e_to));
  const bool error_nil_from  = !mesh.is_valid(e_from);
  const bool error_nil_to    = !mesh.is_valid(e_to);

  if ( error_rank || error_nil_from || error_nil_to ) {
    std::ostringstream msg ;

    msg << "Could not " << action << " relation from entity "
        << mesh.entity_key(e_from) << " to entity " << mesh.entity_key(e_to) << "\n";

    STK_ThrowErrorMsgIf( error_nil_from  || error_nil_to,
                     msg.str() << ", entity was destroyed");
    STK_ThrowErrorMsgIf( error_rank, msg.str() <<
                     "A relation must be from higher to lower ranking entity");
  }
}

bool is_valid_relation(const BulkData& mesh,
                       Entity e_from,
                       Entity e_to,
                       EntityRank e_to_rank,
                       ConnectivityOrdinal ord)
{
  const MeshIndex& meshIndex = mesh.mesh_index(e_from);
  const Bucket* bPtr = meshIndex.bucket;
  const unsigned bOrd = meshIndex.bucket_ordinal;
  const unsigned num = bPtr->num_connectivity(bOrd, e_to_rank);
  if (num > 0) {
    const Entity* conn = bPtr->begin(bOrd, e_to_rank);
    const ConnectivityOrdinal* ords = bPtr->begin_ordinals(bOrd, e_to_rank);
    for(unsigned i=0; i<num; ++i) {
      if (ords[i] == ord && conn[i] == e_to) {
        return true;
      }
    }
  }
  return false;
}

bool is_good_rank_and_id(const MetaData& meta,
                         EntityRank rank,
                         EntityId id)
{
  const EntityRank rank_count = meta.entity_rank_count();
  const bool ok_id   = EntityKey::is_valid_id(id);
  const bool ok_rank = rank < rank_count &&
                 !(rank == stk::topology::FACE_RANK && meta.spatial_dimension() == 2);

  return ok_id && ok_rank;
}

EntityId get_global_max_id_in_use(const BulkData& mesh,
                                  EntityRank rank,
                                  const std::vector<Entity::entity_value_type>& deletedEntitiesCurModCycle,
                                  const std::vector<EntityId>& reservedIds)
{
  EntityId localMax = stk::mesh::get_max_id_on_local_proc(mesh, rank);
  EntityId localMaxReserved = reservedIds.empty() ? 0 : *std::max_element(reservedIds.begin(), reservedIds.end());
  localMax = std::max(localMax, localMaxReserved);

  for (Entity::entity_value_type local_offset : deletedEntitiesCurModCycle) {
    stk::mesh::Entity entity(local_offset);
    if ( mesh.is_valid(entity) && mesh.entity_rank(entity) == rank ) {
      localMax = std::max(localMax, mesh.entity_key(entity).id());
    }
  }

  return stk::get_global_max(mesh.parallel(), localMax);
}

void check_declare_element_side_inputs(const BulkData & mesh,
                                       const Entity elem,
                                       const unsigned localSideId)
{
  stk::topology elem_top = mesh.bucket(elem).topology();

  STK_ThrowErrorMsgIf( elem_top == stk::topology::INVALID_TOPOLOGY,
          "Element[" << mesh.identifier(elem) << "] has no defined topology");

  stk::topology invalid = stk::topology::INVALID_TOPOLOGY;
  stk::topology side_top =
          ((elem_top != stk::topology::INVALID_TOPOLOGY) && (localSideId < elem_top.num_sides()) )
          ? elem_top.side_topology(localSideId) : invalid;

  STK_ThrowErrorMsgIf( elem_top!=stk::topology::INVALID_TOPOLOGY && localSideId >= elem_top.num_sides(),
          "For elem " << mesh.identifier(elem) << ", localSideId " << localSideId
          << ", localSideId exceeds " << elem_top.name()
          << ".num_sides() = " << elem_top.num_sides());

  STK_ThrowErrorMsgIf( side_top == stk::topology::INVALID_TOPOLOGY,
          "For elem " << mesh.identifier(elem) << ", localSideId "
          << localSideId << ", No element topology found");
}

bool connect_edge_or_face_to_elements_impl(stk::mesh::BulkData& bulk, stk::mesh::Entity entity)
{
  stk::mesh::EntityRank entityRank = bulk.entity_rank(entity);
  if (entityRank != stk::topology::EDGE_RANK && entityRank != stk::topology::FACE_RANK) {
    return false;
  }
  const stk::mesh::Entity* nodes = bulk.begin_nodes(entity);
  unsigned numNodes = bulk.num_nodes(entity);
  stk::mesh::EntityVector elems;
  stk::mesh::impl::find_entities_these_nodes_have_in_common(bulk, stk::topology::ELEM_RANK, numNodes, nodes, elems);
  if (elems.empty()) {
    return false;
  }

  stk::mesh::EntityVector entityNodes(bulk.begin_nodes(entity), bulk.end_nodes(entity));
  stk::topology entityTopo = bulk.bucket(entity).topology();
  for(stk::mesh::Entity elem : elems) {
    stk::mesh::OrdinalAndPermutation ordinalAndPerm = stk::mesh::get_ordinal_and_permutation(bulk, elem,
                                                                                            entityRank, entityNodes);

    if(ordinalAndPerm.first == stk::mesh::INVALID_CONNECTIVITY_ORDINAL) { return false; }

    stk::mesh::impl::connect_element_to_entity(bulk, elem, entity, ordinalAndPerm.first, stk::mesh::PartVector{}, entityTopo);
  }
  return true;
}

bool connect_edge_to_elements(stk::mesh::BulkData& bulk, stk::mesh::Entity edge)
{
  return connect_edge_or_face_to_elements_impl(bulk, edge);
}

void connect_face_to_elements(stk::mesh::BulkData& bulk, stk::mesh::Entity face)
{
  STK_ThrowRequireMsg(connect_edge_or_face_to_elements_impl(bulk, face),
                  "Face with id: " << bulk.identifier(face) << " has no valid connectivity to elements");
}

bool has_upward_recv_ghost_connectivity(const stk::mesh::BulkData &bulk,
                                        const stk::mesh::Ghosting& ghosting,
                                        stk::mesh::Entity entity)
{
  if(!bulk.is_valid(entity))
    return false;

  const stk::mesh::EntityRank entityRank = bulk.entity_rank(entity);
  const stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(bulk.mesh_meta_data().entity_rank_count());
  for(stk::mesh::EntityRank conRank = static_cast<stk::mesh::EntityRank>(entityRank + 1); conRank <= endRank; ++conRank)
  {
    unsigned numConnected = bulk.num_connectivity(entity, conRank);
    if(numConnected > 0) {
      const stk::mesh::Entity* conn = bulk.begin(entity, conRank);
      for(unsigned i=0; i<numConnected; ++i) {
        if (bulk.in_receive_ghost(ghosting, conn[i])) {
          return true;
        }
      }
    }
  }

  return false;
}

bool has_upward_send_ghost_connectivity(const stk::mesh::BulkData &bulk,
                                        const stk::mesh::Ghosting& ghosting,
                                        int proc,
                                        stk::mesh::Entity entity)
{
  if(!bulk.is_valid(entity))
    return false;

  const stk::mesh::EntityRank entityRank = bulk.entity_rank(entity);
  const stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(bulk.mesh_meta_data().entity_rank_count());
  for(stk::mesh::EntityRank conRank = static_cast<stk::mesh::EntityRank>(entityRank + 1); conRank <= endRank; ++conRank)
  {
    unsigned numConnected = bulk.num_connectivity(entity, conRank);
    if(numConnected > 0) {
      const stk::mesh::Entity* conn = bulk.begin(entity, conRank);
      for(unsigned i=0; i<numConnected; ++i) {
        if (bulk.in_send_ghost(ghosting, bulk.entity_key(conn[i]), proc)) {
          return true;
        }
      }
    }
  }

  return false;
}

bool has_upward_connectivity(const stk::mesh::BulkData &bulk, stk::mesh::Entity entity)
{
  if(!bulk.is_valid(entity))
    return false;

  const stk::mesh::EntityRank entityRank = bulk.entity_rank(entity);
  const stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(bulk.mesh_meta_data().entity_rank_count());
  for(stk::mesh::EntityRank conRank = static_cast<stk::mesh::EntityRank>(entityRank + 1); conRank <= endRank; ++conRank)
  {
    unsigned numConnected = bulk.num_connectivity(entity, conRank);
    if(numConnected > 0)
      return true;
  }

  return false;
}

bool can_destroy_entity(const stk::mesh::BulkData &bulk, stk::mesh::Entity entity)
{
  return bulk.is_valid(entity) && !impl::has_upward_connectivity(bulk, entity);
}

void destroy_upward_connected_aura_entities(stk::mesh::BulkData &bulk,
                                            stk::mesh::Entity connectedEntity,
                                            EntityVector& scratchSpace)
{
  impl::StoreEntity storeEntity(bulk);
  impl::VisitUpwardClosure(bulk, connectedEntity, storeEntity);

  storeEntity.store_visited_entities_in_vec(scratchSpace);
  stk::util::sort_and_unique(scratchSpace, EntityLess(bulk));

  for(unsigned i=0; i<scratchSpace.size(); ++i) {
    int reverseIdx = scratchSpace.size() - 1 - i;
    Entity upwardEntity = scratchSpace[reverseIdx];

    if (bulk.is_valid(upwardEntity) && bulk.bucket(upwardEntity).in_aura()) {
      bulk.destroy_entity(upwardEntity);
    }
  }
}

void print_upward_connected_entities(stk::mesh::BulkData& bulk,
                                     stk::mesh::Entity entity,
                                     std::ostream& os)
{
  impl::StoreEntity storeEntity(bulk);
  impl::VisitUpwardClosure(bulk, entity, storeEntity);

  EntityVector scratchSpace;
  storeEntity.store_visited_entities_in_vec(scratchSpace);
  stk::util::sort_and_unique(scratchSpace, EntityLess(bulk));

  os << "upw-conn-entities of " << bulk.entity_key(entity) << ": ";
  for(unsigned i=0; i<scratchSpace.size(); ++i) {
    Entity upwardEntity = scratchSpace[i];
    if (upwardEntity != entity) {
      const bool owned = bulk.bucket(upwardEntity).owned();
      const bool shrd = bulk.bucket(upwardEntity).shared();
      const bool recvAura = bulk.bucket(upwardEntity).in_aura();
      const bool recvCG = bulk.in_receive_custom_ghost(bulk.entity_key(upwardEntity));
      os << bulk.entity_key(upwardEntity) << "{"<<(owned?"Owned":"")<<(shrd?"Shrd":"")<<(recvAura?"RcvAu":"")<<(recvCG?"RcvCG":"")
         << bulk.state(upwardEntity) << "} ";
    }
  }
  os << std::endl;
}

} // namespace impl
} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

