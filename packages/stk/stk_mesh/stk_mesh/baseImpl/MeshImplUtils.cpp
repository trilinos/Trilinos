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
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>

#include <vector>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {
namespace impl {

void find_entities_these_nodes_have_in_common(BulkData& mesh, stk::mesh::EntityRank rank, unsigned numNodes, const Entity* nodes, std::vector<Entity>& entity_vector)
{
  entity_vector.clear();
  std::vector<Entity> tmp;
  for(unsigned i=0; i<numNodes; ++i) {
    const Entity* entities = mesh.begin(nodes[i],rank);
    unsigned numEntities = mesh.num_connectivity(nodes[i],rank);
    tmp.assign(entities, entities+numEntities);
    std::sort(tmp.begin(), tmp.end());
    if (i==0) {
      entity_vector.assign(tmp.begin(), tmp.end());
    }
    else {
       std::vector<Entity> intersect;
       std::back_insert_iterator<std::vector<Entity> > intersect_itr(intersect);
       std::set_intersection(entity_vector.begin(), entity_vector.end(),
                             tmp.begin(), tmp.end(),
                             intersect_itr);
       entity_vector.swap(intersect);
    }
  }
}

void find_elements_these_nodes_have_in_common(BulkData& mesh, unsigned numNodes, const Entity* nodes, std::vector<Entity>& entity_vector)
{
    find_entities_these_nodes_have_in_common(mesh,stk::topology::ELEMENT_RANK,numNodes,nodes,entity_vector);
}

void find_faces_these_nodes_have_in_common(BulkData& mesh, unsigned numNodes, const Entity* nodes, std::vector<Entity>& entity_vector)
{
    find_entities_these_nodes_have_in_common(mesh,stk::topology::FACE_RANK,numNodes,nodes,entity_vector);
}

bool do_these_nodes_have_any_shell_elements_in_common(BulkData& mesh, unsigned numNodes, const Entity* nodes)
{
  std::vector<Entity> elems;
  find_elements_these_nodes_have_in_common(mesh, numNodes, nodes, elems);
  bool found_shell = false;
  for (unsigned count = 0; count < elems.size(); ++count) {
      if (mesh.bucket(elems[count]).topology().is_shell()) {
          found_shell = true;
      }
  }
  return found_shell;
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

void find_face_nodes_for_side(BulkData& mesh, Entity element, int side_ordinal, EntityVector & permuted_face_nodes)
{
    stk::topology elemTopology = mesh.bucket(element).topology();
    stk::topology faceTopology = elemTopology.face_topology(side_ordinal);
    const size_t num_elem_nodes = elemTopology.num_nodes();

    std::vector<EntityId> elem_node_ids(num_elem_nodes);
    Entity const *elem_nodes = mesh.begin_nodes(element);
    ThrowRequire(mesh.num_nodes(element) == num_elem_nodes);
    for (size_t n=0; n<num_elem_nodes; ++n) {
        elem_node_ids[n] = mesh.identifier(elem_nodes[n]);
    }

    // Use node identifier instead of node local_offset for cross-processor consistency.
    typedef std::vector<EntityId>  EntityIdVector;
    EntityIdVector side_node_ids(faceTopology.num_nodes());
    elemTopology.face_nodes(elem_node_ids, side_ordinal, side_node_ids.begin());
    unsigned smallest_permutation;
    permuted_face_nodes.resize(faceTopology.num_nodes());
    //if this is a shell OR these nodes are connected to a shell
    EntityVector side_nodes(faceTopology.num_nodes());
    for (unsigned count=0 ; count<faceTopology.num_nodes() ; ++count) {
        side_nodes[count] = mesh.get_entity(stk::topology::NODE_RANK,side_node_ids[count]);
    }
    bool is_connected_to_shell = stk::mesh::impl::do_these_nodes_have_any_shell_elements_in_common(mesh,faceTopology.num_nodes(),&side_nodes[0]);

    if (elemTopology.is_shell() || is_connected_to_shell) {

        EntityIdVector element_node_id_vector(faceTopology.num_nodes());
        EntityIdVector element_node_ordinal_vector(faceTopology.num_nodes());
        EntityVector element_node_vector(faceTopology.num_nodes());
        elemTopology.face_node_ordinals(side_ordinal, &element_node_ordinal_vector[0]);
        for (unsigned count = 0; count < faceTopology.num_nodes(); ++count) {
            element_node_vector[count] = mesh.begin_nodes(element)[element_node_ordinal_vector[count]];
            element_node_id_vector[count] = mesh.identifier(element_node_vector[count]);
        }
        smallest_permutation = faceTopology.lexicographical_smallest_permutation_preserve_polarity(side_node_ids, element_node_id_vector);
        faceTopology.permutation_nodes(&element_node_vector[0], smallest_permutation, permuted_face_nodes.begin());
    }
    else {
        smallest_permutation = faceTopology.lexicographical_smallest_permutation(side_node_ids);
        EntityVector face_nodes(faceTopology.num_nodes());
        elemTopology.face_nodes(elem_nodes, side_ordinal, face_nodes.begin());
        faceTopology.permutation_nodes(face_nodes, smallest_permutation, permuted_face_nodes.begin());
    }
}


stk::mesh::Entity get_or_create_face_at_element_side(stk::mesh::BulkData & bulk,
                                                     stk::mesh::Entity elem,
                                                     int side_ordinal,
                                                     int new_face_global_id,
                                                     stk::mesh::Part & part)
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
        new_face = stk::mesh::declare_element_side(bulk, new_face_global_id, elem, side_ordinal, &part);
    } else {
        stk::mesh::PartVector add_parts(1, &part);
        bulk.change_entity_parts(new_face, add_parts );
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
    elem_topology.face_nodes(bulk.begin_nodes(elem_with_face),elem_with_face_side_ordinal,&side_nodes[0]);
    std::vector<stk::mesh::Entity> common_elements;
    stk::mesh::impl::find_elements_these_nodes_have_in_common(bulk,num_side_nodes,&side_nodes[0],common_elements);

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
                other_elem_topology.face_nodes(bulk.begin_nodes(other_elem),other_elem_side,&other_elem_side_nodes[0]);
                if (should_face_be_connected_to_element_side(side_nodes,other_elem_side_nodes,other_elem_side_topology,element_shell_status[count])) {
                    stk::mesh::declare_element_side(bulk, other_elem, face, other_elem_side);
                    break;
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

