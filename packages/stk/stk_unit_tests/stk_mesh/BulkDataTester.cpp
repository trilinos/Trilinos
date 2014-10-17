#include "BulkDataTester.hpp"
#include <stk_mesh/base/BulkData.hpp>
#include <stddef.h>                     // for size_t, NULL
#include <string.h>                     // for memcpy, strcmp
#include <algorithm>                    // for sort, lower_bound, unique, etc
#include <boost/foreach.hpp>            // for auto_any_base, etc
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <sstream>
#include <fstream>
#include <iterator>                     // for back_insert_iterator, etc
#include <set>                          // for set, set<>::iterator, etc
#include <stk_mesh/base/Bucket.hpp>     // for Bucket, BucketIdComparator, etc
#include <stk_mesh/base/GetEntities.hpp>  // for get_selected_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, print_entity_key, etc
#include <stk_mesh/baseImpl/EntityRepository.hpp>  // for EntityRepository, etc
#include <stk_mesh/baseImpl/Partition.hpp>  // for Partition
#include <stk_util/environment/ReportHandler.hpp>  // for ThrowRequireMsg, etc
#include <stk_util/parallel/ParallelComm.hpp>  // for CommBuffer, CommAll, etc
#include <stk_util/parallel/ParallelReduce.hpp>  // for Reduce, all_reduce, etc
#include <stk_util/util/StaticAssert.hpp>  // for StaticAssert, etc
#include <stk_util/util/string_case_compare.hpp>
#include <string>                       // for char_traits, string, etc
#include <utility>                      // for pair, make_pair, swap
#include <vector>                       // for vector, etc
#include "boost/mpl/bool.hpp"           // for bool_
#include "boost/mpl/bool_fwd.hpp"       // for false_
#include "boost/unordered/detail/buckets.hpp"  // for iterator, etc
#include "stk_mesh/base/ConnectivityMap.hpp"  // for ConnectivityMap
#include "stk_mesh/base/Entity.hpp"     // for Entity, operator<<, etc
#include "stk_mesh/base/EntityCommDatabase.hpp"  // for pack_entity_info, etc
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey, etc
#include "stk_mesh/base/FieldBase.hpp"  // for FieldBase, FieldMetaData, etc
#include "stk_mesh/base/FieldDataManager.hpp"  // for FieldDataManager, etc
#include "stk_mesh/base/Ghosting.hpp"   // for Ghosting
#include "stk_mesh/base/Part.hpp"       // for Part, remove, etc
#include "stk_mesh/base/Relation.hpp"   // for Relation, etc
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_mesh/base/Trace.hpp"      // for DiagIfWatching, Trace_, etc
#include "stk_mesh/base/Types.hpp"      // for EntityProc, EntityRank, etc
#include "stk_mesh/baseImpl/BucketRepository.hpp"  // for BucketRepository
#include "stk_mesh/baseImpl/FieldRepository.hpp"  // for FieldVector
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/parallel/DistributedIndex.hpp"  // for DistributedIndex, etc
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine, etc
#include "stk_util/util/NamedPair.hpp"
#include "stk_util/util/PairIter.hpp"   // for PairIter
#include "stk_util/util/SameType.hpp"   // for SameType, etc
#include "stk_util/util/TrackingAllocator.hpp"  // for tracking_allocator

using namespace stk::mesh;
using namespace stk;

void BulkDataTester::change_entity_owner_exp( const std::vector<EntityProc> & arg_change,
                                    bool regenerate_aura,
                                    modification_optimization mod_optimization )
{
    INCREMENT_MODIFICATION_COUNTER(PUBLIC, CHANGE_ENTITY_OWNER);
    internal_change_entity_owner_exp(arg_change, regenerate_aura, mod_optimization);
}

void BulkDataTester::internal_change_entity_owner_exp( const std::vector<EntityProc> & arg_change,
                                             bool regenerate_aura,
                                             modification_optimization mod_optimization )
{
  INCREMENT_MODIFICATION_COUNTER(INTERNAL, CHANGE_ENTITY_OWNER);

  Trace_("stk::mesh::BulkData::change_entity_owner");
  DiagIf(LOG_ENTITY, "arg_change: " << arg_change);

  const bool modStatus = modification_begin("change_entity_owner");
  ThrowRequireMsg(modStatus, "BulkData::change_entity_owner() must not be called from within a modification cycle.");

  require_ok_to_modify();

  const MetaData  & meta = m_mesh_meta_data ;
  const int       p_rank = m_parallel_rank ;
  const int       p_size = m_parallel_size ;
  ParallelMachine p_comm = m_parallel_machine ;

  //------------------------------
  // Verify the input changes, generate a clean local change list, and
  // generate the remote change list so that all processes know about
  // pending changes.

  std::vector<EntityProc> local_change( arg_change );

  // Parallel synchronous clean up and verify the requested changes:
  impl::internal_clean_and_verify_parallel_change( *this , local_change );

  //----------------------------------------
  // Parallel synchronous determination of changing shared and ghosted.

  // The two vectors below will contain changes to ghosted and shared
  // entities on this process coming from change-entity-owner requests
  // on other processes.
  std::vector<EntityProc> ghosted_change ;
  std::vector<EntityProc> shared_change ;

  impl::internal_generate_parallel_change_lists( *this , local_change ,
                            shared_change , ghosted_change );

  //------------------------------
  // Have enough information to delete all effected ghosts.
  // If the closure of a ghost contains a changing entity
  // then that ghost must be deleted.
  // Request that all ghost entities in the closure of the ghost be deleted.

  StoreEntityProcInSet store_entity_proc_in_set(*this);

  // Compute the closure of all the locally changing entities
  for ( std::vector<EntityProc>::iterator
        i = local_change.begin() ; i != local_change.end() ; ++i ) {
      store_entity_proc_in_set.target = i->second;
      impl::VisitClosureGeneral(*this,i->first,store_entity_proc_in_set,store_entity_proc_in_set);
  }
  std::set<EntityProc,EntityLess> & send_closure = store_entity_proc_in_set.entity_proc_set;


  // Calculate all the ghosts that are impacted by the set of ownership
  // changes. We look at ghosted, shared, and local changes looking for ghosts
  // that are either in the closure of the changing entity, or have the
  // changing entity in their closure. All modified ghosts will be removed.
  {
      impl::OnlyVisitGhostsOnce only_visit_ghosts_once(*this);
      StoreEntityKeyInSet store_entity_key(*this);
      for ( std::vector<EntityProc>::const_iterator i = ghosted_change.begin() ; i != ghosted_change.end() ; ++i) {
          impl::VisitAuraClosureGeneral(*this,i->first,store_entity_key,only_visit_ghosts_once);
      }
      for ( std::vector<EntityProc>::const_iterator i = shared_change.begin() ; i != shared_change.end() ; ++i) {
          impl::VisitAuraClosureGeneral(*this,i->first,store_entity_key,only_visit_ghosts_once);
      }
      for ( std::set<EntityProc,EntityLess>::const_iterator i = send_closure.begin() ; i != send_closure.end() ; ++i) {
          impl::VisitAuraClosureGeneral(*this,i->first,store_entity_key,only_visit_ghosts_once);
      }

    std::set<EntityKey> & modified_ghosts = store_entity_key.entity_key_set;

    // The ghosted change list will become invalid
    ghosted_change.clear();

    std::vector<EntityProc> empty_add ;
    std::vector<EntityKey>  remove_modified_ghosts( modified_ghosts.begin() ,
                                                    modified_ghosts.end() );

    // Skip 'm_ghosting[0]' which is the shared subset.
    for ( std::vector<Ghosting*>::iterator
          ig = m_ghosting.begin() + 1; ig != m_ghosting.end() ; ++ig ) {
      // parallel synchronous:
      internal_change_ghosting( **ig , empty_add , remove_modified_ghosts );
    }
  }

  //------------------------------
  // Consistently change the owner on all processes.
  // 1) The local_change list is giving away ownership.
  // 2) The shared_change may or may not be receiving ownership

  {
    PartVector owned;
    owned.push_back(& meta.locally_owned_part());

    for ( std::vector<EntityProc>::iterator
          i = local_change.begin() ; i != local_change.end() ; ++i ) {
      // Giving ownership, change the parts first and then
      // the owner rank to pass the ownership test.
      Entity entity = i->first;

      internal_verify_and_change_entity_parts( entity , PartVector() , owned );

      const bool changed = this->only_call_from_fmwk_set_parallel_owner_rank( entity, i->second );
      if (changed) {
        internal_change_owner_in_comm_data(entity_key(entity), i->second);
      }
    }

    for ( std::vector<EntityProc>::iterator
          i = shared_change.begin() ; i != shared_change.end() ; ++i ) {
      Entity entity = i->first;
      const bool changed = this->only_call_from_fmwk_set_parallel_owner_rank( entity, i->second );
      if (changed) {
        internal_change_owner_in_comm_data(entity_key(entity), i->second);
      }
      if ( p_rank == i->second ) { // I receive ownership
          internal_verify_and_change_entity_parts( entity , owned , PartVector() );
      }
    }
  }


  //------------------------------
  // Send entities, along with their closure, to the new owner processes
  {
    std::ostringstream error_msg ;
    int error_count = 0 ;

    CommAll comm( p_comm );

    EntityVector unique_list_of_send_closure;
    unique_list_of_send_closure.reserve(send_closure.size());

    for ( std::set<EntityProc,EntityLess>::iterator
          i = send_closure.begin() ; i != send_closure.end() ; ++i ) {
      CommBuffer & buffer = comm.send_buffer( i->second );
      Entity entity = i->first;
      pack_entity_info(*this, buffer , entity );
      pack_field_values(*this, buffer , entity );

      if (unique_list_of_send_closure.empty() || entity_key(unique_list_of_send_closure.back()) != entity_key(entity)) {
        unique_list_of_send_closure.push_back(entity);
      }
    }

    BABBLE_STK_PARALLEL_COMM(p_comm, "          change_entity_owner calling allocate_buffers");
    comm.allocate_buffers( p_size / 4 );

    for ( std::set<EntityProc,EntityLess>::iterator
          i = send_closure.begin() ; i != send_closure.end() ; ++i ) {
      CommBuffer & buffer = comm.send_buffer( i->second );
      Entity entity = i->first;
      pack_entity_info(*this, buffer , entity );
      pack_field_values(*this, buffer , entity );
    }

    BABBLE_STK_PARALLEL_COMM(p_comm, "          change_entity_owner calling communicate");
    comm.communicate();

    for ( int p = 0 ; p < p_size ; ++p ) {
      CommBuffer & buf = comm.recv_buffer(p);
      while ( buf.remaining() ) {
        PartVector parts ;
        std::vector<Relation> relations ;
        EntityKey key ;
        int owner = ~0u ;

        unpack_entity_info( buf, *this, key, owner, parts, relations );

        // Received entity information will be correct,
        // modulo the owned and shared parts

        remove( parts , meta.globally_shared_part() );

        if ( owner == p_rank ) {
          // Must have the locally_owned_part
          insert( parts , meta.locally_owned_part() );
        }
        else {
          // Must not have the locally_owned_part
          remove( parts , meta.locally_owned_part() );
        }

        std::pair<Entity ,bool> result = internal_create_entity( key );

        Entity entity = result.first;

        // The entity was copied and not created.

        internal_change_entity_parts( entity , parts , PartVector() );

        log_created_parallel_copy( entity );

        const bool changed = this->only_call_from_fmwk_set_parallel_owner_rank( entity, owner );
        if (changed) {
          internal_change_owner_in_comm_data(entity_key(entity), owner);
        }

        internal_declare_relation( entity , relations );

        if ( ! unpack_field_values(*this, buf , entity , error_msg ) ) {
          ++error_count ;
        }
      }
    }

    all_reduce( p_comm , ReduceSum<1>( & error_count ) );
    ThrowErrorMsgIf( error_count, error_msg.str() );

    // Any entity that I sent and is not in an owned closure is deleted.
    // The owned closure will be effected by received entities, so can
    // only clean up after the newly owned entities have been received.
    // Destroy backwards so as not to invalidate closures in the process.

    {
      for ( EntityVector::reverse_iterator i = unique_list_of_send_closure.rbegin() ;
            i != unique_list_of_send_closure.rend() ;
            ++i) {
        if ( ! stk::mesh::impl::member_of_owned_closure(*this, *i , p_rank ) ) {
          ThrowRequireMsg( internal_destroy_entity( *i ),
                           "Failed to destroy entity " << identifier(*i) );
        }
      }
    }

    send_closure.clear(); // Has been invalidated
  }
}

void BulkDataTester::internal_update_distributed_index_exp(
        std::vector<Entity> & shared_new )
{
    Trace_("stk::mesh::BulkData::internal_update_distributed_index");
    BABBLE_STK_PARALLEL_COMM(m_parallel_machine, "      entered internal_update_distributed_index");

    // Scoop up a list of all nodes that have had their sharing information
    // provided directly by the end-user.  These need special handling so
    // that DistributedIndex isn't used to query sharing information.
    //
    std::vector<EntityKey> shared_nodes;
    impl::gather_shared_nodes(*this, shared_nodes);

    std::vector<shared_entity_type> shared_edges;
    impl::markEntitiesForResolvingSharingInfoUsingNodes(*this, stk::topology::EDGE_RANK, shared_edges);

    std::vector<shared_entity_type> shared_faces;
    impl::markEntitiesForResolvingSharingInfoUsingNodes(*this, stk::topology::FACE_RANK, shared_faces);

    // Generate a list of all entities that have been created or modified
    // locally, who must have their sharing information in the comm lists
    // updated by making queries to DistributedIndex.  This list *does not*
    // include nodes that have had their sharing information explicitly
    // provided or edges and faces that will have their sharing information
    // resolved through node sharing.
    //
    parallel::DistributedIndex::KeyTypeVector local_created_or_modified;
    fillLocallyCreatedOrModifiedEntities(local_created_or_modified);

    // Update distributed index's ownership and sharing information with
    // locally-owned/shared entities that have been created or modified.
    // We need to manually insert the explicitly-shared nodes into this list
    // so that DistributedIndex is notified about IDs that are now in use
    // (since it manages available IDs).  These nodes will not have their
    // sharing info updated through queries to DistributedIndex.
    //
    parallel::DistributedIndex::KeyTypeVector keys_requiring_update(local_created_or_modified);
    keys_requiring_update.reserve(keys_requiring_update.size() + shared_nodes.size() + shared_edges.size() + shared_faces.size() );
    keys_requiring_update.insert(keys_requiring_update.end(), shared_nodes.begin(), shared_nodes.end());
    for (size_t i = 0; i < shared_edges.size(); ++i)
    {
      keys_requiring_update.push_back(shared_edges[i].global_key);
    }
    for (size_t i = 0; i < shared_faces.size(); ++i)
    {
      keys_requiring_update.push_back(shared_faces[i].global_key);
    }

    parallel::DistributedIndex::KeyTypeVector::const_iterator begin = keys_requiring_update.begin();
    parallel::DistributedIndex::KeyTypeVector::const_iterator end = keys_requiring_update.end();
    m_entities_index.update_keys( begin, end );

    // Update the sharing and ownership information in the comm maps for
    // locally-created or modified entities through queries to DistributedIndex.
    //

    std::vector<EntityKey> entity_keys;  // List of shared and modified entities

    update_shared_entities_global_ids( shared_edges );
    update_shared_entities_global_ids( shared_faces );

    // Add the explicitly-shared nodes to our list of shared and modified entities.
    //
    entity_keys.reserve(entity_keys.size() + shared_nodes.size());
    entity_keys.insert(entity_keys.end(), shared_nodes.begin(), shared_nodes.end());

    // Add the edges that had their sharing information updated through
    // node sharing to our list of shared and modified entities.
    //
    for (size_t i=0; i<shared_edges.size(); ++i)
    {
        Entity entity = get_entity(shared_edges[i].global_key);
        if ( is_entity_marked(entity) == BulkData::IS_SHARED )
        {
            entity_keys.push_back(shared_edges[i].global_key);
        }
    }

    for (size_t i=0; i<shared_faces.size(); ++i)
    {
        Entity entity = get_entity(shared_faces[i].global_key);
        if ( is_entity_marked(entity) == BulkData::IS_SHARED )
        {
            entity_keys.push_back(shared_faces[i].global_key);
        }
    }

    // Reset our marking array for all entities now that all sharing information
    // has been properly resolved.
    //
    std::fill(m_mark_entity.begin(), m_mark_entity.end(), static_cast<int>(BulkData::NOT_MARKED));

    // Return a list of all shared and modified entities, sorted by their entity key.
    //
    std::sort(entity_keys.begin(), entity_keys.end());
    shared_new.clear();
    shared_new.resize(entity_keys.size());
    for (size_t i=0; i<entity_keys.size(); ++i)
    {
        shared_new[i] = get_entity(entity_keys[i]);
    }
}

void BulkDataTester::internal_resolve_parallel_create_exp()
{
  Trace_("stk::mesh::BulkData::internal_resolve_parallel_create");
  BABBLE_STK_PARALLEL_COMM(m_parallel_machine, "  entered internal_resolve_parallel_create");

  ThrowRequireMsg(parallel_size() > 1, "Do not call this in serial");
  std::vector<Entity> shared_modified ;

  // Update the parallel index and
  // output shared and modified entities.
  internal_update_distributed_index_exp(shared_modified );

  // ------------------------------------------------------------
  // Claim ownership on all shared_modified entities that I own
  // and which were not created in this modification cycle. All
  // sharing procs will need to be informed of this claim.

  resolve_ownership_of_modified_entities( shared_modified );

  // ------------------------------------------------------------
  // Update shared created entities.
  // - Revise ownership to selected processor
  // - Update sharing.
  // - Work backward so the 'in_owned_closure' function
  //   can evaluate related higher ranking entities.

  move_entities_to_proper_part_ownership( shared_modified );

  my_update_comm_list( shared_modified );
}

bool BulkDataTester::internal_modification_end_for_change_entity_owner_exp( bool regenerate_aura, modification_optimization opt )
{
  Trace_("stk::mesh::BulkData::internal_modification_end");

  // The two states are MODIFIABLE and SYNCHRONiZED
  if ( m_sync_state == SYNCHRONIZED ) { return false ; }

  ThrowAssertMsg(impl::check_for_connected_nodes(*this)==0, "BulkData::modification_end ERROR, all entities with rank higher than node are required to have connected nodes.");

  ThrowAssertMsg(add_fmwk_data() || impl::check_no_shared_elements_or_higher(*this)==0, "BulkData::modification_end ERROR, Sharing of entities with rank ELEMENT_RANK or higher is not allowed.");

  if (parallel_size() > 1) {
    // Resolve modification or deletion of shared entities
    // which can cause deletion of ghost entities.

    //passes np1 segfaults (np2-np4) if not called
    internal_resolve_shared_modify_delete();

    // Resolve modification or deletion of ghost entities
    // by destroying ghost entities that have been touched.

    //doesn't fail if not called (np1-np4)
    internal_resolve_ghosted_modify_delete();

    //passes np1 segfaults np2-np4 if not called
    update_comm_list_based_on_changes_in_comm_map();

    // Resolve creation of entities: discover sharing and set unique ownership.
    //passes np1 fails on 10 tests np2, 3 tests np3, segfaults np4
    internal_resolve_parallel_create_exp();

    // Resolve part membership for shared entities.
    // This occurs after resolving creation so created and shared
    // entities are resolved along with previously existing shared entities.

    //doesn't fail if not called (np1-np4)
    internal_resolve_shared_membership();

    // Regenerate the ghosting aura around all shared mesh entities.
    if ( regenerate_aura )
    {
      //passes on np1 hangs on np2-np4 if not called
      internal_regenerate_aura();
    }

    // ------------------------------
    // Verify parallel consistency of mesh entities.
    // Unique ownership, communication lists, sharing part membership,
    // application part membership consistency.
#ifndef NDEBUG
    std::ostringstream msg ;
    bool is_consistent = true;
    is_consistent = impl::comm_mesh_verify_parallel_consistency( *this , msg );
    ThrowErrorMsgIf( !is_consistent, msg.str() );
#endif
  }
  else {
      std::vector<Entity> shared_modified ;
      //doesn't fail if not called (np1-np4)
      internal_update_distributed_index( shared_modified );
  }

  // ------------------------------
  // Now sort the bucket entities.
  // This does not change the entities, relations, or field data.
  // However, it insures that the ordering of entities and buckets
  // is independent of the order in which a set of changes were
  // performed.
  //
  //optimize_buckets combines multiple buckets in a bucket-family into
  //a single larger bucket, and also does a sort.
  //If optimize_buckets has not been requested, still do the sort.

  if ( opt == MOD_END_COMPRESS_AND_SORT ) {
    //doesn't fail if not called (np1-np4)
    m_bucket_repository.optimize_buckets();
  }
  else {
    //doesn't fail if not called (np1-np4)
    m_bucket_repository.internal_sort_bucket_entities();
  }

  // ------------------------------

  //doesn't fail if not called (np1-np4)
  m_bucket_repository.internal_modification_end();

  //doesn't fail if not called (np1-np4)
  internal_update_fast_comm_maps();

  m_sync_state = SYNCHRONIZED ;
  m_add_node_sharing_called = false;

  //doesn't fail if not called (np1-np4)
  update_deleted_entities_container();

  return true ;
}

