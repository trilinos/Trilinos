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
#include "stk_mesh/base/Types.hpp"      // for EntityProc, EntityRank, etc
#include <stk_mesh/base/HostMeshManager.hpp>
#include <stk_mesh/base/DeviceMeshManager.hpp>
#include "stk_mesh/baseImpl/BucketRepository.hpp"  // for BucketRepository
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_mesh/baseImpl/MeshCommImplUtils.hpp"
#include "stk_mesh/baseImpl/PartVectorUtils.hpp"
#include "stk_mesh/baseImpl/MeshPrintUtils.hpp"
#include "stk_mesh/baseImpl/MeshModification.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/diag/StringUtil.hpp"
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine, etc
#include "stk_util/util/NamedPair.hpp"
#include "stk_util/util/PairIter.hpp"   // for PairIter
#include "stk_util/util/SameType.hpp"   // for SameType, etc
#include "stk_util/util/SortAndUnique.hpp"
#include <algorithm>                    // for sort, lower_bound, unique, etc
#include <fstream>
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <iterator>                     // for back_insert_iterator, etc
#include <set>                          // for set, set<>::iterator, etc
#include <sstream>
#include <stddef.h>                     // for size_t, NULL
#include <stk_mesh/base/Bucket.hpp>     // for Bucket, BucketIdComparator, etc
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FaceCreator.hpp>
#include <stk_mesh/base/GetEntities.hpp>  // for get_selected_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/SideSetEntry.hpp>
#include <stk_mesh/baseImpl/ElementTopologyDeletions.hpp>
#include <stk_mesh/baseImpl/EntityRepository.hpp>  // for EntityRepository, etc
#include <stk_mesh/baseImpl/check_comm_list.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraphUpdater.hpp>
#include <stk_mesh/baseImpl/elementGraph/SideConnector.hpp>   // for SideConnector
#include <stk_mesh/baseImpl/elementGraph/SideSharingUsingGraph.hpp>
#include <stk_util/parallel/CommSparse.hpp>  // for CommSparse
#include <stk_util/parallel/GenerateParallelUniqueIDs.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>  // for Reduce, all_reduce, etc
#include <stk_util/util/ReportHandler.hpp>  // for ThrowRequireMsg, etc
#include <stk_util/util/StaticAssert.hpp>  // for StaticAssert, etc
#include <stk_util/util/string_case_compare.hpp>
#include <string.h>                     // for memcpy, strcmp
#include <string>                       // for char_traits, string, etc
#include <utility>                      // for pair, make_pair, swap
#include <vector>                       // for vector, etc

namespace stk {
namespace mesh {

namespace impl {
int Counter::counter = 0;
}

// Static constant on BulkData:
const uint16_t BulkData::orphaned_node_marking = 25000;

///////////////////////////////////////////// Functions for creating entities

void BulkData::fillEntityCommInfoForEntity(stk::mesh::Ghosting &ghost_id, stk::mesh::BulkData &mesh, std::vector<stk::mesh::EntityKey> nodes, EntityCommInfoVector &sharing_processors)
{
  size_t num_nodes = nodes.size();
  PairIterEntityComm initial_shared = mesh.internal_entity_comm_map(nodes[0],ghost_id);
  sharing_processors.assign(initial_shared.first, initial_shared.second);

  for(size_t i = 1; i < num_nodes; ++i)
  {
    PairIterEntityComm tmp_shared  = mesh.internal_entity_comm_map(nodes[i], ghost_id);
    EntityCommInfoVector new_shared_vec;

    std::set_intersection( sharing_processors.begin(), sharing_processors.end(),
                           tmp_shared.first, tmp_shared.second,
                           std::back_inserter(new_shared_vec) );
    sharing_processors = new_shared_vec;
    if (sharing_processors.empty())
    {
      break;
    }
  }
}

void BulkData::fillVectorOfSharedEntitiesByProcessor(std::vector<shared_entity_type> & potentially_shared_sides,
        std::vector<std::vector<shared_entity_type> > &shared_entities_by_proc )
{
    for(std::vector<shared_entity_type>::const_iterator itr = potentially_shared_sides.begin(),
            end = potentially_shared_sides.end(); itr != end; ++itr)
    {
        std::vector<int> sharing_processors;
        shared_procs_intersection(itr->nodes, sharing_processors);

        for(int proc_id : sharing_processors)
        {
            if(proc_id != parallel_rank())
                shared_entities_by_proc[proc_id].push_back(*itr);
        }
    }
}

void BulkData::update_owner_global_key_and_sharing_proc(stk::mesh::EntityKey global_key_other_proc, shared_entity_type& shared_entity_this_proc, int proc_id) const
{
    shared_entity_this_proc.sharing_procs.push_back(proc_id);
    if(proc_id < this->parallel_rank())
        shared_entity_this_proc.global_key = global_key_other_proc;
}

void BulkData::update_shared_entity_this_proc(EntityKey global_key_other_proc, shared_entity_type& shared_entity_this_proc, int proc_id)
{
    update_owner_global_key_and_sharing_proc(global_key_other_proc, shared_entity_this_proc, proc_id);
    this->internal_mark_entity(shared_entity_this_proc.entity, BulkData::IS_SHARED);
}

void BulkData::check_if_entity_from_other_proc_exists_on_this_proc_and_update_info_if_shared(std::vector<shared_entity_type>& shared_entities_this_proc, int proc_id, const shared_entity_type &shared_entity_from_other_proc)
{
    std::vector<shared_entity_type>::iterator shared_itr = std::lower_bound(shared_entities_this_proc.begin(), shared_entities_this_proc.end(), shared_entity_from_other_proc);
    size_t shared_itr_index = shared_itr - shared_entities_this_proc.begin();
    bool entitiesAreTheSame = impl::is_received_entity_in_local_shared_entity_list(this->use_entity_ids_for_resolving_sharing(),
           shared_itr, shared_entities_this_proc, shared_entity_from_other_proc);

    if( entitiesAreTheSame )
    {
        update_shared_entity_this_proc(shared_entity_from_other_proc.global_key, shared_entities_this_proc[shared_itr_index], proc_id);
    }
}

void removeEntitiesNotSelected(stk::mesh::BulkData &mesh, stk::mesh::Selector selected, stk::mesh::EntityVector &entities)
{
    if(selected != stk::mesh::Selector(mesh.mesh_meta_data().universal_part()))
    {
        stk::mesh::EntityVector filteredEntities;
        filteredEntities.reserve(entities.size());
        for(size_t i=0; i<entities.size(); i++)
        {
            if(selected(mesh.bucket(entities[i])))
            {
                filteredEntities.push_back(entities[i]);
            }
        }
        entities.swap(filteredEntities);
    }
}

void BulkData::mark_shared_sides_and_fill_list_of_sides_not_on_boundary(std::vector<shared_entity_type>& shared_entity_map, int proc_id, shared_entity_type &sentity,
        std::vector<stk::mesh::EntityKeyProc> &entities_to_send_data, const stk::mesh::Selector *only_consider_second_element_from_this_selector)
{
    std::vector<shared_entity_type>::iterator shared_itr = std::lower_bound(shared_entity_map.begin(), shared_entity_map.end(), sentity);
    bool entitiesAreTheSame = impl::is_received_entity_in_local_shared_entity_list(this->use_entity_ids_for_resolving_sharing(),
           shared_itr, shared_entity_map, sentity);

    if( entitiesAreTheSame )
    {
        Entity entity = this->get_entity(shared_itr->local_key);
        this->internal_mark_entity(entity, BulkData::IS_SHARED);
    }
    else
    {
        stk::mesh::EntityVector common_elements;
        stk::mesh::EntityVector nodes(sentity.nodes.size());
        for (size_t i=0;i<sentity.nodes.size();++i)
        {
            nodes[i] = this->get_entity(sentity.nodes[i]);
        }

        stk::mesh::impl::find_locally_owned_elements_these_nodes_have_in_common(*this, nodes.size(), nodes.data(), common_elements);

        if (only_consider_second_element_from_this_selector != NULL)
        {
            removeEntitiesNotSelected(*this, *only_consider_second_element_from_this_selector, common_elements);
        }

        if ( common_elements.size() > 0 )
        {
            entities_to_send_data.emplace_back(sentity.global_key, proc_id);
        }
    }
}

void BulkData::unpack_shared_entities(stk::CommSparse &comm, std::vector< std::pair<int, shared_entity_type> > &shared_entities_and_proc)
{
    for(int ip = this->parallel_size() - 1; ip >= 0; --ip)
    {
        if(ip != this->parallel_rank())
        {
            CommBuffer & buf = comm.recv_buffer(ip);
            while(buf.remaining())
            {
                shared_entity_type sentity(stk::mesh::EntityKey(), stk::mesh::Entity(), stk::topology::INVALID_TOPOLOGY);

                buf.unpack<stk::topology::topology_t>(sentity.topology);
                stk::topology entity_topology(sentity.topology);
                size_t num_nodes_on_entity = entity_topology.num_nodes();
                sentity.nodes.resize(num_nodes_on_entity);
                for (size_t i = 0; i < num_nodes_on_entity; ++i )
                {
                    buf.unpack<EntityKey>(sentity.nodes[i]);
                }
                buf.unpack<EntityKey>(sentity.global_key);

                shared_entities_and_proc.emplace_back(ip, sentity);
            }
        }
    }
}

void BulkData::unpackEntityFromOtherProcAndUpdateInfoIfSharedLocally(stk::CommSparse &comm, std::vector<shared_entity_type> & shared_entities_this_proc)
{
    std::vector< std::pair<int, shared_entity_type> > shared_entities_received_from_other_procs;
    this->unpack_shared_entities(comm, shared_entities_received_from_other_procs);

    for(size_t i=0;i<shared_entities_received_from_other_procs.size();++i)
    {
        this->check_if_entity_from_other_proc_exists_on_this_proc_and_update_info_if_shared(shared_entities_this_proc, shared_entities_received_from_other_procs[i].first, shared_entities_received_from_other_procs[i].second);
    }
}

void BulkData::change_entity_key_to_match_owner(const std::vector<shared_entity_type> & potentially_shared_sides)
{
    for(size_t i = 0, e = potentially_shared_sides.size(); i < e; ++i)
    {
        Entity entity = potentially_shared_sides[i].entity;
        if(potentially_shared_sides[i].global_key != potentially_shared_sides[i].local_key)
            internal_change_entity_key(potentially_shared_sides[i].local_key, potentially_shared_sides[i].global_key, entity);
    }
}

void BulkData::insert_sharing_info_into_comm_map(const std::vector<shared_entity_type> & potentially_shared_sides)
{
    for(size_t i = 0, e = potentially_shared_sides.size(); i < e; ++i)
    {
        Entity entity = potentially_shared_sides[i].entity;
        for(size_t j = 0; j < potentially_shared_sides[i].sharing_procs.size(); j++)
            entity_comm_map_insert(entity, EntityCommInfo(stk::mesh::BulkData::SHARED, potentially_shared_sides[i].sharing_procs[j]));
    }
}

void BulkData::change_entity_key_and_update_sharing_info(std::vector<shared_entity_type> & potentially_shared_sides)
{
   change_entity_key_to_match_owner(potentially_shared_sides);
   insert_sharing_info_into_comm_map(potentially_shared_sides);
}

void BulkData::set_common_entity_key_and_fix_ordering_of_nodes_and_update_comm_map(std::vector<shared_entity_type> & potentially_shared_sides_this_proc)
{
    std::sort(potentially_shared_sides_this_proc.begin(), potentially_shared_sides_this_proc.end());

    std::vector<std::vector<shared_entity_type> > shared_entities_by_processor(parallel_size());
    fillVectorOfSharedEntitiesByProcessor(potentially_shared_sides_this_proc, shared_entities_by_processor);

    stk::CommSparse comm(parallel());
    impl::communicate_shared_entity_info(*this, comm, shared_entities_by_processor);
    unpackEntityFromOtherProcAndUpdateInfoIfSharedLocally(comm, potentially_shared_sides_this_proc);
    change_entity_key_and_update_sharing_info(potentially_shared_sides_this_proc);
}

void BulkData::internal_update_sharing_comm_map_and_fill_list_modified_shared_entities_of_node_rank(stk::mesh::EntityVector& shared_new)
{
    this->gather_shared_nodes(shared_new);
}

void BulkData::internal_update_sharing_comm_map_and_fill_list_modified_shared_entities_of_rank(stk::mesh::EntityRank entityRank, std::vector<Entity> &entities)
{
    entities.clear();
    if (entityRank == stk::topology::NODE_RANK) {
        internal_update_sharing_comm_map_and_fill_list_modified_shared_entities_of_node_rank(entities);
    }
    else {
        fill_shared_entities_of_rank_while_updating_sharing_info(entityRank, entities);
    }
    std::sort(entities.begin(), entities.end(), EntityLess(*this));
}

void pack_entity_keys_to_send(stk::CommSparse &comm, const std::vector<stk::mesh::EntityKeyProc> &entities_to_send_data)
{
    for(size_t i=0;i<entities_to_send_data.size();++i)
    {
        stk::mesh::EntityKey entityKeyToSend = entities_to_send_data[i].first;
        int destinationProc = entities_to_send_data[i].second;
        comm.send_buffer(destinationProc).pack(entityKeyToSend);
    }
}

void unpack_entity_keys_from_procs(stk::CommSparse &comm, std::vector<stk::mesh::EntityKey> &receivedEntityKeys)
{
    for(int procId = comm.parallel_size() - 1; procId >= 0; --procId)
    {
        if(procId != comm.parallel_rank())
        {
            CommBuffer & buf = comm.recv_buffer(procId);
            while(buf.remaining())
            {
                stk::mesh::EntityKey entityKey;
                buf.unpack<stk::mesh::EntityKey>(entityKey);
                receivedEntityKeys.push_back(entityKey);
            }
        }
    }
}

void BulkData::find_and_delete_internal_faces(stk::mesh::EntityRank entityRank, const stk::mesh::Selector *only_consider_second_element_from_this_selector)
{
    std::vector<shared_entity_type> shared_entities;
    this->markEntitiesForResolvingSharingInfoUsingNodes(entityRank, false, shared_entities);
    std::sort(shared_entities.begin(), shared_entities.end());

    // shared_edges[0] will contain all the edges this processor shares with processor 0
    std::vector<std::vector<shared_entity_type> > shared_entities_to_each_proc(parallel_size());
    fillVectorOfSharedEntitiesByProcessor(shared_entities, shared_entities_to_each_proc);

    stk::CommSparse comm(parallel());
    impl::communicate_shared_entity_info(*this, comm, shared_entities_to_each_proc);

    std::vector< std::pair<int, shared_entity_type> > shared_entities_and_proc;
    this->unpack_shared_entities(comm, shared_entities_and_proc);

    std::vector<stk::mesh::EntityKeyProc> entities_to_send_data;
    for(size_t i=0;i<shared_entities_and_proc.size();++i)
    {
        this->mark_shared_sides_and_fill_list_of_sides_not_on_boundary(shared_entities, shared_entities_and_proc[i].first, shared_entities_and_proc[i].second,
                entities_to_send_data, only_consider_second_element_from_this_selector);
    }

    std::vector<Entity> entities;
    for (size_t i=0; i<shared_entities.size(); ++i)
    {
        Entity entity = shared_entities[i].entity;
        if ( internal_is_entity_marked(entity) == BulkData::IS_SHARED && state(entity) == Created)
        {
            entities.push_back(entity);
        }
    }

    stk::CommSparse commForInternalSides(parallel());
    pack_entity_keys_to_send(commForInternalSides, entities_to_send_data);
    commForInternalSides.allocate_buffers();
    pack_entity_keys_to_send(commForInternalSides, entities_to_send_data);
    commForInternalSides.communicate();

    std::vector<stk::mesh::EntityKey> receivedEntityKeys;
    unpack_entity_keys_from_procs(commForInternalSides, receivedEntityKeys);
    for (size_t i=0; i<receivedEntityKeys.size(); ++i)
    {
        stk::mesh::Entity tempEntity = this->get_entity(receivedEntityKeys[i]);
        if ( state(tempEntity) == Created)
        {
            entities.push_back(tempEntity);
        }
    }

    stk::mesh::impl::delete_entities_and_upward_relations(*this, entities);
}

//----------------------------------------------------------------------
BulkData::BulkData(MetaData & mesh_meta_data,
                   ParallelMachine parallel,
                   enum AutomaticAuraOption auto_aura_option,
#ifdef SIERRA_MIGRATION
                   bool add_fmwk_data,
#endif
                   FieldDataManager *field_data_manager,
                   unsigned bucket_capacity)
  :
#ifdef SIERRA_MIGRATION
    m_check_invalid_rels(true),
#endif
    m_entity_comm_map(),
    m_ghosting(),
    m_mesh_meta_data( mesh_meta_data ),
    m_mark_entity(),
    m_add_node_sharing_called(false),
    m_closure_count(),
    m_mesh_indexes(),
    m_entity_repo(new impl::EntityRepository()),
    m_entity_comm_list(),
    m_entitycomm(),
    m_owner(),
    m_comm_list_updater(m_entity_comm_list, m_entitycomm),
    m_deleted_entities_current_modification_cycle(),
    m_ghost_reuse_map(),
    m_entity_keys(),
#ifdef SIERRA_MIGRATION
    m_add_fmwk_data(add_fmwk_data),
    m_fmwk_global_ids(),
    m_fmwk_aux_relations(),
    m_shouldSortFacesByNodeIds(false),
#endif
    m_autoAuraOption(auto_aura_option),
    m_meshModification(*this),
    m_parallel( parallel ),
    m_volatile_fast_shared_comm_map(),
    m_all_sharing_procs(mesh_meta_data.entity_rank_count()),
    m_ghost_parts(),
    m_deleted_entities(),
    m_num_fields(-1), // meta data not necessarily committed yet
    m_keep_fields_updated(true),
    m_local_ids(),
    m_default_field_data_manager(mesh_meta_data.entity_rank_count()),
    m_field_data_manager(field_data_manager),
    m_selector_to_buckets_map(),
    m_bucket_repository(
        *this,
        mesh_meta_data.entity_rank_count(),
        (mesh_meta_data.spatial_dimension() == 2 ? ConnectivityMap::default_map_2d() : ConnectivityMap::default_map()),
        bucket_capacity),
    m_use_identifiers_for_resolving_sharing(false),
    m_modSummary(*this),
    m_meshDiagnosticObserver(std::make_shared<stk::mesh::MeshDiagnosticObserver>(*this)),
    m_sideSetData(*this),
    m_ngpMeshManager(nullptr),
    m_isDeviceMeshRegistered(false),
    m_soloSideIdGenerator(stk::parallel_machine_size(parallel), stk::parallel_machine_rank(parallel), std::numeric_limits<stk::mesh::EntityId>::max()),
    m_supportsLargeIds(false)
{
  try {
     mesh_meta_data.set_mesh_bulk_data(this);
  }
  catch(...) {
      delete m_entity_repo;
      throw;
  }

  m_entity_comm_map.setCommMapChangeListener(&m_comm_list_updater);

  if (m_field_data_manager == NULL)
  {
      m_field_data_manager = &m_default_field_data_manager;
  }

  initialize_arrays();

  m_ghost_parts.clear();
  internal_create_ghosting( "shared" );
  //shared part should reside in m_ghost_parts[0]
  internal_create_ghosting( "shared_aura" );

  impl::check_size_of_types();

  register_observer(m_meshDiagnosticObserver);

//  m_meshDiagnosticObserver->enable_rule(stk::mesh::RULE_3);
//  m_meshDiagnosticObserver->enable_rule(stk::mesh::SOLO_SIDES);

#ifndef NDEBUG
  //m_meshDiagnosticObserver->enable_rule(stk::mesh::RULE_3);
#endif

  m_meshModification.set_sync_state_synchronized();
}

BulkData::~BulkData()
{

#ifdef STK_PROFILE_MEMORY
  ParallelMachine world = MPI_COMM_WORLD; // HACK, but necessary to work with Fmwk
  const int real_rank = parallel_machine_rank(world);
  print_max_stk_memory_usage(world, real_rank, std::cout);
#endif

#ifdef SIERRA_MIGRATION
  for(size_t i=0; i<m_fmwk_aux_relations.size(); ++i) {
    delete m_fmwk_aux_relations[i];
  }
#endif

  while ( ! m_ghosting.empty() ) {
    delete m_ghosting.back();
    m_ghosting.pop_back();
  }

  mesh_meta_data().set_mesh_bulk_data(NULL);

  delete m_elemElemGraph;
  delete m_entity_repo;

  delete m_ngpMeshManager;
}

NgpMesh &
BulkData::get_updated_ngp_mesh() const
{
  update_ngp_mesh();

  return m_ngpMeshManager->get_mesh();
}

void
BulkData::update_ngp_mesh() const
{
  if (m_ngpMeshManager == nullptr) {
#ifdef STK_USE_DEVICE_MESH
    m_ngpMeshManager = new DeviceMeshManager(*this);
#else
    m_ngpMeshManager = new HostMeshManager(*this);
#endif
  }

  m_ngpMeshManager->update_mesh();
}

void
BulkData::register_device_mesh() const
{
  ThrowRequireMsg(m_isDeviceMeshRegistered==false, "Cannot register more than one DeviceMesh against a BulkData");
  m_isDeviceMeshRegistered = true;
}

void
BulkData::unregister_device_mesh() const
{
  m_isDeviceMeshRegistered = false;

  const stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(mesh_meta_data().entity_rank_count());
  for (stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank < endRank; ++rank) {
    const stk::mesh::BucketVector & stkBuckets = buckets(rank);
    for (Bucket * bucket : stkBuckets) {
      bucket->set_ngp_bucket_id(INVALID_BUCKET_ID);
    }
  }
}

void BulkData::update_deleted_entities_container()
{
  while(!m_deleted_entities_current_modification_cycle.empty()) {
    Entity::entity_value_type entity_offset = m_deleted_entities_current_modification_cycle.front();
    m_deleted_entities_current_modification_cycle.pop_front();
    m_deleted_entities.push_front(entity_offset);
  }

  // Reclaim offsets for deleted ghosts that were not regenerated
  for (auto keyAndOffset : m_ghost_reuse_map) {
    m_deleted_entities.push_front(keyAndOffset.second);
  }

  m_ghost_reuse_map.clear();
}

size_t BulkData::total_field_data_footprint(EntityRank rank) const
{
  const std::vector< FieldBase * > & field_set = mesh_meta_data().get_fields();

  size_t retval = 0;
  for ( int i = 0; i < m_num_fields; ++i) {
    const FieldBase  & field = * field_set[i];
    retval += total_field_data_footprint(field, rank);
  }

  return retval;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void BulkData::require_ok_to_modify() const
{
  ThrowRequireMsg( !this->in_synchronized_state(),
                   "NOT in the ok-to-modify state" );
}

void BulkData::require_entity_owner( const Entity entity ,
                                     int owner ) const
{
  if (parallel_size() > 1 && bucket_ptr(entity) != NULL) {
    const bool error_not_owner = owner != parallel_owner_rank(entity) ;

    ThrowRequireMsg( !error_not_owner,
                     "P" << parallel_rank() << " " << entity_key(entity) << " owner is " <<
                     parallel_owner_rank(entity) << ", expected " << owner);
  }
}

void BulkData::require_good_rank_and_id(EntityRank ent_rank, EntityId ent_id) const
{
  const size_t rank_count = m_mesh_meta_data.entity_rank_count();
  const bool ok_id   = EntityKey::is_valid_id(ent_id);
  const bool ok_rank = ent_rank < rank_count && !(ent_rank == stk::topology::FACE_RANK && mesh_meta_data().spatial_dimension() == 2);

  ThrowRequireMsg( ok_rank,
                   "Bad key rank: " << ent_rank << " for id " << ent_id );

  ThrowRequireMsg( ok_id, "Bad id : " << ent_id);
}

bool BulkData::is_good_rank_and_id(EntityRank ent_rank, EntityId ent_id) const
{
  const size_t rank_count = m_mesh_meta_data.entity_rank_count();
  const bool ok_id   = EntityKey::is_valid_id(ent_id);
  const bool ok_rank = ent_rank < rank_count && !(ent_rank == stk::topology::FACE_RANK && mesh_meta_data().spatial_dimension() == 2);

  return ok_id && ok_rank;
}

void BulkData::mark_entity_and_upward_related_entities_as_modified(Entity entity)
{
  impl::OnlyVisitUnchanged ovu(*this);
  BulkData& mesh = *this;

  auto markAsModified = [&](Entity ent) { mesh.set_state(ent, Modified); };

  auto onlyVisitUnchanged = [&](Entity ent) { return mesh.state(ent) == Unchanged; };

  impl::VisitUpwardClosureGeneral(*this, entity, markAsModified, onlyVisitUnchanged);
}

size_t BulkData::count_relations(Entity entity) const
{
  const MeshIndex &mesh_idx = mesh_index(entity);

  const EntityRank end_rank = static_cast<EntityRank>(m_mesh_meta_data.entity_rank_count());
  size_t count = 0;
  for (EntityRank irank = stk::topology::BEGIN_RANK; irank < end_rank; ++irank)
  {
    count += mesh_idx.bucket->num_connectivity(mesh_idx.bucket_ordinal, irank);
  }
  return count;
}

bool BulkData::has_no_relations(Entity entity) const
{
  const MeshIndex &mesh_idx = mesh_index(entity);

  const EntityRank end_rank = static_cast<EntityRank>(m_mesh_meta_data.entity_rank_count());
  for (EntityRank irank = stk::topology::BEGIN_RANK; irank < end_rank; ++irank)
  {
    if (mesh_idx.bucket->num_connectivity(mesh_idx.bucket_ordinal, irank) > 0)
    {
      return false;
    }
  }
  return true;
}

unsigned BulkData::count_valid_connectivity(Entity entity, EntityRank rank) const
{
  if (bucket(entity).connectivity_type(rank) == FIXED_CONNECTIVITY) {

    m_check_invalid_rels = false;
    Entity const *rel_iter = begin(entity, rank);
    Entity const *rel_end = end(entity, rank);
    m_check_invalid_rels = true;

    unsigned count = 0;
    for (; rel_iter != rel_end; ++rel_iter)
    {
      if (rel_iter->is_local_offset_valid())
      {
        ++count;
      }
    }
    return count;
  }
  else {
    return bucket(entity).num_connectivity(bucket_ordinal(entity), rank);
  }
}

Entity BulkData::generate_new_entity(unsigned preferred_offset)
{
  Entity::entity_value_type new_local_offset = m_mesh_indexes.size();

  if (preferred_offset != 0) {
    new_local_offset = preferred_offset;
  }
  else if (!m_deleted_entities.empty()) {
    new_local_offset = m_deleted_entities.front();
    m_deleted_entities.pop_front();
  }

  MeshIndex mesh_index = {NULL, 0};
  EntityKey invalid_key;

  if (new_local_offset == m_mesh_indexes.size()) {
    m_mesh_indexes.push_back(mesh_index);
    m_entity_keys.push_back(invalid_key);
    m_entitycomm.push_back(nullptr);
    m_owner.push_back(parallel_rank());
    m_meshModification.add_created_entity_state();
    m_mark_entity.push_back(NOT_MARKED);
    m_closure_count.push_back(static_cast<uint16_t>(0));
    m_local_ids.push_back(stk::mesh::GetInvalidLocalId());

#ifdef SIERRA_MIGRATION
    if (m_add_fmwk_data) {
      m_fmwk_aux_relations.push_back(NULL);
      m_fmwk_global_ids.push_back(0);
    }
#endif
  }
  else {
    //re-claiming space from a previously-deleted entity:

    m_mesh_indexes[new_local_offset] = mesh_index;
    m_entity_keys[new_local_offset] = invalid_key;
    m_entitycomm[new_local_offset] = nullptr;
    m_owner[new_local_offset] = parallel_rank();
    m_mark_entity[new_local_offset] = NOT_MARKED;
    m_meshModification.mark_entity_as_created(new_local_offset);
    m_closure_count[new_local_offset] = static_cast<uint16_t>(0);
    m_local_ids[new_local_offset] = stk::mesh::GetInvalidLocalId();

#ifdef SIERRA_MIGRATION
    if (m_add_fmwk_data) {
      //bulk-data allocated aux-relation vector, so delete it here.
      delete m_fmwk_aux_relations[new_local_offset];
      m_fmwk_aux_relations[new_local_offset] = NULL;
      m_fmwk_global_ids[new_local_offset] = 0;
    }
#endif
  }

  return Entity(new_local_offset);
}

void BulkData::initialize_arrays()
{
  ThrowRequireMsg((m_mesh_indexes.size() == 0) && (m_entity_keys.size() == 0),
                   "BulkData::initialize_arrays() called by something other than constructor");

  MeshIndex mesh_index = {NULL, 0};
  m_mesh_indexes.push_back(mesh_index);

  EntityKey invalid_key;
  m_entity_keys.push_back(invalid_key);
  m_entitycomm.push_back(nullptr);
  m_owner.push_back(parallel_rank());

  m_mark_entity.push_back(NOT_MARKED);
  m_closure_count.push_back(static_cast<uint16_t>(0));
  m_local_ids.push_back(stk::mesh::GetInvalidLocalId());

#ifdef SIERRA_MIGRATION
  if (m_add_fmwk_data) {
    m_fmwk_aux_relations.push_back(NULL);
    m_fmwk_global_ids.push_back(0);
  }
#endif
}

Entity BulkData::declare_node(EntityId id)
{
    return declare_entity(stk::topology::NODE_RANK, id);
}
Entity BulkData::declare_element(EntityId id)
{
    return declare_entity(stk::topology::ELEM_RANK, id);
}
Entity BulkData::declare_edge(EntityId id)
{
    ThrowRequireMsg(mesh_meta_data().side_rank() != stk::topology::EDGE_RANK, "Program Error!");
    return declare_entity(stk::topology::EDGE_RANK, id);
}
Entity BulkData::declare_constraint(EntityId id)
{
    return declare_entity(stk::topology::CONSTRAINT_RANK, id);
}

Entity BulkData::declare_entity( EntityRank ent_rank , EntityId ent_id)
{
    ThrowRequireMsg(ent_rank != mesh_meta_data().side_rank(),"declare_entity for side is not supported. Please use declare_element_side() functionality.");
    return internal_declare_entity(ent_rank, ent_id, ConstPartVector{&mesh_meta_data().universal_part()});
}

template<typename PARTVECTOR>
Entity BulkData::declare_node(EntityId id, const PARTVECTOR& parts)
{
    return declare_entity(stk::topology::NODE_RANK, id, parts);
}

template
Entity BulkData::declare_node(EntityId id, const PartVector& parts);
template
Entity BulkData::declare_node(EntityId id, const ConstPartVector& parts);

template<typename PARTVECTOR>
Entity BulkData::declare_edge(EntityId id, const PARTVECTOR& parts)
{
    ThrowRequireMsg(mesh_meta_data().side_rank() != stk::topology::EDGE_RANK, "Program Error!");
    return declare_entity(stk::topology::EDGE_RANK, id, parts);
}

template
Entity BulkData::declare_edge(EntityId id, const PartVector& parts);
template
Entity BulkData::declare_edge(EntityId id, const ConstPartVector& parts);

template<typename PARTVECTOR>
Entity BulkData::declare_element(EntityId id, const PARTVECTOR& parts)
{
    return declare_entity(stk::topology::ELEM_RANK, id, parts);
}

template
Entity BulkData::declare_element(EntityId id, const PartVector& parts);
template
Entity BulkData::declare_element(EntityId id, const ConstPartVector& parts);

template<typename PARTVECTOR>
Entity BulkData::declare_constraint(EntityId id, const PARTVECTOR& parts)
{
    return declare_entity(stk::topology::CONSTRAINT_RANK, id, parts);
}

template
Entity BulkData::declare_constraint(EntityId id, const PartVector& parts);
template
Entity BulkData::declare_constraint(EntityId id, const ConstPartVector& parts);

Entity BulkData::declare_entity( EntityRank ent_rank , EntityId ent_id , Part& part)
{
    return declare_entity( ent_rank, ent_id, ConstPartVector{&part});
}

template<typename PARTVECTOR>
Entity BulkData::declare_entity( EntityRank ent_rank , EntityId ent_id ,
                                 const PARTVECTOR & parts )
{
    return internal_declare_entity(ent_rank, ent_id, parts);
}

template
Entity BulkData::declare_entity(EntityRank ent_rank, EntityId ent_id, const PartVector& parts);
template
Entity BulkData::declare_entity(EntityRank ent_rank, EntityId ent_id, const ConstPartVector& parts);

Entity BulkData::declare_solo_side( const PartVector & parts )
{
    EntityId id = m_soloSideIdGenerator.get_solo_side_id();
    return internal_declare_entity(mesh_meta_data().side_rank(), id, parts);
}

Entity BulkData::declare_solo_side( EntityId ent_id, const PartVector & parts )
{
    return internal_declare_entity(mesh_meta_data().side_rank(), ent_id, parts);
}

stk::mesh::EntityId BulkData::select_side_id(Entity elem, unsigned sideOrdinal)
{
    stk::mesh::EntityId globalSideId = impl::side_id_formula(identifier(elem), sideOrdinal);
    if(has_face_adjacent_element_graph())
    {
        stk::mesh::ElemElemGraph &graph = get_face_adjacent_element_graph();
        stk::mesh::SideIdChooser sideIdChooser = graph.get_side_id_chooser();
        globalSideId = sideIdChooser.get_chosen_side_id(elem, sideOrdinal);
    }
    return globalSideId;
}

template<typename PARTVECTOR>
PARTVECTOR add_root_topology_part(const PARTVECTOR &parts, stk::mesh::Part &rootTopoPart)
{
    PARTVECTOR initialParts(parts.size() + 1);
    initialParts = parts;
    initialParts.push_back(&rootTopoPart);
    return initialParts;
}

void use_graph_to_connect_side(stk::mesh::ElemElemGraph &graph, Entity side, Entity elem, unsigned localSideId)
{
    stk::mesh::SideNodeConnector sideNodeConnector = graph.get_side_node_connector();
    sideNodeConnector.connect_side_to_nodes(side, elem, localSideId);

    stk::mesh::SideConnector sideConnector = graph.get_side_connector();
    sideConnector.connect_side_to_all_elements(side, elem, localSideId);
}

template<typename PARTVECTOR>
Entity BulkData::create_and_connect_side(const stk::mesh::EntityId globalSideId,
                                         Entity elem,
                                         const unsigned localSideId,
                                         const PARTVECTOR& parts)
{
    stk::topology sideTop = bucket(elem).topology().side_topology(localSideId);
    Entity side = internal_declare_entity(sideTop.rank(), globalSideId, PARTVECTOR{});
    if(has_face_adjacent_element_graph())
    {
        change_entity_parts(side, add_root_topology_part(parts, mesh_meta_data().get_topology_root_part(sideTop)));
        use_graph_to_connect_side(get_face_adjacent_element_graph(), side, elem, localSideId);
    }
    else
    {
        impl::connect_element_to_entity(*this, elem, side, localSideId, parts, sideTop);
    }
    return side;
}

template
Entity BulkData::create_and_connect_side(const stk::mesh::EntityId globalSideId,
                                         Entity elem, const unsigned localSideId,
                                         const PartVector& parts);

template
Entity BulkData::create_and_connect_side(const stk::mesh::EntityId globalSideId,
                                         Entity elem, const unsigned localSideId,
                                         const ConstPartVector& parts);

template<typename PARTVECTOR>
Entity BulkData::declare_element_side(Entity elem, const unsigned sideOrd, const PARTVECTOR& parts)
{
    stk::mesh::EntityId chosenId = select_side_id(elem, sideOrd);

    return declare_element_side_with_id(chosenId, elem, sideOrd, parts);
}

template Entity BulkData::declare_element_side<stk::mesh::PartVector>(Entity, const unsigned, const stk::mesh::PartVector&);
template Entity BulkData::declare_element_side<stk::mesh::ConstPartVector>(Entity, const unsigned, const stk::mesh::ConstPartVector&);

template<typename PARTVECTOR>
Entity BulkData::declare_element_side_with_id(const stk::mesh::EntityId globalSideId,
                                      Entity elem,
                                      const unsigned sideOrd,
                                      const PARTVECTOR& parts)
{
    impl::check_declare_element_side_inputs(*this, elem, sideOrd);

    Entity side = stk::mesh::get_side_entity_for_elem_side_pair(*this, elem, sideOrd);
    if(is_valid(side))
    {
        if(bucket(side).owned())
            change_entity_parts(side, parts, PARTVECTOR{});
    }
    else
    {
        if(has_face_adjacent_element_graph())
        {
            ThrowRequireMsg(!is_valid(get_entity(mesh_meta_data().side_rank(), globalSideId)),
                        "Conflicting id in declare_elem_side, found pre-existing side with id " << globalSideId
                        << " not connected to requested element " << identifier(elem) << ", side " << sideOrd);
        }
        side = create_and_connect_side(globalSideId, elem, sideOrd, parts);
    }
    return side;
}

template Entity BulkData::declare_element_side_with_id<stk::mesh::PartVector>(const stk::mesh::EntityId, Entity, const unsigned, const stk::mesh::PartVector&);
template Entity BulkData::declare_element_side_with_id<stk::mesh::ConstPartVector>(const stk::mesh::EntityId, Entity, const unsigned, const stk::mesh::ConstPartVector&);


//----------------------------------------------------------------------

namespace {

// A method for quickly finding an entity within a comm list
const EntityCommListInfo& find_entity(const BulkData& mesh,
                               const EntityCommListInfoVector& entities,
                               const EntityKey& key)
{
  EntityCommListInfoVector::const_iterator lb_itr = std::lower_bound(entities.begin(), entities.end(), key);
  ThrowAssertMsg(lb_itr != entities.end() && lb_itr->key == key,
                 "proc " << mesh.parallel_rank() << " Cannot find entity-key " << key << " in comm-list" );
  return *lb_itr;
}

}

void BulkData::entity_comm_list_insert(Entity node)
{
  stk::mesh::EntityKey key = entity_key(node);
  EntityCommListInfoVector::iterator lb_itr = std::lower_bound(m_entity_comm_list.begin(), m_entity_comm_list.end(), key);
  if(lb_itr == m_entity_comm_list.end() || lb_itr->key != key)
  {
    const EntityComm* entity_comm = m_entity_comm_map.entity_comm(key);
    EntityCommListInfo comm_info = {key, node, entity_comm};
    m_entity_comm_list.insert(lb_itr, comm_info);
  }
}

void BulkData::add_node_sharing( Entity node, int sharing_proc )
{
  // Only valid to specify sharing information for non-deleted nodes
  ThrowRequire(entity_rank(node) == stk::topology::NODE_RANK);
  ThrowRequire(state(node) != Deleted);

  protect_orphaned_node(node);

  m_add_node_sharing_called = true;

  if (state(node) == Unchanged)
  {
      mark_entity_and_upward_related_entities_as_modified(node);
  }

  internal_mark_entity(node, IS_SHARED);
  entity_comm_map_insert(node, EntityCommInfo(stk::mesh::BulkData::SHARED, sharing_proc));
  entity_comm_list_insert(node);
}

EntityId BulkData::get_solo_side_id()
{
    return m_soloSideIdGenerator.get_solo_side_id();
}

template<typename PARTVECTOR>
void BulkData::internal_verify_and_change_entity_parts( Entity entity,
                                                        const PARTVECTOR & add_parts ,
                                                        const PARTVECTOR & remove_parts)
{
    require_ok_to_modify();

#ifdef SIERRA_MIGRATION
    if(!m_add_fmwk_data)
    {
        require_entity_owner(entity, parallel_rank());
    }
#endif //SIERRA_MIGRATION

    OrdinalVector addPartsAndSupersets;
    impl::fill_add_parts_and_supersets(add_parts, addPartsAndSupersets);

    OrdinalVector removePartsAndSubsetsMinusPartsInAddPartsList;
    impl::fill_remove_parts_and_subsets_minus_parts_in_add_parts_list(remove_parts,
                                                      addPartsAndSupersets,
                                                      bucket(entity),
                                                      removePartsAndSubsetsMinusPartsInAddPartsList);

    OrdinalVector scratchOrdinalVec, scratchSpace;
    internal_change_entity_parts(entity,
                                 addPartsAndSupersets,
                                 removePartsAndSubsetsMinusPartsInAddPartsList,
                                 scratchOrdinalVec, scratchSpace);
}

template void BulkData::internal_verify_and_change_entity_parts(Entity, const PartVector&, const PartVector&);
template void BulkData::internal_verify_and_change_entity_parts(Entity, const ConstPartVector&, const ConstPartVector&);

template<typename PARTVECTOR>
void BulkData::internal_verify_and_change_entity_parts( const EntityVector& entities,
                                                        const PARTVECTOR & add_parts ,
                                                        const PARTVECTOR & remove_parts)
{
    require_ok_to_modify();

    OrdinalVector addPartsAndSupersets;
    OrdinalVector removePartsAndSubsetsMinusPartsInAddPartsList;
    OrdinalVector scratchOrdinalVec, scratchSpace;

    for(Entity entity : entities) {
#ifdef SIERRA_MIGRATION
      if(!m_add_fmwk_data)
      {
          require_entity_owner(entity, parallel_rank());
      }
#endif //SIERRA_MIGRATION

      addPartsAndSupersets.clear();
      impl::fill_add_parts_and_supersets(add_parts, addPartsAndSupersets);

      impl::fill_remove_parts_and_subsets_minus_parts_in_add_parts_list(remove_parts,
                                                        addPartsAndSupersets,
                                                        bucket(entity),
                                                        removePartsAndSubsetsMinusPartsInAddPartsList);

      internal_change_entity_parts(entity,
                                   addPartsAndSupersets,
                                   removePartsAndSubsetsMinusPartsInAddPartsList,
                                   scratchOrdinalVec, scratchSpace);
    }
}

template void BulkData::internal_verify_and_change_entity_parts(const EntityVector&, const PartVector&, const PartVector&);
template void BulkData::internal_verify_and_change_entity_parts(const EntityVector&, const ConstPartVector&, const ConstPartVector&);

template<typename PARTVECTOR>
Entity BulkData::internal_declare_entity( EntityRank ent_rank , EntityId ent_id ,
                                 const PARTVECTOR & parts )
{
  m_check_invalid_rels = false;

  require_ok_to_modify();

  require_good_rank_and_id(ent_rank, ent_id);

  EntityKey key( ent_rank , ent_id );

  std::pair< Entity , bool > result = internal_get_or_create_entity_with_notification( key );

  Entity declared_entity = result.first;

  if ( !result.second) {
    // An existing entity, the owner must match.
    require_entity_owner( declared_entity , parallel_rank() );
  }

  //------------------------------

  Part * const owns = & m_mesh_meta_data.locally_owned_part();

  PARTVECTOR rem ;
  PARTVECTOR add( parts );
  add.push_back( owns );

  internal_verify_and_change_entity_parts( declared_entity , add , rem );

  if ( result.second ) {
    this->internal_set_owner(declared_entity, parallel_rank());
  }

  m_check_invalid_rels = true;

  return declared_entity ;
}

template Entity BulkData::internal_declare_entity(EntityRank ent_rank, EntityId ent_id, const PartVector& parts);
template Entity BulkData::internal_declare_entity(EntityRank ent_rank, EntityId ent_id, const ConstPartVector& parts);

void BulkData::clone_solo_side_id_generator(const stk::mesh::BulkData &oldBulk)
{
    m_soloSideIdGenerator = oldBulk.m_soloSideIdGenerator;
}

bool entity_is_purely_local(const BulkData& mesh, Entity entity)
{
    const Bucket& bucket = mesh.bucket(entity);
    return bucket.owned() && !bucket.shared()
            && !bucket.in_aura() && !mesh.in_send_ghost(mesh.entity_key(entity));
}

void require_fmwk_or_entity_purely_local(const BulkData& mesh, Entity entity, const std::string& caller)
{
    ThrowRequireMsg(mesh.add_fmwk_data() || entity_is_purely_local(mesh, entity),
                   "Error, "<<caller<< " requires that stk-mesh is running under fmwk, or entity "<<mesh.entity_key(entity)<<" must be purely local.");
}

void BulkData::change_entity_id( EntityId id, Entity entity)
{
// THIS ThrowAssertMsg IS ONLY MACRO CONTROLLED TO ALLOW EXPERIMENTATION WITH
// Fmwk USING stk_parallel.  WHEN stk parallel IS USED WITHN Fmwk, THIS ASSERTION
// IS VIOLATED.
#ifndef SIERRA_MIGRATION
  ThrowAssertMsg(parallel_size() == 1,
                 "change_entity_id only supported in serial");
#endif

  require_fmwk_or_entity_purely_local(*this, entity, "BulkData::change_entity_id");

  EntityRank e_rank = entity_rank(entity);

  require_ok_to_modify();
  m_modSummary.track_change_entity_id(id, entity);

  require_good_rank_and_id(e_rank, id);

  EntityKey new_key(e_rank,id);
  EntityKey old_key = entity_key(entity);

  internal_change_entity_key(old_key, new_key, entity);
}

void BulkData::internal_change_entity_key( EntityKey old_key, EntityKey new_key, Entity entity)
{
    m_entity_repo->update_entity_key(new_key, old_key, entity);
    set_entity_key(entity, new_key);
    m_bucket_repository.set_needs_to_be_sorted(this->bucket(entity), true);
}

//----------------------------------------------------------------------

bool BulkData::destroy_entity(Entity entity, bool wasGhost)
{
    return internal_destroy_entity_with_notification(entity, wasGhost);
}

bool BulkData::internal_destroy_entity_with_notification(Entity entity, bool wasGhost)
{
    if (is_valid(entity)) {
        notifier.notify_entity_deleted(entity);
    }
    return internal_destroy_entity(entity, wasGhost);
}

bool BulkData::internal_destroy_entity(Entity entity, bool wasGhost)
{
  require_ok_to_modify();

  m_check_invalid_rels = false;
  if (!is_valid(entity)) {
    m_check_invalid_rels = true;
    return false;
  }

  const bool ghost = wasGhost || in_receive_ghost(entity);
  const stk::mesh::EntityRank erank = entity_rank(entity);

  const stk::mesh::EntityRank end_rank = static_cast<stk::mesh::EntityRank>(m_mesh_meta_data.entity_rank_count());
  for (stk::mesh::EntityRank irank = static_cast<stk::mesh::EntityRank>(erank + 1); irank != end_rank; ++irank) {
    if (num_connectivity(entity, irank) > 0) {
      m_check_invalid_rels = true;
      return false;
    }
  }

  m_modSummary.track_destroy_entity(entity);

  //-----------------------------------------------------------------------------
  // Immediately remove it from relations and buckets.
  // Postpone deletion until modification_end to be sure that
  // 1) No attempt is made to re-create it.
  // 2) Parallel index is cleaned up.
  // 3) Parallel sharing is cleaned up.
  // 4) Parallel ghosting is cleaned up.
  //
  // Must clean up the parallel lists before fully deleting the entity.

  // It is important that relations be destroyed from highest to lowest rank so
  // that the back relations are destroyed first.

  for (stk::mesh::EntityRank irank = end_rank; irank != stk::topology::BEGIN_RANK; )
  {
    --irank;

    int num_conn     = num_connectivity(entity, irank);
    Entity const* rel_entities = begin(entity, irank);
    stk::mesh::ConnectivityOrdinal const* rel_ordinals = begin_ordinals(entity, irank);

    for (int j = num_conn; j > 0; )
    {
      --j;
      if (is_valid(rel_entities[j])) {
        internal_destroy_relation(entity, rel_entities[j], rel_ordinals[j]);
        m_check_invalid_rels = false;
      }
    }
  }

  // If this is a ghosted entity, store key->local_offset so that local_offset can be
  // reused if the entity is recreated in the next aura-regen. This will prevent clients
  // from having their handles to ghosted entities go invalid when the ghost is refreshed.
  const stk::mesh::EntityKey key = entity_key(entity);
  if ( ghost ) {
    m_ghost_reuse_map[key] = entity.local_offset();
  }

  // Need to invalidate Entity handles in comm-list
  stk::mesh::EntityCommListInfoVector::iterator lb_itr =
    std::lower_bound(m_entity_comm_list.begin(), m_entity_comm_list.end(), key);
  if (lb_itr != m_entity_comm_list.end() && lb_itr->key == key) {
    lb_itr->entity = stk::mesh::Entity();
  }

  remove_entity_callback(erank, bucket(entity).bucket_id(), bucket_ordinal(entity));

  m_bucket_repository.remove_entity(mesh_index(entity));

  record_entity_deletion(entity);

  if ( !ghost ) {
    m_deleted_entities_current_modification_cycle.push_front(entity.local_offset());
  }

  m_check_invalid_rels = true;
  return true ;
}

void BulkData::record_entity_deletion(Entity entity)
{
    const EntityKey key = entity_key(entity);
    set_mesh_index(entity, 0, 0);
    m_entity_repo->destroy_entity(key, entity);
    notifier.notify_local_entities_created_or_deleted(key.rank());
    notifier.notify_local_buckets_changed(key.rank());
    m_meshModification.mark_entity_as_deleted(entity.local_offset());
    m_mark_entity[entity.local_offset()] = NOT_MARKED;
    m_closure_count[entity.local_offset()] = static_cast<uint16_t>(0u);
}

size_t get_max_num_ids_needed_across_all_procs(const stk::mesh::BulkData& bulkData, size_t numIdsNeededThisProc)
{
    size_t maxNumNeeded = 0;
    ThrowRequireMsg(MPI_Allreduce(&numIdsNeededThisProc, &maxNumNeeded, 1, sierra::MPI::Datatype<size_t>::type(), MPI_MAX, bulkData.parallel()) == MPI_SUCCESS,
            "Program error (MPI_Allreduce failure). Please contact sierra-help@sandia.gov for support.");
    return maxNumNeeded;
}

//----------------------------------------------------------------------

std::vector<uint64_t> BulkData::internal_get_ids_in_use(stk::topology::rank_t rank, const std::vector<stk::mesh::EntityId>& reserved_ids) const
{
    std::vector<uint64_t> ids_in_use;
    ids_in_use.reserve(m_entity_keys.size() + m_deleted_entities_current_modification_cycle.size());

    for (size_t i=0; i<m_entity_keys.size(); ++i)
    {
        if ( stk::mesh::EntityKey::is_valid_id(m_entity_keys[i].id()) && m_entity_keys[i].rank() == rank )
        {
            ids_in_use.push_back(m_entity_keys[i].id());
        }
    }

    for (Entity::entity_value_type local_offset : m_deleted_entities_current_modification_cycle)
    {
        stk::mesh::Entity entity;
        entity.set_local_offset(local_offset);
        if ( is_valid(entity) && entity_rank(entity) == rank )
        {
            ids_in_use.push_back(entity_key(entity).id());
        }
    }

    ids_in_use.insert(ids_in_use.end(), reserved_ids.begin(), reserved_ids.end());

    stk::util::sort_and_unique(ids_in_use);
    return ids_in_use;
}

uint64_t  BulkData::get_max_allowed_id() const {
#ifdef SIERRA_MIGRATION
  if(m_add_fmwk_data) {
    return std::numeric_limits<FmwkId>::max();
  } else {
    return stk::mesh::EntityKey::MAX_ID;
  }
#else
  return stk::mesh::EntityKey::MAX_ID;
#endif
}

void BulkData::generate_new_ids_given_reserved_ids(stk::topology::rank_t rank, size_t numIdsNeeded, const std::vector<stk::mesh::EntityId>& reserved_ids, std::vector<stk::mesh::EntityId>& requestedIds) const
{
    size_t maxNumNeeded = get_max_num_ids_needed_across_all_procs(*this, numIdsNeeded);
    if ( maxNumNeeded == 0 ) return;
    std::vector<uint64_t> ids_in_use = this->internal_get_ids_in_use(rank, reserved_ids);

    uint64_t maxAllowedId = get_max_allowed_id();

    requestedIds = generate_parallel_unique_ids(maxAllowedId, ids_in_use, numIdsNeeded, this->parallel());
}

void BulkData::generate_new_ids(stk::topology::rank_t rank, size_t numIdsNeeded, std::vector<stk::mesh::EntityId>& requestedIds)
{
    size_t maxNumNeeded = get_max_num_ids_needed_across_all_procs(*this, numIdsNeeded);
    if ( maxNumNeeded == 0 ) return;

    if(rank == mesh_meta_data().side_rank())
    {
        requestedIds.resize(numIdsNeeded);
        for(size_t i = 0; i < numIdsNeeded; i++)
            requestedIds[i] = m_soloSideIdGenerator.get_solo_side_id();
    }
    else
    {
        std::vector<uint64_t> ids_in_use = internal_get_ids_in_use(rank);
        uint64_t maxAllowedId = get_max_allowed_id();
        requestedIds = generate_parallel_unique_ids(maxAllowedId, ids_in_use, numIdsNeeded, parallel());
    }
}

void BulkData::generate_new_entities(const std::vector<size_t>& requests,
                                 std::vector<Entity>& requested_entities)
// requests = number of nodes needed, number of elements needed, etc.
{
    size_t numRanks = requests.size();

    std::vector< std::vector<EntityId> > requestedIds(numRanks);

    for (size_t i=0;i<numRanks;++i)
    {
        stk::topology::rank_t rank = static_cast<stk::topology::rank_t>(i);
        generate_new_ids(rank, requests[i], requestedIds[i]);
    }

    //generating 'owned' entities
    Part * const owns = &m_mesh_meta_data.locally_owned_part();

    PartVector addParts;
    addParts.push_back(owns);

    requested_entities.clear();

    for (size_t i=0;i<numRanks;++i)
    {
        stk::topology::rank_t rank = static_cast<stk::topology::rank_t>(i);
        std::vector<Entity> new_entities;
        declare_entities(rank, requestedIds[i], addParts, new_entities);
        requested_entities.insert(requested_entities.end(), new_entities.begin(), new_entities.end());
    }
}

std::pair<Entity, bool> BulkData::internal_create_entity(EntityKey key, size_t preferred_offset)
{
    m_modSummary.track_declare_entity(key.rank(), key.id(), stk::mesh::PartVector());

    std::pair<entity_iterator ,bool> entityBoolPair = m_entity_repo->internal_create_entity(key);

    if(entityBoolPair.second)
    {
        entityBoolPair.first->second = this->generate_new_entity(preferred_offset);
        this->set_entity_key(entityBoolPair.first->second, key);
        notifier.notify_local_entities_created_or_deleted(key.rank());
        notifier.notify_local_buckets_changed(key.rank());
    }

    return std::make_pair(entityBoolPair.first->second, entityBoolPair.second);
}

std::pair<Entity, bool> BulkData::internal_get_or_create_entity_with_notification(EntityKey key, size_t preferred_offset)
{
    std::pair<Entity, bool> entityBoolPair = internal_create_entity(key, preferred_offset);
    if (entityBoolPair.second) {
        notifier.notify_entity_added(entityBoolPair.first);
    }
    return entityBoolPair;
}

template<typename IDVECTOR>
void BulkData::declare_entities(stk::topology::rank_t rank, const IDVECTOR& newIds,
       const PartVector &parts, std::vector<Entity>& requested_entities)
{
    requested_entities.resize(newIds.size());
    if (newIds.empty()) {
        return;
    }

    OrdinalVector partsAndSupersets;
    impl::fill_add_parts_and_supersets(parts, partsAndSupersets);
    Part * const ownedPart = & m_mesh_meta_data.locally_owned_part();
    impl::fill_add_parts_and_supersets(ConstPartVector{ownedPart}, partsAndSupersets);

    std::vector<std::pair<stk::mesh::EntityId,unsigned>> sortedIds(newIds.size());
    for(unsigned i=0; i<newIds.size(); ++i) {
        sortedIds[i].first = newIds[i];
        sortedIds[i].second = i;
        requested_entities[i] = this->generate_new_entity();
    }
    std::sort(sortedIds.begin(), sortedIds.end());

    OrdinalVector rem, parts_removed;
    OrdinalVector newBucketPartList;
    OrdinalVector inducible_parts_added, inducible_parts_removed, scratchSpace;

    internal_fill_new_part_list_and_removed_part_list(requested_entities[0], partsAndSupersets, rem,
                                                      newBucketPartList, parts_removed);
    internal_determine_inducible_parts(requested_entities[0], partsAndSupersets, parts_removed,
                                       inducible_parts_added, inducible_parts_removed);

    notifier.notify_local_entities_created_or_deleted(rank);
    notifier.notify_local_buckets_changed(rank);

    for(size_t i=0;i<sortedIds.size();++i)
    {
        EntityKey key(rank, sortedIds[i].first);
        require_good_rank_and_id(key.rank(), key.id());

        m_modSummary.track_declare_entity(key.rank(), key.id(), stk::mesh::PartVector());
    
        std::pair<entity_iterator ,bool> entityBoolPair = m_entity_repo->internal_create_entity(key);
    
        ThrowErrorMsgIf( ! entityBoolPair.second,
                "Generated id " << key.id() << " of rank " << key.rank() << " which was already used.");

        entityBoolPair.first->second = requested_entities[sortedIds[i].second];
        Entity new_entity = entityBoolPair.first->second;
        this->set_entity_key(new_entity, key);
        notifier.notify_entity_added(new_entity);

        notifier.notify_entity_parts_added(new_entity, partsAndSupersets);

        internal_adjust_closure_count(new_entity, partsAndSupersets, rem);
        internal_move_entity_to_new_bucket(new_entity, newBucketPartList, scratchSpace);
        internal_propagate_induced_part_changes_to_downward_connected_entities(new_entity, inducible_parts_added, inducible_parts_removed);

        this->internal_set_owner(new_entity, parallel_rank());
    }
    //now clear the created-entity cache...
    begin_entities(rank);
}

template
void BulkData::declare_entities(stk::topology::rank_t rank, const std::vector<int>& newIds, const PartVector& parts, std::vector<Entity>& requested_entities);
template
void BulkData::declare_entities(stk::topology::rank_t rank, const std::vector<int64_t>& newIds, const PartVector& parts, std::vector<Entity>& requested_entities);

bool BulkData::in_shared(EntityKey key, int proc) const
{
  PairIterEntityComm sharing = internal_entity_comm_map_shared(key);
  for ( ; !sharing.empty(); ++sharing ) {
    if ( proc == sharing->proc ) {
      return true ;
    }
  }
  return false ;
}

bool BulkData::in_shared(Entity entity, int proc) const
{
  const EntityComm* entityComm = m_entitycomm[entity.local_offset()];
  if (entityComm != nullptr && entityComm->isShared) {
    const EntityCommInfoVector& vec = entityComm->comm_map;
    EntityCommInfoVector::const_iterator it = vec.begin();
    EntityCommInfoVector::const_iterator end = vec.end();
    
    while(it!=end && it->ghost_id == SHARED) {
      if (it->proc == proc) {
        return true;
      }
      ++it;
    }
  }
  return false ;
}

bool BulkData::is_aura_ghosted_onto_another_proc( EntityKey key ) const
{
  if (key.is_valid()) {
    const int proc = parallel_rank();
    const int owner_rank = parallel_owner_rank(get_entity(key));
    if ( proc == owner_rank )
    {
        for ( PairIterEntityComm ec = internal_entity_comm_map(key); ! ec.empty() ; ++ec ) {
          if ( ec->ghost_id == BulkData::AURA &&
               ec->proc     != proc ) {
            return true;
          }
        }
    }
  }
  return false;
}

bool BulkData::in_send_ghost( EntityKey key , int proc ) const
{
  const int owner_rank = parallel_owner_rank(get_entity(key));
  for ( PairIterEntityComm ec = internal_entity_comm_map(key); ! ec.empty() ; ++ec ) {
    if ( ec->ghost_id != BulkData::SHARED &&
         ec->proc     != owner_rank &&
         ec->proc     == proc ) {
      return true;
    }
  }
  return false;
}

bool BulkData::in_ghost( const Ghosting & ghost , EntityKey key , int proc ) const
{
  // Ghost communication from owner.
  EntityCommInfo tmp( ghost.ordinal() , proc );

  PairIterEntityComm ec = internal_entity_comm_map(key);
  EntityCommInfoVector::const_iterator i =
    std::lower_bound( ec.begin(), ec.end() , tmp );

  return i != ec.end() && tmp == *i ;
}

bool BulkData::in_ghost( const Ghosting & ghost , Entity entity , int proc ) const
{
  if (m_entitycomm[entity.local_offset()] == nullptr) {
    return false;
  }

  const EntityCommInfoVector& vec = m_entitycomm[entity.local_offset()]->comm_map;
  
  EntityCommInfoVector::const_iterator i = vec.begin();
  EntityCommInfoVector::const_iterator end = vec.end();
  
  while(i!=end && i->ghost_id < ghost.ordinal()) {
    ++i;
  }

  while(i!=end && i->ghost_id == ghost.ordinal()) {
    if (i->proc == proc) {
      return true;
    }
    ++i;
  }

  return false;
}

bool BulkData::in_ghost( const Ghosting & ghost , Entity entity ) const
{
  if (m_entitycomm[entity.local_offset()] == nullptr) {
    return false;
  }

  const EntityCommInfoVector& vec = m_entitycomm[entity.local_offset()]->comm_map;
  
  EntityCommInfoVector::const_iterator i = vec.begin();
  EntityCommInfoVector::const_iterator end = vec.end();
  
  while(i!=end && i->ghost_id < ghost.ordinal()) {
    ++i;
  }

  if (i!=end && i->ghost_id == ghost.ordinal()) {
    return true;
  }

  return false;
}

bool BulkData::in_send_ghost( const Ghosting & ghost , EntityKey key , int proc ) const
{
  bool ret_val = false;
  const int owner_rank = parallel_owner_rank(get_entity(key));

  if (owner_rank == parallel_rank())
  {
      EntityCommInfo tmp( ghost.ordinal() , proc );

      PairIterEntityComm ec = internal_entity_comm_map(key);
      EntityCommInfoVector::const_iterator i =
        std::lower_bound( ec.begin(), ec.end() , tmp );

      ret_val = i != ec.end() && tmp == *i ;
  }

  return ret_val;
}

bool BulkData::is_communicated_with_proc(Entity entity, int proc) const
{
  if (m_entitycomm[entity.local_offset()] == nullptr) {
    return false;
  }

  const EntityCommInfoVector& vec = m_entitycomm[entity.local_offset()]->comm_map;
  
  EntityCommInfoVector::const_iterator i = vec.begin();
  EntityCommInfoVector::const_iterator end = vec.end();
  
  while(i != end) {
    if (i->proc == proc) {
      return true;
    }
    ++i;
  }

  return false;
}

void fill_sorted_procs(const PairIterEntityComm& ec, std::vector<int>& procs)
{
    procs.clear();
    int n = ec.size();
    for (int i=0; i<n; ++i) {
      procs.push_back( ec[i].proc );
    }
    stk::util::sort_and_unique(procs);
}

void fill_ghosting_procs(const PairIterEntityComm& ec, unsigned ghost_id, std::vector<int>& procs)
{
    procs.clear();
    int n = ec.size();
    for (int i=0; i<n; ++i) {
        if (ghost_id == ec[i].ghost_id) {
            procs.push_back( ec[i].proc );
        }
    }
}

void BulkData::comm_procs( EntityKey key, std::vector<int> & procs ) const
{
  ThrowAssertMsg(is_valid(get_entity(key)),
                  "BulkData::comm_procs ERROR, input key "<<key<<" not a valid entity. Contact sierra-help@sandia.gov");

#ifndef NDEBUG
  EntityCommListInfoVector::const_iterator lb_itr = std::lower_bound(m_entity_comm_list.begin(),
                                                                     m_entity_comm_list.end(),
                                                                     key);
  if (lb_itr != m_entity_comm_list.end() && lb_itr->key == key) {
      ThrowAssertMsg( lb_itr->entity != Entity(),
                      "comm-list contains invalid entity for key "<<key<<". Contact sierra-help@sandia.gov");
  }
#endif
  fill_sorted_procs(internal_entity_comm_map(key), procs);
}

void BulkData::comm_shared_procs(EntityKey key, std::vector<int> & procs ) const
{
    procs.clear();
    const EntityComm* entityComm = m_entity_comm_map.entity_comm(key);
    if (entityComm != nullptr && !entityComm->comm_map.empty()) {
       const EntityCommInfoVector& comm_map = entityComm->comm_map;
       for(const EntityCommInfo& commInfo : comm_map) {
           if (commInfo.ghost_id == BulkData::SHARED) {
               procs.push_back(commInfo.proc);
           }
           else {
               break;
           }
       }
    }
}

void BulkData::comm_shared_procs(Entity entity, std::vector<int> & procs ) const
{
    procs.clear();
    const EntityComm* entityComm = m_entitycomm[entity.local_offset()];
    if (entityComm != nullptr && !entityComm->comm_map.empty()) {
       const EntityCommInfoVector& comm_map = entityComm->comm_map;
       for(const EntityCommInfo& commInfo : comm_map) {
           if (commInfo.ghost_id == BulkData::SHARED) {
               procs.push_back(commInfo.proc);
           }
           else {
               break;
           }
       }
    }
}

void BulkData::shared_procs_intersection(const std::vector<EntityKey> & keys, std::vector<int> & procs ) const
{
  procs.clear();
  int num = keys.size();
  std::vector<int> procs_tmp;
  std::vector<int> result;
  for (int i = 0; i < num; ++i)
  {
    comm_shared_procs(keys[i], procs_tmp);

    if (i == 0)
      procs.swap(procs_tmp);
    else
    {
      // subsequent loops keep the intersection
      result.clear();
      std::back_insert_iterator<std::vector<int> > result_itr(result);
      std::set_intersection(procs.begin(),
                            procs.end(),
                            procs_tmp.begin(),
                            procs_tmp.end(),
                            result_itr,
                            std::less<int>());
      procs.swap(result);
    }
  }
}

void BulkData::shared_procs_intersection(const EntityVector& entities,
                                         std::vector<int> & procs ) const
{
  procs.clear();
  int num = entities.size();
  std::vector<int> procs_tmp;
  std::vector<int> result;
  for (int i = 0; i < num; ++i)
  {
    comm_shared_procs(entities[i], procs_tmp);

    if (i == 0)
      procs.swap(procs_tmp);
    else
    {
      // subsequent loops keep the intersection
      result.clear();
      std::back_insert_iterator<std::vector<int> > result_itr(result);
      std::set_intersection(procs.begin(),
                            procs.end(),
                            procs_tmp.begin(),
                            procs_tmp.end(),
                            result_itr,
                            std::less<int>());
      procs.swap(result);
    }
  }
}

void BulkData::comm_procs( const Ghosting & ghost ,
                           EntityKey key, std::vector<int> & procs ) const
{
  fill_ghosting_procs(internal_entity_comm_map(key), ghost.ordinal(), procs);
}

void BulkData::internal_set_owner(Entity entity, int new_owner)
{
  m_owner[entity.local_offset()] = new_owner;
}

void BulkData::internal_sync_comm_list_owners()
{
}

void BulkData::deactivate_field_updating()
{
  if (m_num_fields > -1) {
    //if fields have already been allocated, then we can't deactivate the updating
    //of field-data.
    m_keep_fields_updated = true;
    return;
  }

  m_keep_fields_updated = false;
}

void BulkData::allocate_field_data()
{
  if (m_keep_fields_updated == true) {
    //fields are already allocated, nothing to do here.
    return;
  }

  //temporary (hopefully) kludge:
  //calling the buckets(rank) getter causes partitions/buckets to potentially
  //be reorganized (including deleting buckets) and so we need to do it
  //before flipping the m_keep_fields_updated flag...
  for(EntityRank rank = stk::topology::NODE_RANK; rank < mesh_meta_data().entity_rank_count(); ++rank) {
    this->buckets(rank);
  }

  m_keep_fields_updated = true;
  //now loop over all buckets and call the 'new_bucket_callback' method which
  //will allocate field-data for that bucket.

  const std::vector< FieldBase * > & field_set = mesh_meta_data().get_fields();
  if (m_num_fields == -1) {
    // hasn't been set yet
    m_num_fields = field_set.size();
  }

  for(EntityRank rank = stk::topology::NODE_RANK; rank < mesh_meta_data().entity_rank_count(); ++rank) {
      const std::vector<Bucket*>& buckets = this->buckets(rank);
      m_field_data_manager->allocate_field_data(rank, buckets, field_set);
  }
}

void BulkData::reallocate_field_data(stk::mesh::FieldBase & field)
{
  const std::vector<FieldBase *> & field_set = mesh_meta_data().get_fields();
  for(EntityRank rank = stk::topology::NODE_RANK; rank < mesh_meta_data().entity_rank_count(); ++rank) {
    const std::vector<Bucket*>& buckets = this->buckets(rank);
    m_field_data_manager->reallocate_field_data(rank, buckets, field, field_set);
  }
}

void BulkData::register_observer(std::shared_ptr<ModificationObserver> observer) const
{
    notifier.register_observer(observer);
}

void BulkData::unregister_observer(std::shared_ptr<ModificationObserver> observer) const
{
    notifier.unregister_observer(observer);
}

void BulkData::new_bucket_caching(EntityRank rank, Bucket* new_bucket)
{
    // update selector map
    if (new_bucket != NULL) {
      for (SelectorBucketMap::iterator itr = m_selector_to_buckets_map.begin(), end = m_selector_to_buckets_map.end();
           itr != end; ++itr) {
        Selector const& sel = itr->first.second;
        const EntityRank map_rank = itr->first.first;
        if (map_rank == rank && sel(*new_bucket)) {
          BucketVector & cached_buckets = itr->second;
          BucketVector::iterator lb_itr = std::lower_bound(cached_buckets.begin(), cached_buckets.end(), new_bucket, BucketIdComparator());
          cached_buckets.insert(lb_itr, new_bucket);
        }
      }
    }
}

void BulkData::new_bucket_callback(EntityRank rank, const PartVector& superset_parts, size_t capacity, Bucket* new_bucket)
{
  this->new_bucket_caching(rank, new_bucket);

  if (!m_keep_fields_updated) {
    return;
  }

  const std::vector< FieldBase * > & field_set = mesh_meta_data().get_fields();

  if (m_num_fields == -1) {
    // hasn't been set yet
    m_num_fields = field_set.size();
  }

  m_field_data_manager->allocate_bucket_field_data(rank, field_set, superset_parts, capacity);
}

//
//  Copy fields from src to dst entity.  If the field size of either entity is zero, do nothing.  If the field
//  size of of both entities are non-zero, then the sizes must match
//

void BulkData::copy_entity_fields_callback(EntityRank dst_rank, unsigned dst_bucket_id, Bucket::size_type dst_bucket_ord,
                                           unsigned src_bucket_id, Bucket::size_type src_bucket_ord, const std::vector<FieldBase*>* field_set)
{
    //
    //  If field set is passed in copy only the defined fields.  Also assume the fields are valid for the bucket
    //
    if(field_set)
    {
        for(int i = 0, iend = field_set->size(); i < iend; ++i)
        {
            const int src_size = (*field_set)[i]->get_meta_data_for_field()[src_bucket_id].m_bytes_per_entity;
            unsigned char * const src = (*field_set)[i]->get_meta_data_for_field()[src_bucket_id].m_data;
            unsigned char * const dst = (*field_set)[i]->get_meta_data_for_field()[dst_bucket_id].m_data;

            ThrowAssert(src_size == (*field_set)[i]->get_meta_data_for_field()[dst_bucket_id].m_bytes_per_entity);

            std::memcpy(dst + src_size * dst_bucket_ord,
                    src + src_size * src_bucket_ord,
                    src_size);
        }
    }
    else
    {
        if(!m_keep_fields_updated)
        {
            return;
        }

        const std::vector<FieldBase *>& allFields = mesh_meta_data().get_fields(dst_rank);
        for(int i = 0, iend = allFields.size(); i < iend; ++i)
        {
            const int src_size = allFields[i]->get_meta_data_for_field()[src_bucket_id].m_bytes_per_entity;
            if(src_size == 0)
            {
                continue;
            }

            unsigned char * const src = allFields[i]->get_meta_data_for_field()[src_bucket_id].m_data;
            const int dst_size = allFields[i]->get_meta_data_for_field()[dst_bucket_id].m_bytes_per_entity;

            if(dst_size)
            {
                unsigned char * const dst = allFields[i]->get_meta_data_for_field()[dst_bucket_id].m_data;
                ThrowAssertMsg( dst_size == src_size,
                        "Incompatible field sizes: " << dst_size << " != " << src_size);

                std::memcpy(dst + dst_size * dst_bucket_ord,
                        src + src_size * src_bucket_ord,
                        dst_size);
            }
        }
    }
}

void BulkData::add_entity_callback(EntityRank rank, unsigned bucket_id, Bucket::size_type bucket_ord)
{
    if(!m_keep_fields_updated)
    {
        return;
    }
    const std::vector<FieldBase *> &fields = mesh_meta_data().get_fields();
    m_field_data_manager->add_field_data_for_entity(fields, rank, bucket_id, bucket_ord);
}

void BulkData::remove_entity_field_data_callback(EntityRank rank, unsigned bucket_id, Bucket::size_type bucket_ord)
{
    if (!m_keep_fields_updated) {
      return;
    }
    const std::vector<FieldBase *> &fields = mesh_meta_data().get_fields();
    m_field_data_manager->remove_field_data_for_entity(rank, bucket_id, bucket_ord, fields);
}

void BulkData::remove_entity_callback(EntityRank rank, unsigned bucket_id, Bucket::size_type bucket_ord)
{
}

void BulkData::destroy_bucket_callback(EntityRank rank, Bucket const& dying_bucket, unsigned capacity)
{
  // Remove destroyed bucket out of memoized get_buckets result, but
  // don't bother if the mesh is being destructed.
  const unsigned bucket_id = dying_bucket.bucket_id();

  if (!m_bucket_repository.being_destroyed()) {
    for (SelectorBucketMap::iterator itr = m_selector_to_buckets_map.begin(), end = m_selector_to_buckets_map.end();
         itr != end; ++itr) {
      Selector const& sel = itr->first.second;
      const EntityRank map_rank = itr->first.first;
      if (map_rank == rank && sel(dying_bucket)) {
        BucketVector & cached_buckets = itr->second;
        BucketVector::iterator lb_itr = std::lower_bound(cached_buckets.begin(), cached_buckets.end(), bucket_id, BucketIdComparator());
        ThrowAssertMsg(lb_itr != cached_buckets.end() && (*lb_itr)->bucket_id() == bucket_id,
                       "Error, bucket id " << bucket_id << ":\n " << dying_bucket << "\nWas selected by selector " << sel << " but was not in the cache");
        cached_buckets.erase(lb_itr);
      }
    }
  }

  if (!m_keep_fields_updated) {
    return;
  }

  const std::vector<FieldBase*>&  fields = mesh_meta_data().get_fields();
  m_field_data_manager->deallocate_bucket_field_data(rank, bucket_id, capacity, fields);
}

void BulkData::update_field_data_states(FieldBase* field)
{
  const int numStates = field->number_of_states();
  if (numStates > 1) {
#ifdef STK_DEBUG_FIELD_SYNC
    ThrowRequireMsg(!field->has_device_field_debug_data(), "Field state rotation currently unsupported for Field sync debugging");
#endif

    field->rotate_multistate_data();

    unsigned fieldOrdinal = field->mesh_meta_data_ordinal();
    for (int s = 1; s < numStates; ++s) {
      m_field_data_manager->swap_fields(fieldOrdinal+s, fieldOrdinal);
    }
  }
}

void BulkData::update_field_data_states()
{
  const std::vector<FieldBase*> & fields = mesh_meta_data().get_fields();

  for (int i = 0; i < m_num_fields; ) {
    FieldBase* field = fields[i];
    const int numStates = field->number_of_states();
    if (numStates > 1) {
      update_field_data_states(field);
    }
    i += numStates ;
  }
}

void sync_to_host_and_mark_modified(MetaData& meta)
{
  const std::vector<FieldBase*>& fields = meta.get_fields();
  for(FieldBase* field : fields) {
    if (field->number_of_states() > 1) {
      field->sync_to_host();
      field->modify_on_host();
    }
  }
}

const_entity_iterator BulkData::begin_entities(EntityRank ent_rank) const
{
  return m_entity_repo->begin_rank(ent_rank);
}

const_entity_iterator BulkData::end_entities(EntityRank ent_rank) const
{
  return m_entity_repo->end_rank(ent_rank);
}

Entity BulkData::get_entity( EntityRank ent_rank , EntityId entity_id ) const
{
  if (!is_good_rank_and_id(ent_rank, entity_id)) {
      return Entity();
  }
  return m_entity_repo->get_entity( EntityKey(ent_rank, entity_id));
}

Entity BulkData::get_entity( const EntityKey key ) const
{
  return m_entity_repo->get_entity(key);
}

void BulkData::reorder_buckets_callback(EntityRank rank, const std::vector<unsigned>& reorderedBucketIds)
{
  for (SelectorBucketMap::iterator itr = m_selector_to_buckets_map.begin(), end = m_selector_to_buckets_map.end();
       itr != end; ++itr) {
    BucketVector& cached_buckets = itr->second;
    std::sort(cached_buckets.begin(), cached_buckets.end(), BucketIdComparator());
  }

  if (!m_keep_fields_updated) {
    return;
  }

  const std::vector<FieldBase*>  fields = mesh_meta_data().get_fields();
  m_field_data_manager->reorder_bucket_field_data(rank, fields, reorderedBucketIds);
}

void BulkData::dump_mesh_bucket_info(std::ostream& out, Bucket* bucket) const
{
    impl::print_bucket_parts(*this, bucket, out);

    for (size_t b_ord = 0, b_end = bucket->size(); b_ord < b_end; ++b_ord) {
        MeshIndex meshIndex(bucket, b_ord);
        impl::print_entity_offset_and_state(*this, meshIndex, out);
        impl::print_entity_connectivity(*this, meshIndex, out);
        if (m_num_fields > 0) {
            impl::print_field_data_for_entity(*this, meshIndex, out);
        }
        impl::print_comm_data_for_entity(*this, entity_key((*bucket)[b_ord]), out);
    }
}

void BulkData::dump_all_mesh_info(std::ostream& out) const
{
    // Dump output for metadata first
    m_mesh_meta_data.dump_all_meta_info(out);

    out << "BulkData "
            << " info...(ptr=" << this << ")\n";

    // Iterate all buckets for all ranks...
    const std::vector<std::string> & rank_names = m_mesh_meta_data.entity_rank_names();
    for (size_t i = 0, e = rank_names.size(); i < e; ++i) {
        EntityRank rank = static_cast<EntityRank>(i);
        out << "  All " << rank_names[i] << " entities:" << std::endl;

        const BucketVector& buckets = this->buckets(rank);
        for(Bucket* bucket : buckets) {
            dump_mesh_bucket_info(out, bucket);
        }
    }
}

void BulkData::reserve_relation(Entity entity, const unsigned num)
{
  if (num == 0 && aux_relations(entity).empty()) {
    RelationVector tmp;
    aux_relations(entity).swap(tmp); // clear memory of m_relations.
  }
  else {
    aux_relations(entity).reserve(num);
  }
}

void BulkData::erase_and_clear_if_empty(Entity entity, RelationIterator rel_itr)
{
  ThrowAssert(!impl::internal_is_handled_generically(rel_itr->getRelationType()));

  RelationVector& aux_rels = aux_relations(entity);
  aux_rels.erase(aux_rels.begin() + (rel_itr - aux_rels.begin())); // Need to convert to non-const iterator

  if (aux_rels.empty()) {
    reserve_relation(entity, 0);
  }
}

void BulkData::internal_verify_initialization_invariant(Entity entity)
{
#ifndef NDEBUG
  int my_global_id = global_id(entity);
  EntityKey my_key = entity_key(entity);
#endif
  ThrowAssert ( !(my_global_id < 0 && my_key.id() == static_cast<EntityId>(my_global_id)) &&
                !(my_global_id > 0 && my_key.id() != static_cast<EntityId>(my_global_id)) );
}

BucketVector const& BulkData::get_buckets(EntityRank rank, Selector const& selector) const
{
  std::pair<EntityRank, Selector> search_item = std::make_pair(rank, selector);
  SelectorBucketMap::iterator fitr =
    m_selector_to_buckets_map.find(search_item);
  if (fitr != m_selector_to_buckets_map.end()) {
    BucketVector const& rv = fitr->second;
    return rv;
  }
  else {
    BucketVector const& all_buckets_for_rank = buckets(rank); // lots of potential side effects! Need to happen before map insertion
    std::pair<SelectorBucketMap::iterator, bool> insert_rv =
      m_selector_to_buckets_map.emplace(std::make_pair(rank, selector), BucketVector() );
    ThrowAssertMsg(insert_rv.second, "Should not have already been in map");
    BucketVector& map_buckets = insert_rv.first->second;
    for (size_t i = 0, e = all_buckets_for_rank.size(); i < e; ++i) {
      if (selector(*all_buckets_for_rank[i])) {
        map_buckets.push_back(all_buckets_for_rank[i]);
      }
    }

    return reinterpret_cast<BucketVector const&>(map_buckets);
  }
}

void BulkData::get_buckets(EntityRank rank, Selector const& selector, BucketVector & output_buckets) const
{
  output_buckets.clear();

  BucketVector const& all_buckets_for_rank = buckets(rank);
  for (size_t i = 0, e = all_buckets_for_rank.size(); i < e; ++i) {
    if (selector(*all_buckets_for_rank[i])) {
      output_buckets.push_back(all_buckets_for_rank[i]);
    }
  }
}

void BulkData::get_entities(EntityRank rank, Selector const& selector, EntityVector& output_entities) const {
  output_entities.clear();
  const stk::mesh::BucketVector &bucket_ptrs = get_buckets(rank, selector);
  for(size_t ib=0, ib_end=bucket_ptrs.size(); ib<ib_end; ++ib) {
     const stk::mesh::Bucket *bucket = bucket_ptrs[ib];
     for(size_t iobj=0, iobj_end=bucket->size(); iobj<iobj_end; ++iobj) {
       output_entities.push_back((*bucket)[iobj]);
     }
  }
}

void BulkData::require_valid_relation( const char action[] ,
                                       const BulkData & mesh ,
                                       const Entity e_from ,
                                       const Entity e_to )
{
  const bool error_rank      = !(mesh.entity_rank(e_from) > mesh.entity_rank(e_to));
  const bool error_nil_from  = !mesh.is_valid(e_from);
  const bool error_nil_to    = !mesh.is_valid(e_to);

  if ( error_rank || error_nil_from || error_nil_to ) {
    std::ostringstream msg ;

    msg << "Could not " << action << " relation from entity "
        << mesh.entity_key(e_from) << " to entity " << mesh.entity_key(e_to) << "\n";

    ThrowErrorMsgIf( error_nil_from  || error_nil_to,
                     msg.str() << ", entity was destroyed");
    ThrowErrorMsgIf( error_rank, msg.str() <<
                     "A relation must be from higher to lower ranking entity");
  }
}


//----------------------------------------------------------------------
bool BulkData::internal_declare_relation(Entity e_from, Entity e_to,
                                         RelationIdentifier local_id,
                                         Permutation permut)
{
  m_modSummary.track_declare_relation(e_from, e_to, local_id, permut);

  const MeshIndex& idx = mesh_index(e_from);

  ThrowAssertMsg(local_id <= INVALID_CONNECTIVITY_ORDINAL, "local_id = " << local_id << ", max = " << (uint32_t)INVALID_CONNECTIVITY_ORDINAL);
  bool modified = idx.bucket->declare_relation(idx.bucket_ordinal, e_to, static_cast<ConnectivityOrdinal>(local_id), permut);

  if (modified)
  {
      if ( idx.bucket->owned() ) // owned entity with relation to node, true shared
      {
          unprotect_orphaned_node(e_to);
      }
      else if ( idx.bucket->in_aura() && bucket(e_to).owned() ) // aura with relation to owned node, mostly true shared
      {
          unprotect_orphaned_node(e_to);
      }

      if (idx.bucket->owned() && (idx.bucket->entity_rank() > entity_rank(e_to)) )
      {
          ++m_closure_count[e_to.local_offset()];
      }

  }
  return modified;
}

bool BulkData::check_permutation(Entity entity, Entity rel_entity, unsigned rel_ordinal, Permutation expected) const
{
    const stk::topology &entity_topo = mesh_index(entity).bucket->topology();
    const stk::topology &rel_topo    = mesh_index(rel_entity).bucket->topology();
    Entity const *entity_nodes     = begin_nodes(entity);
    Entity const *rel_entity_nodes = begin_nodes(rel_entity);

    Permutation computed_permutation = find_permutation(entity_topo, entity_nodes,
                                                        rel_topo, rel_entity_nodes, rel_ordinal);

    return computed_permutation == expected;
}

Permutation BulkData::find_permutation( const stk::topology &hr_entity_topo,
                              Entity const *hr_entity_nodes,
                              const stk::topology &side_topo,
                              Entity const *side_nodes,
                              unsigned side_ordinal) const
{
    Entity expected_nodes[100];
    switch (side_topo.rank())
    {
    case stk::topology::EDGE_RANK:
        hr_entity_topo.edge_nodes(hr_entity_nodes, side_ordinal, expected_nodes);
        break;
    case stk::topology::FACE_RANK:
        hr_entity_topo.face_nodes(hr_entity_nodes, side_ordinal, expected_nodes);
        break;
    default:
        return INVALID_PERMUTATION;
    }

    Permutation retval = INVALID_PERMUTATION;

    int permuted[100];
    const int nv = side_topo.num_nodes();
    const int np = side_topo.num_permutations() ;
    int p = 0 ;
    for ( ; p < np ; ++p ) {
      side_topo.permutation_node_ordinals(p, permuted);

      // ALAN: can we replace this with equivalent? method on topology
      int j = 0 ;
      for ( ; j < nv && side_nodes[j] == expected_nodes[permuted[j]] ; ++j );

      if ( nv == j )
      {
          retval = static_cast<Permutation>(p);
          break;
      }
    }

    return retval;
}


void BulkData::declare_relation( Entity e_from ,
                                 Entity e_to ,
                                 const RelationIdentifier local_id ,
                                 Permutation permut)
{
  OrdinalVector ordinal_scratch, scratch2, scratch3;
  internal_declare_relation(e_from, e_to, local_id, permut, ordinal_scratch, scratch2, scratch3);
}

void BulkData::declare_relation( Entity e_from , const std::vector<Entity>& to_entities)
{
#ifndef NDEBUG
  stk::topology topo = bucket(e_from).topology();
  ThrowAssertMsg(!to_entities.empty(), "ERROR, BulkData::declare_relation given empty to_entities vector.");
  EntityRank rank = entity_rank(to_entities[0]);
  ThrowAssertMsg(to_entities.size() == topo.num_sub_topology(rank), "ERROR, BulkData::declare_relation given wrong number of downward relations.");
  for(unsigned i=1; i<to_entities.size(); ++i) {
    ThrowAssertMsg(entity_rank(to_entities[i]) == rank, "ERROR, BulkData::declare_relation: downward relations must all have the same rank.");
  }
#endif

  OrdinalVector scratch1, scratch2, scratch3;
  for(unsigned i=0; i<to_entities.size(); ++i) {
    internal_declare_relation(e_from, to_entities[i], static_cast<RelationIdentifier>(i),
                              INVALID_PERMUTATION, scratch1, scratch2, scratch3);
  }
}

void BulkData::declare_relation( Entity e_from ,
                                 Entity e_to ,
                                 const RelationIdentifier local_id ,
                                 Permutation permut,
                                 OrdinalVector& scratch1,
                                 OrdinalVector& scratch2,
                                 OrdinalVector& scratch3)
{
    internal_declare_relation(e_from, e_to, local_id, permut, scratch1, scratch2, scratch3);
}

void BulkData::internal_declare_relation( Entity e_from ,
                                 Entity e_to ,
                                 const RelationIdentifier local_id ,
                                 Permutation permut,
                                 OrdinalVector& scratch1,
                                 OrdinalVector& scratch2,
                                 OrdinalVector& scratch3)
{
    require_ok_to_modify();

    require_valid_relation("declare", *this, e_from, e_to);

    // TODO: Don't throw if exact relation already exists, that should be a no-op.
    // Should be an exact match if relation of local_id already exists (e_to should be the same).
    bool is_new_relation = internal_declare_relation(e_from, e_to, local_id, permut);

    if(is_new_relation)
    {
        internal_declare_relation(e_to, e_from, local_id, permut);
    }

    // It is critical that the modification be done AFTER the relations are
    // added so that the propagation can happen correctly.
    if(is_new_relation)
    {
        this->mark_entity_and_upward_related_entities_as_modified(e_to);
        this->mark_entity_and_upward_related_entities_as_modified(e_from);
    }

    // Deduce and set new part memberships:
    scratch1.clear();

    impl::get_part_ordinals_to_induce_on_lower_ranks(*this, e_from, entity_rank(e_to), scratch1);

    OrdinalVector emptyParts;
    internal_change_entity_parts(e_to, scratch1, emptyParts, scratch2, scratch3);
}

//----------------------------------------------------------------------

void BulkData::internal_declare_relation( Entity entity ,
                                 const std::vector<Relation> & rel )
{
  OrdinalVector scratch1;
  internal_declare_relation(entity, rel, scratch1);
}

void BulkData::internal_declare_relation( Entity entity ,
                                 const std::vector<Relation> & rel,
                                 OrdinalVector& ordinal_scratch)
{
  require_ok_to_modify();

  stk::mesh::EntityRank erank = entity_rank(entity);

  OrdinalVector scratch2, scratch3;

  std::vector<Relation>::const_iterator i ;
  for ( i = rel.begin() ; i != rel.end() ; ++i ) {
    Entity e = i->entity();
    const unsigned n = i->relation_ordinal();
    const Permutation permut = i->getPermutation();
    if ( entity_rank(e) < erank ) {
      internal_declare_relation( entity , e , n, permut, ordinal_scratch, scratch2, scratch3);
    }
    else if ( erank < entity_rank(e) ) {
      internal_declare_relation( e , entity , n, permut, ordinal_scratch, scratch2, scratch3);
    }
    else {
      ThrowErrorMsg("Given entities of the same entity rank. entity is " <<
                    identifier(entity));
    }
  }
}

//----------------------------------------------------------------------

bool BulkData::destroy_relation( Entity e_from ,
                                 Entity e_to,
                                 const RelationIdentifier local_id )
{
    return internal_destroy_relation(e_from, e_to,  local_id);
}

bool BulkData::internal_destroy_relation( Entity e_from ,
                                 Entity e_to,
                                 const RelationIdentifier local_id )
{
  require_ok_to_modify();
  m_modSummary.track_destroy_relation(e_from, e_to, local_id);

  require_valid_relation( "destroy" , *this , e_from , e_to );

  const EntityRank end_rank = static_cast<EntityRank>(m_mesh_meta_data.entity_rank_count());
  const EntityRank e_to_entity_rank = entity_rank(e_to);

  //------------------------------
  // When removing a relationship may need to
  // remove part membership and set field relation pointer to NULL

  m_check_invalid_rels = false; // OK to have gaps when deleting
  OrdinalVector scratchOrdinalVec, scratchSpace;

  if ( parallel_size() < 2 || internal_entity_comm_map_shared(e_to).empty() ) {

    //------------------------------
    // 'keep' contains the parts deduced from kept relations
    // 'del'  contains the parts deduced from deleted relations
    //        that are not in 'keep'
    // Only remove these part memberships the entity is not shared.
    // If the entity is shared then wait until modificaton_end_synchronize.
    //------------------------------

    OrdinalVector del, keep, empty;

    // For all relations that are *not* being deleted, add induced parts for
    // these relations to the 'keep' vector
    {
      Entity const* rel_entities = NULL;
      ConnectivityOrdinal const* rel_ordinals = NULL;
      int num_rels = 0;
      for (EntityRank irank = static_cast<EntityRank>(e_to_entity_rank + 1); irank < end_rank; ++irank) {
        num_rels     = num_connectivity(e_to, irank);
        rel_entities = begin(e_to, irank);
        rel_ordinals = begin_ordinals(e_to, irank);

        for (int j = 0; j < num_rels; ++j) {
          ThrowAssertMsg(is_valid(rel_entities[j]), "Error, entity " << e_to.local_offset() << " with key " << entity_key(e_to) << " has invalid back-relation for ordinal: "
                         << (uint32_t)rel_ordinals[j] << " to rank: " << irank << ", target entity is: " << rel_entities[j].local_offset());
          if ( !(rel_entities[j] == e_from && rel_ordinals[j] == static_cast<ConnectivityOrdinal>(local_id) ) )
          {
            impl::get_part_ordinals_to_induce_on_lower_ranks(*this, rel_entities[j], e_to_entity_rank, keep);
          }
        }
      }
    }

    // Find the relation this is being deleted and add the parts that are
    // induced from that relation (and that are not in 'keep') to 'del'
    {
      size_t num_rels = num_connectivity(e_from, e_to_entity_rank);
      Entity const *rel_entities = begin(e_from, e_to_entity_rank);
      ConnectivityOrdinal const *rel_ordinals = begin_ordinals(e_from, e_to_entity_rank);
      for (size_t j = 0; j < num_rels; ++j)
      {
        if ( rel_entities[j] == e_to && rel_ordinals[j] == static_cast<ConnectivityOrdinal>(local_id) )
        {
          impl::get_part_ordinals_to_induce_on_lower_ranks_except_for_omits(*this, e_from, keep, e_to_entity_rank, del);
          break; // at most 1 relation can match our specification
        }
      }
    }

    if ( !del.empty() ) {
      internal_change_entity_parts( e_to , empty, del, scratchOrdinalVec, scratchSpace);
    }
  }

  //delete relations from the entities
  bool caused_change_fwd = bucket(e_from).destroy_relation(e_from, e_to, local_id);

  if (caused_change_fwd && bucket(e_from).owned() && (entity_rank(e_from) > entity_rank(e_to)) ) {
    --m_closure_count[e_to.local_offset()];
  }


  // Relationships should always be symmetrical
  if ( caused_change_fwd ) {
    bool caused_change = bucket(e_to).destroy_relation(e_to, e_from, local_id);
    if (caused_change && bucket(e_to).owned() && (entity_rank(e_to) > entity_rank(e_from)) ) {
      --m_closure_count[e_from.local_offset()];
    }
  }

  // It is critical that the modification be done AFTER the relations are
  // changed so that the propagation can happen correctly.
  if ( caused_change_fwd ) {
    this->mark_entity_and_upward_related_entities_as_modified(e_to);
    this->mark_entity_and_upward_related_entities_as_modified(e_from);
  }

  m_check_invalid_rels = true;

  return caused_change_fwd;
}

//----------------------------------------------------------------------
//ParallelVerify
//----------------------------------------------------------------------

bool BulkData::is_entity_in_sharing_comm_map(stk::mesh::Entity entity)
{
    EntityKey entityKey = this->entity_key(entity);
    bool is_entity_in_shared_comm_map = !this->internal_entity_comm_map(entityKey, this->shared_ghosting()).empty();
    return is_entity_in_shared_comm_map;
}

void BulkData::erase_sharing_info_using_key(EntityKey key, stk::mesh::BulkData::GhostingId ghostingId)
{
    this->entity_comm_map_erase(key, *this->ghostings()[ghostingId]);
}

void BulkData::add_sharing_info(stk::mesh::Entity entity, stk::mesh::BulkData::GhostingId ghostingId, int sharingProc)
{
    this->entity_comm_map_insert(entity, stk::mesh::EntityCommInfo(ghostingId, sharingProc));
}

void BulkData::get_entities_that_have_sharing(std::vector<stk::mesh::Entity> &entitiesThatHaveSharingInfo,
        stk::mesh::EntityToDependentProcessorsMap &entityKeySharing)
{
    if(parallel_size() > 1)
    {
        // this communicates states of the entities to all procs so that entity states are consistent
        stk::mesh::EntityVector entitiesNoLongerShared;
        m_meshModification.internal_resolve_shared_modify_delete(entitiesNoLongerShared);
    }

    int myProcId = this->parallel_rank();
    stk::CommSparse commStage1(this->parallel());

    size_t numEntitiesThatHaveSharingInfo = 0;

    for(int phase = 0; phase < 2; ++phase)
    {
        if(phase == 1)
        {
            entitiesThatHaveSharingInfo.resize(numEntitiesThatHaveSharingInfo);
            numEntitiesThatHaveSharingInfo = 0;
        }

        for(stk::mesh::EntityRank irank = stk::topology::NODE_RANK; irank <= stk::topology::FACE_RANK; ++irank)
        {
            stk::mesh::BucketVector buckets_of_rank = this->buckets(irank);
            for(size_t bucket_i = 0; bucket_i != buckets_of_rank.size(); ++bucket_i)
            {
                stk::mesh::Bucket & bucket = *buckets_of_rank[bucket_i];
                for(size_t entity_i = 0; entity_i != bucket.size(); ++entity_i)
                {
                    stk::mesh::Entity entity = bucket[entity_i];
                    if((is_valid(entity) && this->state(entity) == stk::mesh::Modified))
                    {
                        if(phase == 0 && this->is_entity_in_sharing_comm_map(entity))
                        {
                            numEntitiesThatHaveSharingInfo++;
                        }
                        else if(phase == 1 && this->is_entity_in_sharing_comm_map(entity))
                        {
                            entitiesThatHaveSharingInfo[numEntitiesThatHaveSharingInfo] = entity;
                            numEntitiesThatHaveSharingInfo++;
                        }

                        int procThatOwnsEntity = this->parallel_owner_rank(entity);
                        bool anotherProcOwnsThisEntity = procThatOwnsEntity != myProcId;
                        bool entityIsNotGhosted = this->owned_closure(entity);

                        if(anotherProcOwnsThisEntity && entityIsNotGhosted)
                        {
                            stk::mesh::EntityKey entityKey = this->entity_key(entity);
                            commStage1.send_buffer(procThatOwnsEntity).pack<stk::mesh::EntityKey>(entityKey);
                        }
                    }
                }
            }
        }

        if(phase == 0)
        {
            commStage1.allocate_buffers();
        }
        else
        {
            commStage1.communicate();
        }
    }

    for(int procIndex = 0; procIndex < this->parallel_size(); procIndex++)
    {
        if(myProcId == procIndex)
            continue;
        stk::CommBuffer & dataFromAnotherProc = commStage1.recv_buffer(procIndex);
        EntityKey key;
        int sharingProc = procIndex;
        while(dataFromAnotherProc.remaining())
        {
            dataFromAnotherProc.unpack<stk::mesh::EntityKey>(key);
            ThrowAssertMsg(this->is_valid(this->get_entity(key)) && this->parallel_owner_rank(this->get_entity(key)) == myProcId, "Entitykey " << key << " is not owned by receiving processor " << myProcId);
            entityKeySharing[key].insert(sharingProc);
        }
    }
}

void extractEntityToMapInfoIntoVectorOfEntityKeyAndProcPairs(const int myProcId, stk::mesh::EntityToDependentProcessorsMap &entityKeySharing, std::vector<std::pair<stk::mesh::EntityKey, int> >& sharedEntities)
{
    stk::mesh::EntityToDependentProcessorsMap::iterator iter = entityKeySharing.begin();
    for(; iter != entityKeySharing.end(); iter++)
    {
        std::vector<int> sharingProcs(iter->second.begin(), iter->second.end());
        iter->second.insert(myProcId);
        for(size_t j = 0; j < sharingProcs.size(); j++)
        {
            sharedEntities.emplace_back(iter->first, sharingProcs[j]);
        }
    }
}

void BulkData::get_locally_modified_shared_entities(stk::mesh::EntityToDependentProcessorsMap &entityKeySharing, std::vector<std::pair<stk::mesh::EntityKey, int> >& sharedEntities)
{
    extractEntityToMapInfoIntoVectorOfEntityKeyAndProcPairs(this->parallel_rank(), entityKeySharing, sharedEntities);

    stk::CommSparse commStage2(this->parallel());
    impl::communicateSharingInfoToProcsThatShareEntity(this->parallel_size(), this->parallel_rank(), commStage2, entityKeySharing);
    impl::unpackCommunicationsAndStoreSharedEntityToProcPair(this->parallel_size(), this->parallel_rank(), commStage2, sharedEntities);
}

void BulkData::erase_all_sharing_for_invalid_entities_on_comm_map()
{
    for(size_t i=0; i<this->internal_comm_list().size(); ++i)
    {
        stk::mesh::EntityKey key = this->internal_comm_list()[i].key;
        stk::mesh::Entity entity = this->get_entity(key);
        if( !this->is_valid(entity) )
        {
            this->erase_sharing_info_using_key(key, stk::mesh::BulkData::SHARED);
        }
    }
}

void BulkData::fill_entities_that_have_lost_sharing_info(const std::vector<std::pair<stk::mesh::EntityKey, int> > &sharedEntities,
        const std::vector<stk::mesh::Entity>& entitiesThatUsedToHaveSharingInfoBeforeCEO, std::vector<stk::mesh::Entity>& modifiedEntitiesForWhichCommMapsNeedUpdating)
{
    std::set<stk::mesh::EntityKey> keysThatNeedToHaveCorrectSharingInfo;
    for (size_t i=0;i<sharedEntities.size();i++)
    {
        keysThatNeedToHaveCorrectSharingInfo.insert(sharedEntities[i].first);
    }

    for (size_t i=0;i<entitiesThatUsedToHaveSharingInfoBeforeCEO.size();i++)
    {
        stk::mesh::EntityKey entityKey = this->entity_key(entitiesThatUsedToHaveSharingInfoBeforeCEO[i]);
        const bool keyNotInSharedKeysList = keysThatNeedToHaveCorrectSharingInfo.find(entityKey) == keysThatNeedToHaveCorrectSharingInfo.end();
        if ( keyNotInSharedKeysList )
        {
            modifiedEntitiesForWhichCommMapsNeedUpdating.push_back(entitiesThatUsedToHaveSharingInfoBeforeCEO[i]);
            this->erase_sharing_info_using_key(entityKey, stk::mesh::BulkData::SHARED);
        }
    }
}

void BulkData::fill_modified_entities_and_add_sharing_comm_map_info_for_shared_entities(const std::vector<std::pair<stk::mesh::EntityKey, int> > &sharedEntities,
        const std::vector<stk::mesh::Entity>& entitiesThatUsedToHaveSharingInfoBeforeCEO, std::vector<stk::mesh::Entity>& modifiedEntitiesForWhichCommMapsNeedUpdating)
{
    erase_all_sharing_for_invalid_entities_on_comm_map();
    fill_entities_that_have_lost_sharing_info(sharedEntities, entitiesThatUsedToHaveSharingInfoBeforeCEO, modifiedEntitiesForWhichCommMapsNeedUpdating);

    for (size_t i=0;i<sharedEntities.size();i++)
    {
        stk::mesh::EntityKey key = sharedEntities[i].first;
        this->erase_sharing_info_using_key(key, stk::mesh::BulkData::SHARED);
    }

    for(size_t i = 0; i < sharedEntities.size(); i++)
    {
        stk::mesh::Entity entity = this->get_entity(sharedEntities[i].first);
        this->add_sharing_info(entity, stk::mesh::BulkData::SHARED, sharedEntities[i].second);
        modifiedEntitiesForWhichCommMapsNeedUpdating.push_back(entity);
    }
}

void BulkData::resolve_entity_ownership_and_part_membership_and_comm_list(std::vector<stk::mesh::Entity>& modifiedEntities)
{
    this->resolve_ownership_of_modified_entities(modifiedEntities);
    this->move_entities_to_proper_part_ownership(modifiedEntities);
    this->update_comm_list_based_on_changes_in_comm_map();
    this->add_comm_list_entries_for_entities(modifiedEntities);
}

void BulkData::update_sharing_after_change_entity_owner()
{
    std::vector<stk::mesh::Entity> entitiesThatHaveSharingInfo;
    stk::mesh::EntityToDependentProcessorsMap ownerReceiviesInfoOnOtherProcessorsThatShareEntitiesThisProcOwns;

    get_entities_that_have_sharing(entitiesThatHaveSharingInfo, ownerReceiviesInfoOnOtherProcessorsThatShareEntitiesThisProcOwns);

    std::vector<std::pair<stk::mesh::EntityKey, int> > sharedEntities;
    get_locally_modified_shared_entities(ownerReceiviesInfoOnOtherProcessorsThatShareEntitiesThisProcOwns, sharedEntities);

    std::vector<stk::mesh::Entity> modifiedEntities;
    fill_modified_entities_and_add_sharing_comm_map_info_for_shared_entities(sharedEntities, entitiesThatHaveSharingInfo, modifiedEntities);

    resolve_entity_ownership_and_part_membership_and_comm_list(modifiedEntities);
}



void BulkData::internal_change_entity_owner( const std::vector<EntityProc> & arg_change,
                                             impl::MeshModification::modification_optimization mod_optimization )
{
  require_ok_to_modify();
  m_modSummary.track_change_entity_owner(arg_change);

  const MetaData  & meta = m_mesh_meta_data ;
  const int       p_rank = parallel_rank() ;
  const int       p_size = parallel_size() ;
  ParallelMachine p_comm = parallel() ;

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
      StoreEntityInSet store_entity;
      for ( std::vector<EntityProc>::const_iterator i = ghosted_change.begin() ; i != ghosted_change.end() ; ++i) {
          impl::VisitAuraClosureGeneral(*this,i->first,store_entity,only_visit_ghosts_once);
      }
      for ( std::vector<EntityProc>::const_iterator i = shared_change.begin() ; i != shared_change.end() ; ++i) {
          impl::VisitAuraClosureGeneral(*this,i->first,store_entity,only_visit_ghosts_once);
      }
      for ( std::set<EntityProc,EntityLess>::const_iterator i = send_closure.begin() ; i != send_closure.end() ; ++i) {
          impl::VisitAuraClosureGeneral(*this,i->first,store_entity,only_visit_ghosts_once);
      }

    std::set<Entity> & modified_ghosts = store_entity.entity_set;

    std::set<Entity>::iterator iter = modified_ghosts.begin();
    std::vector<Entity> keysThatAreBothSharedAndCustomGhosted;

    for (;iter!=modified_ghosts.end();++iter)
    {
        if ( in_shared(*iter) )
        {
            keysThatAreBothSharedAndCustomGhosted.push_back(*iter);
        }
    }

    for(size_t i=0;i<keysThatAreBothSharedAndCustomGhosted.size();++i)
    {
        modified_ghosts.erase(keysThatAreBothSharedAndCustomGhosted[i]);
        entity_comm_map_clear_ghosting(entity_key(keysThatAreBothSharedAndCustomGhosted[i]));
    }

    // The ghosted change list will become invalid
    ghosted_change.clear();

    std::vector<EntityProc> empty_add ;
    std::vector<Entity>  remove_modified_ghosts( modified_ghosts.begin(), modified_ghosts.end() );

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
    ConstPartVector owned;
    owned.push_back(& meta.locally_owned_part());

    for ( std::vector<EntityProc>::iterator
          i = local_change.begin() ; i != local_change.end() ; ++i ) {
      // Giving ownership, change the parts first and then
      // the owner rank to pass the ownership test.
      Entity entity = i->first;

      internal_verify_and_change_entity_parts( entity , ConstPartVector() , owned );

      internal_set_owner(entity, i->second);
    }

    for ( std::vector<EntityProc>::iterator
          i = shared_change.begin() ; i != shared_change.end() ; ++i ) {
      Entity entity = i->first;
      internal_set_owner(entity, i->second);
      if ( p_rank == i->second ) { // I received ownership
          internal_verify_and_change_entity_parts( entity , owned , ConstPartVector() );
      }
    }
  }


  //------------------------------
  // Send entities, along with their closure, to the new owner processes
  {
    std::ostringstream error_msg ;
    int error_count = 0 ;

    stk::CommSparse comm( p_comm );

    EntityVector unique_list_of_send_closure;
    unique_list_of_send_closure.reserve(send_closure.size());

    for ( std::set<EntityProc,EntityLess>::iterator
          i = send_closure.begin() ; i != send_closure.end() ; ++i ) {
      CommBuffer & buffer = comm.send_buffer( i->second );
      Entity entity = i->first;
      pack_entity_info(*this, buffer , entity );
      if (!is_communicated_with_proc(entity, i->second) ||
          std::binary_search(local_change.begin(), local_change.end(), *i, EntityLess(*this))) {
        buffer.pack<int>(1);
        pack_field_values(*this, buffer , entity );
      }
      else {
        buffer.pack<int>(0);
      }
      pack_sideset_info(*this, buffer , entity );

      if (unique_list_of_send_closure.empty() || entity_key(unique_list_of_send_closure.back()) != entity_key(entity)) {
        unique_list_of_send_closure.push_back(entity);
      }
    }

    comm.allocate_buffers();

    for ( std::set<EntityProc,EntityLess>::iterator
          i = send_closure.begin() ; i != send_closure.end() ; ++i ) {
      CommBuffer & buffer = comm.send_buffer( i->second );
      Entity entity = i->first;
      pack_entity_info(*this, buffer , entity );
      if (!is_communicated_with_proc(entity, i->second) ||
          std::binary_search(local_change.begin(), local_change.end(), *i, EntityLess(*this))) {
        buffer.pack<int>(1);
        pack_field_values(*this, buffer , entity );
      }
      else {
        buffer.pack<int>(0);
      }
      pack_sideset_info(*this, buffer , entity );
    }

    comm.communicate();

    for ( std::set<EntityProc,EntityLess>::iterator
          i = send_closure.begin() ; i != send_closure.end() ; ++i ) {
      Entity entity = i->first;
      remove_element_entries_from_sidesets(*this, entity);
    }

    OrdinalVector partOrdinals;
    OrdinalVector scratchOrdinalVec, scratchSpace;
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
        partOrdinals.clear();
        for(const stk::mesh::Part* part : parts) {
            partOrdinals.push_back(part->mesh_meta_data_ordinal());
        }

        internal_change_entity_parts( entity , partOrdinals , OrdinalVector(), scratchOrdinalVec, scratchSpace );

        log_created_parallel_copy( entity );

        internal_set_owner(entity, owner);

        internal_declare_relation( entity , relations );

        int shouldUnpackFieldValues = 0;
        buf.unpack<int>(shouldUnpackFieldValues);
        if ( shouldUnpackFieldValues==1 ) {
          if ( ! unpack_field_values(*this, buf , entity , error_msg ) ) {
            ++error_count ;
          }
        }

        unpack_sideset_info( buf, *this, entity);
      }
    }

#ifndef NDEBUG
    all_reduce( p_comm , ReduceSum<1>( & error_count ) );
    ThrowErrorMsgIf( error_count, error_msg.str() );
#endif

    // Any entity that I sent and is not in an owned closure is deleted.
    // The owned closure will be effected by received entities, so can
    // only clean up after the newly owned entities have been received.
    // Destroy backwards so as not to invalidate closures in the process.

    {
        for ( EntityVector::reverse_iterator i = unique_list_of_send_closure.rbegin() ; i != unique_list_of_send_closure.rend() ; ++i) {
            stk::mesh::Entity entity = *i;
            if ( ! this->owned_closure(entity) ) {
                ThrowRequireMsg( internal_destroy_entity( entity ), "Failed to destroy entity " << identifier(entity) );
            }
        }
    }

    send_closure.clear(); // Has been invalidated
  }
}

Ghosting & BulkData::create_ghosting( const std::string & name )
{
    return internal_create_ghosting(name);
}

Ghosting & BulkData::internal_create_ghosting( const std::string & name )
{
  require_ok_to_modify();

#ifndef NDEBUG
  // Verify name is the same on all processors,
  // if not then throw an exception on all processors.
  if (parallel_size() > 1) {
    CommBroadcast bc( parallel() , 0 );

    if ( bc.parallel_rank() == 0 ) {
      bc.send_buffer().skip<char>( name.size() + 1 );
    }

    bc.allocate_buffer();

    if ( bc.parallel_rank() == 0 ) {
      bc.send_buffer().pack<char>( name.c_str() , name.size() + 1 );
    }

    bc.communicate();

    const char * const bc_name =
      reinterpret_cast<const char *>( bc.recv_buffer().buffer() );

    int error = 0 != std::strcmp( bc_name , name.c_str() );

    all_reduce( parallel() , ReduceMax<1>( & error ) );

    ThrowErrorMsgIf( error, "Parallel name inconsistency");
  }
#endif

  for(Ghosting* ghosting : ghostings()) {
    if (ghosting->name() == name) {
      return *ghosting;
    }
  }

  Ghosting * const g =
    new Ghosting( *this , name , m_ghosting.size() );

  m_ghosting.push_back( g );

  if (m_ghost_parts.size() == 0) {
    ThrowRequireMsg(equal_case(std::string("shared"), name), "Expect shared to be the first ghosting created.");
    m_ghost_parts.push_back(&mesh_meta_data().globally_shared_part());
  }
  else if (m_ghost_parts.size() == 1) {
    ThrowRequireMsg(equal_case(std::string("shared_aura"), name), "Expect aura to be the second ghosting created.");
    Part & aura_part = mesh_meta_data().aura_part();
    aura_part.entity_membership_is_parallel_consistent(false);
    m_ghost_parts.push_back(&aura_part);
  }
  else {
    std::ostringstream oss;
    oss << "custom_ghosting_" << m_ghost_parts.size();
    std::string ghostPartName = stk::mesh::impl::convert_to_internal_name(oss.str());
    Part& ghost_part = mesh_meta_data().declare_part(ghostPartName);
    ghost_part.entity_membership_is_parallel_consistent(false);
    m_ghost_parts.push_back(&ghost_part);
  }

  ThrowRequireMsg(m_ghost_parts.size() == m_ghosting.size(), "m_ghost_parts.size()="<<m_ghost_parts.size()<<", must be same as m_ghosting.size()="<<m_ghosting.size());

  return *g ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void BulkData::destroy_ghosting( Ghosting& ghost_layer )
{
  std::vector<EntityKey> receive_list;
  ghost_layer.receive_list(receive_list);
  internal_verify_inputs_and_change_ghosting(ghost_layer, std::vector<stk::mesh::EntityProc>(), receive_list);
}

//----------------------------------------------------------------------

void BulkData::destroy_all_ghosting()
{
  require_ok_to_modify();

  // Iterate backwards so as not to invalidate a closure.

  for ( EntityCommListInfoVector::reverse_iterator
        i =  m_entity_comm_list.rbegin() ;
        i != m_entity_comm_list.rend() ; ++i) {

    if ( in_receive_ghost( i->key ) ) {
      entity_comm_map_clear_ghosting( i->key );
      bool entity_is_not_shared = !bucket(i->entity).shared();
      internal_destroy_entity_with_notification( i->entity );
      if(entity_is_not_shared)
      {
        i->key = EntityKey();
        i->entity_comm = NULL;
      }
    }
    else {
      entity_comm_map_clear_ghosting(i->key);
      if ( internal_entity_comm_map(i->key).empty() ) {
        i->key = EntityKey();
        i->entity_comm = NULL;
      }
    }
  }

  delete_unneeded_entries_from_the_comm_list();
}

//----------------------------------------------------------------------

void BulkData::change_ghosting(
  Ghosting & ghosts ,
  const std::vector<EntityProc> & add_send ,
  const std::vector<EntityKey> & remove_receive )
{
    internal_verify_inputs_and_change_ghosting(ghosts, add_send, remove_receive);
}

void BulkData::verify_and_filter_add_send(Ghosting & ghosts, const std::vector<EntityProc> & add_send, bool &need_to_change_ghosting,
                          bool &add_send_is_owned, std::vector <EntityProc> &filtered_add_send )
{
    filtered_add_send.reserve(add_send.size());

    for ( size_t i = 0; add_send_is_owned && i < add_send.size() ; ++i ) {
      add_send_is_owned = parallel_owner_rank(add_send[i].first) == parallel_rank();
      const bool ghosting_to_myself = parallel_rank() == add_send[i].second;
      const bool already_ghosted_to_proc = in_send_ghost(ghosts, entity_key(add_send[i].first), add_send[i].second);
      const bool need_to_send_ghost = !ghosting_to_myself && !already_ghosted_to_proc;
      if (need_to_send_ghost)
      {
          filtered_add_send.push_back(add_send[i]);
          need_to_change_ghosting = true;
      }
    }
}

void BulkData::verify_remove_receive(Ghosting & ghosts, const std::vector<Entity> & remove_receive, bool &need_to_change_ghosting, bool &remove_receive_are_part_of_this_ghosting)
{
    for ( size_t i = 0; remove_receive_are_part_of_this_ghosting && i < remove_receive.size() ; ++i ) {
      remove_receive_are_part_of_this_ghosting = in_receive_ghost( ghosts , remove_receive[i] );
      need_to_change_ghosting = true;
    }
}

bool BulkData::check_errors_and_determine_if_ghosting_needed_in_parallel(const stk::mesh::Ghosting &ghosts,
                                        bool add_send_is_owned,
                                        bool remove_receive_are_part_of_this_ghosting,
                                        bool need_to_change_ghosting,
                                        const std::vector<EntityProc> & add_send,
                                        const std::vector<Entity> & remove_receive)
{
    const bool ok_mesh  = &ghosts.mesh() == this;
    const bool is_custom_ghost = BulkData::AURA < ghosts.ordinal();
    int ok = ok_mesh && is_custom_ghost && add_send_is_owned && remove_receive_are_part_of_this_ghosting;
    int statuses[2];
    statuses[0] = ok;
    statuses[1] = need_to_change_ghosting ? 0 : 1;

    all_reduce( parallel() , ReduceMin<2>( statuses ) );

    ok = statuses[0];
    bool something_wrong_on_any_proc = (0 == ok);
    if ( something_wrong_on_any_proc ) {
      std::ostringstream msg ;
      msg << "For ghosts " << ghosts.name() << ", " ;
      if ( ! ok_mesh )  { msg << " : Mesh does not own this ghosting" ; }
      if ( ! is_custom_ghost ) { msg << " : Cannot modify this ghosting" ; }
      if ( ! add_send_is_owned ) {
        msg << " : Not owned add {" ;
        for ( std::vector<EntityProc>::const_iterator
              i = add_send.begin() ; i != add_send.end() ; ++i ) {
          if ( parallel_owner_rank(i->first) != parallel_rank() ) {
            msg << " " << identifier(i->first);
          }
        }
        msg << " }" ;
      }
      if ( ! remove_receive_are_part_of_this_ghosting ) {
        msg << " : Not in ghost receive {" ;
        for (stk::mesh::Entity entity : remove_receive) {
          if ( ! in_receive_ghost( ghosts , entity ) ) {
            msg << " " << identifier(entity);
          }
        }
      }

      ThrowErrorMsg( msg.str() );
    }

    bool anyProcsHaveNewGhosts = statuses[1] == 0;
    return anyProcsHaveNewGhosts;
}

bool BulkData::inputs_ok_and_need_ghosting(Ghosting & ghosts ,
                             const std::vector<EntityProc> & add_send ,
                             const std::vector<Entity> & remove_receive,
                             std::vector<EntityProc> &filtered_add_send)
{
    bool add_send_is_owned    = true ;
    bool need_to_change_ghosting = false;

    verify_and_filter_add_send(ghosts, add_send, need_to_change_ghosting, add_send_is_owned, filtered_add_send );

    bool remove_receive_are_part_of_this_ghosting = true;
    verify_remove_receive(ghosts, remove_receive, need_to_change_ghosting, remove_receive_are_part_of_this_ghosting);

    bool anyProcsHaveNewGhosts = check_errors_and_determine_if_ghosting_needed_in_parallel(ghosts, add_send_is_owned,
        remove_receive_are_part_of_this_ghosting, need_to_change_ghosting, add_send, remove_receive);

    return anyProcsHaveNewGhosts;
}

void BulkData::batch_add_to_ghosting(Ghosting &ghosting, const EntityProcVec &entitiesAndDestinationProcs)
{
    std::vector<stk::mesh::EntityProc> filtered_add_send;
    std::vector<stk::mesh::Entity> empty_vector;

    if (inputs_ok_and_need_ghosting(ghosting , entitiesAndDestinationProcs , empty_vector, filtered_add_send))
    {
        internal_batch_add_to_ghosting(ghosting, filtered_add_send);
    }
}

void BulkData::internal_batch_add_to_ghosting(Ghosting &ghosting, const EntityProcVec &entitiesAndDestinationProcs)
{
    bool starting_modification = modification_begin();
    ThrowRequireMsg(starting_modification, "ERROR: BulkData already being modified,\n"
                    <<"BulkData::batch_add_to_ghosting(...) can not be called within an outer modification scope.");
    internal_add_to_ghosting( ghosting , entitiesAndDestinationProcs );
    internal_modification_end_for_change_ghosting();
}

void BulkData::internal_verify_inputs_and_change_ghosting(
  Ghosting & ghosts ,
  const std::vector<EntityProc> & add_send ,
  const std::vector<EntityKey> & remove_receive )
{
  require_ok_to_modify();

  std::vector<EntityProc> filtered_add_send;
  std::vector<Entity> remove_receive_entities(remove_receive.size());
  for(unsigned i=0; i<remove_receive.size(); ++i) {
    remove_receive_entities[i] = get_entity(remove_receive[i]);
  }
  bool needToDoGhosting = inputs_ok_and_need_ghosting(ghosts, add_send , remove_receive_entities, filtered_add_send);

  if(needToDoGhosting)
  {
    internal_change_ghosting( ghosts , filtered_add_send , remove_receive_entities );
  }
}

//----------------------------------------------------------------------

void BulkData::ghost_entities_and_fields(Ghosting & ghosting, const std::set<EntityProc , EntityLess>& sendGhosts, bool isFullRegen)
{
    //------------------------------------
    // Push newly ghosted entities to the receivers and update the comm list.

    const size_t record_entity_comm_size_before_changing_it = m_entity_comm_list.size();
    const int p_size = parallel_size() ;

    stk::CommSparse commSparse( parallel() );
    for ( int phase = 0; phase < 2; ++phase ) {
      for ( std::set< EntityProc , EntityLess >::iterator
              j = sendGhosts.begin(); j != sendGhosts.end() ; ++j ) {

        Entity entity = j->first;
        const int proc = j->second;

        if ( isFullRegen || !in_ghost(ghosting , entity, proc) ) {
          // Not already being sent , must send it.
          CommBuffer & buf = commSparse.send_buffer( proc );
          buf.pack<unsigned>( entity_rank(entity) );
          pack_entity_info(*this, buf , entity );
          pack_field_values(*this, buf , entity );

          if (phase == 1) {
            std::pair<EntityComm*,bool> result = entity_comm_map_insert(entity, EntityCommInfo(ghosting.ordinal(), proc));
            const EntityComm* entity_comm = result.first;
            EntityCommListInfo comm_info = {entity_key(entity), entity,
                                            entity_comm};
            m_entity_comm_list.push_back( comm_info );
          }
        }
      }

      if (phase == 0) {
        commSparse.allocate_buffers();
      }
      else {
        commSparse.communicate();
      }
    }

    std::ostringstream error_msg ;
    int error_count = 0 ;
    OrdinalVector ordinal_scratch, empty, partOrdinals, scratchSpace;
    PartVector parts ;
    std::vector<Relation> relations ;

    const MetaData & meta = m_mesh_meta_data ;
    const unsigned rank_count = meta.entity_rank_count();

    // Unpacking must proceed in entity-rank order so that higher ranking
    // entities that have relations to lower ranking entities will have
    // the lower ranking entities unpacked first.  The higher and lower
    // ranking entities may be owned by different processes,
    // as such unpacking must be performed in rank order.

    for ( unsigned rank = 0 ; rank < rank_count ; ++rank ) {
      for ( int p = 0 ; p < p_size ; ++p ) {
        CommBuffer & buf = commSparse.recv_buffer(p);
        while ( buf.remaining() ) {
          // Only unpack if of the current entity rank.
          // If not the current entity rank, break the iteration
          // until a subsequent entity rank iteration.
          {
            unsigned this_rank = ~0u ;
            buf.peek<unsigned>( this_rank );

            if ( this_rank != rank ) break ;

            buf.unpack<unsigned>( this_rank );
          }

          parts.clear();
          relations.clear();
          EntityKey key ;
          int owner = ~0u ;

          unpack_entity_info( buf, *this, key, owner, parts, relations );

          if (owner != this->parallel_rank()) {
            // We will also add the entity to the part corresponding to the 'ghosts' ghosting.
            stk::mesh::Part& ghost_part = *m_ghost_parts[ghosting.ordinal()];
            insert( parts, ghost_part );
          }

          GhostReuseMap::iterator f_itr = m_ghost_reuse_map.find(key);
          const size_t use_this_offset = f_itr == m_ghost_reuse_map.end() ? 0 : f_itr->second;
          if (use_this_offset != 0) {
            m_ghost_reuse_map.erase(f_itr);
          }

          std::pair<Entity ,bool> result = internal_get_or_create_entity_with_notification( key, use_this_offset );

          Entity entity = result.first;
          const bool created   = result.second ;

          require_entity_owner( entity , owner );

          partOrdinals.clear();
          for(const stk::mesh::Part* part : parts) {
              partOrdinals.push_back(part->mesh_meta_data_ordinal());
          }

          internal_change_entity_parts( entity , partOrdinals , empty, ordinal_scratch, scratchSpace );

          if ( created ) {
            log_created_parallel_copy( entity );
          }

          const EntityCommInfo tmp( ghosting.ordinal() , owner );

          std::pair<EntityComm*, bool> insertResult = entity_comm_map_insert(entity, tmp);
          if ( insertResult.second ) {
            const EntityComm* entity_comm = insertResult.first;
            EntityCommListInfo comm_info = {entity_key(entity), entity,
                                            entity_comm};
            m_entity_comm_list.push_back( comm_info );
          }

          //now, change owner. (needed to wait until comm-info was created)
          internal_set_owner( entity, owner);

          internal_declare_relation( entity , relations, ordinal_scratch);

          if ( ! unpack_field_values(*this, buf , entity , error_msg ) ) {
            ++error_count ;
          }
        }
      }
    }

#ifndef NDEBUG
    if (parallel_size() > 1) {
      all_reduce( parallel() , ReduceSum<1>( & error_count ) );
    }
    ThrowErrorMsgIf( error_count, error_msg.str() );
#endif

    if ( record_entity_comm_size_before_changing_it < m_entity_comm_list.size() ) {
      // Added new ghosting entities to the list,
      // must now sort and merge.

      EntityCommListInfoVector::iterator i = m_entity_comm_list.begin();
      i += record_entity_comm_size_before_changing_it ;
      std::sort( i , m_entity_comm_list.end() );
      std::inplace_merge( m_entity_comm_list.begin() , i ,
                          m_entity_comm_list.end() );
      m_entity_comm_list.erase( std::unique( m_entity_comm_list.begin() , m_entity_comm_list.end() ) ,
                                m_entity_comm_list.end() );
    }
}

void BulkData::conditionally_add_entity_to_ghosting_set(const stk::mesh::Ghosting &ghosting, stk::mesh::Entity entity, int toProc, std::set <stk::mesh::EntityProc , stk::mesh::EntityLess > &entitiesWithClosure)
{
    const bool notOwnedByRecvGhostProc = toProc != parallel_owner_rank(entity);
    const bool entityIsShared = in_shared( entity_key(entity) , toProc );
    const bool alreadyGhostedToProc = in_send_ghost(ghosting, entity_key(entity), toProc);
    const bool ghostingIsAura = ghosting.ordinal() == BulkData::AURA;

    bool shouldAddToGhostingSet = notOwnedByRecvGhostProc && !alreadyGhostedToProc;
    if (ghostingIsAura && entityIsShared) {
        shouldAddToGhostingSet = false;
    }

    if (shouldAddToGhostingSet)
    {
        entitiesWithClosure.insert(EntityProc(entity , toProc));
    }
}

void BulkData::add_closure_entities(const stk::mesh::Ghosting &ghosting, const stk::mesh::EntityProcVec& entities, std::set <stk::mesh::EntityProc , stk::mesh::EntityLess > &entitiesWithClosure)
{
    for ( std::vector< EntityProc >::const_iterator
          i = entities.begin() ; i != entities.end() ; ++i )
    {
        if(is_valid(i->first))
        {
            conditionally_add_entity_to_ghosting_set(ghosting, i->first, i->second, entitiesWithClosure);

            stk::mesh::EntityRank entityRank = entity_rank(i->first);
            for(stk::mesh::EntityRank irank = stk::topology::NODE_RANK; irank < entityRank; ++irank)
            {
                unsigned numEntities = num_connectivity(i->first, irank);
                const stk::mesh::Entity* connectedEntities = begin(i->first, irank);
                for(unsigned entityIndex = 0; entityIndex < numEntities; ++entityIndex)
                {
                    if(is_valid(connectedEntities[entityIndex]))
                    {
                        conditionally_add_entity_to_ghosting_set(ghosting, connectedEntities[entityIndex], i->second, entitiesWithClosure);
                    }
                }
            }
        }
    }
}

void BulkData::internal_add_to_ghosting(
        Ghosting &ghosting,
        const std::vector<EntityProc> &add_send)
{
    m_modSummary.track_add_to_ghosting(ghosting, add_send);

    //------------------------------------
    // Copy ghosting lists into more efficiently edited container.
    // The send and receive lists must be in entity rank-order.
    std::set<EntityProc , EntityLess> entitiesToGhostOntoOtherProcessors(EntityLess(*this));

    add_closure_entities(ghosting, add_send, entitiesToGhostOntoOtherProcessors);

    stk::mesh::impl::move_unowned_entities_for_owner_to_ghost(*this, entitiesToGhostOntoOtherProcessors);

    ghost_entities_and_fields(ghosting, entitiesToGhostOntoOtherProcessors);
}

void BulkData::generate_ghosting_receive_list(const stk::mesh::Ghosting &ghosting, const std::vector <Entity> &remove_receive,
    std::vector<Entity> &recvGhosts, std::vector<bool>& ghostStatus)
{
    // Remove any entities that are in the remove list.

    for ( Entity rmEntity : remove_receive) {
      if (is_valid(rmEntity) && in_receive_ghost(ghosting, rmEntity)) {
        ghostStatus[rmEntity.local_offset()] = true;
      }
    }

    // Iterate over all entities with communication information, adding
    // the entity if it's a ghost on this process. recvGhosts will contain
    // all received-ghosts on this process by the end of the loop.
    for ( EntityCommListInfoVector::const_iterator
          i = internal_comm_list().begin() ; i != internal_comm_list().end() ; ++i ) {
      const bool inRemoveReceive = ghostStatus[i->entity.local_offset()];
      if ( is_valid(i->entity) && !inRemoveReceive && in_receive_ghost(ghosting, i->entity) ) {
        recvGhosts.push_back(i->entity);
        ghostStatus[i->entity.local_offset()] = true;
      }
      else {
        ghostStatus[i->entity.local_offset()] = false;
      }
    }

    //Add back in the closure-entities of each recv-ghost, if those closure-entities were
    //removed due to being in the remove_receive list
    //
    impl::OnlyRecvGhosts org(*this,ghosting,ghostStatus);
    impl::VecPushBack vpb(recvGhosts, ghostStatus);

    unsigned len = recvGhosts.size();
    for (unsigned ii=0; ii<len; ++ii) {
      Entity e = recvGhosts[ii];
      const unsigned erank = entity_rank(e);

      for (EntityRank irank = stk::topology::BEGIN_RANK; irank < erank; ++irank) {
        Entity const *rels_i = begin(e, irank);
        Entity const *rels_e = end(e, irank);
        for (; rels_i != rels_e; ++rels_i)
        {
          if (erank > stk::topology::ELEM_RANK) {
            impl::VisitClosureGeneral(*this, *rels_i, vpb, org);
          }
          else {
            if ( is_valid(*rels_i) &&
                 in_receive_ghost( ghosting , *rels_i ) && !ghostStatus[(*rels_i).local_offset()] )
            {
              recvGhosts.push_back(*rels_i);
              ghostStatus[(*rels_i).local_offset()] = true;
            }
          }
        }
      }
    }
}

void BulkData::delete_unneeded_entries_from_the_comm_list()
{
    EntityCommListInfoVector::iterator i =
      std::remove_if( m_entity_comm_list.begin() ,
                      m_entity_comm_list.end() , IsInvalid() );
    m_entity_comm_list.erase( i , m_entity_comm_list.end() );
}

void BulkData::internal_change_ghosting( Ghosting & ghosting,
                                         EntityProcMapping& entityProcMapping)
{
  //------------------------------------
  // Copy ghosting lists into more efficiently edited container.
  // The send and receive lists must be in entity rank-order.

  std::vector<EntityProc> add_send;
  entityProcMapping.fill_vec(add_send);

  //------------------------------------
  // Add the specified entities and their closure to entityProcMapping

  impl::StoreInEntityProcMapping siepm(*this, entityProcMapping);
  EntityProcMapping epm(this->get_size_of_entity_index_space());
  impl::OnlyGhostsEPM og(*this, epm);
  for ( const EntityProc& entityProc : add_send ) {
      og.proc = entityProc.second;
      siepm.proc = entityProc.second;
      impl::VisitClosureGeneral(*this,entityProc.first,siepm,og);
  }

  entityProcMapping.fill_vec(add_send);

  // Synchronize the send and receive list.
  // If the send list contains a not-owned entity
  // inform the owner and receiver to add that entity
  // to their ghost send and receive lists.

  std::vector<bool> ghostStatus(get_size_of_entity_index_space(), false);

  stk::mesh::impl::comm_sync_aura_send_recv(*this, add_send,
                                            entityProcMapping, ghostStatus );

  //------------------------------------
  // Remove the ghost entities that will not remain.
  // If the last reference to the receive ghost entity then delete it.

  OrdinalVector addParts;
  OrdinalVector removeParts(1, m_ghost_parts[ghosting.ordinal()]->mesh_meta_data_ordinal());
  OrdinalVector scratchOrdinalVec, scratchSpace;
  bool removed = false ;

  std::vector<EntityCommInfo> comm_ghost ;
  for ( EntityCommListInfoVector::reverse_iterator
        i = m_entity_comm_list.rbegin() ; i != m_entity_comm_list.rend() ; ++i) {

    if (!i->entity_comm || !i->entity_comm->isGhost) {
      continue;
    }

    const bool is_owner = parallel_owner_rank(i->entity) == parallel_rank() ;
    const bool remove_recv = ( ! is_owner ) &&
                             !ghostStatus[i->entity.local_offset()] && in_receive_ghost(ghosting, i->entity);

    if(is_valid(i->entity))
    {
      if ( is_owner ) {
        // Is owner, potentially removing ghost-sends
        // Have to make a copy

          const PairIterEntityComm ec = internal_entity_comm_map(i->entity, ghosting);
          comm_ghost.assign( ec.first , ec.second );
  
          for ( ; ! comm_ghost.empty() ; comm_ghost.pop_back() ) {
            const EntityCommInfo tmp = comm_ghost.back();

            if (!entityProcMapping.find(i->entity, tmp.proc) ) {
              entity_comm_map_erase(i->key, tmp);
            }
            else {
              entityProcMapping.eraseEntityProc(i->entity, tmp.proc);
            }
          }
      }
      else if ( remove_recv ) {
          entity_comm_map_erase(i->key, ghosting);
          internal_change_entity_parts(i->entity, addParts, removeParts, scratchOrdinalVec, scratchSpace);
      }

      if ( internal_entity_comm_map(i->entity).empty() ) {
        removed = true ;
        i->key = EntityKey(); // No longer communicated
        if ( remove_recv ) {
          ThrowRequireMsg( internal_destroy_entity_with_notification( i->entity, remove_recv ),
                           "P[" << this->parallel_rank() << "]: FAILED attempt to destroy entity: "
                           << entity_key(i->entity) );
        }
      }
    }
  }

  // if an entry in the comm_list has the EntityKey() value, it is invalid,
  // and removed from the comm_list

  if ( removed ) {
    delete_unneeded_entries_from_the_comm_list();
  }

  std::set<EntityProc , EntityLess> sendGhosts(EntityLess(*this));
  entityProcMapping.fill_set(sendGhosts);

  const bool isFullRegen = true;
  ghost_entities_and_fields(ghosting, sendGhosts, isFullRegen);
}

void BulkData::internal_change_ghosting(
  Ghosting & ghosting ,
  const std::vector<EntityProc> & add_send ,
  const std::vector<Entity> & remove_receive,
  bool is_full_regen)
{
  m_modSummary.track_change_ghosting(ghosting, add_send, remove_receive);

  //------------------------------------
  // Copy ghosting lists into more efficiently edited container.
  // The send and receive lists must be in entity rank-order.

  std::set<EntityProc , EntityLess> sendGhosts(EntityLess(*this));
  std::vector<Entity>               recvGhosts;
  std::vector<bool> ghostStatus(get_size_of_entity_index_space(), false);

  //------------------------------------
  // Insert the current ghost receives and then remove from that list.

  // This if-check is an optimization; if doing a full regen then we are
  // removing all ghosting information and recvGhosts should be left empty.
  if ( !is_full_regen ) {
    generate_ghosting_receive_list(ghosting, remove_receive, recvGhosts, ghostStatus);

    //  Initialize sendGhosts from recvGhosts
    stk::mesh::impl::send_entity_keys_to_owners(*this , recvGhosts , sendGhosts);
  }

  //------------------------------------
  // Add the specified entities and their closure to sendGhosts

  impl::StoreInEntityProcSet sieps(*this, sendGhosts);
  impl::OnlyGhosts og(*this);
  for ( const EntityProc& entityProc : add_send ) {
      og.proc = entityProc.second;
      sieps.proc = entityProc.second;
      impl::VisitClosureGeneral(*this,entityProc.first,sieps,og);
  }

  // Synchronize the send and receive list.
  // If the send list contains a not-owned entity
  // inform the owner and receiver to add that entity
  // to their ghost send and receive lists.

  stk::mesh::impl::comm_sync_send_recv( *this , sendGhosts , recvGhosts, ghostStatus );

  // The sendGhosts list is now parallel complete and parallel accurate
  // recvGhosts has those ghost entities that are to be kept.
  //------------------------------------
  // Remove the ghost entities that will not remain.
  // If the last reference to the receive ghost entity then delete it.

  OrdinalVector addParts;
  OrdinalVector removeParts(1, m_ghost_parts[ghosting.ordinal()]->mesh_meta_data_ordinal());
  OrdinalVector scratchOrdinalVec, scratchSpace;
  bool removed = false ;

  std::vector<EntityCommInfo> comm_ghost ;
  for ( EntityCommListInfoVector::reverse_iterator
        i = m_entity_comm_list.rbegin() ; i != m_entity_comm_list.rend() ; ++i) {

    if (!i->entity_comm || !i->entity_comm->isGhost) {
      continue;
    }

    const bool is_owner = parallel_owner_rank(i->entity) == parallel_rank() ;
    const bool remove_recv = ( ! is_owner ) &&
                             !ghostStatus[i->entity.local_offset()] && in_receive_ghost(ghosting, i->entity);

    if(is_valid(i->entity))
    {
      if ( is_owner ) {
        // Is owner, potentially removing ghost-sends
        // Have to make a copy

          const PairIterEntityComm ec = internal_entity_comm_map(i->entity, ghosting);
          comm_ghost.assign( ec.first , ec.second );
  
          for ( ; ! comm_ghost.empty() ; comm_ghost.pop_back() ) {
            const EntityCommInfo tmp = comm_ghost.back();

            std::set<EntityProc,EntityLess>::iterator iter =
                sendGhosts.find(EntityProc(i->entity, tmp.proc));

            if ( iter == sendGhosts.end() ) {
              entity_comm_map_erase(i->key, tmp);
            }
            else {
              sendGhosts.erase(iter);
            }
          }
      }
      else if ( remove_recv ) {
          entity_comm_map_erase(i->key, ghosting);
          internal_change_entity_parts(i->entity, addParts, removeParts, scratchOrdinalVec, scratchSpace);
      }

      if ( internal_entity_comm_map(i->entity).empty() ) {
        removed = true ;
        i->key = EntityKey(); // No longer communicated
        if ( remove_recv ) {
          ThrowRequireMsg( internal_destroy_entity_with_notification( i->entity, remove_recv ),
                           "P[" << this->parallel_rank() << "]: FAILED attempt to destroy entity: "
                           << entity_key(i->entity) );
        }
      }
    }
  }

  // if an entry in the comm_list has the EntityKey() value, it is invalid,
  // and removed from the comm_list

  if ( removed ) {
    delete_unneeded_entries_from_the_comm_list();
  }

  ghost_entities_and_fields(ghosting, sendGhosts, is_full_regen);
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void BulkData::fill_list_of_entities_to_send_for_aura_ghosting(EntityProcMapping& send)
{
    const EntityRank end_rank = static_cast<EntityRank>(m_mesh_meta_data.entity_rank_count());

    // Iterate over all entities with communication info, get the sharing
    // comm info for each entity, and ensure that upwardly related
    // entities to the shared entity are ghosted on the sharing proc.
    Entity const* rels = NULL;
    int num_rels = 0;
    for ( EntityCommListInfoVector::const_iterator
        i = internal_comm_list().begin() ; i != internal_comm_list().end() ; ++i )
    {
      if (i->entity_comm && i->entity_comm->isShared) {
        const EntityRank erank = static_cast<EntityRank>(i->key.rank());

        const PairIterEntityComm shared = shared_comm_info_range(i->entity_comm->comm_map);

        for ( size_t j = 0 ; j < shared.size() ; ++j ) {

          const int share_proc = shared[j].proc ;

          for (EntityRank k_rank = static_cast<EntityRank>(erank + 1); k_rank < end_rank; ++k_rank)
          {
            num_rels = num_connectivity(i->entity, k_rank);
            rels     = begin(i->entity, k_rank);

            for (int r = 0; r < num_rels; ++r)
            {
              if (is_valid(rels[r])) {
                stk::mesh::impl::insert_upward_relations(*this, rels[r], erank, share_proc, send);
              }
            }
          }
        }
      }
    }
}

void BulkData::internal_regenerate_aura()
{
  require_ok_to_modify();

  EntityProcMapping entityProcMapping(get_size_of_entity_index_space());
  fill_list_of_entities_to_send_for_aura_ghosting(entityProcMapping);

  internal_change_ghosting(aura_ghosting() , entityProcMapping);
}

int BulkData::determine_new_owner( Entity entity ) const
{
  // We will decide the new owner by looking at all the processes sharing
  // this entity. The new owner will be the sharing process with lowest rank.

  // The local process is a candidate only if the entity is not destroyed.
  int new_owner = is_valid(entity) ? parallel_rank() : ~0u;

  for ( PairIterEntityComm
        share = internal_entity_comm_map_shared(entity_key(entity)); ! share.empty() ; ++share ) {
    if ( share->proc < parallel_size() &&
         ( share->proc < new_owner || parallel_size() <= new_owner ) ) {
      new_owner = share->proc ;
    }
  }

  return new_owner ;
}

bool BulkData::pack_entity_modification( const bool packShared , stk::CommSparse & comm )
{
  bool flag = false;
  bool packGhosted = packShared == false;

  const EntityCommListInfoVector & entityCommList = this->internal_comm_list();

  for ( EntityCommListInfoVector::const_iterator
        i = entityCommList.begin() ; i != entityCommList.end() ; ++i ) {
    if (i->entity_comm != nullptr) {
      Entity entity = i->entity;
      EntityState status = this->is_valid(entity) ? this->state(entity) : Deleted;

      if ( status == Modified || status == Deleted ) {
        int owned_closure_int = owned_closure(entity) ? 1 : 0;
  
        for ( PairIterEntityComm ec(i->entity_comm->comm_map); ! ec.empty() ; ++ec )
        {
          if ( ( packGhosted && ec->ghost_id > BulkData::SHARED ) || ( packShared && ec->ghost_id == BulkData::SHARED ) )
          {
            comm.send_buffer( ec->proc )
                .pack<EntityKey>( i->key )
                .pack<EntityState>( status )
                .pack<int>(owned_closure_int);
  
            flag = true ;
          }
        }
      }
    }
  }

  return flag ;
}

void BulkData::communicate_entity_modification( const bool shared , std::vector<EntityParallelState > & data )
{
  stk::CommSparse comm( this->parallel() );
  const int p_size = comm.parallel_size();

  // Sizing send buffers:
  pack_entity_modification(shared , comm);

  const bool needToSendOrRecv = comm.allocate_buffers();
  if ( needToSendOrRecv )
  {

    // Packing send buffers:
    pack_entity_modification(shared , comm);

    comm.communicate();

    const EntityCommListInfoVector & entityCommList = this->internal_comm_list();
    for ( int procNumber = 0 ; procNumber < p_size ; ++procNumber ) {
      CommBuffer & buf = comm.recv_buffer( procNumber );
      EntityKey key;
      EntityState state;
      int remote_owned_closure_int;
      bool remote_owned_closure;

      while ( buf.remaining() ) {

        buf.unpack<EntityKey>( key )
           .unpack<EntityState>( state )
           .unpack<int>( remote_owned_closure_int);
        remote_owned_closure = ((remote_owned_closure_int==1)?true:false);

        // search through entity_comm, should only receive info on entities
        // that are communicated.
        EntityCommListInfo info = find_entity(*this, entityCommList, key);
        EntityParallelState parallel_state = {procNumber, state, info, remote_owned_closure, this};
        data.push_back( parallel_state );
      }
    }
  }

  std::sort( data.begin() , data.end() );
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

// Postconditions:
//  * Comm lists for shared entities are up-to-date.
//  * shared_new contains all entities that were modified/created on a
//    different process

void BulkData::extract_entity_from_shared_entity_type(const std::vector<shared_entity_type>& shared_entities, std::vector<Entity>& shared_new)
{
    for (size_t i=0; i<shared_entities.size(); ++i)
    {
        Entity entity = this->get_entity(shared_entities[i].global_key);
        if ( internal_is_entity_marked(entity) == BulkData::IS_SHARED )
        {
            shared_new.push_back(entity);
        }
    }
}

stk::mesh::Permutation get_permutation(stk::mesh::BulkData &bulk, stk::mesh::Entity element, const stk::mesh::EntityVector& nodes)
{
    std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> ordinalAndPermutation =
                  stk::mesh::get_ordinal_and_permutation(bulk, element, bulk.mesh_meta_data().side_rank(), nodes);
    return ordinalAndPermutation.second;
}

unsigned get_index_of_side_in_element_bucket(stk::mesh::BulkData &bulk, stk::mesh::Entity element, stk::mesh::Entity side)
{
    unsigned elements_edge_offset = stk::mesh::INVALID_CONNECTIVITY_ORDINAL;
    unsigned num_edges_or_faces = bulk.num_connectivity(element, bulk.entity_rank(side));
    const stk::mesh::Entity* entities = bulk.begin(element, bulk.entity_rank(side));
    for(unsigned j=0;j<num_edges_or_faces;++j)
    {
        if (entities[j]==side)
        {
            elements_edge_offset = static_cast<stk::mesh::ConnectivityOrdinal>(j);
            break;
        }
    }
    return elements_edge_offset;
}

void BulkData::change_connectivity_for_edge_or_face(stk::mesh::Entity side, const std::vector<stk::mesh::EntityKey>& node_keys)
{
    stk::mesh::EntityVector nodes = impl::convert_keys_to_entities(*this, node_keys);

    unsigned num_elems = this->num_elements(side);
    const stk::mesh::Entity *elements = this->begin_elements(side);
    for (unsigned i=0;i<num_elems;++i)
    {
        if(bucket(elements[i]).owned())
        {
            stk::mesh::Bucket& bucket_edge = bucket(side);
            bucket_edge.change_existing_connectivity(bucket_ordinal(side), nodes.data());

            stk::mesh::Permutation new_permutation = get_permutation(*this, elements[i], nodes);
            ThrowRequireWithSierraHelpMsg(new_permutation!=stk::mesh::INVALID_PERMUTATION);

            unsigned edges_element_offset = static_cast<stk::mesh::ConnectivityOrdinal>(i);
            bucket_edge.change_existing_permutation_for_connected_element(bucket_ordinal(side), edges_element_offset, new_permutation);

            unsigned elements_edge_offset = get_index_of_side_in_element_bucket(*this, elements[i], side);
            stk::mesh::Bucket& bucket_elem = bucket(elements[i]);
            bucket_elem.change_existing_permutation_for_connected_edge(bucket_ordinal(elements[i]), elements_edge_offset, new_permutation);
        }
    }
}

void BulkData::update_side_elem_permutations(Entity side)
{
    const stk::mesh::Entity* sideNodes = begin_nodes(side);

    unsigned numElems = num_elements(side);
    const stk::mesh::Entity* elems = begin_elements(side);
    const stk::mesh::ConnectivityOrdinal* side_elem_ords = begin_element_ordinals(side);
    const stk::mesh::Permutation* side_elem_perms = begin_element_permutations(side);
    for(unsigned i=0; i<numElems; ++i) {
        stk::EquivalentPermutation equiv = stk::mesh::side_equivalent(*this, elems[i], side_elem_ords[i], sideNodes);
        ThrowRequireMsg(equiv.is_equivalent,"ERROR, side not equivalent to element-side");
        stk::mesh::Permutation correctPerm = static_cast<stk::mesh::Permutation>(equiv.permutation_number);
        internal_declare_relation(elems[i], side, side_elem_ords[i], correctPerm);
        ThrowRequireMsg(side_elem_perms[i] == correctPerm, "ERROR, failed to set permutation");
    }
}

void BulkData::resolve_parallel_side_connections(std::vector<SideSharingData>& sideSharingDataToSend,
                                                 std::vector<SideSharingData>& sideSharingDataReceived)
{
    std::vector<stk::mesh::impl::IdViaSidePair> idAndSides;

    stk::mesh::EntityVector sharedSides = fill_shared_entities_that_need_fixing(*this);

    fill_sharing_data(*this, *m_elemElemGraph, sharedSides, sideSharingDataToSend, idAndSides);

    stk::CommSparse comm(parallel());
    allocate_and_send(comm, sideSharingDataToSend, idAndSides);
    unpack_data(comm, parallel_rank(), parallel_size(), sideSharingDataReceived);
    for(SideSharingData &sideSharingData : sideSharingDataReceived)
    {
        stk::mesh::Entity element = get_entity(stk::topology::ELEM_RANK, sideSharingData.elementAndSide.id);
        int sideOrdinal = sideSharingData.elementAndSide.side;

        stk::mesh::Entity side = stk::mesh::get_side_entity_for_elem_side_pair(*this, element, sideOrdinal);
        if(!is_valid(side))
        {
            stk::mesh::PartVector addParts;
            stk::mesh::impl::convert_part_ordinals_to_parts(mesh_meta_data(), sideSharingData.partOrdinals, addParts);
            side = declare_element_side(element, sideOrdinal, addParts);
        }
        else
        {
            if (parallel_rank() != sideSharingData.owningProc)
            {
                stk::mesh::EntityKey newKey(mesh_meta_data().side_rank(), sideSharingData.chosenSideId);
                stk::mesh::EntityKey existingKey = entity_key(side);
                if(newKey != existingKey)
                    internal_change_entity_key(existingKey, newKey, side);
            }
        }
        sideSharingData.side = side;
    }
}

void BulkData::add_comm_map_for_sharing(const std::vector<SideSharingData>& sidesSharingData, stk::mesh::EntityVector& shared_entities)
{
    for(const SideSharingData &sideSharingData : sidesSharingData)
    {
        shared_entities.push_back(sideSharingData.side);

        for(int sharingProc : sideSharingData.allSharingProcs)
        {
            if(sharingProc != parallel_rank())
                entity_comm_map_insert(sideSharingData.side, stk::mesh::EntityCommInfo(stk::mesh::BulkData::SHARED, sharingProc));
        }
    }
}

void BulkData::use_elem_elem_graph_to_determine_shared_entities(std::vector<stk::mesh::Entity>& shared_entities)
{
    std::vector<SideSharingData> sideSharingDataReceived;
    std::vector<SideSharingData> sideSharingDataToSend;

    resolve_parallel_side_connections(sideSharingDataToSend, sideSharingDataReceived);

    add_comm_map_for_sharing(sideSharingDataToSend, shared_entities);
    add_comm_map_for_sharing(sideSharingDataReceived, shared_entities);

    std::sort(shared_entities.begin(), shared_entities.end(), stk::mesh::EntityLess(*this));
}

void BulkData::use_nodes_to_resolve_sharing(stk::mesh::EntityRank rank, std::vector<Entity>& shared_new, bool onlyConsiderSoloSides)
{
    std::vector<shared_entity_type> potentially_shared_entities;
    markEntitiesForResolvingSharingInfoUsingNodes(rank, onlyConsiderSoloSides, potentially_shared_entities);
    set_common_entity_key_and_fix_ordering_of_nodes_and_update_comm_map( potentially_shared_entities );
    extract_entity_from_shared_entity_type(potentially_shared_entities, shared_new);
}

void BulkData::fill_shared_entities_of_rank_while_updating_sharing_info(stk::mesh::EntityRank rank, std::vector<Entity> &shared_new)
{
    if (this->has_face_adjacent_element_graph() && rank == mesh_meta_data().side_rank())
    {
        use_elem_elem_graph_to_determine_shared_entities(shared_new);
        const bool onlyConsiderSoloSides = true;
        use_nodes_to_resolve_sharing(rank, shared_new, onlyConsiderSoloSides);
    }
    else
    {
        use_nodes_to_resolve_sharing(rank, shared_new);
    }
}

void BulkData::internal_update_sharing_comm_map_and_fill_list_modified_shared_entities(
        std::vector<Entity> & shared_new )
{
    std::vector<Entity> shared_nodes;
    this->gather_shared_nodes(shared_nodes);

    shared_new.clear();
    shared_new.reserve(shared_nodes.size()*2);
    shared_new.insert(shared_new.end(), shared_nodes.begin(), shared_nodes.end());

    this->fill_shared_entities_of_rank_while_updating_sharing_info(stk::topology::EDGE_RANK, shared_new);
    this->fill_shared_entities_of_rank_while_updating_sharing_info(stk::topology::FACE_RANK, shared_new);

    std::fill(m_mark_entity.begin(), m_mark_entity.end(), BulkData::NOT_MARKED);
}



//----------------------------------------------------------------------

void BulkData::internal_establish_new_owner(stk::mesh::Entity entity)
{
    const int new_owner = determine_new_owner(entity);

    internal_set_owner(entity, new_owner);
}

void BulkData::internal_update_parts_for_shared_entity(stk::mesh::Entity entity, const bool is_entity_shared, const bool did_i_just_become_owner)
{
    OrdinalVector parts_to_add_entity_to , parts_to_remove_entity_from, scratchOrdinalVec, scratchSpace;

    if ( !is_entity_shared ) {
      parts_to_remove_entity_from.push_back(m_mesh_meta_data.globally_shared_part().mesh_meta_data_ordinal());
    }

    if ( did_i_just_become_owner ) {
      parts_to_add_entity_to.push_back(m_mesh_meta_data.locally_owned_part().mesh_meta_data_ordinal());
    }

    if ( ! parts_to_add_entity_to.empty() || ! parts_to_remove_entity_from.empty() ) {
      internal_change_entity_parts( entity , parts_to_add_entity_to , parts_to_remove_entity_from, scratchOrdinalVec, scratchSpace );
    }
}

void BulkData::filter_upward_ghost_relations(const Entity entity, std::function<void(Entity)> filter)
{
  EntityRank rank = entity_rank(entity);
  EntityRank endRank = static_cast<EntityRank>(mesh_meta_data().entity_rank_count());

  for(EntityRank iRank = static_cast<EntityRank>(rank+1); iRank < endRank; iRank++) {
    unsigned numRelations = num_connectivity(entity, iRank);
    const Entity* relations = begin(entity, iRank);

    for(unsigned i = 0; i < numRelations; i++) {
      filter(relations[i]);
    }
  }
}

EntityVector BulkData::get_upward_recv_ghost_relations(const Entity entity)
{
  EntityVector ghosts;
  filter_upward_ghost_relations(entity, [&](Entity relation) {
    if(in_receive_ghost(entity_key(relation))) {
      ghosts.push_back(relation);
    }
  });

  return ghosts;
}

EntityVector BulkData::get_upward_send_ghost_relations(const Entity entity)
{
  EntityVector ghosts;
  filter_upward_ghost_relations(entity, [&](Entity relation) {
    if(in_send_ghost(entity_key(relation))) {
      ghosts.push_back(relation);
    }
  });

  return ghosts;
}

void BulkData::add_entity_to_same_ghosting(Entity entity, Entity connectedGhost) {
  for(PairIterEntityComm ec(internal_entity_comm_map(entity_key(connectedGhost))); ! ec.empty(); ++ec) {
    if (ec->ghost_id > BulkData::AURA) {
      entity_comm_map_insert(entity, EntityCommInfo(ec->ghost_id, ec->proc));
      entity_comm_list_insert(entity);
    }
  }
}

void BulkData::internal_resolve_formerly_shared_entities(const EntityVector& entitiesNoLongerShared)
{
  for(Entity entity : entitiesNoLongerShared) {
    EntityVector ghostRelations = get_upward_send_ghost_relations(entity);

    for(Entity ghost : ghostRelations) {
      add_entity_to_same_ghosting(entity, ghost);
    }
  }
}

//----------------------------------------------------------------------
// Resolve modifications for ghosted entities:
// If a ghosted entity is modified or destroyed on the owning
// process then the ghosted entity must be destroyed.
//
// Post condition:
//  Ghosted entities of modified or deleted entities are destroyed.
//  Ghosted communication lists are cleared to reflect all deletions.

void BulkData::internal_resolve_ghosted_modify_delete(const stk::mesh::EntityVector& entitiesNoLongerShared)
{
  ThrowRequireMsg(parallel_size() > 1, "Do not call this in serial");
  // Resolve modifications for ghosted entities:

  std::vector<EntityParallelState > remotely_modified_ghosted_entities ;
  internal_resolve_formerly_shared_entities(entitiesNoLongerShared);

  // Communicate entity modification state for ghost entities
  const bool communicate_shared = false ;
  communicate_entity_modification( communicate_shared , remotely_modified_ghosted_entities );

  const size_t ghosting_count = m_ghosting.size();
  const size_t ghosting_count_minus_shared = ghosting_count - 1;

  // We iterate backwards over remote_mod to ensure that we hit the
  // higher-ranking entities first. This is important because higher-ranking
  // entities like element must be deleted before the nodes they have are
  // deleted.
  for ( std::vector<EntityParallelState>::reverse_iterator
        i = remotely_modified_ghosted_entities.rbegin(); i != remotely_modified_ghosted_entities.rend() ; ++i )
  {
    Entity entity                 = i->comm_info.entity;
    const EntityKey key           = i->comm_info.key;
    const int      remote_proc    = i->from_proc;
    const bool     local_owner    = parallel_owner_rank(entity) == parallel_rank() ;
    const bool remotely_destroyed = Deleted == i->state ;
    const bool remote_proc_is_owner = remote_proc == parallel_owner_rank(entity);
    const bool isAlreadyDestroyed  = !is_valid(entity);
    if ( local_owner ) { // Sending to 'remote_proc' for ghosting

      if ( remotely_destroyed ) {

        // remove from ghost-send list

        for ( size_t j = ghosting_count_minus_shared ; j>=1 ; --j) {
          entity_comm_map_erase( key, EntityCommInfo( j , remote_proc ) );
        }
      }

      // Remotely modified ghosts are ignored
    }
    else if (remote_proc_is_owner) { // Receiving from 'remote_proc' for ghosting

      const bool hasBeenPromotedToSharedOrOwned = this->owned_closure(entity);
      bool isAuraGhost = false;
      bool isCustomGhost = false;
      PairIterEntityComm pairIterEntityComm = internal_entity_comm_map(key);
      if(pairIterEntityComm.empty()) {
        if(std::binary_search(entitiesNoLongerShared.begin(), entitiesNoLongerShared.end(), entity)) {
          EntityVector ghosts = get_upward_recv_ghost_relations(entity);

          for(Entity ghost : ghosts) {
            add_entity_to_same_ghosting(entity, ghost);
          }
        }
      } else {
        for(unsigned j=0; j<pairIterEntityComm.size(); ++j)
        {
          if (pairIterEntityComm[j].ghost_id == AURA)
          {
            isAuraGhost = true;
          }
          else if (pairIterEntityComm[j].ghost_id > AURA)
          {
            isCustomGhost = true;
          }
        }
      }

      if ( isAuraGhost ) {
        entity_comm_map_erase(key, aura_ghosting());
      }

      if(!isAlreadyDestroyed)
      {
        const bool wasDestroyedByOwner = remotely_destroyed;
        const bool shouldDestroyGhost = wasDestroyedByOwner || (isAuraGhost && !isCustomGhost && !hasBeenPromotedToSharedOrOwned);
        const bool shouldRemoveFromGhosting = remotely_destroyed && !isAuraGhost && hasBeenPromotedToSharedOrOwned;

        if (shouldRemoveFromGhosting) {
            for ( size_t j = ghosting_count_minus_shared ; j >=1 ; --j ) {
                entity_comm_map_erase( key, *m_ghosting[j] );
            }
        }

        if ( shouldDestroyGhost )
        {
          const bool was_ghost = true;
          internal_destroy_entity_with_notification(entity, was_ghost);
        }

        entity_comm_list_insert(entity);
      }
    }
  } // end loop on remote mod

  // Erase all ghosting communication lists for:
  // 1) Destroyed entities.
  // 2) Owned and modified entities.

  for ( EntityCommListInfoVector::const_reverse_iterator
        i = internal_comm_list().rbegin() ; i != internal_comm_list().rend() ; ++i) {

    Entity entity = i->entity;

    const bool locally_destroyed = !is_valid(entity);
    const bool locally_owned_and_modified = locally_destroyed ? false :
      (Modified == state(entity) && (parallel_rank() == parallel_owner_rank(entity)));

    if ( locally_destroyed ) {
      for ( size_t j = ghosting_count_minus_shared ; j >=1 ; --j ) {
        entity_comm_map_erase( i->key, *m_ghosting[j] );
      }
    }
    else if ( locally_owned_and_modified ) {
      entity_comm_map_erase( i->key, aura_ghosting() );
    }
  }
}

void BulkData::resolve_ownership_of_modified_entities( const std::vector<Entity> &shared_modified )
{
    stk::CommSparse comm_sparse( parallel() );

    for ( int phase = 0; phase < 2; ++phase ) {
        for (Entity entity : shared_modified) {
            if ( parallel_owner_rank(entity) == parallel_rank() && state(entity) != Created ) {
                for ( PairIterEntityComm jc = internal_entity_comm_map_shared(entity_key(entity)) ; ! jc.empty() ; ++jc ) {
                    comm_sparse.send_buffer( jc->proc ) .pack<EntityKey>( entity_key(entity) );
                }
            }
        }

        if (phase == 0) { //allocation phase
            comm_sparse.allocate_buffers();
        }
        else { // communication phase
            comm_sparse.communicate();
        }
    }

    for ( int receive_proc = 0 ; receive_proc < parallel_size() ; ++receive_proc ) {
        CommBuffer & buf = comm_sparse.recv_buffer( receive_proc );
        EntityKey key ;
        while ( buf.remaining() ) {
            buf.unpack<EntityKey>( key );
            Entity entity = get_entity( key );

            // Set owner, will correct part membership later
            internal_set_owner(entity, receive_proc);
        }
    }
}

void BulkData::move_entities_to_proper_part_ownership( const std::vector<Entity> &shared_modified )
{
    std::ostringstream error_msg;
    int error_flag = 0;

    OrdinalVector shared_part, owned_part, empty;
    OrdinalVector scratchOrdinalVec, scratchSpace;
    shared_part.push_back(m_mesh_meta_data.globally_shared_part().mesh_meta_data_ordinal());
    owned_part.push_back(m_mesh_meta_data.locally_owned_part().mesh_meta_data_ordinal());

    std::vector<Entity>::const_reverse_iterator iend = shared_modified.rend();
    for(std::vector<Entity>::const_reverse_iterator i = shared_modified.rbegin(); i != iend; ++i)
    {
        Entity entity = *i;

        if(parallel_owner_rank(entity) == parallel_rank() && state(entity) == Created)
        {
            // Created and not claimed by an existing owner

            const int new_owner = determine_new_owner(entity);

            internal_set_owner(entity, new_owner);
        }

        if(parallel_owner_rank(entity) != parallel_rank())
        {
            // Do not own it and still have it.
            // Remove the locally owned, add the globally_shared
            internal_change_entity_parts(entity, shared_part /*add*/, owned_part /*remove*/, scratchOrdinalVec, scratchSpace);
        }
        else if(!internal_entity_comm_map_shared(entity_key(entity)).empty())
        {
            // Own it and has sharing information.
            // Add the globally_shared
            unprotect_orphaned_node(entity);
            internal_change_entity_parts(entity, shared_part /*add*/, empty /*remove*/, scratchOrdinalVec, scratchSpace);
        }
        else
        {
            // Own it and does not have sharing information.
            // Remove the globally_shared
            internal_change_entity_parts(entity, empty /*add*/, shared_part /*remove*/, scratchOrdinalVec, scratchSpace);
        }

        // Newly created shared entity had better be in the owned closure
        bool isEntityGhost = !this->owned_closure(entity) && state(entity)==Created;
        if(isEntityGhost)
        {
            if(0 == error_flag)
            {
                error_flag = 1;
                error_msg << "\nP" << parallel_rank() << ": " << " FAILED\n"
                          << "  The following entities were declared on multiple processors,\n"
                          << "  cannot be parallel-shared, and were declared with"
                          << "  parallel-ghosting information. {\n";
            }
            error_msg << "    " << impl::print_entity_key(m_mesh_meta_data, entity_key(entity));
            error_msg << " also declared on";
            for(PairIterEntityComm ec = internal_entity_comm_map_shared(entity_key(entity)); !ec.empty(); ++ec)
            {
                error_msg << " P" << ec->proc;
            }
            error_msg << "\n";
        }
    }

    // Parallel-consistent error checking of above loop
    if(error_flag)
    {
        error_msg << "}\n";
    }
    all_reduce(parallel(), ReduceMax<1>(&error_flag));
    ThrowErrorMsgIf(error_flag, error_msg.str());
}


void BulkData::add_comm_list_entries_for_entities(const std::vector<stk::mesh::Entity>& sharedModifiedEntities)
{
    // ------------------------------------------------------------
    // Update m_entity_comm_list based on shared_modified

    const size_t n_old = m_entity_comm_list.size();

    m_entity_comm_list.reserve(m_entity_comm_list.size() + sharedModifiedEntities.size());
    for (size_t i = 0, e = sharedModifiedEntities.size(); i < e; ++i)
    {
      Entity entity = sharedModifiedEntities[i];
      EntityKey key = entity_key(entity);
      const EntityComm* entity_comm = m_entity_comm_map.entity_comm(key);
      EntityCommListInfo comm_info = {key, entity, entity_comm};
      m_entity_comm_list.push_back(comm_info);
    }

    std::sort( m_entity_comm_list.begin() + n_old , m_entity_comm_list.end() );

    std::inplace_merge( m_entity_comm_list.begin() ,
                        m_entity_comm_list.begin() + n_old ,
                        m_entity_comm_list.end() );

    EntityCommListInfoVector::iterator iter =
      std::unique( m_entity_comm_list.begin() , m_entity_comm_list.end() );

    m_entity_comm_list.erase( iter , m_entity_comm_list.end() );
}
//----------------------------------------------------------------------

// Postconditions:
//  * All shared entities have parallel-consistent owner
//  * Part membership of shared entities is up-to-date
//  * m_entity_comm is up-to-date
void BulkData::internal_resolve_parallel_create_nodes()
{
    std::vector<stk::mesh::EntityRank> ranks={stk::topology::NODE_RANK};
    this->internal_resolve_parallel_create(ranks);
}

void BulkData::internal_resolve_parallel_create_edges_and_faces()
{
    std::vector<stk::mesh::EntityRank> ranks;
    for(EntityRank rank=stk::topology::EDGE_RANK; rank<=mesh_meta_data().side_rank(); rank++)
        ranks.push_back(rank);
    this->internal_resolve_parallel_create(ranks);
}

void BulkData::internal_resolve_parallel_create()
{
    this->internal_resolve_parallel_create_nodes();
    this->internal_resolve_parallel_create_edges_and_faces();
}

void BulkData::internal_resolve_parallel_create(const std::vector<stk::mesh::EntityRank>& ranks)
{
  ThrowRequireMsg(parallel_size() > 1, "Do not call this in serial");

  for(EntityRank rank : ranks)
  {
      std::vector<Entity> shared_modified;

      internal_update_sharing_comm_map_and_fill_list_modified_shared_entities_of_rank(rank, shared_modified);

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

      add_comm_list_entries_for_entities( shared_modified );
  }
  update_comm_list_based_on_changes_in_comm_map();
  std::fill(m_mark_entity.begin(), m_mark_entity.end(), BulkData::NOT_MARKED);
}

bool BulkData::modification_end_for_entity_creation( const std::vector<EntityRank> & entity_rank_vector, impl::MeshModification::modification_optimization opt)
{
  notifier.notify_started_modification_end();

  bool return_value = internal_modification_end_for_entity_creation( entity_rank_vector, opt );
  return return_value;
}

void BulkData::update_comm_list_based_on_changes_in_comm_map()
// Resolution of shared and ghost modifications can empty
// the communication information for entities.
// If there is no communication information then the
// entity must be removed from the communication list.
{
  bool changed = false ;
  for (EntityCommListInfo& entityCommInfo : m_entity_comm_list) {
      if (entityCommInfo.entity_comm == NULL || entityCommInfo.entity_comm->comm_map.empty() ) {
          entityCommInfo.key = EntityKey();
          changed = true;
      }
  }

  if ( changed ) {
      delete_unneeded_entries_from_the_comm_list();
  }
}

void BulkData::notify_finished_mod_end()
{
    notifier.notify_finished_modification_end(parallel());
}

void BulkData::internal_modification_end_for_change_ghosting()
{
    internal_resolve_send_ghost_membership();

    m_bucket_repository.internal_default_sort_bucket_entities();

    m_modSummary.write_summary(m_meshModification.synchronized_count());
    if(parallel_size() > 1)
    {
        check_mesh_consistency();
    }

    m_bucket_repository.internal_modification_end();
    internal_update_fast_comm_maps();

    m_meshModification.set_sync_state_synchronized();
    notify_finished_mod_end();
}

bool BulkData::internal_modification_end_for_change_parts()
{
    int global_shared_modified = 0;
    int local_shared_modified = 0;
    if (m_meshModification.did_any_shared_entity_change_parts())
    {
        local_shared_modified = 1;
    }

    stk::all_reduce_max(parallel(), &local_shared_modified, &global_shared_modified, 1);

    if (this->parallel_size() > 1 && global_shared_modified > 0) {
      stk::mesh::EntityVector entitiesNoLongerShared;
      m_meshModification.internal_resolve_shared_modify_delete(entitiesNoLongerShared);
      internal_resolve_shared_membership(entitiesNoLongerShared);
    }
    m_modSummary.write_summary(m_meshModification.synchronized_count());

    m_bucket_repository.internal_default_sort_bucket_entities();

    m_bucket_repository.internal_modification_end();
    internal_update_fast_comm_maps();

    m_meshModification.set_sync_state_synchronized();
    notify_finished_mod_end();
    return true;
}

bool BulkData::internal_modification_end_for_change_entity_owner( impl::MeshModification::modification_optimization opt )
{
  // The two states are MODIFIABLE and SYNCHRONiZED
  if ( this->in_synchronized_state() ) { return false ; }

  ThrowAssertMsg(impl::check_for_connected_nodes(*this)==0, "BulkData::modification_end ERROR, all entities with rank higher than node are required to have connected nodes.");

  ThrowAssertMsg(add_fmwk_data() || impl::check_no_shared_elements_or_higher(*this)==0, "BulkData::modification_end ERROR, Sharing of entities with rank ELEMENT_RANK or higher is not allowed.");

  if (parallel_size() > 1)
  {
      if ( m_autoAuraOption == AUTO_AURA )
      {
          internal_regenerate_aura();
      }
      internal_resolve_send_ghost_membership();
      m_modSummary.write_summary(m_meshModification.synchronized_count());
      check_mesh_consistency();
  }
  else
  {
      std::vector<Entity> shared_modified ;
      internal_update_sharing_comm_map_and_fill_list_modified_shared_entities( shared_modified );
      m_modSummary.write_summary(m_meshModification.synchronized_count());
  }

  this->internal_finish_modification_end(opt);

  return true ;
}

void BulkData::check_mesh_consistency()
{
    // ------------------------------
    // Verify parallel consistency of mesh entities.
    // Unique ownership, communication lists, sharing part membership,
    // application part membership consistency.

#ifndef NDEBUG
    ThrowErrorMsgIf(!stk::mesh::impl::check_permutations_on_all(*this), "Permutation checks failed.");
    std::ostringstream msg ;
    bool is_consistent = true;
    is_consistent = comm_mesh_verify_parallel_consistency( msg );
    std::string error_msg = msg.str();
    ThrowErrorMsgIf( !is_consistent, error_msg );
#endif
}


//////////////////////////////////// Free functions to help with modification end (exp) for edges

void fillElementsConnectedToNodes(stk::mesh::BulkData &stkMeshBulkData, const stk::mesh::Entity* nodes, size_t numNodes,
                                  std::vector<stk::mesh::Entity> & elementsConnectedToNodes)
{
    stk::mesh::Entity const * elemStartNode = stkMeshBulkData.begin_elements(nodes[0]);
    stk::mesh::Entity const * elemEndNode = stkMeshBulkData.end_elements(nodes[0]);
    elementsConnectedToNodes.assign(elemStartNode, elemEndNode);
    std::sort(elementsConnectedToNodes.begin(), elementsConnectedToNodes.end());
    std::vector<stk::mesh::Entity>::iterator intersectionEnd = elementsConnectedToNodes.end();
    std::vector<stk::mesh::Entity> elems;
    std::vector<stk::mesh::Entity> elementsConnectedToNodesTemp;
    for(size_t i = 1; i < numNodes; ++i)
    {
        elementsConnectedToNodesTemp.clear();
        elementsConnectedToNodesTemp.assign(elementsConnectedToNodes.begin(), intersectionEnd);
        elemStartNode = stkMeshBulkData.begin_elements(nodes[i]);
        elemEndNode = stkMeshBulkData.end_elements(nodes[i]);
        elems.assign(elemStartNode, elemEndNode);
        std::sort(elems.begin(), elems.end());

        intersectionEnd = std::set_intersection(elems.begin(), elems.end(),
                elementsConnectedToNodesTemp.begin(), elementsConnectedToNodesTemp.end(),
                elementsConnectedToNodes.begin());
        if(intersectionEnd == elementsConnectedToNodes.begin())
            break; // Empty set
    }

    elementsConnectedToNodes.resize(intersectionEnd - elementsConnectedToNodes.begin());
}

void fillFacesConnectedToNodes(stk::mesh::BulkData &stkMeshBulkData, const stk::mesh::Entity* nodes, size_t numNodes,
              std::vector<stk::mesh::Entity> & facesConnectedToNodes)
{
    facesConnectedToNodes.clear();

    //TODO: fill faces for all edge nodes
    ThrowRequireMsg(numNodes==2, "fillFacesConnectedToNodes ERROR, num-edge-nodes must be 2 currently.");
    if ( stkMeshBulkData.num_faces(nodes[0]) > 0 && stkMeshBulkData.num_faces(nodes[1]) > 0 )
    {
        facesConnectedToNodes.resize(20);
        stk::mesh::Entity const * facesStartNode1 = stkMeshBulkData.begin_faces(nodes[0]);
        stk::mesh::Entity const * facesEndNode1 = stkMeshBulkData.end_faces(nodes[0]);
        stk::mesh::Entity const * facesStartNode2 = stkMeshBulkData.begin_faces(nodes[1]);
        stk::mesh::Entity const * facesEndNode2 = stkMeshBulkData.end_faces(nodes[1]);

        std::vector<stk::mesh::Entity> faces1(facesStartNode1, facesEndNode1);
        std::sort(faces1.begin(), faces1.end());
        std::vector<stk::mesh::Entity> faces2(facesStartNode2, facesEndNode2);
        std::sort(faces2.begin(), faces2.end());

        std::vector<stk::mesh::Entity>::iterator iter = std::set_intersection( faces1.begin(), faces1.end(),
              faces2.begin(), faces2.end(), facesConnectedToNodes.begin());

        facesConnectedToNodes.resize(iter-facesConnectedToNodes.begin());
    }
}

bool doesEdgeNeedGhostingCommunication(stk::mesh::BulkData &stkMeshBulkData, std::vector<stk::mesh::Entity>& connectedEntities)
{
    bool communicate_edge_for_ghosting = false;
    for (size_t j=0;j<connectedEntities.size();j++)
    {
        bool isEntityGhostedToAnotherProc = stkMeshBulkData.is_aura_ghosted_onto_another_proc(stkMeshBulkData.entity_key(connectedEntities[j]));
        if ( isEntityGhostedToAnotherProc )
        {
            communicate_edge_for_ghosting = true;
        }
        else
        {
            connectedEntities[j] = Entity();
        }
    }
    return communicate_edge_for_ghosting;
}

int get_element_ordinal_for_connected_face(stk::mesh::BulkData &stkMeshBulkData,
                                           const stk::mesh::Entity localElement,
                                           const stk::mesh::Entity sideEntity)
{
    const Entity * sides = stkMeshBulkData.begin(localElement, stkMeshBulkData.mesh_meta_data().side_rank());
    ConnectivityOrdinal const * ordinals = stkMeshBulkData.begin_ordinals(localElement, stkMeshBulkData.mesh_meta_data().side_rank());
    const unsigned numSides = stkMeshBulkData.num_sides(localElement);

    for(unsigned k=0; k<numSides; ++k)
    {
        if(sideEntity == sides[k])
            return ordinals[k];
    }

    return -1;
}

std::vector<bool> get_allowable_ghost_connections(stk::mesh::BulkData &stkMeshBulkData,
                                                  std::vector<stk::mesh::Entity> &entitiesConnectedToNodes,
                                                  const stk::mesh::Entity sideEntity)
{
    std::vector<bool> isAllowableConnection(entitiesConnectedToNodes.size(), false);

    for (size_t j=0; j<entitiesConnectedToNodes.size();j++)
    {
        bool isEntityGhostedOntoThisProc = stkMeshBulkData.in_receive_ghost(stkMeshBulkData.aura_ghosting(), stkMeshBulkData.entity_key(entitiesConnectedToNodes[j]));
        if ( isEntityGhostedOntoThisProc)
            isAllowableConnection[j] = true;
    }

    return isAllowableConnection;
}

void connectGhostedEntitiesToEntity(stk::mesh::BulkData &stkMeshBulkData,
                                    std::vector<stk::mesh::Entity> &entitiesConnectedToNodes,
                                    stk::mesh::Entity entity,
                                    const stk::mesh::Entity* nodes)
{
    std::vector<bool> isAllowableGhostConnection = get_allowable_ghost_connections(stkMeshBulkData, entitiesConnectedToNodes, entity);

    for (size_t j=0; j<entitiesConnectedToNodes.size();j++)
    {
        if ( isAllowableGhostConnection[j])
        {
            stk::mesh::impl::connectUpwardEntityToEntity(stkMeshBulkData, entitiesConnectedToNodes[j], entity, nodes);
        }
    }
}

void BulkData::determineEntitiesThatNeedGhosting(stk::mesh::BulkData &stkMeshBulkData,
                                                 stk::mesh::Entity edge,
                                                 std::vector<stk::mesh::Entity>& entitiesConnectedToNodes,
                                                 const stk::mesh::Entity* nodes,
                                                 std::set<EntityProc, EntityLess> &addGhostedEntities)
{
    // Grab all the entities attached to the 2 nodes
    // If the entity is ghosted and the edge is owned, then the edge needs to be ghosted.
    bool doesEdgeNeedToBeGhosted = doesEdgeNeedGhostingCommunication(stkMeshBulkData, entitiesConnectedToNodes);
    if ( doesEdgeNeedToBeGhosted )
    {
        for (size_t j=0;j<entitiesConnectedToNodes.size();j++)
        {
            if ( entitiesConnectedToNodes[j] != Entity() )
            {
                PairIterEntityComm ghosted = stkMeshBulkData.internal_entity_comm_map( stkMeshBulkData.entity_key(entitiesConnectedToNodes[j]) , stkMeshBulkData.aura_ghosting());
                for (PairIterEntityComm ec = ghosted; !ec.empty(); ++ec)
                {
                    if ( ec->proc != stkMeshBulkData.parallel_rank() )
                    {
                        bool isEdgeSharedWithOtherProc = stkMeshBulkData.in_shared(stkMeshBulkData.entity_key(edge), ec->proc);
                        if ( !isEdgeSharedWithOtherProc )
                        {
                            addGhostedEntities.insert(EntityProc(edge, ec->proc));
                        }
                    }
                }
            }
        }
    }
}

void BulkData::find_upward_connected_entities_to_ghost_onto_other_processors(stk::mesh::BulkData &mesh,
                                                                             std::set<EntityProc, EntityLess> &entitiesToGhostOntoOtherProcessors,
                                                                             EntityRank entity_rank,
                                                                             stk::mesh::Selector selected,
                                                                             bool connectFacesToPreexistingGhosts)
{
    const stk::mesh::BucketVector& entity_buckets = mesh.buckets(entity_rank);
    bool isedge = (entity_rank == stk::topology::EDGE_RANK && mesh_meta_data().spatial_dimension() == 3);

    std::vector<stk::mesh::Entity> facesConnectedToNodes;
    std::vector<stk::mesh::Entity> elementsConnectedToNodes;
    for(size_t bucketIndex = 0; bucketIndex < entity_buckets.size(); bucketIndex++)
    {
        const stk::mesh::Bucket& bucket = *entity_buckets[bucketIndex];
        size_t numNodes = bucket.topology().num_nodes();
        for(size_t entityIndex = 0; entityIndex < bucket.size(); entityIndex++)
        {
            Entity entity = bucket[entityIndex];
            if ( mesh.state(entity) != Unchanged )
            {
                const stk::mesh::Entity* nodes = mesh.begin_nodes(entity);

                if(isedge)
                {
                    fillFacesConnectedToNodes(mesh, nodes, numNodes, facesConnectedToNodes);
                    removeEntitiesNotSelected(mesh, selected, facesConnectedToNodes);
                    if(connectFacesToPreexistingGhosts)
                        connectGhostedEntitiesToEntity(mesh, facesConnectedToNodes, entity, nodes);
                }

                fillElementsConnectedToNodes(mesh, nodes, numNodes, elementsConnectedToNodes);
                removeEntitiesNotSelected(mesh, selected, elementsConnectedToNodes);
                if(connectFacesToPreexistingGhosts)
                    connectGhostedEntitiesToEntity(mesh, elementsConnectedToNodes, entity, nodes);

                if ( bucket.owned() || bucket.shared() )
                {
                    if (isedge)
                    {
                      determineEntitiesThatNeedGhosting(mesh, entity, facesConnectedToNodes, nodes, entitiesToGhostOntoOtherProcessors);
                    }
                    determineEntitiesThatNeedGhosting(mesh, entity, elementsConnectedToNodes, nodes, entitiesToGhostOntoOtherProcessors);
                }
            }
        }
    }

    std::set< EntityKey > entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs;
    stk::mesh::impl::comm_sync_send_recv(mesh, entitiesToGhostOntoOtherProcessors, entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs);
}

void BulkData::internal_finish_modification_end(impl::MeshModification::modification_optimization opt)
{
    if(opt == impl::MeshModification::MOD_END_SORT)
    {
        m_bucket_repository.internal_default_sort_bucket_entities();
    }

    m_bucket_repository.internal_modification_end();

    internal_update_fast_comm_maps();

    m_meshModification.set_sync_state_synchronized();
    m_add_node_sharing_called = false;

    update_deleted_entities_container();

#ifdef STK_DEBUG_FIELD_SYNC
    for (FieldBase * stkField : mesh_meta_data().get_fields()) {
      if (stkField->has_device_field_debug_data()) {
        stkField->reset_debug_state();
        stkField->get_debug_ngp_field()->update_debug_storage(synchronized_count());
      }
    }
#endif

    notify_finished_mod_end();
}

bool BulkData::internal_modification_end_for_skin_mesh( EntityRank entity_rank, impl::MeshModification::modification_optimization opt, stk::mesh::Selector selectedToSkin,
        const Selector * only_consider_second_element_from_this_selector)
{
  // The two states are MODIFIABLE and SYNCHRONiZED
  if ( this->in_synchronized_state() ) { return false ; }

  ThrowAssertMsg(impl::check_for_connected_nodes(*this)==0, "BulkData::modification_end ERROR, all entities with rank higher than node are required to have connected nodes.");

  if (parallel_size() > 1)
  {
      if ( !this->is_automatic_aura_on())
      {
          find_and_delete_internal_faces(entity_rank, only_consider_second_element_from_this_selector);
      }

      stk::mesh::EntityVector entitiesNoLongerShared;
      m_meshModification.internal_resolve_shared_modify_delete(entitiesNoLongerShared);
      this->internal_resolve_shared_membership(entitiesNoLongerShared);

      if(is_automatic_aura_on())
      {
          bool connectFacesToPreexistingGhosts = true;
          resolve_incremental_ghosting_for_entity_creation_or_skin_mesh(entity_rank, selectedToSkin, connectFacesToPreexistingGhosts);
      }
  }
  m_modSummary.write_summary(m_meshModification.synchronized_count());
  if (parallel_size() > 1)
  {
      check_mesh_consistency();
  }

  this->internal_finish_modification_end(opt);

  return true ;
}

void BulkData::resolve_incremental_ghosting_for_entity_creation_or_skin_mesh(EntityRank entity_rank, stk::mesh::Selector selectedToSkin, bool connectFacesToPreexistingGhosts)
{
    std::set<EntityProc, EntityLess> entitiesToGhostOntoOtherProcessors(EntityLess(*this));
    find_upward_connected_entities_to_ghost_onto_other_processors(*this, entitiesToGhostOntoOtherProcessors, entity_rank, selectedToSkin, connectFacesToPreexistingGhosts);

    ghost_entities_and_fields(aura_ghosting(), entitiesToGhostOntoOtherProcessors);
}

bool BulkData::internal_modification_end_for_entity_creation( const std::vector<EntityRank> & entity_rank_vector, impl::MeshModification::modification_optimization opt )
{
  // The two states are MODIFIABLE and SYNCHRONiZED
  if (in_synchronized_state()) { return false; }

  ThrowAssertMsg(impl::check_for_connected_nodes(*this)==0, "BulkData::modification_end ERROR, all entities with rank higher than node are required to have connected nodes.");

  if (parallel_size() > 1) {
      internal_resolve_parallel_create(entity_rank_vector);
      stk::mesh::EntityVector entitiesNoLongerShared;
      internal_resolve_shared_membership(entitiesNoLongerShared);

      if (is_automatic_aura_on()) {
          bool connectFacesToPreexistingGhosts = true;
          for (size_t rank_idx=0; rank_idx < entity_rank_vector.size(); ++rank_idx) {
              EntityRank entity_rank = entity_rank_vector[rank_idx];
              resolve_incremental_ghosting_for_entity_creation_or_skin_mesh(entity_rank, mesh_meta_data().universal_part(), connectFacesToPreexistingGhosts);
          }
      }

      m_modSummary.write_summary(m_meshModification.synchronized_count());
      check_mesh_consistency();
  }
  else {
      std::vector<Entity> shared_modified;
      for (size_t rank_idx=0; rank_idx < entity_rank_vector.size(); ++rank_idx) {
          EntityRank entity_rank = entity_rank_vector[rank_idx];
          internal_update_sharing_comm_map_and_fill_list_modified_shared_entities_of_rank( entity_rank, shared_modified );
      }
      m_modSummary.write_summary(m_meshModification.synchronized_count());
  }

  internal_finish_modification_end(opt);

  return true;
}

bool is_modified_or_created(const BulkData& bulkData, Entity entity)
{
    return bulkData.state(entity)==stk::mesh::Modified || bulkData.state(entity)==stk    ::mesh::Created;
}

void BulkData::fill_entity_procs_for_owned_modified_or_created(std::vector<EntityProc> & send_list ) const
{
    const int p_rank = this->parallel_rank();
    for(size_t i=0; i < m_entity_comm_list.size(); ++i)
    {
        const int owner = parallel_owner_rank(m_entity_comm_list[i].entity);
        stk::mesh::Entity entity = m_entity_comm_list[i].entity;

        if(owner == p_rank && is_modified_or_created(*this, entity))
        {
            const EntityComm* entity_comm = m_entity_comm_list[i].entity_comm;
            for(PairIterEntityComm ec(entity_comm->comm_map); !ec.empty(); ++ec)
            {
                EntityProc tmp(entity, ec->proc);
                send_list.push_back(tmp);
            }
        }
    }

    stk::util::sort_and_unique(send_list, EntityLess(*this));
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

struct PartStorage
{
   OrdinalVector induced_part_ordinals;
   OrdinalVector removeParts;
   OrdinalVector scratch;
   OrdinalVector scratch2;
};


void BulkData::remove_unneeded_induced_parts(stk::mesh::Entity entity, const EntityCommInfoVector& entity_comm_info,
        PartStorage& part_storage, stk::CommSparse& comm)
{
    part_storage.induced_part_ordinals.clear();
    induced_part_membership(*this, entity, part_storage.induced_part_ordinals);

    impl::unpack_induced_parts_from_sharers(part_storage.induced_part_ordinals, entity_comm_info, comm, entity_key(entity));
    impl::filter_out_unneeded_induced_parts(*this, entity, part_storage.induced_part_ordinals, part_storage.removeParts);

    internal_change_entity_parts(entity, part_storage.induced_part_ordinals, part_storage.removeParts, part_storage.scratch, part_storage.scratch2);
}

void BulkData::internal_resolve_shared_membership(const stk::mesh::EntityVector & entitiesNoLongerShared)
{
    ThrowRequireMsg(parallel_size() > 1, "Do not call this in serial");

    ParallelMachine p_comm = parallel();
    const int p_rank = parallel_rank();

    stk::CommSparse comm(p_comm);
    impl::pack_and_send_induced_parts_from_sharers_to_owners(*this, comm, m_entity_comm_list);

    PartStorage part_storage;

#ifndef NDEBUG
    bool allOk = true;
    try
    {
#endif
        for(const EntityCommListInfo& info : m_entity_comm_list)
        {
            bool i_own_this_entity_in_comm_list = is_valid(info.entity) &&
                                        (parallel_owner_rank(info.entity) == p_rank);
            if( i_own_this_entity_in_comm_list )
            {
                remove_unneeded_induced_parts(info.entity, info.entity_comm->comm_map, part_storage,  comm);
            }
        }
#ifndef NDEBUG
    }
    catch(...) { allOk = false; }

    int numericAllOk = allOk;
    int minAllOk = 0;
    stk::all_reduce_min(this->parallel(), &numericAllOk, &minAllOk, 1);
    allOk = minAllOk == 1;
    ThrowRequireWithSierraHelpMsg(allOk);
#endif

    OrdinalVector scratch, scratch2;
    for (stk::mesh::Entity entity : entitiesNoLongerShared) {
      part_storage.induced_part_ordinals.clear();
      part_storage.removeParts.clear();
      induced_part_membership(*this, entity, part_storage.induced_part_ordinals);
      impl::filter_out_unneeded_induced_parts(*this, entity, part_storage.induced_part_ordinals, part_storage.removeParts);
      internal_change_entity_parts(entity, part_storage.induced_part_ordinals, part_storage.removeParts, scratch, scratch2);
    }

    std::vector<EntityProc> send_list;
    fill_entity_procs_for_owned_modified_or_created(send_list);
    internal_send_part_memberships_from_owner(send_list);
}

bool is_less_than_element_rank(const BulkData& bulkData, Entity entity)
{
    return bulkData.entity_rank(entity) <= bulkData.mesh_meta_data().side_rank();
}

void BulkData::internal_resolve_shared_part_membership_for_element_death()
{
    ThrowRequireMsg(parallel_size() > 1, "Do not call this in serial");

    ParallelMachine p_comm = parallel();
    const int p_rank = parallel_rank();

    stk::CommSparse comm(p_comm);

    const EntityCommListInfoVector& entityCommList = m_entity_comm_list;
    pack_and_communicate(comm, [this, &comm, &entityCommList]() {
        impl::pack_induced_memberships_for_entities_less_than_element_rank(*this, comm, entityCommList);
    });

    OrdinalVector induced_part_ordinals;
    PartStorage part_storage;

    for(EntityCommListInfoVector::iterator i = m_entity_comm_list.begin(); i != m_entity_comm_list.end(); ++i)
    {
        stk::mesh::Entity entity = i->entity;

        if(is_less_than_element_rank(*this, entity) && is_modified_or_created(*this, entity))
        {
            bool i_own_this_entity_in_comm_list = parallel_owner_rank(i->entity) == p_rank;
            if( i_own_this_entity_in_comm_list )
            {
                remove_unneeded_induced_parts(entity, i->entity_comm->comm_map, part_storage,  comm);
            }
        }
    }

    std::vector<EntityProc> send_list;
    fill_entity_procs_for_owned_modified_or_created(send_list);
    internal_send_part_memberships_from_owner(send_list);
}

void BulkData::internal_send_part_memberships_from_owner(const std::vector<EntityProc> &send_list)
{
    ParallelMachine p_comm = parallel();
    const int p_size = parallel_size();
    const PartVector & all_parts = m_mesh_meta_data.get_parts();

    stk::CommSparse comm(p_comm);

    pack_and_communicate(comm, [this, &comm, &send_list]() {
      impl::pack_part_memberships(*this, comm, send_list);
    });

    OrdinalVector owner_parts, remove_parts, parts_removed;
    OrdinalVector scratchOrdinalVec, scratchSpace;

    const MetaData & meta = m_mesh_meta_data;
    for(int p = 0; p < p_size; ++p)
    {
        CommBuffer & buf = comm.recv_buffer(p);
        while(buf.remaining())
        {
            EntityKey key;
            buf.unpack<EntityKey>(key);
            unsigned count = 0;
            buf.unpack<unsigned>(count);
            owner_parts.clear();
            for(unsigned j = 0; j < count; ++j)
            {
                unsigned part_ord = 0;
                buf.unpack<unsigned>(part_ord);
                if (all_parts[part_ord]->entity_membership_is_parallel_consistent()) {
                    stk::util::insert_keep_sorted_and_unique(part_ord, owner_parts);
                }
            }

            // Any current part that is not a member of owners_parts
            // must be removed.

            remove_parts.clear();

            Entity const entity = find_entity(*this, m_entity_comm_list, key).entity;

            const PartVector& current_parts = this->bucket(entity).supersets();

            for(PartVector::const_iterator
            ip = current_parts.begin(); ip != current_parts.end(); ++ip)
            {
                Part * const part = *ip;
                const unsigned part_ord = part->mesh_meta_data_ordinal();

                if(meta.universal_part().mesh_meta_data_ordinal() != part_ord &&
                   meta.locally_owned_part().mesh_meta_data_ordinal() != part_ord &&
                   meta.globally_shared_part().mesh_meta_data_ordinal() != part_ord &&
                   !contain(m_ghost_parts, *part) &&
                   !contains_ordinal(owner_parts.begin(), owner_parts.end(), part_ord))
                {
                    remove_parts.push_back(part_ord);
                }
            }

            internal_change_entity_parts_without_propogating_to_downward_connected_entities(entity, owner_parts, remove_parts, parts_removed, scratchOrdinalVec, scratchSpace);
        }
    }
}

void BulkData::internal_resolve_send_ghost_membership()
{
    // This virtual method can be removed when we no longer need the
    // StkTransitionBulkData derived class in Framework.
}

void add_bucket_and_ord(unsigned bucket_id, unsigned bucket_ord, BucketIndices& bktIndices)
{
  if (bktIndices.bucket_info.empty() || bktIndices.bucket_info.back().bucket_id != bucket_id) {
    bktIndices.bucket_info.push_back(BucketInfo(bucket_id, 0u));
  }

  bktIndices.bucket_info.back().num_entities_this_bucket += 1;
  bktIndices.ords.push_back(bucket_ord);
}

void BulkData::internal_update_fast_comm_maps()
{
    const EntityRank num_ranks = static_cast<EntityRank>(m_mesh_meta_data.entity_rank_count());
    m_all_sharing_procs.resize(num_ranks);
    m_volatile_fast_shared_comm_map.resize(num_ranks);
    if (parallel_size() > 1) {
        EntityCommListInfoVector& all_comm = m_entity_comm_list;

        // Flush previous map
        for (EntityRank r = stk::topology::BEGIN_RANK; r < num_ranks; ++r) {
            m_all_sharing_procs[r].clear();
            m_volatile_fast_shared_comm_map[r].resize(parallel_size());
            for (int proc = 0; proc < parallel_size(); ++proc) {
                m_volatile_fast_shared_comm_map[r][proc].bucket_info.clear();
                m_volatile_fast_shared_comm_map[r][proc].ords.clear();
            }
        }

        // Assemble map, find all shared entities and pack into volatile fast map
        std::vector<std::vector<unsigned> > shared_entity_counts(num_ranks);
        for (EntityRank r = stk::topology::BEGIN_RANK; r < num_ranks; ++r) {
            shared_entity_counts[r].assign(parallel_size(), 0);
        }

        for (size_t i = 0, ie = all_comm.size(); i < ie; ++i) {
            EntityKey const key   = all_comm[i].key;
            EntityRank const rank = key.rank();

            if (all_comm[i].entity_comm != nullptr) {
                PairIterEntityComm ec(all_comm[i].entity_comm->comm_map);
                for(; !ec.empty() && ec->ghost_id == BulkData::SHARED; ++ec) {
                    shared_entity_counts[rank][ec->proc]++;
                }
            }
        }

        for (EntityRank r = stk::topology::BEGIN_RANK; r < num_ranks; ++r) {
            for(int p=0; p<parallel_size(); ++p) {
                if (shared_entity_counts[r][p] > 0) {
                    m_volatile_fast_shared_comm_map[r][p].ords.reserve(shared_entity_counts[r][p]);
                }
            }
        }

        for (size_t i = 0, ie = all_comm.size(); i < ie; ++i) {
            Entity const e        = all_comm[i].entity;
            EntityKey const key   = all_comm[i].key;
            MeshIndex const& idx  = mesh_index(e);

            EntityRank const rank = key.rank();

            unsigned bucket_id  = idx.bucket->bucket_id();
            unsigned bucket_ord = idx.bucket_ordinal;

            if (all_comm[i].entity_comm != nullptr) {
                PairIterEntityComm ec(all_comm[i].entity_comm->comm_map);
                for(; !ec.empty() && ec->ghost_id == BulkData::SHARED; ++ec) {
                    add_bucket_and_ord(bucket_id, bucket_ord, m_volatile_fast_shared_comm_map[rank][ec->proc]);
                    stk::util::insert_keep_sorted_and_unique(ec->proc, m_all_sharing_procs[rank]);
                }
            }
        }
    }
}
 
template<typename PARTVECTOR>
void internal_throw_error_if_manipulating_internal_part_memberships(const PARTVECTOR & parts)
{
    for(size_t i=0; i<parts.size(); i++)
    {
        ThrowErrorMsgIf( stk::mesh::is_auto_declared_part(*parts[i]) && !stk::mesh::is_topology_root_part(*parts[i]),
                         "Cannot add or remove entity from internally managed part " << parts[i]->name() );
    }
}

template<typename PARTVECTOR>
void BulkData::change_entity_parts( Entity entity,
    const PARTVECTOR & add_parts ,
    const PARTVECTOR & remove_parts)
{
    bool stkMeshRunningUnderFramework = m_add_fmwk_data;
    if(!stkMeshRunningUnderFramework)
    {
        internal_throw_error_if_manipulating_internal_part_memberships(add_parts);
        internal_throw_error_if_manipulating_internal_part_memberships(remove_parts);
    }
    internal_verify_and_change_entity_parts(entity, add_parts, remove_parts);
}

template void BulkData::change_entity_parts(Entity, const PartVector&, const PartVector&);
template void BulkData::change_entity_parts(Entity, const ConstPartVector&, const ConstPartVector&);

template<typename PARTVECTOR>
void BulkData::change_entity_parts( const EntityVector& entities,
    const PARTVECTOR & add_parts ,
    const PARTVECTOR & remove_parts)
{
    bool stkMeshRunningUnderFramework = m_add_fmwk_data;
    if(!stkMeshRunningUnderFramework)
    {
        internal_throw_error_if_manipulating_internal_part_memberships(add_parts);
        internal_throw_error_if_manipulating_internal_part_memberships(remove_parts);
    }
    internal_verify_and_change_entity_parts(entities, add_parts, remove_parts);
} 
  
template void BulkData::change_entity_parts(const EntityVector&, const PartVector&, const PartVector&);
template void BulkData::change_entity_parts(const EntityVector&, const ConstPartVector&, const ConstPartVector&);

void BulkData::batch_change_entity_parts( const stk::mesh::EntityVector& entities,
                          const std::vector<PartVector>& add_parts,
                          const std::vector<PartVector>& remove_parts)
{
    bool stkMeshRunningUnderFramework = m_add_fmwk_data;
    if(!stkMeshRunningUnderFramework)
    {
        for(size_t i=0; i<add_parts.size(); i++)
        {
            internal_throw_error_if_manipulating_internal_part_memberships(add_parts[i]);
        }
        for(size_t i=0; i<remove_parts.size(); i++)
        {
            internal_throw_error_if_manipulating_internal_part_memberships(remove_parts[i]);
        }
    }

    bool starting_modification = modification_begin();
    ThrowRequireMsg(starting_modification, "ERROR: BulkData already being modified,\n"
                    <<"BulkData::change_entity_parts(vector-of-entities) can not be called within an outer modification scope.");

    for(size_t i=0; i<entities.size(); ++i) {
      internal_verify_and_change_entity_parts(entities[i], add_parts[i], remove_parts[i]);
    }

    internal_modification_end_for_change_parts();
}

void BulkData::batch_change_entity_parts( const stk::mesh::EntityVector& entities,
                          const PartVector& add_parts,
                          const PartVector& remove_parts)
{
    bool stkMeshRunningUnderFramework = m_add_fmwk_data;
    if(!stkMeshRunningUnderFramework)
    {
        internal_throw_error_if_manipulating_internal_part_memberships(add_parts);
        internal_throw_error_if_manipulating_internal_part_memberships(remove_parts);
    }

    bool starting_modification = modification_begin();
    ThrowRequireMsg(starting_modification, "ERROR: BulkData already being modified,\n"
                    <<"BulkData::change_entity_parts(vector-of-entities) can not be called within an outer modification scope.");

    for(size_t i=0; i<entities.size(); ++i) {
      internal_verify_and_change_entity_parts(entities[i], add_parts, remove_parts);
    }

    internal_modification_end_for_change_parts();
}

void BulkData::internal_adjust_entity_and_downward_connectivity_closure_count(stk::mesh::Entity entity, stk::mesh::Bucket *bucket_old, int closureCountAdjustment)
{
    m_closure_count[entity.local_offset()] += closureCountAdjustment;

    // update downward connectivity closure count
    if(bucket_old)
    {
        for(EntityRank rank = stk::topology::NODE_RANK, end_rank = bucket_old->entity_rank(); rank < end_rank; ++rank)
        {
            unsigned num = num_connectivity(entity, rank);
            Entity const * entities = begin(entity, rank);
            for(unsigned i = 0; i < num; ++i)
            {
                m_closure_count[entities[i].local_offset()] += closureCountAdjustment;
            }
        }
    }
}

void BulkData::internal_change_entity_parts_without_propogating_to_downward_connected_entities(Entity entity, const OrdinalVector& add_parts, const OrdinalVector& remove_parts, OrdinalVector& parts_removed, OrdinalVector& newBucketPartList, OrdinalVector& scratchSpace)
{
    internal_adjust_closure_count(entity, add_parts, remove_parts);

    internal_fill_new_part_list_and_removed_part_list(entity, add_parts, remove_parts, newBucketPartList, parts_removed);
    internal_move_entity_to_new_bucket(entity, newBucketPartList, scratchSpace);
}

void BulkData::internal_determine_inducible_parts(Entity entity, const OrdinalVector& add_parts, const OrdinalVector& parts_removed, OrdinalVector& inducible_parts_added, OrdinalVector& inducible_parts_removed)
{
    EntityRank e_rank = entity_rank(entity);

    impl::fill_inducible_parts_from_list(mesh_meta_data(), add_parts, e_rank, inducible_parts_added);
    impl::fill_inducible_parts_from_list(mesh_meta_data(), parts_removed, e_rank, inducible_parts_removed);
}

void BulkData::internal_determine_inducible_parts_and_propagate_to_downward_connected_entities(Entity entity, const OrdinalVector& add_parts, const OrdinalVector& parts_removed)
{
    OrdinalVector inducible_parts_added, inducible_parts_removed;
    internal_determine_inducible_parts(entity, add_parts, parts_removed, inducible_parts_added, inducible_parts_removed);
    internal_propagate_induced_part_changes_to_downward_connected_entities(entity, inducible_parts_added, inducible_parts_removed);
}

void BulkData::internal_change_entity_parts(
  Entity entity ,
  const OrdinalVector& add_parts ,
  const OrdinalVector& remove_parts,
  OrdinalVector& scratchOrdinalVec,
  OrdinalVector& scratchSpace)
{
    require_ok_to_modify();
//    m_modSummary.track_change_entity_parts(entity, add_parts, remove_parts);

    Bucket * const bucket_old = bucket_ptr(entity);
    bool needToChangeParts = bucket_old == NULL
            || !bucket_old->member_all(add_parts)
            || bucket_old->member_any(remove_parts);
    if(needToChangeParts)
    {
        if (!add_parts.empty())
        {
            notifier.notify_entity_parts_added(entity, add_parts);
        }
        if (!remove_parts.empty())
        {
            notifier.notify_entity_parts_removed(entity, remove_parts);
        }
        OrdinalVector parts_removed;
        internal_change_entity_parts_without_propogating_to_downward_connected_entities(entity, add_parts, remove_parts, parts_removed, scratchOrdinalVec, scratchSpace);
        internal_determine_inducible_parts_and_propagate_to_downward_connected_entities(entity, add_parts, parts_removed);
    }
}

void BulkData::internal_move_entity_to_new_bucket(stk::mesh::Entity entity, const OrdinalVector &newBucketPartList, OrdinalVector& /*scratchSpace*/)
{
    const MeshIndex &meshIndex = mesh_index(entity);
    Bucket *bucketOld = meshIndex.bucket;
    if (bucketOld != nullptr)
    {
        bool isEntityCommunicated = (bucketOld->shared() || in_send_ghost(entity) || in_receive_ghost(entity));
        if (!m_meshModification.did_any_shared_entity_change_parts() && isEntityCommunicated)
        {
            m_meshModification.set_shared_entity_changed_parts();
        }

        m_bucket_repository.change_entity_part_membership(meshIndex, newBucketPartList);
    }
    else
    {
        EntityRank rank = entity_rank(entity);
        m_bucket_repository.add_entity_with_part_memberships(entity, rank, newBucketPartList);
    }

    notifier.notify_local_buckets_changed(entity_rank(entity));

    mark_entity_and_upward_related_entities_as_modified(entity);
}

void BulkData::internal_fill_new_part_list_and_removed_part_list(stk::mesh::Entity entity,
                                                                 const OrdinalVector & add_parts,
                                                                 const OrdinalVector & remove_parts,
                                                                 OrdinalVector &newBucketPartList,
                                                                 OrdinalVector &parts_removed)
{
    newBucketPartList.clear();
    parts_removed.clear();

    Bucket *bucket_old = bucket_ptr( entity );
    if(bucket_old != NULL)
    {
        const std::pair<const unsigned *, const unsigned*> oldEntityPartMembership =
                bucket_old->superset_part_ordinals();

        const size_t num_bucket_parts = bucket_old->supersets().size();
        newBucketPartList.reserve(num_bucket_parts + add_parts.size());
        newBucketPartList.assign(oldEntityPartMembership.first, oldEntityPartMembership.second);

        if(!remove_parts.empty())
        {
            parts_removed.reserve(remove_parts.size());
            impl::filter_out(newBucketPartList, remove_parts, parts_removed);
        }
    }
    else
    {
        newBucketPartList.reserve(add_parts.size());
    }

    if(!add_parts.empty())
    {
        impl::merge_in(newBucketPartList, add_parts);
    }

    if ( newBucketPartList.empty() )
    {
      newBucketPartList.push_back( m_mesh_meta_data.universal_part().mesh_meta_data_ordinal() );
    }
    stk::util::sort_and_unique(newBucketPartList);
}

void BulkData::internal_adjust_closure_count(Entity entity,
                                             const OrdinalVector& add_parts,
                                             const OrdinalVector& remove_parts)
{
    Bucket *bucket = bucket_ptr( entity );

    const unsigned locally_owned_ordinal = m_mesh_meta_data.locally_owned_part().mesh_meta_data_ordinal();

    bool isInOwnedPart = bucket && bucket->owned();

    bool add_to_locally_owned_part = !isInOwnedPart
                                     && contains_ordinal(add_parts.begin(),
                                                         add_parts.end(),
                                                         locally_owned_ordinal);

    bool remove_from_locally_owned_part = isInOwnedPart
                                          && contains_ordinal(remove_parts.begin(),
                                                              remove_parts.end(),
                                                              locally_owned_ordinal);
    if(add_to_locally_owned_part)
    {
        unprotect_orphaned_node(entity);
        int incrementClosureCount = 1;
        internal_adjust_entity_and_downward_connectivity_closure_count(entity, bucket, incrementClosureCount);
    }
    else if(remove_from_locally_owned_part)
    {
        int decrementClosureCount = -1;
        internal_adjust_entity_and_downward_connectivity_closure_count(entity, bucket, decrementClosureCount);
    }
}

void BulkData::internal_fill_parts_to_actually_remove(const OrdinalVector & parts_to_remove_assuming_not_induced_from_other_entities,
                                                      OrdinalVector &partsThatShouldStillBeInduced,
                                                      OrdinalVector &parts_to_actually_remove)
{
    for(size_t k = 0; k < parts_to_remove_assuming_not_induced_from_other_entities.size(); k++)
    {
        if(!contains_ordinal(partsThatShouldStillBeInduced.begin(),
                             partsThatShouldStillBeInduced.end(),
                             parts_to_remove_assuming_not_induced_from_other_entities[k]))
        {
            stk::util::insert_keep_sorted_and_unique(parts_to_remove_assuming_not_induced_from_other_entities[k], parts_to_actually_remove);
        }
    }
}

//----------------------------------------------------------------------
// Deduce propagation of part membership changes to a 'from' entity
// to the related 'to' entities.  There can be both additions and
// removals.

void BulkData::internal_propagate_induced_part_changes_to_downward_connected_entities(
  Entity entity,
  const OrdinalVector & addParts,
  const OrdinalVector & parts_to_remove_assuming_not_induced_from_other_entities )
{
    m_check_invalid_rels = false;

    const EntityRank erank = entity_rank(entity);

    OrdinalVector partsThatShouldStillBeInduced;
    OrdinalVector parts_to_actually_remove, emptyParts;
    OrdinalVector scratchOrdinalVec, scratchSpace;
    for (EntityRank irank = stk::topology::BEGIN_RANK; irank < erank; ++irank)
    {
        size_t num_rels = num_connectivity(entity, irank);
        if (num_rels > 0) {
            Entity const *rel_entities = begin(entity, irank);
            for (size_t j = 0; j < num_rels; ++j)
            {
                Entity e_to = rel_entities[j];

                if (e_to == Entity::InvalidEntity)
                {
                    continue;
                }

                parts_to_actually_remove.clear();

                const bool remote_changes_needed = !( parallel_size() == 1 || !bucket(e_to).shared() );
                if (remote_changes_needed)
                {
                    Bucket *bucket_old = bucket_ptr(e_to);
                    if ( !m_meshModification.did_any_shared_entity_change_parts() && bucket_old && (bucket_old->shared()  || this->in_send_ghost(entity_key(entity)) || this->in_receive_ghost(entity_key(entity)) ))
                    {
                        m_meshModification.set_shared_entity_changed_parts();
                    }

                    // Don't remove parts until modification_end to avoid losing field data with bucket churn.
                    mark_entity_and_upward_related_entities_as_modified(e_to);
                }
                else
                {
                    if(!parts_to_remove_assuming_not_induced_from_other_entities.empty())
                    {
                        partsThatShouldStillBeInduced.clear();
                        internal_insert_all_parts_induced_from_higher_rank_entities_to_vector(entity,
                                                                                              e_to,
                                                                                              partsThatShouldStillBeInduced);
                        internal_fill_parts_to_actually_remove(parts_to_remove_assuming_not_induced_from_other_entities,
                                                               partsThatShouldStillBeInduced,
                                                               parts_to_actually_remove);
                    }
                }
                m_modSummary.track_induced_parts(entity, e_to, addParts, parts_to_actually_remove);
                internal_change_entity_parts( e_to , addParts , parts_to_actually_remove, scratchOrdinalVec, scratchSpace );
            }
        }
    }
    m_check_invalid_rels = true;
}

void BulkData::internal_insert_all_parts_induced_from_higher_rank_entities_to_vector(stk::mesh::Entity entity,
                                                                                     stk::mesh::Entity e_to,
                                                                                     OrdinalVector &to_add)
{
    EntityRank e_to_rank = entity_rank(e_to);
    EntityRank start_rank = static_cast<EntityRank>(e_to_rank + 1);
    EntityRank end_rank = static_cast<EntityRank>(m_mesh_meta_data.entity_rank_count());
    for(EntityRank to_rel_rank_i = start_rank; to_rel_rank_i < end_rank; ++to_rel_rank_i)
    {
        int num_upward_rels = num_connectivity(e_to, to_rel_rank_i);
        Entity const* upward_rel_entities = begin(e_to, to_rel_rank_i);

        for (int k = 0; k < num_upward_rels; ++k)
        {
            if (entity != upward_rel_entities[k])  // Already did this entity
            {
                // Relation from to_rel->entity() to e_to
                impl::get_part_ordinals_to_induce_on_lower_ranks(*this, upward_rel_entities[k], e_to_rank, to_add );
            }
        }
    }
}

// TODO Change the methods below to requirements (private, const invariant checkers)

// Do not allow any of the induced part memberships to explicitly
// appear in the add or remove parts lists.
// 1) Intersection part
// 3) Part that does not match the entity rank.

void BulkData::internal_verify_change_parts( const MetaData   & meta ,
                                             const Entity entity ,
                                             const PartVector & parts ) const
{
  const std::vector<std::string> & rank_names = meta.entity_rank_names();
  const EntityRank erank = entity_rank(entity);

  bool ok = true ;
  std::ostringstream msg ;

  for (PartVector::const_iterator
        i = parts.begin() ; i != parts.end() ; ++i ) {

    const Part & p = **i;
    const unsigned part_rank = p.primary_entity_rank();

    bool rank_ok = (part_rank == InvalidEntityRank || part_rank == erank);

    if ( !rank_ok ) {
      if ( ok ) {
        msg << "change parts for entity " << identifier(entity);
        msg << " , { " ;
      }
      else {
        msg << " , " ;
      }
      ok = false ;

      msg << p.name() << "[" ;
      if ( part_rank < rank_names.size() ) {
        msg << rank_names[ part_rank ];
      }
      else {
        msg << part_rank ;
      }
      msg << "] " ;
      if ( !rank_ok )         { msg << "is_bad_rank " ; }
    }
  }

  ThrowErrorMsgIf( !ok, msg.str() << "}" );
}

void BulkData::sortNodesIfNeeded(std::vector<stk::mesh::EntityKey>& nodes)
{
    std::sort(nodes.begin(),nodes.end());
}

//----------------------------------------------------------------------
void BulkData::markEntitiesForResolvingSharingInfoUsingNodes(stk::mesh::EntityRank entityRank, bool onlyConsiderSoloSides, std::vector<shared_entity_type>& shared_entities)
{
    const stk::mesh::BucketVector& entity_buckets = this->buckets(entityRank);
    const bool add_node_sharing_called = this->internal_add_node_sharing_called();

    for(size_t bucketIndex = 0; bucketIndex < entity_buckets.size(); bucketIndex++)
    {
        const stk::mesh::Bucket& bucket = *entity_buckets[bucketIndex];
        stk::topology topology = bucket.topology();
        for(size_t entityIndex = 0; entityIndex < bucket.size(); entityIndex++)
        {
            Entity entity = bucket[entityIndex];
            const unsigned num_nodes_on_entity = bucket.num_nodes(entityIndex);

            if (!add_node_sharing_called && this->state(entity) == stk::mesh::Unchanged)
            {
              // No nodes newly shared and entity has not had nodes added, so entity cannot become shared.
              continue;
            }

            const bool isSide = entityRank==this->mesh_meta_data().side_rank();
            const bool isSoloSide = isSide && (this->num_elements(entity)==0);
            if (isSide && onlyConsiderSoloSides && !isSoloSide) {
                continue;
            }

            if ( num_nodes_on_entity > 1 )
            {
                if(owned_closure(entity))
                {
                    Entity const * nodes = bucket.begin_nodes(entityIndex);

                    //do we need to do some sorting operation here?
                    //sort entity nodes into lexicographical smallest permutation?


                    bool shared_entity = true;
                    for(size_t n = 0; n < num_nodes_on_entity; ++n)
                    {
                        Entity node = nodes[n];
                        shared_entity = shared_entity && (this->bucket(node).shared() || (this->internal_is_entity_marked(node) == BulkData::IS_SHARED));
                    }

                    if(shared_entity)
                    {
                        shared_entity_type sentity(this->entity_key(entity), entity, topology);
                        sentity.nodes.resize(num_nodes_on_entity);
                        for(size_t n = 0; n < num_nodes_on_entity; ++n)
                        {
                            sentity.nodes[n]=this->entity_key(nodes[n]);
                        }
                        //Sort will have to go away
                        this->sortNodesIfNeeded(sentity.nodes);
                        shared_entities.push_back(sentity);
                        this->internal_mark_entity(entity, BulkData::POSSIBLY_SHARED);
                    }
                }
            }
        }
    }
}

void BulkData::gather_shared_nodes(std::vector<Entity> & shared_nodes)
{
    const stk::mesh::BucketVector & node_buckets = this->buckets(stk::topology::NODE_RANK);

    for(size_t nodeIndex = 0; nodeIndex < node_buckets.size(); ++nodeIndex)
    {
        const stk::mesh::Bucket & bucket = *node_buckets[nodeIndex];
        for(size_t entityIndex = 0; entityIndex < bucket.size(); ++entityIndex)
        {
            Entity node = bucket[entityIndex];
            if (this->internal_is_entity_marked(node) == BulkData::IS_SHARED)
            {
                shared_nodes.push_back(node);
            }
        }
    }
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

void BulkData::unpack_not_owned_verify_compare_comm_info( CommBuffer&            buf,
                                                Entity                 entity,
                                                EntityKey &            recv_entity_key,
                                                int       &            recv_owner_rank,
                                                unsigned  &            recv_comm_count,
                                                PartVector&    recv_parts,
                                                std::vector<Relation>& recv_relations,
                                                std::vector<int>    &  recv_comm,
                                                bool&                  bad_comm)
{
  if (!in_shared(entity) && !in_ghost(aura_ghosting(), entity)) {
    //don't do consistency-check for custom-ghosts
    bad_comm = false;
    return;
  }

  const EntityKey key = entity_key(entity);
  const PairIterEntityComm ec = internal_entity_comm_map(key);
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
  if ( in_shared( key ) ) {
    // only packed shared size, so only compare with shared here
    bad_comm = ec_idx_shared.size() != recv_comm.size();
    if ( ! bad_comm ) {
      size_t j = 0 ;
      for ( ; j < ec_idx_shared.size() &&
              ec[ec_idx_shared[j]].ghost_id == BulkData::SHARED &&
              ec[ec_idx_shared[j]].proc   == recv_comm[j] ; ++j );
      bad_comm = j != ec_idx_shared.size() ;

      // unpack count of additional ghosts
      bad_comm = bad_comm || check_tag(*this, buf, PACK_TAG_GHOST_COUNT_AFTER_SHARED);
      buf.unpack<unsigned>( ghost_after_shared_count);
    }
  }

  if ( ! bad_comm ) {

    if (ghost_after_shared_count) {
      bad_comm = bad_comm || check_tag(*this, buf, PACK_TAG_ENTITY_GHOST);
      unpack_entity_info( buf , *this ,
                          recv_entity_key , recv_owner_rank ,
                          recv_parts , recv_relations );

      bad_comm = bad_comm || check_tag(*this, buf, PACK_TAG_GHOST_COUNT);
      buf.unpack<unsigned>(recv_comm_count);
      recv_comm.resize( recv_comm_count);
      buf.unpack<int>( recv_comm.data() , recv_comm_count);
    }

    if ( !in_shared( key ) || ghost_after_shared_count) {
        // recv_comm contains ghost_ids for ghosted entities
        bad_comm = !impl::all_ghost_ids_are_found_in_comm_data(ec, parallel_owner_rank(entity), recv_comm);
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

        ThrowAssertMsg((rels_itr != rels_end && ords_itr == nullptr) == false, "Relations found without ordinals");


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

std::string bool_str(bool flag)
{
  return flag ? "true" : "false";
}

void BulkData::unpack_not_owned_verify_report_errors(Entity entity,
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
  const int p_rank = parallel_rank();

  const Ordinal bucketOrdinal = bucket_ordinal(entity);
  const EntityRank erank = entity_rank(entity);
  const EntityKey key = entity_key(entity);

  error_log << __FILE__ << ":" << __LINE__ << ": ";
  error_log << "P" << p_rank << ": " ;
  error_log << key;
  error_log << " owner(P" << parallel_owner_rank(entity) << ") shared: " << bool_str(bucket(entity).shared()) << " in aura: " << bool_str(bucket(entity).in_aura()) << " ";

  if ( bad_key || bad_own ) {
    error_log << " != received " ;
    error_log << recv_entity_key;
    error_log << " owner(" << recv_owner_rank
              << ")" << std::endl ;
  }
  else if ( bad_comm ) {
    const PairIterEntityComm ec = internal_entity_comm_map(key);
    if ( in_shared( key ) ) {
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
        error_log << ",p" << parallel_owner_rank(entity);
        error_log << ")" ;
      }
      error_log << " )" << std::endl ;
    }
  }
  else if ( bad_part ) {
    error_log << " Comparing parts from this processor(" << parallel_rank() << ") against processor (" << recv_owner_rank << ")" << std::endl;

    std::set<std::string> thisProcExtraParts;
    std::set<std::string> otherProcExtraParts;
    impl::fill_part_list_differences(*this, entity, recv_parts, thisProcExtraParts, otherProcExtraParts);

    if ( !thisProcExtraParts.empty() )
    {
        error_log << "\tParts on this proc, not on other proc:" << std::endl;
        std::set<std::string>::iterator iter = thisProcExtraParts.begin();
        for (;iter!=thisProcExtraParts.end();++iter)
        {
            error_log << "\t\t" << *iter << std::endl;
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
    const Bucket & entityBucket = bucket(entity);
    for (EntityRank irank = stk::topology::BEGIN_RANK;
         irank < erank; ++irank)
    {
      error_log << " " << irank << ": ";
      Entity const *ir_itr = entityBucket.begin(bucketOrdinal, irank);
      Entity const *ir_end = entityBucket.end(bucketOrdinal, irank);
      for ( ; ir_itr != ir_end; ++ir_itr ) {
        error_log << identifier(*ir_itr)<<" " ;
        if (irank != stk::topology::NODE_RANK) { 
          Entity const * nodes_begin = begin_nodes(*ir_itr);
          Entity const * nodes_end   = end_nodes(*ir_itr);
          error_log << "nodes (";
          for (Entity const* nodeId = nodes_begin; nodeId != nodes_end; ++nodeId)
          {
              error_log << identifier(*nodeId) << ", ";
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
      error_log << identifier(jr->entity()) << " ";
      if (jr->entity_rank() != stk::topology::NODE_RANK) {
        Entity const * nodes_begin = begin_nodes(jr->entity());
        Entity const * nodes_end   = end_nodes(jr->entity());
        error_log << " nodes (";
        for (Entity const* nodeId = nodes_begin; nodeId != nodes_end; ++nodeId)
        {
            error_log << identifier(*nodeId) << ", ";
        }
        error_log << ") ";
      }
    }
    error_log << " )" << std::endl ;
  }
}


//----------------------------------------------------------------------------
// Unpacking all of my not-owned entities.
bool BulkData::unpack_not_owned_verify( CommSparse & commSparse , std::ostream & error_log )
{
  const int               p_rank = parallel_rank();
  const EntityCommListInfoVector & entity_comm = internal_comm_list();

  bool result = true ;

  EntityKey             recv_entity_key ;
  int                   recv_owner_rank = 0 ;
  unsigned              recv_comm_count = 0 ;
  PartVector    recv_parts ;
  std::vector<Relation> recv_relations ;
  std::vector<int>      recv_comm ;

  for ( EntityCommListInfoVector::const_iterator
        i = entity_comm.begin() ; i != entity_comm.end() ; ++i ) {

    EntityKey key = i->key;
    Entity entity = i->entity;
    ThrowRequire( entity_key(entity) == key );


    if ( parallel_owner_rank(i->entity) != p_rank ) {

      bool broken_tag = false;
      CommBuffer & buf = commSparse.recv_buffer( parallel_owner_rank(i->entity));

      broken_tag = broken_tag || check_tag(*this, buf, PACK_TAG_ENTITY_SHARED, PACK_TAG_ENTITY_GHOST);
      if(!broken_tag)
      {
          unpack_entity_info( buf , *this,
                              recv_entity_key , recv_owner_rank ,
                              recv_parts , recv_relations );

          if (in_shared(key)) {
              broken_tag = broken_tag || check_tag(*this, buf, PACK_TAG_SHARED_COUNT);
          }
          else {
              broken_tag = broken_tag || check_tag(*this, buf, PACK_TAG_GHOST_COUNT);
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

      const bool bad_key = key                              != recv_entity_key ;
      const bool bad_own = parallel_owner_rank(entity) != recv_owner_rank ;
      bool bad_part = false ;
      bool bad_rel  = false ;
      bool bad_comm = false ;

      bool broken = broken_tag || bad_key || bad_own;

      // Compare communication information:

      if ( ! broken ) {
        unpack_not_owned_verify_compare_comm_info( buf,
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
        unpack_not_owned_verify_compare_parts(*this,
                                              entity,
                                              recv_parts,
                                              bad_part);
        broken = bad_part;
      }

      // Compare the closure relations:
      if ( ! broken )
      {
        unpack_not_owned_verify_compare_closure_relations( *this,
                                                           entity,
                                                           recv_relations,
                                                           bad_rel );
        broken = bad_rel;

      }

      // The rest of this code is just error handling
      if ( broken ) {
        unpack_not_owned_verify_report_errors(entity,
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

void BulkData::pack_owned_verify( CommSparse & commSparse )
{
  const EntityCommListInfoVector & entity_comm = internal_comm_list();
  const int p_rank = commSparse.parallel_rank();

  for ( EntityCommListInfoVector::const_iterator
        i = entity_comm.begin() ; i != entity_comm.end() ; ++i ) {

    if ( parallel_owner_rank(i->entity) == p_rank ) {

      std::vector<int> share_procs ;
      std::vector<int> ghost_procs ;

      const PairIterEntityComm comm = internal_entity_comm_map(i->key);

      for ( size_t j = 0 ; j < comm.size() ; ++j ) {
        if ( comm[j].ghost_id == stk::mesh::BulkData::SHARED ) {
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

        pack_entity_info(*this, buf , i->entity );

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
        pack_entity_info(*this, buf , i->entity );

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

void printConnectivityOfRank(BulkData& M, Entity entity, stk::topology::rank_t connectedRank, std::ostream & error_log)
{
    if (M.is_valid(entity)) {
    error_log << connectedRank << "-connectivity(";
    const Entity* connectedEntities = M.begin(entity, connectedRank);
    unsigned numConnected = M.num_connectivity(entity, connectedRank);
    for(unsigned i=0; i<numConnected; ++i) {
        if (M.is_valid(connectedEntities[i])) {
            error_log<<"{"<<M.identifier(connectedEntities[i])<<",topo="<<M.bucket(connectedEntities[i]).topology()
                  <<",owned="<<M.bucket(connectedEntities[i]).owned()<<",shared="<<M.bucket(connectedEntities[i]).shared()
                  <<",in_aura="<<M.bucket(connectedEntities[i]).in_aura()
                  <<",custom-recv-ghost="<<M.in_receive_custom_ghost(M.entity_key(connectedEntities[i]))
                  <<"}";
        }
        else {
            error_log << "{invalid entity!}";
        }
    }
    error_log<<"), ";
    }
    else {
        error_log << "invalid entity!";
    }
}

bool BulkData::verify_parallel_attributes_for_bucket( Bucket const& bucket, std::ostream & error_log)
{
  const int p_rank = parallel_rank();

  bool result = true;

  Part & owns_part = m_mesh_meta_data.locally_owned_part();
  Part & shares_part = m_mesh_meta_data.globally_shared_part();

  const bool has_owns_part   = has_superset( bucket , owns_part );
  const bool has_shares_part = has_superset( bucket , shares_part );

  const Bucket::iterator j_end = bucket.end();
  Bucket::iterator j           = bucket.begin();

  while ( j != j_end ) {
    Entity entity = *j ; ++j ;

    bool this_result = true;

    const EntityKey key = entity_key(entity);
    const int      p_owner    = parallel_owner_rank(entity);
    const bool     ordered    = impl::is_comm_ordered(internal_entity_comm_map(key));
    const bool     shares     = in_shared( key );
    const bool     recv_aura = in_receive_ghost( aura_ghosting(), key );
    const bool     recv_any_ghost = in_receive_ghost( key );
    const bool     custom_ghost = recv_any_ghost && !recv_aura;
    const bool     send_ghost = in_send_ghost( key );
    const bool     ownedClosure = owned_closure(entity);

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
                << "in comm-list: " << (shares?"true":"false") << ", entity key " << entity_key(entity) <<" "<<this->bucket(entity).topology() << std::endl;
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
      error_log << "problem: entity "<<key<<" (with topology "<<bucket.topology()<<") is both recv aura ghost and in owned_closure;"<<std::endl;
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
      comm_shared_procs(entity_key(entity),shared_procs);
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
      error_log << "P" << parallel_rank() << ": " << " entity " << this->entity_key(entity) << " " << this->bucket(entity).topology();
      error_log << " details: owner(" << p_owner<<"), shared=" << (bucket.shared() ? "true" : "false");
      error_log << ", aura=" << (bucket.in_aura() ? "true" : "false");
      error_log <<", custom-recv-ghost="<<this->in_receive_custom_ghost(this->entity_key(entity));
      error_log << std::endl;

      for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<mesh_meta_data().entity_rank_count(); rank++)
      {
          printConnectivityOfRank(*this, entity, rank, error_log);
      }

      error_log<<"comm(";
      PairIterEntityComm ip = internal_entity_comm_map(entity_key(entity));
      for ( ; ! ip.empty() ; ++ip ) {
        error_log << " ghost_id=" << ip->ghost_id << ":proc=" << ip->proc ;
      }
      error_log << " )" << std::endl ;
    }
  }

  return result;
}

bool BulkData::verify_parallel_attributes( std::ostream & error_log )
{
  bool result = true ;

  const EntityRank entityRankEnd = static_cast<EntityRank>(m_mesh_meta_data.entity_rank_count());

  for ( EntityRank rank = stk::topology::NODE_RANK ; rank < entityRankEnd ; ++rank ) {
    const BucketVector & all_buckets = buckets(rank);

    for(const Bucket* bucketptr : all_buckets)
    {
      result = result && verify_parallel_attributes_for_bucket(*bucketptr, error_log);
    }
  }

  bool isGloballyConsistentCommList = impl::is_comm_list_globally_consistent(*this, this->m_entity_comm_list, error_log);
  result = result && isGloballyConsistentCommList;

  return result ;
}

bool BulkData::comm_mesh_verify_parallel_consistency(std::ostream & error_log )
{
  int verified_ok = 1 ;

  // Verify consistency of parallel attributes

  verified_ok = verify_parallel_attributes( error_log );

  if (parallel_size() > 1) {
    all_reduce( parallel() , ReduceMin<1>( & verified_ok ) );
  }

  // Verify entities against owner.

  if ( verified_ok ) {
    CommSparse comm( parallel() );

    pack_owned_verify( comm );

    comm.allocate_buffers();

    pack_owned_verify(comm);

    comm.communicate();

    verified_ok = unpack_not_owned_verify( comm , error_log );

    if (parallel_size() > 1) {
      all_reduce( parallel() , ReduceMin<1>( & verified_ok ) );
    }
  }

  return verified_ok == 1 ;
}

// Enforce that shared entities must be in the owned closure:

void BulkData::destroy_dependent_ghosts( Entity entity, EntityProcVec& entitiesToRemoveFromSharing )
{
  EntityRank entity_rank = this->entity_rank(entity);
  const EntityRank end_rank = static_cast<EntityRank>(this->mesh_meta_data().entity_rank_count());

  for (EntityRank irank = static_cast<EntityRank>(end_rank - 1); irank > entity_rank; --irank)
  {
    int num_rels = this->num_connectivity(entity, irank);
    const Entity* rels     = this->begin(entity, irank);

    for (int r = num_rels - 1; r >= 0; --r)
    {
      Entity e = rels[r];

      bool upwardRelationOfEntityIsInClosure = this->owned_closure(e);
      ThrowRequireMsg( !upwardRelationOfEntityIsInClosure, this->entity_rank(e) << " with id " << this->identifier(e) << " should not be in closure." );

      // Recursion
      if (this->is_valid(e) && this->bucket(e).in_aura())
      {
          this->destroy_dependent_ghosts( e, entitiesToRemoveFromSharing );
      }
    }
  }

  const bool successfully_destroyed_entity = this->destroy_entity(entity);
  if (!successfully_destroyed_entity)
  {
      std::vector<int> sharing_procs;
      comm_shared_procs(entity_key(entity), sharing_procs);
      for(int p : sharing_procs) {
          entitiesToRemoveFromSharing.emplace_back(entity, p);
      }
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

void BulkData::delete_shared_entities_which_are_no_longer_in_owned_closure(EntityProcVec& entitiesToRemoveFromSharing)
{
  for ( EntityCommListInfoVector::const_reverse_iterator
        i =  internal_comm_list().rbegin() ;
        i != internal_comm_list().rend() ; ++i)
  {
    Entity entity = i->entity;
    if (is_valid(entity) && !owned_closure(entity)) {
      if ( in_shared(entity) )
      {
        destroy_dependent_ghosts( entity, entitiesToRemoveFromSharing );
      }
    }
  }
}

void BulkData::remove_entities_from_sharing(const EntityProcVec& entitiesToRemoveFromSharing, stk::mesh::EntityVector & entitiesNoLongerShared)
{
  entitiesNoLongerShared.clear();
  OrdinalVector scratchOrdinalVec, scratchSpace;
  for(const EntityProc& entityAndProc : entitiesToRemoveFromSharing) {
      EntityKey key = this->entity_key(entityAndProc.first);
      this->entity_comm_map_erase(key,EntityCommInfo(BulkData::SHARED, entityAndProc.second));
      if (!this->in_shared(key)) {
          entitiesNoLongerShared.push_back(entityAndProc.first);
          this->internal_change_entity_parts(entityAndProc.first,{},{this->mesh_meta_data().globally_shared_part().mesh_meta_data_ordinal()}, scratchOrdinalVec, scratchSpace);
          this->internal_mark_entity(entityAndProc.first, NOT_SHARED);
      }
  }
  stk::util::sort_and_unique(entitiesNoLongerShared);
}

namespace
{
bool is_node_connected_to_active_element_locally(const stk::mesh::BulkData &mesh, stk::mesh::Entity node, const stk::mesh::Part &activePart)
{
    bool activeNode = false;
    const int numElements = mesh.num_elements(node);
    const stk::mesh::Entity * elements = mesh.begin_elements(node);
    for (int elementI=0 ; elementI<numElements ; ++elementI)
    {
        stk::mesh::Entity connectedElement = elements[elementI];
        stk::mesh::Bucket &connectedElementBucket = mesh.bucket(connectedElement);
        if (connectedElementBucket.owned() && connectedElementBucket.member(activePart))
        {
            activeNode = true;
            break;
        }
    }
    return activeNode;
}
} //emtpy namespace

void BulkData::delete_sides_on_all_procs(const stk::mesh::EntityVector& deletedSides)
{
    stk::CommSparse comm(this->parallel());
    std::vector<int> procs;
    for(int phase = 0; phase < 2; ++phase)
    {
        for(size_t i = 0; i < deletedSides.size(); ++i)
        {
            stk::mesh::Entity side = deletedSides[i];
            stk::mesh::EntityKey key = this->entity_key(side);
            const bool is_comm_entity_and_locally_owned = this->parallel_owner_rank(side) == this->parallel_rank();
            if(is_comm_entity_and_locally_owned)
            {
                procs.clear();
                for ( PairIterEntityComm ec = internal_entity_comm_map(key); ! ec.empty() ; ++ec )
                {
                    procs.push_back( ec->proc );
                }
                stk::util::sort_and_unique(procs);

                for(size_t proc_index = 0; proc_index < procs.size(); ++proc_index)
                {
                    const int proc = procs[proc_index];
                    stk::CommBuffer & buf = comm.send_buffer(proc);
                    buf.pack<stk::mesh::EntityKey>(entity_key(side));
                }

                if(phase == 1)
                {
                    this->entity_comm_map_clear(this->entity_key(side));
                }
            }
        }

        if(phase == 0)
        {
            comm.allocate_buffers();
        }
        else
        {
            comm.communicate();
        }
    }

    stk::mesh::EntityVector recvSidesToDelete;
    for(int p = 0; p < this->parallel_size(); ++p)
    {
        stk::CommBuffer & buf = comm.recv_buffer(p);
        while(buf.remaining())
        {
            stk::mesh::EntityKey key;
            buf.unpack<stk::mesh::EntityKey>(key);
            this->entity_comm_map_clear(key);
            stk::mesh::Entity side = this->get_entity(key);
            if(this->is_valid(side))
            {
                recvSidesToDelete.push_back(side);
            }
        }
    }

    stk::mesh::impl::delete_entities_and_upward_relations(*this, recvSidesToDelete);
    this->update_comm_list_based_on_changes_in_comm_map();
}

void BulkData::set_shared_owned_parts_and_ownership_on_comm_data(const std::vector<sharing_info>& shared_modified)
{
    stk::mesh::OrdinalVector shared_part, owned_part, empty;
    shared_part.push_back(m_mesh_meta_data.globally_shared_part().mesh_meta_data_ordinal());
    owned_part.push_back(m_mesh_meta_data.locally_owned_part().mesh_meta_data_ordinal());

    stk::mesh::EntityVector modified_entities(shared_modified.size());
    OrdinalVector scratchOrdinalVec, scratchSpace;
    for(size_t i = 0; i < shared_modified.size(); ++i)
    {
        stk::mesh::Entity entity = shared_modified[i].m_entity;
        int sharing_proc = shared_modified[i].m_sharing_proc;
        entity_comm_map_insert(entity, stk::mesh::EntityCommInfo(stk::mesh::BulkData::SHARED, sharing_proc));
        int old_owner = parallel_owner_rank(entity);
        int owning_proc = shared_modified[i].m_owner;
        if(old_owner != owning_proc)
        {
            internal_set_owner(entity, owning_proc);
            internal_change_entity_parts(entity, shared_part /*add*/, owned_part /*remove*/, scratchOrdinalVec, scratchSpace);
        }
        else
        {
            internal_change_entity_parts(entity, shared_part /*add*/, empty /*remove*/, scratchOrdinalVec, scratchSpace);
        }
        modified_entities[i] = entity;
    }

    stk::util::sort_and_unique(modified_entities, stk::mesh::EntityLess(*this));

    add_comm_list_entries_for_entities(modified_entities);
}


void BulkData::make_mesh_parallel_consistent_after_element_death(const std::vector<sharing_info> &shared_modified,
                                                                 const stk::mesh::EntityVector &deletedSides,
                                                                 stk::mesh::ElemElemGraph &elementGraph,
                                                                 const stk::mesh::EntityVector &killedElements,
                                                                 stk::mesh::Part &activePart)
{
    if(!in_synchronized_state())
    {
        for(size_t i = 0; i < deletedSides.size(); ++i)
        {
            ThrowAssertMsg(entity_rank(deletedSides[i]) == mesh_meta_data().side_rank(), "ERROR, modification_end_for_face_deletion only handles faces");
        }

        ThrowAssertMsg(stk::mesh::impl::check_for_connected_nodes(*this)==0, "BulkData::modification_end ERROR, all entities with rank higher than node are required to have connected nodes.");

        ThrowAssertMsg(add_fmwk_data() || stk::mesh::impl::check_no_shared_elements_or_higher(*this)==0,
                "BulkData::modification_end ERROR, Sharing of entities with rank ELEMENT_RANK or higher is not allowed.");

        if(parallel_size() > 1)
        {
            delete_sides_on_all_procs(deletedSides);
            set_shared_owned_parts_and_ownership_on_comm_data(shared_modified);
        }

        de_induce_parts_from_nodes(killedElements, activePart);
        remove_boundary_faces_from_part(elementGraph, killedElements, activePart);

        if(parallel_size() > 1)
        {
            bool connectFacesToPreexistingGhosts = true;
            internal_resolve_sharing_and_ghosting_for_sides(connectFacesToPreexistingGhosts);
        }

        m_modSummary.write_summary(m_meshModification.synchronized_count(), false);

        if(parallel_size() > 1)
            check_mesh_consistency();

        internal_finish_modification_end(impl::MeshModification::MOD_END_SORT);
    }
}

void BulkData::internal_resolve_sharing_and_ghosting_for_sides(bool connectFacesToPreexistingGhosts)
{
    stk::mesh::EntityVector entitiesNoLongerShared;
    m_meshModification.internal_resolve_shared_modify_delete(entitiesNoLongerShared);
    internal_resolve_shared_part_membership_for_element_death();
    if(is_automatic_aura_on())
        resolve_incremental_ghosting_for_entity_creation_or_skin_mesh(mesh_meta_data().side_rank(),
                                                                      mesh_meta_data().universal_part(),
                                                                      connectFacesToPreexistingGhosts);
}

void BulkData::make_mesh_parallel_consistent_after_skinning(const std::vector<sharing_info>& sharedModified)
{
    if(!in_synchronized_state())
    {
        ThrowAssertMsg(stk::mesh::impl::check_for_connected_nodes(*this)==0, "BulkData::modification_end ERROR, all entities with rank higher than node are required to have connected nodes.");

        if(parallel_size() > 1)
        {
            set_shared_owned_parts_and_ownership_on_comm_data(sharedModified);

            bool connectFacesToPreexistingGhosts = false;
            internal_resolve_sharing_and_ghosting_for_sides(connectFacesToPreexistingGhosts);
        }

        m_modSummary.write_summary(m_meshModification.synchronized_count(), false);

        if (parallel_size() > 1)
          check_mesh_consistency();

        internal_finish_modification_end(impl::MeshModification::MOD_END_SORT);
    }
}

void BulkData::remove_boundary_faces_from_part(stk::mesh::ElemElemGraph &graph, const stk::mesh::EntityVector & deactivatedElements, const stk::mesh::Part & activePart)
{
    stk::mesh::EntityVector sidesToRemoveFromPart;
    OrdinalVector scratchOrdinalVec, scratchSpace;
    for (stk::mesh::Entity element : deactivatedElements)
    {
        size_t numSides = this->num_connectivity(element, mesh_meta_data().side_rank());
        const stk::mesh::Entity * sides = this->begin(element, mesh_meta_data().side_rank());
        const stk::mesh::ConnectivityOrdinal * sideOrdinals = this->begin_ordinals(element, mesh_meta_data().side_rank());
        for (size_t sideI=0 ; sideI<numSides ; ++sideI)
        {
            stk::mesh::Entity side = sides[sideI];
            stk::mesh::ConnectivityOrdinal sideOrdinal = sideOrdinals[sideI];
            if(!graph.is_connected_to_other_element_via_side_ordinal(element, sideOrdinal))
            {
                sidesToRemoveFromPart.push_back(side);
            }
            // find if this element is connected to any other element through this sideOrdinal
            // if not, deactivate it.
            // deactiveate the skin of the killed elements
            // We can deactivate all sides at this point because we already deleted the death-created-sides
        }
    }

    const stk::mesh::OrdinalVector rm_parts(1, activePart.mesh_meta_data_ordinal());

    for (stk::mesh::Entity side : sidesToRemoveFromPart)
    {
        this->internal_change_entity_parts(side, {}, rm_parts, scratchOrdinalVec, scratchSpace);
    }

    std::vector<int> commProcs;
    stk::CommSparse comm(this->parallel());
    pack_and_communicate(comm,
        [this,&comm,&sidesToRemoveFromPart,&commProcs]()
        {
            for (stk::mesh::Entity side : sidesToRemoveFromPart)
            {
                const stk::mesh::EntityKey entityKey = this->entity_key(side);
                this->comm_procs(entityKey, commProcs);
                for (int otherProc : commProcs)
                {
                    comm.send_buffer(otherProc).pack<stk::mesh::EntityId>(entityKey.id());
                }
            }
        }
    );
    unpack_communications(comm,
        [this,&comm,&rm_parts, &scratchOrdinalVec, &scratchSpace](int procId)
        {
            stk::mesh::EntityId sideId;
            comm.recv_buffer(procId).unpack<stk::mesh::EntityId>(sideId);
            stk::mesh::Entity side = this->get_entity(mesh_meta_data().side_rank(), sideId);
            ThrowAssertMsg(this->is_valid(side),"Error in communication for de-imprinting the active part on nodes of killed elements in element death!");
            this->internal_change_entity_parts(side, {}, rm_parts, scratchOrdinalVec, scratchSpace);
        }
    );
}

std::ostream &operator<<(std::ostream &out, const stk::mesh::PartVector &partVector)
{
    out << "{ ";
      for(Part* part : partVector) {
        out << part->name() << " ";
      }
    out << "}";
    return out;
}

stk::mesh::EntityVector BulkData::get_nodes_to_deactivate(const stk::mesh::EntityVector & deactivatedElements, const stk::mesh::Part & activePart) const
{
    stk::mesh::EntityVector nodesToDeactivate;

    stk::mesh::EntityVector potentiallyDeactivatedNodes;
    for (stk::mesh::Entity element : deactivatedElements)
    {
        const int numNodes = this->num_nodes(element);
        const stk::mesh::Entity * nodes = this->begin_nodes(element);
        for (int nodeI=0 ; nodeI<numNodes ; ++nodeI)
        {
            potentiallyDeactivatedNodes.push_back(nodes[nodeI]);
        }
    }
    stk::util::sort_and_unique(potentiallyDeactivatedNodes);

    stk::mesh::EntityVector nodesToCommunicate;
    for (stk::mesh::Entity node : potentiallyDeactivatedNodes)
    {
        if (this->bucket(node).owned() || this->bucket(node).shared())
        {
            bool activeNode = is_node_connected_to_active_element_locally(*this, node, activePart);
            if (!activeNode)
            {
                if (this->bucket(node).shared())
                {
                    nodesToCommunicate.push_back(node);
                }
                else
                {
                    nodesToDeactivate.push_back(node);
                }
            }
        }
    }

    std::vector<int> sharedProcs;
    stk::CommSparse inquiryComm(this->parallel());
    pack_and_communicate(inquiryComm,
        [this,&inquiryComm,&nodesToCommunicate,&sharedProcs]()
        {
            for (stk::mesh::Entity node : nodesToCommunicate)
            {
                const stk::mesh::EntityKey nodeKey = this->entity_key(node);
                this->comm_shared_procs(nodeKey,sharedProcs);
                for (int otherProc : sharedProcs)
                {
                    inquiryComm.send_buffer(otherProc).pack<stk::mesh::EntityId>(nodeKey.id());
                }
            }
        }
    );
    stk::mesh::EntityVector incomingNodes;
    unpack_communications(inquiryComm,
        [this,&inquiryComm,&incomingNodes](int procId)
        {
            stk::mesh::EntityId nodeId;
            inquiryComm.recv_buffer(procId).unpack<stk::mesh::EntityId>(nodeId);
            stk::mesh::Entity node = this->get_entity(stk::topology::NODE_RANK, nodeId);
            ThrowAssertMsg(this->is_valid(node),"Error in communication for de-imprinting the active part on nodes of killed elements in element death!");
            incomingNodes.push_back(node);
        }
    );

    std::map<stk::mesh::Entity,bool> nodeToActiveStatusMap;
    stk::CommSparse answerComm(this->parallel());
    pack_and_communicate(answerComm,
        [this,&answerComm,&incomingNodes,&nodeToActiveStatusMap,&activePart]()
        {
            for (stk::mesh::Entity incomingNode : incomingNodes)
            {
                std::vector<int> sharingProcs;
                this->comm_shared_procs(this->entity_key(incomingNode),sharingProcs);
                bool activeStatus = is_node_connected_to_active_element_locally(*this, incomingNode, activePart);
                for (int otherProc : sharingProcs)
                {
                    answerComm.send_buffer(otherProc).pack<stk::mesh::EntityId>(this->identifier(incomingNode));
                    answerComm.send_buffer(otherProc).pack<bool>(activeStatus);
                }
                auto nodeLocationInMap = nodeToActiveStatusMap.find(incomingNode);
                if (nodeLocationInMap == nodeToActiveStatusMap.end())
                {
                    nodeToActiveStatusMap.emplace(incomingNode, activeStatus);
                }
                else
                {
                    nodeLocationInMap->second = nodeLocationInMap->second || activeStatus;
                }
            }
        }
    );

    unpack_communications(answerComm,
        [this,&answerComm,&nodeToActiveStatusMap](int procId)
        {
            stk::mesh::EntityId nodeId;
            answerComm.recv_buffer(procId).unpack<stk::mesh::EntityId>(nodeId);
            bool activeStatus = false;
            answerComm.recv_buffer(procId).unpack<bool>(activeStatus);
            stk::mesh::Entity node = this->get_entity(stk::topology::NODE_RANK,nodeId);
            ThrowAssertMsg(this->is_valid(node),"Error in communication for de-imprinting the active part on nodes of killed elements in element death!");
            auto nodeLocationInMap = nodeToActiveStatusMap.find(node);
            if (nodeLocationInMap == nodeToActiveStatusMap.end())
            {
                nodeToActiveStatusMap.emplace(node, activeStatus);
            }
            else
            {
                nodeLocationInMap->second = nodeLocationInMap->second || activeStatus;
            }
        }
    );

    for (auto nodeActiveStatusPair : nodeToActiveStatusMap)
    {
        stk::mesh::Entity node = nodeActiveStatusPair.first;
        bool nodeIsActiveOnAnyOtherProcessors = nodeActiveStatusPair.second;
        if (!nodeIsActiveOnAnyOtherProcessors)
        {
            nodesToDeactivate.push_back(node);
        }
    }

    return nodesToDeactivate;
}

void BulkData::de_induce_parts_from_nodes(const stk::mesh::EntityVector & deactivatedElements, stk::mesh::Part & activePart)
{
    stk::mesh::EntityVector nodesToDeactivate = get_nodes_to_deactivate(deactivatedElements, activePart);
    OrdinalVector scratchOrdinalVec, scratchSpace;
    for (stk::mesh::Entity nodeToDeactivate : nodesToDeactivate)
    {
        this->internal_change_entity_parts(nodeToDeactivate,{}, {activePart.mesh_meta_data_ordinal()}, scratchOrdinalVec, scratchSpace);
    }
}

unsigned BulkData::num_sides(Entity entity) const
{
    return num_connectivity(entity, mesh_meta_data().side_rank());
}

void BulkData::sort_entities(const stk::mesh::EntitySorterBase& sorter)
{
    ThrowRequireMsg(synchronized_count()>0,"Error, sort_entities must be called after at least one modification cycle.");
    ThrowRequireMsg(in_synchronized_state(), "Error, sort_entities cannot be called from inside a modification cycle.");
    m_bucket_repository.internal_custom_sort_bucket_entities(sorter);

    m_bucket_repository.internal_modification_end();
    internal_update_fast_comm_maps();

    if(parallel_size() > 1)
        check_mesh_consistency();
}

void BulkData::enable_mesh_diagnostic_rule(stk::mesh::MeshDiagnosticFlag flag)
{
    m_meshDiagnosticObserver->enable_rule(flag);
}

unsigned BulkData::get_mesh_diagnostic_error_count() const
{
    return m_meshDiagnosticObserver->get_number_of_errors();
}

void BulkData::throw_on_mesh_diagnostic_error()
{
    m_meshDiagnosticObserver->throw_if_errors_exist();
}

void BulkData::initialize_face_adjacent_element_graph()
{
    if (m_elemElemGraph == nullptr)
    {
        m_elemElemGraph = new ElemElemGraph(*this);
        m_elemElemGraphUpdater = std::make_shared<ElemElemGraphUpdater>(*this,*m_elemElemGraph);
        register_observer(m_elemElemGraphUpdater);
    }
}

void BulkData::delete_face_adjacent_element_graph()
{
    unregister_observer(m_elemElemGraphUpdater);
    delete m_elemElemGraph; m_elemElemGraph = nullptr;
}

stk::mesh::ElemElemGraph& BulkData::get_face_adjacent_element_graph()
{
    ThrowRequireMsg(m_elemElemGraph != nullptr, "Error, Please call initialize_face_adjacent_element_graph before calling get_face_adjacent_element_graph!");
    return *m_elemElemGraph;
}

const stk::mesh::ElemElemGraph& BulkData::get_face_adjacent_element_graph() const
{
    ThrowRequireMsg(m_elemElemGraph != nullptr, "Error, Please call initialize_face_adjacent_element_graph before calling get_face_adjacent_element_graph!");
    return *m_elemElemGraph;
}

bool BulkData::has_face_adjacent_element_graph() const
{
    return m_elemElemGraph != nullptr;
}

void BulkData::destroy_elements_of_topology(stk::topology topologyToDelete)
{
    impl::ElementTopologyDeletions topoBasedBucketDestroyer(*this, topologyToDelete);
    topoBasedBucketDestroyer.find_buckets_to_destroy_and_relations_to_sever();
    break_boundary_relations_and_delete_buckets(topoBasedBucketDestroyer.get_relations_to_sever(), topoBasedBucketDestroyer.get_buckets_to_delete());
}

void BulkData::break_boundary_relations_and_delete_buckets(const std::vector<impl::RelationEntityToNode> & relationsToDestroy, const stk::mesh::BucketVector & bucketsToDelete)
{
    modification_begin();
    for(const impl::RelationEntityToNode & relation : relationsToDestroy)
        destroy_relation(relation.entity, relation.node, relation.ordinal);
    delete_buckets(bucketsToDelete);
    modification_end();
}

void BulkData::delete_buckets(const stk::mesh::BucketVector & buckets)
{
    for(stk::mesh::Bucket * bucket : buckets)
    {
        mark_entities_as_deleted(bucket);
        m_bucket_repository.delete_bucket(bucket);
    }
}

void BulkData::mark_entities_as_deleted(stk::mesh::Bucket * bucket)
{
    for(Entity e : *bucket)
    {
        notifier.notify_entity_deleted(e);
        record_entity_deletion(e);
    }
}

void BulkData::create_side_entities(const SideSet &sideSet, const stk::mesh::PartVector& parts)
{
    if(has_face_adjacent_element_graph())
        FaceCreator(*this, *m_elemElemGraph).create_side_entities_given_sideset(sideSet, parts);
}

void BulkData::dump_mesh_per_proc(const std::string& fileNamePrefix) const
{
    std::ostringstream oss;
    oss << fileNamePrefix << "." << parallel_rank();
    std::ofstream out(oss.str());
    dump_all_mesh_info(out);
    out.close();
}

bool BulkData::does_sideset_exist(const stk::mesh::Part &part) const
{
    return m_sideSetData.does_sideset_exist(part);
}

SideSet& BulkData::create_sideset(const stk::mesh::Part &part, bool fromInput)
{
    return m_sideSetData.create_sideset(part, fromInput);
}

const SideSet& BulkData::get_sideset(const stk::mesh::Part &part) const
{
    return m_sideSetData.get_sideset(part);
}

SideSet& BulkData::get_sideset(const stk::mesh::Part &part)
{
    return m_sideSetData.get_sideset(part);
}

size_t BulkData::get_number_of_sidesets() const
{
    return m_sideSetData.size();
}

bool BulkData::was_mesh_modified_since_sideset_creation()
{
    return m_sideSetData.was_mesh_modified_since_sideset_creation();
}

void BulkData::synchronize_sideset_sync_count()
{
    m_sideSetData.set_sideset_sync_count(this->synchronized_count());
}

void BulkData::clear_sidesets()
{
    m_sideSetData.clear_sidesets();
}

void BulkData::clear_sideset(const stk::mesh::Part &part)
{
    m_sideSetData.clear_sideset(part);
}

std::vector<SideSet *> BulkData::get_sidesets()
{
    return m_sideSetData.get_sidesets();
}

std::vector<const SideSet *> BulkData::get_sidesets() const
{
    return m_sideSetData.get_sidesets();
}

#ifdef SIERRA_MIGRATION
EntityLess::EntityLess(const BulkData& mesh)
  : m_mesh(&mesh),
    m_shouldSortFacesByNodeIds(mesh.should_sort_faces_by_node_ids()),
    m_sideRank(mesh.mesh_meta_data().side_rank())
{}
#endif

} // namespace mesh
} // namespace stk
