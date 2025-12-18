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

#include "stk_mesh/base/Types.hpp"      // for EntityProc, EntityRank, etc
#include <stk_mesh/base/BulkData.hpp>
#include "stk_mesh/base/Entity.hpp"     // for Entity, operator<<, etc
#include "stk_mesh/base/EntityCommDatabase.hpp"  // for pack_entity_info, etc
#include "stk_mesh/baseImpl/EntityCommHelpers.hpp"
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey, etc
#include "stk_mesh/base/EntityLess.hpp"
#include "stk_mesh/base/FieldBase.hpp"  // for FieldBase, FieldMetaData, etc
#include "stk_mesh/base/FieldDataManager.hpp"  // for FieldDataManager, etc
#include "stk_mesh/base/Ghosting.hpp"   // for Ghosting
#include "stk_mesh/base/Part.hpp"       // for Part, remove, etc
#include "stk_mesh/base/Relation.hpp"   // for Relation, etc
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_mesh/base/DestroyRelations.hpp"
#include "stk_mesh/base/ForEachEntity.hpp"
#include "stk_mesh/baseImpl/AuraGhosting.hpp"
#include "stk_mesh/baseImpl/BucketRepository.hpp"  // for BucketRepository
#include "stk_mesh/baseImpl/Visitors.hpp"
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_mesh/baseImpl/ElemDeathImpl.hpp"
#include "stk_mesh/baseImpl/MeshCommImplUtils.hpp"
#include "stk_mesh/baseImpl/MeshCommVerify.hpp"
#include "stk_mesh/baseImpl/PartVectorUtils.hpp"
#include "stk_mesh/baseImpl/MeshModification.hpp"
#include "stk_mesh/baseImpl/CommEntityMods.hpp"
#include <stk_mesh/baseImpl/SideSetPartImpl.hpp>
#include "stk_mesh/baseImpl/ConnectEdgesImpl.hpp"
#include "stk_mesh/baseImpl/Partition.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine, etc
#include "stk_util/util/PairIter.hpp"   // for PairIter
#include "stk_util/util/SortAndUnique.hpp"
#include "stk_util/util/GetEnv.hpp"
#include <algorithm>                    // for sort, lower_bound, unique, etc
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <iterator>                     // for back_insert_iterator, etc
#include <set>                          // for set, set<>::iterator, etc
#include <sstream>
#include <stddef.h>                     // for size_t
#include <stk_mesh/base/Bucket.hpp>     // for Bucket, BucketIdComparator, etc
#include <stk_mesh/base/FaceCreator.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/SideSetEntry.hpp>
#include <stk_mesh/base/SideSetHelper.hpp>
#include <stk_mesh/base/NgpMeshBase.hpp>
#include <stk_mesh/baseImpl/ElementTopologyDeletions.hpp>
#include <stk_mesh/baseImpl/EntityKeyMapping.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraphUpdater.hpp>
#include <stk_mesh/baseImpl/elementGraph/SideConnector.hpp>   // for SideConnector
#include <stk_mesh/baseImpl/elementGraph/SideSharingUsingGraph.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/CommSparse.hpp>  // for CommSparse
#include <stk_util/parallel/GenerateParallelUniqueIDs.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>  // for Reduce, all_reduce, etc
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <stk_util/util/ReportHandler.hpp>  // for ThrowRequireMsg, etc
#include <stk_util/util/StaticAssert.hpp>  // for StaticAssert, etc
#include <stk_util/util/string_case_compare.hpp>
#include <cstring>                      // for strcmp
#include <string>                       // for char_traits, string, etc
#include <utility>                      // for pair, make_pair, swap
#include <vector>                       // for vector, etc

namespace stk {
namespace mesh {

// Static constant on BulkData:
const uint16_t BulkData::orphaned_node_marking = 25000;

///////////////////////////////////////////// Functions for creating entities

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

void removeEntitiesNotSelected(stk::mesh::BulkData &mesh, const stk::mesh::Selector& selected, stk::mesh::EntityVector &entities)
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

        if (only_consider_second_element_from_this_selector != nullptr)
        {
            removeEntitiesNotSelected(*this, *only_consider_second_element_from_this_selector, common_elements);
        }

        if ( common_elements.size() > 0 )
        {
            entities_to_send_data.emplace_back(sentity.global_key, proc_id);
        }
    }
}

void BulkData::unpackEntityFromOtherProcAndUpdateInfoIfSharedLocally(stk::CommSparse &comm, std::vector<shared_entity_type> & shared_entities_this_proc)
{
    std::vector< std::pair<int, shared_entity_type> > shared_entities_received_from_other_procs;
    impl::unpack_shared_entities(*this, comm, shared_entities_received_from_other_procs);

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
    if (entityRank == stk::topology::NODE_RANK) {
        internal_update_sharing_comm_map_and_fill_list_modified_shared_entities_of_node_rank(entities);
    }
    else {
        fill_shared_entities_of_rank_while_updating_sharing_info(entityRank, entities);
    }
    std::sort(entities.begin(), entities.end(), EntityLess(*this));
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
    impl::unpack_shared_entities(*this, comm, shared_entities_and_proc);

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
    impl::pack_entity_keys_to_send(commForInternalSides, entities_to_send_data);
    commForInternalSides.allocate_buffers();
    impl::pack_entity_keys_to_send(commForInternalSides, entities_to_send_data);
    commForInternalSides.communicate();

    std::vector<stk::mesh::EntityKey> receivedEntityKeys;
    impl::unpack_entity_keys_from_procs(commForInternalSides, receivedEntityKeys);
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
BulkData::BulkData(std::shared_ptr<MetaData> mesh_meta_data,
                   ParallelMachine parallel,
                   enum AutomaticAuraOption auto_aura_option,
#ifdef SIERRA_MIGRATION
                   bool add_fmwk_data,
#endif
                   std::unique_ptr<FieldDataManager> field_data_manager,
                   unsigned initialBucketCapacity,
                   unsigned maximumBucketCapacity,
                   std::shared_ptr<impl::AuraGhosting> auraGhosting,
                   bool createUpwardConnectivity)
  :
#ifdef SIERRA_MIGRATION
    m_check_invalid_rels(true),
#endif
    m_createUpwardConnectivity(createUpwardConnectivity),
    m_auraGhosting((auraGhosting!=nullptr ? auraGhosting : std::make_shared<impl::AuraGhosting>())),
    m_entity_comm_map(),
    m_ghosting(),
    m_meta_data(mesh_meta_data),
    m_mark_entity(),
    m_add_node_sharing_called(false),
    m_closure_count(),
    m_mesh_indexes(),
    m_entityKeyMapping(new impl::EntityKeyMapping()),
    m_entitycomm(),
    m_owner(),
    m_comm_list_updater(m_entity_comm_map.comm_list(), m_entitycomm, m_removedGhosts),
    m_entity_keys(),
    m_field_data_manager(std::move(field_data_manager)),
    m_bucket_repository(*this, mesh_meta_data->entity_rank_count(), initialBucketCapacity, maximumBucketCapacity),
#ifdef SIERRA_MIGRATION
    m_add_fmwk_data(add_fmwk_data),
    m_fmwk_global_ids(),
    m_shouldSortFacesByNodeIds(false),
#endif
    m_autoAuraOption(auto_aura_option),
    m_turningOffAutoAura(false),
    m_meshModification(*this),
    m_parallel( parallel ),
    m_ngpMeshHostDataBase(),
    m_all_sharing_procs(mesh_meta_data->entity_rank_count()),
    m_all_sharing_procs_sync_count(0),
    m_ghost_parts(),
    m_num_fields(-1), // meta data not necessarily committed yet
    m_keep_fields_updated(true),
    m_local_ids(),
    m_selector_to_buckets_maps(stk::topology::NUM_RANKS),
    m_use_identifiers_for_resolving_sharing(false),
    m_modSummary(*this),
    m_meshDiagnosticObserver(std::make_shared<stk::mesh::MeshDiagnosticObserver>(*this)),
    m_sideSetData(*this),
    m_ngpMeshBase(nullptr),
    m_isDeviceMeshRegistered(false),
    m_ngpFieldSyncBufferModCount(0),
    m_deviceFastMeshIndicesSynchronizedCount(0),
    m_runConsistencyCheck(false),
    m_symmetricGhostInfo(false),
    m_maintainLocalIds(false),
    m_soloSideIdGenerator(stk::parallel_machine_size(parallel), stk::parallel_machine_rank(parallel), std::numeric_limits<stk::mesh::EntityId>::max()),
    m_supportsLargeIds(false)
{
  mesh_meta_data->set_mesh_bulk_data(this);

  m_entity_comm_map.setCommMapChangeListener(&m_comm_list_updater);

  if (not m_field_data_manager) {
    m_field_data_manager = std::make_unique<FieldDataManager>(mesh_meta_data->entity_rank_count());
  }

  initialize_arrays();

  m_ghost_parts.clear();
  internal_create_ghosting( "shared" );
  //shared part should reside in m_ghost_parts[0]
  internal_create_ghosting( "shared_aura" );

  impl::check_size_of_types();

  register_observer(m_meshDiagnosticObserver);

  init_mesh_consistency_check_mode();

  m_meshModification.set_sync_state_synchronized();
}

BulkData::~BulkData()
{
  while ( ! m_ghosting.empty() ) {
    delete m_ghosting.back();
    m_ghosting.pop_back();
  }

  if (this == &mesh_meta_data().mesh_bulk_data()) {
    mesh_meta_data().set_mesh_bulk_data(nullptr);
  }

  delete m_elemElemGraph;

  delete m_ngpMeshBase;
}

void
BulkData::register_device_mesh() const
{
  STK_ThrowRequireMsg(m_isDeviceMeshRegistered==false, "Cannot register more than one DeviceMesh against a BulkData");
  m_isDeviceMeshRegistered = true;
}

void BulkData::set_automatic_aura_option(AutomaticAuraOption auraOption, bool applyImmediately)
{
  STK_ThrowRequireMsg(in_synchronized_state(),"set_automatic_aura_option currently can only be used when the mesh is not already being modified.");

  if (auraOption != m_autoAuraOption) {
    m_autoAuraOption = auraOption;

    if (applyImmediately) {
      modification_begin();
      if (m_autoAuraOption == BulkData::NO_AUTO_AURA) {
        m_auraGhosting->remove_aura(*this);
      }
      modification_end();
    }
    else if (m_autoAuraOption == BulkData::NO_AUTO_AURA) {
      m_turningOffAutoAura = true;
    }
  }
}

size_t BulkData::ngp_mesh_synchronized_count() const
{
  return get_ngp_mesh()->synchronized_count();
}


//----------------------------------------------------------------------
//----------------------------------------------------------------------

void BulkData::require_ok_to_modify() const
{
  STK_ThrowRequireMsg( !this->in_synchronized_state(), "NOT in the ok-to-modify state" );
}

void BulkData::require_entity_owner( const Entity entity ,
                                     int owner ) const
{
  if (parallel_size() > 1 && bucket_ptr(entity) != nullptr) {
    const bool error_not_owner = owner != parallel_owner_rank(entity) ;

    STK_ThrowRequireMsg( !error_not_owner,
                     "P" << parallel_rank() << " " << entity_key(entity) << " owner is " <<
                     parallel_owner_rank(entity) << ", expected " << owner);
  }
}

void BulkData::require_good_rank_and_id(EntityRank ent_rank, EntityId ent_id) const
{
  const EntityRank rank_count = mesh_meta_data().entity_rank_count();
  const bool ok_id   = EntityKey::is_valid_id(ent_id);
  const bool ok_rank = ent_rank < rank_count && !(ent_rank == stk::topology::FACE_RANK && mesh_meta_data().spatial_dimension() == 2);

  STK_ThrowRequireMsg( ok_rank, "Bad key rank: " << ent_rank << " for id " << ent_id );

  STK_ThrowRequireMsg( ok_id, "Bad id : " << ent_id);
}

void BulkData::mark_entity_and_upward_related_entities_as_modified(Entity entity, bool markUpDownClosureIfShared)
{
  BulkData& mesh = *this;

  auto markAsModified = [&](Entity ent) { mesh.set_state(ent, Modified); };

  auto onlyVisitUnchanged = [&](Entity ent) { return mesh.state(ent) == Unchanged; };

  if (markUpDownClosureIfShared && in_shared(entity)) {
    impl::VisitUpDownClosureGeneral(*this, entity, markAsModified, onlyVisitUnchanged);
  }
  else {
    const EntityRank entityRank = entity_rank(entity);
    const EntityRank endRank = static_cast<EntityRank>(mesh_meta_data().entity_rank_count());
    if (mesh.state(entity) == Unchanged) {
      impl::VisitUpwardClosureGeneral(mesh, entity, entityRank, endRank, markAsModified, onlyVisitUnchanged);
    }
    else if (mesh.state(entity) == Modified) {
      const EntityRank beginRank = static_cast<EntityRank>(entityRank+1);

      for(EntityRank rank=beginRank; rank<endRank; ++rank) {
        const ConnectedEntities connected = get_connected_entities(entity, rank);
        if (connected.size() > 0) {
          impl::VisitUpwardClosureGeneral(mesh, connected.data(), connected.data()+connected.size(), rank, endRank, markAsModified, onlyVisitUnchanged);
        }
      }
    }
  }
}

size_t BulkData::count_relations(Entity entity, bool onlyDownwardRelations) const
{
  const MeshIndex &mesh_idx = mesh_index(entity);

  const EntityRank end_rank = onlyDownwardRelations ? entity_rank(entity) : static_cast<EntityRank>(mesh_meta_data().entity_rank_count());
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

  const EntityRank end_rank = static_cast<EntityRank>(mesh_meta_data().entity_rank_count());
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
  else {
    Entity::entity_value_type local_offset = m_meshModification.get_deleted_entity_cache().get_entity_for_reuse();
    if (local_offset != Entity::InvalidEntity)
    {
      new_local_offset = local_offset;
    }
  }

  MeshIndex mesh_index = {nullptr, 0};
  EntityKey invalid_key;

  if (new_local_offset == m_mesh_indexes.size()) {
    m_mesh_indexes.push_back(mesh_index);
    m_entity_keys.push_back(invalid_key);
    m_entitycomm.push_back(-1);
    m_owner.push_back(parallel_rank());
    m_meshModification.add_created_entity_state();
    m_mark_entity.push_back(NOT_MARKED);
    m_closure_count.push_back(static_cast<uint16_t>(0));
    m_local_ids.push_back(stk::mesh::GetInvalidLocalId());

#ifdef SIERRA_MIGRATION
    if (add_fmwk_data()) {
      m_fmwk_global_ids.push_back(0);
    }
#endif
  }
  else {
    //re-claiming space from a previously-deleted entity:

    m_mesh_indexes[new_local_offset] = mesh_index;
    m_entity_keys[new_local_offset] = invalid_key;
    m_entitycomm[new_local_offset] = -1;
    m_owner[new_local_offset] = parallel_rank();
    m_mark_entity[new_local_offset] = NOT_MARKED;
    m_meshModification.mark_entity_as_created(new_local_offset);
    m_closure_count[new_local_offset] = static_cast<uint16_t>(0);
    m_local_ids[new_local_offset] = stk::mesh::GetInvalidLocalId();

#ifdef SIERRA_MIGRATION
    if (add_fmwk_data()) {
      m_fmwk_global_ids[new_local_offset] = 0;
    }
#endif
  }

  return Entity(new_local_offset);
}

void BulkData::initialize_arrays()
{
  STK_ThrowRequireMsg((m_mesh_indexes.size() == 0) && (m_entity_keys.size() == 0),
                   "BulkData::initialize_arrays() called by something other than constructor");

  MeshIndex mesh_index = {nullptr, 0};
  m_mesh_indexes.push_back(mesh_index);

  EntityKey invalid_key;
  m_entity_keys.push_back(invalid_key);
  m_entitycomm.push_back(-1);
  m_owner.push_back(parallel_rank());

  m_mark_entity.push_back(NOT_MARKED);
  m_closure_count.push_back(static_cast<uint16_t>(0));
  m_local_ids.push_back(stk::mesh::GetInvalidLocalId());

#ifdef SIERRA_MIGRATION
  if (add_fmwk_data()) {
    m_fmwk_global_ids.push_back(0);
  }
#endif
}

#ifdef SIERRA_MIGRATION
void
BulkData::initialize_global_ids()
{
  m_fmwk_global_ids.clear();
  m_fmwk_global_ids.resize(m_entity_keys.size(), 0);
  for (stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank < m_meta_data->entity_rank_count(); ++rank) {
    for (const stk::mesh::Bucket * bucket : buckets(rank)) {
      for (stk::mesh::Entity entity : *bucket) {
        m_fmwk_global_ids[entity.local_offset()] = identifier(entity);
      }
    }
  }
}
#endif

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
    STK_ThrowRequireMsg(mesh_meta_data().side_rank() != stk::topology::EDGE_RANK, "Program Error!");
    return declare_entity(stk::topology::EDGE_RANK, id);
}

Entity BulkData::declare_entity( EntityRank ent_rank , EntityId ent_id)
{
    STK_ThrowRequireMsg(ent_rank != mesh_meta_data().side_rank(),"declare_entity for side is not supported. Please use declare_element_side() functionality.");
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
    STK_ThrowRequireMsg(mesh_meta_data().side_rank() != stk::topology::EDGE_RANK, "Program Error!");
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
  stk::topology elemTop = bucket(elem).topology();
  stk::topology sideTop = elemTop.side_topology(localSideId);
  Entity side = internal_declare_entity(sideTop.rank(), globalSideId, add_root_topology_part(parts, mesh_meta_data().get_topology_root_part(sideTop)));

  if (has_face_adjacent_element_graph()) {
    use_graph_to_connect_side(get_face_adjacent_element_graph(), side, elem, localSideId);
  }
  else {
    impl::connect_element_to_side(*this, elem, side, localSideId, parts, sideTop);

    impl::connect_side_to_other_elements(*this, side, elem, localSideId);
  }
  if(this->num_edges(elem) != 0) {
    impl::connect_face_to_edges(*this, side);
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
    if(is_valid(side)) {
        if(bucket(side).owned()) {
            change_entity_parts(side, parts, PARTVECTOR{});
        }
    }
    else {
        stk::topology elemTop = bucket(elem).topology();
        stk::topology sideTop = elemTop.side_topology(sideOrd);
        EntityKey sideKey(sideTop.rank(), globalSideId);

        std::pair<Entity,bool> result = internal_get_or_create_entity_with_notification(sideKey);
        side = result.first;
        const bool newlyCreated = result.second;

        if (newlyCreated) {
          PARTVECTOR allParts = add_root_topology_part(parts, mesh_meta_data().get_topology_root_part(sideTop));
          allParts.push_back(&mesh_meta_data().locally_owned_part());

          OrdinalVector scratch1, scratch2;
          internal_verify_and_change_entity_parts(side, allParts, PARTVECTOR{}, scratch1, scratch2);

          internal_set_owner(side, parallel_rank());
        }

        if (has_face_adjacent_element_graph()) {
          use_graph_to_connect_side(get_face_adjacent_element_graph(), side, elem, sideOrd);
        }
        else {
          impl::connect_element_to_side(*this, elem, side, sideOrd, parts, sideTop);

          impl::connect_side_to_other_elements(*this, side, elem, sideOrd);
        }

        if(this->num_edges(elem) != 0) {
          impl::connect_face_to_edges(*this, side);
        }
    }
    return side;
}

template Entity BulkData::declare_element_side_with_id<stk::mesh::PartVector>(const stk::mesh::EntityId, Entity, const unsigned, const stk::mesh::PartVector&);
template Entity BulkData::declare_element_side_with_id<stk::mesh::ConstPartVector>(const stk::mesh::EntityId, Entity, const unsigned, const stk::mesh::ConstPartVector&);


void BulkData::entity_comm_list_insert(Entity node)
{
  EntityCommListInfoVector& commList = m_entity_comm_map.comm_list();
  stk::mesh::EntityKey key = entity_key(node);
  EntityCommListInfoVector::iterator lb_itr = std::lower_bound(commList.begin(), commList.end(), key);
  if(lb_itr == commList.end() || lb_itr->key != key)
  {
    const int entityCommIndex = m_entity_comm_map.entity_comm(key);
    EntityCommListInfo comm_info = {key, node, entityCommIndex};
    commList.insert(lb_itr, comm_info);
  }
}

void BulkData::add_node_sharing(Entity node, int sharing_proc)
{
  STK_ThrowRequireMsg(sharing_proc != parallel_rank() && sharing_proc < parallel_size(),
      "add_node_sharing: illegal sharing_proc="<<sharing_proc
       <<", must be < numProcs("<<parallel_size()<<") and != localProc("<<parallel_rank()<<")");
  // Only valid to specify sharing information for non-deleted nodes
  STK_ThrowRequire(entity_rank(node) == stk::topology::NODE_RANK);
  STK_ThrowRequire(state(node) != Deleted);

  protect_orphaned_node(node);

  m_add_node_sharing_called = true;

  if (state(node) == Unchanged) {
      mark_entity_and_upward_related_entities_as_modified(node);
  }

  internal_mark_entity(node, IS_SHARED);
  entity_comm_map_insert(node, EntityCommInfo(stk::mesh::BulkData::SHARED, sharing_proc));
  entity_comm_map_erase(entity_key(node), EntityCommInfo(stk::mesh::BulkData::AURA, sharing_proc));
  entity_comm_list_insert(node);
}

void BulkData::add_node_sharing(const EntityProcVec& nodesAndProcs)
{
  EntityCommListInfoVector newCommListEntries;
  newCommListEntries.reserve(nodesAndProcs.size());

  for(const EntityProc& nodeAndSharingProc : nodesAndProcs) {
    Entity node = nodeAndSharingProc.first;
    int sharing_proc = nodeAndSharingProc.second;

    // Only valid to specify sharing information for non-deleted nodes
    STK_ThrowRequire(entity_rank(node) == stk::topology::NODE_RANK);
    STK_ThrowRequire(state(node) != Deleted);

    protect_orphaned_node(node);

    if (state(node) == Unchanged) {
        mark_entity_and_upward_related_entities_as_modified(node);
    }

    internal_mark_entity(node, IS_SHARED);
    std::pair<int,bool> result = entity_comm_map_insert(node, EntityCommInfo(stk::mesh::BulkData::SHARED, sharing_proc));
    entity_comm_map_erase(entity_key(node), EntityCommInfo(stk::mesh::BulkData::AURA, sharing_proc));

    const bool inserted = result.second;
    if (inserted) {
      EntityKey key = entity_key(node);
      EntityCommListInfo info = {key, node, result.first};
      newCommListEntries.push_back(info);
    }
  }

  internal_add_comm_list_entries(newCommListEntries);
  m_add_node_sharing_called = true;
}

EntityId BulkData::get_solo_side_id()
{
    return m_soloSideIdGenerator.get_solo_side_id();
}

template<typename PARTVECTOR>
void BulkData::internal_verify_and_change_entity_parts( Entity entity,
                                                        const PARTVECTOR & add_parts ,
                                                        const PARTVECTOR & remove_parts,
                                                        OrdinalVector & scratchOrdinalVec,
                                                        OrdinalVector & scratchSpace)
{
    require_ok_to_modify();

    OrdinalVector addPartsAndSupersets;
    impl::fill_add_parts_and_supersets(add_parts, addPartsAndSupersets);

    OrdinalVector removePartsAndSubsetsMinusPartsInAddPartsList;
    impl::fill_remove_parts_and_subsets_minus_parts_in_add_parts_list(remove_parts,
                                                      addPartsAndSupersets,
                                                      bucket_ptr(entity),
                                                      removePartsAndSubsetsMinusPartsInAddPartsList);

    scratchOrdinalVec.clear();
    scratchSpace.clear();
    internal_change_entity_parts(entity,
                                 addPartsAndSupersets,
                                 removePartsAndSubsetsMinusPartsInAddPartsList,
                                 scratchOrdinalVec, scratchSpace);
}

template void BulkData::internal_verify_and_change_entity_parts(Entity, const PartVector&, const PartVector&, OrdinalVector&, OrdinalVector&);
template void BulkData::internal_verify_and_change_entity_parts(Entity, const ConstPartVector&, const ConstPartVector&, OrdinalVector&, OrdinalVector&);

template<typename PARTVECTOR>
void BulkData::internal_verify_and_change_entity_parts( const EntityVector& entities,
                                                        const PARTVECTOR & add_parts ,
                                                        const PARTVECTOR & remove_parts)
{
    require_ok_to_modify();

    OrdinalVector addPartsAndSupersets;
    OrdinalVector removePartsAndSubsetsMinusPartsInAddPartsList;
    OrdinalVector scratchOrdinalVec, scratchSpace;

    impl::fill_add_parts_and_supersets(add_parts, addPartsAndSupersets);

    for(Entity entity : entities) {
      impl::fill_remove_parts_and_subsets_minus_parts_in_add_parts_list(remove_parts,
                                                        addPartsAndSupersets,
                                                        bucket_ptr(entity),
                                                        removePartsAndSubsetsMinusPartsInAddPartsList);

      internal_change_entity_parts(entity, addPartsAndSupersets,
                                   removePartsAndSubsetsMinusPartsInAddPartsList,
                                   scratchOrdinalVec, scratchSpace);
    }
}

template void BulkData::internal_verify_and_change_entity_parts(const EntityVector&, const PartVector&, const PartVector&);
template void BulkData::internal_verify_and_change_entity_parts(const EntityVector&, const ConstPartVector&, const ConstPartVector&);

template<typename PARTVECTOR>
void BulkData::internal_verify_and_change_entity_parts( const Selector& selector,
                                                        const EntityRank rank,
                                                        const PARTVECTOR & add_parts ,
                                                        const PARTVECTOR & remove_parts)
{
  require_ok_to_modify();

  OrdinalVector addPartsAndSupersets;
  OrdinalVector removePartsAndSubsetsMinusPartsInAddPartsList;
  OrdinalVector scratchOrdinalVec, scratchSpace;

  impl::fill_add_parts_and_supersets(add_parts, addPartsAndSupersets);

  impl::fill_remove_parts_and_subsets_minus_parts_in_add_parts_list(remove_parts,
                                                    addPartsAndSupersets,
                                                    get_buckets(rank, selector),
                                                    removePartsAndSubsetsMinusPartsInAddPartsList);

  internal_change_entity_parts(selector,
                               rank,
                               addPartsAndSupersets,
                               removePartsAndSubsetsMinusPartsInAddPartsList,
                               scratchOrdinalVec, scratchSpace);
}

template void BulkData::internal_verify_and_change_entity_parts(const Selector&, const EntityRank, const PartVector&, const PartVector&);
template void BulkData::internal_verify_and_change_entity_parts(const Selector&, const EntityRank, const ConstPartVector&, const ConstPartVector&);


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

  Part * const owns = & mesh_meta_data().locally_owned_part();

  PARTVECTOR rem ;
  PARTVECTOR add( parts );
  add.push_back( owns );

  OrdinalVector scratchOrdinalVec, scratchSpace;
  internal_verify_and_change_entity_parts( declared_entity , add , rem, scratchOrdinalVec, scratchSpace );

  if ( result.second ) {
    this->internal_set_owner(declared_entity, parallel_rank());
  }

  m_check_invalid_rels = true;

  return declared_entity ;
}

template Entity BulkData::internal_declare_entity(EntityRank ent_rank, EntityId ent_id, const PartVector& parts);
template Entity BulkData::internal_declare_entity(EntityRank ent_rank, EntityId ent_id, const ConstPartVector& parts);

bool entity_is_purely_local(const BulkData& mesh, Entity entity)
{
    const Bucket& bucket = mesh.bucket(entity);
    return bucket.owned() && !bucket.shared()
            && !bucket.in_aura() && !mesh.in_send_ghost(entity);
}

void require_fmwk_or_entity_purely_local(const BulkData& mesh, Entity entity, const std::string& caller)
{
    STK_ThrowRequireMsg(mesh.add_fmwk_data() || entity_is_purely_local(mesh, entity),
                   "Error, "<<caller<< " requires that stk-mesh is running under fmwk, or entity "<<mesh.entity_key(entity)<<" must be purely local.");
}

void BulkData::change_entity_id( EntityId id, Entity entity)
{
// THIS STK_ThrowAssertMsg IS ONLY MACRO CONTROLLED TO ALLOW EXPERIMENTATION WITH
// Fmwk USING stk_parallel.  WHEN stk parallel IS USED WITHIN Fmwk, THIS ASSERTION
// IS VIOLATED.
#ifndef SIERRA_MIGRATION
  STK_ThrowAssertMsg(parallel_size() == 1,
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
    m_entityKeyMapping->update_entity_key(new_key, old_key, entity);
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
  if(impl::can_destroy_entity(*this, entity)) {
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

  if(impl::has_upward_connectivity(*this, entity)) {
    m_check_invalid_rels = true;
    return false;
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

  const stk::mesh::EntityRank end_rank = static_cast<stk::mesh::EntityRank>(mesh_meta_data().entity_rank_count());
  for (stk::mesh::EntityRank irank = end_rank; irank != stk::topology::BEGIN_RANK; ) {
    --irank;
    destroy_relations(*this, entity, irank);
  }

  // Need to invalidate Entity handles in comm-list
  const stk::mesh::EntityKey key = entity_key(entity);
  EntityCommListInfoVector& commList = m_entity_comm_map.comm_list();
  EntityCommListInfoVector::iterator lb_itr =
    std::lower_bound(commList.begin(), commList.end(), key);
  if (lb_itr != commList.end() && lb_itr->key == key) {
    lb_itr->entity = stk::mesh::Entity();
  }

  remove_entity_callback(erank, bucket(entity).bucket_id(), bucket_ordinal(entity));

  m_bucket_repository.remove_entity(mesh_index(entity));

  record_entity_deletion(entity, ghost);

  m_check_invalid_rels = true;
  return true ;
}

void BulkData::record_entity_deletion(Entity entity, bool isGhost)
{
    const EntityKey key = entity_key(entity);
    set_mesh_index(entity, 0, 0);
    m_entityKeyMapping->destroy_entity(key, entity);
    notifier.notify_local_entities_created_or_deleted(key.rank());
    notifier.notify_local_buckets_changed(key.rank());
    m_meshModification.mark_entity_as_deleted(entity, isGhost);
    m_mark_entity[entity.local_offset()] = NOT_MARKED;
    m_closure_count[entity.local_offset()] = static_cast<uint16_t>(0u);
}

size_t get_max_num_ids_needed_across_all_procs(const stk::mesh::BulkData& bulkData, size_t numIdsNeededThisProc)
{
    size_t maxNumNeeded = numIdsNeededThisProc;
    STK_ThrowRequireMsg(MPI_Allreduce(MPI_IN_PLACE, &maxNumNeeded, 1, sierra::MPI::Datatype<size_t>::type(), MPI_MAX, bulkData.parallel()) == MPI_SUCCESS,
            "Program error (MPI_Allreduce failure). Please contact sierra-help@sandia.gov for support.");
    return maxNumNeeded;
}

//----------------------------------------------------------------------

std::vector<uint64_t> BulkData::internal_get_ids_in_use(stk::topology::rank_t rank, const std::vector<stk::mesh::EntityId>& reserved_ids) const
{
  const_entity_iterator beg = begin_entities(rank);
  const_entity_iterator end = end_entities(rank);
  std::vector<uint64_t> ids_in_use;
  ids_in_use.reserve(std::distance(beg,end) + m_meshModification.get_deleted_entity_cache().get_deleted_entities_current_mod_cycle().size() + reserved_ids.size());

  for(const_entity_iterator i = beg; i!=end; ++i) {
    ids_in_use.push_back(i->first.id());
  }

  std::vector<uint64_t> reserved_and_deleted_ids(reserved_ids.begin(), reserved_ids.end());
  for (Entity::entity_value_type local_offset : m_meshModification.get_deleted_entity_cache().get_deleted_entities_current_mod_cycle()) {
    stk::mesh::Entity entity;
    entity.set_local_offset(local_offset);
    if ((entity_rank(entity) == rank) && (is_valid(entity) || state(entity)==Deleted)) {
      reserved_and_deleted_ids.push_back(identifier(entity));
    }
  }

  if (!reserved_and_deleted_ids.empty()) {
    stk::util::sort_and_unique(reserved_and_deleted_ids);
    stk::util::insert_keep_sorted_and_unique(reserved_and_deleted_ids, ids_in_use);
  }
  return ids_in_use;
}

uint64_t  BulkData::get_max_allowed_id() const {
  if(add_fmwk_data()) {
#ifdef SIERRA_MIGRATION
    return std::numeric_limits<FmwkId>::max();
#else
    return stk::mesh::EntityKey::MAX_ID;
#endif
  }
  return stk::mesh::EntityKey::MAX_ID;
}

void BulkData::generate_new_ids_given_reserved_ids(stk::topology::rank_t rank, size_t numIdsNeeded, const std::vector<stk::mesh::EntityId>& reserved_ids, std::vector<stk::mesh::EntityId>& requestedIds) const
{
    size_t maxNumNeeded = get_max_num_ids_needed_across_all_procs(*this, numIdsNeeded);
    if ( maxNumNeeded == 0 ) return;
    EntityId globalMaxId = impl::get_global_max_id_in_use(*this, rank, m_meshModification.get_deleted_entity_cache().get_deleted_entities_current_mod_cycle(), reserved_ids);
    uint64_t maxAllowedId = get_max_allowed_id();
    uint64_t availableIds = maxAllowedId - globalMaxId;
    uint64_t globalNumIdsRequested = stk::get_global_sum(this->parallel(), numIdsNeeded);
    if (availableIds > globalNumIdsRequested) {
      generate_parallel_ids_above_existing_max(this->parallel(), numIdsNeeded, globalNumIdsRequested, maxNumNeeded, availableIds, globalMaxId, requestedIds);
      return;
    }

    std::vector<uint64_t> ids_in_use = this->internal_get_ids_in_use(rank, reserved_ids);
    generate_parallel_ids_in_gap(this->parallel(), ids_in_use, maxAllowedId, numIdsNeeded, globalNumIdsRequested, requestedIds);
}

void BulkData::generate_new_ids(stk::topology::rank_t rank, size_t numIdsNeeded, std::vector<stk::mesh::EntityId>& requestedIds)
{
  if (rank == mesh_meta_data().side_rank()) {
    requestedIds.resize(numIdsNeeded);
    for(size_t i = 0; i < numIdsNeeded; i++) {
      requestedIds[i] = m_soloSideIdGenerator.get_solo_side_id();
    }
    return;
  }

  uint64_t globalNumIdsRequested = 0;
  uint64_t maxIdsRequested = 0;
  stk::compute_global_sum_and_max(parallel(), numIdsNeeded,
                                       globalNumIdsRequested, maxIdsRequested);
  if ( globalNumIdsRequested == 0 ) return;

  EntityId globalMaxId = impl::get_global_max_id_in_use(*this, rank,
                                   m_meshModification.get_deleted_entity_cache().get_deleted_entities_current_mod_cycle());

  uint64_t maxAllowedId = get_max_allowed_id();
  uint64_t availableIds = maxAllowedId - globalMaxId;
  if (globalNumIdsRequested < availableIds) {
    stk::generate_parallel_ids_above_existing_max(parallel(), numIdsNeeded,
                             globalNumIdsRequested, maxIdsRequested, availableIds,
                             globalMaxId, requestedIds);
    return;
  }

  std::vector<uint64_t> ids_in_use = internal_get_ids_in_use(rank);
  uint64_t localMaxIdInUse = ids_in_use.back();
  requestedIds = generate_parallel_unique_ids(maxAllowedId, ids_in_use, numIdsNeeded, parallel(), localMaxIdInUse);
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
    Part * const owns = &mesh_meta_data().locally_owned_part();

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

    std::pair<entity_iterator ,bool> entityBoolPair = m_entityKeyMapping->internal_create_entity(key);

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
    require_ok_to_modify();

    requested_entities.resize(newIds.size());
    if (newIds.empty()) {
        return;
    }

    OrdinalVector partsAndSupersets;
    impl::fill_add_parts_and_supersets(parts, partsAndSupersets);
    Part * const ownedPart = & mesh_meta_data().locally_owned_part();
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
    OrdinalVector inducible_parts_added, inducible_parts_removed;
    OrdinalVector scratchOrdinalVec, scratchSpace;

    internal_fill_new_part_list_and_removed_part_list(bucket_ptr(requested_entities[0]), partsAndSupersets, rem,
                                                      newBucketPartList, parts_removed);
    internal_determine_inducible_parts(entity_rank(requested_entities[0]), partsAndSupersets, parts_removed,
                                       inducible_parts_added, inducible_parts_removed);

    notifier.notify_local_entities_created_or_deleted(rank);
    notifier.notify_local_buckets_changed(rank);

    for(size_t i=0;i<sortedIds.size();++i)
    {
        EntityKey key(rank, sortedIds[i].first);
        require_good_rank_and_id(key.rank(), key.id());

        m_modSummary.track_declare_entity(key.rank(), key.id(), stk::mesh::PartVector());

        std::pair<entity_iterator ,bool> entityBoolPair = m_entityKeyMapping->internal_create_entity(key);

        STK_ThrowErrorMsgIf( ! entityBoolPair.second,
                "Generated id " << key.id() << " of rank " << key.rank() << " which was already used.");

        entityBoolPair.first->second = requested_entities[sortedIds[i].second];
        Entity new_entity = entityBoolPair.first->second;
        this->set_entity_key(new_entity, key);
        notifier.notify_entity_added(new_entity);

        notifier.notify_entity_parts_added(new_entity, partsAndSupersets);

        internal_adjust_closure_count(new_entity, partsAndSupersets, rem);
        internal_move_entity_to_new_bucket(new_entity, newBucketPartList, scratchSpace);
        internal_propagate_induced_part_changes_to_downward_connected_entities(new_entity, inducible_parts_added, inducible_parts_removed, scratchOrdinalVec, scratchSpace);

        this->internal_set_owner(new_entity, parallel_rank());
    }
    //now clear the created-entity cache...
    begin_entities(rank);
}

template
void BulkData::declare_entities(stk::topology::rank_t rank, const std::vector<int>& newIds, const PartVector& parts, std::vector<Entity>& requested_entities);
template
void BulkData::declare_entities(stk::topology::rank_t rank, const std::vector<int64_t>& newIds, const PartVector& parts, std::vector<Entity>& requested_entities);
template
void BulkData::declare_entities(stk::topology::rank_t rank, const std::vector<unsigned long>& newIds, const PartVector& parts, std::vector<Entity>& requested_entities);
template
void BulkData::declare_entities(stk::topology::rank_t rank, const std::vector<unsigned long long>& newIds, const PartVector& parts, std::vector<Entity>& requested_entities);

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
  const int entityCommIndex = m_entitycomm[entity.local_offset()];
  if (entityCommIndex != -1) {
    PairIterEntityComm commInfo = internal_comm_db().comm(entityCommIndex);

    while(!commInfo.empty() && commInfo->ghost_id == SHARED) {
      if (commInfo->proc == proc) {
        return true;
      }
      ++commInfo;
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
  const EntityCommInfo* i = std::lower_bound( ec.begin(), ec.end() , tmp );

  return i != ec.end() && tmp == *i ;
}

bool BulkData::in_ghost( const Ghosting & ghost , Entity entity , int proc ) const
{
  if (m_entitycomm[entity.local_offset()] == -1) {
    return false;
  }

  PairIterEntityComm commInfo = internal_comm_db().comm(m_entitycomm[entity.local_offset()]);

  while(!commInfo.empty() && commInfo->ghost_id < ghost.ordinal()) {
    ++commInfo;
  }

  while(!commInfo.empty() && commInfo->ghost_id == ghost.ordinal()) {
    if (commInfo->proc == proc) {
      return true;
    }
    ++commInfo;
  }

  return false;
}

bool BulkData::in_ghost( const Ghosting & ghost , Entity entity ) const
{
  if (m_entitycomm[entity.local_offset()] == -1) {
    return false;
  }

  PairIterEntityComm commInfo = internal_comm_db().comm(m_entitycomm[entity.local_offset()]);

  const EntityCommInfo* i = commInfo.begin();
  const EntityCommInfo* end = commInfo.end();

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
  return in_send_ghost(ghost, get_entity(key), proc);
}

bool BulkData::in_send_ghost( const Ghosting & ghost , Entity entity , int proc ) const
{
  bool ret_val = false;
  const int owner_rank = parallel_owner_rank(entity);

  if (owner_rank == parallel_rank()) {
    const int entityCommIndex = m_entitycomm[entity.local_offset()];
    if (entityCommIndex >= 0) {
      EntityCommInfo tmp( ghost.ordinal() , proc );

      PairIterEntityComm ec = internal_comm_db().comm(m_entitycomm[entity.local_offset()]);
      const EntityCommInfo* i = std::lower_bound( ec.begin(), ec.end() , tmp );

      ret_val = i != ec.end() && tmp == *i ;
    }
  }

  return ret_val;
}

bool BulkData::is_communicated_with_proc(Entity entity, int proc) const
{
  if (m_entitycomm[entity.local_offset()] == -1) {
    return false;
  }

  PairIterEntityComm commInfo = internal_comm_db().comm(m_entitycomm[entity.local_offset()]);

  const EntityCommInfo* i = commInfo.begin();
  const EntityCommInfo* end = commInfo.end();

  while(i != end) {
    if (i->proc == proc) {
      return true;
    }
    ++i;
  }

  return false;
}

void BulkData::comm_procs(Entity entity, std::vector<int> & procs ) const
{
  procs.clear();
  const int entityCommIndex = m_entitycomm[entity.local_offset()];
  if (entityCommIndex != -1) {
    impl::fill_sorted_procs(internal_comm_db().comm(entityCommIndex), procs);
  }
}

void BulkData::comm_shared_procs(EntityKey key, std::vector<int> & procs ) const
{
  procs.clear();
  const int entityCommIndex = m_entity_comm_map.entity_comm(key);
  if (entityCommIndex != -1) {
    impl::comm_shared_procs(internal_comm_db().comm(entityCommIndex), procs);
  }
}

void BulkData::comm_shared_procs(Entity entity, std::vector<int> & procs ) const
{
  procs.clear();
  const int entityCommIndex = m_entitycomm[entity.local_offset()];
  if (entityCommIndex != -1) {
    impl::comm_shared_procs(internal_comm_db().comm(entityCommIndex), procs);
  }
}

void BulkData::shared_procs_intersection(const std::vector<EntityKey> & keys, std::vector<int> & procs ) const
{
  confirm_host_mesh_is_synchronized_from_device();
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
  confirm_host_mesh_is_synchronized_from_device();
  procs.clear();
  int num = entities.size();
  for (int i = 0; i < num; ++i) {

    if (i == 0) {
      comm_shared_procs(entities[i], procs);
    }
    else {
      PairIterEntityComm sharing = internal_entity_comm_map_shared(entities[i]);

      auto notFoundInSharing = [&sharing](const int& proc) {
        auto compare = [&proc](const EntityCommInfo& info){ return info.proc == proc; };
        return std::find_if(sharing.begin(), sharing.end(), compare) == sharing.end();
      };

      procs.erase(std::remove_if(procs.begin(), procs.end(), notFoundInSharing),
                  procs.end());
    }

    if (procs.empty()) {
      return;
    }
  }
}

void BulkData::comm_procs( const Ghosting & ghost ,
                           EntityKey key, std::vector<int> & procs ) const
{
  impl::fill_ghosting_procs(internal_entity_comm_map(key), ghost.ordinal(), procs);
}

void BulkData::internal_set_owner(Entity entity, int new_owner)
{
  m_owner[entity.local_offset()] = new_owner;
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

  const unsigned totalNumFields = mesh_meta_data().get_fields().size();
  if (m_num_fields == -1) {  // hasn't been set yet
    m_num_fields = totalNumFields;
  }

  for(EntityRank rank = stk::topology::NODE_RANK; rank < mesh_meta_data().entity_rank_count(); ++rank) {
      const std::vector<Bucket*>& buckets = this->buckets(rank);
      const FieldVector fieldsOfRank = mesh_meta_data().get_fields(rank);
      m_field_data_manager->allocate_field_data(rank, buckets, fieldsOfRank, totalNumFields);
  }
}

void BulkData::reallocate_field_data(stk::mesh::FieldBase& newField)
{
  const FieldVector& allFields = mesh_meta_data().get_fields();
  for (FieldBase* stkField : allFields) {
    stkField->synchronize<stk::mesh::ReadOnly, stk::ngp::HostSpace>();
  }

  const EntityRank rank = newField.entity_rank();
  const FieldVector& fieldsOfRank = mesh_meta_data().get_fields(rank);
  const BucketVector& bucketsOfRank = this->buckets(rank);
  const unsigned totalNumFields = mesh_meta_data().get_fields().size();

  m_field_data_manager->reallocate_field_data(rank, bucketsOfRank, newField, fieldsOfRank, totalNumFields);
  for (Bucket * bucket : bucketsOfRank) {
    bucket->mark_for_modification();
  }
  m_meshModification.increment_sync_count();
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
    if (new_bucket != nullptr) {
      SelectorBucketMap& selectorBucketMap = m_selector_to_buckets_maps[rank];
      for (SelectorBucketMap::iterator itr = selectorBucketMap.begin(), end = selectorBucketMap.end();
           itr != end; ++itr) {
        Selector const& sel = itr->first;
        if (sel(*new_bucket)) {
          BucketVector & cached_buckets = itr->second;
          BucketVector::iterator lb_itr = std::lower_bound(cached_buckets.begin(), cached_buckets.end(), new_bucket, BucketIdComparator());
          cached_buckets.insert(lb_itr, new_bucket);
        }
      }
    }
}

void BulkData::remove_bucket_caching(EntityRank rank, const Selector selector)
{
  SelectorBucketMap& selectorBucketMap = m_selector_to_buckets_maps[rank];
  for(SelectorBucketMap::iterator iter = selectorBucketMap.begin(), end = selectorBucketMap.end();
      iter != end; ++iter) {

    Selector const& sel = iter->first;
    if (sel == selector) {
      selectorBucketMap.erase(iter);
    }
  }
}

void BulkData::new_bucket_callback(EntityRank rank, const PartVector& superset_parts, size_t capacity, Bucket* new_bucket)
{
  this->new_bucket_caching(rank, new_bucket);

  if (!m_keep_fields_updated) {
    return;
  }

  const unsigned totalNumFields = mesh_meta_data().get_fields().size();

  if (m_num_fields == -1) {  // hasn't been set yet
    m_num_fields = totalNumFields;
  }

  const unsigned newBucketSize = 0;
  const FieldVector fieldsOfRank = mesh_meta_data().get_fields(rank);
  m_field_data_manager->allocate_bucket_field_data(rank, fieldsOfRank, superset_parts, totalNumFields, newBucketSize,
                                                   capacity);
}

//
//  Copy fields from src to dst entity.  If the field size of either entity is zero, do nothing.  If the field
//  size of of both entities are non-zero, then the sizes must match
//

void BulkData::copy_entity_fields_callback(EntityRank dst_rank, unsigned dst_bucket_id, unsigned dst_bucket_ord,
                                           unsigned src_bucket_id, unsigned src_bucket_ord,
                                           const std::vector<FieldBase*>* field_set)
{
  //
  //  If field set is passed in copy only the defined fields.  Also assume the fields are valid for the bucket
  //
  const FastMeshIndex srcFmi{src_bucket_id, src_bucket_ord};
  const FastMeshIndex dstFmi{dst_bucket_id, dst_bucket_ord};

  if (field_set) {
    for (const FieldBase* field : *field_set) {
      auto& srcFieldBytes = field->data_bytes<const std::byte>();
      auto& dstFieldBytes = field->data_bytes<std::byte>();

      if (field->host_data_layout() == Layout::Right) {
        auto srcEntityBytes = srcFieldBytes.entity_bytes<Layout::Right>(srcFmi);
        auto dstEntityBytes = dstFieldBytes.entity_bytes<Layout::Right>(dstFmi);
        for (stk::mesh::ByteIdx byte : srcEntityBytes.bytes()) {
          dstEntityBytes(byte) = srcEntityBytes(byte);
        }
      }
      else if (field->host_data_layout() == Layout::Left) {
        auto srcEntityBytes = srcFieldBytes.entity_bytes<Layout::Left>(srcFmi);
        auto dstEntityBytes = dstFieldBytes.entity_bytes<Layout::Left>(dstFmi);
        for (stk::mesh::ByteIdx byte : srcEntityBytes.bytes()) {
          dstEntityBytes(byte) = srcEntityBytes(byte);
        }
      }
      else {
        STK_ThrowErrorMsg("Unsupported host Field data layout: " << field->host_data_layout());
      }
    }
  }
  else {
    if (not m_keep_fields_updated) {
      return;
    }

    const std::vector<FieldBase *>& allFieldsForRank = mesh_meta_data().get_fields(dst_rank);
    for (const FieldBase* field : allFieldsForRank) {

      const int dstSize = field_bytes_per_entity(*field, dst_bucket_id);

      if (dstSize) {
        const int srcSize = field_bytes_per_entity(*field, src_bucket_id);
        if (srcSize) {
          STK_ThrowAssertMsg(dstSize == srcSize, "Incompatible field sizes: " << srcSize << " != " << dstSize);

          auto& srcFieldBytes = field->data_bytes<const std::byte>();
          auto& dstFieldBytes = field->data_bytes<std::byte>();

          if (field->host_data_layout() == Layout::Right) {
            auto srcEntityBytes = srcFieldBytes.entity_bytes<Layout::Right>(srcFmi);
            auto dstEntityBytes = dstFieldBytes.entity_bytes<Layout::Right>(dstFmi);
            for (stk::mesh::ByteIdx byte : srcEntityBytes.bytes()) {
              dstEntityBytes(byte) = srcEntityBytes(byte);
            }
          }
          else if (field->host_data_layout() == Layout::Left) {
            auto srcEntityBytes = srcFieldBytes.entity_bytes<Layout::Left>(srcFmi);
            auto dstEntityBytes = dstFieldBytes.entity_bytes<Layout::Left>(dstFmi);
            for (stk::mesh::ByteIdx byte : srcEntityBytes.bytes()) {
              dstEntityBytes(byte) = srcEntityBytes(byte);
            }
          }
          else {
            STK_ThrowErrorMsg("Unsupported host Field data layout: " << field->host_data_layout());
          }
        }
        else {
          initialize_field_on_entity(*field, dst_bucket_id, dst_bucket_ord);
        }
      }
    }
  }
}

void BulkData::add_entity_callback(EntityRank rank, unsigned bucketId, unsigned newBucketSize, unsigned bucketCapacity,
                                   unsigned indexInBucket, bool initializeFieldData)
{
  if (not m_keep_fields_updated) {
    return;
  }

  const std::vector<FieldBase*>& fieldsOfRank = mesh_meta_data().get_fields(rank);

  if (m_field_data_manager->get_bucket_capacity(rank, bucketId) < bucketCapacity) {
    m_field_data_manager->grow_bucket_capacity(fieldsOfRank, rank, bucketId, indexInBucket, bucketCapacity);
  }

#ifdef STK_ASAN_FIELD_ACCESS
  bool mustInitializeRigorously = true;
#else
  bool mustInitializeRigorously = initializeFieldData;
#endif

  if (mustInitializeRigorously) {
    m_field_data_manager->initialize_entity_field_data(fieldsOfRank, rank, bucketId, indexInBucket, newBucketSize);
  }
  else {
    m_field_data_manager->add_field_data_for_entity(fieldsOfRank, rank, bucketId, indexInBucket, newBucketSize);
  }
}

void BulkData::remove_entity_field_data_callback(EntityRank rank, unsigned bucket_id, unsigned newBucketSize,
                                                 unsigned bucket_ord)
{
    if (!m_keep_fields_updated) {
      return;
    }
    const std::vector<FieldBase*>& fieldsOfRank = mesh_meta_data().get_fields(rank);
    m_field_data_manager->remove_field_data_for_entity(rank, bucket_id, bucket_ord, newBucketSize, fieldsOfRank);
}

void BulkData::remove_entity_callback(EntityRank /*rank*/, unsigned /*bucket_id*/, unsigned /*bucket_ord*/)
{
}

void BulkData::destroy_bucket_callback(EntityRank rank, Bucket const& dying_bucket, unsigned capacity)
{
  // Remove destroyed bucket out of memoized get_buckets result, but
  // don't bother if the mesh is being destructed.
  const unsigned bucket_id = dying_bucket.bucket_id();

  if (!m_bucket_repository.being_destroyed()) {
    SelectorBucketMap& selectorBucketMap = m_selector_to_buckets_maps[rank];
    for (SelectorBucketMap::iterator itr = selectorBucketMap.begin(), end = selectorBucketMap.end();
         itr != end; ++itr) {
      Selector const& sel = itr->first;
      if (!itr->second.empty() && sel(dying_bucket)) {
        BucketVector & cached_buckets = itr->second;
        BucketVector::iterator lb_itr = std::lower_bound(cached_buckets.begin(), cached_buckets.end(), bucket_id, BucketIdComparator());
        STK_ThrowAssertMsg(lb_itr != cached_buckets.end() && (*lb_itr)->bucket_id() == bucket_id,
                       "Error, bucket id " << bucket_id << ":\n " << dying_bucket << "\nWas selected by selector " << sel << " but was not in the cache");
        cached_buckets.erase(lb_itr);
      }
    }
  }

  if (!m_keep_fields_updated) {
    return;
  }

  const std::vector<FieldBase*>& fieldsOfRank = mesh_meta_data().get_fields(rank);
  m_field_data_manager->deallocate_bucket_field_data(rank, bucket_id, capacity, fieldsOfRank);
}

void BulkData::reset_empty_field_data_callback(EntityRank rank, unsigned bucketId, unsigned bucketSize,
                                               unsigned bucketCapacity, const FieldVector& fieldsOfRank)
{
  m_field_data_manager->reset_empty_field_data(rank, bucketId, bucketSize, bucketCapacity, fieldsOfRank);
}

void BulkData::update_field_data_states(FieldBase* field, bool rotateNgpFieldViews)
{
  const int numStates = field->number_of_states();
  if (numStates > 1) {
    field->rotate_multistate_data(rotateNgpFieldViews);
  }
}

void BulkData::update_field_data_states(bool rotateNgpFieldViews)
{
  const std::vector<FieldBase*> & fields = mesh_meta_data().get_fields();

  for (int i = 0; i < m_num_fields; ) {
    FieldBase* field = fields[i];
    const int numStates = field->number_of_states();
    if (numStates > 1) {
      update_field_data_states(field, rotateNgpFieldViews);
    }
    i += numStates ;
  }
}

const_entity_iterator BulkData::begin_entities(EntityRank ent_rank) const
{
  confirm_host_mesh_is_synchronized_from_device();
  return m_entityKeyMapping->begin_rank(ent_rank);
}

const_entity_iterator BulkData::end_entities(EntityRank ent_rank) const
{
  return m_entityKeyMapping->end_rank(ent_rank);
}

Entity BulkData::get_entity( EntityRank ent_rank , EntityId entity_id ) const
{
  confirm_host_mesh_is_synchronized_from_device();
  if (!impl::is_good_rank_and_id(mesh_meta_data(), ent_rank, entity_id)) {
      return Entity();
  }
  return m_entityKeyMapping->get_entity( EntityKey(ent_rank, entity_id));
}

Entity BulkData::get_entity( const EntityKey key ) const
{
  confirm_host_mesh_is_synchronized_from_device();
  return m_entityKeyMapping->get_entity(key);
}

void BulkData::reorder_buckets_callback(EntityRank rank, const std::vector<unsigned>& reorderedBucketIds)
{
  SelectorBucketMap& selectorBucketMap = m_selector_to_buckets_maps[rank];
  for (SelectorBucketMap::iterator itr = selectorBucketMap.begin(), end = selectorBucketMap.end();
       itr != end; ++itr) {
    BucketVector& cached_buckets = itr->second;
    std::sort(cached_buckets.begin(), cached_buckets.end(), BucketIdComparator());
  }

  if (!m_keep_fields_updated) {
    return;
  }

  const FieldVector fieldsOfRank = mesh_meta_data().get_fields(rank);
  m_field_data_manager->reorder_bucket_field_data(rank, fieldsOfRank, reorderedBucketIds);
}

BucketVector const& BulkData::get_buckets(EntityRank rank, Selector const& selector) const
{
  confirm_host_mesh_is_synchronized_from_device();

  if (rank == stk::topology::INVALID_RANK) {
    static BucketVector empty;
    return empty;
  }

  STK_ThrowAssertMsg(static_cast<size_t>(rank) < m_selector_to_buckets_maps.size(), "BulkData::get_buckets, EntityRank ("<<rank<<") out of range.");

  auto& selectorBucketMap = m_selector_to_buckets_maps[rank];

  auto [it, inserted] = selectorBucketMap.try_emplace(selector, BucketVector{});
  auto& mapBuckets = it->second;

  if (inserted) {
    // Only on the first time we see this selector do the work:
    // 1) gather the full bucket list (this may have side-effects, so do it first)
    auto const& allBuckets = buckets(rank);

    for (Bucket* b : allBuckets) {
      if (selector(*b)) {
        mapBuckets.push_back(b);
      }
    }

    if (this->should_sort_buckets_by_first_entity_identifier()) {
      std::sort(mapBuckets.begin(),
                mapBuckets.end(),
                BucketIdComparator());
    }
  }

  return mapBuckets;
}

void BulkData::get_entities(EntityRank rank, Selector const& selector, EntityVector& output_entities) const {
    confirm_host_mesh_is_synchronized_from_device();
    output_entities.clear();
    const stk::mesh::BucketVector &bucket_ptrs = get_buckets(rank, selector);

    // Reserve enough space in the output_entities vector to avoid multiple reallocations
    size_t total_size = 0;
    for (const auto* bucket : bucket_ptrs) {
        total_size += bucket->size();
    }
    output_entities.reserve(total_size);

    for (const auto* bucket : bucket_ptrs) {
        output_entities.insert(output_entities.end(), bucket->begin(), bucket->end());
    }
}

//----------------------------------------------------------------------
bool BulkData::internal_declare_relation(Entity e_from, Entity e_to,
                                         RelationIdentifier local_id,
                                         Permutation permut)
{
  m_modSummary.track_declare_relation(e_from, e_to, local_id, permut);

  const MeshIndex& idx = mesh_index(e_from);

  STK_ThrowAssertMsg(local_id <= INVALID_CONNECTIVITY_ORDINAL, "local_id = " << local_id << ", max = " << (uint32_t)INVALID_CONNECTIVITY_ORDINAL);
  bool modified = idx.bucket->declare_relation(idx.bucket_ordinal, e_to, static_cast<ConnectivityOrdinal>(local_id), permut);

  if (modified)
  {
      notifier.notify_relation_declared(e_from, e_to, local_id);

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
  STK_ThrowAssertMsg(!to_entities.empty(), "ERROR, BulkData::declare_relation given empty to_entities vector.");
  EntityRank rank = entity_rank(to_entities[0]);
  STK_ThrowAssertMsg(to_entities.size() == topo.num_sub_topology(rank), "ERROR, BulkData::declare_relation given wrong number of downward relations.");
  for(unsigned i=1; i<to_entities.size(); ++i) {
    STK_ThrowAssertMsg(entity_rank(to_entities[i]) == rank, "ERROR, BulkData::declare_relation: downward relations must all have the same rank.");
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

    impl::require_valid_relation("declare", *this, e_from, e_to);

    // TODO: Don't throw if exact relation already exists, that should be a no-op.
    // Should be an exact match if relation of local_id already exists (e_to should be the same).
    bool is_new_relation = internal_declare_relation(e_from, e_to, local_id, permut);

    if(is_new_relation && m_createUpwardConnectivity)
    {
        internal_declare_relation(e_to, e_from, local_id, permut);
    }

    // It is critical that the modification be done AFTER the relations are
    // added so that the propagation can happen correctly.
    if(is_new_relation)
    {
        this->mark_entity_and_upward_related_entities_as_modified(e_to, false);
        this->mark_entity_and_upward_related_entities_as_modified(e_from, false);
    }

    // Deduce and set new part memberships:
    scratch1.clear();

    const Bucket& bucketFrom = bucket(e_from);
    impl::get_part_ordinals_to_induce_on_lower_ranks(*this, bucketFrom, entity_rank(e_to), scratch1);

    OrdinalVector emptyParts;
    internal_change_entity_parts(e_to, scratch1, emptyParts, scratch2, scratch3);
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
      STK_ThrowErrorMsg("declare_relation given entities of the same entity rank ("<<erank<<"). "
             <<entity_key(entity)<<"("<<bucket(entity).topology()<<") <--> "
             <<entity_key(e)<<"("<<bucket(e).topology()<<")");
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

  notifier.notify_relation_destroyed(e_from, e_to, local_id);

  impl::require_valid_relation( "destroy" , *this , e_from , e_to );

  const EntityRank e_to_entity_rank = entity_rank(e_to);
  const EntityRank e_from_entity_rank = entity_rank(e_from);

  m_check_invalid_rels = false; // OK to have gaps when deleting

  // When removing a relationship may need to remove part membership

  const stk::mesh::Bucket& bucketFrom = bucket(e_from);

  if ( parallel_size() < 2 || (internal_entity_comm_map_shared(e_to).empty() && bucketFrom.owned())) {
    // Only remove induced part memberships if the entity is not shared and if
    // the 'from' entity is owned.
    // (If the 'to' entity is shared then wait until modificaton_end.)

    // We will remove parts from the 'to' entity, that were induced from the
    // upward-connected 'from' entity, but only if no other upward-connected
    // entities are in those parts.

    const stk::mesh::PartVector& parts = bucketFrom.supersets();
    OrdinalVector del;

    for(const stk::mesh::Part* part : parts) {
      if (part->should_induce(e_from_entity_rank)) {
        const bool needToKeepPart = impl::part_is_on_upward_entity_except(*this, e_to, e_from_entity_rank, *part, e_from);
        if (!needToKeepPart) {
          del.push_back(part->mesh_meta_data_ordinal());
        }
      }
    }

    if ( !del.empty() ) {
      OrdinalVector empty, scratchOrdinalVec, scratchSpace;
      internal_change_entity_parts( e_to , empty, del, scratchOrdinalVec, scratchSpace);
    }
  }

  if (impl::is_valid_relation(*this, e_from, e_to, e_to_entity_rank, local_id)) {
    const bool markAuraClosureIfNotRecvGhost = false; //experiment more with this in follow-on commit: !bucket(e_from).in_aura();
    this->mark_entity_and_upward_related_entities_as_modified(e_to, markAuraClosureIfNotRecvGhost);
  }

  //delete relations from the entities
  bool caused_change_fwd = bucket(e_from).destroy_relation(e_from, e_to, local_id);

  if (caused_change_fwd && bucket(e_from).owned() && (e_from_entity_rank > e_to_entity_rank) ) {
    --m_closure_count[e_to.local_offset()];
  }


  // Relationships should always be symmetrical
  if ( caused_change_fwd ) {
    bool caused_change = bucket(e_to).destroy_relation(e_to, e_from, local_id);
    if (caused_change && bucket(e_to).owned() && (e_to_entity_rank > e_from_entity_rank) ) {
      --m_closure_count[e_from.local_offset()];
    }
  }

  m_check_invalid_rels = true;

  return caused_change_fwd;
}

void BulkData::add_sharing_info(stk::mesh::Entity entity, stk::mesh::BulkData::GhostingId ghostingId, int sharingProc)
{
    this->entity_comm_map_insert(entity, stk::mesh::EntityCommInfo(ghostingId, sharingProc));
}

bool is_modified_or_created(const BulkData& bulkData, Entity entity)
{
    return bulkData.state(entity)==stk::mesh::Modified || bulkData.state(entity)==stk::mesh::Created;
}

void BulkData::get_entities_that_have_sharing(std::vector<stk::mesh::Entity> &entitiesThatHaveSharingInfo,
        stk::mesh::EntityToDependentProcessorsMap &entityKeySharing)
{
    if(parallel_size() > 1)
    {
        // this communicates states of the entities to all procs so that entity states are consistent
        stk::mesh::EntityVector entitiesNoLongerShared;
        stk::mesh::EntityProcVec entitiesToRemoveFromSharing;
        m_meshModification.delete_shared_entities_which_are_no_longer_in_owned_closure(entitiesToRemoveFromSharing);

        impl::CommEntityMods commEntityMods(*this, internal_comm_db(), internal_comm_list());
        commEntityMods.communicate(impl::CommEntityMods::PACK_SHARED);
        m_meshModification.internal_resolve_shared_modify_delete(commEntityMods.get_shared_mods(), entitiesToRemoveFromSharing, entitiesNoLongerShared);
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
                    if(is_valid(entity) && is_modified_or_created(*this, entity))
                    {
                        if(phase == 0 && this->in_shared(entity))
                        {
                            numEntitiesThatHaveSharingInfo++;
                        }
                        else if(phase == 1 && this->in_shared(entity))
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
            STK_ThrowAssertMsg(this->is_valid(this->get_entity(key)) && this->parallel_owner_rank(this->get_entity(key)) == myProcId, "Entitykey " << key << " is not owned by receiving processor " << myProcId);
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
            this->entity_comm_map_erase(key, shared_ghosting());
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
            this->entity_comm_map_erase(entityKey, shared_ghosting());
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
        this->entity_comm_map_erase(key, shared_ghosting());
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

    stk::pack_and_communicate(bc, [&bc,&name](){
      if ( bc.parallel_rank() == 0 ) {
        bc.send_buffer().pack<char>( name.c_str() , name.size() + 1 );
      }
    });

    const char * const bc_name =
      reinterpret_cast<const char *>( bc.recv_buffer().buffer() );

    int error = 0 != std::strcmp( bc_name , name.c_str() );

    all_reduce( parallel() , ReduceMax<1>( & error ) );

    STK_ThrowErrorMsgIf( error, "Parallel name inconsistency");
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
    STK_ThrowRequireMsg(equal_case(std::string("shared"), name), "Expect shared to be the first ghosting created.");
    m_ghost_parts.push_back(&mesh_meta_data().globally_shared_part());
  }
  else if (m_ghost_parts.size() == 1) {
    STK_ThrowRequireMsg(equal_case(std::string("shared_aura"), name), "Expect aura to be the second ghosting created.");
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

  STK_ThrowRequireMsg(m_ghost_parts.size() == m_ghosting.size(), "m_ghost_parts.size()="<<m_ghost_parts.size()<<", must be same as m_ghosting.size()="<<m_ghosting.size());

  return *g ;
}

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

  EntityCommListInfoVector& commList = m_entity_comm_map.comm_list();
  for ( EntityCommListInfoVector::reverse_iterator
        i =  commList.rbegin() ;
        i != commList.rend() ; ++i) {

    if ( in_receive_ghost( i->key ) ) {
      entity_comm_map_clear_ghosting( i->key );
      bool entity_is_not_shared = !bucket(i->entity).shared();
      internal_destroy_entity_with_notification( i->entity );
      if(entity_is_not_shared)
      {
        i->key = EntityKey();
        i->entity_comm = -1;
      }
    }
    else {
      entity_comm_map_clear_ghosting(i->key);
      if ( internal_entity_comm_map(i->key).empty() ) {
        i->key = EntityKey();
        i->entity_comm = -1;
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
      const bool already_ghosted_to_proc = in_send_ghost(ghosts, add_send[i].first, add_send[i].second);
      const bool need_to_send_ghost = !ghosting_to_myself && !already_ghosted_to_proc;
      if (need_to_send_ghost)
      {
          filtered_add_send.push_back(add_send[i]);
          need_to_change_ghosting = true;
      }
    }
}

void BulkData::verify_and_filter_remove_receive(Ghosting & ghosts, const std::vector<Entity> & remove_receive, bool &need_to_change_ghosting, std::vector<Entity> & filtered_remove_receive)
{
  for (Entity entity : remove_receive) {
    if (in_receive_ghost( ghosts , entity)) {
      filtered_remove_receive.push_back(entity);
      need_to_change_ghosting = true;
    }
  }
}

bool BulkData::check_errors_and_determine_if_ghosting_needed_in_parallel(const stk::mesh::Ghosting &ghosts,
                                        bool add_send_is_owned,
                                        bool need_to_change_ghosting,
                                        const std::vector<EntityProc> & add_send)
{
    const bool ok_mesh  = &ghosts.mesh() == this;
    const bool is_custom_ghost = BulkData::AURA < ghosts.ordinal();
    int ok = ok_mesh && is_custom_ghost && add_send_is_owned;
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

      STK_ThrowErrorMsg( msg.str() );
    }

    bool anyProcsHaveNewGhosts = statuses[1] == 0;
    return anyProcsHaveNewGhosts;
}

bool BulkData::inputs_ok_and_need_ghosting(Ghosting & ghosts ,
                             const std::vector<EntityProc> & add_send ,
                             const std::vector<Entity> & remove_receive,
                             std::vector<EntityProc> &filtered_add_send,
                             std::vector<Entity> & filtered_remove_receive)
{
    bool add_send_is_owned    = true ;
    bool need_to_change_ghosting = false;

    verify_and_filter_add_send(ghosts, add_send, need_to_change_ghosting, add_send_is_owned, filtered_add_send );

    verify_and_filter_remove_receive(ghosts, remove_receive, need_to_change_ghosting, filtered_remove_receive);

    bool anyProcsHaveNewGhosts = check_errors_and_determine_if_ghosting_needed_in_parallel(ghosts, add_send_is_owned, need_to_change_ghosting, add_send);

    return anyProcsHaveNewGhosts;
}

bool BulkData::batch_add_to_ghosting(Ghosting &ghosting, const EntityProcVec &entitiesAndDestinationProcs)
{
  std::vector<stk::mesh::EntityProc> filtered_add_send;
  std::vector<stk::mesh::Entity> empty_vector;

  bool needToGhostEntities = false;
  if (inputs_ok_and_need_ghosting(ghosting , entitiesAndDestinationProcs , empty_vector, filtered_add_send, empty_vector))
  {
    needToGhostEntities = true;
    internal_batch_add_to_ghosting(ghosting, filtered_add_send);
  }
  return needToGhostEntities;
}

void BulkData::internal_batch_add_to_ghosting(Ghosting &ghosting, const EntityProcVec &entitiesAndDestinationProcs)
{
    bool starting_modification = modification_begin();
    STK_ThrowRequireMsg(starting_modification, "ERROR: BulkData already being modified,\n"
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
  std::vector<Entity> filtered_remove_receive_entities(remove_receive.size());
  for(unsigned i=0; i<remove_receive.size(); ++i) {
    remove_receive_entities[i] = get_entity(remove_receive[i]);
  }
  filtered_remove_receive_entities.clear();
  bool needToDoGhosting = inputs_ok_and_need_ghosting(ghosts, add_send , remove_receive_entities, filtered_add_send, filtered_remove_receive_entities);

  if(needToDoGhosting)
  {
    internal_change_ghosting( ghosts , filtered_add_send , filtered_remove_receive_entities );
  }
}

//----------------------------------------------------------------------

void BulkData::ghost_entities_and_fields(Ghosting & ghosting,
                                         EntityProcVec&& sendGhosts,
                                         bool isFullRegen,
                                         const std::vector<EntityProc>& removedSendGhosts)
{
    //------------------------------------
    // Push newly ghosted entities to the receivers and update the comm list.

    const int p_size = parallel_size() ;
    const bool onlyPackDownwardRelations = isFullRegen ? true : false;
    EntityCommListInfoVector newCommListEntries;
    const unsigned arbitraryInitialCapacity = 512;
    newCommListEntries.reserve(arbitraryInitialCapacity);

    stk::CommSparse commSparse( parallel() );
    for ( int phase = 0; phase < 2; ++phase ) {
      Entity prevEntity;
      for (const EntityProc& entProc : sendGhosts) {
        Entity entity = entProc.first;
        const int proc = entProc.second;

        if ( isFullRegen || !in_ghost(ghosting , entity, proc) ) {
          // Not already being sent , must send it.
          CommBuffer & buf = commSparse.send_buffer( proc );
          buf.pack<unsigned>( entity_rank(entity) );
          unsigned flag = 1;
          buf.pack<unsigned>(flag);
          impl::pack_entity_info(*this, buf, entity, onlyPackDownwardRelations);
          impl::pack_field_values(*this, buf , entity );

          if (phase == 1) {
            std::pair<int,bool> result = entity_comm_map_insert(entity, EntityCommInfo(ghosting.ordinal(), proc));
            if(result.second && entity != prevEntity) {
              const int entityCommIndex = result.first;
              EntityCommListInfo comm_info = {entity_key(entity), entity, entityCommIndex};
              newCommListEntries.push_back(comm_info);
              prevEntity = entity;
            }
          }
        }
      }

      for(const EntityProc& ep : removedSendGhosts) {
        CommBuffer& buf = commSparse.send_buffer(ep.second);
        buf.pack<unsigned>(entity_rank(ep.first));
        unsigned flag = 0;
        buf.pack<unsigned>(flag);
        buf.pack<EntityKey>(entity_key(ep.first));
      }

      if (phase == 0) {
        commSparse.allocate_buffers();
      }
      else {
        const bool deallocateSendBuffers = true;
        commSparse.communicate(deallocateSendBuffers);
      }
    }

    {
      std::vector<EntityProc>().swap(sendGhosts);
    }

    OrdinalVector ordinal_scratch, removeParts, partOrdinals, scratchSpace, scratch3;
    PartVector parts ;
    std::vector<Relation> relations ;
    std::vector<EntityProc> removedRecvGhosts;

    const MetaData & meta = mesh_meta_data() ;
    const unsigned rank_count = meta.entity_rank_count();

    // Unpacking must proceed in entity-rank order so that higher ranking
    // entities that have relations to lower ranking entities will have
    // the lower ranking entities unpacked first.  The higher and lower
    // ranking entities may be owned by different processes,
    // as such unpacking must be performed in rank order.

    std::ostringstream error_msg ;
    int error_count = 0 ;

    for ( unsigned rank = 0 ; rank < rank_count ; ++rank ) {
      for ( int p = 0 ; p < p_size ; ++p ) {
        CommBuffer & buf = commSparse.recv_buffer(p);
        while ( buf.remaining() ) {
          // Only unpack if of the current entity rank.
          // If not the current entity rank, break the iteration
          // until a subsequent entity rank iteration.
          {
            unsigned rankAndFlag[2] = {~0u,~0u};
            buf.peek<unsigned>( rankAndFlag, 2 );

            if ( rankAndFlag[1] == 1 && rankAndFlag[0] != rank ) break ;

            if (rankAndFlag[1] == 0) {
              while(buf.remaining()) {
                buf.unpack<unsigned>( rankAndFlag[0] );
                buf.unpack<unsigned>( rankAndFlag[1] );

                STK_ThrowAssert(rankAndFlag[1] == 0);
                EntityKey key;
                buf.unpack<EntityKey>(key);
                Entity rmEnt = get_entity(key);
                if (!is_valid(rmEnt)) {
                  continue;
                }
                removedRecvGhosts.push_back(EntityProc(rmEnt,p));
              }
              break;
            }

            buf.unpack<unsigned>( rankAndFlag[0] );
            buf.unpack<unsigned>( rankAndFlag[1] );
          }

          parts.clear();
          removeParts.clear();
          relations.clear();
          EntityKey key ;
          int owner = ~0u ;

          impl::unpack_entity_info( buf, *this, key, owner, parts, relations );

          if (owner != this->parallel_rank()) {
            // We will also add the entity to the part corresponding to the 'ghosts' ghosting.
            stk::mesh::Part& ghost_part = *m_ghost_parts[ghosting.ordinal()];
            insert( parts, ghost_part );
          }

          auto& ghost_reuse_map = m_meshModification.get_deleted_entity_cache().get_ghost_reuse_map();
          GhostReuseMap::iterator f_itr = ghost_reuse_map.find(key);
          const size_t use_this_offset = f_itr == ghost_reuse_map.end() ? 0 : f_itr->second;
          if (use_this_offset != 0) {
            ghost_reuse_map.erase(f_itr);
          }

          std::pair<Entity ,bool> result = internal_get_or_create_entity_with_notification( key, use_this_offset );

          Entity entity = result.first;
          const bool created   = result.second ;

          if (!created && owner != parallel_owner_rank(entity)) {
            internal_set_owner(entity, owner);
          }

          partOrdinals.clear();
          for(const stk::mesh::Part* part : parts) {
              partOrdinals.push_back(part->mesh_meta_data_ordinal());
          }

          if (!created) {
            const PartVector& currentParts = bucket(entity).supersets();
            for(const Part* currentPart : currentParts) {
              Ordinal currentPartOrdinal = currentPart->mesh_meta_data_ordinal();
              bool excludedPart = ( (currentPartOrdinal == meta.locally_owned_part().mesh_meta_data_ordinal()) ||
                   (currentPartOrdinal == meta.globally_shared_part().mesh_meta_data_ordinal()) ||
                   (!meta.get_parts()[currentPartOrdinal]->entity_membership_is_parallel_consistent() ));

              if (!excludedPart && !contains_ordinal(partOrdinals, currentPartOrdinal)) {
                removeParts.push_back(currentPartOrdinal);
              }
            }
          }

          internal_change_entity_parts_without_propagating_to_downward_connected_entities(entity, partOrdinals, removeParts, ordinal_scratch, scratchSpace, scratch3);

          if ( created ) {
            log_created_parallel_copy( entity );
          }

          const EntityCommInfo tmp( ghosting.ordinal() , owner );

          std::pair<int, bool> insertResult = entity_comm_map_insert(entity, tmp);
          if ( insertResult.second ) {
            const int entityCommIndex = insertResult.first;
            EntityCommListInfo comm_info = {entity_key(entity), entity, entityCommIndex};
            newCommListEntries.push_back(comm_info);
          }

          //now, change owner. (needed to wait until comm-info was created)
          internal_set_owner( entity, owner);

          internal_declare_relation( entity , relations, ordinal_scratch);

          if ( ! impl::unpack_field_values(*this, buf , entity , error_msg ) ) {
            ++error_count ;
          }
        }
      }
    }

#ifndef NDEBUG
    if (parallel_size() > 1) {
      all_reduce( parallel() , ReduceSum<1>( & error_count ) );
    }
#endif
    STK_ThrowRequireMsg(error_count==0, error_msg.str() );

    internal_add_comm_list_entries(newCommListEntries);

    OrdinalVector addParts, scratchOrdinalVec;
    removeParts = {ghosting_part(ghosting).mesh_meta_data_ordinal()};
    EntityCommListInfoVector& commList = m_entity_comm_map.comm_list();

    stk::util::sort_and_unique(removedRecvGhosts, EntityLess(*this));
    for(unsigned i=0; i<removedRecvGhosts.size(); ++i) {
      const unsigned reverseIdx = removedRecvGhosts.size() - i - 1;
      const EntityKey key = entity_key(removedRecvGhosts[reverseIdx].first);
      const int proc = removedRecvGhosts[reverseIdx].second;
      const bool removed = entity_comm_map_erase(key, EntityCommInfo(ghosting.ordinal(), proc));
      if (removed) {
        EntityCommListInfoVector::iterator itr = std::lower_bound(commList.begin(), commList.end(), key);
        if (itr != commList.end() && itr->key == key) {
          const int owner = parallel_owner_rank(itr->entity);
          if (owner != parallel_rank()) {
            if ( internal_entity_comm_map(itr->entity).empty() ) {
              if (is_valid(itr->entity)) {
                internal_destroy_entity_with_notification(itr->entity, true);
              }
            }
            else {
              internal_change_entity_parts(itr->entity, addParts, removeParts, scratchOrdinalVec, scratchSpace);
            }
          }
        }
      }
      else {
        Entity recvGhostEntity = get_entity(key);
        if (is_valid(recvGhostEntity) && bucket(recvGhostEntity).member(ghosting_part(ghosting))
            && !bucket(recvGhostEntity).owned()) {
          if ( internal_entity_comm_map(recvGhostEntity).empty() ) {
            const bool destroyed = internal_destroy_entity_with_notification(recvGhostEntity, true);
            if (destroyed) {
              EntityCommListInfoVector::iterator itr = std::lower_bound(commList.begin(), commList.end(), key);
              if (itr != commList.end() && itr->key == key) {
                itr->key = EntityKey();
              }
            }
          }
          else {
            internal_change_entity_parts(recvGhostEntity, addParts, removeParts, scratchOrdinalVec, scratchSpace);
          }
        }
      }
    }
    delete_unneeded_entries_from_the_comm_list();
}

void BulkData::conditionally_add_entity_to_ghosting_set(const Ghosting &ghosting, Entity entity, int toProc, EntityProcVec& entitiesWithClosure)
{
    const bool notOwnedByRecvGhostProc = toProc != parallel_owner_rank(entity);
    const bool entityIsShared = in_shared( entity_key(entity) , toProc );
    const bool alreadyGhostedToProc = in_send_ghost(ghosting, entity, toProc);
    const bool ghostingIsAura = ghosting.ordinal() == BulkData::AURA;

    bool shouldAddToGhostingSet = notOwnedByRecvGhostProc && !alreadyGhostedToProc;
    if (ghostingIsAura && entityIsShared) {
        shouldAddToGhostingSet = false;
    }

    if (shouldAddToGhostingSet)
    {
        entitiesWithClosure.push_back(EntityProc(entity , toProc));
    }
}

void BulkData::add_closure_entities(const Ghosting &ghosting, const EntityProcVec& entities, EntityProcVec& entitiesWithClosure)
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

    std::vector<EntityProc> entitiesToGhostOntoOtherProcessors;
    entitiesToGhostOntoOtherProcessors.reserve(add_send.size());

    add_closure_entities(ghosting, add_send, entitiesToGhostOntoOtherProcessors);

    stk::mesh::impl::move_unowned_entities_for_owner_to_ghost(*this, entitiesToGhostOntoOtherProcessors);

    stk::util::sort_and_unique(entitiesToGhostOntoOtherProcessors, EntityLess(*this));
    ghost_entities_and_fields(ghosting, std::move(entitiesToGhostOntoOtherProcessors));
}

void BulkData::filter_ghosting_remove_receives(const stk::mesh::Ghosting &ghosting,
                                               const std::vector <Entity> &remove_receive,
                                               std::vector<Entity> &removeRecvGhosts,
                                               std::vector<bool>& ghostStatus)
{
  for ( Entity rmEntity : remove_receive) {
    if (is_valid(rmEntity) && in_receive_ghost(ghosting, rmEntity)) {
      ghostStatus[rmEntity.local_offset()] = true;
    }
  }

  // Iterate over all entities with communication information, adding
  // the entity if it's a ghost on this process. recvGhosts will contain
  // all received-ghosts on this process by the end of the loop.
  EntityVector recvGhosts;
  recvGhosts.reserve(std::max(64u, static_cast<unsigned>(internal_comm_list().size()/4)));
  for ( const EntityCommListInfo& info : internal_comm_list()) {
    if (info.entity_comm != -1) {
      const bool inRemoveReceive = ghostStatus[info.entity.local_offset()];
      if ( is_valid(info.entity) && !inRemoveReceive && in_receive_ghost(ghosting, info.entity) ) {
        recvGhosts.push_back(info.entity);
        ghostStatus[info.entity.local_offset()] = true;
      }
      else {
        ghostStatus[info.entity.local_offset()] = false;
      }
    }
  }

  //Add back in the closure-entities of each recv-ghost, if those closure-entities were
  //removed due to being in the remove_receive list
  impl::OnlyRecvGhosts org(*this,ghosting,ghostStatus);
  impl::VecPushBack vpb(recvGhosts, ghostStatus);

  unsigned len = recvGhosts.size();
  for (unsigned ii=0; ii<len; ++ii) {
    Entity e = recvGhosts[ii];
    const EntityRank erank = entity_rank(e);

    for (EntityRank irank = stk::topology::BEGIN_RANK; irank < erank; ++irank) {
      const ConnectedEntities rels = get_connected_entities(e, irank);
      for (unsigned j=0; j<rels.size(); ++j) {
        if (irank > stk::topology::ELEM_RANK) {
          impl::VisitClosureGeneral(*this, rels[j], irank, vpb, org);
        }
        else {
          if ( is_valid(rels[j]) &&
               in_receive_ghost( ghosting , rels[j] ) && !ghostStatus[rels[j].local_offset()] )
          {
            recvGhosts.push_back(rels[j]);
            ghostStatus[rels[j].local_offset()] = true;
          }
        }
      }
    }
  }

  removeRecvGhosts.clear();
  removeRecvGhosts.reserve(remove_receive.size());
  for(Entity entity : remove_receive) {
    if (in_receive_ghost(ghosting, entity) && !ghostStatus[entity.local_offset()]) {
      removeRecvGhosts.push_back(entity);
    }
  }
}

void BulkData::delete_unneeded_entries_from_the_comm_list()
{
    EntityCommListInfoVector& commList = m_entity_comm_map.comm_list();
    EntityCommListInfoVector::iterator i =
      std::remove_if( commList.begin() ,
                      commList.end() , IsInvalid() );
    commList.erase( i , commList.end() );
}

void BulkData::internal_change_ghosting(
  Ghosting & ghosting ,
  const std::vector<EntityProc> & add_send ,
  const std::vector<Entity> & remove_receive,
  bool add_send_is_globally_empty)
{
  m_modSummary.track_change_ghosting(ghosting, add_send, remove_receive);

  // put add_send entities and their closure in newSendGhosts
  EntityProcVec newSendGhosts;
  impl::StoreInEntityProcVec storeInVec(*this, newSendGhosts);
  impl::OnlyGhosts og(*this);
  for ( const EntityProc& entityProc : add_send ) {
      storeInVec.proc = entityProc.second;
      og.proc = entityProc.second;
      impl::VisitClosureGeneral(*this,entityProc.first,entity_rank(entityProc.first),storeInVec,og);
  }
  stk::util::sort_and_unique(newSendGhosts, EntityLess(*this));

  //remove newSendGhosts that are already in comm-list:
  for (EntityProc& sendGhost : newSendGhosts) {
    if (in_send_ghost(ghosting, sendGhost.first, sendGhost.second)) {
      sendGhost.first = Entity();
    }
  }
  auto shouldRemove = [&](const EntityProc& ep) {
    return ep.first.local_offset() == 0;
  };
  newSendGhosts.erase(
       std::remove_if(newSendGhosts.begin(), newSendGhosts.end(), shouldRemove),
       newSendGhosts.end());

  std::vector<Entity> removeRecvGhosts;
  std::vector<bool> ghostStatus(get_size_of_entity_index_space(), false);
  filter_ghosting_remove_receives(ghosting, remove_receive, removeRecvGhosts, ghostStatus);

  std::set<EntityKeyProc> removeSendGhosts;

  stk::mesh::impl::comm_sync_send_recv(*this , removeRecvGhosts, newSendGhosts, removeSendGhosts);

  if (!add_send_is_globally_empty) {
    stk::util::sort_and_unique(newSendGhosts, EntityLess(*this));
    ghost_entities_and_fields(ghosting, std::move(newSendGhosts), false);
  }
  else {
    STK_ThrowRequireMsg(newSendGhosts.empty(), "internal_change_ghosting: add_send_is_globally_empty but newSendGhosts not empty");
  }

  const unsigned ghostingPartOrdinal = ghosting_part(ghosting).mesh_meta_data_ordinal();
  OrdinalVector removeGhostingPart = {ghostingPartOrdinal};
  const unsigned ownedPartOrdinal = mesh_meta_data().locally_owned_part().mesh_meta_data_ordinal();
  OrdinalVector addOwnedPart = {ownedPartOrdinal};
  OrdinalVector addParts;
  OrdinalVector scratchOrdinalVec, scratchSpace;

  for (std::set<EntityKeyProc>::reverse_iterator reverseIterator = removeSendGhosts.rbegin(); reverseIterator != removeSendGhosts.rend(); ++reverseIterator) {
    const EntityKey key = reverseIterator->first;
    Entity entity = get_entity(key);
    const int proc = reverseIterator->second;
    if (!is_valid(entity) ||
        impl::has_upward_send_ghost_connectivity(*this, ghosting, proc, get_entity(key))) {
      continue;
    }

    entity_comm_map_erase(key, EntityCommInfo(ghosting.ordinal(), proc));
    internal_change_entity_parts(get_entity(key), {}, removeGhostingPart, scratchOrdinalVec, scratchSpace);
  }

  const EntityCommDatabase& commDB = internal_comm_db();
  for(EntityCommListInfo& info : m_entity_comm_map.comm_list()) {
    if (info.entity_comm == -1 || commDB.comm(info.entity_comm).empty()) {
      info.key = EntityKey();
    }
  }
  delete_unneeded_entries_from_the_comm_list();

  stk::util::sort_and_unique(removeRecvGhosts, EntityLess(*this));

  for (unsigned i=0; i<removeRecvGhosts.size(); ++i) {
    const unsigned reverseIdx = removeRecvGhosts.size() - i - 1;
    Entity rmEntity = removeRecvGhosts[reverseIdx];
    const EntityKey key = entity_key(rmEntity);
    if (impl::has_upward_recv_ghost_connectivity(*this, ghosting, rmEntity)) {
      continue;
    }

    entity_comm_map_erase(key, ghosting);

    if (internal_entity_comm_map(rmEntity).empty()) {
      const bool destroyed = internal_destroy_entity_with_notification(rmEntity, true);
      if (!destroyed && owned_closure(rmEntity)) {
        internal_set_owner(rmEntity, parallel_rank());
        internal_change_entity_parts(rmEntity, addOwnedPart, removeGhostingPart, scratchOrdinalVec, scratchSpace);
      }
    }
    else {
      internal_change_entity_parts(rmEntity, addParts, removeGhostingPart, scratchOrdinalVec, scratchSpace);
    }
  }

  for(EntityCommListInfo& info : m_entity_comm_map.comm_list()) {
    if (info.entity_comm == -1 || commDB.comm(info.entity_comm).empty()) {
      info.key = EntityKey();
    }
  }
  delete_unneeded_entries_from_the_comm_list();
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void BulkData::internal_regenerate_aura()
{
  require_ok_to_modify();

  m_auraGhosting->generate_aura(*this);
}

void BulkData::internal_remove_aura()
{
  require_ok_to_modify();

  if (m_autoAuraOption == BulkData::NO_AUTO_AURA) {
    m_auraGhosting->remove_aura(*this);
    m_turningOffAutoAura = false;
  }
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
            STK_ThrowRequireWithSierraHelpMsg(new_permutation!=stk::mesh::INVALID_PERMUTATION);

            unsigned edges_element_offset = static_cast<stk::mesh::ConnectivityOrdinal>(i);
            bucket_edge.change_existing_permutation_for_connected_element(bucket_ordinal(side), edges_element_offset, new_permutation);

            unsigned elements_edge_offset = get_index_of_side_in_element_bucket(*this, elements[i], side);
            stk::mesh::Bucket& bucket_elem = bucket(elements[i]);
            bucket_elem.change_existing_permutation_for_connected_edge(bucket_ordinal(elements[i]), elements_edge_offset, new_permutation);
        }
    }
}

void BulkData::resolve_parallel_side_connections(stk::mesh::EntityRank rank,
                                                 std::vector<SideSharingData>& sideSharingDataToSend,
                                                 std::vector<SideSharingData>& sideSharingDataReceived)
{
    std::vector<stk::mesh::impl::IdViaSidePair> idAndSides;

    stk::mesh::EntityVector sharedSides = fill_shared_entities_that_need_fixing(*this, rank);

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
                stk::mesh::EntityKey newKey(rank, sideSharingData.chosenSideId);
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

void BulkData::use_elem_elem_graph_to_determine_shared_entities(std::vector<stk::mesh::Entity>& shared_entities,
                                                                stk::mesh::EntityRank rank)
{
    std::vector<SideSharingData> sideSharingDataReceived;
    std::vector<SideSharingData> sideSharingDataToSend;

    resolve_parallel_side_connections(rank, sideSharingDataToSend, sideSharingDataReceived);

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
    if (this->has_face_adjacent_element_graph() && (rank == stk::topology::EDGE_RANK || rank == stk::topology::FACE_RANK))
    {
        use_elem_elem_graph_to_determine_shared_entities(shared_new, rank);
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
  if (parallel_size() > 1) {
    std::vector<Entity> shared_nodes;
    this->gather_shared_nodes(shared_nodes);

    shared_new.clear();
    shared_new.reserve(shared_nodes.size()*2);
    shared_new.insert(shared_new.end(), shared_nodes.begin(), shared_nodes.end());

    this->fill_shared_entities_of_rank_while_updating_sharing_info(stk::topology::EDGE_RANK, shared_new);
    this->fill_shared_entities_of_rank_while_updating_sharing_info(stk::topology::FACE_RANK, shared_new);

  }
  std::fill(m_mark_entity.begin(), m_mark_entity.end(), BulkData::NOT_MARKED);
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
    if(in_receive_ghost(relation)) {
      ghosts.push_back(relation);
    }
  });

  return ghosts;
}

EntityVector BulkData::get_upward_send_ghost_relations(const Entity entity)
{
  EntityVector ghosts;
  filter_upward_ghost_relations(entity, [&](Entity relation) {
    if(in_send_ghost(relation)) {
      ghosts.push_back(relation);
    }
  });

  return ghosts;
}

void BulkData::resolve_ownership_of_modified_entities( const std::vector<Entity> &shared_modified )
{
  const BulkData& bulk = *this;
  stk::CommSparse comm_sparse( parallel() );
  const bool anythingToUnpack =
    stk::pack_and_communicate(comm_sparse, [&comm_sparse, &bulk, &shared_modified]() {
      for (Entity entity : shared_modified) {
        if (bulk.parallel_owner_rank(entity) == bulk.parallel_rank() &&
            bulk.state(entity) != Created )
        {
          PairIterEntityComm jc = bulk.internal_entity_comm_map_shared(entity);
          for (; !jc.empty(); ++jc) {
            comm_sparse.send_buffer(jc->proc).pack<EntityKey>( bulk.entity_key(entity) );
          }
        }
      }
    });

  if (anythingToUnpack) {
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
}

void BulkData::move_entities_to_proper_part_ownership( const std::vector<Entity> &shared_modified )
{
    std::ostringstream error_msg;
    int error_flag = 0;

    OrdinalVector shared_part, owned_part, empty;
    OrdinalVector scratchOrdinalVec, scratchSpace;
    shared_part.push_back(mesh_meta_data().globally_shared_part().mesh_meta_data_ordinal());
    owned_part.push_back(mesh_meta_data().locally_owned_part().mesh_meta_data_ordinal());

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
        else if(!internal_entity_comm_map_shared(entity).empty())
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
            error_msg << "    " << entity_key(entity);
            error_msg << " also declared on";
            for(PairIterEntityComm ec = internal_entity_comm_map_shared(entity); !ec.empty(); ++ec)
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
    STK_ThrowErrorMsgIf(error_flag, error_msg.str());
}


void BulkData::add_comm_list_entries_for_entities(const std::vector<stk::mesh::Entity>& sharedModifiedEntities)
{
    // ------------------------------------------------------------
    // Update comm-list based on shared_modified

    EntityCommListInfoVector& commList = m_entity_comm_map.comm_list();
    const size_t n_old = commList.size();

    commList.reserve(commList.size() + sharedModifiedEntities.size());
    for (size_t i = 0, e = sharedModifiedEntities.size(); i < e; ++i)
    {
      Entity entity = sharedModifiedEntities[i];
      EntityKey key = entity_key(entity);
      const int entityCommIndex = m_entity_comm_map.entity_comm(key);
      EntityCommListInfo comm_info = {key, entity, entityCommIndex};
      commList.push_back(comm_info);
    }

    std::sort( commList.begin() + n_old , commList.end() );

    std::inplace_merge( commList.begin() ,
                        commList.begin() + n_old ,
                        commList.end() );

    EntityCommListInfoVector::iterator iter =
      std::unique( commList.begin() , commList.end() );

    commList.erase( iter , commList.end() );
}

void BulkData::internal_add_comm_list_entries(EntityCommListInfoVector& newCommListEntries)
{
  EntityCommListInfoVector& commList = m_entity_comm_map.comm_list();
  if (!newCommListEntries.empty()) {
    stk::util::sort_and_unique(newCommListEntries);
    if (newCommListEntries.capacity() >= newCommListEntries.size()+commList.size()) {
      commList.swap(newCommListEntries);
    }
    stk::util::insert_keep_sorted_and_unique(newCommListEntries, commList);
  }
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
    std::vector<stk::mesh::EntityRank> ranks;
    for(EntityRank rank=stk::topology::NODE_RANK; rank<=mesh_meta_data().side_rank(); rank++)
        ranks.push_back(rank);
    this->internal_resolve_parallel_create(ranks);
}

void BulkData::internal_resolve_parallel_create(const std::vector<stk::mesh::EntityRank>& ranks)
{
  STK_ThrowRequireMsg(parallel_size() > 1, "Do not call this in serial");

  std::vector<Entity> shared_modified;
  for (EntityRank rank : ranks) {
    internal_update_sharing_comm_map_and_fill_list_modified_shared_entities_of_rank(rank, shared_modified);
  }

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

  update_comm_list_based_on_changes_in_comm_map();
  std::fill(m_mark_entity.begin(), m_mark_entity.end(), BulkData::NOT_MARKED);
}

bool BulkData::modification_end_after_node_sharing_resolution()
{
  notifier.notify_started_modification_end();
  return m_meshModification.modification_end_after_node_sharing_resolution();
}

bool BulkData::modification_end_for_entity_creation( const std::vector<EntityRank> & entity_rank_vector, ModEndOptimizationFlag opt)
{
  notifier.notify_started_modification_end();

  bool return_value = internal_modification_end_for_entity_creation( entity_rank_vector, opt );
  return return_value;
}

void BulkData::modification_end_for_sync_to_host(ModEndOptimizationFlag opt)
{
  for (EntityRank rank=stk::topology::BEGIN_RANK; rank <= stk::topology::END_RANK; ++rank)
  {
    notifier.notify_local_buckets_changed(rank);
  }

  notifier.notify_started_modification_end();

  if (opt == ModEndOptimizationFlag::MOD_END_SORT)
  {
    m_bucket_repository.internal_default_sort_bucket_entities(should_sort_faces_by_node_ids());
  }

  m_bucket_repository.internal_modification_end();

  m_meshModification.set_sync_state_synchronized();
  for (SelectorBucketMap& selectorBucketMap : m_selector_to_buckets_maps)
  {
    selectorBucketMap.clear();
  }

  if (get_maintain_local_ids())
  {
    impl::set_local_ids(*this);
  }

  notify_finished_mod_end();

  if(parallel_size() > 1)
  {
    STK_ThrowRequireMsg(false, "parallel mesh sync to host not supported");
    check_mesh_consistency();
  }
}

void BulkData::update_comm_list_based_on_changes_in_comm_map()
// Resolution of shared and ghost modifications can empty
// the communication information for entities.
// If there is no communication information then the
// entity must be removed from the communication list.
{
  bool changed = false ;
  for (EntityCommListInfo& entityCommInfo : m_entity_comm_map.comm_list()) {
      if (entityCommInfo.entity_comm == -1) {
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
  bool anyModOnAnyProc = notifier.notify_finished_modification_end(parallel());
  if (!anyModOnAnyProc) {
    m_meshModification.set_sync_count(synchronized_count()-1);
  }
}

void BulkData::internal_modification_end_for_change_ghosting()
{
    internal_resolve_send_ghost_membership();

    m_modSummary.write_summary(m_meshModification.synchronized_count());

    internal_finish_modification_end(ModEndOptimizationFlag::MOD_END_SORT);
}

bool BulkData::internal_modification_end_for_change_parts(ModEndOptimizationFlag opt)
{
    notifier.notify_started_modification_end();

    int global_shared_modified = 0;
    int local_shared_modified = 0;
    if (m_meshModification.did_any_shared_entity_change_parts())
    {
        local_shared_modified = 1;
    }

    stk::all_reduce_max(parallel(), &local_shared_modified, &global_shared_modified, 1);

    if (this->parallel_size() > 1 && global_shared_modified > 0) {
      stk::mesh::EntityVector entitiesNoLongerShared;
      stk::mesh::EntityProcVec entitiesToRemoveFromSharing;
      m_meshModification.delete_shared_entities_which_are_no_longer_in_owned_closure(entitiesToRemoveFromSharing);

      impl::CommEntityMods commEntityMods(*this, internal_comm_db(), internal_comm_list());
      commEntityMods.communicate(impl::CommEntityMods::PACK_SHARED);
      m_meshModification.internal_resolve_shared_modify_delete(commEntityMods.get_shared_mods(), entitiesToRemoveFromSharing, entitiesNoLongerShared);
      internal_resolve_shared_membership(entitiesNoLongerShared);
    }
    m_modSummary.write_summary(m_meshModification.synchronized_count());

    constexpr bool resetSymGhostInfo = false;
    internal_finish_modification_end(opt, resetSymGhostInfo);

    return true;
}

bool BulkData::internal_modification_end_for_change_entity_owner(ModEndOptimizationFlag opt )
{
  // The two states are MODIFIABLE and SYNCHRONiZED
  if ( this->in_synchronized_state() ) { return false ; }

  STK_ThrowAssertMsg(impl::check_for_connected_nodes(*this)==0, "BulkData::modification_end ERROR, all entities with rank higher than node are required to have connected nodes.");

  STK_ThrowAssertMsg(add_fmwk_data() || impl::check_no_shared_elements_or_higher(*this)==0, "BulkData::modification_end ERROR, Sharing of entities with rank ELEMENT_RANK or higher is not allowed.");

  if (parallel_size() > 1)
  {
      if ( m_autoAuraOption == AUTO_AURA )
      {
          internal_regenerate_aura();
      }
      else if (m_turningOffAutoAura) {
          internal_remove_aura();
      }

      internal_resolve_send_ghost_membership();
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

  if(m_runConsistencyCheck) {
    STK_ThrowErrorMsgIf(!stk::mesh::impl::check_permutations_on_all(*this), "Permutation checks failed.");
    std::ostringstream msg ;
    bool is_consistent = impl::comm_mesh_verify_parallel_consistency(*this, internal_comm_db(), internal_comm_list(), [&](Entity entity){return internal_entity_comm_map(entity);}, msg );
    std::string error_msg = msg.str();
    STK_ThrowErrorMsgIf( !is_consistent, error_msg );
  }
}


//////////////////////////////////// Free functions to help with modification end (exp) for edges

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

std::vector<bool> get_allowable_ghost_connections(stk::mesh::BulkData &stkMeshBulkData,
                                                  std::vector<stk::mesh::Entity> &entitiesConnectedToNodes,
                                                  const stk::mesh::Entity /*sideEntity*/)
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

void BulkData::determineEntitiesThatNeedGhosting(stk::mesh::Entity edge,
                                                 std::vector<stk::mesh::Entity>& entitiesConnectedToNodes,
                                                 const stk::mesh::Entity* /*nodes*/,
                                                 EntityProcVec& addGhostedEntities)
{
    // Grab all the entities attached to the 2 nodes
    // If the entity is ghosted and the edge is owned, then the edge needs to be ghosted.
    bool doesEdgeNeedToBeGhosted = doesEdgeNeedGhostingCommunication(*this, entitiesConnectedToNodes);
    if ( doesEdgeNeedToBeGhosted )
    {
        for (size_t j=0;j<entitiesConnectedToNodes.size();j++)
        {
            if ( entitiesConnectedToNodes[j] != Entity() )
            {
                PairIterEntityComm ghosted = internal_entity_comm_map( entitiesConnectedToNodes[j] , aura_ghosting());
                for (PairIterEntityComm ec = ghosted; !ec.empty(); ++ec)
                {
                    if ( ec->proc != parallel_rank() )
                    {
                        bool isEdgeSharedWithOtherProc = in_shared(entity_key(edge), ec->proc);
                        if ( !isEdgeSharedWithOtherProc )
                        {
                            addGhostedEntities.push_back(EntityProc(edge, ec->proc));
                        }
                    }
                }
            }
        }
    }
}

void BulkData::find_upward_connected_entities_to_ghost_onto_other_processors(EntityProcVec& entitiesToGhostOntoOtherProcessors,
                                                                             EntityRank entity_rank,
                                                                             const stk::mesh::Selector& selected,
                                                                             bool connectFacesToPreexistingGhosts)
{
    if(entity_rank == stk::topology::NODE_RANK) { return; }

    const stk::mesh::BucketVector& entity_buckets = buckets(entity_rank);
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
            if ( is_valid(entity) && state(entity) != Unchanged )
            {
                const stk::mesh::Entity* nodes = begin_nodes(entity);

                if(isedge)
                {
                    impl::find_entities_these_nodes_have_in_common(*this, stk::topology::FACE_RANK,
                          numNodes, nodes, facesConnectedToNodes);
                    removeEntitiesNotSelected(*this, selected, facesConnectedToNodes);
                    if(connectFacesToPreexistingGhosts)
                        connectGhostedEntitiesToEntity(*this, facesConnectedToNodes, entity, nodes);
                }

                impl::find_entities_these_nodes_have_in_common(*this, stk::topology::ELEM_RANK,
                            numNodes, nodes, elementsConnectedToNodes);
                removeEntitiesNotSelected(*this, selected, elementsConnectedToNodes);
                if(connectFacesToPreexistingGhosts)
                    connectGhostedEntitiesToEntity(*this, elementsConnectedToNodes, entity, nodes);

                if ( bucket.owned() || bucket.shared() )
                {
                    if (isedge)
                    {
                      determineEntitiesThatNeedGhosting(entity, facesConnectedToNodes, nodes, entitiesToGhostOntoOtherProcessors);
                    }
                    determineEntitiesThatNeedGhosting(entity, elementsConnectedToNodes, nodes, entitiesToGhostOntoOtherProcessors);
                }
            }
        }
    }
    stk::util::sort_and_unique(entitiesToGhostOntoOtherProcessors, EntityLess(*this));

    std::set< EntityKey > entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs;
    stk::mesh::impl::comm_sync_send_recv(*this, entitiesToGhostOntoOtherProcessors, entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs);
}

void BulkData::make_ghost_info_symmetric()
{
  remove_symmetric_ghost_info();

  const EntityCommDatabase& commDB = m_entity_comm_map;
  const EntityCommListInfoVector& all_comm = commDB.comm_list();

  stk::CommSparse commSparse(parallel());
  stk::pack_and_communicate(commSparse,
    [this, &commSparse, &commDB, all_comm]() {
      for(size_t i=0; i<all_comm.size(); ++i) {

        if (this->parallel_owner_rank(all_comm[i].entity) == this->parallel_rank()) {
          PairIterEntityComm ec = commDB.comm(all_comm[i].entity_comm);
          PairIterEntityComm ecSh = shared_comm_info_range(ec);

          const bool sharedAndGhostedOrMultipleGhosted =
                 (ecSh.size() > 0 && ecSh.size() < ec.size()) ||
                 (ecSh.size() == 0 && ec.size() > 1);

          if (sharedAndGhostedOrMultipleGhosted) {
            for(unsigned ii=0; ii<ec.size(); ++ii) {
              const int proc_ii = ec[ii].proc;
              commSparse.send_buffer(proc_ii).pack<EntityKey>(all_comm[i].key);
              commSparse.send_buffer(proc_ii).pack<unsigned>(ec.size()-1);

              for(unsigned jj=0; jj<ec.size(); ++jj) {
                if (ii != jj) {
                  const int proc_jj = ec[jj].proc;
                  commSparse.send_buffer(proc_ii).pack<unsigned>(ec[jj].ghost_id);
                  commSparse.send_buffer(proc_ii).pack<int>(proc_jj);
                }
              }
            }
          }
        }
      }
    });

  EntityCommListInfoVector newCommListEntries;

  for(int p=0; p<parallel_size(); ++p) {
    CommBuffer& buf = commSparse.recv_buffer(p);
    while(buf.remaining()) {
      EntityKey key;
      buf.unpack<EntityKey>(key);

      unsigned numOtherProcs = 0;
      buf.unpack<unsigned>(numOtherProcs);

      for(unsigned i=0; i<numOtherProcs; ++i) {
        unsigned ghostId = 0;
        buf.unpack<unsigned>(ghostId);
        int otherProc = -1;
        buf.unpack<int>(otherProc);

        if (otherProc != parallel_rank()) {
          EntityCommInfo entityCommInfo(SYMM_INFO+ghostId, otherProc);

          constexpr int ownerArgIsNotUsed = -1;
          m_entity_comm_map.insert(key, entityCommInfo, ownerArgIsNotUsed);
        }
      }
    }
  }
}

void BulkData::remove_symmetric_ghost_info()
{
  EntityCommDatabase& commDB = m_entity_comm_map;
  const EntityCommListInfoVector& all_comm = commDB.comm_list();
  for(size_t i=0; i<all_comm.size(); ++i) {
    if (all_comm[i].entity_comm != -1) {
      commDB.erase(all_comm[i].entity_comm, all_comm[i].key, SYMM_INFO);
    }
  }
}


void BulkData::internal_finish_modification_end(ModEndOptimizationFlag opt,
                                                bool resetSymGhostInfo)
{
  if(opt == impl::MeshModification::MOD_END_SORT) {
    m_bucket_repository.internal_default_sort_bucket_entities(should_sort_faces_by_node_ids());
  }

  m_bucket_repository.internal_modification_end();

  m_meshModification.set_sync_state_synchronized();
  m_add_node_sharing_called = false;

  m_meshModification.get_deleted_entity_cache().update_deleted_entities_container();

  for(SelectorBucketMap& selectorBucketMap : m_selector_to_buckets_maps) {
    for (SelectorBucketMap::iterator itr = selectorBucketMap.begin(), end = selectorBucketMap.end(); itr != end; ++itr) {
      if (itr->second.empty()) {
        itr = selectorBucketMap.erase(itr);
        if (itr == end) {
          break;
        }
      }
    }
  }

  if (resetSymGhostInfo && has_symmetric_ghost_info() && parallel_size() > 2) {
    make_ghost_info_symmetric();
  }

  if (get_maintain_local_ids()) {
    impl::set_local_ids(*this);
  }

  notify_finished_mod_end();

  if(parallel_size() > 1) {
    check_mesh_consistency();
  }
}

bool BulkData::internal_modification_end_for_skin_mesh( EntityRank entity_rank, ModEndOptimizationFlag opt, const stk::mesh::Selector& selectedToSkin,
        const Selector * only_consider_second_element_from_this_selector)
{
  // The two states are MODIFIABLE and SYNCHRONiZED
  if ( this->in_synchronized_state() ) { return false ; }

  STK_ThrowAssertMsg(impl::check_for_connected_nodes(*this)==0, "BulkData::modification_end ERROR, all entities with rank higher than node are required to have connected nodes.");

  if (parallel_size() > 1)
  {
      if ( !this->is_automatic_aura_on())
      {
          find_and_delete_internal_faces(entity_rank, only_consider_second_element_from_this_selector);
      }

      stk::mesh::EntityVector entitiesNoLongerShared;
      stk::mesh::EntityProcVec entitiesToRemoveFromSharing;
      m_meshModification.delete_shared_entities_which_are_no_longer_in_owned_closure(entitiesToRemoveFromSharing);

      impl::CommEntityMods commEntityMods(*this, internal_comm_db(), internal_comm_list());
      commEntityMods.communicate(impl::CommEntityMods::PACK_SHARED);
      m_meshModification.internal_resolve_shared_modify_delete(commEntityMods.get_shared_mods(), entitiesToRemoveFromSharing, entitiesNoLongerShared);
      this->internal_resolve_shared_membership(entitiesNoLongerShared);

      if(is_automatic_aura_on())
      {
          bool connectFacesToPreexistingGhosts = true;
          resolve_incremental_ghosting_for_entity_creation_or_skin_mesh(entity_rank, selectedToSkin, connectFacesToPreexistingGhosts);
      }
  }
  m_modSummary.write_summary(m_meshModification.synchronized_count());

  this->internal_finish_modification_end(opt);

  return true ;
}

void BulkData::resolve_incremental_ghosting_for_entity_creation_or_skin_mesh(EntityRank entity_rank, const stk::mesh::Selector& selectedToSkin, bool connectFacesToPreexistingGhosts)
{
    EntityProcVec sendGhosts;
    find_upward_connected_entities_to_ghost_onto_other_processors(sendGhosts, entity_rank, selectedToSkin, connectFacesToPreexistingGhosts);

    stk::util::sort_and_unique(sendGhosts, EntityLess(*this));
    ghost_entities_and_fields(aura_ghosting(), std::move(sendGhosts));
}

bool BulkData::internal_modification_end_for_entity_creation( const std::vector<EntityRank> & entity_rank_vector, ModEndOptimizationFlag opt )
{
  // The two states are MODIFIABLE and SYNCHRONiZED
  if (in_synchronized_state()) { return false; }

  STK_ThrowAssertMsg(impl::check_for_connected_nodes(*this)==0, "BulkData::modification_end ERROR, all entities with rank higher than node are required to have connected nodes.");

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

void BulkData::fill_entity_procs_for_owned_modified_or_created(std::vector<EntityProc> & send_list ) const
{
    const int p_rank = this->parallel_rank();
    const EntityCommDatabase& commDB = internal_comm_db();
    for(const EntityCommListInfo& info : commDB.comm_list())
    {
        stk::mesh::Entity entity = info.entity;
        const int owner = parallel_owner_rank(entity);

        if(owner == p_rank && is_modified_or_created(*this, entity))
        {
            const int entityCommIndex = info.entity_comm;
            for(PairIterEntityComm ec = commDB.comm(entityCommIndex); !ec.empty(); ++ec)
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


void BulkData::remove_unneeded_induced_parts(stk::mesh::Entity entity, PairIterEntityComm entity_comm_info,
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
    STK_ThrowRequireMsg(parallel_size() > 1, "Do not call this in serial");

    EntityProcVec auraEntitiesToResend;
    EntityVector newlyShared;
    {
      stk::CommSparse comm(parallel());
      const bool anythingToUnpack = impl::pack_and_send_modified_shared_entity_states(
                                      comm, *this, m_entity_comm_map.comm_list());

      if (anythingToUnpack) {
        for(int p=0; p<parallel_size(); ++p) {
          stk::CommBuffer& buf = comm.recv_buffer(p);
          while(buf.remaining()) {
            EntityKey key;
            EntityState remoteState;
            buf.unpack<EntityKey>(key);
            buf.unpack<EntityState>(remoteState);
            Entity entity = get_entity(key);
            int isShared = 0;
            buf.unpack<int>(isShared);
            if (isShared == 1) {
              auto rslt = entity_comm_map_insert(entity, EntityCommInfo(SHARED, p));
              if (rslt.second){
                newlyShared.push_back(entity);
              }
            }
            int numSharingProcs = 0;
            buf.unpack<int>(numSharingProcs);
            for(int sp=0; sp<numSharingProcs; ++sp) {
              int sharingProc = -1;
              buf.unpack<int>(sharingProc);
              entity_comm_map_insert(entity, EntityCommInfo(SHARED, sharingProc));
            }

            const EntityState localState = state(entity);
            if (!in_shared(entity,p)) {
              if (localState == Unchanged) {
                stk::util::insert_keep_sorted_and_unique(EntityProc(entity,p),
                                        auraEntitiesToResend, EntityLess(*this));
              }
              continue;
            }

            if (localState == Unchanged && remoteState != Unchanged) {
              set_state(entity, Modified);
            }
          }
        }
      }
    }

    ParallelMachine p_comm = parallel();
    const int p_rank = parallel_rank();

    stk::CommSparse comm(p_comm);
    const EntityCommDatabase& commDB = internal_comm_db();
    const EntityCommListInfoVector& commList = commDB.comm_list();
    impl::pack_and_send_induced_parts_from_sharers_to_owners(*this, comm, commList);

    PartStorage part_storage;

#ifndef NDEBUG
    bool localOk = true;
#endif
    std::string errorMsg;
    try
    {
      std::vector<bool> shouldProcess(commList.size(), false);
      for(unsigned i=0; i<commList.size(); ++i) {
        const EntityCommListInfo& info = commList[i];
        const int owner = parallel_owner_rank(info.entity);
        bool i_own_this_entity = (owner == p_rank);

        if (i_own_this_entity && state(info.entity) != Unchanged) {
          shouldProcess[i] = true;
        }
      }
      for (unsigned i=0; i<commList.size(); ++i) {
        if (shouldProcess[i]) {
          const EntityCommListInfo& info = commList[i];
          remove_unneeded_induced_parts(info.entity, commDB.comm(info.entity_comm), part_storage, comm);
        }
      }
    }
    catch(std::exception& e) {
      errorMsg = "P"+std::to_string(parallel_rank())
               + " stk::mesh::BulkData::internal_resolve_shared_membership exception: "
               + e.what();
#ifndef NDEBUG
      localOk = false;
#else
      std::ostringstream os;
      os<<errorMsg<<std::endl;
      os<<"Run a debug build to get parallel-consistent catch and shutdown."<<std::endl;
      std::cerr<<os.str();
#endif
    }

#ifndef NDEBUG
    if (!stk::is_true_on_all_procs(this->parallel(), localOk)) {
      std::string msg = localOk ? "error on another proc" : errorMsg;
      STK_ThrowErrorMsg(msg);
    }
#endif

    OrdinalVector scratch, scratch2;
    part_storage.induced_part_ordinals.clear();
    part_storage.removeParts.clear();
    part_storage.induced_part_ordinals.push_back(mesh_meta_data().globally_shared_part().mesh_meta_data_ordinal());
    OrdinalVector& sharedPart = part_storage.induced_part_ordinals;

    for(Entity ent : newlyShared) {
      internal_change_entity_parts(ent, sharedPart, part_storage.removeParts, scratch, scratch2);
    }

    for (stk::mesh::Entity entity : entitiesNoLongerShared) {
      part_storage.induced_part_ordinals.clear();
      part_storage.removeParts.clear();
      induced_part_membership(*this, entity, part_storage.induced_part_ordinals);
      impl::filter_out_unneeded_induced_parts(*this, entity, part_storage.induced_part_ordinals, part_storage.removeParts);
      internal_change_entity_parts(entity, part_storage.induced_part_ordinals, part_storage.removeParts, scratch, scratch2);
    }

    std::vector<EntityProc> send_list;
    fill_entity_procs_for_owned_modified_or_created(send_list);
    stk::util::insert_keep_sorted_and_unique(auraEntitiesToResend,
             send_list, EntityLess(*this));
    internal_send_part_memberships_from_owner(send_list);
}

bool is_less_than_element_rank(const BulkData& bulkData, Entity entity)
{
    return bulkData.entity_rank(entity) <= bulkData.mesh_meta_data().side_rank();
}

void BulkData::internal_resolve_shared_part_membership_for_element_death()
{
    STK_ThrowRequireMsg(parallel_size() > 1, "Do not call this in serial");

    ParallelMachine p_comm = parallel();
    const int p_rank = parallel_rank();

    stk::CommSparse comm(p_comm);

    const EntityCommDatabase& commDB = internal_comm_db();
    const EntityCommListInfoVector& entityCommList = commDB.comm_list();
    pack_and_communicate(comm, [this, &comm, &entityCommList]() {
        impl::pack_induced_memberships_for_entities_less_than_element_rank(*this, comm, entityCommList);
    });

    OrdinalVector induced_part_ordinals;
    PartStorage part_storage;

    for(const EntityCommListInfo& info : entityCommList) {
      stk::mesh::Entity entity = info.entity;

      if(is_less_than_element_rank(*this, entity) && is_modified_or_created(*this, entity))
      {
        bool i_own_this_entity_in_comm_list = parallel_owner_rank(entity) == p_rank;
        if( i_own_this_entity_in_comm_list )
        {
          remove_unneeded_induced_parts(entity, commDB.comm(info.entity_comm), part_storage, comm);
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
    const PartVector & all_parts = mesh_meta_data().get_parts();

    stk::CommSparse comm(p_comm);

    pack_and_communicate(comm, [this, &comm, &send_list]() {
      impl::pack_part_memberships(*this, comm, send_list);
    });

    OrdinalVector owner_parts, remove_parts, parts_removed;
    OrdinalVector scratchOrdinalVec, scratchSpace;

    const MetaData & meta = mesh_meta_data();
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

            Entity const entity = impl::find_entity(*this, m_entity_comm_map.comm_list(), key).entity;

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

            internal_change_entity_parts_without_propagating_to_downward_connected_entities_with_notification(entity, owner_parts, remove_parts, parts_removed, scratchOrdinalVec, scratchSpace);
        }
    }
}

void BulkData::internal_resolve_send_ghost_membership()
{
    // This virtual method can be removed when we no longer need the
    // StkTransitionBulkData derived class in Framework.
}

const std::vector<int>& BulkData::all_sharing_procs(stk::mesh::EntityRank rank) const
{
  confirm_host_mesh_is_synchronized_from_device();
  internal_update_all_sharing_procs();
  return m_all_sharing_procs[rank];
}

void BulkData::internal_update_all_sharing_procs() const
{
  if (m_all_sharing_procs_sync_count < synchronized_count()) {
    const EntityRank numRanks = static_cast<EntityRank>(mesh_meta_data().entity_rank_count());
    m_all_sharing_procs.resize(numRanks);
    if (parallel_size() > 1) {
      for (EntityRank r = stk::topology::BEGIN_RANK; r < numRanks; ++r) {
        m_all_sharing_procs[r].clear();
      }

      const EntityCommDatabase& commDB = internal_comm_db();
      const EntityCommListInfoVector& all_comm = commDB.comm_list();
      for (const EntityCommListInfo& info : all_comm) {
        const Entity entity = info.entity;
        const Bucket& bkt = bucket(entity);

        if (bkt.shared() && info.entity_comm != -1) {
          const EntityRank rank = bkt.entity_rank();

          PairIterEntityComm commInfo = commDB.comm(info.entity_comm);
          const unsigned len = commInfo.size();
          unsigned i = 0;
          while (i < len && commInfo[i].ghost_id == BulkData::SHARED) {
            stk::util::insert_keep_sorted_and_unique(commInfo[i++].proc, m_all_sharing_procs[rank]);
          }
        }
      }
    }
    m_all_sharing_procs_sync_count = synchronized_count();
  }
}

template<typename PARTVECTOR>
void internal_throw_error_if_manipulating_internal_part_memberships(const PARTVECTOR & parts)
{
    for(size_t i=0; i<parts.size(); i++)
    {
        STK_ThrowErrorMsgIf( stk::mesh::is_auto_declared_part(*parts[i]) && !stk::mesh::is_topology_root_part(*parts[i]),
                         "Cannot add or remove entity from internally managed part " << parts[i]->name() );
    }
}

template<typename PARTVECTOR>
void BulkData::change_entity_parts( Entity entity,
    const PARTVECTOR & add_parts ,
    const PARTVECTOR & remove_parts)
{
    const bool stkMeshRunningUnderFramework = add_fmwk_data();
    if(!stkMeshRunningUnderFramework)
    {
        internal_throw_error_if_manipulating_internal_part_memberships(add_parts);
        internal_throw_error_if_manipulating_internal_part_memberships(remove_parts);

        require_entity_owner(entity, parallel_rank());
    }
    OrdinalVector scratchOrdinalVec, scratchSpace;
    internal_verify_and_change_entity_parts(entity, add_parts, remove_parts,
                        scratchOrdinalVec, scratchSpace);
}

template void BulkData::change_entity_parts(Entity, const PartVector&, const PartVector&);
template void BulkData::change_entity_parts(Entity, const ConstPartVector&, const ConstPartVector&);

template<typename PARTVECTOR>
void BulkData::change_entity_parts( const EntityVector& entities,
    const PARTVECTOR & add_parts ,
    const PARTVECTOR & remove_parts)
{
    const bool stkMeshRunningUnderFramework = add_fmwk_data();
    if(!stkMeshRunningUnderFramework) {
      internal_throw_error_if_manipulating_internal_part_memberships(add_parts);
      internal_throw_error_if_manipulating_internal_part_memberships(remove_parts);
      for(Entity entity : entities) {
        require_entity_owner(entity, parallel_rank());
      }
    }
    internal_verify_and_change_entity_parts(entities, add_parts, remove_parts);
}

template void BulkData::change_entity_parts(const EntityVector&, const PartVector&, const PartVector&);
template void BulkData::change_entity_parts(const EntityVector&, const ConstPartVector&, const ConstPartVector&);

void BulkData::batch_change_entity_parts( const stk::mesh::EntityVector& entities,
                          const std::vector<PartVector>& add_parts,
                          const std::vector<PartVector>& remove_parts,
                          ModEndOptimizationFlag opt)
{
    const bool stkMeshRunningUnderFramework = add_fmwk_data();
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

    constexpr bool resetSymGhostInfo = false;
    bool starting_modification = internal_modification_begin("batch_change_entity_parts", resetSymGhostInfo);
    STK_ThrowRequireMsg(starting_modification, "ERROR: BulkData already being modified,\n"
                    <<"BulkData::change_entity_parts(vector-of-entities) can not be called within an outer modification scope.");

    if (opt == ModEndOptimizationFlag::MOD_END_SORT) {
      m_bucket_repository.set_remove_mode_tracking();
    }

    OrdinalVector scratchOrdinalVec, scratchSpace;
    for(size_t i=0; i<entities.size(); ++i) {
      if (!stkMeshRunningUnderFramework) {
        require_entity_owner(entities[i], parallel_rank());
      }
      internal_verify_and_change_entity_parts(entities[i], add_parts[i], remove_parts[i],
                       scratchOrdinalVec, scratchSpace);
    }

    internal_modification_end_for_change_parts(opt);
    if (opt == ModEndOptimizationFlag::MOD_END_SORT) {
      m_bucket_repository.set_remove_mode_fill_and_sort();
    }
}

void BulkData::batch_change_entity_parts(const stk::mesh::EntityVector& entities,
                                         const PartVector& add_parts,
                                         const PartVector& remove_parts,
                                         ModEndOptimizationFlag opt)
{
    const bool stkMeshRunningUnderFramework = add_fmwk_data();
    if(!stkMeshRunningUnderFramework) {
      internal_throw_error_if_manipulating_internal_part_memberships(add_parts);
      internal_throw_error_if_manipulating_internal_part_memberships(remove_parts);
      for(Entity entity : entities) {
        require_entity_owner(entity, parallel_rank());
      }
    }

    constexpr bool resetSymGhostInfo = false;
    bool starting_modification = internal_modification_begin("batch_change_entity_parts", resetSymGhostInfo);
    STK_ThrowRequireMsg(starting_modification, "ERROR: BulkData already being modified,\n"
                    <<"BulkData::change_entity_parts(vector-of-entities) can not be called within an outer modification scope.");

    if (opt == ModEndOptimizationFlag::MOD_END_SORT) {
      m_bucket_repository.set_remove_mode_tracking();
    }

    internal_verify_and_change_entity_parts(entities, add_parts, remove_parts);

    internal_modification_end_for_change_parts(opt);
    if (opt == ModEndOptimizationFlag::MOD_END_SORT) {
      m_bucket_repository.set_remove_mode_fill_and_sort();
    }
}

void BulkData::change_entity_parts(const Selector& selector,
                                   EntityRank rank,
                                   const PartVector& add_parts,
                                   const PartVector& remove_parts)
{
    if(m_runConsistencyCheck) {
      impl::check_matching_selectors_and_parts_across_procs(selector, add_parts, remove_parts, parallel());
    }

    bool stkMeshRunningUnderFramework = add_fmwk_data();
    if(!stkMeshRunningUnderFramework)
    {
        internal_throw_error_if_manipulating_internal_part_memberships(add_parts);
        internal_throw_error_if_manipulating_internal_part_memberships(remove_parts);
    }

    internal_verify_and_change_entity_parts(selector, rank, add_parts, remove_parts);
}

void BulkData::batch_change_entity_parts(const Selector& selector,
                                         EntityRank rank,
                                         const PartVector& add_parts,
                                         const PartVector& remove_parts,
                                         ModEndOptimizationFlag opt)
{
    constexpr bool resetSymGhostInfo = false;
    bool starting_modification = internal_modification_begin("batch_change_entity_parts", resetSymGhostInfo);
    STK_ThrowRequireMsg(starting_modification, "ERROR: BulkData already being modified,\n"
                    <<"BulkData::change_entity_parts(vector-of-entities) can not be called within an outer modification scope.");

    change_entity_parts(selector, rank, add_parts, remove_parts);

    internal_modification_end_for_change_parts(opt);
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

void BulkData::internal_change_entity_parts_without_propagating_to_downward_connected_entities(Entity entity, const OrdinalVector& add_parts, const OrdinalVector& remove_parts,
                                                                                               OrdinalVector& ranked_parts_removed, OrdinalVector& newBucketPartList, OrdinalVector& scratchSpace)
{
    internal_adjust_closure_count(entity, add_parts, remove_parts);

    internal_fill_new_part_list_and_removed_part_list(bucket_ptr(entity), add_parts, remove_parts, newBucketPartList, ranked_parts_removed);
    internal_move_entity_to_new_bucket(entity, newBucketPartList, scratchSpace);
}

void BulkData::internal_change_bucket_parts_without_propagating_to_downward_connected_entities(Bucket* bucket, EntityRank rank, const OrdinalVector& add_parts, const OrdinalVector& remove_parts, OrdinalVector& ranked_parts_removed, OrdinalVector& newBucketPartList)
{
  impl::Partition* originalPartition = bucket->getPartition();

  internal_adjust_closure_count(bucket, add_parts, remove_parts);
  internal_fill_new_part_list_and_removed_part_list(bucket, add_parts, remove_parts, newBucketPartList, ranked_parts_removed);

  impl::Partition* partition = m_bucket_repository.get_partition(rank, newBucketPartList);

  if(partition == nullptr) {
    bucket->reset_bucket_parts(newBucketPartList);
    originalPartition->reset_partition_key(bucket->key_vector());
  } else {
    if(originalPartition->get_legacy_partition_id() < partition->get_legacy_partition_id() ||
       partition->get_legacy_partition_id() < originalPartition->get_legacy_partition_id()) {

      originalPartition->remove_bucket(bucket);
      bucket->reset_bucket_parts(newBucketPartList);
      partition->add_bucket(bucket);
    } else {
      bucket->reset_bucket_parts(newBucketPartList);
    }
  }
  m_bucket_repository.set_need_sync_from_partitions(rank);
}

void BulkData::internal_change_entity_parts_without_propagating_to_downward_connected_entities_with_notification(Entity entity, const OrdinalVector& add_parts, const OrdinalVector& remove_parts,
                                                                                                                 OrdinalVector& ranked_parts_removed, OrdinalVector& newBucketPartList, OrdinalVector& scratchSpace)
{
  if (!remove_parts.empty()) {
      notifier.notify_entity_parts_removed(entity, remove_parts);
  }

  internal_change_entity_parts_without_propagating_to_downward_connected_entities(entity, add_parts, remove_parts, ranked_parts_removed, newBucketPartList, scratchSpace);

  if (!add_parts.empty()) {
      notifier.notify_entity_parts_added(entity, add_parts);
  }
}

void BulkData::internal_determine_inducible_parts(EntityRank e_rank, const OrdinalVector& add_parts, const OrdinalVector& parts_removed, OrdinalVector& inducible_parts_added, OrdinalVector& inducible_parts_removed)
{
    impl::fill_inducible_parts_from_list(mesh_meta_data(), add_parts, e_rank, inducible_parts_added);
    impl::fill_inducible_parts_from_list(mesh_meta_data(), parts_removed, e_rank, inducible_parts_removed);
}

void BulkData::internal_determine_inducible_parts_and_propagate_to_downward_connected_entities(Entity entity, const OrdinalVector& add_parts, const OrdinalVector& parts_removed, OrdinalVector& scratchOrdinalVec, OrdinalVector& scratchSpace)
{
  if (impl::are_any_parts_ranked(mesh_meta_data(), add_parts) || !parts_removed.empty()) {
    OrdinalVector inducible_parts_added, inducible_parts_removed;
    internal_determine_inducible_parts(entity_rank(entity), add_parts, parts_removed, inducible_parts_added, inducible_parts_removed);
    internal_propagate_induced_part_changes_to_downward_connected_entities(entity, inducible_parts_added, inducible_parts_removed, scratchOrdinalVec, scratchSpace);
  }
}

void BulkData::internal_determine_inducible_parts_and_propagate_to_downward_connected_entities(Bucket* bucket, const OrdinalVector& add_parts, const OrdinalVector& parts_removed, OrdinalVector& scratchOrdinalVec, OrdinalVector& scratchSpace)
{
  if (impl::are_any_parts_ranked(mesh_meta_data(), add_parts) || !parts_removed.empty()) {
    OrdinalVector inducible_parts_added, inducible_parts_removed;
    internal_determine_inducible_parts(bucket->entity_rank(), add_parts, parts_removed, inducible_parts_added, inducible_parts_removed);
    internal_propagate_induced_part_changes_to_downward_connected_entities(bucket, inducible_parts_added, inducible_parts_removed, scratchOrdinalVec, scratchSpace);
  }
}

bool need_to_change_parts(const Bucket* bkt,
                          const OrdinalVector& addParts,
                          const OrdinalVector& removeParts)
{
  return bkt == nullptr
      || !bkt->member_all(addParts)
      || bkt->member_any(removeParts);
}

void BulkData::internal_change_entity_parts(
  Entity entity ,
  const OrdinalVector& add_parts ,
  const OrdinalVector& remove_parts,
  OrdinalVector& scratchOrdinalVec,
  OrdinalVector& scratchSpace)
{
  require_ok_to_modify();

  Bucket * const bucket_old = bucket_ptr(entity);
  if(need_to_change_parts(bucket_old, add_parts, remove_parts)) {
    if (!remove_parts.empty()) {
        notifier.notify_entity_parts_removed(entity, remove_parts);
    }
    OrdinalVector parts_removed;
    internal_change_entity_parts_without_propagating_to_downward_connected_entities(entity, add_parts, remove_parts, parts_removed, scratchOrdinalVec, scratchSpace);
    internal_determine_inducible_parts_and_propagate_to_downward_connected_entities(entity, add_parts, parts_removed, scratchOrdinalVec, scratchSpace);
    if (!add_parts.empty()) {
      notifier.notify_entity_parts_added(entity, add_parts);
    }
  }
}

void BulkData::internal_change_entity_parts(
  const Selector& selector,
  const EntityRank rank,
  const OrdinalVector& add_parts,
  const OrdinalVector& remove_parts,
  OrdinalVector& scratchOrdinalVec,
  OrdinalVector& scratchSpace)
{
  require_ok_to_modify();

  const BucketVector& buckets = get_buckets(rank, selector);

  for(auto bucket : buckets) {
    bool needToChangeParts = !bucket->member_all(add_parts) || bucket->member_any(remove_parts);

    if(needToChangeParts) {
      notifier.notify_local_buckets_changed(rank);
      Bucket* modifiableBucket = const_cast<Bucket*>(bucket);
      OrdinalVector parts_removed;

      internal_change_bucket_parts_without_propagating_to_downward_connected_entities(modifiableBucket, rank, add_parts, remove_parts, parts_removed, scratchOrdinalVec);
      internal_determine_inducible_parts_and_propagate_to_downward_connected_entities(modifiableBucket, add_parts, parts_removed, scratchOrdinalVec, scratchSpace);
    }
  }

  for(SelectorBucketMap& selectorBucketMap : m_selector_to_buckets_maps) {
    selectorBucketMap.clear();
  }
}

void BulkData::internal_move_entity_to_new_bucket(stk::mesh::Entity entity, const OrdinalVector &newBucketPartList, OrdinalVector& /*scratchSpace*/)
{
    const MeshIndex &meshIndex = mesh_index(entity);
    Bucket *bucketOld = meshIndex.bucket;
    if (bucketOld != nullptr)
    {
        if (!m_meshModification.did_any_shared_entity_change_parts() && parallel_size() > 1) {
          const bool isEntityCommunicated = (bucketOld->shared() || in_send_ghost(entity) || in_receive_ghost(entity));
          if (isEntityCommunicated) {
            m_meshModification.set_shared_entity_changed_parts();
          }
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

void BulkData::internal_fill_new_part_list_and_removed_part_list(stk::mesh::Bucket* bucket_old,
                                                                 const OrdinalVector & add_parts,
                                                                 const OrdinalVector & remove_parts,
                                                                 OrdinalVector &newBucketPartList,
                                                                 OrdinalVector &ranked_parts_removed)
{
    newBucketPartList.clear();
    ranked_parts_removed.clear();

    if(bucket_old != nullptr)
    {
        const std::pair<const unsigned *, const unsigned*> oldEntityPartMembership =
                bucket_old->superset_part_ordinals();

        const size_t num_bucket_parts = bucket_old->supersets().size();
        newBucketPartList.reserve(num_bucket_parts + add_parts.size());
        newBucketPartList.assign(oldEntityPartMembership.first, oldEntityPartMembership.second);

        if(!remove_parts.empty())
        {
            const bool trackRemovedParts = impl::are_any_parts_ranked(mesh_meta_data(), remove_parts);
            ranked_parts_removed.reserve(remove_parts.size());
            impl::filter_out(newBucketPartList, remove_parts, ranked_parts_removed, trackRemovedParts);
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
      newBucketPartList.push_back( mesh_meta_data().universal_part().mesh_meta_data_ordinal() );
    }
    stk::util::sort_and_unique(newBucketPartList);
}

void BulkData::internal_adjust_closure_count(Entity entity,
                                             const OrdinalVector& add_parts,
                                             const OrdinalVector& remove_parts)
{
    Bucket *bucket = bucket_ptr( entity );

    const unsigned locally_owned_ordinal = mesh_meta_data().locally_owned_part().mesh_meta_data_ordinal();

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

void BulkData::internal_adjust_closure_count(Bucket* bucket,
                                             const OrdinalVector& add_parts,
                                             const OrdinalVector& remove_parts)
{
  const unsigned locally_owned_ordinal = mesh_meta_data().locally_owned_part().mesh_meta_data_ordinal();

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
    for(Entity entity : *bucket) {
      unprotect_orphaned_node(entity);
      int incrementClosureCount = 1;
      internal_adjust_entity_and_downward_connectivity_closure_count(entity, bucket, incrementClosureCount);
    }
  }
  else if(remove_from_locally_owned_part)
  {
    for(Entity entity : *bucket) {
      int decrementClosureCount = -1;
      internal_adjust_entity_and_downward_connectivity_closure_count(entity, bucket, decrementClosureCount);
    }
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
  const OrdinalVector & parts_to_remove_assuming_not_induced_from_other_entities,
  OrdinalVector & scratchOrdinalVec,
  OrdinalVector & scratchSpace )
{
    m_check_invalid_rels = false;

    const EntityRank erank = entity_rank(entity);

    OrdinalVector parts_to_actually_remove, emptyParts;
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

                if(!parts_to_remove_assuming_not_induced_from_other_entities.empty())
                {
                    scratchSpace.clear();
                    internal_insert_all_parts_induced_from_higher_rank_entities_to_vector(entity,
                                                                                          e_to,
                                                                                          scratchSpace);
                    internal_fill_parts_to_actually_remove(parts_to_remove_assuming_not_induced_from_other_entities,
                                                           scratchSpace,
                                                           parts_to_actually_remove);
                }

                if (!parts_to_actually_remove.empty())
                {
                    const bool remote_changes_needed = !( parallel_size() == 1 || !bucket(e_to).shared() );
                    if (remote_changes_needed)
                    {
                        Bucket *bucket_old = bucket_ptr(e_to);
                        if (bucket_old && (bucket_old->shared()  || this->in_send_ghost(entity) || this->in_receive_ghost(entity) ))
                        {
                            m_meshModification.set_shared_entity_changed_parts();
                            // Don't remove parts until modification_end to avoid losing field data with bucket churn.
                            parts_to_actually_remove.clear();
                            mark_entity_and_upward_related_entities_as_modified(e_to);
                        }
                    }
                }
                m_modSummary.track_induced_parts(entity, e_to, addParts, parts_to_actually_remove);
                internal_change_entity_parts( e_to , addParts , parts_to_actually_remove, scratchOrdinalVec, scratchSpace );
            }
        }
    }
    m_check_invalid_rels = true;
}

void BulkData::internal_propagate_induced_part_changes_to_downward_connected_entities(
  Bucket* bucket,
  const OrdinalVector & addParts,
  const OrdinalVector & removeParts,
  OrdinalVector& scratchOrdinalPartsRemoved,
  OrdinalVector& scratchOrdinalVec)
{
  EntityRank rank = bucket->entity_rank();

  OrdinalVector partsRemoved;
  EntityVector subs;
  for(EntityRank subRank = stk::topology::BEGIN_RANK; subRank < rank; ++subRank) {
    Bucket* prevBkt = nullptr;
    for(unsigned i=0; i<bucket->size(); ++i) {
      subs.assign(bucket->begin(i, subRank), bucket->end(i, subRank));
      for(stk::mesh::Entity sub : subs) {
        Bucket* subBkt = bucket_ptr(sub);
        if(subBkt == prevBkt || need_to_change_parts(subBkt, addParts, removeParts)) {
          prevBkt = subBkt;
          internal_change_entity_parts_without_propagating_to_downward_connected_entities(sub, addParts, removeParts, partsRemoved, scratchOrdinalPartsRemoved, scratchOrdinalVec);
        }
      }
    }
  }
}

void BulkData::internal_insert_all_parts_induced_from_higher_rank_entities_to_vector(stk::mesh::Entity entity,
                                                                                     stk::mesh::Entity e_to,
                                                                                     OrdinalVector &to_add)
{
    EntityRank e_to_rank = entity_rank(e_to);
    EntityRank start_rank = static_cast<EntityRank>(e_to_rank + 1);
    EntityRank end_rank = static_cast<EntityRank>(mesh_meta_data().entity_rank_count());
    for(EntityRank to_rel_rank_i = start_rank; to_rel_rank_i < end_rank; ++to_rel_rank_i)
    {
        int num_upward_rels = num_connectivity(e_to, to_rel_rank_i);
        Entity const* upward_rel_entities = begin(e_to, to_rel_rank_i);

        const Bucket* prevBucketPtr = nullptr;
        for (int k = 0; k < num_upward_rels; ++k)
        {
            if (entity != upward_rel_entities[k])  // Already did this entity
            {
              STK_ThrowAssertMsg(is_valid(upward_rel_entities[k]), "invalid entity "<<entity_key(upward_rel_entities[k])<<", connected to "<<entity_key(e_to));
              const Bucket* curBucketPtr = bucket_ptr(upward_rel_entities[k]);
              STK_ThrowAssertMsg(curBucketPtr != nullptr, "found null bucket for entity "<<entity_key(upward_rel_entities[k]));
              if (prevBucketPtr != curBucketPtr) {
                prevBucketPtr = curBucketPtr;
                impl::get_part_ordinals_to_induce_on_lower_ranks(*this, *curBucketPtr, e_to_rank, to_add );
              }
            }
        }
    }
}

void BulkData::sortNodesIfNeeded(std::vector<stk::mesh::EntityKey>& nodes)
{
    std::sort(nodes.begin(),nodes.end());
}

inline bool determine_if_is_side(const BulkData& bulk, Entity entity, EntityRank entityRank)
{
  if(entityRank == bulk.mesh_meta_data().side_rank()) {
    return true;
  }

  if(stk::topology::EDGE_RANK == entityRank) {
    // Check if edge is connected to a shell
    auto elements = bulk.begin_elements(entity);
    unsigned numElements = bulk.num_elements(entity);

    for(unsigned i=0; i<numElements; i++) {
      if(bulk.bucket(elements[i]).topology().is_shell()) {
        return true;
      }
    }
  }

  return false;
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
            if (!is_valid(entity)) {
              continue;
            }
            const unsigned num_nodes_on_entity = bucket.num_nodes(entityIndex);

            if (!add_node_sharing_called && this->state(entity) == stk::mesh::Unchanged)
            {
              // No nodes newly shared and entity has not had nodes added, so entity cannot become shared.
              continue;
            }

            const bool isSide = determine_if_is_side(*this, entity, entityRank);
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

void BulkData::remove_entities_from_sharing(const EntityProcVec& entitiesToRemoveFromSharing, EntityVector & entitiesNoLongerShared)
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
                for ( PairIterEntityComm ec = internal_entity_comm_map(side); ! ec.empty() ; ++ec )
                {
                    procs.push_back( ec->proc );
                }
                stk::util::sort_and_unique(procs);

                for(size_t proc_index = 0; proc_index < procs.size(); ++proc_index)
                {
                    const int proc = procs[proc_index];
                    stk::CommBuffer & buf = comm.send_buffer(proc);
                    buf.pack<stk::mesh::EntityKey>(key);
                }

                if(phase == 1)
                {
                    this->entity_comm_map_clear(key);
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
    stk::mesh::OrdinalVector shared_part, owned_part, owned_and_shared, empty;
    shared_part.push_back(mesh_meta_data().globally_shared_part().mesh_meta_data_ordinal());
    owned_part.push_back(mesh_meta_data().locally_owned_part().mesh_meta_data_ordinal());
    owned_and_shared.push_back(mesh_meta_data().locally_owned_part().mesh_meta_data_ordinal());
    owned_and_shared.push_back(mesh_meta_data().globally_shared_part().mesh_meta_data_ordinal());

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
            if (owning_proc == parallel_rank()) {
              internal_change_entity_parts(entity, owned_and_shared /*add*/, empty /*remove*/, scratchOrdinalVec, scratchSpace);
            }
            else {
              internal_change_entity_parts(entity, shared_part /*add*/, owned_part /*remove*/, scratchOrdinalVec, scratchSpace);
            }
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
                                                                 stk::mesh::Part &activePart,
                                                                 ModEndOptimizationFlag opt)
{
    if(!in_synchronized_state())
    {
        STK_ThrowAssertMsg(stk::mesh::impl::check_for_connected_nodes(*this)==0, "BulkData::modification_end ERROR, all entities with rank higher than node are required to have connected nodes.");

        STK_ThrowAssertMsg(add_fmwk_data() || stk::mesh::impl::check_no_shared_elements_or_higher(*this)==0,
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

        internal_finish_modification_end(opt);
    }
}

void BulkData::internal_resolve_sharing_and_ghosting_for_sides(bool connectFacesToPreexistingGhosts)
{
    stk::mesh::EntityVector entitiesNoLongerShared;
    stk::mesh::EntityProcVec entitiesToRemoveFromSharing;
    m_meshModification.delete_shared_entities_which_are_no_longer_in_owned_closure(entitiesToRemoveFromSharing);

    impl::CommEntityMods commEntityMods(*this, internal_comm_db(), internal_comm_list());
    commEntityMods.communicate(impl::CommEntityMods::PACK_SHARED);
    m_meshModification.internal_resolve_shared_modify_delete(commEntityMods.get_shared_mods(), entitiesToRemoveFromSharing, entitiesNoLongerShared);
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
        notifier.notify_started_modification_end();

        STK_ThrowAssertMsg(stk::mesh::impl::check_for_connected_nodes(*this)==0, "BulkData::modification_end ERROR, all entities with rank higher than node are required to have connected nodes.");

        if(parallel_size() > 1)
        {
            set_shared_owned_parts_and_ownership_on_comm_data(sharedModified);

            bool connectFacesToPreexistingGhosts = false;
            internal_resolve_sharing_and_ghosting_for_sides(connectFacesToPreexistingGhosts);
        }

        m_modSummary.write_summary(m_meshModification.synchronized_count(), false);

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
                this->comm_procs(side, commProcs);
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
            STK_ThrowAssertMsg(this->is_valid(side),"Error in communication for de-imprinting the active part on nodes of killed elements in element death!");
            this->internal_change_entity_parts(side, {}, rm_parts, scratchOrdinalVec, scratchSpace);
        }
    );
}

void BulkData::init_mesh_consistency_check_mode()
{
  bool checkByDefault = false;
#ifndef NDEBUG
  checkByDefault = true;
#endif
  m_runConsistencyCheck = get_env_var_as_bool("STK_MESH_RUN_CONSISTENCY_CHECK", checkByDefault);
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

void BulkData::de_induce_parts_from_nodes(const stk::mesh::EntityVector & deactivatedElements, stk::mesh::Part & activePart)
{
    stk::mesh::EntityVector nodesToDeactivate = impl::get_nodes_to_deactivate(*this, deactivatedElements, activePart);
    OrdinalVector scratchOrdinalVec, scratchSpace;
    for (stk::mesh::Entity nodeToDeactivate : nodesToDeactivate)
    {
        this->internal_change_entity_parts(nodeToDeactivate,{}, {activePart.mesh_meta_data_ordinal()}, scratchOrdinalVec, scratchSpace);
    }
}

unsigned BulkData::num_sides(Entity entity) const
{
  if (bucket(entity).topology().has_mixed_rank_sides()) {
    auto num_connected_edges = num_connectivity(entity, stk::topology::EDGE_RANK);
    auto num_connected_faces = num_connectivity(entity, stk::topology::FACE_RANK);

    return num_connected_edges + num_connected_faces;
  } else {
    return num_connectivity(entity, mesh_meta_data().side_rank());
  }
}

bool BulkData::modification_begin(const std::string description)
{
  return internal_modification_begin(description);
}

bool BulkData::modification_begin_for_sync_to_host(const std::string description)
{
  return internal_modification_begin(description, true, true);
}

bool BulkData::internal_modification_begin(const std::string& description, bool resetSymGhostInfo, bool isSyncToHost)
{
  ProfilingBlock block("mod begin:" + description);
  if (!isSyncToHost)
  {
    confirm_host_mesh_is_synchronized_from_device();
    if(m_meshModification.in_modifiable_state()) {
      return false;
    }
  }
  notifier.notify_modification_begin();
  m_lastModificationDescription = description;
  return m_meshModification.modification_begin(description, resetSymGhostInfo, isSyncToHost);
}

bool BulkData::modification_end(ModEndOptimizationFlag modEndOpt)
{
  bool endStatus = false;
  {
      ProfilingBlock block("mod end begin:"+m_lastModificationDescription);
      notifier.notify_started_modification_end();
  }
  {
      ProfilingBlock block("mod end end:"+m_lastModificationDescription);
      endStatus = m_meshModification.modification_end(modEndOpt);
  }
  return endStatus;
}

void BulkData::sort_entities(const stk::mesh::EntitySorterBase& sorter)
{
    STK_ThrowRequireMsg(synchronized_count()>0,"Error, sort_entities must be called after at least one modification cycle.");
    STK_ThrowRequireMsg(in_synchronized_state(), "Error, sort_entities cannot be called from inside a modification cycle.");
    m_bucket_repository.internal_custom_sort_bucket_entities(sorter);

    m_bucket_repository.internal_modification_end();

    if(parallel_size() > 1) {
        check_mesh_consistency();
    }
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

bool BulkData::initialize_face_adjacent_element_graph()
{
    if (m_elemElemGraph == nullptr && m_createUpwardConnectivity)
    {
        m_elemElemGraph = new ElemElemGraph(*this);
        m_elemElemGraphUpdater = std::make_shared<ElemElemGraphUpdater>(*this,*m_elemElemGraph);
        register_observer(m_elemElemGraphUpdater);
        return true;
    }
    return false;
}

void BulkData::delete_face_adjacent_element_graph()
{
    unregister_observer(m_elemElemGraphUpdater);
    delete m_elemElemGraph; m_elemElemGraph = nullptr;
}

stk::mesh::ElemElemGraph& BulkData::get_face_adjacent_element_graph()
{
    STK_ThrowRequireMsg(m_elemElemGraph != nullptr, "Error, Please call initialize_face_adjacent_element_graph before calling get_face_adjacent_element_graph!");
    return *m_elemElemGraph;
}

const stk::mesh::ElemElemGraph& BulkData::get_face_adjacent_element_graph() const
{
    STK_ThrowRequireMsg(m_elemElemGraph != nullptr, "Error, Please call initialize_face_adjacent_element_graph before calling get_face_adjacent_element_graph!");
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
    for(const impl::RelationEntityToNode & relation : relationsToDestroy) {
        destroy_relation(relation.entity, relation.node, relation.ordinal);
    }
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
        record_entity_deletion(e, false);  // the only other user of record_entity_deletion adds the
                                           // entity to the m_deleted_entities_current_modification_cycle if
                                           // it is not a ghost.  Not sure why this doesn't.
    }
}

void
BulkData::internal_check_unpopulated_relations([[maybe_unused]] Entity entity, [[maybe_unused]] EntityRank rank) const
{
#if !defined(NDEBUG) && !defined(__HIP_DEVICE_COMPILE__)
  if (m_check_invalid_rels) {
    const MeshIndex &mesh_idx = mesh_index(entity);
    const Bucket &b = *mesh_idx.bucket;
    const unsigned bucket_ord = mesh_idx.bucket_ordinal;
    STK_ThrowAssertMsg(count_valid_connectivity(entity, rank) == b.num_connectivity(bucket_ord, rank),
                   count_valid_connectivity(entity,rank) << " = count_valid_connectivity("<<entity_key(entity)<<","<<rank<<") != b.num_connectivity("<<bucket_ord<<","<<rank<<") = " << b.num_connectivity(bucket_ord,rank);
                  );

  }
#endif
}

void
BulkData::log_created_parallel_copy(Entity entity)
{
  if (state(entity) == Unchanged) {
    set_state(entity, Modified);
  }
}

bool
BulkData::is_valid_connectivity(Entity entity, EntityRank rank) const
{
  if (!is_valid(entity)) return false;
  if (bucket_ptr(entity) == NULL) return false;
  internal_check_unpopulated_relations(entity, rank);
  return true;
}

void
BulkData::copy_entity_fields(Entity src, Entity dst)
{
  if (src == dst) return;

  //TODO fix const correctness for src
  MeshIndex & src_mesh_idx = mesh_index(src);
  MeshIndex & dst_mesh_idx = mesh_index(dst);

  copy_entity_fields_callback(dst_mesh_idx.bucket->entity_rank(),
                              dst_mesh_idx.bucket->bucket_id(),
                              dst_mesh_idx.bucket_ordinal,
                              src_mesh_idx.bucket->bucket_id(),
                              src_mesh_idx.bucket_ordinal);
}

void BulkData::create_side_entities(const SideSet &sideSet, const stk::mesh::PartVector& parts)
{
    if(has_face_adjacent_element_graph())
        FaceCreator(*this, *m_elemElemGraph).create_side_entities_given_sideset(sideSet, parts);
}

bool BulkData::does_sideset_exist(const stk::mesh::Part &part) const
{
    return m_sideSetData.does_sideset_exist(part);
}

SideSet& BulkData::create_sideset(const stk::mesh::Part &part, bool fromInput)
{
  if(!m_sideSetData.does_sideset_exist(part)) {
    impl::check_sideset_part_constraints(*this, part);
  }

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

EntityRank BulkData::get_entity_rank_count() const
{
  return mesh_meta_data().entity_rank_count();
}

void BulkData::confirm_host_mesh_is_synchronized_from_device(const char * fileName, int lineNumber) const
{
#ifdef STK_USE_DEVICE_MESH
  STK_ThrowRequireMsg((not get_ngp_mesh()) || (not get_ngp_mesh()->need_update_bulk_data()),
                      std::string(fileName) + ":" + std::to_string(lineNumber) +
                      " Accessing host-side BulkData or Field data after a device-side mesh modification without "
                      "calling NgpMesh::need_update_bulk_data()");
#endif
}

namespace impl {

void set_ngp_mesh(const BulkData & bulk, NgpMeshBase * ngpMesh) {
  bulk.set_ngp_mesh(ngpMesh);
  bulk.m_meshModification.set_last_device_synchronized_count(ngpMesh->synchronized_count());
}
}

} // namespace mesh
} // namespace stk
