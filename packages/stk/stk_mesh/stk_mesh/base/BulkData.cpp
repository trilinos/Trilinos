/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.               */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

/**
 * @author H. Carter Edwards
 */

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <stk_util/environment/ReportHandler.hpp>

#include <stk_util/util/StaticAssert.hpp>

#include <stk_util/diag/Trace.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/FindRestriction.hpp>
#include <stk_mesh/baseImpl/EntityRepository.hpp>

#include <boost/foreach.hpp>

namespace stk {
namespace mesh {

namespace impl {
int Counter::counter = 0;

#ifdef STK_MESH_ANALYZE_DYN_CONN
std::vector<DynConnData> DynConnMetrics::m_data;
#endif
}

namespace {

parallel::DistributedIndex::KeySpanVector
convert_entity_keys_to_spans( const MetaData & meta )
{
  // Make sure the distributed index can handle the EntityKey

  enum { OK = StaticAssert<
                SameType< uint64_t,
                          parallel::DistributedIndex::KeyType >::value >::OK };

  // Default constructed EntityKey has all bits set.

  const EntityKey invalid_key ;
  const EntityId  min_id = 1 ;
  const EntityId  max_id = invalid_key.id();

  const size_t rank_count = meta.entity_rank_count();

  parallel::DistributedIndex::KeySpanVector spans( rank_count );

  for ( size_t rank = 0 ; rank < rank_count ; ++rank ) {
    EntityKey key_min( rank , min_id );
    EntityKey key_max( rank , max_id );
    spans[rank].first  = key_min;
    spans[rank].second = key_max;
  }

  return spans ;
}

}

//----------------------------------------------------------------------

BulkData::BulkData( MetaData & mesh_meta_data ,
                    ParallelMachine parallel
#ifdef SIERRA_MIGRATION
                    , bool add_fmwk_data
#endif
                    , ConnectivityMap const* arg_connectivity_map
                    )
  : m_entities_index( parallel, convert_entity_keys_to_spans(mesh_meta_data) ),
    m_entity_repo(*this),
    m_entity_comm_list(),
    m_ghosting(),
    m_deleted_entities(),
    m_deleted_entities_current_modification_cycle(),
    m_mesh_meta_data( mesh_meta_data ),
    m_parallel_machine( parallel ),
    m_parallel_size( parallel_machine_size( parallel ) ),
    m_parallel_rank( parallel_machine_rank( parallel ) ),
    m_sync_count( 0 ),
    m_sync_state( MODIFIABLE ),
    m_meta_data_verified( false ),
    m_mesh_finalized(false),
#ifdef SIERRA_MIGRATION
    m_add_fmwk_data(add_fmwk_data),
    m_fmwk_bulk_ptr(NULL),
    m_check_invalid_rels(true),
#endif
    m_num_fields(-1), // meta data not necessarily committed yet
    m_keep_fields_updated(true),
    m_mesh_indexes(),
    m_entity_keys(),
    m_entity_states(),
    m_closure_count(),
    m_entity_sync_counts(),
    m_local_ids(),
#ifdef SIERRA_MIGRATION
    m_fmwk_aux_relations(),
    m_fmwk_global_ids(),
    m_fmwk_shared_attrs(),
    m_fmwk_connect_counts(),
#endif
    m_field_raw_data(mesh_meta_data.entity_rank_count()),
    m_selector_to_buckets_map(),
#ifdef GATHER_GET_BUCKETS_METRICS
    m_selector_to_count_map(),
    m_num_memoized_get_buckets_calls(0),
    m_num_non_memoized_get_buckets_calls(0),
    m_num_buckets_inserted_in_cache(0),
    m_num_buckets_removed_from_cache(0),
    m_num_modifications(0),
#endif
    m_bucket_repository(
        *this,
        mesh_meta_data.entity_rank_count(),
        arg_connectivity_map != NULL ? *arg_connectivity_map :
           (mesh_meta_data.spatial_dimension() == 2 ? ConnectivityMap::default_map_2d() : ConnectivityMap::default_map())
/*           (mesh_meta_data.spatial_dimension() == 2 ? ConnectivityMap::fixed_edges_map_2d() : ConnectivityMap::fixed_edges_map()) */
                        )
{
  initialize_arrays();

  create_ghosting( "shared" );
  create_ghosting( "shared_aura" );

  m_sync_state = SYNCHRONIZED ;
}

#ifdef GATHER_GET_BUCKETS_METRICS
void BulkData::gather_and_print_get_buckets_metrics()
{
  ParallelMachine world = MPI_COMM_WORLD; // HACK, but necessary to work with Fmwk
  CommAll comm_all(world);

  int real_rank = parallel_machine_rank(world);
  int real_size = parallel_machine_size(world);

  for (int phase = 0; phase < 2; ++phase) {
    CommBuffer& proc_0_buff = comm_all.send_buffer(0); // everything to rank 0

    proc_0_buff.pack<size_t>(m_num_buckets_inserted_in_cache);
    proc_0_buff.pack<size_t>(m_num_buckets_removed_from_cache);

    for (EntityRank r = 0; r < m_mesh_meta_data.entity_rank_count(); ++r) {
      proc_0_buff.pack<size_t>(buckets(r).size());
    }

    proc_0_buff.pack<size_t>(m_num_memoized_get_buckets_calls);
    proc_0_buff.pack<size_t>(m_num_non_memoized_get_buckets_calls);

    proc_0_buff.pack<size_t>(m_selector_to_count_map.size());
    for (SelectorCountMap::const_iterator itr = m_selector_to_count_map.begin(), end = m_selector_to_count_map.end(); itr != end; ++itr) {
      const EntityRank field_array_rank = itr->first.first;
      const size_t cache_hits = itr->second.first;
      const size_t bucket_trav_saved = itr->second.second;

      std::ostringstream out;
      out << itr->first.second;
      const std::string sel_str = out.str();

      proc_0_buff.pack<EntityRank>(field_array_rank);
      proc_0_buff.pack<size_t>(sel_str.size() + 1);
      proc_0_buff.pack<char>(sel_str.c_str(), sel_str.size() + 1);
      proc_0_buff.pack<size_t>(cache_hits);
      proc_0_buff.pack<size_t>(bucket_trav_saved);
    }

    if (phase == 0) { //allocation phase
      comm_all.allocate_buffers( m_parallel_size / 4 );
    }
    else { // communication phase
      comm_all.communicate();
    }
  }

  if (real_rank == 0) {
    std::ostringstream out;
    out << "get_buckets memoization_data:\n";
    out << "  Num modifications: " << m_num_modifications << "\n";

    typedef std::map< std::pair<EntityRank, std::string>, std::pair<size_t, size_t> > GlobalUsageMap;
    GlobalUsageMap global_usage_map;

    size_t global_num_memoized_get_buckets_calls_sum = 0;
    size_t global_num_non_memoized_get_buckets_calls_sum = 0;
    size_t global_num_cache_hits = 0;
    size_t global_num_bucket_trav_saved = 0;
    size_t global_num_buckets_inserted_in_cache = 0;
    size_t global_num_buckets_removed_from_cache = 0;

    std::vector<size_t> bucket_count(m_mesh_meta_data.entity_rank_count(), 0u);

    const size_t MAX_TEXT_LEN = 32768;
    char sel_text[ MAX_TEXT_LEN ];

    for ( int p = 0 ; p < real_size ; ++p ) {
      CommBuffer & buf = comm_all.recv_buffer(p);

      size_t num_buckets_inserted_in_cache = 0;
      size_t num_buckets_removed_from_cache = 0;

      buf.unpack<size_t>(num_buckets_inserted_in_cache);
      buf.unpack<size_t>(num_buckets_removed_from_cache);

      global_num_buckets_inserted_in_cache  += num_buckets_inserted_in_cache;
      global_num_buckets_removed_from_cache += num_buckets_removed_from_cache;

      for (EntityRank r = 0; r < m_mesh_meta_data.entity_rank_count(); ++r) {
        size_t count = 0;
        buf.unpack<size_t>(count);
        bucket_count[r] += count;
      }

      size_t num_memoized_get_buckets_calls = 0;
      size_t num_non_memoized_get_buckets_calls = 0;
      size_t map_size = 0;

      buf.unpack<size_t>(num_memoized_get_buckets_calls);
      buf.unpack<size_t>(num_non_memoized_get_buckets_calls);
      buf.unpack<size_t>(map_size);

      global_num_memoized_get_buckets_calls_sum     += num_memoized_get_buckets_calls;
      global_num_non_memoized_get_buckets_calls_sum += num_non_memoized_get_buckets_calls;

      for (size_t i = 0; i < map_size; ++i) {
        EntityRank field_array_rank = 0;
        size_t str_size = 0;
        size_t cache_hits = 0;
        size_t bucket_trav_saved = 0;

        buf.unpack<EntityRank>(field_array_rank);
        buf.unpack<size_t>(str_size);
        ThrowRequire(str_size < MAX_TEXT_LEN);
        buf.unpack<char>(sel_text, str_size);
        buf.unpack<size_t>(cache_hits);
        buf.unpack<size_t>(bucket_trav_saved);

        global_num_cache_hits += cache_hits;
        global_num_bucket_trav_saved += bucket_trav_saved;

        const std::string sel_str = sel_text;
        std::pair<EntityRank, std::string> search_key = std::make_pair(field_array_rank, sel_str);
        GlobalUsageMap::iterator f_itr = global_usage_map.find(search_key);
        if (f_itr == global_usage_map.end()) {
          global_usage_map[search_key] = std::make_pair(cache_hits, bucket_trav_saved);
        }
        else {
          f_itr->second.first  += cache_hits;
          f_itr->second.second += bucket_trav_saved;
        }
      }
    }

    size_t total_buckets = 0;
    out << "  Bucket counts by rank:\n";
    for (EntityRank r = 0; r < m_mesh_meta_data.entity_rank_count(); ++r) {
      out << "    " << r << ": " << bucket_count[r] << "\n";
      total_buckets += bucket_count[r];
    }
    out << "  Total buckets: " << total_buckets << "\n";
    out << "  Num buckets incrementally inserted in cache: " << global_num_buckets_inserted_in_cache << "\n";
    out << "  Num buckets incrementally removed from cache: " << global_num_buckets_removed_from_cache << "\n";
    out << "  Num memoized get_buckets calls: " << global_num_memoized_get_buckets_calls_sum << "\n";
    out << "  Num non-memoized get_buckets calls: " << global_num_non_memoized_get_buckets_calls_sum << "\n";
    out << "  Num cache hits: " << global_num_cache_hits << "\n";
    out << "  Num bucket traversals saved: " << global_num_bucket_trav_saved << "\n";
    out << "  Total number of selectors used: " << global_usage_map.size() << "\n";
    out << "  Selector usage counts:\n";
    for (GlobalUsageMap::const_iterator itr = global_usage_map.begin(), end = global_usage_map.end(); itr != end; ++itr) {
      out << "  (" << itr->first.first << ", " << itr->first.second << ") -> (" << itr->second.first << ", " << itr->second.second << ")\n";
    }
    std::cout << out.str() << std::endl;
  }
}
#endif

BulkData::~BulkData()
{
#ifdef STK_PROFILE_MEMORY
  print_max_stk_memory_usage(parallel(), parallel_rank(), std::cout);
#endif
#ifdef STK_MESH_ANALYZE_DYN_CONN
  print_dynamic_connectivity_profile(parallel(), parallel_rank(), std::cout);
#endif

#ifdef SIERRA_MIGRATION
  for(size_t i=0; i<m_fmwk_aux_relations.size(); ++i) {
    delete m_fmwk_aux_relations[i];
  }
#endif

#ifdef GATHER_GET_BUCKETS_METRICS
  gather_and_print_get_buckets_metrics();
#endif

  while ( ! m_ghosting.empty() ) {
    delete m_ghosting.back();
    m_ghosting.pop_back();
  }
}

void BulkData::get_selected_nodes(stk::mesh::Selector selector, stk::mesh::EntityVector& nodes)
{
    get_selected_entities( selector, this->buckets(stk::topology::NODE_RANK), nodes );
}

void BulkData::update_deleted_entities_container()
{
  //Question: should the m_deleted_entities container be sorted and uniqued?
  //I.e., should we guard against the same entity being deleted in consecutive modification cycles?

  while(!m_deleted_entities_current_modification_cycle.empty()) {
    size_t entity_offset = m_deleted_entities_current_modification_cycle.front();
    m_deleted_entities_current_modification_cycle.pop_front();
    m_deleted_entities.push_front(entity_offset);
  }
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
  ThrowRequireMsg( m_sync_state != SYNCHRONIZED,
                   "NOT in the ok-to-modify state" );
}

void BulkData::require_entity_owner( const Entity entity ,
                                     int owner ) const
{
  if (parallel_size() > 1 && bucket_ptr(entity) != NULL) {
    const bool error_not_owner = owner != parallel_owner_rank(entity) ;

    ThrowRequireMsg( !error_not_owner,
                     "Entity " << identifier(entity) << " owner is " <<
                     parallel_owner_rank(entity) << ", expected " << owner);
  }
}

void BulkData::require_good_rank_and_id(EntityRank ent_rank, EntityId ent_id) const
{
  const size_t rank_count = m_mesh_meta_data.entity_rank_count();
  const bool ok_id   = EntityKey::is_valid_id(ent_id);
  const bool ok_rank = ent_rank < rank_count && !(ent_rank == MetaData::FACE_RANK && mesh_meta_data().spatial_dimension() == 2);

  ThrowRequireMsg( ok_rank,
                   "Bad key rank: " << ent_rank << " for id " << ent_id );

  ThrowRequireMsg( ok_id, "Bad id : " << ent_id);
}

bool BulkData::is_good_rank_and_id(EntityRank ent_rank, EntityId ent_id) const
{
  const size_t rank_count = m_mesh_meta_data.entity_rank_count();
  const bool ok_id   = EntityKey::is_valid_id(ent_id);
  const bool ok_rank = ent_rank < rank_count && !(ent_rank == MetaData::FACE_RANK && mesh_meta_data().spatial_dimension() == 2);

  return ok_id && ok_rank;
}

void BulkData::require_metadata_committed()
{
  if (!m_mesh_meta_data.is_commit()) {
    m_mesh_meta_data.commit();
  }
}

//----------------------------------------------------------------------

bool BulkData::modification_begin()
{
  Trace_("stk::mesh::BulkData::modification_begin");

  parallel_machine_barrier( m_parallel_machine );

  ThrowRequireMsg( m_mesh_finalized == false, "Unable to modifiy, BulkData has been finalized.");

  if (m_sync_count == 0) {
    m_mesh_meta_data.set_mesh_on_fields(this);
  }

  if ( m_sync_state == MODIFIABLE && m_mesh_finalized == false ) return false ;

  if ( ! m_meta_data_verified ) {
    require_metadata_committed();

    if (parallel_size() > 1) {
      verify_parallel_consistency( m_mesh_meta_data , m_parallel_machine );
    }

    m_meta_data_verified = true ;
  }
  else {
    ++m_sync_count ;

    //Set all entity states to 'Unchanged',
    for ( impl::EntityRepository::const_iterator
            i = m_entity_repo.begin() ; i != m_entity_repo.end() ; ++i )
    {
      m_entity_states[i->second.local_offset()] = Unchanged;
    }
  }

  // // It might be overkill to call this on every modification cycle.

  m_sync_state = MODIFIABLE ;

  return true ;
}

void BulkData::modified(Entity entity)
{
  TraceIfWatching("stk::mesh::BulkData::log_modified_and_propagate", LOG_ENTITY, entity_key(entity));

  // If already in modified state, return
  EntityState entity_state = this->state(entity);
  if (entity_state != Unchanged) {
    return;
  }

  // mark this entity as modified
  this->set_state(entity, Modified);

  // recurse on related entities w/ higher rank
  // outer loop iterates backwards as an optimization to reduce function call depth.
  EntityVector temp_entities;
  Entity const* rels_i = NULL;
  int num_rels = 0;
  EntityRank rank_of_original_entity = entity_rank(entity);
  for (EntityRank irank = m_mesh_meta_data.entity_rank_count() - 1;
        irank > rank_of_original_entity;
        --irank)
  {
    if (connectivity_map().valid(rank_of_original_entity, irank)) {
      num_rels = num_connectivity(entity, irank);
      rels_i   = begin(entity, irank);
    }
    else {
      num_rels = get_connectivity(*this, entity, irank, temp_entities);
      rels_i   = &*temp_entities.begin();
    }

    for (int i = 0; i < num_rels; ++i)
    {
      Entity other_entity = rels_i[i];
      if ( this->state(other_entity) == Unchanged ) {
        this->modified(other_entity);
      }
    }
  }
}

size_t BulkData::count_relations(Entity entity) const
{
  const MeshIndex &mesh_idx = mesh_index(entity);

  const EntityRank end_rank = m_mesh_meta_data.entity_rank_count();
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

  const EntityRank end_rank = m_mesh_meta_data.entity_rank_count();
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

unsigned BulkData::count_valid_connectivity(Entity entity) const
{
  unsigned count = 0;
  const EntityRank end_rank = m_mesh_meta_data.entity_rank_count();
  for (EntityRank irank = stk::topology::BEGIN_RANK; irank < end_rank; ++irank)
  {
    count += count_valid_connectivity(entity, irank);
  }
  return count;
}

size_t BulkData::generate_next_local_offset()
{
  size_t new_local_offset = m_mesh_indexes.size();

  if (!m_deleted_entities.empty()) {
    new_local_offset = m_deleted_entities.front();
    m_deleted_entities.pop_front();
  }

  MeshIndex mesh_index = {NULL, 0};
  EntityKey invalid_key;

  if (new_local_offset == m_mesh_indexes.size()) {
    m_mesh_indexes.push_back(mesh_index);
    m_entity_keys.push_back(invalid_key);
    m_entity_states.push_back(Created);
    m_closure_count.push_back(static_cast<uint16_t>(0));
    m_entity_sync_counts.push_back(0);
    m_local_ids.push_back(stk::mesh::GetInvalidLocalId());

#ifdef SIERRA_MIGRATION
    if (m_add_fmwk_data) {
      m_fmwk_aux_relations.push_back(NULL);
      m_fmwk_global_ids.push_back(0);
      m_fmwk_shared_attrs.push_back(NULL);
      m_fmwk_connect_counts.push_back(0);
    }
#endif
  }
  else {
    //re-claiming space from a previously-deleted entity:

    m_mesh_indexes[new_local_offset] = mesh_index;
    m_entity_keys[new_local_offset] = invalid_key;
    m_entity_states[new_local_offset] = Created;
    m_closure_count[new_local_offset] = static_cast<uint16_t>(0);
    m_entity_sync_counts[new_local_offset] = 0;
    m_local_ids[new_local_offset] = stk::mesh::GetInvalidLocalId();

#ifdef SIERRA_MIGRATION
    if (m_add_fmwk_data) {
      //bulk-data allocated aux-relation vector, so delete it here.
      delete m_fmwk_aux_relations[new_local_offset];
      m_fmwk_aux_relations[new_local_offset] = NULL;
      m_fmwk_global_ids[new_local_offset] = 0;
      //don't delete shared-attr, it was allocated by fmwk.
      m_fmwk_shared_attrs[new_local_offset] = NULL;
      m_fmwk_connect_counts[new_local_offset] = 0;
    }
#endif
  }

  return new_local_offset;
}

void BulkData::initialize_arrays()
{
  ThrowRequireMsg((m_mesh_indexes.size() == 0) && (m_entity_keys.size() == 0)
                        && (m_entity_states.size() == 0) && (m_entity_sync_counts.size() == 0),
                   "BulkData::initialize_arrays() called by something other than constructor");

  MeshIndex mesh_index = {NULL, 0};
  m_mesh_indexes.push_back(mesh_index);

  EntityKey invalid_key;
  m_entity_keys.push_back(invalid_key);

  m_entity_states.push_back(Deleted);
  m_closure_count.push_back(static_cast<uint16_t>(0));
  m_entity_sync_counts.push_back(0);
  m_local_ids.push_back(stk::mesh::GetInvalidLocalId());

#ifdef SIERRA_MIGRATION
  if (m_add_fmwk_data) {
    m_fmwk_aux_relations.push_back(NULL);
    m_fmwk_global_ids.push_back(0);
    m_fmwk_shared_attrs.push_back(NULL);
    m_fmwk_connect_counts.push_back(0);
  }
#endif
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// The add_parts must be full ordered and consistent,
// i.e. no bad parts, all supersets included, and
// owner & used parts match the owner value.

//----------------------------------------------------------------------

Entity BulkData::declare_entity( EntityRank ent_rank , EntityId ent_id ,
                                 const PartVector & parts )
{
  m_check_invalid_rels = false;

  require_ok_to_modify();

  require_good_rank_and_id(ent_rank, ent_id);

  EntityKey key( ent_rank , ent_id );
  TraceIfWatching("stk::mesh::BulkData::declare_entity", LOG_ENTITY, key);
  DiagIfWatching(LOG_ENTITY, key, "declaring entity with parts " << parts);

  std::pair< Entity , bool > result = m_entity_repo.internal_create_entity( key );

  Entity declared_entity = result.first;

  if ( !result.second) {
    // An existing entity, the owner must match.
    require_entity_owner( declared_entity , m_parallel_rank );
    DiagIfWatching(LOG_ENTITY, key, "existing entity: " << declared_entity);
  }

  //------------------------------

  Part * const owns = & m_mesh_meta_data.locally_owned_part();

  PartVector rem ;
  PartVector add( parts );
  add.push_back( owns );

  change_entity_parts( declared_entity , add , rem );

  if ( result.second ) {
    this->set_parallel_owner_rank(declared_entity, m_parallel_rank);
    set_synchronized_count(declared_entity, m_sync_count);
    DiagIfWatching(LOG_ENTITY, key, "new entity: " << declared_entity);
  }

  m_check_invalid_rels = true;

  return declared_entity ;
}

Entity BulkData::declare_entity( EntityRank ent_rank , EntityId ent_id)
{
    Part& universal = mesh_meta_data().universal_part();
    return declare_entity(ent_rank, ent_id, universal);
}

void BulkData::change_entity_id( EntityId id, Entity entity)
{
// When stk parallel is used within Fmwk, this assertion is violated
#ifndef SIERRA_MIGRATION
  ThrowAssertMsg(parallel_size() == 1,
                 "change_entity_id only supported in serial");
#endif

  EntityRank e_rank = entity_rank(entity);

  require_ok_to_modify();
  require_good_rank_and_id(e_rank, id);

  EntityKey new_key(e_rank,id);
  EntityKey old_key = entity_key(entity);

  internal_change_entity_key(old_key, new_key, entity);
}

void BulkData::internal_change_entity_key( EntityKey old_key, EntityKey new_key, Entity entity)
{
  m_entity_repo.update_entity_key(new_key, old_key, entity);
  set_entity_key(entity, new_key);
}

//----------------------------------------------------------------------

bool BulkData::destroy_entity( Entity entity )
{
  TraceIfWatching("stk::mesh::BulkData::destroy_entity", LOG_ENTITY, entity_key(entity));
  DiagIfWatching(LOG_ENTITY, entity_key(entity), "entity state: " << entity);

  require_ok_to_modify();

  m_check_invalid_rels = false;

  if (!is_valid(entity)) {
    m_check_invalid_rels = true;
    return false;
  }

  const EntityRank erank = entity_rank(entity);
  const EntityRank end_rank = m_mesh_meta_data.entity_rank_count();
  for (EntityRank irank = erank + 1; irank != end_rank; ++irank) {
    if (num_connectivity(entity, irank) > 0) {
      m_check_invalid_rels = true;
      return false;
    }
  }

  //------------------------------
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
  EntityVector temp_entities;
  std::vector<ConnectivityOrdinal> temp_ordinals;
  Entity const* rel_entities = NULL;
  int num_conn = 0;
  ConnectivityOrdinal const* rel_ordinals;
  for (EntityRank irank = end_rank; irank != stk::topology::BEGIN_RANK; )
  {
    --irank;

    if (connectivity_map().valid(erank, irank)) {
      num_conn     = num_connectivity(entity, irank);
      rel_entities = begin(entity, irank);
      rel_ordinals = begin_ordinals(entity, irank);
    }
    else {
      num_conn     = get_connectivity(*this, entity, irank, temp_entities, temp_ordinals);
      rel_entities = &*temp_entities.begin();
      rel_ordinals = &*temp_ordinals.begin();
    }

    for (int j = num_conn; j > 0; )
    {
      --j;
      if (is_valid(rel_entities[j])) {
        destroy_relation(entity, rel_entities[j], rel_ordinals[j]);
      }
    }
  }

  // We need to save these items and call remove_entity AFTER the call to
  // destroy_later because remove_entity may destroy the bucket
  // which would cause problems in m_entity_repo.destroy_later because it
  // makes references to the entity's original bucket.

  // Need to invalidate Entity handles in comm-list
  EntityCommListInfoVector::iterator lb_itr =
    std::lower_bound(m_entity_comm_list.begin(), m_entity_comm_list.end(), entity_key(entity));
  if (lb_itr != m_entity_comm_list.end() && lb_itr->key == entity_key(entity)) {
    lb_itr->entity = Entity();
  }

  remove_entity_callback(erank, bucket(entity).bucket_id(), bucket_ordinal(entity));

  m_entities_index.register_removed_key( entity_key(entity) );

  bucket(entity).getPartition()->remove(entity);
  m_entity_repo.destroy_entity(entity_key(entity), entity );
  m_entity_states[entity.local_offset()] = Deleted;
  m_closure_count[entity.local_offset()] = static_cast<uint16_t>(0u);
  m_deleted_entities_current_modification_cycle.push_front(entity.local_offset());

  m_check_invalid_rels = true;
  return true ;
}

//----------------------------------------------------------------------

void BulkData::generate_new_entities(const std::vector<size_t>& requests,
                                 std::vector<Entity>& requested_entities)
{
  Trace_("stk::mesh::BulkData::generate_new_entities");

  typedef stk::parallel::DistributedIndex::KeyType       KeyType;
  typedef stk::parallel::DistributedIndex::KeyTypeVector KeyTypeVector;
  typedef std::vector< KeyTypeVector > RequestKeyVector;

  RequestKeyVector requested_key_types;

  m_entities_index.generate_new_keys(requests, requested_key_types);

  //generating 'owned' entities
  Part * const owns = & m_mesh_meta_data.locally_owned_part();

  std::vector<Part*> rem ;
  std::vector<Part*> add;
  add.push_back( owns );

  requested_entities.clear();
  unsigned cnt=0;
  for (RequestKeyVector::const_iterator itr = requested_key_types.begin(); itr != requested_key_types.end(); ++itr) {
    const KeyTypeVector& key_types = *itr;
    for (KeyTypeVector::const_iterator
        kitr = key_types.begin(); kitr != key_types.end(); ++kitr) {
      ++cnt;
    }
  }
  requested_entities.reserve(cnt);

  for (RequestKeyVector::const_iterator itr = requested_key_types.begin(); itr != requested_key_types.end(); ++itr) {
    const KeyTypeVector & key_types = *itr;
    for ( KeyTypeVector::const_iterator
        kitr = key_types.begin(); kitr != key_types.end(); ++kitr) {
      EntityKey key( static_cast<EntityKey::entity_key_t>((*kitr) ));
      require_good_rank_and_id(key.rank(), key.id());
      std::pair<Entity , bool> result = m_entity_repo.internal_create_entity(key);

      //if an entity is declared with the declare_entity function in
      //the same modification cycle as the generate_new_entities
      //function, and it happens to generate a key that was declared
      //previously in the same cycle it is an error
      ThrowErrorMsgIf( ! result.second,
                       "Generated id " << key.id() <<
                       " which was already used in this modification cycle.");

      // A new application-created entity

      Entity new_entity = result.first;

      //add entity to 'owned' part
      change_entity_parts( new_entity , add , rem );
      requested_entities.push_back(new_entity);

      this->set_parallel_owner_rank( new_entity, m_parallel_rank);
      set_synchronized_count( new_entity, m_sync_count);
    }
  }
}

bool BulkData::in_shared(EntityKey key, int proc) const
{
  PairIterEntityComm sharing = entity_comm_sharing(key);
  for ( ; !sharing.empty(); ++sharing ) {
    if ( proc == sharing->proc ) {
      return true ;
    }
  }
  return false ;
}

bool BulkData::in_send_ghost( EntityKey key , int proc ) const
{
  const int owner_rank = entity_comm_owner(key);
  for ( PairIterEntityComm ec = entity_comm(key); ! ec.empty() ; ++ec ) {
    if ( ec->ghost_id != 0 &&
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

  PairIterEntityComm ec = entity_comm(key);
  EntityCommInfoVector::const_iterator i =
    std::lower_bound( ec.begin(), ec.end() , tmp );

  return i != ec.end() && tmp == *i ;
}

void BulkData::comm_procs( EntityKey key, std::vector<int> & procs ) const
{
  procs.clear();
  for ( PairIterEntityComm ec = entity_comm(key); ! ec.empty() ; ++ec ) {
    procs.push_back( ec->proc );
  }
  std::sort( procs.begin() , procs.end() );
  std::vector<int>::iterator
    i = std::unique( procs.begin() , procs.end() );
  procs.erase( i , procs.end() );
}
void BulkData::comm_shared_procs( EntityKey key, std::vector<int> & procs ) const
{
  procs.clear();
  for ( PairIterEntityComm ec = entity_comm_sharing(key); ! ec.empty() ; ++ec ) {
    procs.push_back( ec->proc );
  }
  std::sort( procs.begin() , procs.end() );
  std::vector<int>::iterator
    i = std::unique( procs.begin() , procs.end() );
  procs.erase( i , procs.end() );
}

void BulkData::shared_procs_intersection( std::vector<EntityKey> & keys, std::vector<int> & procs ) const
{

  procs.clear();
  int num = keys.size();
  std::vector<int> procs_tmp;
  for (int i = 0; i < num; ++i)
  {

    comm_shared_procs(keys[i], procs_tmp);

    if (i == 0)
      procs.swap(procs_tmp);
    else
    {
      // subsequent loops keep the intersection
      std::vector<int> result;
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
  procs.clear();
  for ( PairIterEntityComm ec = entity_comm(key); ! ec.empty() ; ++ec ) {
    if ( ec->ghost_id == ghost.ordinal() ) {
      procs.push_back( ec->proc );
    }
  }
}

void BulkData::internal_change_owner_in_comm_data(const EntityKey& key, int new_owner)
{
  const bool changed = m_entity_comm_map.change_owner_rank(key, new_owner);
  if (changed) {
    EntityCommListInfoVector::iterator lb_itr = std::lower_bound(m_entity_comm_list.begin(),
                                                                        m_entity_comm_list.end(),
                                                                        key);
    if (lb_itr != m_entity_comm_list.end() && lb_itr->key == key) {
      lb_itr->owner = new_owner;
    }
  }
}

void BulkData::internal_sync_comm_list_owners()
{
  for (size_t i = 0, e = m_entity_comm_list.size(); i < e; ++i) {
    m_entity_comm_list[i].owner = parallel_owner_rank(m_entity_comm_list[i].entity);
  }
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

  for(EntityRank rank = stk::topology::NODE_RANK; rank < mesh_meta_data().entity_rank_count(); ++rank) {
    const std::vector<Bucket*>& buckets = this->buckets(rank);
    for(size_t i=0; i<buckets.size(); ++i) {
      new_bucket_callback(rank, buckets[i]->supersets(), buckets[i]->capacity(), NULL);
    }
  }
}

void BulkData::new_bucket_callback(EntityRank rank, const PartVector& superset_parts, size_t capacity, Bucket* new_bucket)
{
  // update selector map
  if (new_bucket != NULL) {
    for (SelectorBucketMap::iterator itr = m_selector_to_buckets_map.begin(), end = m_selector_to_buckets_map.end();
         itr != end; ++itr) {
      Selector const& sel = itr->first.second;
      const EntityRank map_rank = itr->first.first;
      if (map_rank == rank && sel(*new_bucket)) {
#ifdef GATHER_GET_BUCKETS_METRICS
        ++m_num_buckets_inserted_in_cache;
#endif
        TrackedBucketVector& cached_buckets = itr->second;
        TrackedBucketVector::iterator lb_itr = std::lower_bound(cached_buckets.begin(), cached_buckets.end(), new_bucket, BucketIdComparator());
        cached_buckets.insert(lb_itr, new_bucket);
      }
    }
  }

  if (!m_keep_fields_updated) {
    return;
  }

  const std::vector< FieldBase * > & field_set = mesh_meta_data().get_fields();

  if (m_num_fields == -1) {
    // hasn't been set yet
    m_num_fields = field_set.size();

    for (int i = 0; i < m_num_fields; ++i) {
      FieldBase& field = * field_set[i];
      field.get_meta_data_for_field().resize(1);
    }

  }

  // Sizing loop
  size_t total_field_data_size = 0;
  for (int i = 0; i < m_num_fields; ++i) {
    FieldMetaData field_meta_data = {0, NULL, NULL};

    const FieldBase  & field = * field_set[i];
    if (field.entity_rank() == rank)
    {
        unsigned num_bytes_per_entity = 0;

        const FieldBase::Restriction & restriction =
          find_and_check_restriction(field, rank, superset_parts);

        if ( restriction.dimension() > 0 ) { // Exists

          const unsigned type_stride = field.data_traits().stride_of ;
          const unsigned field_rank  = field.field_array_rank();

          num_bytes_per_entity = type_stride *
            ( field_rank ? restriction.stride( field_rank - 1 ) : 1 );

          if (num_bytes_per_entity > 0) {
            field_meta_data.m_size   = num_bytes_per_entity;
            field_meta_data.m_stride = &restriction.stride(0); // JGF: why is this a pointer?

            total_field_data_size += num_bytes_per_entity * capacity;
          }
        }
        field_set[i]->get_meta_data_for_field()[0].push_back(field_meta_data);
    }
  }

  // Allocate all field data for this bucket
  if (total_field_data_size > 0) {
    unsigned char* all_data = field_data_allocator().allocate(total_field_data_size);
    m_field_raw_data[rank].push_back(all_data);

    // Set data ptrs in field meta datas
    size_t current_field_offset = 0;
    for ( int i = 0; i < m_num_fields; ++i ) {
      const FieldBase  & field = * field_set[i];
      if (field.entity_rank() == rank)
      {
          FieldMetaData& field_meta_data = const_cast<FieldMetaData&>(field.get_meta_data_for_field()[0].back());

          if (field_meta_data.m_size > 0) {
            field_meta_data.m_data = all_data + current_field_offset;
            current_field_offset += field_meta_data.m_size * capacity;

            // initialize field data
            const unsigned char* init_val = reinterpret_cast<const unsigned char*>(field.get_initial_value());
            if (init_val != NULL) {
              for (size_t j = 0; j < capacity; ++j) {
                std::memcpy( field_meta_data.m_data + j * field_meta_data.m_size, init_val, field_meta_data.m_size );
              }
            }
            else {
              std::memset( field_meta_data.m_data, 0, capacity * field_meta_data.m_size );
            }
          }
      }
    }
  }
  else {
    m_field_raw_data[rank].push_back(NULL);
  }
}

//
//  Copy fields from src to dst entity.  If the field size of either entity is zero, do nothing.  If the field
//  size of of both entities are non-zero, then the sizes must match
//

void BulkData::copy_entity_fields_callback(EntityRank dst_rank, unsigned dst_bucket_id, Bucket::size_type dst_bucket_ord,
                                           EntityRank src_rank, unsigned src_bucket_id, Bucket::size_type src_bucket_ord)
{
  if (!m_keep_fields_updated) {
    return;
  }

  const std::vector< FieldBase * > & field_set = mesh_meta_data().get_fields();
  for (int i = 0; i < m_num_fields; ++i) {
    if (field_set[i]->entity_rank() == src_rank && field_set[i]->entity_rank() == dst_rank)
    {
        const int src_size        = field_set[i]->get_meta_data_for_field()[0][src_bucket_id].m_size;
        if (src_size == 0) {
          continue;
        }


        unsigned char * const src = field_set[i]->get_meta_data_for_field()[0][src_bucket_id].m_data;
        const int dst_size        = field_set[i]->get_meta_data_for_field()[0][dst_bucket_id].m_size;

        if ( dst_size ) {
          unsigned char * const dst = field_set[i]->get_meta_data_for_field()[0][dst_bucket_id].m_data;
          ThrowAssertMsg( dst_size == src_size,
                          "Incompatible field sizes: " << dst_size << " != " << src_size );

          std::memcpy( dst + dst_size * dst_bucket_ord,
                       src + src_size * src_bucket_ord,
                       dst_size );
        }
    }
  }
}

// More specific version, of copy_entity_fields, here assume know entities are of the same rank

void BulkData::copy_entity_fields_callback_same_rank(EntityRank rank,
                                                     unsigned dst_bucket_id, Bucket::size_type dst_bucket_ord,
                                                     unsigned src_bucket_id, Bucket::size_type src_bucket_ord)
{
  if (!m_keep_fields_updated || m_num_fields == 0) {
    return;
  }


  const std::vector< FieldBase * > & field_set = mesh_meta_data().get_fields();

  for (int i = 0; i < m_num_fields; ++i) {
    if (field_set[i]->entity_rank() == rank)
    {
        const FieldMetaDataVector& metaVec = field_set[i]->get_meta_data_for_field()[0];

        const FieldMetaData& srcMeta = metaVec[src_bucket_id];

        const int src_size        = srcMeta.m_size;
        if (src_size == 0) {
          continue;
        }



        const FieldMetaData& dstMeta = metaVec[dst_bucket_id];

        const int dst_size = dstMeta.m_size;

        ThrowAssertMsg( dst_size == src_size || dst_size == 0, "Incompatible field sizes: " << dst_size << " != " << src_size );



        std::memcpy( dstMeta.m_data + dst_size * dst_bucket_ord,
                     srcMeta.m_data + src_size * src_bucket_ord,
                     dst_size );
    }
  }
}

void BulkData::remove_entity_callback(EntityRank rank, unsigned bucket_id, Bucket::size_type bucket_ord)
{
  if (!m_keep_fields_updated) {
    return;
  }


  const std::vector< FieldBase * > & field_set = mesh_meta_data().get_fields();
  for ( int i = 0; i < m_num_fields; ++i) {
    const FieldBase  & field      = *field_set[i];
    if (field.entity_rank() == rank)
    {
        FieldMetaData field_meta_data = field_set[i]->get_meta_data_for_field()[0][bucket_id];
        const int num_bytes_per_entity = field_meta_data.m_size;

        if (num_bytes_per_entity > 0) {
          // reset field data
          const unsigned char* init_val = reinterpret_cast<const unsigned char*>(field.get_initial_value());
          if (init_val != NULL) {
            std::memcpy( field_meta_data.m_data + bucket_ord * num_bytes_per_entity, init_val, num_bytes_per_entity );
          }
          else {
            std::memset( field_meta_data.m_data + bucket_ord * num_bytes_per_entity, 0, num_bytes_per_entity );
          }
        }
    }
  }
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
#ifdef GATHER_GET_BUCKETS_METRICS
        ++m_num_buckets_removed_from_cache;
#endif
        TrackedBucketVector & cached_buckets = itr->second;
        TrackedBucketVector::iterator lb_itr = std::lower_bound(cached_buckets.begin(), cached_buckets.end(), bucket_id, BucketIdComparator());
        ThrowAssertMsg(lb_itr != cached_buckets.end() && (*lb_itr)->bucket_id() == bucket_id,
                       "Error, bucket id " << bucket_id << ":\n " << dying_bucket << "\nWas selected by selector " << sel << " but was not in the cache");
        cached_buckets.erase(lb_itr);
      }
    }
  }

  if (!m_keep_fields_updated) {
    return;
  }
  const std::vector<FieldBase*>  field_set = mesh_meta_data().get_fields();

  if(field_set.size() == 0) return;

  if (m_field_raw_data[rank][bucket_id] != NULL) {
    size_t bytes_to_delete = 0;
    for (unsigned int i = 0; i < field_set.size(); ++i) {
      if(field_set[i] == NULL || field_set[i]->entity_rank() != rank) continue;
      FieldMetaData& field_data = field_set[i]->get_meta_data_for_field()[0][bucket_id];
      if (field_data.m_data != NULL) {
        bytes_to_delete += field_data.m_size * capacity;
        field_data.m_size = 0;
        field_data.m_data = NULL;
      }
    }
    field_data_allocator().deallocate(m_field_raw_data[rank][bucket_id], bytes_to_delete);
    m_field_raw_data[rank][bucket_id] = NULL;
  }
}

void BulkData::update_field_data_states()
{
  const std::vector<FieldBase*> & field_set = mesh_meta_data().get_fields();

    for ( int i = 0 ; i < m_num_fields ; ) {
      const FieldBase & field = * field_set[i];
      const int outer_idx = i;
      const int num_state = field.number_of_states();
      i += num_state ;




      if (num_state > 1) {
        for ( int b = 0, be = field_set[outer_idx]->get_meta_data_for_field()[0].size(); b < be; ++b) {
          if ( field_set[outer_idx]->get_meta_data_for_field()[0][b].m_size > 0 ) {
            unsigned char* data_last = field_set[outer_idx]->get_meta_data_for_field()[0][b].m_data;
            for ( int s = 1; s < num_state; ++s ) {
              std::swap(field_set[outer_idx+s]->get_meta_data_for_field()[0][b].m_data, data_last);
            }
            field_set[outer_idx]->get_meta_data_for_field()[0][b].m_data = data_last;
          }
        }
      }
    }
}

void BulkData::reorder_buckets_callback(EntityRank rank, const std::vector<unsigned>& id_map)
{
  for (SelectorBucketMap::iterator itr = m_selector_to_buckets_map.begin(), end = m_selector_to_buckets_map.end();
       itr != end; ++itr) {
    TrackedBucketVector& cached_buckets = itr->second;
    std::sort(cached_buckets.begin(), cached_buckets.end(), BucketIdComparator());
  }

  if (!m_keep_fields_updated) {
    return;
  }

  std::vector<unsigned char*> field_raw_data(id_map.size());
  for (unsigned m = 0, e = id_map.size(); m < e; ++m) {
    field_raw_data[m] = m_field_raw_data[rank][id_map[m]];
  }
  m_field_raw_data[rank].swap(field_raw_data);

  const std::vector<FieldBase*> & field_set = mesh_meta_data().get_fields();
  for ( int i = 0 ; i < m_num_fields ; ++i ) {
    if (field_set[i]->entity_rank() == rank)
    {
        FieldMetaDataVector new_fields(id_map.size());
        for ( unsigned m = 0, e = id_map.size(); m < e; ++m ) {
          new_fields[m] = field_set[i]->get_meta_data_for_field()[0][id_map[m]];
        }
        new_fields.swap(field_set[i]->get_meta_data_for_field()[0]);
    }
  }
}

void BulkData::dump_all_mesh_info(std::ostream& out) const
{
  // Dump output for metadata first
  m_mesh_meta_data.dump_all_meta_info(out);

  out << "BulkData "
      << " info...\n";

  const FieldVector& all_fields = m_mesh_meta_data.get_fields();

  // Iterate all buckets for all ranks...
  const std::vector<std::string> & rank_names = m_mesh_meta_data.entity_rank_names();
  for (size_t i = 0, e = rank_names.size(); i < e; ++i) {
    EntityRank rank = i;
    out << "  All " << rank_names[i] << " entities:" << std::endl;

    const std::vector<Bucket*>& buckets = this->buckets(rank);
    BOOST_FOREACH(Bucket* bucket, buckets) {
      out << "    Found bucket " << bucket->bucket_id() << " with superset parts: { ";
      PartVector supersets;
      bucket->supersets(supersets);
      BOOST_FOREACH(Part* part, supersets) {
        out << part->name() << " ";
      }
      out << "}" << std::endl;

      for (size_t b_ord = 0, b_end = bucket->size(); b_ord < b_end; ++b_ord) {
        Entity entity = (*bucket)[b_ord];
        out << "      " << print_entity_key(m_mesh_meta_data, entity_key(entity)) << "(offset: " << entity.local_offset() << ")" << std::endl;

        // Print connectivity
        for (size_t r = 0, re = rank_names.size(); r < re; ++r) {
          if (connectivity_map().valid(rank, r)) {
            out << "        Connectivity to " << rank_names[r] << std::endl;
            Entity const* entities = bucket->begin(b_ord, r);
            ConnectivityOrdinal const* ordinals = bucket->begin_ordinals(b_ord, r);
            const int num_conn         = bucket->num_connectivity(b_ord, r);
            for (int c_itr = 0; c_itr < num_conn; ++c_itr) {
              out << "          " << print_entity_key(m_mesh_meta_data, entity_key(entities[c_itr])) << "[" << ordinals[c_itr] << "]" << std::endl;
            }
          }
        }

        // Print field data
        if (m_num_fields > 0) {
          BOOST_FOREACH(FieldBase* field, all_fields) {

            if(field->entity_rank() != bucket->entity_rank()) continue;

            FieldMetaData field_meta_data = field->get_meta_data_for_field()[0][bucket->bucket_id()];

            unsigned data_size = field_meta_data.m_size;
            if (data_size > 0) { // entity has this field?
              void* data = field_meta_data.m_data + field_meta_data.m_size * b_ord;
              out << "        For field: " << *field << ", has data: ";
              field->print_data(out, data, data_size);
              out << std::endl;
            }
          }
        }
      }

    }
  }
}

RelationIterator BulkData::find_aux_relation(Entity entity, const Relation& relation) const
{
  // Extremely hacky: It would be better to set up the < operator for relations so that lower_bound
  // can return the desired iterator, but any sane definition would probably force a change in
  // relation ordering and that's more than I'm willing to take on now.
  //
  // The current semantics for relation-searching is as follows:
  // Ordered based on derived_type, relation_type, and ordinal in descending precedence
  //   If multiple relations have the same derived_type, relation_type, and ordinal, a linear
  //   scan takes place looking for a matching meshobj. If no such meshobj was found, then
  //   we are left with an iterator pointing to the first relation with a different derived_type,
  //   relation_type, or ordinal. To sum up, the result of the search can either be equivalent to
  //   lower_bound OR upper_bound depending upon the state of the relations... YUCK!

  ThrowAssert(!impl::internal_is_handled_generically(relation.getRelationType()));
  const RelationVector& aux_rels = aux_relations(entity);

  for (RelationIterator rel = aux_rels.begin(); rel != aux_rels.end(); ++rel) {
    if (same_specification(*rel, relation) && rel->entity() != relation.entity()) {
      return rel;
    }
  }

  return sierra::Fmwk::INVALID_RELATION_ITR;
}

void BulkData::set_relation_orientation(Entity from, Entity to, ConnectivityOrdinal to_ord, unsigned to_orientation)
{
  const EntityRank from_rank = entity_rank(from);
  const EntityRank to_rank   = entity_rank(to);

  Entity const*              fwd_rels  = begin(from, to_rank);
  ConnectivityOrdinal const* fwd_ords  = begin_ordinals(from, to_rank);
  Permutation *              fwd_perms = const_cast<Permutation*>(begin_permutations(from, to_rank));
  const int                  num_fwd   = num_connectivity(from, to_rank);

  Entity const*              back_rels  = begin(to, from_rank);
  ConnectivityOrdinal const* back_ords  = begin_ordinals(to, from_rank);
  Permutation *              back_perms = const_cast<Permutation*>(begin_permutations(to, from_rank));
  const int                  num_back   = num_connectivity(to,from_rank);

  // Find and change fwd connectivity
  for (int i = 0; i < num_fwd; ++i, ++fwd_rels, ++fwd_ords, ++fwd_perms) {
    // Allow clients to make changes to orientation
    // Orientations do not affect Relation ordering, so this is safe.
    if (*fwd_rels == to && *fwd_ords == to_ord) {
      *fwd_perms = static_cast<Permutation>(to_orientation);
    }
  }

  // Find and change back connectivity
  for (int i = 0; i < num_back; ++i, ++back_rels, ++back_ords, ++back_perms) {
    // Allow clients to make changes to orientation
    // Orientations do not affect Relation ordering, so this is safe.
    if (*back_rels == from && *back_ords == to_ord) {
      *back_perms = static_cast<Permutation>(to_orientation);
    }
  }
}

namespace {

void find_potential_upward_entities( const BulkData & mesh, Entity node, EntityRank rank, EntityVector & entities)
{
  entities.clear();

  // NOTE: it's adequate to just look at one node

  if ( mesh.connectivity_map().valid(stk::topology::NODE_RANK, rank) ) {
    entities.assign(mesh.begin(node, rank), mesh.end(node, rank));
  }
  else {
    // Have to go up to the elements, then back down to the target rank
    EntityVector elements;
    find_potential_upward_entities(mesh, node, stk::topology::ELEMENT_RANK, elements);

    for (unsigned i = 0, e = elements.size(); i < e; ++i) {
      entities.insert(entities.end(), mesh.begin(elements[i], rank), mesh.end(elements[i], rank));
    }

    std::sort(entities.begin(), entities.end());
    EntityVector::iterator new_end = std::unique(entities.begin(), entities.end());
    entities.erase(new_end, entities.end());
  }
}

template <bool GatherOrdinals, bool GatherPermutations>
size_t get_connectivity_impl(const BulkData & mesh,
                             Entity entity,
                             EntityRank to_rank,
                             EntityVector & entity_scratch_storage,
                             std::vector<ConnectivityOrdinal>* ordinal_scratch_storage = NULL,
                             std::vector<Permutation>* permutation_scratch_storage = NULL)
{
  ThrowAssert( !GatherOrdinals || ordinal_scratch_storage != NULL );
  ThrowAssert( !GatherPermutations || permutation_scratch_storage != NULL );

  const EntityRank source_rank = mesh.entity_rank(entity);

  // if the connectivity stored on the mesh, we shouldn't be calling this
  ThrowAssert( !mesh.connectivity_map().valid(source_rank, to_rank) );
  if ( source_rank >= to_rank || (mesh.mesh_meta_data().spatial_dimension() == 2 && to_rank == stk::topology::FACE_RANK) ) {
    // There is no possible connectivity
    return 0;
  }
  // go through the nodes, keeping in mind that we may be a node
  else {
    entity_scratch_storage.clear();
    if (GatherOrdinals)     { ordinal_scratch_storage->clear(); }
    if (GatherPermutations) { permutation_scratch_storage->clear(); }

    unsigned num_nodes;
    Entity const * nodes;
    if (source_rank == stk::topology::NODE_RANK) {
      num_nodes = 1;
      nodes     = &entity;
    }
    else {
      num_nodes = mesh.num_nodes(entity);
      nodes     = mesh.begin_nodes(entity);
    }

    if (num_nodes > 0) {
      EntityVector potential_upward_entities;
      find_potential_upward_entities( mesh, nodes[0], to_rank, potential_upward_entities);

      for (size_t i = 0u, e = potential_upward_entities.size(); i < e; ++i) {
        Entity potential_upward_entity = potential_upward_entities[i];
        Entity const * potential_sources = mesh.begin(potential_upward_entity, source_rank);
        ConnectivityOrdinal const * potential_ordinals = GatherOrdinals ? mesh.begin_ordinals(potential_upward_entity, source_rank) : NULL;
        Permutation const * potential_permutations = GatherPermutations ? mesh.begin_permutations(potential_upward_entity, source_rank) : NULL;
        const unsigned num_sources = mesh.num_connectivity(potential_upward_entity, source_rank);

        for (unsigned is=0u; is < num_sources; ++is) {
          if ( potential_sources[is] == entity) {
            entity_scratch_storage.push_back(potential_upward_entity);
            if (GatherOrdinals) { ordinal_scratch_storage->push_back(potential_ordinals[is]); }
            if (GatherPermutations && potential_permutations != NULL) { permutation_scratch_storage->push_back(potential_permutations[is]); }
          }
        }
      }
    }

    return entity_scratch_storage.size();
  }
}

static const bool GATHER_ORDINALS     = true;
static const bool GATHER_PERMUTATIONS = true;
static const bool DONT_GATHER_ORDINALS     = false;
static const bool DONT_GATHER_PERMUTATIONS = false;


} // namespace

size_t get_connectivity( const BulkData & mesh,
                         Entity entity,
                         EntityRank to_rank,
                         EntityVector & entity_scratch_storage)
{
  return get_connectivity_impl<DONT_GATHER_ORDINALS, DONT_GATHER_PERMUTATIONS>(mesh, entity, to_rank, entity_scratch_storage);
}

size_t get_connectivity( const BulkData & mesh,
                         Entity entity,
                         EntityRank to_rank,
                         EntityVector & entity_scratch_storage,
                         std::vector<ConnectivityOrdinal> & ordinals )
{
  return get_connectivity_impl<GATHER_ORDINALS, DONT_GATHER_PERMUTATIONS>(mesh, entity, to_rank, entity_scratch_storage, &ordinals);
}

size_t get_connectivity( const BulkData & mesh,
                         Entity entity,
                         EntityRank to_rank,
                         EntityVector & entity_scratch_storage,
                         std::vector<Permutation> & permutations )
{
  std::vector<ConnectivityOrdinal>* ignore = NULL;
  return get_connectivity_impl<DONT_GATHER_ORDINALS, GATHER_PERMUTATIONS>(mesh, entity, to_rank, entity_scratch_storage, ignore, &permutations);
}

size_t get_connectivity( const BulkData & mesh,
                         Entity entity,
                         EntityRank to_rank,
                         EntityVector & entity_scratch_storage,
                         std::vector<ConnectivityOrdinal> & ordinals,
                         std::vector<Permutation> & permutations )
{
  return get_connectivity_impl<GATHER_ORDINALS, GATHER_PERMUTATIONS>(mesh, entity, to_rank, entity_scratch_storage, &ordinals, &permutations);
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
#ifdef GATHER_GET_BUCKETS_METRICS
  ++m_num_memoized_get_buckets_calls;
#endif

  std::pair<EntityRank, Selector> search_item = std::make_pair(rank, selector);
  SelectorBucketMap::iterator fitr =
    m_selector_to_buckets_map.find(search_item);
  if (fitr != m_selector_to_buckets_map.end()) {
#ifdef GATHER_GET_BUCKETS_METRICS
    std::pair<size_t, size_t> & data = m_selector_to_count_map[search_item];
    ++data.first;
    data.second += buckets(field_array_rank).size();
#endif
    TrackedBucketVector const& rv = fitr->second;
    return reinterpret_cast<BucketVector const&>(rv);
  }
  else {
    BucketVector const& all_buckets_for_rank = buckets(rank); // lots of potential side effects! Need to happen before map insertion
    std::pair<SelectorBucketMap::iterator, bool> insert_rv =
      m_selector_to_buckets_map.insert(std::make_pair( std::make_pair(rank, selector), TrackedBucketVector() ));
    ThrowAssertMsg(insert_rv.second, "Should not have already been in map");
    TrackedBucketVector& map_buckets = insert_rv.first->second;
    for (size_t i = 0, e = all_buckets_for_rank.size(); i < e; ++i) {
      if (selector(*all_buckets_for_rank[i])) {
        map_buckets.push_back(all_buckets_for_rank[i]);
      }
    }

#ifdef GATHER_GET_BUCKETS_METRICS
    m_selector_to_count_map[search_item] = std::make_pair(0, 0); // (hits, bucket traversals saved)
#endif

    return reinterpret_cast<BucketVector const&>(map_buckets);
  }
}

void BulkData::get_buckets(EntityRank rank, Selector const& selector, BucketVector & output_buckets) const
{
#ifdef GATHER_GET_BUCKETS_METRICS
  ++m_num_non_memoized_get_buckets_calls;
#endif

  output_buckets.clear();

  BucketVector const& all_buckets_for_rank = buckets(rank);
  for (size_t i = 0, e = all_buckets_for_rank.size(); i < e; ++i) {
    if (selector(*all_buckets_for_rank[i])) {
      output_buckets.push_back(all_buckets_for_rank[i]);
    }
  }
}

} // namespace mesh
} // namespace stk
