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
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/parallel/DistributedIndex.hpp"  // for DistributedIndex, etc
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine, etc
#include "stk_util/util/NamedPair.hpp"
#include "stk_util/util/PairIter.hpp"   // for PairIter
#include "stk_util/util/SameType.hpp"   // for SameType, etc
#include "stk_util/util/TrackingAllocator.hpp"  // for tracking_allocator

namespace stk {
namespace mesh {

namespace impl {
int Counter::counter = 0;

#ifdef STK_MESH_ANALYZE_DYN_CONN
std::vector<DynConnData> DynConnMetrics::m_data;
#endif
}


///////////////////////////////////////////// Functions for creating entities

struct shared_entity_type
{
  stk::topology::topology_t topology;
  std::vector<EntityKey>    nodes;
  EntityKey                 local_key;
  EntityKey                 global_key;
  std::vector<int>          sharing_procs;

  friend inline bool operator < (shared_entity_type const& l, shared_entity_type const& r)
  {
    if (l.topology < r.topology)   return true;
    if (l.topology > r.topology)   return false;
    return l.nodes < r.nodes;
  }

  friend inline bool operator == (shared_entity_type const& l, shared_entity_type const& r)
  {

    bool sameTopologyAndNodes =  (l.topology == r.topology) && (l.nodes==r.nodes);
    return sameTopologyAndNodes;
  }
};

void fillEntityCommInfoForEntity(stk::mesh::Ghosting &ghost_id, stk::mesh::BulkData &mesh, std::vector<stk::mesh::EntityKey> nodes, EntityCommInfoVector &sharing_processors)
{
  size_t num_nodes = nodes.size();
  PairIterEntityComm initial_shared = mesh.entity_comm_map(nodes[0],ghost_id);
  sharing_processors.assign(initial_shared.first, initial_shared.second);

  for(size_t i = 1; i < num_nodes; ++i)
  {
    PairIterEntityComm tmp_shared  = mesh.entity_comm_map(nodes[i], ghost_id);
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

void fillSharedEntities(stk::mesh::Ghosting& ghost_id, stk::mesh::BulkData &mesh, std::vector<shared_entity_type> & shared_entity_map, std::vector<std::vector<shared_entity_type> > &shared_entities )
{
    for(std::vector<shared_entity_type>::const_iterator itr = shared_entity_map.begin(),
            end = shared_entity_map.end(); itr != end; ++itr)
    {
        EntityCommInfoVector sharing_processors;
        fillEntityCommInfoForEntity(ghost_id, mesh, itr->nodes, sharing_processors);

        for(EntityCommInfoVector::const_iterator comm_itr = sharing_processors.begin(),
                comm_end = sharing_processors.end(); comm_itr != comm_end; ++comm_itr)
        {
            if(comm_itr->proc != mesh.parallel_rank())
                shared_entities[comm_itr->proc].push_back(*itr);
        }
    }
}

void communicateSharedEntityInfo(stk::mesh::BulkData &mesh, CommAll &comm, std::vector<std::vector<shared_entity_type> > &shared_entities)
{
    for(int allocation_pass = 0; allocation_pass < 2; ++allocation_pass)
    {
        for(int proc = 0, parallel_size = mesh.parallel_size(); proc < parallel_size; ++proc)
        {
            if(proc != mesh.parallel_rank())
            {
                for(size_t e = 0, num_shared = shared_entities[proc].size(); e < num_shared; ++e)
                {
                    shared_entity_type const & sentity = shared_entities[proc][e];
                    size_t num_nodes_on_entity = sentity.nodes.size();
                    comm.send_buffer(proc).pack<stk::topology::topology_t>(sentity.topology);
                    comm.send_buffer(proc).pack<size_t>(num_nodes_on_entity);  // TEMPORARY HACK: DO NOT PUSH
                    for (size_t i = 0; i < num_nodes_on_entity; ++i )
                    {
                        comm.send_buffer(proc).pack(sentity.nodes[i]);
                    }
                    comm.send_buffer(proc).pack(sentity.local_key);
                }
            }
        }

        if(allocation_pass == 0)
        {
            comm.allocate_buffers(mesh.parallel_size() / 4, 0);
        }
    }
    comm.communicate();
}

void unpackEntityInfromFromOtherProcsAndMarkEntitiesAsSharedAndTrackProcessorsThatNeedAlsoHaveEntity(CommAll &comm, stk::mesh::BulkData& mesh, std::vector<shared_entity_type> & shared_entity_map)
{
    for(int ip = mesh.parallel_size() - 1; ip >= 0; --ip)
    {
        if(ip != mesh.parallel_rank())
        {
            CommBuffer & buf = comm.recv_buffer(ip);
            while(buf.remaining())
            {
                shared_entity_type sentity;

                buf.unpack<stk::topology::topology_t>(sentity.topology);
                //stk::topology entity_topology(sentity.topology);
                //size_t num_nodes_on_entity = entity_topology.num_nodes();
                size_t num_nodes_on_entity = 0;            // TEMPORARY HACK:
                buf.unpack<size_t>(num_nodes_on_entity);   // DO NOT PUSH
                sentity.nodes.resize(num_nodes_on_entity);
                for (size_t i = 0; i < num_nodes_on_entity; ++i )
                {
                    buf.unpack<EntityKey>(sentity.nodes[i]);
                }
                buf.unpack<EntityKey>(sentity.global_key);

                std::vector<shared_entity_type>::iterator shared_itr = std::lower_bound(shared_entity_map.begin(), shared_entity_map.end(), sentity);

                //update the global global_key
                if(shared_itr != shared_entity_map.end() && *shared_itr == sentity)
                {
                    Entity entity = mesh.get_entity(shared_itr->local_key);
                    shared_itr->sharing_procs.push_back(ip);
                    if(ip < mesh.parallel_rank())
                    {
                        shared_itr->global_key = sentity.global_key;
                    }
                    mesh.markEntity(entity, BulkData::IS_SHARED);
                }
            }
        }
    }
}

void resolveUniqueIdForSharedEntityAndCreateCommMapInfoForSharingProcs(stk::mesh::BulkData &mesh, std::vector<shared_entity_type> & shared_entity_map)
{
   for(size_t i = 0, e = shared_entity_map.size(); i < e; ++i)
   {
       Entity entity = mesh.get_entity(shared_entity_map[i].local_key);
       if(shared_entity_map[i].global_key != shared_entity_map[i].local_key)
       {
           mesh.internal_change_entity_key(shared_entity_map[i].local_key, shared_entity_map[i].global_key, entity);
       }
       for(size_t j = 0; j < shared_entity_map[i].sharing_procs.size(); j++)
       {
           mesh.entity_comm_map_insert(entity, EntityCommInfo(stk::mesh::BulkData::SHARED, shared_entity_map[i].sharing_procs[j]));
       }
   }
}

void update_shared_entities_global_ids( BulkData & mesh, std::vector<shared_entity_type> & shared_entity_map)
{
    std::sort(shared_entity_map.begin(), shared_entity_map.end());

    // shared_edges[0] will contain all the edges this processor shares with processor 0
    std::vector<std::vector<shared_entity_type> > shared_entities(mesh.parallel_size());
    fillSharedEntities(mesh.shared_ghosting(), mesh, shared_entity_map, shared_entities);

    CommAll comm(mesh.parallel());
    communicateSharedEntityInfo(mesh, comm, shared_entities);
    unpackEntityInfromFromOtherProcsAndMarkEntitiesAsSharedAndTrackProcessorsThatNeedAlsoHaveEntity(comm, mesh, shared_entity_map);
    resolveUniqueIdForSharedEntityAndCreateCommMapInfoForSharingProcs(mesh, shared_entity_map);
}

void markEntitiesForResolvingSharingInfoUsingNodes(stk::mesh::BulkData &mesh, stk::mesh::EntityRank entityRank, std::vector<shared_entity_type>& shared_entities)
{
    const stk::mesh::BucketVector& entity_buckets = mesh.buckets(entityRank);

    for(size_t bucketIndex = 0; bucketIndex < entity_buckets.size(); bucketIndex++)
    {
        const stk::mesh::Bucket& bucket = *entity_buckets[bucketIndex];
        stk::topology topology = bucket.topology();
        for(size_t entityIndex = 0; entityIndex < bucket.size(); entityIndex++)
        {
            Entity entity = bucket[entityIndex];
            const unsigned num_nodes_on_entity = bucket.num_nodes(entityIndex);

            if ( num_nodes_on_entity > 1 )
            {
                if(stk::mesh::in_owned_closure(mesh, entity, mesh.parallel_rank()))
                {
                    EntityVector entity_nodes(num_nodes_on_entity);

                    Entity const * nodes = bucket.begin_nodes(entityIndex);
                    entity_nodes.assign(nodes, nodes+num_nodes_on_entity);

                    //do we need to do some sorting operation here?
                    //sort entity nodes into lexicographical smallest permutation?


                    bool shared_entity = true;
                    for(size_t n = 0; n < num_nodes_on_entity; ++n)
                    {
                        Entity node = entity_nodes[n];
                        shared_entity = shared_entity && (mesh.bucket(node).shared() || (mesh.isEntityMarked(node) == BulkData::IS_SHARED));
                    }

                    if(shared_entity)
                    {
                        std::sort(entity_nodes.begin(),entity_nodes.end(),EntityLess(mesh));
                        shared_entity_type sentity;
                        sentity.topology = topology;
                        sentity.nodes.resize(num_nodes_on_entity);
                        for(size_t n = 0; n < num_nodes_on_entity; ++n)
                        {
                            sentity.nodes[n]=mesh.entity_key(entity_nodes[n]);
                        }
                        const EntityKey &entity_key = mesh.entity_key(entity);
                        sentity.local_key = entity_key;
                        sentity.global_key = entity_key;
                        shared_entities.push_back(sentity);
                        mesh.markEntity(entity, BulkData::POSSIBLY_SHARED);
                    }
                }
            }
        }
    }
}

void gather_shared_nodes(stk::mesh::BulkData & mesh, std::vector<EntityKey> & shared_nodes)
{
    const stk::mesh::BucketVector & node_buckets = mesh.buckets(stk::topology::NODE_RANK);

    for(size_t nodeIndex = 0; nodeIndex < node_buckets.size(); ++nodeIndex)
    {
        const stk::mesh::Bucket & bucket = *node_buckets[nodeIndex];
        for(size_t entityIndex = 0; entityIndex < bucket.size(); ++entityIndex)
        {
            Entity node = bucket[entityIndex];
            if (mesh.isEntityMarked(node) == BulkData::IS_SHARED)
            {
                shared_nodes.push_back(mesh.entity_key(node));
            }
        }
    }
}

void resolve_entity_sharing(stk::mesh::BulkData &stkMeshBulkData, stk::mesh::EntityRank entityRank, std::vector<EntityKey> &entity_keys)
{
    std::vector<shared_entity_type> shared_entities;
    markEntitiesForResolvingSharingInfoUsingNodes(stkMeshBulkData, entityRank, shared_entities);
    update_shared_entities_global_ids( stkMeshBulkData, shared_entities );

    for (size_t i=0; i<shared_entities.size();i++)
    {
        Entity entity = stkMeshBulkData.get_entity(shared_entities[i].global_key);
        if ( stkMeshBulkData.isEntityMarked(entity) == BulkData::IS_SHARED )
        {
            entity_keys.push_back(shared_entities[i].global_key);
        }
        stkMeshBulkData.markEntity(entity, BulkData::NOT_MARKED);
    }
    std::sort(entity_keys.begin(), entity_keys.end());
}

/////////////////////////////////////// End functions for create edges

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

//----------------------------------------------------------------------

#ifdef STK_MESH_MODIFICATION_COUNTERS
unsigned BulkData::m_num_bulk_data_counter = 0;
#endif

BulkData::BulkData( MetaData & mesh_meta_data ,
                    ParallelMachine parallel
#ifdef SIERRA_MIGRATION
                    , bool add_fmwk_data
#endif
                    , ConnectivityMap const* arg_connectivity_map
                    , FieldDataManager *field_data_manager
                    )
  : m_entities_index( parallel, convert_entity_keys_to_spans(mesh_meta_data) ),
    m_entity_repo(*this),
    m_entity_comm_list(),
    m_ghosting(),
    m_ghost_parts(),
    m_deleted_entities(),
    m_deleted_entities_current_modification_cycle(),
    m_ghost_reuse_map(),
    m_mesh_meta_data( mesh_meta_data ),
    m_parallel_machine( parallel ),
    m_parallel_size( parallel_machine_size( parallel ) ),
    m_parallel_rank( parallel_machine_rank( parallel ) ),
    m_sync_count( 0 ),
    m_sync_state( MODIFIABLE ),
    m_meta_data_verified( false ),
    m_mesh_finalized(false),
    m_modification_begin_description("UNSET"),
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
    m_mark_entity(),
    m_closure_count(),
    m_entity_sync_counts(),
    m_local_ids(),
#ifdef SIERRA_MIGRATION
    m_fmwk_aux_relations(),
    m_fmwk_global_ids(),
    m_fmwk_shared_attrs(),
    m_fmwk_connect_counts(),
#endif
    m_default_field_data_manager(mesh_meta_data.entity_rank_count()),
    m_field_data_manager(field_data_manager),
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
#ifdef STK_MESH_MODIFICATION_COUNTERS
,
    m_modification_counters()
#endif
{
#ifdef STK_MESH_MODIFICATION_COUNTERS
  m_num_bulk_data_counter++;
  std::ofstream outfile(create_modification_counts_filename().c_str());
  write_modification_labels_to_stream(outfile);
  outfile.close();
  reset_modification_counters();
#endif

  mesh_meta_data.set_mesh_bulk_data(this);

  if (m_field_data_manager == NULL)
  {
      m_field_data_manager = &m_default_field_data_manager;
  }

  initialize_arrays();

  m_ghost_parts.clear();
  internal_create_ghosting( "shared" );
  //shared part should reside in m_ghost_parts[0]
  internal_create_ghosting( "shared_aura" );

  m_sync_state = SYNCHRONIZED ;
}

void BulkData::reset_modification_counters()
{
#ifdef STK_MESH_MODIFICATION_COUNTERS
    for(unsigned j=0; j<static_cast<unsigned>(NumMethodTypes); ++j)
    {
        for(unsigned i=0; i<static_cast<unsigned>(NumModificationTypes); ++i)
        {
            m_modification_counters[j][i] = 0;
        }
        for(unsigned i=0; i<static_cast<unsigned>(NumEntityModificationTypes); ++i)
        {
            for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<mesh_meta_data().entity_rank_count(); rank++)
            {
                m_entity_modification_counters[j][rank][i] = 0;
            }
        }
    }
#endif
}

std::string BulkData::convert_label_for_method_type(const std::string &label, enum PublicOrInternalMethod methodType)
{
    std::string newLabel = label;
    if(methodType == INTERNAL)
    {
        newLabel = "INTERNAL-" + label;
    }
    return newLabel;
}

void BulkData::write_modification_entry_label(std::ostream& out, const std::string& label, enum PublicOrInternalMethod methodType)
{
    out << convert_label_for_method_type(label, methodType) << ", ";
}

void BulkData::write_entity_modification_entry_label(std::ostream& out, const std::string& label, enum PublicOrInternalMethod methodType)
{
#ifdef STK_MESH_MODIFICATION_COUNTERS
    for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<mesh_meta_data().entity_rank_count(); rank++)
    {
        out << convert_label_for_method_type(label, methodType) <<"["<<rank<<"], ";
    }
    out << convert_label_for_method_type(label, methodType) << ", ";
#endif
}

void BulkData::write_modification_labels_to_stream_for_method_type(std::ostream& out, enum PublicOrInternalMethod methodType)
{
#ifdef STK_MESH_MODIFICATION_COUNTERS
    write_entity_modification_entry_label(out, "declare_entity", methodType);
    write_entity_modification_entry_label(out, "destroy_entity", methodType);
    write_entity_modification_entry_label(out, "change_entity_id", methodType);
    write_entity_modification_entry_label(out, "change_entity_parts", methodType);
    write_modification_entry_label(out, "change_entity_owner", methodType);
    write_modification_entry_label(out, "create_ghosting", methodType);
    write_modification_entry_label(out, "change_ghosting", methodType);
    write_modification_entry_label(out, "destroy_ghosting", methodType);
    write_modification_entry_label(out, "destroy_all_ghosting", methodType);
    write_modification_entry_label(out, "declare_relation", methodType);
    out << convert_label_for_method_type("destroy_relation", methodType);
#endif
}

void BulkData::write_modification_labels_to_stream(std::ostream& out)
{
#ifdef STK_MESH_MODIFICATION_COUNTERS
    out << "modification cycle, ";
    write_modification_labels_to_stream_for_method_type(out, PUBLIC);
    out << ", ";
    write_modification_labels_to_stream_for_method_type(out, INTERNAL);
    out << std::endl;
#endif
}

void BulkData::write_entity_modification_entry(std::ostream& out,
                                               enum PublicOrInternalMethod methodType,
                                               EntityModificationTypes entityModification)
{
#ifdef STK_MESH_MODIFICATION_COUNTERS
    unsigned sum = 0;
    for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<mesh_meta_data().entity_rank_count(); rank++)
    {
        out << m_entity_modification_counters[methodType][rank][entityModification]<<", ";
        sum += m_entity_modification_counters[methodType][rank][entityModification];
    }
    out << sum << ", ";
#endif
}

void BulkData::write_modification_counts_to_stream_for_method_type(std::ostream& out, enum PublicOrInternalMethod methodType)
{
#ifdef STK_MESH_MODIFICATION_COUNTERS
    write_entity_modification_entry(out, methodType, DECLARE_ENTITY);
    write_entity_modification_entry(out, methodType, DESTROY_ENTITY);
    write_entity_modification_entry(out, methodType, CHANGE_ENTITY_ID);
    write_entity_modification_entry(out, methodType, CHANGE_ENTITY_PARTS);
    out << m_modification_counters[methodType][CHANGE_ENTITY_OWNER]<<", ";
    out << m_modification_counters[methodType][CREATE_GHOSTING]<<", ";
    out << m_modification_counters[methodType][CHANGE_GHOSTING]<<", ";
    out << m_modification_counters[methodType][DESTROY_GHOSTING]<<", ";
    out << m_modification_counters[methodType][DESTROY_ALL_GHOSTING]<<", ";
    out << m_modification_counters[methodType][DECLARE_RELATION]<<", ";
    out << m_modification_counters[methodType][DESTROY_RELATION];
#endif
}

void BulkData::write_modification_counts_to_stream(std::ostream& out)
{
#ifdef STK_MESH_MODIFICATION_COUNTERS
    out << synchronized_count()<<", ";
    write_modification_counts_to_stream_for_method_type(out, PUBLIC);
    out << ", ";
    write_modification_counts_to_stream_for_method_type(out, INTERNAL);
    out << std::endl;
#endif
}

std::string BulkData::create_modification_counts_filename() const
{
    std::string fileName;
#ifdef STK_MESH_MODIFICATION_COUNTERS
    std::ostringstream oss;
    int numProcs = parallel_machine_size(MPI_COMM_WORLD);
    int procId = parallel_machine_rank(MPI_COMM_WORLD);
    oss<<"modification_counts_"<<m_num_bulk_data_counter<<"_np"<<numProcs<<"."<<procId<<".csv";
    fileName = oss.str();
#endif
    return fileName;
}

void BulkData::write_modification_counts()
{
#ifdef STK_MESH_MODIFICATION_COUNTERS
    std::ofstream outfile(create_modification_counts_filename().c_str(), std::ios::app);
    write_modification_counts_to_stream(outfile);
    reset_modification_counters();
#endif
}

#ifdef GATHER_GET_BUCKETS_METRICS
void BulkData::gather_and_print_get_buckets_metrics() const
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
      const EntityRank rank = itr->first.first;
      const size_t cache_hits = itr->second.first;
      const size_t bucket_trav_saved = itr->second.second;

      std::ostringstream out;
      out << itr->first.second;
      const std::string sel_str = out.str();

      proc_0_buff.pack<EntityRank>(rank);
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
        EntityRank rank = 0;
        size_t str_size = 0;
        size_t cache_hits = 0;
        size_t bucket_trav_saved = 0;

        buf.unpack<EntityRank>(rank);
        buf.unpack<size_t>(str_size);
        ThrowRequire(str_size < MAX_TEXT_LEN);
        buf.unpack<char>(sel_text, str_size);
        buf.unpack<size_t>(cache_hits);
        buf.unpack<size_t>(bucket_trav_saved);

        global_num_cache_hits += cache_hits;
        global_num_bucket_trav_saved += bucket_trav_saved;

        const std::string sel_str = sel_text;
        std::pair<EntityRank, std::string> search_key = std::make_pair(rank, sel_str);
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

void BulkData::gather_and_print_mesh_partitioning() const
{
  ParallelMachine world = MPI_COMM_WORLD; // HACK, but necessary to work with Fmwk
  CommAll comm_all(world);

  int real_rank = parallel_machine_rank(world);
  int real_size = parallel_machine_size(world);

  for (int phase = 0; phase < 2; ++phase) {
    CommBuffer& proc_0_buff = comm_all.send_buffer(0); // everything to rank 0

    const EntityRank rank_count = static_cast<EntityRank>(m_mesh_meta_data.entity_rank_count());

    for (EntityRank r = stk::topology::NODE_RANK; r < rank_count; ++r) {
      std::vector<impl::Partition*> partitions = m_bucket_repository.get_partitions(r);
      proc_0_buff.pack<size_t>(partitions.size());

      for (size_t p = 0, pe = partitions.size(); p < pe; ++p) {
        const std::vector<PartOrdinal> & parts = partitions[p]->get_legacy_partition_id();
        proc_0_buff.pack<size_t>(partitions[p]->size());

        std::ostringstream out;
        for (size_t i = 1, ie = parts.size() - 1; i < ie; ++i) {
          Part& part = m_mesh_meta_data.get_part(parts[i]);
          out << part.name() << "[" << part.primary_entity_rank() << "] ";
        }

        const std::string partition_str = out.str();

        proc_0_buff.pack<size_t>(partition_str.size() + 1);
        proc_0_buff.pack<char>(partition_str.c_str(), partition_str.size() + 1);
      }
    }

    if (phase == 0) { //allocation phase
      BABBLE_STK_PARALLEL_COMM(world, "          gather_and_print_mesh_partitioning calling allocate_buffers");
      comm_all.allocate_buffers( m_parallel_size / 4 );
    }
    else { // communication phase
      BABBLE_STK_PARALLEL_COMM(world, "          gather_and_pring_mesh_partitioning calling allocate_buffers");
      comm_all.communicate();
    }
  }

  if (real_rank == 0) {
    std::ostringstream out;
    out << "mesh partitioning data:\n";

    typedef std::map< std::pair<EntityRank, std::string>, size_t> GlobalPartitionMap;
    GlobalPartitionMap global_partition_map;

    const size_t MAX_TEXT_LEN = 32768;
    char partition_text[ MAX_TEXT_LEN ];

    for ( int p = 0 ; p < real_size ; ++p ) {
      CommBuffer & buf = comm_all.recv_buffer(p);

      const EntityRank rank_count = static_cast<EntityRank>(m_mesh_meta_data.entity_rank_count());

      for (EntityRank r = stk::topology::NODE_RANK; r < rank_count; ++r) {
        size_t num_partitions_for_this_rank = 0;
        buf.unpack<size_t>(num_partitions_for_this_rank);

        for (size_t partition = 0; partition < num_partitions_for_this_rank; ++partition) {
          size_t num_entities_in_partition = 0;
          buf.unpack<size_t>(num_entities_in_partition);

          size_t str_size = 0;
          buf.unpack<size_t>(str_size);
          ThrowRequire(str_size < MAX_TEXT_LEN);
          buf.unpack<char>(partition_text, str_size);

          const std::string partition_str = partition_text;
          std::pair<EntityRank, std::string> search_key = std::make_pair(r, partition_str);
          GlobalPartitionMap::iterator f_itr = global_partition_map.find(search_key);
          if (f_itr == global_partition_map.end()) {
            global_partition_map[search_key] = num_entities_in_partition;
          }
          else {
            f_itr->second += num_entities_in_partition;
          }
        }
      }
    }

    for (GlobalPartitionMap::const_iterator itr = global_partition_map.begin(), end = global_partition_map.end(); itr != end; ++itr) {
      out << "  (" << itr->first.first << ", " << itr->first.second << ") -> " << itr->second << "\n";
    }
    std::cout << out.str() << std::endl;
  }
}

BulkData::~BulkData()
{

#ifdef STK_PROFILE_MEMORY
  ParallelMachine world = MPI_COMM_WORLD; // HACK, but necessary to work with Fmwk
  const int real_rank = parallel_machine_rank(world);
  print_max_stk_memory_usage(world, real_rank, std::cout);
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

#ifdef PRINT_MESH_PARTITIONING
  gather_and_print_mesh_partitioning();
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

  // Reclaim offsets for deleted ghosted that were not regenerated
  for (GhostReuseMap::iterator m_itr = m_ghost_reuse_map.begin(), m_end = m_ghost_reuse_map.end(); m_itr != m_end; ++m_itr) {
    m_deleted_entities.push_front(m_itr->second);
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

void BulkData::require_metadata_committed()
{
  if (!m_mesh_meta_data.is_commit()) {
    m_mesh_meta_data.commit();
  }
}

//----------------------------------------------------------------------

const std::string& BulkData::get_modification_begin_description() const {
  return m_modification_begin_description;
}

//----------------------------------------------------------------------

bool BulkData::modification_begin(const std::string description)
{
  Trace_("stk::mesh::BulkData::modification_begin");

  parallel_machine_barrier( m_parallel_machine );

  ThrowRequireMsg( m_mesh_finalized == false, "Unable to modifiy, BulkData has been finalized.");

  if (m_sync_count == 0) {
    m_mesh_meta_data.set_mesh_on_fields(this);
  }

  if ( m_sync_state == MODIFIABLE && m_mesh_finalized == false ) return false ;


  m_modification_begin_description = description;


  if ( ! m_meta_data_verified ) {
    require_metadata_committed();

    if (parallel_size() > 1) {
      verify_parallel_consistency( m_mesh_meta_data , m_parallel_machine );
    }

    m_meta_data_verified = true ;
  }
  else {
    ++m_sync_count ;

    for(unsigned i=0, iend=m_entity_states.size(); i<iend; ++i) {
      if(m_entity_states[i] != Deleted) {
        m_entity_states[i] = Unchanged;
      }
    }
  }

  // // It might be overkill to call this on every modification cycle.

  m_sync_state = MODIFIABLE ;

  return true ;
}

void BulkData::mark_entity_and_upward_related_entities_as_modified(Entity entity)
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
  for (EntityRank irank = static_cast<EntityRank>(m_mesh_meta_data.entity_rank_count() - 1);
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
        this->mark_entity_and_upward_related_entities_as_modified(other_entity);
      }
    }
  }
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

unsigned BulkData::count_valid_connectivity(Entity entity) const
{
  unsigned count = 0;
  const EntityRank end_rank = static_cast<EntityRank>(m_mesh_meta_data.entity_rank_count());
  for (EntityRank irank = stk::topology::BEGIN_RANK; irank < end_rank; ++irank)
  {
    count += count_valid_connectivity(entity, irank);
  }
  return count;
}

size_t BulkData::generate_next_local_offset(size_t preferred_offset)
{
  size_t new_local_offset = m_mesh_indexes.size();

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
    m_entity_states.push_back(Created);
    m_mark_entity.push_back(NOT_MARKED);
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
    m_mark_entity[new_local_offset] = NOT_MARKED;
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
  m_mark_entity.push_back(NOT_MARKED);
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

Entity BulkData::declare_entity( EntityRank ent_rank , EntityId ent_id)
{
    INCREMENT_ENTITY_MODIFICATION_COUNTER(PUBLIC, ent_rank, DECLARE_ENTITY);
    PartVector parts(1, &mesh_meta_data().universal_part());
    return internal_declare_entity(ent_rank, ent_id, parts);
}

Entity BulkData::declare_entity( EntityRank ent_rank , EntityId ent_id , Part& part)
{
    INCREMENT_ENTITY_MODIFICATION_COUNTER(PUBLIC, ent_rank, DECLARE_ENTITY);
    PartVector parts(1, &part);
    return internal_declare_entity( ent_rank, ent_id, parts);
}

Entity BulkData::declare_entity( EntityRank ent_rank , EntityId ent_id ,
                                 const PartVector & parts )
{
    INCREMENT_ENTITY_MODIFICATION_COUNTER(PUBLIC, ent_rank, DECLARE_ENTITY);
    return internal_declare_entity(ent_rank, ent_id, parts);
}

void BulkData::add_node_sharing( Entity node, int sharing_proc )
{
  // Only valid to specify sharing information for newly-created nodes
  ThrowRequire(entity_rank(node) == stk::topology::NODE_RANK);
  ThrowRequire(state(node) == Created);

  markEntity(node, IS_SHARED);
  entity_comm_map_insert(node, EntityCommInfo(stk::mesh::BulkData::SHARED, sharing_proc));
}

Entity BulkData::internal_declare_entity( EntityRank ent_rank , EntityId ent_id ,
                                 const PartVector & parts )
{
    m_check_invalid_rels = false;

  require_ok_to_modify();

  require_good_rank_and_id(ent_rank, ent_id);

  EntityKey key( ent_rank , ent_id );
  TraceIfWatching("stk::mesh::BulkData::internal_declare_entity", LOG_ENTITY, key);
  DiagIfWatching(LOG_ENTITY, key, "declaring entity with parts " << parts);

  std::pair< Entity , bool > result = internal_create_entity( key );

  Entity declared_entity = result.first;

  if ( !result.second) {
    // An existing entity, the owner must match.
    require_entity_owner( declared_entity , m_parallel_rank );
    DiagIfWatching(LOG_ENTITY, key, "existing entity: " << entity_key(declared_entity));
  }

  //------------------------------

  Part * const owns = & m_mesh_meta_data.locally_owned_part();

  PartVector rem ;
  PartVector add( parts );
  add.push_back( owns );

  internal_verify_and_change_entity_parts( declared_entity , add , rem );

  if ( result.second ) {
    this->set_parallel_owner_rank(declared_entity, m_parallel_rank);
    set_synchronized_count(declared_entity, m_sync_count);
    DiagIfWatching(LOG_ENTITY, key, "new entity: " << entity_key(declared_entity));
  }

  m_check_invalid_rels = true;

  return declared_entity ;
}

void BulkData::change_entity_id( EntityId id, Entity entity)
{
// When stk parallel is used within Fmwk, this assertion is violated
#ifndef SIERRA_MIGRATION
  ThrowAssertMsg(parallel_size() == 1,
                 "change_entity_id only supported in serial");
#endif

  EntityRank e_rank = entity_rank(entity);
  INCREMENT_ENTITY_MODIFICATION_COUNTER(PUBLIC, e_rank, CHANGE_ENTITY_ID);

  require_ok_to_modify();
  require_good_rank_and_id(e_rank, id);

  EntityKey new_key(e_rank,id);
  EntityKey old_key = entity_key(entity);

  internal_change_entity_key(old_key, new_key, entity);
}

void BulkData::internal_change_entity_key( EntityKey old_key, EntityKey new_key, Entity entity)
{
  INCREMENT_ENTITY_MODIFICATION_COUNTER(INTERNAL, old_key.rank(), CHANGE_ENTITY_ID);

  m_entity_repo.update_entity_key(new_key, old_key, entity);
  set_entity_key(entity, new_key);
  this->bucket(entity).getPartition()->set_flag_needs_to_be_sorted(true);
}

//----------------------------------------------------------------------

bool BulkData::destroy_entity( Entity entity, bool was_ghost )
{
    INCREMENT_ENTITY_MODIFICATION_COUNTER(PUBLIC, entity_rank(entity), DESTROY_ENTITY);
    return internal_destroy_entity(entity, was_ghost);
}

bool BulkData::internal_destroy_entity( Entity entity, bool was_ghost )
{
  const EntityKey key = entity_key(entity);
  TraceIfWatching("stk::mesh::BulkData::destroy_entity", LOG_ENTITY, key);
  DiagIfWatching(LOG_ENTITY, key, "entity state: " << key);

  require_ok_to_modify();

  m_check_invalid_rels = false;

  if (!is_valid(entity)) {
    m_check_invalid_rels = true;
    return false;
  }

  const bool ghost = was_ghost || in_receive_ghost(key);
  const EntityRank erank = entity_rank(entity);

  INCREMENT_ENTITY_MODIFICATION_COUNTER(INTERNAL, erank, DESTROY_ENTITY);

  const EntityRank end_rank = static_cast<EntityRank>(m_mesh_meta_data.entity_rank_count());
  for (EntityRank irank = static_cast<EntityRank>(erank + 1); irank != end_rank; ++irank) {
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
        internal_destroy_relation(entity, rel_entities[j], rel_ordinals[j]);
      }
    }
  }

  // If this is a ghosted entity, store key->local_offset so that local_offset can be
  // reused if the entity is recreated in the next aura-regen. This will prevent clients
  // from having their handles to ghosted entities go invalid when the ghost is refreshed.
  if ( ghost ) {
    DiagIfWatching(LOG_ENTITY, key, "Adding to ghost reuse map");
    m_ghost_reuse_map[key] = entity.local_offset();
  }

  // We need to save these items and call remove_entity AFTER the call to
  // destroy_later because remove_entity may destroy the bucket
  // which would cause problems in m_entity_repo.destroy_later because it
  // makes references to the entity's original bucket.

  // Need to invalidate Entity handles in comm-list
  EntityCommListInfoVector::iterator lb_itr =
    std::lower_bound(m_entity_comm_list.begin(), m_entity_comm_list.end(), key);
  if (lb_itr != m_entity_comm_list.end() && lb_itr->key == key) {
    lb_itr->entity = Entity();
  }

  remove_entity_callback(erank, bucket(entity).bucket_id(), bucket_ordinal(entity));

  m_entities_index.register_removed_key( key );

  bucket(entity).getPartition()->remove(entity);
  m_entity_repo.destroy_entity(key, entity );
  m_entity_states[entity.local_offset()] = Deleted;
  m_mark_entity[entity.local_offset()] = NOT_MARKED;
  m_closure_count[entity.local_offset()] = static_cast<uint16_t>(0u);
  if ( !ghost ) {
    m_deleted_entities_current_modification_cycle.push_front(entity.local_offset());
  }

  m_check_invalid_rels = true;
  return true ;
}

//----------------------------------------------------------------------

void BulkData::generate_new_entities(const std::vector<size_t>& requests,
                                 std::vector<Entity>& requested_entities)
// requests = number of nodes needed, number of elements needed, etc.
{
    Trace_("stk::mesh::BulkData::generate_new_entities");

    typedef stk::parallel::DistributedIndex::KeyType KeyType;
    typedef stk::parallel::DistributedIndex::KeyTypeVector KeyTypeVector;
    typedef std::vector<KeyTypeVector> RequestKeyVector;

    RequestKeyVector requested_key_types;

    m_entities_index.generate_new_keys(requests, requested_key_types);

    //generating 'owned' entities
    Part * const owns = &m_mesh_meta_data.locally_owned_part();

    std::vector<Part*> rem;
    std::vector<Part*> add;
    add.push_back(owns);

    requested_entities.clear();
    unsigned total_number_of_key_types_requested = 0;
    for(RequestKeyVector::const_iterator itr = requested_key_types.begin(); itr != requested_key_types.end(); ++itr)
    {
        const KeyTypeVector& key_types = *itr;
        for(KeyTypeVector::const_iterator kitr = key_types.begin(); kitr != key_types.end(); ++kitr)
        {
            ++total_number_of_key_types_requested;
        }
    }
    requested_entities.reserve(total_number_of_key_types_requested);
    addMeshEntities(requested_key_types, rem, add, requested_entities);
}

std::pair<Entity, bool> BulkData::internal_create_entity(EntityKey key, size_t preferred_offset)
{
    INCREMENT_ENTITY_MODIFICATION_COUNTER(INTERNAL, key.rank(), DECLARE_ENTITY);
    return m_entity_repo.internal_create_entity(key, preferred_offset);
}

void BulkData::addMeshEntities(const std::vector< stk::parallel::DistributedIndex::KeyTypeVector >& requested_key_types,
       const std::vector<Part*> &rem, const std::vector<Part*> &add, std::vector<Entity>& requested_entities)
{
    typedef stk::parallel::DistributedIndex::KeyType       KeyType;
    typedef stk::parallel::DistributedIndex::KeyTypeVector KeyTypeVector;
    typedef std::vector< KeyTypeVector > RequestKeyVector;

    for(RequestKeyVector::const_iterator itr = requested_key_types.begin(); itr != requested_key_types.end(); ++itr)
    {
        const KeyTypeVector & key_types = *itr;
        for(KeyTypeVector::const_iterator kitr = key_types.begin(); kitr != key_types.end(); ++kitr)
        {
            EntityKey key(static_cast<EntityKey::entity_key_t>((*kitr)));
            require_good_rank_and_id(key.rank(), key.id());
            std::pair<Entity, bool> result = internal_create_entity(key);

            //if an entity is declared with the declare_entity function in
            //the same modification cycle as the generate_new_entities
            //function, and it happens to generate a key that was declared
            //previously in the same cycle it is an error
            ThrowErrorMsgIf( ! result.second,
                    "Generated id " << key.id() << " of rank " << key.rank() <<
                    " which was already used in this modification cycle.");

            // A new application-created entity

            Entity new_entity = result.first;

            //add entity to 'owned' part
            internal_verify_and_change_entity_parts(new_entity, add, rem);
            requested_entities.push_back(new_entity);

            this->set_parallel_owner_rank(new_entity, m_parallel_rank);
            set_synchronized_count(new_entity, m_sync_count);
        }
    }
}


bool BulkData::in_shared(EntityKey key, int proc) const
{
  PairIterEntityComm sharing = entity_comm_map_shared(key);
  for ( ; !sharing.empty(); ++sharing ) {
    if ( proc == sharing->proc ) {
      return true ;
    }
  }
  return false ;
}

bool BulkData::is_ghosted_onto_proc( EntityKey key, int otherProc ) const
{
  const int proc = parallel_rank();
  const int owner_rank = entity_comm_map_owner(key);
  if ( proc == owner_rank )
  {
      for ( PairIterEntityComm ec = entity_comm_map(key); ! ec.empty() ; ++ec ) {
        if ( ec->ghost_id == 1 &&
             ec->proc  == otherProc &&
             ec->proc     != proc ) {
          return true;
        }
      }
  }
  return false;
}

bool BulkData::is_ghosted_onto_another_proc( EntityKey key ) const
{
  const int proc = parallel_rank();
  const int owner_rank = entity_comm_map_owner(key);
  if ( proc == owner_rank )
  {
      for ( PairIterEntityComm ec = entity_comm_map(key); ! ec.empty() ; ++ec ) {
        if ( ec->ghost_id == 1 &&
             ec->proc     != proc ) {
          return true;
        }
      }
  }
  return false;
}

bool BulkData::in_send_ghost( EntityKey key , int proc ) const
{
  const int owner_rank = entity_comm_map_owner(key);
  for ( PairIterEntityComm ec = entity_comm_map(key); ! ec.empty() ; ++ec ) {
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

  PairIterEntityComm ec = entity_comm_map(key);
  EntityCommInfoVector::const_iterator i =
    std::lower_bound( ec.begin(), ec.end() , tmp );

  return i != ec.end() && tmp == *i ;
}

void BulkData::comm_procs( EntityKey key, std::vector<int> & procs ) const
{
  procs.clear();
  for ( PairIterEntityComm ec = entity_comm_map(key); ! ec.empty() ; ++ec ) {
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
  for ( PairIterEntityComm ec = entity_comm_map_shared(key); ! ec.empty() ; ++ec ) {
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
  for ( PairIterEntityComm ec = entity_comm_map(key); ! ec.empty() ; ++ec ) {
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

void BulkData::new_bucket_caching(EntityRank rank, Bucket* new_bucket)
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

        const std::vector<FieldBase *>& allFields = mesh_meta_data().get_fields((stk::topology::rank_t) dst_rank);
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
  if (!m_keep_fields_updated) {
    return;
  }
  const std::vector<FieldBase *> &fields = mesh_meta_data().get_fields();
  m_field_data_manager->reinitialize_removed_entity_field_data(rank, bucket_id, bucket_ord, fields);
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

  const std::vector<FieldBase*>  fields = mesh_meta_data().get_fields();
  m_field_data_manager->deallocate_bucket_field_data(rank, bucket_id, capacity, fields);
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
        for ( int b = 0, be = field_set[outer_idx]->get_meta_data_for_field().size(); b < be; ++b) {
          if ( field_set[outer_idx]->get_meta_data_for_field()[b].m_bytes_per_entity > 0 ) {
            unsigned char* data_last = field_set[outer_idx]->get_meta_data_for_field()[b].m_data;
            for ( int s = 1; s < num_state; ++s ) {
              std::swap(field_set[outer_idx+s]->get_meta_data_for_field()[b].m_data, data_last);
            }
            field_set[outer_idx]->get_meta_data_for_field()[b].m_data = data_last;
          }
        }
      }

      for ( int s = 1; s < num_state; ++s )
      {
          m_field_data_manager->swap_fields(outer_idx+s, outer_idx);
      }

    }
}

void BulkData::reorder_buckets_callback(EntityRank rank, const std::vector<unsigned>& reorderedBucketIds)
{
  for (SelectorBucketMap::iterator itr = m_selector_to_buckets_map.begin(), end = m_selector_to_buckets_map.end();
       itr != end; ++itr) {
    TrackedBucketVector& cached_buckets = itr->second;
    std::sort(cached_buckets.begin(), cached_buckets.end(), BucketIdComparator());
  }

  if (!m_keep_fields_updated) {
    return;
  }

  const std::vector<FieldBase*>  fields = mesh_meta_data().get_fields();
  m_field_data_manager->reorder_bucket_field_data(rank, fields, reorderedBucketIds);
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
    EntityRank rank = static_cast<EntityRank>(i);
    out << "  All " << rank_names[i] << " entities:" << std::endl;

    const BucketVector& buckets = this->buckets(rank);
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
        for (EntityRank r = stk::topology::NODE_RANK, re = static_cast<EntityRank>(rank_names.size()); r < re; ++r) {
          if (connectivity_map().valid(static_cast<EntityRank>(rank), r)) {
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

            if(static_cast<unsigned>(field->entity_rank()) != bucket->entity_rank()) continue;

            FieldMetaData field_meta_data = field->get_meta_data_for_field()[bucket->bucket_id()];

            unsigned data_size = field_meta_data.m_bytes_per_entity;
            if (data_size > 0) { // entity has this field?
              void* data = field_meta_data.m_data + field_meta_data.m_bytes_per_entity * b_ord;
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
    ThrowRequireMsg(mesh.connectivity_map().valid(stk::topology::NODE_RANK, stk::topology::ELEMENT_RANK), "ERROR, stk::mesh::get_connectivity was called but upward node->element connections are disabled. For get_connectivity to find non-stored connectivities, node->element connections must be enabled.");

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
    data.second += buckets(rank).size();
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


//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------

//----------------------------------------------------------------------

void BulkData::require_valid_relation( const char action[] ,
                                       const BulkData & mesh ,
                                       const Entity e_from ,
                                       const Entity e_to )
{
  const bool error_type      = mesh.entity_rank(e_from) <= mesh.entity_rank(e_to);
  const bool error_nil_from  = !mesh.is_valid(e_from);
  const bool error_nil_to    = !mesh.is_valid(e_to);

  if ( error_type || error_nil_from || error_nil_to ) {
    std::ostringstream msg ;

    msg << "Could not " << action << " relation from entity "
        << print_entity_key(MetaData::get(mesh), mesh.entity_key(e_from)) << " to entity "
        << print_entity_key(MetaData::get(mesh), mesh.entity_key(e_to)) << "\n";

    ThrowErrorMsgIf( error_nil_from  || error_nil_to,
                     msg.str() << ", entity was destroyed");
    ThrowErrorMsgIf( error_type, msg.str() <<
                     "A relation must be from higher to lower ranking entity");
  }
}


//----------------------------------------------------------------------
bool BulkData::internal_declare_relation(Entity e_from, Entity e_to,
                                         RelationIdentifier local_id,
                                         unsigned sync_count, bool is_back_relation,
                                         Permutation permut)
{
    INCREMENT_MODIFICATION_COUNTER(INTERNAL, DECLARE_RELATION);

  TraceIfWatching("stk::mesh::BuilkData::internal_declare_relation", LOG_ENTITY, entity_key(e_from));

  const MeshIndex& idx = mesh_index(e_from);

  bool modified = idx.bucket->declare_relation(idx.bucket_ordinal, e_to, static_cast<ConnectivityOrdinal>(local_id),
                                               permut);

  if (modified) {
    set_synchronized_count( e_from, sync_count );
  }
  return modified;
}

void BulkData::declare_relation( Entity e_from ,
                                 Entity e_to ,
                                 const RelationIdentifier local_id ,
                                 Permutation permut)
{
  INCREMENT_MODIFICATION_COUNTER(PUBLIC, DECLARE_RELATION);
  OrdinalVector ordinal_scratch;
  PartVector part_scratch;
  internal_declare_relation(e_from, e_to, local_id, permut, ordinal_scratch, part_scratch);
}

void BulkData::declare_relation( Entity e_from ,
                                 Entity e_to ,
                                 const RelationIdentifier local_id ,
                                 Permutation permut,
                                 OrdinalVector& ordinal_scratch,
                                 PartVector& part_scratch)
{
    INCREMENT_MODIFICATION_COUNTER(PUBLIC, DECLARE_RELATION);
    internal_declare_relation(e_from, e_to, local_id, permut, ordinal_scratch, part_scratch);
}

void BulkData::internal_declare_relation( Entity e_from ,
                                 Entity e_to ,
                                 const RelationIdentifier local_id ,
                                 Permutation permut,
                                 OrdinalVector& ordinal_scratch,
                                 PartVector& part_scratch)
{
  TraceIfWatching("stk::mesh::BulkData::declare_relation", LOG_ENTITY, entity_key(e_from));
  TraceIfWatchingDec("stk::mesh::BulkData::declare_relation", LOG_ENTITY, entity_key(e_to), 1);
  DiagIfWatching(LOG_ENTITY, entity_key(e_from),
                 "from: " << entity_key(e_from) << ";  " <<
                 "to: " << entity_key(e_to) << ";  " <<
                 "id: " << local_id);
  DiagIfWatching(LOG_ENTITY, entity_key(e_to),
                 "from: " << entity_key(e_from) << ";  " <<
                 "to: " << entity_key(e_to) << ";  " <<
                 "id: " << local_id);

  require_ok_to_modify();

  require_valid_relation( "declare" , *this , e_from , e_to );

  // TODO: Don't throw if exact relation already exists, that should be a no-op.
  // Should be an exact match if relation of local_id already exists (e_to should be the same).
  bool is_converse = false;
  bool caused_change_fwd = internal_declare_relation(e_from, e_to, local_id, m_sync_count,
                                                     is_converse, permut);

  //TODO: check connectivity map
  // Relationships should always be symmetrical
  if ( caused_change_fwd ) {

    const bool higher_order_relation = stk::topology::ELEMENT_RANK < entity_rank(e_from);
    if (    higher_order_relation
         || m_bucket_repository.connectivity_map()(entity_rank(e_to), entity_rank(e_from)) != stk::mesh::INVALID_CONNECTIVITY_TYPE
       )
    {
      // the setup for the converse relationship works slightly differently
      is_converse = true;
      internal_declare_relation(e_to, e_from, local_id, m_sync_count, is_converse, permut );
    }
  }

  // It is critical that the modification be done AFTER the relations are
  // added so that the propagation can happen correctly.
  if ( caused_change_fwd ) {
    this->mark_entity_and_upward_related_entities_as_modified(e_to);
    this->mark_entity_and_upward_related_entities_as_modified(e_from);
  }

  OrdinalVector empty ;

  // Deduce and set new part memberships:
  ordinal_scratch.clear();

  induced_part_membership(*this, e_from, empty, entity_rank(e_to), ordinal_scratch );

  PartVector emptyParts;
  part_scratch.clear();
  for(unsigned ipart=0; ipart<ordinal_scratch.size(); ++ipart) {
    part_scratch.push_back(&m_mesh_meta_data.get_part(ordinal_scratch[ipart]));
  }

  internal_change_entity_parts( e_to , part_scratch , emptyParts );
}

//----------------------------------------------------------------------

void BulkData::internal_declare_relation( Entity entity ,
                                 const std::vector<Relation> & rel )
{
  OrdinalVector ordinal_scratch;
  PartVector part_scratch;
  internal_declare_relation(entity, rel, ordinal_scratch, part_scratch);
}

void BulkData::internal_declare_relation( Entity entity ,
                                 const std::vector<Relation> & rel,
                                 OrdinalVector& ordinal_scratch,
                                 PartVector& part_scratch )
{
  require_ok_to_modify();

  const unsigned erank = entity_rank(entity);

  std::vector<Relation>::const_iterator i ;
  for ( i = rel.begin() ; i != rel.end() ; ++i ) {
    Entity e = i->entity();
    const unsigned n = i->relation_ordinal();
    const Permutation permut = static_cast<Permutation>(i->getOrientation());
    if ( entity_rank(e) < erank ) {
      internal_declare_relation( entity , e , n, permut, ordinal_scratch, part_scratch );
    }
    else if ( erank < entity_rank(e) ) {
      internal_declare_relation( e , entity , n, permut, ordinal_scratch, part_scratch );
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
    INCREMENT_MODIFICATION_COUNTER(PUBLIC, DESTROY_RELATION);
    return internal_destroy_relation(e_from, e_to,  local_id);
}

bool BulkData::internal_destroy_relation( Entity e_from ,
                                 Entity e_to,
                                 const RelationIdentifier local_id )
{
    INCREMENT_MODIFICATION_COUNTER(INTERNAL, DESTROY_RELATION);

  TraceIfWatching("stk::mesh::BulkData::destroy_relation", LOG_ENTITY, entity_key(e_from));
  TraceIfWatchingDec("stk::mesh::BulkData::destroy_relation", LOG_ENTITY, entity_key(e_to), 1);
  DiagIfWatching(LOG_ENTITY, entity_key(e_from),
                 "from: " << entity_key(e_from) << ";  " <<
                 "to: " << entity_key(e_to) << ";  " <<
                 "id: " << local_id);
  DiagIfWatching(LOG_ENTITY, entity_key(e_to),
                 "from: " << entity_key(e_from) << ";  " <<
                 "to: " << entity_key(e_to) << ";  " <<
                 "id: " << local_id);

  require_ok_to_modify();

  require_valid_relation( "destroy" , *this , e_from , e_to );

  const EntityRank end_rank = static_cast<EntityRank>(m_mesh_meta_data.entity_rank_count());
  const EntityRank e_to_entity_rank = entity_rank(e_to);
  const EntityRank e_from_entity_rank = entity_rank(e_from);

  //------------------------------
  // When removing a relationship may need to
  // remove part membership and set field relation pointer to NULL

  m_check_invalid_rels = false; // OK to have gaps when deleting

  if ( parallel_size() < 2 || entity_comm_map_shared(entity_key(e_to)).empty() ) {

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
      EntityVector temp_entities;
      std::vector<ConnectivityOrdinal> temp_ordinals;
      Entity const* rel_entities = NULL;
      ConnectivityOrdinal const* rel_ordinals = NULL;
      int num_rels = 0;
      for (EntityRank irank = static_cast<EntityRank>(e_to_entity_rank + 1); irank < end_rank; ++irank) {
        if (connectivity_map().valid(e_to_entity_rank, irank)) {
          num_rels     = num_connectivity(e_to, irank);
          rel_entities = begin(e_to, irank);
          rel_ordinals = begin_ordinals(e_to, irank);
        }
        else {
          num_rels     = get_connectivity(*this, e_to, irank, temp_entities, temp_ordinals);
          rel_entities = &*temp_entities.begin();
          rel_ordinals = &*temp_ordinals.begin();
        }

        for (int j = 0; j < num_rels; ++j) {
          ThrowAssertMsg(is_valid(rel_entities[j]), "Error, entity " << e_to.local_offset() << " has invalid back-relation for ordinal: "
                         << rel_ordinals[j] << " to rank: " << irank << ", target entity is: " << rel_entities[j].local_offset());
          if ( !(rel_entities[j] == e_from && rel_ordinals[j] == static_cast<ConnectivityOrdinal>(local_id) ) )
          {
            induced_part_membership(*this, rel_entities[j], empty, e_to_entity_rank, keep,
                                    false /*Do not look at supersets*/);
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
          induced_part_membership(*this, e_from, keep, e_to_entity_rank, del,
                                  false /*Do not look at supersets*/);
          break; // at most 1 relation can match our specification
        }
      }
    }

    if ( !del.empty() ) {

      PartVector delParts, emptyParts;
      delParts.reserve(del.size());
      for(unsigned ipart=0; ipart<del.size(); ++ipart) {
        delParts.push_back(&m_mesh_meta_data.get_part(del[ipart]));
      }

      internal_change_entity_parts( e_to , emptyParts , delParts );
    }
  }

  //delete relations from the entities
  bool caused_change_fwd = bucket(e_from).destroy_relation(e_from, e_to, local_id);

  // Relationships should always be symmetrical
  if ( caused_change_fwd &&
       (e_to_entity_rank > stk::topology::ELEMENT_RANK || e_from_entity_rank > stk::topology::ELEMENT_RANK ||
        connectivity_map().valid(entity_rank(e_to), entity_rank(e_from))) ) {
    bucket(e_to).destroy_relation(e_to, e_from, local_id);
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
//----------------------------------------------------------------------

//----------------------------------------------------------------------
//ParallelVerify
//----------------------------------------------------------------------

namespace {

//#define DEBUG_PRINT_COMM_LIST
//#define DEBUG_PRINT_COMM_LIST_UNPACK

#ifdef DEBUG_PRINT_COMM_LIST

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

      for ( PairIterEntityComm ec = mesh.entity_comm_map(i->key); ! ec.empty() ; ++ec ) {
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

    BABBLE_STK_PARALLEL_COMM(M.parallel(), "          comm_mesh_verify_parallel_consistency calling allocate_buffers");
    all.allocate_buffers( all.parallel_size() / 4 );

    pack_owned_verify( all , M );

    BABBLE_STK_PARALLEL_COMM(M.parallel(), "          comm_mesh_verify_parallel_consistency calling communicate");
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
  const PairIterEntityComm ec = bulk.entity_comm_map(bulk.entity_key(entity));
  const size_t n = ec.size();
  for ( size_t i = 1 ; i < n ; ++i ) {
    if ( ! ( ec[i-1] < ec[i] ) ) {
      return false ;
    }
  }
  return true ;
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
                << "shares: " << shares << " has entity key " << M.entity_key(entity).m_value << std::endl;
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
      error_log << "problem: entity is both recv ghost and in owned_closure;"<<std::endl;
      this_result = false ;
    }
    if ( ! owned_closure && ! recv_ghost ) {
      error_log << "problem: entity is neither recv_ghost nor in owned_closure;"<<std::endl;
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
      PairIterEntityComm ip = M.entity_comm_map_shared(M.entity_key(entity));
      for ( ; ! ip.empty() && p_owner != ip->proc ; ++ip );
      if ( ip.empty() ) {
        error_log << "problem: entity shared-not-owned, but entity_comm_map_shared returns empty;" << std::endl;
        this_result = false ;
      }
    }

    if ( shares || recv_ghost || send_ghost ) { ++comm_count ; }

    if ( ! this_result ) {
      result = false ;
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
          if ( comm[kk].ghost_id == 1 && comm[kk].proc == share_proc )
          {
              std::cerr << "[" << mesh.parallel_rank() << "] bad news for entity " << i->key << " with id: " << mesh.identifier(i->entity) << std::endl;
              stk::mesh::Entity const * nodes = mesh.begin_nodes(i->entity);
              stk::mesh::Entity const * nodes_end = mesh.end_nodes(i->entity);
              std::cerr << "( ";
              for (; nodes != nodes_end; nodes++)
              {
                  std::cerr << mesh.identifier(*nodes) << " ";
              }
              std::cerr << ") " << std::endl;
          }

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

} // namespace<>


//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------
//Owner
//----------------------------------------------------------------------

//----------------------------------------------------------------------

namespace {

// Given an entity, if it's a ghost, insert the closure of the entity
// into work_list.
void insert_closure_ghost(const BulkData& mesh, Entity const entity ,
                           const int proc_local ,
                           std::set<EntityKey> & work_list )
{
  bool entityIsAGhost = ! in_owned_closure( mesh, entity , proc_local );
  if ( entityIsAGhost ) {
    // This entity is a ghost, put it on the work_list
    // along with all ghosts in its closure

    const bool was_inserted = work_list.insert(mesh.entity_key(entity)).second;

    if ( was_inserted ) {
      // This ghost entity is new to the list, traverse its closure.

      const unsigned erank = mesh.entity_rank(entity);

      // Recurse over downward relations
      for (EntityRank irank = stk::topology::BEGIN_RANK; irank < erank; ++irank)
      {
        Entity const *irels_j = mesh.begin(entity, irank);
        Entity const *irels_e = mesh.end(entity, irank);
        for (; irels_j != irels_e; ++irels_j)
        {
          insert_closure_ghost(mesh, *irels_j, proc_local, work_list);
        }
      }
    }
  }
}

// Given an entity, insert the closures of every entity that has this entity
// in its closure. Only ghosts will be inserted.
void insert_transitive_ghost(const BulkData& mesh, Entity const entity ,
                              const int proc_local ,
                              std::set<EntityKey> & work_list )
{
  insert_closure_ghost(mesh, entity , proc_local , work_list );

  // Transitive:
  // If this entity is a member of another entity's closure
  // then that other entity is part of the traversal.

  const EntityRank erank = mesh.entity_rank(entity);

  // Recurse over upward relations
  const EntityRank end_rank = static_cast<EntityRank>(mesh.mesh_meta_data().entity_rank_count());
  EntityVector temp_entities;
  Entity const* rels = NULL;
  int num_rels = 0;
  for (EntityRank irank = static_cast<EntityRank>(erank + 1); irank < end_rank; ++irank)
  {
    if (mesh.connectivity_map().valid(erank, irank)) {
      num_rels = mesh.num_connectivity(entity, irank);
      rels     = mesh.begin(entity, irank);
    }
    else {
      num_rels = get_connectivity( mesh, entity, irank, temp_entities);
      rels     = &*temp_entities.begin();
    }

    for (int r = 0; r < num_rels; ++r)
    {
      insert_transitive_ghost(mesh, rels[r] , proc_local , work_list );
    }
  }
}

//----------------------------------------------------------------------

// Add EntityProc pairs to send_list for every entity in the closure of the
// entity in send_entry. All these entities will be sent to the same proc as
// the original send_entry.
void insert_closure_send(
  const BulkData &mesh,
  const EntityProc                  send_entry ,
  std::set<EntityProc,EntityLess> & send_list )
{
  ThrowRequireMsg( mesh.is_valid(send_entry.first),
                   "Cannot send destroyed entity");

  std::pair< std::set<EntityProc,EntityLess>::iterator , bool >
    result = send_list.insert( send_entry );

  if ( result.second ) {
    // First time this entity was inserted into the send_list.

    const unsigned erank  = mesh.entity_rank(send_entry.first);
    const Bucket &ebucket = mesh.bucket(send_entry.first);
    const Ordinal ebordinal = mesh.bucket_ordinal(send_entry.first);

    // Recurse over downward relations
    for (EntityRank irank = stk::topology::BEGIN_RANK; irank < erank; ++irank)
    {
      Entity const *rels_itr = ebucket.begin(ebordinal, irank);
      Entity const *rels_end= ebucket.end(ebordinal, irank);
      for (; rels_itr != rels_end; ++rels_itr)
      {
        const EntityProc rel_send_entry( *rels_itr, send_entry.second );
        insert_closure_send(mesh, rel_send_entry , send_list );
      }
    }
  }
}

//----------------------------------------------------------------------

bool member_of_owned_closure(const BulkData& mesh, const Entity e , const int p_rank )
{
  if (p_rank == mesh.parallel_owner_rank(e)) {
    return true;
  }

  const EntityRank erank = mesh.entity_rank(e);
  const EntityRank end_rank = static_cast<EntityRank>(mesh.mesh_meta_data().entity_rank_count());

  // Any higher ranking entities locally owned?
  EntityVector temp_entities;
  Entity const* rels = NULL;
  int num_rels = 0;
  for (EntityRank irank = static_cast<EntityRank>(end_rank - 1); irank > erank; --irank)
  {
    if (mesh.connectivity_map().valid(erank, irank)) {
      num_rels = mesh.num_connectivity(e, irank);
      rels     = mesh.begin(e, irank);
    }
    else {
      num_rels = get_connectivity( mesh, e, irank, temp_entities);
      rels     = &*temp_entities.begin();
    }

    for (int r = 0; r < num_rels; ++r) {
      if (p_rank == mesh.parallel_owner_rank(rels[r]) ||  member_of_owned_closure(mesh, rels[r], p_rank) ) {
        return true;
      }
    }
  }

  return false;
}

//----------------------------------------------------------------------

// Given a vector of local ownership changes, remove duplicates and
// sanity check.
void clean_and_verify_parallel_change(
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

void generate_parallel_change( const BulkData & mesh ,
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
      BABBLE_STK_PARALLEL_COMM(mesh.parallel(), "          generate_parallel_change calling allocate_buffers");
      comm.allocate_buffers( p_size / 4 , 0 );
    }
    else { // communication phase
      BABBLE_STK_PARALLEL_COMM(mesh.parallel(), "          generate_parallel_change calling communicate");
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

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void BulkData::change_entity_owner( const std::vector<EntityProc> & arg_change )
{
    INCREMENT_MODIFICATION_COUNTER(PUBLIC, CHANGE_ENTITY_OWNER);
    internal_change_entity_owner(arg_change);
}

void BulkData::internal_change_entity_owner( const std::vector<EntityProc> & arg_change )
{
    INCREMENT_MODIFICATION_COUNTER(INTERNAL, CHANGE_ENTITY_OWNER);

  Trace_("stk::mesh::BulkData::change_entity_owner");
  DiagIf(LOG_ENTITY, "arg_change: " << arg_change);

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
  clean_and_verify_parallel_change( *this , local_change );

  //----------------------------------------
  // Parallel synchronous determination of changing shared and ghosted.

  // The two vectors below will contain changes to ghosted and shared
  // entities on this process coming from change-entity-owner requests
  // on other processes.
  std::vector<EntityProc> ghosted_change ;
  std::vector<EntityProc> shared_change ;

  generate_parallel_change( *this , local_change ,
                            shared_change , ghosted_change );

  //------------------------------
  // Have enough information to delete all effected ghosts.
  // If the closure of a ghost contains a changing entity
  // then that ghost must be deleted.
  // Request that all ghost entities in the closure of the ghost be deleted.

  typedef std::set<EntityProc,EntityLess> EntityProcSet;

  // Compute the closure of all the locally changing entities
  std::set<EntityProc,EntityLess> send_closure(EntityLess(*this)) ;
  for ( std::vector<EntityProc>::iterator
        i = local_change.begin() ; i != local_change.end() ; ++i ) {
    insert_closure_send(*this, *i , send_closure );
  }

  // Calculate all the ghosts that are impacted by the set of ownership
  // changes. We look at ghosted, shared, and local changes looking for ghosts
  // that are either in the closure of the changing entity, or have the
  // changing entity in their closure. All modified ghosts will be removed.
  {
    std::set<EntityKey> modified_ghosts;

    for ( std::vector<EntityProc>::const_iterator
          i = ghosted_change.begin() ; i != ghosted_change.end() ; ++i ) {
      insert_transitive_ghost(*this, i->first , m_parallel_rank , modified_ghosts );
    }

    for ( std::vector<EntityProc>::const_iterator
          i = shared_change.begin() ; i != shared_change.end() ; ++i ) {
      insert_transitive_ghost(*this, i->first , m_parallel_rank , modified_ghosts );
    }

    for ( EntityProcSet::iterator
          i = send_closure.begin() ; i != send_closure.end() ; ++i ) {
      insert_transitive_ghost(*this, i->first , m_parallel_rank , modified_ghosts );
    }

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

      const bool changed = this->set_parallel_owner_rank( entity, i->second );
      if (changed) {
        internal_change_owner_in_comm_data(entity_key(entity), i->second);
      }
    }

    for ( std::vector<EntityProc>::iterator
          i = shared_change.begin() ; i != shared_change.end() ; ++i ) {
      Entity entity = i->first;
      const bool changed = this->set_parallel_owner_rank( entity, i->second );
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

        const bool changed = this->set_parallel_owner_rank( entity, owner );
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
        if ( ! member_of_owned_closure(*this, *i , p_rank ) ) {
          ThrowRequireMsg( internal_destroy_entity( *i ),
                           "Failed to destroy entity " << identifier(*i) );
        }
      }
    }

    send_closure.clear(); // Has been invalidated
  }
}

//----------------------------------------------------------------------


//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------
//Ghosting
//----------------------------------------------------------------------

//----------------------------------------------------------------------

Ghosting & BulkData::create_ghosting( const std::string & name )
{
    INCREMENT_MODIFICATION_COUNTER(PUBLIC, CREATE_GHOSTING);
    return internal_create_ghosting(name);
}

Ghosting & BulkData::internal_create_ghosting( const std::string & name )
{
    INCREMENT_MODIFICATION_COUNTER(INTERNAL, CREATE_GHOSTING);

  require_ok_to_modify();

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

  Ghosting * const g =
    new Ghosting( *this , name , m_ghosting.size() , m_sync_count );

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

namespace {

void insert_transitive_closure( BulkData& bulk_data,
                                std::set<EntityProc,EntityLess> & new_send ,
                                const EntityProc & entry );

void comm_recv_to_send(
  BulkData & mesh ,
  const std::set< EntityKey > & entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs ,
        std::set< EntityProc , EntityLess > & entitiesToGhostOntoOtherProcessors );

void comm_sync_send_recv(
  BulkData & mesh ,
  std::set< EntityProc , EntityLess > & new_send ,
  std::set< EntityKey > & new_recv );

} // namespace <>

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void BulkData::destroy_ghosting( Ghosting& ghost_layer )
{
  INCREMENT_MODIFICATION_COUNTER(PUBLIC, DESTROY_GHOSTING);
  std::vector<EntityKey> receive_list;
  ghost_layer.receive_list(receive_list);
  internal_verify_inputs_and_change_ghosting(ghost_layer, std::vector<stk::mesh::EntityProc>(), receive_list);
}

//----------------------------------------------------------------------

void BulkData::destroy_all_ghosting()
{
  INCREMENT_MODIFICATION_COUNTER(PUBLIC, DESTROY_ALL_GHOSTING);

  Trace_("stk::mesh::BulkData::destroy_all_ghosting");

  require_ok_to_modify();

  // Clear Ghosting data

  for ( std::vector<Ghosting*>::iterator
        ig = m_ghosting.begin() ; ig != m_ghosting.end() ; ++ig ) {
    Ghosting & gh = **ig ;
    gh.m_sync_count = m_sync_count ;
  }

  // Iterate backwards so as not to invalidate a closure.

  for ( EntityCommListInfoVector::reverse_iterator
        i =  m_entity_comm_list.rbegin() ;
        i != m_entity_comm_list.rend() ; ++i) {

    if ( in_receive_ghost( i->key ) ) {
        entity_comm_map_clear( i->key );
        internal_destroy_entity( i->entity );
      i->key = EntityKey();
      i->entity_comm = NULL;
    }
    else {
        entity_comm_map_clear_ghosting(i->key);
      if ( entity_comm_map(i->key).empty() ) {
        i->key = EntityKey();
        i->entity_comm = NULL;
      }
    }
  }

  EntityCommListInfoVector::iterator i =
    std::remove_if( m_entity_comm_list.begin() ,
                    m_entity_comm_list.end() , IsInvalid() );

  m_entity_comm_list.erase( i , m_entity_comm_list.end() );
}

//----------------------------------------------------------------------

void BulkData::change_ghosting(
  Ghosting & ghosts ,
  const std::vector<EntityProc> & add_send ,
  const std::vector<EntityKey> & remove_receive )
{
    INCREMENT_MODIFICATION_COUNTER(PUBLIC, CHANGE_GHOSTING);
    internal_verify_inputs_and_change_ghosting(ghosts, add_send, remove_receive);
}

void BulkData::internal_verify_inputs_and_change_ghosting(
  Ghosting & ghosts ,
  const std::vector<EntityProc> & add_send ,
  const std::vector<EntityKey> & remove_receive )
{
  Trace_("stk::mesh::BulkData::change_ghosting");

  //----------------------------------------
  // Verify inputs:

  require_ok_to_modify();

  const bool ok_mesh  = & BulkData::get(ghosts) == this ;
  const bool ok_ghost = 1 < ghosts.ordinal();
  bool ok_add    = true ;
  bool ok_remove = true ;

  // Verify all 'add' are locally owned.

  for ( std::vector<EntityProc>::const_iterator
        i = add_send.begin() ; ok_add && i != add_send.end() ; ++i ) {
    ok_add = parallel_owner_rank(i->first) == parallel_rank();
  }

  // Verify all 'remove' are members of the ghosting.

  for ( std::vector<EntityKey>::const_iterator
        i = remove_receive.begin() ;
        ok_remove && i != remove_receive.end() ; ++i ) {
    ok_remove = in_receive_ghost( ghosts , *i );
  }

  int ok = ok_mesh && ok_ghost && ok_add && ok_remove ;

  all_reduce( parallel() , ReduceMin<1>( & ok ) );

  if ( 0 == ok ) {
    std::ostringstream msg ;
    msg << "For ghosts " << ghosts.name() << ", " ;
    if ( ! ok_mesh )  { msg << " : Mesh does not own this ghosting" ; }
    if ( ! ok_ghost ) { msg << " : Cannot modify this ghosting" ; }
    if ( ! ok_add ) {
      msg << " : Not owned add {" ;
      for ( std::vector<EntityProc>::const_iterator
            i = add_send.begin() ; i != add_send.end() ; ++i ) {
        if ( parallel_owner_rank(i->first) != parallel_rank() ) {
          msg << " " << identifier(i->first);
        }
      }
      msg << " }" ;
    }
    if ( ! ok_remove ) {
      msg << " : Not in ghost receive {" ;
      for ( std::vector<EntityKey>::const_iterator
            i = remove_receive.begin() ; i != remove_receive.end() ; ++i ) {
        if ( ! in_receive_ghost( ghosts , *i ) ) {
          msg << " " << i->id();
        }
      }
    }

    ThrowErrorMsg( msg.str() );
  }
  //----------------------------------------
  // Change the ghosting:

  internal_change_ghosting( ghosts , add_send , remove_receive );
}

//----------------------------------------------------------------------

void BulkData::ghost_entities_and_fields(Ghosting & ghosting, const std::set<EntityProc , EntityLess>& entitiesToGhostOntoOtherProcessors)
{
    const size_t record_entity_comm_size_before_changing_it = m_entity_comm_list.size();
    const int p_size = parallel_size() ;

    CommAll comm( m_parallel_machine );

    for ( int phase = 0; phase < 2; ++phase ) {
      for ( std::set< EntityProc , EntityLess >::iterator
              j = entitiesToGhostOntoOtherProcessors.begin(); j != entitiesToGhostOntoOtherProcessors.end() ; ++j ) {

        Entity entity = j->first;
        const int proc = j->second;

        if ( ! in_ghost( ghosting , entity_key(entity) , proc ) ) {
          // Not already being sent , must send it.
          CommBuffer & buf = comm.send_buffer( proc );
          buf.pack<unsigned>( entity_rank(entity) );
          pack_entity_info(*this, buf , entity );
          pack_field_values(*this, buf , entity );

          if (phase == 1) {
            entity_comm_map_insert(entity, EntityCommInfo(ghosting.ordinal(), proc));
            const EntityComm* entity_comm = m_entity_comm_map.entity_comm(entity_key(entity));
            EntityCommListInfo comm_info = {entity_key(entity), entity, parallel_owner_rank(entity), entity_comm};
            m_entity_comm_list.push_back( comm_info );
          }
        }
      }

      if (phase == 0) {
        comm.allocate_buffers( p_size / 4 );
      }
      else {
        comm.communicate();
      }
    }

    std::ostringstream error_msg ;
    int error_count = 0 ;
    OrdinalVector ordinal_scratch;
    PartVector part_scratch;

    const MetaData & meta = m_mesh_meta_data ;
    const unsigned rank_count = meta.entity_rank_count();

    for ( unsigned rank = 0 ; rank < rank_count ; ++rank ) {
      for ( int p = 0 ; p < p_size ; ++p ) {
        CommBuffer & buf = comm.recv_buffer(p);
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

          PartVector parts ;
          std::vector<Relation> relations ;
          EntityKey key ;
          int owner = ~0u ;

          unpack_entity_info( buf, *this, key, owner, parts, relations );

          // Must not have the locally_owned_part or globally_shared_part

          remove( parts , meta.locally_owned_part() );
          remove( parts , meta.globally_shared_part() );

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

          std::pair<Entity ,bool> result = internal_create_entity( key, use_this_offset );

          Entity entity = result.first;
          const bool created   = result.second ;

          require_entity_owner( entity , owner );

          internal_change_entity_parts( entity , parts , PartVector() );

          if ( created ) {
            log_created_parallel_copy( entity );
            this->set_parallel_owner_rank( entity, owner);
          }

          internal_declare_relation( entity , relations, ordinal_scratch, part_scratch );

          if ( ! unpack_field_values(*this, buf , entity , error_msg ) ) {
            ++error_count ;
          }

          const EntityCommInfo tmp( ghosting.ordinal() , owner );

          if ( entity_comm_map_insert(entity, tmp) ) {
            const EntityComm* entity_comm = m_entity_comm_map.entity_comm(entity_key(entity));
            EntityCommListInfo comm_info = {entity_key(entity), entity, parallel_owner_rank(entity), entity_comm};
            m_entity_comm_list.push_back( comm_info );
          }
        }
      }
    }

    if (parallel_size() > 1) {
      all_reduce( m_parallel_machine , ReduceSum<1>( & error_count ) );
    }

    ThrowErrorMsgIf( error_count, error_msg.str() );

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

      internal_sync_comm_list_owners();
    }
    ghosting.m_sync_count = m_sync_count ;
}

void BulkData::internal_change_ghosting(
  Ghosting & ghosting ,
  const std::vector<EntityProc> & add_send ,
  const std::vector<EntityKey> & remove_receive,
  bool is_full_regen)
{
  Trace_("stk::mesh::BulkData::internal_change_ghosting");

  BABBLE_STK_PARALLEL_COMM(m_parallel_machine, "      entered internal_change_ghosting");

  INCREMENT_MODIFICATION_COUNTER(INTERNAL, CHANGE_GHOSTING);

  //------------------------------------
  // Copy ghosting lists into more efficiently edited container.
  // The send and receive lists must be in entity rank-order.

  std::set<EntityProc , EntityLess> entitiesToGhostOntoOtherProcessors(EntityLess(*this));
  std::set<EntityKey>               entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs;

  //------------------------------------
  // Insert the current ghost receives and then remove from that list.

  // This if-check is an optimization; if doing a full regen
  // then we are removing all ghosting information and new_recv should
  // be left empty.
  if ( !is_full_regen ) {

    // Iterate over all entities with communication information, adding
    // the entity if it's a ghost on this process. new_recv will contain
    // all ghosts on this process by the end of the loop.
    for ( EntityCommListInfoVector::const_iterator
          i = comm_list().begin() ; i != comm_list().end() ; ++i ) {
      if ( in_receive_ghost( ghosting , i->key ) ) {
        entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs.insert( i->key );
      }
    }

    // Remove any entities that are in the remove list.

    for ( std::vector<EntityKey>::const_iterator
          i = remove_receive.begin() ; i != remove_receive.end() ; ++i ) {
      entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs.erase( *i );
    }

    // Keep the closure of the remaining received ghosts.
    // Working from highest-to-lowest key (rank entity type)
    // results in insertion of the closure because
    // inserted entities will get looped over after they are inserted.

    // Insertion will not invalidate the associative container's iterator.

    for ( std::set<EntityKey>::reverse_iterator
          i = entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs.rbegin() ; i != entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs.rend() ; ++i) {
      const unsigned erank = i->rank();

      Entity e = get_entity(*i); // Could be performance issue? Not if you're just doing full regens

      for (EntityRank irank = stk::topology::BEGIN_RANK; irank < erank; ++irank) {
        Entity const *rels_i = begin(e, irank);
        Entity const *rels_e = end(e, irank);
        for (; rels_i != rels_e; ++rels_i)
        {
          if ( is_valid(*rels_i) && in_receive_ghost( ghosting , entity_key(*rels_i) ) )
          {
            entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs.insert( entity_key(*rels_i) );
          }
        }
      }
    }
  }

  //  Initialize the new_send from the new_recv
  comm_recv_to_send( *this , entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs , entitiesToGhostOntoOtherProcessors );

  //------------------------------------
  // Add the specified entities and their closure to the send ghosting

  for ( std::vector< EntityProc >::const_iterator
        i = add_send.begin() ; i != add_send.end() ; ++i ) {
    insert_transitive_closure( *this, entitiesToGhostOntoOtherProcessors , *i );
  }

  // Synchronize the send and receive list.
  // If the send list contains a not-owned entity
  // inform the owner and receiver to add that entity
  // to their ghost send and receive lists.

  comm_sync_send_recv( *this , entitiesToGhostOntoOtherProcessors , entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs );

  // The new_send list is now parallel complete and parallel accurate
  // The new_recv has those ghost entities that are to be kept.
  //------------------------------------
  // Remove the ghost entities that will not remain.
  // If the last reference to the receive ghost entity then delete it.

  bool removed = false ;

  for ( EntityCommListInfoVector::reverse_iterator
        i = m_entity_comm_list.rbegin() ; i != m_entity_comm_list.rend() ; ++i) {

    const bool is_owner = i->owner == m_parallel_rank ;
    const bool remove_recv = ( ! is_owner ) &&
                             0 == entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs.count(i->key);

    if ( is_owner ) {
      // Is owner, potentially removing ghost-sends
      // Have to make a copy

      std::vector<EntityCommInfo> comm_ghost ;
      const PairIterEntityComm ec = entity_comm_map(i->key, ghosting);
      comm_ghost.assign( ec.first , ec.second );

      for ( ; ! comm_ghost.empty() ; comm_ghost.pop_back() ) {
        const EntityCommInfo tmp = comm_ghost.back();

        if ( 0 == entitiesToGhostOntoOtherProcessors.count( EntityProc( i->entity , tmp.proc ) ) ) {
          entity_comm_map_erase(i->key, tmp);
        }
      }
    }
    else if ( remove_recv ) {
        entity_comm_map_erase(i->key, ghosting);
    }

    if ( entity_comm_map(i->key).empty() ) {
      removed = true ;
      i->key = EntityKey(); // No longer communicated
      if ( remove_recv ) {
        ThrowRequireMsg( internal_destroy_entity( i->entity, remove_recv ),
                         " FAILED attempt to destroy entity: " << identifier(i->entity) );
      }
    }
  }

  if ( removed ) {
    EntityCommListInfoVector::iterator i =
      std::remove_if( m_entity_comm_list.begin() ,
                      m_entity_comm_list.end() , IsInvalid() );
    m_entity_comm_list.erase( i , m_entity_comm_list.end() );
  }

  //------------------------------------
  // Push newly ghosted entities to the receivers and update the comm list.
  // Unpacking must proceed in entity-rank order so that higher ranking
  // entities that have relations to lower ranking entities will have
  // the lower ranking entities unpacked first.  The higher and lower
  // ranking entities may be owned by different processes,
  // as such unpacking must be performed in rank order.

  ghost_entities_and_fields(ghosting, entitiesToGhostOntoOtherProcessors);
}

//----------------------------------------------------------------------

namespace {

void insert_transitive_closure( BulkData& bulk_data,
                                std::set<EntityProc,EntityLess> & new_send ,
                                const EntityProc & entry )
{
  // Do not insert if I can determine that this entity is already
  // owned or shared by the receiving processor.

  bool iDoNotOwnEntity = entry.second != bulk_data.parallel_owner_rank(entry.first);
  bool entityIsShared = bulk_data.in_shared( bulk_data.entity_key(entry.first) , entry.second );
  if (  iDoNotOwnEntity && !entityIsShared  ) {

    std::pair< std::set<EntityProc,EntityLess>::iterator , bool >
      result = new_send.insert( entry );

    if ( result.second ) {
      // A new insertion, must also insert the closure

      const unsigned erank = bulk_data.entity_rank(entry.first);

      for (EntityRank irank = stk::topology::BEGIN_RANK; irank < erank; ++irank)
      {
        Entity const *irels_j = bulk_data.begin(entry.first, irank);;
        Entity const *irels_e = bulk_data.end(entry.first, irank);
        for (; irels_j != irels_e; ++irels_j)
        {
          if (bulk_data.is_valid(*irels_j)) {
            EntityProc tmp( *irels_j , entry.second );
            insert_transitive_closure( bulk_data, new_send , tmp );
          }
        }
      }
    }
  }
}

// Fill a new send list from the receive list.
void comm_recv_to_send(
  BulkData & mesh ,
  const std::set< EntityKey > & entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs,
        std::set< EntityProc , EntityLess > & entitiesToGhostOntoOtherProcessors )
{
  const int parallel_size = mesh.parallel_size();

  CommAll all( mesh.parallel() );

  // For all entity keys in new recv, send the entity key to the owning proc
  for ( int phase = 0; phase < 2; ++phase) {
    for ( std::set<EntityKey>::const_iterator
            i = entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs.begin() ; i != entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs.end() ; ++i ) {
      Entity e = mesh.get_entity(*i); // Performance issue? Not if you're just doing full regens
      const int owner = mesh.parallel_owner_rank(e);
      const EntityKey key = mesh.entity_key(e);
      all.send_buffer( owner ).pack<EntityKey>( key );
    }
    if (phase == 0) { //allocation phase
      BABBLE_STK_PARALLEL_COMM(mesh.parallel(), "          comm_recv_to_send calling allocate_buffers");
      all.allocate_buffers( parallel_size / 4 , false /* Not symmetric */ );
    }
    else { //communication phase
      BABBLE_STK_PARALLEL_COMM(mesh.parallel(), "          comm_recv_to_send calling communicate");
      all.communicate();
    }
  }

  // Insert into entitiesToGhostOntoOtherProcessors that there is an entity on another proc
  for ( int proc_rank = 0 ; proc_rank < parallel_size ; ++proc_rank ) {
    CommBuffer & buf = all.recv_buffer(proc_rank);
    while ( buf.remaining() ) {
      EntityKey key ;
      buf.unpack<EntityKey>( key );
      EntityProc tmp( mesh.get_entity( key ) , proc_rank );
      entitiesToGhostOntoOtherProcessors.insert( tmp );
    }
  }
}

// Synchronize the send list to the receive list.

void comm_sync_send_recv(
  BulkData & mesh ,
  std::set< EntityProc , EntityLess > & entitiesToGhostOntoOtherProcessors ,
  std::set< EntityKey > & entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs )
{
  const int parallel_rank = mesh.parallel_rank();
  const int parallel_size = mesh.parallel_size();

  CommAll all( mesh.parallel() );

  // Communication sizing:

  for ( std::set< EntityProc , EntityLess >::iterator
        i = entitiesToGhostOntoOtherProcessors.begin() ; i != entitiesToGhostOntoOtherProcessors.end() ; ++i ) {
    const int owner = mesh.parallel_owner_rank(i->first);
    all.send_buffer( i->second ).skip<EntityKey>(1).skip<int>(1);
    if ( owner != parallel_rank ) {
      all.send_buffer( owner ).skip<EntityKey>(1).skip<int>(1);
    }
  }

  BABBLE_STK_PARALLEL_COMM(mesh.parallel(), "          comm_sync_send_recv calling allocate_buffers");
  all.allocate_buffers( parallel_size / 4 , false /* Not symmetric */ );

  // Loop thru all entities in entitiesToGhostOntoOtherProcessors, send the entity key to the sharing/ghosting proc
  // Also, if the owner of the entity is NOT me, also send the entity key to the owing proc

  // Communication packing (with message content comments):
  for ( std::set< EntityProc , EntityLess >::iterator
        i = entitiesToGhostOntoOtherProcessors.begin() ; i != entitiesToGhostOntoOtherProcessors.end() ; ) {
    const int owner = mesh.parallel_owner_rank(i->first);

    // Inform receiver of ghosting, the receiver does not own
    // and does not share this entity.
    // The ghost either already exists or is a to-be-done new ghost.
    // This status will be resolved on the final communication pass
    // when new ghosts are packed and sent.

    const EntityKey entity_key = mesh.entity_key(i->first);
    const int proc = i->second;

    all.send_buffer( i->second ).pack(entity_key).pack(proc);

    if ( owner != parallel_rank ) {
      // I am not the owner of this entity.
      // Inform the owner of this ghosting need.
      all.send_buffer( owner ).pack(entity_key).pack(proc);

      // Erase it from my processor's ghosting responsibility:
      // The iterator passed to the erase method will be invalidated.
      std::set< EntityProc , EntityLess >::iterator jrem = i ; ++i ;
      entitiesToGhostOntoOtherProcessors.erase( jrem );
    }
    else {
      ++i ;
    }
  }

  BABBLE_STK_PARALLEL_COMM(mesh.parallel(), "          comm_sync_send_recv calling communicate");
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
        ThrowRequireMsg(mesh.is_valid(e),
            "Unknown entity key: " <<
            MetaData::get(mesh).entity_rank_name(entity_key.rank()) <<
            "[" << entity_key.id() << "]");
        EntityProc tmp( e , proc );
        entitiesToGhostOntoOtherProcessors.insert( tmp );
      }
      else if ( mesh.is_valid(e) ) {
        //  I am the receiver for this ghost.
        //  If I already have it add it to the receive list,
        //  otherwise don't worry about it - I will receive
        //  it in the final new-ghosting communication.
        entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs.insert( mesh.entity_key(e) );
      }
    }
  }
}

void insert_upward_relations(const BulkData& bulk_data, Entity rel_entity,
                             const EntityRank rank_of_orig_entity,
                             const int my_rank,
                             const int share_proc,
                             std::vector<EntityProc>& send)
{
  EntityRank rel_entity_rank = bulk_data.entity_rank(rel_entity);
  ThrowAssert(rel_entity_rank > rank_of_orig_entity);

  // If related entity is higher rank, I own it, and it is not
  // already shared by proc, ghost it to the sharing processor.
  if ( bulk_data.parallel_owner_rank(rel_entity) == my_rank &&
       ! bulk_data.in_shared( bulk_data.entity_key(rel_entity) , share_proc ) ) {

    EntityProc entry( rel_entity , share_proc );
    send.push_back( entry );

    // There may be even higher-ranking entities that need to be ghosted, so we must recurse
    const EntityRank end_rank = static_cast<EntityRank>(bulk_data.mesh_meta_data().entity_rank_count());
    EntityVector temp_entities;
    Entity const* rels = NULL;
    int num_rels = 0;
    for (EntityRank irank = static_cast<EntityRank>(rel_entity_rank + 1); irank < end_rank; ++irank)
    {
      if (bulk_data.connectivity_map().valid(rel_entity_rank, irank)) {
        num_rels = bulk_data.num_connectivity(rel_entity, irank);
        rels     = bulk_data.begin(rel_entity, irank);
      }
      else {
        num_rels = get_connectivity(bulk_data, rel_entity, irank, temp_entities);
        rels     = &*temp_entities.begin();
      }

      for (int r = 0; r < num_rels; ++r)
      {
        Entity const rel_of_rel_entity = rels[r];
        if (bulk_data.is_valid(rel_of_rel_entity)) {
          insert_upward_relations(bulk_data, rel_of_rel_entity, rel_entity_rank, my_rank, share_proc, send);
        }
      }
    }
  }
}

} // namespace <>

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void BulkData::internal_regenerate_aura()
{
  Trace_("stk::mesh::BulkData::internal_regenerate_shared_aura");
  BABBLE_STK_PARALLEL_COMM(m_parallel_machine, "  entered internal_regenerate_shared_aura");

  require_ok_to_modify();

  std::vector<EntityProc> send ;

  const EntityRank end_rank = static_cast<EntityRank>(m_mesh_meta_data.entity_rank_count());

  // Iterate over all entities with communication info, get the sharing
  // comm info for each entity, and ensure that upwardly related
  // entities to the shared entity are ghosted on the sharing proc.
  EntityVector temp_entities;
  Entity const* rels = NULL;
  int num_rels = 0;
  for ( EntityCommListInfoVector::const_iterator
      i = comm_list().begin() ; i != comm_list().end() ; ++i ) {

    const EntityRank erank = static_cast<EntityRank>(i->key.rank());

    const PairIterEntityComm aura = entity_comm_map_shared(i->key);

    for ( size_t j = 0 ; j < aura.size() ; ++j ) {

      const int share_proc = aura[j].proc ;

      for (EntityRank k_rank = static_cast<EntityRank>(erank + 1); k_rank < end_rank; ++k_rank)
      {
        if (connectivity_map().valid(erank, k_rank)) {
          num_rels = num_connectivity(i->entity, k_rank);
          rels     = begin(i->entity, k_rank);
        }
        else {
          num_rels = get_connectivity(*this, i->entity, k_rank, temp_entities);
          rels     = &*temp_entities.begin();
        }

        for (int r = 0; r < num_rels; ++r)
        {
          if (is_valid(rels[r])) {
            insert_upward_relations(*this, rels[r], erank, m_parallel_rank, share_proc, send);
          }
        }
      }
    }
  }

  // Add new aura, remove all of the old aura.
  // The change_ghosting figures out what to actually delete and add.
  internal_change_ghosting( aura_ghosting() , send , std::vector<EntityKey>(), true /*full regen*/ );

  // BABBLE_STK_PARALLEL_COMM(m_parallel_machine, "  exiting internal_regenerate_shared_aura");

}


//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------
//EndSync
//----------------------------------------------------------------------


//----------------------------------------------------------------------

int BulkData::determine_new_owner( Entity entity ) const
{
  // We will decide the new owner by looking at all the processes sharing
  // this entity. The new owner will be the sharing process with lowest rank.

  // The local process is a candidate only if the entity is not destroyed.
  int new_owner = is_valid(entity) ? m_parallel_rank : ~0u;

  for ( PairIterEntityComm
        share = entity_comm_map_shared(entity_key(entity)); ! share.empty() ; ++share ) {
    if ( share->proc < m_parallel_size &&
         ( share->proc < new_owner || m_parallel_size <= new_owner ) ) {
      new_owner = share->proc ;
    }
  }

  return new_owner ;
}

//----------------------------------------------------------------------

namespace {

// A method for quickly finding an entity within a comm list
EntityCommListInfo find_entity(const BulkData& mesh,
                               const EntityCommListInfoVector& entities,
                               const EntityKey& key)
{
  EntityCommListInfoVector::const_iterator lb_itr = std::lower_bound(entities.begin(), entities.end(), key);
  ThrowAssertMsg(lb_itr != entities.end() && lb_itr->key == key,
                 "Cannot find id: " << key.id() << " in comm-list" );
  return *lb_itr;
}

struct EntityParallelState {
  int                 from_proc;
  EntityState         state;
  EntityCommListInfo  comm_info;
  const BulkData* mesh;

  bool operator<(const EntityParallelState& rhs) const
  { return EntityLess(*mesh)(comm_info.entity, rhs.comm_info.entity); }
};


bool pack_entity_modification( const BulkData & mesh ,
                               const bool packShared ,
                               CommAll & comm )
{
  bool flag = false;
  bool packGhosted = packShared == false;

  const EntityCommListInfoVector & entityCommList = mesh.comm_list();

  for ( EntityCommListInfoVector::const_iterator
        i = entityCommList.begin() ; i != entityCommList.end() ; ++i ) {

    Entity entity = i->entity;
    EntityState status = mesh.is_valid(entity) ? mesh.state(entity) : Deleted;

    if ( status == Modified || status == Deleted ) {

      for ( PairIterEntityComm ec = mesh.entity_comm_map(i->key); ! ec.empty() ; ++ec )
      {
        if ( ( packGhosted && ec->ghost_id > 0 ) || ( packShared && ec->ghost_id == 0 ) )
        {
          comm.send_buffer( ec->proc )
              .pack<EntityKey>( i->key )
              .pack<EntityState>( status );

          flag = true ;
        }
      }
    }
  }

  return flag ;
}

void communicate_entity_modification( const BulkData & mesh ,
                                      const bool shared ,
                                      std::vector<EntityParallelState > & data )
{
  BABBLE_STK_PARALLEL_COMM(mesh.parallel(), "      entered communicate_entity_modification");
  CommAll comm( mesh.parallel() );
  const int p_size = comm.parallel_size();

  // Sizing send buffers:
  const bool local_mod = pack_entity_modification( mesh , shared , comm );

  // Allocation of send and receive buffers:
  BABBLE_STK_PARALLEL_COMM(mesh.parallel(), "          communicate_entity_modification calling allocate_buffers")
  const bool global_mod = comm.allocate_buffers( comm.parallel_size() / 4 , false , local_mod );

  if ( global_mod )
  {
    const EntityCommListInfoVector & entityCommList = mesh.comm_list();

    // Packing send buffers:
    pack_entity_modification( mesh , shared , comm );

    BABBLE_STK_PARALLEL_COMM(mesh.parallel(), "          communicate_entity_modification calling communicate")
    comm.communicate();

    for ( int procNumber = 0 ; procNumber < p_size ; ++procNumber ) {
      CommBuffer & buf = comm.recv_buffer( procNumber );
      EntityKey key;
      EntityState state;

      while ( buf.remaining() ) {

        buf.unpack<EntityKey>( key )
           .unpack<EntityState>( state );

        // search through entity_comm, should only receive info on entities
        // that are communicated.
        EntityCommListInfo info = find_entity(mesh, entityCommList, key);
        EntityParallelState parallel_state = {procNumber, state, info, &mesh};
        data.push_back( parallel_state );
      }
    }
  }

  std::sort( data.begin() , data.end() );

  // BABBLE_STK_PARALLEL_COMM(mesh.parallel(), "     exiting communicate_entity_modification");
}

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

// Postconditions:
//  * DistributedIndex is updated based on entity creation/deletions in the
//    last modification cycle.
//  * Comm lists for shared entities are up-to-date.
//  * shared_new contains all entities that were modified/created on a
//    different process

void BulkData::internal_update_distributed_index(
        stk::mesh::EntityRank entityRank, std::vector<Entity> & shared_new )
{
  Trace_("stk::mesh::BulkData::internal_update_distributed_index");
  std::vector<EntityKey> entity_keys;
  //resolve_edge_sharing(*this, entity_keys);
  resolve_entity_sharing(*this, entityRank, entity_keys);

  shared_new.clear();
  shared_new.resize(entity_keys.size());
  for (size_t i=0;i<entity_keys.size();i++)
  {
      shared_new[i] = get_entity(entity_keys[i]);
  }
}

void BulkData::fillLocallyCreatedOrModifiedEntities(parallel::DistributedIndex::KeyTypeVector &local_created_or_modified)
{
  size_t num_created_or_modified = 0;

  for ( impl::EntityRepository::const_iterator
        i = m_entity_repo.begin() ; i != m_entity_repo.end() ; ++i )
  {
    Entity entity = i->second ;

    if ( state(entity) != Unchanged &&  isEntityMarked(entity)==BulkData::NOT_MARKED &&
         in_owned_closure( *this, entity , m_parallel_rank ) )
    {
      // Has been changed and is in owned closure, may be shared
      ++num_created_or_modified;
    }
  }

  local_created_or_modified.reserve(num_created_or_modified);

  for ( impl::EntityRepository::const_iterator
        i = m_entity_repo.begin() ; i != m_entity_repo.end() ; ++i )
  {
    Entity entity = i->second ;

    if ( state(entity) != Unchanged &&  isEntityMarked(entity)==BulkData::NOT_MARKED &&
         in_owned_closure(*this, entity , m_parallel_rank ) )
    {
      // Has been changed and is in owned closure, may be shared
      local_created_or_modified.push_back( entity_key(entity) );
    }
  }
}

void BulkData::internal_update_distributed_index(
        std::vector<Entity> & shared_new )
{
    Trace_("stk::mesh::BulkData::internal_update_distributed_index");
    BABBLE_STK_PARALLEL_COMM(m_parallel_machine, "      entered internal_update_distributed_index");

    // Scoop up a list of all nodes that have had their sharing information
    // provided directly by the end-user.  These need special handling so
    // that DistributedIndex isn't used to query sharing information.
    //
    std::vector<EntityKey> shared_nodes;
    gather_shared_nodes(*this, shared_nodes);

    std::vector<shared_entity_type> shared_edges;
    markEntitiesForResolvingSharingInfoUsingNodes(*this, stk::topology::EDGE_RANK, shared_edges);

    std::vector<shared_entity_type> shared_faces;
    markEntitiesForResolvingSharingInfoUsingNodes(*this, stk::topology::FACE_RANK, shared_faces);

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
    if (parallel_size() > 1)
    {
        // Retrieve data regarding which processes use the local_created_or_modified
        // including this process.
        parallel::DistributedIndex::KeyProcVector global_created_or_modified ;
        m_entities_index.query_to_usage( local_created_or_modified ,
                                 global_created_or_modified );

        //------------------------------
        // Take the usage data and update the sharing comm lists
        {
            Entity entity = Entity();

            // Iterate over all global modifications to this entity, this vector is
            // sorted, so we're guaranteed that all modifications to a particular
            // entities will be adjacent in this vector.

            for ( parallel::DistributedIndex::KeyProcVector::iterator
                    i =  global_created_or_modified.begin() ;
                    i != global_created_or_modified.end() ; ++i )
            {
                EntityKey key( static_cast<EntityKey::entity_key_t>(i->first) );
                int modifying_proc = i->second;

                if ( m_parallel_rank != modifying_proc )
                {
                    // Another process also created or updated this entity.
                    // Only want to look up entities at most once
                    if ( !is_valid(entity) || entity_key(entity) != key )
                    {
                        // Have not looked this entity up by key
                        entity = get_entity( key );
                        entity_keys.push_back(key);
                    }

                    // Add the other_process to the entity's sharing info.
                    entity_comm_map_insert(entity, EntityCommInfo( stk::mesh::BulkData::SHARED, modifying_proc ));
                }
            }
        }
    }

    update_shared_entities_global_ids( *this, shared_edges );
    update_shared_entities_global_ids( *this, shared_faces );

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
        if ( isEntityMarked(entity) == BulkData::IS_SHARED )
        {
            entity_keys.push_back(shared_edges[i].global_key);
        }
    }

    for (size_t i=0; i<shared_faces.size(); ++i)
    {
        Entity entity = get_entity(shared_faces[i].global_key);
        if ( isEntityMarked(entity) == BulkData::IS_SHARED )
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

//----------------------------------------------------------------------
//----------------------------------------------------------------------

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

      bool upwardRelationOfEntityIsInClosure = in_owned_closure(mesh, e , mesh.parallel_rank());
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
    bool isSharedEntity = !mesh.entity_comm_map_shared(i->key).empty();
    bool isNotInOwnedClosure = !in_owned_closure(mesh, entity , mesh.parallel_rank() );
    bool entityIsSharedButNotInClosure =  entityisValid && isSharedEntity && isNotInOwnedClosure;

    if ( entityIsSharedButNotInClosure )
    {
      destroy_dependent_ghosts( mesh , entity );
    }
  }
}

void BulkData::internal_establish_new_owner(stk::mesh::Entity entity)
{
    const int new_owner = determine_new_owner(entity);

    const bool changed = this->set_parallel_owner_rank(entity, new_owner);
    if(changed)
    {
        internal_change_owner_in_comm_data(entity_key(entity), new_owner);
    }
    set_synchronized_count(entity, m_sync_count);
}

void BulkData::internal_update_parts_for_shared_entity(stk::mesh::Entity entity, const bool is_entity_shared, const bool did_i_just_become_owner)
{
    PartVector parts_to_add_entity_to , parts_to_remove_entity_from ;

    if ( !is_entity_shared ) {
      parts_to_remove_entity_from.push_back(& m_mesh_meta_data.globally_shared_part());
    }

    if ( did_i_just_become_owner ) {
      parts_to_add_entity_to.push_back( & m_mesh_meta_data.locally_owned_part() );
    }

    if ( ! parts_to_add_entity_to.empty() || ! parts_to_remove_entity_from.empty() ) {
      internal_change_entity_parts( entity , parts_to_add_entity_to , parts_to_remove_entity_from );
    }
}

// Resolve modifications for shared entities:
// If not locally destroyed and remotely modified
// then set to locally modified.
// If remotely destroyed then determine the new owner.
//
// Post condition:
//  Shared entities are in-sync with respect to modification state.
//  Shared communication lists are updated to reflect all deletions.
//  Ownership has been re-assigned as necessary for deletion
//  of shared entities.

void BulkData::internal_resolve_shared_modify_delete()
{
  Trace_("stk::mesh::BulkData::internal_resolve_shared_modify_delete");
  BABBLE_STK_PARALLEL_COMM(m_parallel_machine, "  entered internal_resolve_shared_modify_delete");

  ThrowRequireMsg(parallel_size() > 1, "Do not call this in serial");

  delete_shared_entities_which_are_no_longer_in_owned_closure( *this );

  std::vector< EntityParallelState > remotely_modified_shared_entities ;

  // Communicate entity modification state for shared entities
  // the resulting vector is sorted by entity and process.
  const bool communicate_shared = true ;
  communicate_entity_modification( *this , communicate_shared , remotely_modified_shared_entities );

  // We iterate backwards over remote_mod to ensure that we hit the
  // higher-ranking entities first.
  for ( std::vector<EntityParallelState>::reverse_iterator
        i = remotely_modified_shared_entities.rbegin(); i != remotely_modified_shared_entities.rend() ; ) {

    Entity entity                = i->comm_info.entity;
    EntityKey key                = i->comm_info.key;
    int owner                    = i->comm_info.owner;
    const bool locally_destroyed = !is_valid(entity);
    bool remote_owner_destroyed  = false;

    // Iterate over all of this entity's remote changes
    for ( ; i != remotely_modified_shared_entities.rend() && i->comm_info.entity == entity ; ++i ) {

      const int remote_proc    = i->from_proc ;
      const bool remotely_destroyed = Deleted == i->state ;

      // When a shared entity is remotely modified or destroyed
      // then the local copy is also modified.  This modification
      // status is applied to all related higher ranking entities.

      if ( ! locally_destroyed ) {
        this->mark_entity_and_upward_related_entities_as_modified( entity );
      }

      // A shared entity is being deleted on the remote process.
      // Remove it from the sharing communication list.
      // Ownership changes are processed later, but we'll need
      // to know if the remote owner destroyed the entity in order
      // to correctly resolve ownership (it is not sufficient to just
      // look at the comm list of the entity since there is no
      // guarantee that the comm list is correct or up-to-date).

      if ( remotely_destroyed ) {
        entity_comm_map_erase( key, EntityCommInfo(SHARED,remote_proc) );

        // check if owner is destroying
        if ( owner == remote_proc ) {
          remote_owner_destroyed = true ;
        }
      }
    }

    // Have now processed all remote changes knowledge for this entity.

    if(!locally_destroyed)
    {
        const bool am_i_old_local_owner = m_parallel_rank == parallel_owner_rank(entity);

        if ( remote_owner_destroyed ) {
            internal_establish_new_owner(entity);
        }

        const bool am_i_new_local_owner = m_parallel_rank == parallel_owner_rank(entity);
        const bool did_i_just_become_owner = ( ! am_i_old_local_owner && am_i_new_local_owner );

        const bool is_entity_shared = !entity_comm_map_shared(key).empty();
        internal_update_parts_for_shared_entity(entity, is_entity_shared, did_i_just_become_owner);
    }
  } // remote mod loop

  // Erase all sharing communication lists for Destroyed entities:
  for ( EntityCommListInfoVector::const_reverse_iterator
        i = comm_list().rbegin() ; i != comm_list().rend() ; ++i) {
    if ( !is_valid(i->entity) ) {
        entity_comm_map_erase( i->key, shared_ghosting() );
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

void BulkData::internal_resolve_ghosted_modify_delete()
{
  Trace_("stk::mesh::BulkData::internal_resolve_ghosted_modify_delete");
  BABBLE_STK_PARALLEL_COMM(m_parallel_machine, "  entered internal_resolve_ghosted_modify_delete");


  ThrowRequireMsg(parallel_size() > 1, "Do not call this in serial");
  // Resolve modifications for ghosted entities:

  std::vector<EntityParallelState > remotely_modified_ghosted_entities ;

  // Communicate entity modification state for ghost entities
  const bool communicate_shared = false ;
  communicate_entity_modification( *this , communicate_shared , remotely_modified_ghosted_entities );

  const size_t ghosting_count = m_ghosting.size();

  std::vector< int > ghosting_change_flags( ghosting_count , 0 );

  // We iterate backwards over remote_mod to ensure that we hit the
  // higher-ranking entities first. This is important because higher-ranking
  // entities like element must be deleted before the nodes they have are
  // deleted.
  for ( std::vector<EntityParallelState>::reverse_iterator
        i = remotely_modified_ghosted_entities.rbegin(); i != remotely_modified_ghosted_entities.rend() ; ++i ) {
    Entity entity                 = i->comm_info.entity;
    const EntityKey key           = i->comm_info.key;
    const int      remote_proc    = i->from_proc;
    const bool     local_owner    = i->comm_info.owner == m_parallel_rank ;
    const bool remotely_destroyed = Deleted == i->state ;
    const bool locally_destroyed  = !is_valid(entity);

    if ( local_owner ) { // Sending to 'remote_proc' for ghosting

      if ( remotely_destroyed ) {

        // remove from ghost-send list

        for ( size_t j = ghosting_count ; j-- ; ) {
          if ( entity_comm_map_erase( key, EntityCommInfo( j , remote_proc ) ) ) {
            ghosting_change_flags[ j ] = true ;
          }
        }
      }

      // Remotely modified ghosts are ignored

    }
    else { // Receiving from 'remote_proc' for ghosting

      // Owner modified or destroyed, must locally destroy.

        // Manoj: could be error here. I think the loop below
        // should loop thru all ghostings, m_ghosting. But effect might be same.
        // would only impact custom ghostings.
      for ( PairIterEntityComm ec = entity_comm_map(key); !ec.empty() ; ++ec ) {
        ghosting_change_flags[ ec->ghost_id ] = true ;
      }

      // This is a receive ghost so the only communication information
      // is the ghosting information, can clear it all out.
      entity_comm_map_clear(key);

      if ( ! locally_destroyed && ! in_owned_closure(*this, entity , m_parallel_rank ) )
      {
          // If mesh modification causes a ghost entity to become
          // a member of an owned-closure then do not automatically
          // destroy it.  The new sharing status will be resolved
          // in 'internal_resolve_parallel_create'.
          const bool destroy_entity_successful = internal_destroy_entity(entity, true);
          ThrowRequireMsg(destroy_entity_successful,
                          "Could not destroy ghost entity " << identifier(entity));
      }
    }
  } // end loop on remote mod

  // Erase all ghosting communication lists for:
  // 1) Destroyed entities.
  // 2) Owned and modified entities.

  for ( EntityCommListInfoVector::const_reverse_iterator
        i = comm_list().rbegin() ; i != comm_list().rend() ; ++i) {

    Entity entity = i->entity;

    const bool locally_destroyed = !is_valid(entity);
    const bool locally_owned_and_modified = locally_destroyed ? false :
      Modified == state(entity) &&
      m_parallel_rank   == i->owner ;

    if ( locally_destroyed || locally_owned_and_modified ) {

      // m_ghosting[0] is the SHARED communication

      for ( size_t j = ghosting_count ; --j ; ) {
        if ( entity_comm_map_erase( i->key, *m_ghosting[j] ) ) {
          ghosting_change_flags[ j ] = true ;
        }
      }
    }
  }

  std::vector< int > ghosting_change_flags_global( ghosting_count , 0 );

  BABBLE_STK_PARALLEL_COMM(m_parallel_machine, "      internal_resolve_ghosted_modify_delete calling all_reduce_sum");
  all_reduce_sum( m_parallel_machine ,
                  & ghosting_change_flags[0] ,
                  & ghosting_change_flags_global[0] ,
                  ghosting_change_flags.size() );

  for ( unsigned ic = 0 ; ic < ghosting_change_flags_global.size() ; ++ic ) {
    if ( ghosting_change_flags_global[ic] ) {
      m_ghosting[ic]->m_sync_count = m_sync_count ;
    }
  }
  // BABBLE_STK_PARALLEL_COMM(m_parallel_machine, "  exiting internal_resolve_ghosted_modify_delete");
}

void BulkData::resolve_ownership_of_modified_entities( const std::vector<Entity> &shared_modified )
{
    CommAll comm_all( m_parallel_machine );

    for ( int phase = 0; phase < 2; ++phase ) {
        for ( std::vector<Entity>::const_iterator i = shared_modified.begin() ; i != shared_modified.end() ; ++i ) {
            Entity entity = *i ;
            if ( parallel_owner_rank(entity) == m_parallel_rank &&
                   state(entity)  != Created ) {
                for ( PairIterEntityComm jc = entity_comm_map_shared(entity_key(entity)) ; ! jc.empty() ; ++jc ) {
                    comm_all.send_buffer( jc->proc ) .pack<EntityKey>( entity_key(entity) );
                }
            }
        }

        if (phase == 0) { //allocation phase
            BABBLE_STK_PARALLEL_COMM(m_parallel_machine, "      calling allocate_buffers from internal_resolve_parallel_create");
            comm_all.allocate_buffers( m_parallel_size / 4 );
        }
        else { // communication phase
            BABBLE_STK_PARALLEL_COMM(m_parallel_machine, "      calling communicate from internal_resolve_parallel_create");
            comm_all.communicate();
        }
    }

    for ( int p = 0 ; p < m_parallel_size ; ++p ) {
        CommBuffer & buf = comm_all.recv_buffer( p );
        EntityKey key ;
        while ( buf.remaining() ) {
            buf.unpack<EntityKey>( key );
            Entity entity = get_entity( key );

            // Set owner, will correct part membership later
            const bool changed = this->set_parallel_owner_rank( entity, p);
            if (changed) {
                internal_change_owner_in_comm_data(key, p);
            }
        }
    }
}

void BulkData::move_entities_to_proper_part_ownership( const std::vector<Entity> &shared_modified )
{
    std::ostringstream error_msg ;
    int error_flag = 0 ;

    PartVector shared_part , owned_part ;
    shared_part.push_back( & m_mesh_meta_data.globally_shared_part() );
    owned_part.push_back(  & m_mesh_meta_data.locally_owned_part() );

    std::vector<Entity>::const_reverse_iterator iend = shared_modified.rend();
    for ( std::vector<Entity>::const_reverse_iterator
          i = shared_modified.rbegin() ; i != iend ; ++i)
    {
      Entity entity = *i ;

      if ( parallel_owner_rank(entity) == m_parallel_rank &&
           state(entity) == Created )
      {
        // Created and not claimed by an existing owner

        const int new_owner = determine_new_owner( entity );

        const bool changed = this->set_parallel_owner_rank( entity, new_owner);
        if (changed) {
          internal_change_owner_in_comm_data(entity_key(entity), new_owner);
        }
      }

      if ( parallel_owner_rank(entity) != m_parallel_rank ) {
        // Do not own it and still have it.
        // Remove the locally owned, add the globally_shared
        set_synchronized_count( entity, m_sync_count);
        internal_change_entity_parts( entity , shared_part /*add*/, owned_part /*remove*/);
      }
      else if ( ! entity_comm_map_shared(entity_key(entity)).empty() ) {
        // Own it and has sharing information.
        // Add the globally_shared
        internal_change_entity_parts( entity , shared_part /*add*/, PartVector() /*remove*/ );
      }
      else {
        // Own it and does not have sharing information.
        // Remove the globally_shared
        internal_change_entity_parts( entity , PartVector() /*add*/, shared_part /*remove*/);
      }

      // Newly created shared entity had better be in the owned closure
      if ( ! in_owned_closure(*this, entity , m_parallel_rank ) ) {
        if ( 0 == error_flag ) {
          error_flag = 1 ;
          error_msg
            << "\nP" << m_parallel_rank << ": " << " FAILED\n"
            << "  The following entities were declared on multiple processors,\n"
            << "  cannot be parallel-shared, and were declared with"
            << "  parallel-ghosting information. {\n";
        }
        error_msg << "    " << print_entity_key(m_mesh_meta_data, entity_key(entity));
        error_msg << " also declared on" ;
        for ( PairIterEntityComm ec = entity_comm_map_shared(entity_key(entity)); ! ec.empty() ; ++ec ) {
          error_msg << " P" << ec->proc ;
        }
        error_msg << "\n" ;
      }
    }

    // Parallel-consistent error checking of above loop
    if ( error_flag ) { error_msg << "}\n" ; }
    all_reduce( m_parallel_machine , ReduceMax<1>( & error_flag ) );
    ThrowErrorMsgIf( error_flag, error_msg.str() );
}


void BulkData::update_comm_list(const std::vector<stk::mesh::Entity>& shared_modified)
{
    // ------------------------------------------------------------
    // Update m_entity_comm based on shared_modified

    const size_t n_old = m_entity_comm_list.size();

    m_entity_comm_list.reserve(m_entity_comm_list.size() + shared_modified.size());
    for (size_t i = 0, e = shared_modified.size(); i < e; ++i) {
      Entity entity = shared_modified[i];
      EntityCommListInfo new_comm = {entity_key(entity), entity, parallel_owner_rank(entity), NULL};
      m_entity_comm_list.push_back(new_comm);
    }

    std::inplace_merge( m_entity_comm_list.begin() ,
                        m_entity_comm_list.begin() + n_old ,
                        m_entity_comm_list.end() );

    EntityCommListInfoVector::iterator iter =
      std::unique( m_entity_comm_list.begin() , m_entity_comm_list.end() );

    m_entity_comm_list.erase( iter , m_entity_comm_list.end() );

    for(size_t i=0; i<m_entity_comm_list.size(); ++i) {
      EntityKey key = m_entity_comm_list[i].key;
      const EntityComm* entity_comm = m_entity_comm_map.entity_comm(key);
      m_entity_comm_list[i].entity_comm = entity_comm;
    }

    internal_sync_comm_list_owners();
}
//----------------------------------------------------------------------

// Postconditions:
//  * All shared entities have parallel-consistent owner
//  * Part membership of shared entities is up-to-date
//  * m_entity_comm is up-to-date
void BulkData::internal_resolve_parallel_create()
{
  Trace_("stk::mesh::BulkData::internal_resolve_parallel_create");
  BABBLE_STK_PARALLEL_COMM(m_parallel_machine, "  entered internal_resolve_parallel_create");

  ThrowRequireMsg(parallel_size() > 1, "Do not call this in serial");
  std::vector<Entity> shared_modified ;

  // Update the parallel index and
  // output shared and modified entities.
  internal_update_distributed_index(shared_modified );

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

  update_comm_list( shared_modified );
}
//----------------------------------------------------------------------

namespace {

#ifdef STK_VERBOSE_OUTPUT

bool no_buckets(const stk::mesh::BulkData& mesh, const stk::mesh::Part& part)
{
  for(stk::topology::rank_t r = stk::topology::NODE_RANK; r < mesh.mesh_meta_data().entity_rank_count(); ++r) {
    stk::mesh::Selector selector = part;
    const stk::mesh::BucketVector& buckets = mesh.get_buckets(r, selector);
    if (buckets.size() > 0) {
      return false;
    }
  }

  return true;
}

void print_bucket_data(const stk::mesh::BulkData& mesh)
{
  const stk::mesh::PartVector& all_parts = mesh.mesh_meta_data().get_parts();
  for(size_t i=0; i<all_parts.size(); ++i) {
    std::cout << "Part: " << all_parts[i]->name()<<std::endl;
    if (no_buckets(mesh, *all_parts[i])) {
      std::cout<<"\tEmpty"<<std::endl;
      continue;
    }
    for(stk::topology::rank_t r = stk::topology::NODE_RANK; r < mesh.mesh_meta_data().entity_rank_count(); ++r) {
      stk::mesh::Selector selector = *all_parts[i];
      const stk::mesh::BucketVector& buckets = mesh.get_buckets(r, selector);
      std::cout<<"\t"<< buckets.size() << " "<< r << " buckets";
      size_t min_entities = 1000000, max_entities = 0, total_entities = 0;
      double avg_entities = 0.0;
      for(size_t j=0; j<buckets.size(); ++j) {
        total_entities += buckets[j]->size();
        min_entities = std::min(min_entities, buckets[j]->size());
        max_entities = std::max(max_entities, buckets[j]->size());
      }
      avg_entities = buckets.size()>0 ? (1.0*total_entities)/buckets.size() : 0.0;
      if (total_entities == 0) {
        min_entities = 0;
      }
      std::cout << "; min=" << min_entities << ", avg=" << avg_entities << ", max=" << max_entities << ", tot=" << total_entities << std::endl;
    }
  }
}

#endif

}

//----------------------------------------------------------------------

bool BulkData::modification_end( modification_optimization opt)
{
  Trace_("stk::mesh::BulkData::modification_end");

  bool return_value = internal_modification_end( true, opt );

#ifdef GATHER_GET_BUCKETS_METRICS
  ++m_num_modifications;
#endif

#ifdef STK_VERBOSE_OUTPUT
  print_bucket_data(*this);
#endif

  m_modification_begin_description = "NO_MODIFICATION";

  write_modification_counts();

  return return_value;
}

bool BulkData::modification_end_for_entity_creation( EntityRank entity_rank, modification_optimization opt)
{
  Trace_("stk::mesh::BulkData::modification_end");

  bool return_value = internal_modification_end_for_entity_creation( entity_rank, true, opt );

#ifdef GATHER_GET_BUCKETS_METRICS
  ++m_num_modifications;
#endif

#ifdef STK_VERBOSE_OUTPUT
  print_bucket_data(*this);
#endif

  return return_value;
}

#if 0

namespace {

// Very, very handy for debugging parallel resolution...

void print_comm_list( const BulkData & mesh , bool doit )
{
  if ( doit ) {
    std::ostringstream msg ;

    msg << std::endl ;

    for ( EntityCommListInfoVector::const_iterator
          i =  mesh.comm_list().begin() ;
          i != mesh.comm_list().end() ; ++i ) {

      Entity entity = i->entity;
      msg << "P" << mesh.parallel_rank() << ": " ;

      print_entity_key( msg , MetaData::get(mesh) , i->key );

      msg << " owner(" << i->owner << ")" ;

      if ( !entity.is_valid() ) { msg << " del" ; }
      else if ( Modified == entity.state() ) { msg << " mod" ; }
      else { msg << "    " ; }

      for ( PairIterEntityComm ec = mesh.comm_list(i->key); ! ec.empty() ; ++ec ) {
        msg << " (" << ec->ghost_id << "," << ec->proc << ")" ;
      }
      msg << std::endl ;
    }

    std::cout << msg.str();
  }
}

}

#endif

namespace {

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
            std::cerr << "Entity with rank="<<rank<<", identifier="<<mesh.identifier(bucket[j])<<" is connected to an invalid node."<<std::endl;
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

}

void BulkData::update_comm_list_based_on_changes_in_comm_map()
// Resolution of shared and ghost modifications can empty
// the communication information for entities.
// If there is no communication information then the
// entity must be removed from the communication list.
{
  EntityCommListInfoVector::iterator i = m_entity_comm_list.begin();
  bool changed = false ;
  for ( ; i != m_entity_comm_list.end() ; ++i ) {
    const EntityComm* entity_comm = m_entity_comm_map.entity_comm(i->key);
    if ( entity_comm == NULL || entity_comm->comm_map.empty() ) {
      i->key = EntityKey();
      changed = true;
    }
    i->entity_comm = entity_comm;
  }
  if ( changed ) {
    i = std::remove_if( m_entity_comm_list.begin() ,
                        m_entity_comm_list.end() , IsInvalid() );
    m_entity_comm_list.erase( i , m_entity_comm_list.end() );
  }
}

bool BulkData::internal_modification_end( bool regenerate_aura, modification_optimization opt )
{
  Trace_("stk::mesh::BulkData::internal_modification_end");

  // The two states are MODIFIABLE and SYNCHRONiZED
  if ( m_sync_state == SYNCHRONIZED ) { return false ; }

  ThrowAssertMsg(check_for_connected_nodes(*this)==0, "BulkData::modification_end ERROR, all entities with rank higher than node are required to have connected nodes.");

  ThrowAssertMsg(add_fmwk_data() || check_no_shared_elements_or_higher(*this)==0, "BulkData::modification_end ERROR, Sharing of entities with rank ELEMENT_RANK or higher is not allowed.");

  if (parallel_size() > 1) {
    // Resolve modification or deletion of shared entities
    // which can cause deletion of ghost entities.
    internal_resolve_shared_modify_delete();

    // Resolve modification or deletion of ghost entities
    // by destroying ghost entities that have been touched.
    internal_resolve_ghosted_modify_delete();

    update_comm_list_based_on_changes_in_comm_map();

    // Resolve creation of entities: discover sharing and set unique ownership.
    internal_resolve_parallel_create();

    // Resolve part membership for shared entities.
    // This occurs after resolving creation so created and shared
    // entities are resolved along with previously existing shared entities.
    internal_resolve_shared_membership();

    // Regenerate the ghosting aura around all shared mesh entities.
    if ( regenerate_aura ) { internal_regenerate_aura(); }

    // ------------------------------
    // Verify parallel consistency of mesh entities.
    // Unique ownership, communication lists, sharing part membership,
    // application part membership consistency.
#ifndef NDEBUG
    std::ostringstream msg ;
    bool is_consistent = true;
    is_consistent = comm_mesh_verify_parallel_consistency( *this , msg );
    ThrowErrorMsgIf( !is_consistent, msg.str() );
#endif
  }
  else {
      std::vector<Entity> shared_modified ;
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
    m_bucket_repository.optimize_buckets();
  }
  else {
    m_bucket_repository.internal_sort_bucket_entities();
  }

  // ------------------------------

  m_bucket_repository.internal_modification_end();

  internal_update_fast_comm_maps();

  m_sync_state = SYNCHRONIZED ;

  update_deleted_entities_container();

  return true ;
}

//////////////////////////////////// Free funcions to help with modification end (exp) for edges

void fillElementsConnectedToNodes(stk::mesh::BulkData &stkMeshBulkData, std::vector<stk::mesh::Entity> &nodes,
                                  std::vector<stk::mesh::Entity> & elementsConnectedToNodes)
{
    stk::mesh::Entity const * elemStartNode = stkMeshBulkData.begin_elements(nodes[0]);
    stk::mesh::Entity const * elemEndNode = stkMeshBulkData.end_elements(nodes[0]);
    elementsConnectedToNodes.assign(elemStartNode, elemEndNode);
    std::sort(elementsConnectedToNodes.begin(), elementsConnectedToNodes.end());
    std::vector<stk::mesh::Entity>::iterator intersectionEnd = elementsConnectedToNodes.end();

    std::vector<stk::mesh::Entity> elementsConnectedToNodesTemp;
    for (size_t i = 1; i < nodes.size(); ++i)
    {
      elementsConnectedToNodesTemp.clear();
      elementsConnectedToNodesTemp.assign(elementsConnectedToNodes.begin(), intersectionEnd);
      elemStartNode = stkMeshBulkData.begin_elements(nodes[i]);
      elemEndNode = stkMeshBulkData.end_elements(nodes[i]);
      std::vector<stk::mesh::Entity> elems(elemStartNode, elemEndNode);
      std::sort(elems.begin(), elems.end());

      intersectionEnd = std::set_intersection( elems.begin(), elems.end(),
                                               elementsConnectedToNodesTemp.begin(), elementsConnectedToNodesTemp.end(),
                                               elementsConnectedToNodes.begin() );
      if (intersectionEnd == elementsConnectedToNodes.begin()) break;  // Empty set
    }

    elementsConnectedToNodes.resize(intersectionEnd-elementsConnectedToNodes.begin());
}

void fillFacesConnectedToNodes(stk::mesh::BulkData &stkMeshBulkData, std::vector<stk::mesh::Entity> &nodes,
              std::vector<stk::mesh::Entity> & facesConnectedToNodes)
{
    facesConnectedToNodes.clear();
    stk::mesh::EntityRank needConnectivityOfType = stk::topology::FACE_RANK;
    if (stkMeshBulkData.connectivity_map().valid(stk::topology::NODE_RANK, needConnectivityOfType))
    {
        //TODO: fill faces for all edge nodes
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
    else
    {
        stk::mesh::EntityVector entitiesNode1;
        stk::mesh::EntityVector entitiesNode2;
        int numFacesConnectedToNode1 = stk::mesh::get_connectivity(stkMeshBulkData, nodes[0], needConnectivityOfType, entitiesNode1);
        int numFacesConnectedToNode2 = stk::mesh::get_connectivity(stkMeshBulkData, nodes[1], needConnectivityOfType, entitiesNode2);
        if ( numFacesConnectedToNode1 > 0 && numFacesConnectedToNode2 > 0 )
        {
            facesConnectedToNodes.resize(entitiesNode1.size()+entitiesNode2.size());

            std::vector<stk::mesh::Entity>::iterator iter = std::set_intersection( entitiesNode1.begin(), entitiesNode1.end(),
                    entitiesNode2.begin(), entitiesNode2.end(), facesConnectedToNodes.begin());

            facesConnectedToNodes.resize(iter-facesConnectedToNodes.begin());
        }
    }
}

bool doesEdgeNeedGhostingCommunication(stk::mesh::BulkData &stkMeshBulkData, std::vector<stk::mesh::Entity>& connectedEntities)
{
    bool communicate_edge_for_ghosting = false;
    for (size_t j=0;j<connectedEntities.size();j++)
    {
        bool isEntityGhostedToAnotherProc = stkMeshBulkData.is_ghosted_onto_another_proc(stkMeshBulkData.entity_key(connectedEntities[j]));
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


void connectEntityToEdge(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Entity entity,
        stk::mesh::Entity edge, std::vector<stk::mesh::Entity> &nodes)
{
    // get node entity ids
    std::vector<stk::mesh::EntityId> nodeIdsForEdge(2);
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

void connectUpwardEntityToEntity(stk::mesh::BulkData& mesh, stk::mesh::Entity upward_entity,
        stk::mesh::Entity entity, std::vector<stk::mesh::Entity> &nodes)
{
    // get node upwardentity ids
    uint num_nodes = mesh.num_nodes(entity);
    EntityRank entity_rank = mesh.entity_rank(entity);
    std::vector<stk::mesh::EntityId> nodeIdsForEntity(num_nodes);
    for (uint i = 0 ; i < num_nodes ; ++i)
    {
      nodeIdsForEntity[i]= mesh.identifier(nodes[i]);
    }

    // scratch space
    stk::mesh::OrdinalVector ordinal_scratch;
    ordinal_scratch.reserve(64);
    stk::mesh::PartVector part_scratch;
    part_scratch.reserve(64);
    stk::mesh::Permutation perm = static_cast<stk::mesh::Permutation>(0);

    // now what
    stk::topology upward_entity_topology = mesh.bucket(upward_entity).topology();
    std::vector<stk::mesh::EntityId> nodeIds(num_nodes);
    std::vector<stk::mesh::EntityId> upward_entity_node_ids(upward_entity_topology.num_nodes());
    unsigned entity_ordinal = 100000;
    stk::mesh::Entity const * upward_entity_nodes = mesh.begin_nodes(upward_entity);
    for (size_t k=0;k<upward_entity_topology.num_nodes();k++)
    {
        upward_entity_node_ids[k] = mesh.identifier(upward_entity_nodes[k]);
    }

    stk::topology entity_top;
    for (size_t k=0;k<upward_entity_topology.num_sub_topology(entity_rank);k++)
    {
        if(entity_rank == stk::topology::EDGE_RANK)
        {
          upward_entity_topology.edge_nodes(upward_entity_node_ids, k, nodeIds.begin());
          entity_top = upward_entity_topology.edge_topology();
        }
        else
        {
          upward_entity_topology.face_nodes(upward_entity_node_ids, k, nodeIds.begin());
          entity_top = upward_entity_topology.face_topology();
        }
        if ( entity_top.equivalent(nodeIds, nodeIdsForEntity).first )
        {
            entity_ordinal = k;
            break;
        }
    }
    ThrowRequireMsg(entity_ordinal !=100000, "Program error. Contact sierra-help for support.");
    mesh.declare_relation(upward_entity, entity, entity_ordinal, perm, ordinal_scratch, part_scratch);
}

void connectGhostedEntitiesToEdge(stk::mesh::BulkData &stkMeshBulkData, std::vector<stk::mesh::Entity> &entitiesConnectedToNodes, stk::mesh::Entity edge, std::vector<stk::mesh::Entity>& nodes)
{
    for (size_t j=0; j<entitiesConnectedToNodes.size();j++)
    {
        bool isEntityGhostedOntoThisProc = stkMeshBulkData.in_receive_ghost(stkMeshBulkData.aura_ghosting(), stkMeshBulkData.entity_key(entitiesConnectedToNodes[j]));

        if ( isEntityGhostedOntoThisProc )
        {
            connectEntityToEdge(stkMeshBulkData, entitiesConnectedToNodes[j], edge, nodes);
        }
    }
}

void connectGhostedEntitiesToEntity(stk::mesh::BulkData &stkMeshBulkData, std::vector<stk::mesh::Entity> &entitiesConnectedToNodes, stk::mesh::Entity entity, std::vector<stk::mesh::Entity>& nodes)
{
    for (size_t j=0; j<entitiesConnectedToNodes.size();j++)
    {
        bool isEntityGhostedOntoThisProc = stkMeshBulkData.in_receive_ghost(stkMeshBulkData.aura_ghosting(), stkMeshBulkData.entity_key(entitiesConnectedToNodes[j]));

        if ( isEntityGhostedOntoThisProc )
        {
            connectUpwardEntityToEntity(stkMeshBulkData, entitiesConnectedToNodes[j], entity, nodes);
        }
    }
}

void determineEntitiesThatNeedGhosting(stk::mesh::BulkData &stkMeshBulkData, stk::mesh::Entity edge, std::vector<stk::mesh::Entity>& entitiesConnectedToNodes,
        std::vector<stk::mesh::Entity> &nodes, std::set<EntityProc, EntityLess> &addGhostedEntities)
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
                PairIterEntityComm ghosted = stkMeshBulkData.entity_comm_map( stkMeshBulkData.entity_key(entitiesConnectedToNodes[j]) , stkMeshBulkData.aura_ghosting());
                for (PairIterEntityComm ec = ghosted; !ec.empty(); ++ec)
                {
                    if ( ec->proc != stkMeshBulkData.parallel_rank() )
                    {
                        bool isEdgeSharedWithOtherProc = stkMeshBulkData.in_shared(stkMeshBulkData.entity_key(edge), ec->proc);
                        bool isEntityGhostedToAnotherProc = stkMeshBulkData.is_ghosted_onto_proc(stkMeshBulkData.entity_key(entitiesConnectedToNodes[j]), ec->proc);

                        if ( !isEdgeSharedWithOtherProc && isEntityGhostedToAnotherProc )
                        {
                            addGhostedEntities.insert(EntityProc(edge, ec->proc));
                        }
                    }
                }
            }
        }
    }
}

void connect_edges_to_ghosted_entities_and_fill_edges_to_ghost_onto_other_processors(stk::mesh::BulkData &stkMeshBulkata, std::set<EntityProc, EntityLess> &entitiesToGhostOntoOtherProcessors)
{
    const stk::mesh::BucketVector& edge_buckets = stkMeshBulkata.buckets(stk::topology::EDGE_RANK);

    for(size_t bucketIndex = 0; bucketIndex < edge_buckets.size(); bucketIndex++)
    {
        const stk::mesh::Bucket& bucket = *edge_buckets[bucketIndex];
        for(size_t entityIndex = 0; entityIndex < bucket.size(); entityIndex++)
        {
            Entity edge = bucket[entityIndex];
            if ( stkMeshBulkata.state(edge) == Unchanged ) continue;

            std::vector<stk::mesh::Entity> edge_nodes(2);
            stk::mesh::Entity const *  nodes = stkMeshBulkata.begin_nodes(edge);
            edge_nodes[0] = nodes[0];
            edge_nodes[1] = nodes[1];

            std::vector<stk::mesh::Entity> facesConnectedToNodes;
            fillFacesConnectedToNodes(stkMeshBulkata, edge_nodes, facesConnectedToNodes);
            connectGhostedEntitiesToEdge(stkMeshBulkata, facesConnectedToNodes, edge, edge_nodes);

            std::vector<stk::mesh::Entity> elementsConnectedToNodes;
            fillElementsConnectedToNodes(stkMeshBulkata, edge_nodes, elementsConnectedToNodes);
            connectGhostedEntitiesToEdge(stkMeshBulkata, elementsConnectedToNodes, edge, edge_nodes);

            if ( bucket.owned() || bucket.shared() )
            {
                determineEntitiesThatNeedGhosting(stkMeshBulkata, edge, facesConnectedToNodes, edge_nodes, entitiesToGhostOntoOtherProcessors);
                determineEntitiesThatNeedGhosting(stkMeshBulkata, edge, elementsConnectedToNodes, edge_nodes, entitiesToGhostOntoOtherProcessors);
            }
        }
    }

    std::set< EntityKey > entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs;
    comm_sync_send_recv(stkMeshBulkata, entitiesToGhostOntoOtherProcessors, entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs);
}

void find_upward_connected_entities_to_ghost_onto_other_processors(stk::mesh::BulkData &mesh, std::set<EntityProc, EntityLess> &entitiesToGhostOntoOtherProcessors, EntityRank entity_rank)
{
    const stk::mesh::BucketVector& entity_buckets = mesh.buckets(entity_rank);
    bool isedge = (entity_rank == stk::topology::EDGE_RANK);

    for(size_t bucketIndex = 0; bucketIndex < entity_buckets.size(); bucketIndex++)
    {
        const stk::mesh::Bucket& bucket = *entity_buckets[bucketIndex];
        size_t num_nodes_on_entity = bucket.topology().num_nodes();
        std::vector<stk::mesh::Entity> entity_nodes(num_nodes_on_entity);

        for(size_t entityIndex = 0; entityIndex < bucket.size(); entityIndex++)
        {
            Entity entity = bucket[entityIndex];
            if ( mesh.state(entity) == Unchanged ) continue;

            entity_nodes.assign(mesh.begin_nodes(entity), mesh.end_nodes(entity));
            std::vector<stk::mesh::Entity> facesConnectedToNodes;

            if(isedge)
            {
              fillFacesConnectedToNodes(mesh, entity_nodes, facesConnectedToNodes);
              //connectGhostedEntitiesToEdge(mesh, facesConnectedToNodes, entity, entity_nodes);
              connectGhostedEntitiesToEntity(mesh, facesConnectedToNodes, entity, entity_nodes);
            }

            std::vector<stk::mesh::Entity> elementsConnectedToNodes;
            fillElementsConnectedToNodes(mesh, entity_nodes, elementsConnectedToNodes);
            //connectGhostedEntitiesToEdge(mesh, elementsConnectedToNodes, entity, entity_nodes);
            connectGhostedEntitiesToEntity(mesh, elementsConnectedToNodes, entity, entity_nodes);

            if ( bucket.owned() || bucket.shared() )
            {
                if (isedge)
                {
                  determineEntitiesThatNeedGhosting(mesh, entity, facesConnectedToNodes, entity_nodes, entitiesToGhostOntoOtherProcessors);
                }
                determineEntitiesThatNeedGhosting(mesh, entity, elementsConnectedToNodes, entity_nodes, entitiesToGhostOntoOtherProcessors);
            }
        }
    }

    std::set< EntityKey > entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs;
    comm_sync_send_recv(mesh, entitiesToGhostOntoOtherProcessors, entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs);
}

void connect_ghosted_edges_received_to_ghosted_faces_and_elements(stk::mesh::BulkData &stkMeshBulkData)
{
    if (!stkMeshBulkData.connectivity_map().valid(stk::topology::EDGE_RANK, stk::topology::FACE_RANK))
    {
        const stk::mesh::BucketVector& edge_buckets = stkMeshBulkData.buckets(stk::topology::EDGE_RANK);
        for(size_t bucketIndex = 0; bucketIndex < edge_buckets.size(); bucketIndex++)
        {
            const stk::mesh::Bucket& bucket = *edge_buckets[bucketIndex];
            if ( bucket.in_aura() )
            {
                for(size_t entityIndex = 0; entityIndex < bucket.size(); entityIndex++)
                {
                    Entity edge = bucket[entityIndex];
                    if ( stkMeshBulkData.state(edge) == Unchanged ) continue;

                    std::vector<stk::mesh::Entity> edge_nodes(2);
                    stk::mesh::Entity const *  nodes = stkMeshBulkData.begin_nodes(edge);
                    edge_nodes[0] = nodes[0];
                    edge_nodes[1] = nodes[1];

                    std::vector<stk::mesh::Entity> facesConnectedToNodes;
                    fillFacesConnectedToNodes(stkMeshBulkData, edge_nodes, facesConnectedToNodes);
                    connectGhostedEntitiesToEdge(stkMeshBulkData, facesConnectedToNodes, edge, edge_nodes);
                }
            }
        }
    }

    if (!stkMeshBulkData.connectivity_map().valid(stk::topology::EDGE_RANK, stk::topology::ELEMENT_RANK) )
    {
        const stk::mesh::BucketVector& edge_buckets = stkMeshBulkData.buckets(stk::topology::EDGE_RANK);
        for(size_t bucketIndex = 0; bucketIndex < edge_buckets.size(); bucketIndex++)
        {
            const stk::mesh::Bucket& bucket = *edge_buckets[bucketIndex];
            if ( bucket.in_aura() )
            {
                for(size_t entityIndex = 0; entityIndex < bucket.size(); entityIndex++)
                {
                    Entity edge = bucket[entityIndex];
                    if ( stkMeshBulkData.state(edge) == Unchanged ) continue;

                    std::vector<stk::mesh::Entity> edge_nodes(2);
                    stk::mesh::Entity const *  nodes = stkMeshBulkData.begin_nodes(edge);
                    edge_nodes[0] = nodes[0];
                    edge_nodes[1] = nodes[1];

                    std::vector<stk::mesh::Entity> elementsConnectedToNodes;
                    fillElementsConnectedToNodes(stkMeshBulkData, edge_nodes, elementsConnectedToNodes);
                    connectGhostedEntitiesToEdge(stkMeshBulkData, elementsConnectedToNodes, edge, edge_nodes);
                }
            }
        }
    }
}

void connect_ghosted_entities_received_to_ghosted_upwardly_connected_entities(stk::mesh::BulkData &mesh, EntityRank entity_rank)
{
    if (entity_rank == stk::topology::EDGE_RANK && !mesh.connectivity_map().valid(stk::topology::EDGE_RANK, stk::topology::FACE_RANK) )
    {
        const stk::mesh::BucketVector& entity_buckets = mesh.buckets(stk::topology::EDGE_RANK);
        for(size_t bucketIndex = 0; bucketIndex < entity_buckets.size(); bucketIndex++)
        {
            const stk::mesh::Bucket& bucket = *entity_buckets[bucketIndex];
            if ( bucket.in_aura() )
            {
                for(size_t entityIndex = 0; entityIndex < bucket.size(); entityIndex++)
                {
                    Entity entity = bucket[entityIndex];
                    if ( mesh.state(entity) == Unchanged ) continue;

                    std::vector<stk::mesh::Entity> entity_nodes(2);
                    stk::mesh::Entity const *  nodes = mesh.begin_nodes(entity);
                    entity_nodes[0] = nodes[0];
                    entity_nodes[1] = nodes[1];

                    std::vector<stk::mesh::Entity> facesConnectedToNodes;
                    fillFacesConnectedToNodes(mesh, entity_nodes, facesConnectedToNodes);
                    connectGhostedEntitiesToEdge(mesh, facesConnectedToNodes, entity, entity_nodes);
                }
            }
        }
    }

    if (!mesh.connectivity_map().valid(entity_rank, stk::topology::ELEMENT_RANK) )
    {
        const stk::mesh::BucketVector& entity_buckets = mesh.buckets(entity_rank);
        for(size_t bucketIndex = 0; bucketIndex < entity_buckets.size(); bucketIndex++)
        {
            const stk::mesh::Bucket& bucket = *entity_buckets[bucketIndex];
            if ( bucket.in_aura() )
            {
              size_t num_nodes_on_entity = bucket.topology().num_nodes();
              std::vector<stk::mesh::Entity> entity_nodes(num_nodes_on_entity);

                for(size_t entityIndex = 0; entityIndex < bucket.size(); entityIndex++)
                {
                    Entity entity = bucket[entityIndex];
                    if ( mesh.state(entity) == Unchanged ) continue;

                    entity_nodes.assign(mesh.begin_nodes(entity), mesh.end_nodes(entity));

                    std::vector<stk::mesh::Entity> elementsConnectedToNodes;
                    fillElementsConnectedToNodes(mesh, entity_nodes, elementsConnectedToNodes);
                    connectGhostedEntitiesToEdge(mesh, elementsConnectedToNodes, entity, entity_nodes);
                }
            }
        }
    }
}

bool BulkData::internal_modification_end_for_entity_creation( EntityRank entity_rank, bool regenerate_aura, modification_optimization opt )
{
  Trace_("stk::mesh::BulkData::internal_modification_end_for_entity_creation");

  // The two states are MODIFIABLE and SYNCHRONiZED
  if ( m_sync_state == SYNCHRONIZED ) { return false ; }

  ThrowAssertMsg(check_for_connected_nodes(*this)==0, "BulkData::modification_end ERROR, all entities with rank higher than node are required to have connected nodes.");

  if (parallel_size() > 1)
  {
    std::vector<Entity> shared_modified ;

    // Update the parallel index and
    // output shared and modified entities.
    internal_update_distributed_index( entity_rank, shared_modified );

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

    update_comm_list( shared_modified );

    internal_resolve_shared_membership();

    std::set<EntityProc, EntityLess> entitiesToGhostOntoOtherProcessors(EntityLess(*this));

    //connect_edges_to_ghosted_entities_and_fill_edges_to_ghost_onto_other_processors(*this, entitiesToGhostOntoOtherProcessors);
    find_upward_connected_entities_to_ghost_onto_other_processors(*this, entitiesToGhostOntoOtherProcessors, entity_rank);

    ghost_entities_and_fields(aura_ghosting(), entitiesToGhostOntoOtherProcessors);

    //connect_ghosted_edges_received_to_ghosted_faces_and_elements(*this);
    connect_ghosted_entities_received_to_ghosted_upwardly_connected_entities(*this, entity_rank);

#ifndef NDEBUG
    std::ostringstream msg ;
    bool is_consistent = true;
    is_consistent = comm_mesh_verify_parallel_consistency( *this , msg );
    ThrowErrorMsgIf( !is_consistent, msg.str() );
#endif
  }
  else {
      std::vector<Entity> shared_modified ;
      internal_update_distributed_index( entity_rank, shared_modified );
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
    m_bucket_repository.optimize_buckets();
  }
  else {
    m_bucket_repository.internal_sort_bucket_entities();
  }

  // ------------------------------

  m_bucket_repository.internal_modification_end();

  internal_update_fast_comm_maps();

  m_sync_state = SYNCHRONIZED ;

  update_deleted_entities_container();

  return true ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

template <typename T>
T const* get_begin_itr(const Bucket& bucket, Bucket::size_type bucket_ordinal, EntityRank rank);

template <typename T>
T const* get_end_itr(const Bucket& bucket, Bucket::size_type bucket_ordinal, EntityRank rank);


template <>
Entity const* get_begin_itr<Entity>(const Bucket& bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.begin(bucket_ordinal, rank); }

template <>
Entity const* get_end_itr<Entity>(const Bucket& bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.end(bucket_ordinal, rank); }


template <>
ConnectivityOrdinal const* get_begin_itr<ConnectivityOrdinal>(const Bucket& bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.begin_ordinals(bucket_ordinal, rank); }

template <>
ConnectivityOrdinal const* get_end_itr<ConnectivityOrdinal>(const Bucket& bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.end_ordinals(bucket_ordinal, rank); }


template <>
Permutation const* get_begin_itr<Permutation>(const Bucket& bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.begin_permutations(bucket_ordinal, rank); }

template <>
Permutation const* get_end_itr<Permutation>(const Bucket& bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.end_permutations(bucket_ordinal, rank); }

//
// Because the connectivity API in Bucket is not templated...
//

template <typename T>
T const *get_begin_relation_data(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank);

template <typename T>
T const *get_end_relation_data(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank);

template <>
Entity const *get_begin_relation_data<Entity>(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.begin(bucket_ordinal, rank); }

template <>
Entity const *get_end_relation_data<Entity>(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.end(bucket_ordinal, rank); }

template <>
ConnectivityOrdinal const *get_begin_relation_data<ConnectivityOrdinal>(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.begin_ordinals(bucket_ordinal, rank); }

template <>
ConnectivityOrdinal const *get_end_relation_data<ConnectivityOrdinal>(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.end_ordinals(bucket_ordinal, rank); }

template <>
Permutation const *get_begin_relation_data<Permutation>(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.begin_permutations(bucket_ordinal, rank); }

template <>
Permutation const *get_end_relation_data<Permutation>(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.end_permutations(bucket_ordinal, rank); }

template <typename T>
void verify_relation_data(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{
  T const * old_it = get_begin_relation_data<T>(bucket, bucket_ordinal, rank);
  T const * old_it_end = get_end_relation_data<T>(bucket, bucket_ordinal, rank);

  T const * new_data_begin = get_begin_itr<T>(bucket, bucket_ordinal, rank);
  T const * new_data_end   = get_end_itr<T>(bucket, bucket_ordinal, rank);

  ThrowRequire( std::distance(new_data_begin, new_data_end ) ==
                std::distance(old_it, old_it_end ) );

  ThrowRequire( std::distance(new_data_begin, new_data_end) ==
                bucket.num_connectivity(bucket_ordinal, rank) );

  T const * new_it = new_data_begin;
  for (; old_it != old_it_end ; ++old_it , ++new_it) {
    T old_data = *old_it;
    T new_data = *new_it;
    ThrowRequire(old_data == new_data);
  }
}

}

void BulkData::verify_relations(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank) const
{
  verify_relation_data<Entity>(bucket, bucket_ordinal, rank);
  verify_relation_data<ConnectivityOrdinal>(bucket, bucket_ordinal, rank);
  if (bucket.has_permutation(rank)) {
    verify_relation_data<Permutation>(bucket, bucket_ordinal, rank);
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

enum { PART_ORD_UNIVERSAL = 0 };
enum { PART_ORD_OWNED     = 1 };
enum { PART_ORD_SHARED    = 2 };

namespace {

void pack_induced_memberships( BulkData& bulk_data,
                               CommAll & comm ,
                               const EntityCommListInfoVector & entity_comm )
{
  OrdinalVector empty , induced ;
  for ( EntityCommListInfoVector::const_iterator
        i = entity_comm.begin() ; i != entity_comm.end() ; ++i ) {

    if ( bulk_data.in_shared( i->key , i->owner ) ) {
      // Is shared with owner, send to owner.

      empty.clear();
      induced.clear();

      induced_part_membership(bulk_data, i->entity , empty , induced );

      CommBuffer & buf = comm.send_buffer( i->owner );

      unsigned tmp = induced.size();

      buf.pack<unsigned>( tmp );

      for ( OrdinalVector::iterator
            j = induced.begin() ; j != induced.end() ; ++j ) {
        buf.pack<unsigned>( *j );
      }
    }
  }
}

void generate_send_list( BulkData& bulk_data,
                         const size_t sync_count ,
                         const int p_rank ,
                         const EntityCommListInfoVector & entity_comm ,
                               std::vector<EntityProc> & send_list )
{
  for ( EntityCommListInfoVector::const_iterator
        i = entity_comm.begin() ; i != entity_comm.end() ; ++i ) {

    if ( i->owner == p_rank &&
         bulk_data.synchronized_count(i->entity) == sync_count ) {

      for ( PairIterEntityComm ec = bulk_data.entity_comm_map(i->key); ! ec.empty(); ++ec ) {
        EntityProc tmp( i->entity , ec->proc );
        send_list.push_back( tmp );
      }
    }
  }

  {
    std::sort( send_list.begin() , send_list.end() , EntityLess(bulk_data) );
    std::vector<EntityProc>::iterator i =
      std::unique( send_list.begin() , send_list.end() );
    send_list.erase( i , send_list.end() );
  }
}

void pack_part_memberships( BulkData& meshbulk, CommAll & comm ,
                            const std::vector<EntityProc> & send_list )
{
  for ( std::vector<EntityProc>::const_iterator
        i = send_list.begin() ; i != send_list.end() ; ++i ) {

    Entity entity = i->first;

    std::pair<const unsigned *, const unsigned *>
      part_ord = meshbulk.bucket(entity).superset_part_ordinals();

    // I am the owner; therefore, the first three members are
    // universal, uses, and owns.  Don't send them.

    // I am the owner.  The first two memberships are
    // universal_part and locally_owned_part.  The third
    // membership may be globally_shared_part ;

    const unsigned count_all  = part_ord.second - part_ord.first ;
    const unsigned count_skip =
      ( 2 < count_all && part_ord.first[2] == PART_ORD_SHARED ) ? 3 : 2 ;

    const unsigned count_send = count_all - count_skip ;

    const unsigned * const start_send = part_ord.first + count_skip ;

    comm.send_buffer( i->second ).pack<EntityKey>( meshbulk.entity_key(entity) )
                                 .pack<unsigned>( count_send )
                                 .pack<unsigned>( start_send , count_send );
  }
}

}

//  Mesh entity membership changes must be synchronized among
//  processes that share mesh entities and propagated to
//  processes that ghost copies of the mesh entities.
//
//  Precondition: correct shared and ghosting lists.
//
//  Part memberships may have been added or removed
//  either explicitly or indirectly via entity relationships
//  being added or removed.

void BulkData::internal_resolve_shared_membership()
{
    Trace_("stk::mesh::BulkData::internal_resolve_shared_membership");
    BABBLE_STK_PARALLEL_COMM(m_parallel_machine, "  entered internal_resolve_shared_membership");

    ThrowRequireMsg(parallel_size() > 1, "Do not call this in serial");

    const MetaData & meta = m_mesh_meta_data;
    ParallelMachine p_comm = m_parallel_machine;
    const int p_rank = m_parallel_rank;
    const int p_size = m_parallel_size;
    const PartVector & all_parts = meta.get_parts();

    const Part & part_universal = meta.universal_part();
    const Part & part_owned = meta.locally_owned_part();
    const Part & part_shared = meta.globally_shared_part();

    // Quick verification of part ordinal assumptions

    ThrowRequireMsg(PART_ORD_UNIVERSAL == part_universal.mesh_meta_data_ordinal(),
            "Universal part ordinal is wrong, expected "
            << PART_ORD_UNIVERSAL << ", got: "
            << part_universal.mesh_meta_data_ordinal());

    ThrowRequireMsg(PART_ORD_OWNED == part_owned.mesh_meta_data_ordinal(),
            "Owned part ordinal is wrong, expected "
            << PART_ORD_OWNED << ", got: "
            << part_owned.mesh_meta_data_ordinal());

    ThrowRequireMsg(PART_ORD_SHARED == part_shared.mesh_meta_data_ordinal(),
            "Shared part ordinal is wrong, expected "
            << PART_ORD_SHARED << ", got: "
            << part_shared.mesh_meta_data_ordinal());

    //  Shared entities may have been modified due to relationship changes.
    //  Send just the current induced memberships from the sharing to
    //  the owning processes.
    {
        CommAll comm(p_comm);

        pack_induced_memberships(*this, comm, m_entity_comm_list);

        BABBLE_STK_PARALLEL_COMM(p_comm, "          internal_resolve_shared_membership calling allocate_buffers");
        comm.allocate_buffers(p_size / 4);

        pack_induced_memberships(*this, comm, m_entity_comm_list);

        BABBLE_STK_PARALLEL_COMM(p_comm, "          internal_resolve_shared_membership calling communicate");
        comm.communicate();

        OrdinalVector empty, induced_parts, current_parts, remove_parts;
        PartVector inducedParts, removeParts;

        for(EntityCommListInfoVector::iterator i = m_entity_comm_list.begin(); i != m_entity_comm_list.end(); ++i)
        {
            bool i_own_this_entity_in_comm_list = i->owner == p_rank;
            if( i_own_this_entity_in_comm_list )
            {
                // Receiving from all sharing processes

                empty.clear();
                induced_parts.clear();
                current_parts.clear();
                remove_parts.clear();

                induced_part_membership(*this, i->entity, empty, induced_parts);

                for(PairIterEntityComm ec = entity_comm_map_shared(i->key); !ec.empty(); ++ec)
                {

                    CommBuffer & buf = comm.recv_buffer(ec->proc);

                    unsigned count = 0;
                    buf.unpack<unsigned>(count);
                    for(unsigned j = 0; j < count; ++j)
                    {
                        unsigned part_ord = 0;
                        buf.unpack<unsigned>(part_ord);
                        insert_ordinal(induced_parts, part_ord);
                    }
                }

                // Remove any part that is an induced part but is not
                // in the induced parts list.

                this->bucket(i->entity).supersets(current_parts);

                OrdinalVector::const_iterator induced_parts_begin = induced_parts.begin(),
                        induced_parts_end = induced_parts.end();

                for(OrdinalVector::iterator
                p = current_parts.begin(); p != current_parts.end(); ++p)
                {
                    if(meta.get_parts()[*p]->was_induced(i->key.rank()) &&
                            !contains_ordinal(induced_parts_begin, induced_parts_end, *p))
                    {
                        remove_parts.push_back(*p);
                    }
                }

                inducedParts.clear();
                removeParts.clear();

                inducedParts.reserve(induced_parts.size());
                for(unsigned ipart = 0; ipart < induced_parts.size(); ++ipart)
                {
                    inducedParts.push_back(&m_mesh_meta_data.get_part(induced_parts[ipart]));
                }
                removeParts.reserve(remove_parts.size());
                for(unsigned ipart = 0; ipart < remove_parts.size(); ++ipart)
                {
                    removeParts.push_back(&m_mesh_meta_data.get_part(remove_parts[ipart]));
                }

                internal_change_entity_parts(i->entity, inducedParts, removeParts);
            }
        }
    }

    //------------------------------
    // The owners have complete knowledge of memberships.
    // Send membership information to sync the shared and ghosted copies.
    // Only need to do this for entities that have actually changed.

    {
        std::vector<EntityProc> send_list;

        generate_send_list(*this, m_sync_count, p_rank, m_entity_comm_list, send_list);

        CommAll comm(p_comm);

        pack_part_memberships(*this, comm, send_list);

        BABBLE_STK_PARALLEL_COMM(p_comm, "          internal_resolve_shared_membership calling allocate_buffers");
        comm.allocate_buffers(p_size / 4);

        pack_part_memberships(*this, comm, send_list);

        BABBLE_STK_PARALLEL_COMM(p_comm, "          internal_resolve_shared_membership calling communicate");
        comm.communicate();

        for(int p = 0; p < p_size; ++p)
        {
            CommBuffer & buf = comm.recv_buffer(p);
            while(buf.remaining())
            {

                PartVector owner_parts, current_parts, remove_parts;

                EntityKey key;
                buf.unpack<EntityKey>(key);
                unsigned count = 0;
                buf.unpack<unsigned>(count);
                for(unsigned j = 0; j < count; ++j)
                {
                    unsigned part_ord = 0;
                    buf.unpack<unsigned>(part_ord);
                    insert(owner_parts, *all_parts[part_ord]);
                }

                // Any current part that is not a member of owners_parts
                // must be removed.

                Entity const entity = find_entity(*this, m_entity_comm_list, key).entity;

                this->bucket(entity).supersets(current_parts);

                for(PartVector::iterator
                ip = current_parts.begin(); ip != current_parts.end(); ++ip)
                {
                    Part * const part = *ip;
                    const unsigned part_ord = part->mesh_meta_data_ordinal();
                    if(PART_ORD_UNIVERSAL != part_ord &&
                            PART_ORD_OWNED != part_ord &&
                            PART_ORD_SHARED != part_ord &&
                            !contain(m_ghost_parts, *part) &&
                            !contain(owner_parts, *part))
                    {
                        remove_parts.push_back(part);
                    }
                }

                internal_change_entity_parts(entity, owner_parts, remove_parts);
            }
        }
    }
}

void BulkData::internal_update_fast_comm_maps()
{
  if (parallel_size() > 1) {
    EntityCommListInfoVector const& all_comm = comm_list();

    // Flush previous map
    const EntityRank num_ranks = static_cast<EntityRank>(m_mesh_meta_data.entity_rank_count());
    m_volatile_fast_shared_comm_map.resize(num_ranks);
    for (EntityRank r = stk::topology::BEGIN_RANK; r < num_ranks; ++r) {
      m_volatile_fast_shared_comm_map[r].resize(parallel_size());
      for (int proc = 0; proc < parallel_size(); ++proc) {
        m_volatile_fast_shared_comm_map[r][proc].clear();
      }
    }

    // Assemble map, find all shared entities and pack into volatile fast map
    for (size_t i = 0, ie = all_comm.size(); i < ie; ++i) {
      Entity const e        = all_comm[i].entity;
      EntityKey const key   = all_comm[i].key;
      MeshIndex const& idx  = mesh_index(e);
      EntityRank const rank = key.rank();

      FastMeshIndex fast_idx;
      fast_idx.bucket_id  = idx.bucket->bucket_id();
      fast_idx.bucket_ord = idx.bucket_ordinal;

      PairIterEntityComm ec = entity_comm_map(key);
      for (; !ec.empty() && ec->ghost_id == 0; ++ec) {
        m_volatile_fast_shared_comm_map[rank][ec->proc].push_back(fast_idx);
      }
    }

    // Need to shrink-to-fit these vectors?
  }
}



namespace impl {

unsigned get_ordinal(const Part* part)
{ return part->mesh_meta_data_ordinal(); }

unsigned get_ordinal(unsigned ord)
{ return ord; }

const Part& get_part(const Part* part, MetaData& meta)
{ return *part; }

const Part& get_part(unsigned ord, MetaData& meta)
{ return *meta.get_parts()[ord]; }

void filter_out( std::vector<unsigned> & vec ,
                 const std::vector<Part*> & parts ,
                 std::vector<Part*> & removed )
{
  std::vector<unsigned>::iterator i , j ;
  i = j = vec.begin();

  std::vector<Part*>::const_iterator ip = parts.begin() ;

  while ( j != vec.end() && ip != parts.end() ) {
    if      ( get_ordinal(*ip) < *j ) { ++ip ; }
    else if ( *j < get_ordinal(*ip) ) { *i = *j ; ++i ; ++j ; }
    else {
      removed.push_back( *ip );
      ++j ;
      ++ip ;
    }
  }

  if ( i != j ) { vec.erase( i , j ); }
}

void merge_in( std::vector<unsigned> & vec , const std::vector<Part*> & parts )
{
  std::vector<unsigned>::iterator i = vec.begin();
  std::vector<Part*>::const_iterator ip = parts.begin() ;

  for ( ; i != vec.end() && ip != parts.end() ; ++i ) {

    const unsigned ord = get_ordinal(*ip);

    if ( ord <= *i ) {
      if ( ord < *i ) { i = vec.insert( i , ord ); }
      // Now have: ord == *i
      ++ip ;
    }
  }

  for ( ; ip != parts.end() ; ++ip ) {
    const unsigned ord = get_ordinal(*ip);
    vec.push_back( ord );
  }
}

} // namespace impl

void BulkData::change_entity_parts( Entity entity,
    const PartVector & add_parts ,
    const PartVector & remove_parts,
    bool always_propagate_internal_changes)
{
    INCREMENT_ENTITY_MODIFICATION_COUNTER(PUBLIC, entity_rank(entity), DESTROY_ALL_GHOSTING);
    internal_verify_and_change_entity_parts(entity, add_parts, remove_parts, always_propagate_internal_changes);
}

void BulkData::internal_verify_and_change_entity_parts( Entity entity,
                                                        const PartVector & add_parts ,
                                                        const PartVector & remove_parts,
                                                        bool always_propagate_internal_changes)
{
  TraceIfWatching("stk::mesh::BulkData::internal_change_entity_parts", LOG_ENTITY, entity_key(entity));
  DiagIfWatching(LOG_ENTITY, entity_key(entity), "entity state: " << entity_key(entity));

  require_ok_to_modify();

// When stk parallel is used within Fmwk, this assertion (require_entity_owner) is violated
// So, temporarily, don't test this assertion if SIERRA_MIGRATION is defined, and the bulk
// data point is set.  (Any other use case will go ahead and test this assertion.)
#ifdef SIERRA_MIGRATION
  if (NULL == get_fmwk_bulk_data())
    require_entity_owner( entity , m_parallel_rank );
#else
  require_entity_owner( entity , m_parallel_rank );
#endif

  const EntityRank ent_rank = entity_rank(entity);

  const EntityRank undef_rank  = InvalidEntityRank;

  // Transitive addition and removal:
  // 1) Include supersets of add_parts
  // 2) Do not include a remove_part if it appears in the add_parts
  // 3) Include subsets of remove_parts

  // most parts will at least have universal and topology part as supersets
  const unsigned expected_min_num_supersets = 2;

  PartVector a_parts;
  a_parts.reserve( std::distance(add_parts.begin(), add_parts.end()) * (expected_min_num_supersets + 1) );
  for(PartVector::const_iterator add_iter=add_parts.begin(); add_iter!=add_parts.end(); ++add_iter) {
#ifdef FMWK_NO_GLOBALLY_SHARED_ELEMENTS
    ThrowErrorMsgIf(entity_rank == stk::topology::ELEMENT_RANK && **add_iter == mesh_meta_data().globally_shared_part(), "FMWK_NO_GLOBALLY_SHARED_ELEMENTS  Error in BulkData::internal_change_entity_parts, trying to make an element globally shared!");
#endif // FMWK_NO_GLOBALLY_SHARED_ELEMENTS
    a_parts.push_back((*add_iter));
  }
  bool quick_verify_check = true;

  for ( PartVector::const_iterator ia = add_parts.begin(); ia != add_parts.end() ; ++ia ) {
    quick_verify_check = quick_verify_check &&
      internal_quick_verify_change_part(*ia, ent_rank, undef_rank);
    const PartVector& supersets = (*ia)->supersets();
    for(PartVector::const_iterator s_iter=supersets.begin(), s_end=supersets.end();
        s_iter!=s_end; ++s_iter) {
      a_parts.push_back((*s_iter));
    }
  }

  order(a_parts);

  PartVector::const_iterator a_parts_begin = a_parts.begin(),
                                a_parts_end   = a_parts.end();
  PartVector r_parts ;

  for ( PartVector::const_iterator ir = remove_parts.begin(); ir != remove_parts.end() ; ++ir ) {

    // The following guards should be in the public interface to
    // changing parts.  However, internal mechanisms such as changing
    // ownership calls this function to add or remove an entity from
    // the three special parts.  Without refactoring, these guards
    // cannot be put in place.
    /*
    ThrowErrorMsgIf( m_mesh_meta_data.universal_part() == **ir,
                     "Cannot remove entity from universal part" );
    ThrowErrorMsgIf( m_mesh_meta_data.locally_owned_part() == **ir,
                     "Cannot remove entity from locally owned part" );
    ThrowErrorMsgIf( m_mesh_meta_data.globally_shared_part() == **ir,
                     "Cannot remove entity from globally shared part" );
    */

    quick_verify_check = quick_verify_check &&
      internal_quick_verify_change_part(*ir, ent_rank, undef_rank);

    if ( ! contains_ordinal_part( a_parts_begin, a_parts_end , (*ir)->mesh_meta_data_ordinal() ) ) {
      r_parts.push_back( (*ir) );
      for ( PartVector::const_iterator  cur_part = (*ir)->subsets().begin() ;
            cur_part != (*ir)->subsets().end() ;
            ++cur_part )
        if ( bucket(entity).member ( **cur_part ) )
          r_parts.push_back ( (*cur_part) );
    }
  }

  order(r_parts);

  // If it looks like we have a problem, run the full check and we should
  // expect to see an exception thrown; otherwise, only do the full check in
  // debug mode because it incurs significant overhead.
  if ( ! quick_verify_check ) {
    internal_verify_change_parts( m_mesh_meta_data , entity , a_parts );
    internal_verify_change_parts( m_mesh_meta_data , entity , r_parts );
    ThrowRequireMsg(false, "Expected throw from verify methods above.");
  }
  else {
#ifndef NDEBUG
    internal_verify_change_parts( m_mesh_meta_data , entity , a_parts );
    internal_verify_change_parts( m_mesh_meta_data , entity , r_parts );
#endif
  }

  internal_change_entity_parts( entity , a_parts , r_parts , always_propagate_internal_changes );
}

//  The 'add_parts' and 'remove_parts' are complete and disjoint.
//  Changes need to have parallel resolution during
//  modification_end.

void BulkData::internal_change_entity_parts(
  Entity entity ,
  const std::vector<Part*> & add_parts ,
  const std::vector<Part*> & remove_parts,
  bool always_propagate_internal_changes )
{
  TraceIfWatching("stk::mesh::BulkData::internal_change_entity_parts", LOG_ENTITY, entity_key(entity));
  DiagIfWatching(LOG_ENTITY, entity_key(entity), "entity state: " << entity_key(entity));
  DiagIfWatching(LOG_ENTITY, entity_key(entity), "add_parts: " << add_parts);
  DiagIfWatching(LOG_ENTITY, entity_key(entity), "remove_parts: " << remove_parts);

  Bucket * const bucket_old = bucket_ptr( entity );

  if ( bucket_old && bucket_old->member_all( add_parts ) &&
              ! bucket_old->member_any( remove_parts ) ) {
    // Is already a member of all add_parts,
    // is not a member of any remove_parts,
    // thus nothing to do.
    return ;
  }

  INCREMENT_ENTITY_MODIFICATION_COUNTER(INTERNAL, entity_rank(entity), CHANGE_ENTITY_PARTS);

  const unsigned locally_owned_ordinal = m_mesh_meta_data.locally_owned_part().mesh_meta_data_ordinal();

  bool add_to_locally_owned = false;
  for (std::vector<Part*>::const_iterator itr = add_parts.begin(), end_itr = add_parts.end(); itr != end_itr; ++itr) {
    if ( impl::get_ordinal(*itr) == locally_owned_ordinal ) {
      add_to_locally_owned = true;
      break;
    }
  }
  add_to_locally_owned = add_to_locally_owned && (!bucket_old || !bucket_old->owned());


  bool remove_from_locally_owned = false;
  for (std::vector<Part*>::const_iterator itr = remove_parts.begin(), end_itr = remove_parts.end(); itr != end_itr; ++itr) {
    if ( impl::get_ordinal(*itr) == locally_owned_ordinal ) {
      remove_from_locally_owned = true;
      break;
    }
  }
  remove_from_locally_owned = remove_from_locally_owned && (!bucket_old || bucket_old->owned());

  if (add_to_locally_owned) {

    ++m_closure_count[entity.local_offset()];

    // update downward connectivity closure count
    if (bucket_old) {
      for (EntityRank rank = stk::topology::NODE_RANK, end_rank = bucket_old->entity_rank(); rank < end_rank; ++rank) {
        unsigned num = num_connectivity(entity,rank);
        Entity const * entities = begin(entity,rank);
        for (unsigned i =0; i<num; ++i) {
          ++m_closure_count[entities[i].local_offset()];
        }
      }
    }

  }
  else if (remove_from_locally_owned)
  {
    --m_closure_count[entity.local_offset()];

    // update downward connectivity closure count
    if (bucket_old) {
      for (EntityRank rank = stk::topology::NODE_RANK, end_rank = bucket_old->entity_rank(); rank < end_rank; ++rank) {
        unsigned num = num_connectivity(entity,rank);
        Entity const * entities = begin(entity,rank);
        for (unsigned i =0; i<num; ++i) {
          --m_closure_count[entities[i].local_offset()];
        }
      }
    }
  }

  std::vector<Part*> parts_removed ;

  OrdinalVector parts_total ; // The final part list

  //--------------------------------



  if ( bucket_old ) {
    // Keep any of the existing bucket's parts
    // that are not a remove part.
    // This will include the 'intersection' parts.
    //
    // These parts are properly ordered and unique.

    const std::pair<const unsigned *, const unsigned*>
      bucket_parts = bucket_old->superset_part_ordinals();

    const unsigned * parts_begin = bucket_parts.first;
    const unsigned * parts_end   = bucket_parts.second;

    const unsigned num_bucket_parts = parts_end - parts_begin;
    parts_total.reserve( num_bucket_parts + add_parts.size() );
    parts_total.insert( parts_total.begin(), parts_begin , parts_end);

    if ( !remove_parts.empty() ) {
      parts_removed.reserve(remove_parts.size());
      impl::filter_out( parts_total , remove_parts , parts_removed );
    }
  }
  else {
    parts_total.reserve(add_parts.size());
  }

  if ( !add_parts.empty() ) {
    impl::merge_in( parts_total , add_parts );
  }

  if ( parts_total.empty() ) {
    // Always a member of the universal part.
    const unsigned univ_ord =
      m_mesh_meta_data.universal_part().mesh_meta_data_ordinal();
    parts_total.push_back( univ_ord );
  }

  EntityRank e_rank = entity_rank(entity);

  //--------------------------------
  // Move the entity to the new partition.
  stk::mesh::impl::Partition *partition =
          m_bucket_repository.get_or_create_partition(e_rank, parts_total);

  if (bucket_old) {
    bucket_old->getPartition()->move_to(entity, *partition);
  }
  else {
    partition->add(entity);
  }

  ////
  //// SHOULD WE FIND A WAY TO PUT THE REST OF THIS IN Partition::move_to(..)?
  ////

  // Propagate part changes through the entity's relations.
  //(Only propagate part changes for parts which have a primary-entity-rank that matches
  // the entity's rank. Other parts don't get induced...)

  std::vector<Part*> inducable_parts_removed;
  for(std::vector<Part*>::const_iterator pr=parts_removed.begin(), prend=parts_removed.end(); pr!=prend; ++pr) {
    Part const& check_part = impl::get_part(*pr, m_mesh_meta_data);
    if (check_part.should_induce(e_rank)) {
      inducable_parts_removed.push_back(*pr);
    }
  }

  if (always_propagate_internal_changes ||
      !inducable_parts_removed.empty() ) {
    internal_propagate_part_changes( entity , inducable_parts_removed );
  }
}


//----------------------------------------------------------------------
// Deduce propagation of part membership changes to a 'from' entity
// to the related 'to' entities.  There can be both additions and
// removals.

void BulkData::internal_propagate_part_changes(
  Entity entity ,
  const std::vector<Part*> & removed )
{
  TraceIfWatching("stk::mesh::BulkData::internal_propagate_part_changes",
      LOG_ENTITY,
      entity_key(entity));
  DiagIfWatching(LOG_ENTITY, entity_key(entity), "entity state: " << entity_key(entity));
  DiagIfWatching(LOG_ENTITY, entity_key(entity), "Removed: " << removed);

  m_check_invalid_rels = false;

  const EntityRank erank = entity_rank(entity);
  const EntityRank end_rank = static_cast<EntityRank>(m_mesh_meta_data.entity_rank_count());

  OrdinalVector to_del , to_add , empty ;
  EntityVector temp_entities;
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

        to_del.clear();
        to_add.clear();
        empty.clear();

        // Induce part membership from this relationship to
        // pick up any additions.
        induced_part_membership(*this, entity, empty, irank, to_add );

        if ( ! removed.empty() ) {
          // Something was removed from the 'from' entity,
          // deduce what may have to be removed from the 'to' entity.

          // Deduce parts for 'e_to' from all upward relations.
          // Any non-parallel part that I removed that is not deduced for
          // 'e_to' must be removed from 'e_to'

          EntityRank e_to_rank = entity_rank(e_to);

          Entity const* back_rel_entities = NULL;
          int num_back_rels = 0;
          for (EntityRank to_rel_rank_i = static_cast<EntityRank>(e_to_rank + 1); to_rel_rank_i < end_rank; ++to_rel_rank_i)
          {
            if (connectivity_map().valid(e_to_rank, to_rel_rank_i)) {
              num_back_rels     = num_connectivity(e_to, to_rel_rank_i);
              back_rel_entities = begin(e_to, to_rel_rank_i);
            }
            else {
              num_back_rels = get_connectivity(*this, e_to, to_rel_rank_i, temp_entities);
              back_rel_entities = &*temp_entities.begin();
            }

            for (int k = 0; k < num_back_rels; ++k)
            {
              if (entity != back_rel_entities[k])  // Already did this entity
              {
                // Relation from to_rel->entity() to e_to
                induced_part_membership(*this, back_rel_entities[k], empty, e_to_rank, to_add );
              }
            }
          }

          OrdinalVector::const_iterator to_add_begin = to_add.begin(),
            to_add_end   = to_add.end();

          for (std::vector<Part*>::const_iterator
                  k = removed.begin() ; k != removed.end() ; ++k ) {
            if ( ! contains_ordinal( to_add_begin, to_add_end , impl::get_ordinal(*k) ) ) {
              induced_part_membership( impl::get_part(*k, m_mesh_meta_data), erank, irank, to_del );
            }
          }
        }


        PartVector addParts, delParts, emptyParts;

        addParts.reserve(to_add.size());
        for(unsigned ipart=0; ipart<to_add.size(); ++ipart) {
          addParts.push_back(&m_mesh_meta_data.get_part(to_add[ipart]));
        }
        delParts.reserve(to_del.size());
        for(unsigned ipart=0; ipart<to_del.size(); ++ipart) {
          delParts.push_back(&m_mesh_meta_data.get_part(to_del[ipart]));
        }



        if ( parallel_size() < 2 || entity_comm_map_shared(entity_key(e_to)).empty() ) {
          // Entirely local, ok to remove memberships now
          internal_change_entity_parts( e_to , addParts , delParts );
        }
        else {
          // Shared, do not remove memberships now.
          // Wait until modification_end.
          internal_change_entity_parts( e_to , addParts , emptyParts );
        }
      }
    }
  }
  m_check_invalid_rels = true;
}

// TODO Change the methods below to requirements (private, const invariant checkers)

// Do not allow any of the induced part memberships to explicitly
// appear in the add or remove parts lists.
// 1) Intersection part
// 3) Part that does not match the entity rank.

void BulkData::internal_verify_change_parts( const MetaData   & meta ,
                                             const Entity entity ,
                                             const std::vector<Part*> & parts ) const
{
  const std::vector<std::string> & rank_names = meta.entity_rank_names();
  const EntityRank undef_rank  = InvalidEntityRank;
  const EntityRank erank = entity_rank(entity);

  bool ok = true ;
  std::ostringstream msg ;

  for (std::vector<Part*>::const_iterator
        i = parts.begin() ; i != parts.end() ; ++i ) {

    const Part & p = impl::get_part(*i, m_mesh_meta_data);
    const unsigned part_rank = p.primary_entity_rank();

    bool intersection_ok=false, rel_target_ok=false, rank_ok=false;
    internal_basic_part_check(&p, erank, undef_rank, intersection_ok, rel_target_ok, rank_ok);

    if ( !intersection_ok || !rel_target_ok || !rank_ok ) {
      if ( ok ) {
        ok = false ;
        msg << "change parts for entity " << identifier(entity);
        msg << " , { " ;
      }
      else {
        msg << " , " ;
      }

      msg << p.name() << "[" ;
      if ( part_rank < rank_names.size() ) {
        msg << rank_names[ part_rank ];
      }
      else {
        msg << part_rank ;
      }
      msg << "] " ;
      if ( !intersection_ok ) { msg << "is_intersection " ; }
      if ( !rel_target_ok )   { msg << "is_relation_target " ; }
      if ( !rank_ok )         { msg << "is_bad_rank " ; }
    }
  }

  ThrowErrorMsgIf( !ok, msg.str() << "}" );
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


} // namespace mesh
} // namespace stk
