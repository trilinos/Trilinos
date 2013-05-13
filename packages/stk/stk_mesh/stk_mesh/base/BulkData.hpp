/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_BulkData_hpp
#define stk_mesh_BulkData_hpp

//----------------------------------------------------------------------

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/DistributedIndex.hpp>
#include <stk_util/util/TrackingAllocator.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Ghosting.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Trace.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityCommDatabase.hpp>
#include <stk_mesh/base/ConnectivityMap.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

#include <stk_mesh/baseImpl/EntityRepository.hpp>
#include <stk_mesh/baseImpl/BucketRepository.hpp>

#include <algorithm>
#include <map>
#include <limits>

#include <list>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

namespace impl {

class Partition;
class EntityRepository;

template <typename T>
struct IndexToRank{};

template <>
struct IndexToRank<Node>
{
  static const EntityRank rank = stk::topology::NODE_RANK;
};

template <>
struct IndexToRank<Edge>
{
  static const EntityRank rank = stk::topology::EDGE_RANK;
};

template <>
struct IndexToRank<Face>
{
  static const EntityRank rank = stk::topology::FACE_RANK;
};

template <>
struct IndexToRank<Element>
{
  static const EntityRank rank = stk::topology::ELEMENT_RANK;
};

}

struct EntityCommListInfo
{
  EntityKey key;
  Entity    entity; // Might be invalid if entity has been deleted
  int  owner;
};

inline
bool operator<(const EntityKey& key, const EntityCommListInfo& comm)
{ return key < comm.key; }

inline
bool operator<(const EntityCommListInfo& comm, const EntityKey& key)
{ return comm.key < key; }

inline
bool operator<(const EntityCommListInfo& lhs, const EntityCommListInfo& rhs)
{ return lhs.key < rhs.key; }

inline
bool operator==(const EntityCommListInfo& lhs, const EntityCommListInfo& rhs)
{ return lhs.key == rhs.key; }

inline
bool operator!=(const EntityCommListInfo& lhs, const EntityCommListInfo& rhs)
{ return !(lhs == rhs); }

struct IsInvalid
{
  bool operator()(const EntityCommListInfo& comm) const
  {
    return comm.key == EntityKey();
  }
};

inline
bool is_handled_generically(const RelationType relation_type)
{
  return relation_type == RelationType::USES || relation_type == RelationType::USED_BY;
}


/** \addtogroup stk_mesh_module
 *  \{
 */

//----------------------------------------------------------------------
/** \brief  Manager for an integrated collection of
 *          \ref stk::mesh::Entity "entities",
 *          \ref stk::mesh::Relation "entity relations", and
 *          \ref stk::mesh::Bucket "buckets" of
 *          \ref stk_mesh_field_data "field data".
 *
 *  Bulk data should be distributed among all processors.
 */
class BulkData {

  typedef tracking_allocator<unsigned char, FieldBase> field_data_allocator;

public:


  void optimize_buckets_at_modification_end(bool b) { m_optimize_buckets = b; }


#ifdef  STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  inline static BulkData & get( Entity entity);
#endif

  inline static BulkData & get( const Bucket & bucket);
  inline static BulkData & get( const Ghosting & ghost);
  inline static BulkData & get( const impl::BucketRepository & bucket_repo );

  enum BulkDataSyncState { MODIFIABLE = 1 , SYNCHRONIZED = 2 };

  ~BulkData();

  /** \brief  Construct mesh bulk data manager conformal to the given
   *          \ref stk::mesh::MetaData "meta data manager" and will
   *          distribute bulk data over the given parallel machine.
   *
   *  - The maximum number of entities per bucket may be supplied.
   *  - The bulk data is in the synchronized or "locked" state.
   */
  BulkData( MetaData & mesh_meta_data ,
            ParallelMachine parallel ,
            unsigned bucket_max_size = 1000 ,
            bool use_memory_pool = true
#ifdef SIERRA_MIGRATION
            , bool add_fmwk_data = false
#endif
            , bool maintain_fast_indices = true
            , ConnectivityMap* connectivity_map = NULL
            );

  FieldBase* coordinate_field() const { return m_coord_field; }
  void set_coordinate_field(FieldBase* coord_field) { m_coord_field = coord_field; }

  //------------------------------------
  /** \brief  The meta data manager for this bulk data manager. */
  const MetaData & mesh_meta_data() const { return m_mesh_meta_data ; }

  /** \brief  The parallel machine */
  ParallelMachine parallel() const { return m_parallel_machine ; }

  /** \brief  Size of the parallel machine */
  int parallel_size()   const { return m_parallel_size ; }

  /** \brief  Rank of the parallel machine's local processor */
  int parallel_rank()   const { return m_parallel_rank ; }

#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  int bulk_data_id() const { return m_bulk_data_id; }
#endif

  //------------------------------------
  /** \brief  Bulk data has two states:
   *          guaranteed to be parallel synchronized or
   *          modification in progress and may be parallel inconsistent.
   */
  BulkDataSyncState synchronized_state() const { return m_sync_state ; }

  /** \brief  Count of the number of times that the bulk data has been
   *          parallel synchronized.  This count gets updated with
   *          each call to 'modification_end'.
   */
  size_t synchronized_count() const { return m_sync_count ; }

  /** \brief  Begin a modification phase during which the mesh bulk data
   *          could become parallel inconsistent.  This is a parallel
   *          synchronous call.  The first time this method is called
   *          the mesh meta data is verified to be committed and
   *          parallel consistent.  An exception is thrown if this
   *          verification fails.
   *
   *  \return  True if transitioned out of the guaranteed
   *           parallel consistent state to the "ok to modify" state.
   *           False if already in this state.
   */
  bool modification_begin();

  /** \brief  Parallel synchronization of modifications and
   *          transition to the guaranteed parallel consistent state.
   *
   *  Parallel synchronization of accumulated local modifications
   *  is probably an expensive operation.  Operations include:
   *  - Determining ownership and sharing of created entities.
   *  - Synchronizing entity membership in parts for shared entities.
   *  - Refreshing the shared entities ghosting (e.g. aura).
   *  - Updating ghost entities that have change part membership.
   *  - Sorting buckets' entities for a well-defined ordering.
   *
   *  \return  True if transitioned from the "ok to modify" state.
   *           False if already already in this state.
   *
   *  \exception  If modification resolution errors occur then
   *              a parallel-consistent exception will be thrown.
   */
  bool modification_end();

  void verify_relations(const Bucket & bucket, unsigned bucket_ordinal, EntityRank rank) const;

  void verify_SOA_relations() const;

  bool final_modification_end()
  {
    const bool mod_flag =  modification_end();

    //call modification_begin and end one last time to free deleted entities
    modification_begin();
    modification_end();

    m_mesh_finalized = true;

    return mod_flag;
  }

  /** \brief  Give away ownership of entities to other parallel processes.
   *
   *  A parallel-synchronous operation while the mesh is in the
   *  ok-to-modify state.
   *
   *  Each owning process inputs a list of entities and the
   *  new owning process.  Upon completion of the call the owning
   *  processes have the newly owned entities as well as the closure
   *  of those entities (without changing the ownership of the closure
   *  entities).  If a previous owner no longer needs a
   *  changed-owner entity to support the closure of a still-owned
   *  entity then the changed-owner entity is deleted from that process.
   *  All ghosts of all entities effected by the changed ownerships
   *  deleted.
   */
  void change_entity_owner( const std::vector<EntityProc> & arg_change);

  /** \brief  Rotate the field data of multistate fields.
   *
   *  <PRE>
   *  Rotation of states:
   *    StateN   <- StateNP1 (StateOld <- StateNew)
   *    StateNM1 <- StateN   (StateNM1 <- StateOld)
   *    StateNM2 <- StateNM1
   *    StateNM3 <- StateNM2
   *    StateNM3 <- StateNM2
   *  </PRE>
   */
  void update_field_data_states();

  /** \brief  Copy field data from src entity to Dest entity
   *           - Fields that exist on the src that don't exist on the dest will
   *             be ignored
   *           - Fields that exist on the dest that don't exist on the src will
   *             be zeroed or initialized with the Field-specified initial-value.
   */
  void copy_entity_fields( Entity src, Entity dst)
  {
    //TODO fix const correctness for src
    MeshIndex & src_mesh_idx = mesh_index(src);
    MeshIndex & dst_mesh_idx = mesh_index(dst);

    //// Pre-upgrade stk_mesh did not have this restriction, and it was easy enough to remove.
    //    ThrowAssert(src_mesh_idx.bucket->entity_rank() == dst_mesh_idx.bucket->entity_rank());

    copy_entity_fields_callback(dst_mesh_idx.bucket->entity_rank(),
                         dst_mesh_idx.bucket->bucket_id(),
                         dst_mesh_idx.bucket_ordinal,
                         src_mesh_idx.bucket->entity_rank(),
                         src_mesh_idx.bucket->bucket_id(),
                         src_mesh_idx.bucket_ordinal);
  }

  //------------------------------------
  /** \brief  Query all buckets of a given entity rank
   *  Don't call inside BucketRepository member functions!
   */
  const std::vector<Bucket*> & buckets( EntityRank rank ) const
  { return m_bucket_repository.buckets(rank); }

#ifndef SWIG //SRK
  //  Don't call inside BucketRepository member functions!
  AllBucketsRange get_bucket_range() const
  {
    return m_bucket_repository.get_bucket_range();
  }

  //  Don't call inside BucketRepository member functions!
  AllBucketsRange get_bucket_range( EntityRank rank ) const
  {
    return m_bucket_repository.get_bucket_range(rank);
  }
#endif

  /** \brief  Get entity with a given rank and id */
  Entity get_entity( EntityRank ent_rank , EntityId entity_id ) const {
    require_good_rank_and_id(ent_rank, entity_id);
    return m_entity_repo.get_entity( EntityKey(ent_rank, entity_id));
  }

  /** \brief  Get entity with a given key */
  Entity get_entity( const EntityKey key ) const  {
    return m_entity_repo.get_entity(key);
  }

  //------------------------------------
  /** \brief  Create or retrieve a locally owned entity of a
   *          given rank and id.
   *
   *  A parallel-local operation.
   *
   *  The entity is created as locally owned and a member of the input
   *  mesh parts.  The entity a member of the meta data's locally owned
   *  mesh part and the entity's owner_rank() == parallel_rank().
   *
   *  If two or more processes create an entity of the same rank
   *  and identifier then the sharing and ownership of these entities
   *  will be resolved by the call to 'modification_end'.
   */
  Entity declare_entity( EntityRank ent_rank , EntityId ent_id , const PartVector& parts);

  Entity declare_entity( EntityRank ent_rank , EntityId ent_id , Part& part)
  {
    PartVector parts(1, &part);
    return declare_entity( ent_rank, ent_id, parts);
  }

  /** This overloading of declare_entity that doesn't take a part
   * creates the new entity in the 'universal' part.
   */
  Entity declare_entity( EntityRank ent_rank , EntityId ent_id);


  void change_entity_id( EntityId id, Entity entity);

  /** \brief  Change the parallel-locally-owned entity's
   *          part membership by adding and/or removing parts
   *
   *  A parallel-local operation.
   *
   *  If the locally owned entity is shared or ghosted then
   *  the change will be propagated to the sharing or ghosting
   *  processes by modification_end.
   */
  void change_entity_parts( Entity entity,
      const PartVector & add_parts ,
      const PartVector & remove_parts = PartVector() )
  {
    change_entity_parts(entity,
                        add_parts.begin(), add_parts.end(),
                        remove_parts.begin(), remove_parts.end());
  }

  //Optional parameter 'always_propagate_internal_changes' is always true except when this function
  //is being called from the sierra-framework. The fmwk redundantly does its own propagation of the
  //internal part changes (mostly induced-part stuff), so it's a performance optimization to avoid
  //the propagation that stk-mesh does.
  template<typename AddIterator, typename RemoveIterator>
  void change_entity_parts( Entity entity,
                            AddIterator begin_add_parts, AddIterator end_add_parts,
                            RemoveIterator begin_remove_parts, RemoveIterator end_remove_parts,
                            bool always_propagate_internal_changes=true );

  /** \brief  Request the destruction an entity on the local process.
   *
   * \paragraph destroy_requirements  Requirements
   *
   *  An entity cannot be the 'to' member of a relation.
   *  These relations must first be explicitly removed or the
   *  'from' entity be explicitly destroyed.
   *
   * \paragraph destroy_locally_owned  Destroy Locally Owned
   *
   *  Destrunction of entities in the 'locally_owned_part' schedules
   *  all ghost copies of that entity for destruction during
   *  modification_end.  If the entity ia shared with
   *  another process and that process does not also destroy the
   *  entity then ownership of the entity will be transfered to
   *  a sharing process during modification_end.
   *
   * \paragraph destroy_globally_shared  Destroy Locally Used
   *
   *  Entities in the 'globally_shared_part' are deleted
   *  on the local process and removed from the sharing lists on
   *  other processes during modication_end.
   *
   * \paragraph destroy_ghosted  Destroy Ghosted
   *
   *  Entities not in the 'locally_owned_part' and 'globally_shared_part'
   *  are ghosted.
   *  These entities are removed from all ghosting lists
   *  during 'modification_end'.
   *
   *  \return  True if the request for destruction is accepted; i.e.,
   *           if the entity is not the 'to' member of a relation.
   */
  bool destroy_entity( Entity entity );

  //------------------------------------

  /** \brief Generate a set of entites with globally unique id's
   *
   *  Each processor fills a request vector asking for a number of new
   *  entities of the given ranks.
   *
   *  ex. request = { 0, 4,  8}
   *  request 0 entites of rank 0, 4 entites of rank 1, and 8 entites
   *  of rank 2
   */
  void generate_new_entities(const std::vector<size_t>& requests,
      std::vector<Entity>& requested_entities);

  //------------------------------------
  /** \brief  Declare a relation and its converse between
   *          entities in the same mesh.
   *
   *  A parallel-local mesh modificaton operation.
   *
   *  This mapping ( e_from , local_id ) -> e_to  must be unique.
   *
   *  Relations between entities induces part membership as follows.
   *  1) If 'e_from' is a member of 'part' and
   *     part.primary_entity_rank() == entity_rank(e_from)
   *     then 'e_to' has induced membership in 'part'.
   *  2) If there exists a part relation 'part_rel' such that
   *     'e_from' is a member of part_rel.m_root and
   *     the entity relation conforms to the part relation
   *     then 'e_to' has induced membership in part_rel.m_target.
   *
   * Note that relation-declarations must be symmetric across all
   * sharers of the involved entities within a modification cycle.
   */
  void declare_relation( Entity e_from ,
      Entity e_to ,
      const RelationIdentifier local_id,
      Permutation permutation = static_cast<Permutation>(0));

  /** \brief  Declare a collection of relations by simply iterating
   *          the input and calling declare_relation on each entry.
   */
  void declare_relation( Entity entity, const std::vector<Relation> & rel);

  /** \brief  Remove all relations between two entities.
   *
   *  If the relation induced a part membership for 'e_to' and 'e_to'
   *  is not shared with another processor then that part membership
   *  is removed if and only if there does not exist another relation
   *  inducing the same part membership.
   *  If 'e_to' is shared then the check for removing the induced
   *  relatinship does not occur for that entity until the call to
   *  'modification_end'.
   *  The local_id arg is used to differentiate the case when there are
   *  multiple relationships between e_from and e_to.
   *
   *  Returns true if we were able to destroy the relation.
   */
  bool destroy_relation( Entity e_from ,
                         Entity e_to,
                         const RelationIdentifier local_id );


  // Check if entity has a specific relation to an entity of subcell_rank
  bool relation_exist( const Entity entity, EntityRank subcell_rank, RelationIdentifier subcell_id )
  {
    bool found = false;
    Entity const * rel_entity_it = bucket(entity).begin_entities(bucket_ordinal(entity),subcell_rank);
    const unsigned num_rel = bucket(entity).num_connectivity(bucket_ordinal(entity),subcell_rank);
    ConnectivityOrdinal const * rel_ord_it = bucket(entity).begin_ordinals(bucket_ordinal(entity),subcell_rank);

    for (unsigned i=0 ; i < num_rel ; ++i) {
      if (rel_ord_it[i] == static_cast<ConnectivityOrdinal>(subcell_id) && is_valid(rel_entity_it[i])) {
        found = true;
        break;
      }
    }

    return found;
  }

  /** \brief  Determine the polarity of the local side,
   *          more efficient if the local_side_id is known.
   */
  bool element_side_polarity( const Entity elem ,
      const Entity side , unsigned local_side_id ) const
  {
    // 09/14/10:  TODO:  tscoffe:  Will this work in 1D??
    const bool is_side = entity_rank(side) != MetaData::EDGE_RANK;
    const CellTopologyData * const elem_top = get_cell_topology( bucket(elem) ).getCellTopologyData();

    const unsigned side_count = ! elem_top ? 0 : (
        is_side ? elem_top->side_count
            : elem_top->edge_count );

    ThrowErrorMsgIf( elem_top == NULL,
        "For Element[" << identifier(elem) << "], element has no defined topology");

    ThrowErrorMsgIf( static_cast<unsigned>(side_count) <= local_side_id,
        "For Element[" << identifier(elem) << "], " <<
        "side: " << identifier(side) << ", " <<
        "local_side_id = " << local_side_id <<
        " ; unsupported local_side_id");

    const CellTopologyData * const side_top =
        is_side ? elem_top->side[ local_side_id ].topology
            : elem_top->edge[ local_side_id ].topology ;

    const unsigned * const side_map =
        is_side ? elem_top->side[ local_side_id ].node
            : elem_top->edge[ local_side_id ].node ;

    Entity const *elem_nodes = begin_node_entities(elem);
    Entity const *side_nodes = begin_node_entities(side);
    const unsigned n = side_top->node_count;
    bool good = false ;
    for ( unsigned i = 0 ; !good && i < n ; ++i ) {
        good = true;
        for ( unsigned j = 0; good && j < n ; ++j ) {
          good = side_nodes[(j+i)%n] == elem_nodes[ side_map[j] ];
        }
    }
    return good ;
  }

  //------------------------------------
  //------------------------------------
  /** \brief  All entities with communication information. */
  const std::vector<EntityCommListInfo> & comm_list() const
    { return m_entity_comm_list; }

  //------------------------------------
  /** \brief  Query the shared-entity aura.
   *          Is likely to be stale if ownership or sharing has changed
   *          and the 'modification_end' has not been called.
   */
  Ghosting & shared_aura() const { return * m_ghosting[1] ; }

  /** \brief Asymmetric parallel relations for owner-to-ghosted mesh entities.
   *
   *  - A collective parallel operation that must have the
   *    same name on all processors of this distributed mesh.
   */
  Ghosting & create_ghosting( const std::string & name );

  /** \brief  Change the members of a ghosting list on the sending processor.
   *
   *  - A collective parallel operation.
   *  - The ghosting must belong to this mesh.
   *  - Cannot change the 'shared_aura' in this manner.
   *  - Add locally owned entities to the input ghosting on the given
   *    destination processor.  The closure of the input entities
   *    will be ghosted.
   *  - Request removal of ghosted entities on the ghosting processor.
   *    This request will only be honored if the ghosted entity is
   *    not in the closure of another ghosted entity which will remain
   *    in or be added to this ghosting.
   *
   *  To destroy all ghosted entities in a particular ghosting:
   *    mesh.change_ghosting( ghosts , std::vector<EntityProc>() ,
   *                                   ghosts.receive() );
   */
  void change_ghosting( Ghosting & ghosts,
                        const std::vector<EntityProc> & add_send ,
                        const std::vector<EntityKey> & remove_receive );

  /** \brief  Empty every single Ghosting.
   *          Same result, but more efficient than, calling
   *          change_ghosting to remove every single ghosted entity.
   */
  void destroy_all_ghosting();

  /** \brief  Query all ghostings */
  const std::vector<Ghosting*> & ghostings() const { return m_ghosting ; }

  /** \brief  Entity Comm functions that are now moved to BulkData
   */
  PairIterEntityComm entity_comm(const EntityKey & key) const { return m_entity_comm_map.comm(key); }
  PairIterEntityComm entity_comm_sharing(const EntityKey & key) const { return m_entity_comm_map.sharing(key); }
  PairIterEntityComm entity_comm(const EntityKey & key, const Ghosting & sub ) const { return m_entity_comm_map.comm(key,sub); }
  bool entity_comm_insert(Entity entity, const EntityCommInfo & val) { return m_entity_comm_map.insert(entity_key(entity), val, parallel_owner_rank(entity)); }
  bool entity_comm_erase(  const EntityKey & key, const EntityCommInfo & val) { return m_entity_comm_map.erase(key,val); }
  bool entity_comm_erase(  const EntityKey & key, const Ghosting & ghost) { return m_entity_comm_map.erase(key,ghost); }
  void entity_comm_clear_ghosting(const EntityKey & key ) { m_entity_comm_map.comm_clear_ghosting(key); }
  void entity_comm_clear(const EntityKey & key) { m_entity_comm_map.comm_clear(key); }
  int entity_comm_owner(const EntityKey & key) const;

  // Comm-related convenience methods

  bool in_shared(EntityKey key) const { return !entity_comm_sharing(key).empty(); }

  bool in_shared(EntityKey key, int proc) const;

  bool in_receive_ghost( EntityKey key ) const;

  bool in_receive_ghost( const Ghosting & ghost , EntityKey entity ) const;

  bool in_send_ghost( EntityKey key) const;

  bool in_send_ghost( EntityKey key , int proc ) const;

  bool in_ghost( const Ghosting & ghost , EntityKey key , int proc ) const;

  void comm_procs( EntityKey key, std::vector<int> & procs ) const;

  void comm_procs( const Ghosting & ghost ,
                   EntityKey key, std::vector<int> & procs ) const;


  const MeshIndex& mesh_index(Entity entity) const
  {
    return m_mesh_indexes[entity.local_offset()];
  }
  MeshIndex& mesh_index(Entity entity)
  {
    return m_mesh_indexes[entity.local_offset()];
  }

  void set_mesh_index(Entity entity, Bucket * in_bucket, unsigned ordinal )
  {
    // The trace statement forces this method to be defined after Entity
    TraceIfWatching("stk::mesh::BulkData::set_mesh_index", LOG_ENTITY, entity_key(entity));

    if (in_bucket != NULL) {
      ThrowAssertMsg(in_bucket->size() >= ordinal, "Detected bad bucket/ordinal.");
    }
    MeshIndex &mesh_idx = mesh_index(entity);
    mesh_idx.bucket = in_bucket;
    mesh_idx.bucket_ordinal = ordinal;
  }

  ///
  /// Entity Member Functions:
  ///
  EntityId identifier(Entity entity) const
  {
    ThrowAssert(this == &BulkData::get(entity));
    ThrowAssertMsg(entity.local_offset() < m_entity_keys.size(), "For bulk " << m_bulk_data_id << ", Entity offset " << entity.local_offset() << " is beyond max value " << m_entity_keys.size() - 1);
    return m_entity_keys[entity.local_offset()].id();
  }

  EntityRank entity_rank(Entity entity) const
  {
    //Implementation will be replaced with this:
    ThrowAssert(this == &BulkData::get(entity));
    ThrowAssertMsg(entity.local_offset() < m_entity_keys.size(), "For bulk " << m_bulk_data_id << ", Entity offset " << entity.local_offset() << " is beyond max value " << m_entity_keys.size() - 1);
    return m_entity_keys[entity.local_offset()].rank();
  }

  EntityKey entity_key(Entity entity) const
  {
    ThrowAssert(this == &BulkData::get(entity));
    ThrowAssertMsg(entity.local_offset() < m_entity_keys.size(), "For bulk " << m_bulk_data_id << ", Entity offset " << entity.local_offset() << " is beyond max value " << m_entity_keys.size() - 1);
    return m_entity_keys[entity.local_offset()];
  }

  void set_entity_key(Entity entity, EntityKey key)
  {
    ThrowAssert(this == &BulkData::get(entity));
    ThrowAssertMsg(entity.local_offset() < m_entity_keys.size(), "For bulk " << m_bulk_data_id << ", Entity offset " << entity.local_offset() << " is beyond max value " << m_entity_keys.size() - 1);
    m_entity_keys[entity.local_offset()] = key;

    //// Prefer the above because an additional place that can
    //// cause the array size to grow makes the code confusing.
    //// If there are multiple places that can push_back(), they
    //// all would need to check whether it is actually necessary.
    //    if (entity.local_offset() == m_entity_keys.size()) {
    //      m_entity_keys.push_back(key);
    //    }
    //    else {
    //      m_entity_keys[entity.local_offset()] = key;
    //    }
  }

  bool in_index_range(Entity entity) const
  {
    return entity.local_offset() < m_entity_states.size();
  }

  bool is_valid(Entity entity) const
  {
    ThrowAssert(entity.local_offset() < m_entity_states.size());
    return m_entity_states[entity.local_offset()] != Deleted;
  }

  Bucket & bucket(Entity entity) const
  {
    ThrowAssert(entity.local_offset() < m_mesh_indexes.size());
    return *m_mesh_indexes[entity.local_offset()].bucket;
  }

  Bucket * bucket_ptr(Entity entity) const
  {
    if (entity.local_offset() >= m_mesh_indexes.size()) return NULL;
    return m_mesh_indexes[entity.local_offset()].bucket;
  }

  unsigned bucket_ordinal(Entity entity) const
  {
    if (entity.local_offset() >= m_mesh_indexes.size()) return 0;
    return m_mesh_indexes[entity.local_offset()].bucket_ordinal;
  }

  inline EntityState state(Entity entity) const
  {
    return m_entity_states[entity.local_offset()];
  }
  inline void set_state(Entity entity, EntityState entity_state)
  {
    m_entity_states[entity.local_offset()] = entity_state;
  }

  inline int parallel_owner_rank(Entity entity) const
  {
    return bucket(entity).parallel_owner_rank(bucket_ordinal(entity));
  }

  inline bool set_parallel_owner_rank(Entity entity, int in_owner_rank)
  {
     TraceIfWatching("stk::mesh::BulkData::set_entity_owner_rank", LOG_ENTITY, entity_key(entity));
     DiagIfWatching(LOG_ENTITY, entity_key(entity), "new owner: " << in_owner_rank);

     //     bool changed = entity.m_entityImpl->set_owner_rank(in_owner_rank);
     //     if ( changed ) {
     //       this->modified(entity);
     //     }
     //     return changed;
     int &rank = bucket(entity).m_owner_ranks[bucket_ordinal(entity)];
     if ( in_owner_rank != rank ) {
       rank = in_owner_rank;
       modified(entity);
       return true;
     }
     return false;
  }

  inline size_t synchronized_count(Entity entity) const
  {
    return m_entity_sync_counts[entity.local_offset()];
  }
  inline void set_synchronized_count(Entity entity, size_t sync_count)
  {
    m_entity_sync_counts[entity.local_offset()] = sync_count;
  }

  /**
   * Mark this entity as modified (only changes from Unchanged
   * to Modified). Propagates the modification to higher-ranking
   * entities related to this entity. In other words, based on our
   * modification model, all entities that have modified_entity in their
   * closure must also be marked as modified.
   */
  void modified(Entity entity);

  size_t count_relations(Entity entity) const;

  bool has_no_relations(Entity entity) const;

  ////
  //// Rank Id based functions --- only use outside of modification cycle.
  //// Didn't use templates because of specialization issues wrt the Entity
  //// argument versions.
  ////

  Entity get_entity(Node node) const;
  Entity get_entity(Edge edge) const;
  Entity get_entity(Face facee) const;
  Entity get_entity(Element elem) const;

  Node get_node(Entity entity) const;
  Edge get_edge(Entity entity) const;
  Face get_face(Entity entity) const;
  Element get_element(Entity entity) const;

  template <class RankType>
  RankType
  inline get_rank_id(Entity entity)
  {
    MeshIndex mi = m_mesh_indexes[entity.local_offset()];
    ThrowAssert(mi.bucket->entity_rank() == impl::IndexToRank<RankType>::rank);
    return make_rank_id<RankType>(mi.bucket->bucket_id(), mi.bucket_ordinal);
  }

  ////
  //// NEW SAMBA-LIKE RELATIONS GETTERS
  ////

  // TODO: support beyond-element rank (e.g. constaint) connectivity

  Node const* begin_nodes(Entity entity) const;
  Node const* end_nodes  (Entity entity) const;

  Edge const* begin_edges(Entity entity) const;
  Edge const* end_edges  (Entity entity) const;

  Face const* begin_faces(Entity entity) const;
  Face const* end_faces  (Entity entity) const;

  Element const* begin_elements(Entity entity) const;
  Element const* end_elements  (Entity entity) const;

  Entity const* begin_entities(Entity entity, EntityRank rank) const;
  Entity const* begin_node_entities(Entity entity) const;
  Entity const* begin_edge_entities(Entity entity) const;
  Entity const* begin_face_entities(Entity entity) const;
  Entity const* begin_element_entities(Entity entity) const;
  Entity const* begin_other_entities(Entity entity) const
  { return begin_entities(entity, stk::topology::CONSTRAINT_RANK); }

  ConnectivityOrdinal const* begin_ordinals(Entity entity, EntityRank rank) const;
  ConnectivityOrdinal const* begin_node_ordinals(Entity entity) const;
  ConnectivityOrdinal const* begin_edge_ordinals(Entity entity) const;
  ConnectivityOrdinal const* begin_face_ordinals(Entity entity) const;
  ConnectivityOrdinal const* begin_element_ordinals(Entity entity) const;
  ConnectivityOrdinal const* begin_other_ordinals(Entity entity) const
  { return begin_ordinals(entity, stk::topology::CONSTRAINT_RANK); }

  Permutation const* begin_permutations(Entity entity, EntityRank rank) const;
  Permutation const* begin_node_permutations(Entity entity) const;
  Permutation const* begin_edge_permutations(Entity entity) const;
  Permutation const* begin_face_permutations(Entity entity) const;
  Permutation const* begin_element_permutations(Entity entity) const;
  Permutation const* begin_other_permutations(Entity entity) const
  { return begin_permutations(entity, stk::topology::CONSTRAINT_RANK); }

  unsigned num_connectivity(Entity entity, EntityRank rank) const;
  unsigned num_nodes(Entity entity) const;
  unsigned num_edges(Entity entity) const;
  unsigned num_faces(Entity entity) const;
  unsigned num_elements(Entity entity) const;
  unsigned num_other(Entity entity) const
  { return end_other_entities(entity) - begin_other_entities(entity); }

  unsigned count_valid_connectivity(Entity entity, EntityRank rank) const;
  unsigned count_valid_connectivity(Entity entity) const;

  Entity const* end_entities(Entity entity, EntityRank rank) const;
  Entity const* end_node_entities(Entity entity) const;
  Entity const* end_edge_entities(Entity entity) const;
  Entity const* end_face_entities(Entity entity) const;
  Entity const* end_element_entities(Entity entity) const;
  Entity const* end_other_entities(Entity entity) const;

  ConnectivityOrdinal const* end_ordinals(Entity entity, EntityRank rank) const;
  ConnectivityOrdinal const* end_node_ordinals(Entity entity) const;
  ConnectivityOrdinal const* end_edge_ordinals(Entity entity) const;
  ConnectivityOrdinal const* end_face_ordinals(Entity entity) const;
  ConnectivityOrdinal const* end_element_ordinals(Entity entity) const;

  Permutation const* end_permutations(Entity entity, EntityRank rank) const;
  Permutation const* end_node_permutations(Entity entity) const;
  Permutation const* end_edge_permutations(Entity entity) const;
  Permutation const* end_face_permutations(Entity entity) const;
  Permutation const* end_element_permutations(Entity entity) const;

  // Return index (offset) of query ordinal if found, num_connectivity otherwise.
  unsigned find_ordinal(Entity entity, EntityRank rank, ConnectivityOrdinal ordinal) const;

  bool has_permutation(Entity entity, EntityRank rank) const;

  unsigned field_data_size(const FieldBase& f, Entity e) const
  {
    return field_data_size_per_entity(f, bucket(e));
  }

  unsigned field_data_size_per_entity(const FieldBase& f, const Bucket& b) const
  {
    const EntityRank rank   = b.entity_rank();
    return m_field_meta_data[m_num_fields * rank + f.mesh_meta_data_ordinal()][b.bucket_id()].m_size;
  }

  size_t total_field_data_footprint(const FieldBase &f, EntityRank rank) const
  {
    return m_bucket_repository.total_field_data_footprint(f, rank);
  }

  size_t total_field_data_footprint(EntityRank rank) const;

  template<class FieldType>
  typename FieldTraits<FieldType>::data_type*
  field_data(const FieldType & f, const Bucket& b, unsigned bucket_ord) const
  {
    const EntityRank rank         = b.entity_rank();
    // TODO - should be const ref?
    FieldMetaData field_meta_data = m_field_meta_data[m_num_fields * rank + f.mesh_meta_data_ordinal()][b.bucket_id()];
    return reinterpret_cast<typename FieldTraits<FieldType>::data_type*>(field_meta_data.m_data + field_meta_data.m_size * bucket_ord);
  }

  template<class FieldType>
  typename FieldTraits<FieldType>::data_type*
  field_data(const FieldType & f, Entity e) const
  {
    const MeshIndex& mi           = mesh_index(e);
    return field_data(f, *mi.bucket, mi.bucket_ordinal);
  }

  template<class FieldType, typename Index>
  typename FieldTraits<FieldType>::data_type*
  field_data(const FieldType & f, Index e_idx) const
  {
    const EntityRank rank = impl::IndexToRank<Index>::rank;
    const unsigned b_id  = EXTRACT_BUCKET_ID(e_idx);
    const unsigned b_ord = EXTRACT_BUCKET_ORDINAL(e_idx);
    FieldMetaData field_meta_data = m_field_meta_data[m_num_fields * rank + f.mesh_meta_data_ordinal()][b_id];
    return reinterpret_cast<typename FieldTraits<FieldType>::data_type*>(field_meta_data.m_data + field_meta_data.m_size * b_ord);
  }

  const FieldBase::Restriction::size_type * field_data_stride( const FieldBase & field, const Bucket& b ) const
  {
    const EntityRank rank = b.entity_rank();
    return m_field_meta_data[m_num_fields * rank + field.mesh_meta_data_ordinal()][b.bucket_id()].m_stride;
  }

  //reserves space for a new entity, or reclaims space from a previously-deleted entity
  size_t generate_next_local_offset();

#ifdef SIERRA_MIGRATION
  //strictly a transition aid!!! don't add new usage of this!
  void set_fmwk_bulk_data(const void* fmwk_bulk_ptr)
  {
    m_fmwk_bulk_ptr = fmwk_bulk_ptr;
  }

  //strictly a transition aid!!! don't add new usage of this!
  const void* get_fmwk_bulk_data() const
  {
    return m_fmwk_bulk_ptr;
  }

  const RelationVector& aux_relations(Entity entity) const
  {
    ThrowAssert(m_add_fmwk_data);
    if (m_fmwk_aux_relations[entity.local_offset()] == NULL) {
      m_fmwk_aux_relations[entity.local_offset()] = new RelationVector();
    }
    return *m_fmwk_aux_relations[entity.local_offset()];
  }

  RelationVector& aux_relations(Entity entity)
  {
    ThrowAssert(m_add_fmwk_data);
    if (m_fmwk_aux_relations[entity.local_offset()] == NULL) {
      m_fmwk_aux_relations[entity.local_offset()] = new RelationVector();
    }
    return *m_fmwk_aux_relations[entity.local_offset()];
  }

  RelationIterator internal_begin_relation(Entity entity, const Relation::RelationType relation_type) const
  {
    ThrowAssert(m_add_fmwk_data);
    if (is_handled_generically(relation_type)) {
      ThrowErrorMsg("stk::Mesh::BulkData::internal_begin_relation(..) requests native stk::mesh relation type");
      return RelationIterator();
    }
    else {
      return aux_relations(entity).begin();
    }
  }

  RelationIterator internal_end_relation(Entity entity, const Relation::RelationType relation_type) const
  {
    ThrowAssert(m_add_fmwk_data);
    if (is_handled_generically(relation_type)) {
      ThrowErrorMsg("stk::Mesh::BulkData::internal_begin_relation(..) requests native stk::mesh relation type");
      return RelationIterator();
    }
    else {
      return aux_relations(entity).end();
    }
  }

  void compress_relation_capacity(Entity entity)
  {
    RelationVector &rels = aux_relations(entity);
    RelationVector tmp(rels);
    tmp.swap(rels);
  }

  int global_id(size_t local_offset) const
  {
    ThrowAssert(m_add_fmwk_data);
    ThrowAssertMsg(local_offset < m_fmwk_global_ids.size(), "Out of bounds local_offset " << local_offset << ", max is " << m_fmwk_global_ids.size());
    return m_fmwk_global_ids[local_offset];
  }
  int global_id(Entity entity) const
  { return m_fmwk_global_ids[entity.local_offset()]; }

  void set_global_id(Entity entity, int id)
  {
    ThrowAssert(m_add_fmwk_data);
    ThrowAssertMsg(entity.local_offset() < m_fmwk_global_ids.size(), "Out of bounds local_offset " << entity.local_offset() << ", max is " << m_fmwk_global_ids.size());
    m_fmwk_global_ids[entity.local_offset()] = id;
  }

  unsigned local_id(size_t local_offset) const
  {
    ThrowAssert(m_add_fmwk_data);
    ThrowAssertMsg(local_offset < m_fmwk_local_ids.size(), "Out of bounds local_offset " << local_offset << ", max is " << m_fmwk_local_ids.size());
    return m_fmwk_local_ids[local_offset];
  }
  unsigned local_id(Entity entity) const
  { return local_id(entity.local_offset()); }

  void set_local_id(Entity entity, unsigned id)
  {
    ThrowAssert(m_add_fmwk_data);
    ThrowAssertMsg(entity.local_offset() < m_fmwk_local_ids.size(), "Out of bounds local_offset " << entity.local_offset() << ", max is " << m_fmwk_local_ids.size());
    m_fmwk_local_ids[entity.local_offset()] = id;
  }

  const void* get_shared_attr(Entity entity) const
  {
    ThrowAssert(m_add_fmwk_data);
    ThrowAssertMsg(entity.local_offset() < m_fmwk_shared_attrs.size(), "Out of bounds local_offset " << entity.local_offset() << ", max is " << m_fmwk_shared_attrs.size());
    return m_fmwk_shared_attrs[entity.local_offset()];
  }

  template <typename SharedAttr>
  void set_shared_attr(Entity entity, SharedAttr* attr)
  {
    ThrowAssert(m_add_fmwk_data);
    ThrowAssertMsg(entity.local_offset() < m_fmwk_shared_attrs.size(), "Out of bounds local_offset " << entity.local_offset() << ", max is " << m_fmwk_shared_attrs.size());
    m_fmwk_shared_attrs[entity.local_offset()] = attr;
  }

  int get_connect_count(Entity entity) const
  {
    ThrowAssert(m_add_fmwk_data);
    ThrowAssertMsg(entity.local_offset() < m_fmwk_connect_counts.size(), "Out of bounds local_offset " << entity.local_offset() << ", max is " << m_fmwk_connect_counts.size());
    return m_fmwk_connect_counts[entity.local_offset()];
  }

  void set_connect_count(Entity entity, int count)
  {
    ThrowAssert(m_add_fmwk_data);
    m_fmwk_connect_counts[entity.local_offset()] = count;
  }

  bool add_fmwk_data() const { return m_add_fmwk_data; }
#endif

  bool maintain_fast_indices() const { return m_maintain_fast_indices; }

  // Do not call!
  void internal_change_entity_key(EntityKey old_key, EntityKey new_key, Entity entity);

  // Do not call!  Just for a legacy test!
  impl::EntityRepository &get_entity_repository() { return m_entity_repo; }

private:

  void update_deleted_entities_container();

#ifdef  STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  static const int MAX_NUM_BULKDATA = 256;
  static BulkData * the_bulk_data_registry[MAX_NUM_BULKDATA];
#endif

#ifndef DOXYGEN_COMPILE

  // Forbidden
  BulkData();
  BulkData( const BulkData & );
  BulkData & operator = ( const BulkData & );

  //
  // Members
  //

  /** \brief  Parallel index for entity keys */
  parallel::DistributedIndex          m_entities_index;
  impl::EntityRepository              m_entity_repo;

  // Simply a list of data for entities that are being communicated
  std::vector<EntityCommListInfo>     m_entity_comm_list;

  std::vector<Ghosting*>              m_ghosting; /**< Aura is [1] */

  std::list<size_t>     m_deleted_entities;
  std::list<size_t>     m_deleted_entities_current_modification_cycle;
  FieldBase* m_coord_field;

  // Other information:
  MetaData &         m_mesh_meta_data;
  ParallelMachine    m_parallel_machine;
  int                m_parallel_size;
  int                m_parallel_rank;
  size_t             m_sync_count;
  BulkDataSyncState  m_sync_state;
  bool               m_meta_data_verified;
  bool               m_optimize_buckets;
  bool               m_mesh_finalized;
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  int                m_bulk_data_id;
#endif
#ifdef SIERRA_MIGRATION
  bool               m_add_fmwk_data; // flag that will add extra data to buckets to support fmwk
  const void*    m_fmwk_bulk_ptr;

public:
  mutable bool               m_check_invalid_rels; // TODO REMOVE
private:
#endif
  bool               m_maintain_fast_indices; // flag that will tell buckets to allocate+compute fast indices at modification end

  int m_num_fields;

  // The full database of comm info for all communicated entities.
  EntityCommDatabase m_entity_comm_map;

  // Arrays of Entity member-data indexed by entity.local_offset():
  std::vector<MeshIndex> m_mesh_indexes;

  std::vector<EntityKey>   m_entity_keys;
  std::vector<EntityState> m_entity_states;
  std::vector<size_t>      m_entity_sync_counts;

#ifdef SIERRA_MIGRATION
  // Extra data that fmwk needs to have on an entity. These vectors are indexed by local offset.

  mutable std::vector<RelationVector* > m_fmwk_aux_relations;   // Relations that can't be managed by STK such as PARENT/CHILD
  std::vector<int>                      m_fmwk_global_ids;
  std::vector<unsigned>                 m_fmwk_local_ids;
  std::vector<const void*>              m_fmwk_shared_attrs;
  std::vector<unsigned short>           m_fmwk_connect_counts;
#endif

  // There will be one of these per bucket
  struct FieldMetaData
  {
    int m_size; // num bytes per entity, 0 means bucket does not have this field
    const FieldBase::Restriction::size_type * m_stride;
    unsigned char* m_data;
  };
  // Outer index is (m_num_fields * entity rank) + field_ordinal, inner index
  // is bucket id, pair defines num bytes of data per entity and the
  // data for that field on that bucket

#ifdef STK_PROFILE_MEMORY
  typedef std::vector<FieldMetaData, tracking_allocator<FieldMetaData, FieldBase> > FieldMetaDataVector;
  typedef std::vector<FieldMetaDataVector, tracking_allocator<FieldMetaDataVector, FieldBase> > FieldMetaDataVectorVector;
#else
  typedef std::vector<FieldMetaData> FieldMetaDataVector;
  typedef std::vector<FieldMetaDataVector> FieldMetaDataVectorVector;
#endif
  FieldMetaDataVectorVector m_field_meta_data;

  // Outer index is rank, inner is bucket-id. This contains *all* field
  // data for a bucket.
  std::vector<std::vector<unsigned char*> > m_field_raw_data;

  impl::BucketRepository              m_bucket_repository; // needs to be destructed first!

  //
  // Internal methods
  //

  // Field callbacks

  void new_bucket_callback(EntityRank rank, unsigned const* part_ord_begin, unsigned const* part_ord_end, size_t capacity);

  void copy_entity_fields_callback(EntityRank dst_rank, unsigned dst_bucket_id, unsigned dst_bucket_ord,
                            EntityRank src_rank, unsigned src_bucket_id, unsigned src_bucket_ord);

  void destroy_bucket_callback(EntityRank rank, unsigned bucket_id, unsigned capacity);

  // id_map, indexed by new id, maps to old id
  void reorder_buckets_callback(EntityRank rank, const std::vector<unsigned>& id_map);

  void remove_entity_callback(EntityRank rank, unsigned bucket_id, unsigned bucket_ord);

  // Field update helper
  void internal_update_fast_field_data(bool skip_onestate_fields = false);

  // Misc

  void initialize_arrays();

  bool internal_declare_relation(Entity e_from, Entity e_to,
                                 RelationIdentifier local_id,
                                 unsigned sync_count, bool is_back_relation,
                                 Permutation permut);

  /** \brief  The meta data manager for this bulk data manager. */
  MetaData & meta_data() const { return m_mesh_meta_data ; }

  /** methods for managing arrays of entity member-data */

  void log_created_parallel_copy(Entity entity)
  {
    if (state(entity) == Created) {
      set_state(entity, Modified);
    }
  }

  /**
   * For all processors sharing an entity, find one to be the new
   * owner.
   */
  int determine_new_owner( Entity ) const ;

  /*  Entity modification consequences:
   *  1) Change entity relation => update via part relation => change parts
   *  2) Change parts => update forward relations via part relation
   *                  => update via field relation
   */
  template <typename PartT>
  void internal_change_entity_parts( Entity ,
                                     const std::vector<PartT> & add_parts ,
                                     const std::vector<PartT> & remove_parts,
                                     bool always_propagate_internal_changes=true);

  template <typename PartT>
  void internal_propagate_part_changes( Entity entity, const std::vector<PartT> & removed );

  void internal_change_ghosting( Ghosting & ghosts,
                                 const std::vector<EntityProc> & add_send ,
                                 const std::vector<EntityKey> & remove_receive,
                                 bool is_full_regen = false);

  bool internal_modification_end( bool regenerate_aura );
  void internal_resolve_shared_modify_delete();
  void internal_resolve_shared_modify_delete_second_pass();
  void internal_resolve_ghosted_modify_delete();
  void internal_resolve_parallel_create();
  void internal_resolve_shared_membership();

  void internal_update_distributed_index( std::vector<Entity> & shared_new );

  /** \brief  Regenerate the shared-entity aura,
   *          adding and removing ghosted entities as necessary.
   *
   *  - a collective parallel operation.
   */
  void internal_regenerate_shared_aura();

  void internal_basic_part_check(const Part* part,
                                 const unsigned ent_rank,
                                 const unsigned undef_rank,
                                 bool& intersection_ok,
                                 bool& rel_target_ok,
                                 bool& rank_ok) const;

  // Returns false if there is a problem. It is expected that
  // verify_change_parts will be called if quick_verify_change_part detects
  // a problem, therefore we leave the generation of an exception to
  // verify_change_parts. We want this function to be as fast as
  // possible.
  bool internal_quick_verify_change_part(const Part* part,
                                         const unsigned ent_rank,
                                         const unsigned undef_rank) const;

  template <typename PartT>
  void internal_verify_change_parts( const MetaData   & meta ,
                                     const Entity entity ,
                                     const std::vector<PartT> & parts ) const;

  void internal_change_owner_in_comm_data(const EntityKey& key, int new_owner);

  void internal_sync_comm_list_owners();

  //------------------------------------

  /** \name  Invariants/preconditions for MetaData.
   * \{
   */

  /** \brief  All non-const methods assert this */
  void require_ok_to_modify() const ;

  void require_entity_owner( const Entity entity, int owner) const ;

  void require_metadata_committed();

  void require_good_rank_and_id(EntityRank ent_rank, EntityId ent_id) const;

  void require_valid_relation( const char action[] ,
                               const BulkData & mesh ,
                               const Entity e_from ,
                               const Entity e_to );

  /** \} */

  //------------------------------------

  // FIXME: Remove this friend once unit-testing has been refactored
  friend class UnitTestModificationEndWrapper;
  friend class ::stk::mesh::MetaData;
  friend class ::stk::mesh::impl::Partition;
  friend class ::stk::mesh::impl::BucketRepository;
  friend class stk::mesh::Bucket; // for field callbacks

#endif /* DOXYGEN_COMPILE */
};

 /** \brief  Comparison operator for entities compares the entities' keys */
struct EntityLess {
  EntityLess(const BulkData& mesh) : m_mesh(&mesh) {}

  /** \brief  Comparison operator */
  bool operator()(const Entity lhs, const Entity rhs) const
  {
    const EntityKey lhs_key = m_mesh->in_index_range(lhs) ? m_mesh->entity_key(lhs) : EntityKey();
    const EntityKey rhs_key = m_mesh->in_index_range(rhs) ? m_mesh->entity_key(rhs) : EntityKey();
    return lhs_key < rhs_key;
  }

  bool operator()(const Entity lhs, const EntityKey & rhs) const
  {
    const EntityKey lhs_key = m_mesh->in_index_range(lhs) ? m_mesh->entity_key(lhs) : EntityKey();
    return lhs_key < rhs;
  }

  bool operator()( const EntityProc & lhs, const EntityProc & rhs) const
  {
    const EntityKey lhs_key = m_mesh->in_index_range(lhs.first) ? m_mesh->entity_key(lhs.first) : EntityKey() ;
    const EntityKey rhs_key = m_mesh->in_index_range(rhs.first) ? m_mesh->entity_key(rhs.first) : EntityKey() ;
    return lhs_key != rhs_key ? lhs_key < rhs_key : lhs.second < rhs.second ;
  }

  bool operator()( const EntityProc & lhs, const Entity rhs) const
  {
    const EntityKey lhs_key = m_mesh->in_index_range(lhs.first) ? m_mesh->entity_key(lhs.first) : EntityKey();
    const EntityKey rhs_key = m_mesh->in_index_range(rhs)       ? m_mesh->entity_key(rhs)       : EntityKey();
    return lhs_key < rhs_key;
  }

  bool operator()( const EntityProc & lhs, const EntityKey & rhs) const
  {
    const EntityKey lhs_key = m_mesh->in_index_range(lhs.first) ? m_mesh->entity_key(lhs.first) : EntityKey();
    return lhs_key < rhs ;
  }

  EntityLess& operator=(const EntityLess& rhs)
  {
    m_mesh = rhs.m_mesh;
    return *this;
  }

  const BulkData* m_mesh;
}; //struct EntityLess

#ifdef  STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
inline
BulkData & BulkData::get( Entity entity)
{
  return *the_bulk_data_registry[entity.bulk_data_id()];
}
#endif

inline
BulkData & BulkData::get( const Bucket & bucket) {
  return bucket.bulk_data();
}

inline
BulkData & BulkData::get( const Ghosting & ghost) {
  return ghost.bulk_data();
}

inline
BulkData & BulkData::get( const impl::BucketRepository & bucket_repo ) {
  return bucket_repo.mesh();
}

#define RANK_VAL_node stk::topology::NODE_RANK
#define RANK_VAL_edge stk::topology::NODE_RANK
#define RANK_VAL_face stk::topology::NODE_RANK
#define RANK_VAL_element stk::topology::NODE_RANK

//
// Define a begin/end pair
//
#define BEGIN_END_PAIR(rank_name, return_type, data_type)       \
                                                                \
inline                                                          \
return_type const* BulkData::begin_##rank_name##_##data_type(Entity entity) const \
{                                                                       \
  ThrowAssert(is_valid(entity));                                        \
  ThrowAssert(bucket_ptr(entity));                                      \
  ThrowAssert(this == &stk::mesh::BulkData::get(entity));               \
  const MeshIndex &mesh_idx = mesh_index(entity);                       \
  const Bucket &b = *mesh_idx.bucket;                                   \
  unsigned bucket_ord = mesh_idx.bucket_ordinal;                        \
  if (m_check_invalid_rels) {                                           \
    ThrowAssert(count_valid_connectivity(entity, RANK_VAL_##rank_name) == b.num_connectivity(bucket_ord, RANK_VAL_##rank_name)); \
  }                                                                     \
  return b.begin_##rank_name##_##data_type(bucket_ord);                 \
}                                                                       \
                                                                        \
inline                                                                  \
return_type const* BulkData::end_##rank_name##_##data_type(Entity entity) const \
{                                                                       \
  ThrowAssert(is_valid(entity));                                        \
  ThrowAssert(bucket_ptr(entity));                                      \
  ThrowAssert(this == &stk::mesh::BulkData::get(entity));               \
  const MeshIndex &mesh_idx = mesh_index(entity);                       \
  const Bucket &b = *mesh_idx.bucket;                                   \
  unsigned bucket_ord = mesh_idx.bucket_ordinal;                        \
  return b.end_##rank_name##_##data_type(bucket_ord);                   \
}

//
// Define all methods for a rank
//
#define RANK_FUNCTION_DEFS(rank_name, type_name)                       \
                                                                       \
BEGIN_END_PAIR(rank_name, Entity, entities)                            \
                                                                       \
BEGIN_END_PAIR(rank_name, ConnectivityOrdinal, ordinals)               \
                                                                       \
BEGIN_END_PAIR(rank_name, Permutation, permutations)                   \
                                                                       \
inline                                                                 \
Entity BulkData::get_entity(type_name rank_id) const                   \
{                                                                      \
  ThrowAssert(m_sync_state != MODIFIABLE);                             \
  return (*m_bucket_repository.get_bucket(rank_id))[EXTRACT_BUCKET_ORDINAL(rank_id)]; \
}                                                                      \
                                                                       \
inline                                                                 \
type_name BulkData::get_##rank_name(Entity entity) const             \
{                                                                      \
  ThrowAssert(m_sync_state != MODIFIABLE);                             \
  ThrowAssert(is_valid(entity));                                       \
  ThrowAssert(this == &stk::mesh::BulkData::get(entity));              \
  ThrowAssert(bucket_ptr(entity));                                     \
  const MeshIndex &mesh_idx = mesh_index(entity);                      \
  const uint64_t bucket_id  = mesh_idx.bucket->bucket_id();            \
  const uint64_t bucket_ord = mesh_idx.bucket_ordinal;                 \
  return make_rank_id<type_name>(bucket_id,bucket_ord);                \
}                                                                      \
                                                                       \
inline                                                                 \
unsigned BulkData::num_##rank_name##s(Entity entity) const                  \
{                                                                      \
  ThrowAssert(is_valid(entity));                                       \
  ThrowAssert(bucket_ptr(entity));                                     \
  ThrowAssert(this == &stk::mesh::BulkData::get(entity));              \
  const MeshIndex &mesh_idx = mesh_index(entity);                      \
  const Bucket &b = *mesh_idx.bucket;                                  \
  unsigned bucket_ord = mesh_idx.bucket_ordinal;                       \
  return b.num_##rank_name##s(bucket_ord);                             \
}

//
// Define method for runtime rank
//
#define FUNCTION_DEF(begin_str, end_str, return_type)   \
                                                        \
inline                                                  \
return_type const* BulkData::begin_str##_##end_str(Entity entity, EntityRank rank) const \
{                                                                       \
  ThrowAssert(is_valid(entity));                                        \
  ThrowAssert(bucket_ptr(entity));                                      \
  ThrowAssert(this == &stk::mesh::BulkData::get(entity));               \
  const MeshIndex &mesh_idx = mesh_index(entity);                       \
  const Bucket &b = *mesh_idx.bucket;                                   \
  unsigned bucket_ord = mesh_idx.bucket_ordinal;                        \
  if (m_check_invalid_rels) {                                           \
    ThrowAssert(count_valid_connectivity(entity, rank) == b.num_connectivity(bucket_ord, rank)); \
  }                                                                     \
  return b.begin_str##_##end_str(bucket_ord, rank);                 \
}

#define FAST_RANK_FUNCTION_DEFS(rank_name, type_name)                  \
                                                                       \
inline                                                                 \
type_name const* BulkData::begin_##rank_name##s(Entity entity) const    \
{                                                                       \
  ThrowAssert(m_sync_state != MODIFIABLE);                              \
  ThrowAssert(is_valid(entity));                                        \
  ThrowAssert(bucket_ptr(entity));                                      \
  ThrowAssert(this == &stk::mesh::BulkData::get(entity));               \
  const MeshIndex &mesh_idx = mesh_index(entity);                       \
  const Bucket &b = *mesh_idx.bucket;                                   \
  unsigned bucket_ord = mesh_idx.bucket_ordinal;                        \
  return b.begin_##rank_name##s(bucket_ord);                            \
}                                                                       \
                                                                        \
inline                                                                  \
type_name const* BulkData::end_##rank_name##s(Entity entity) const      \
{                                                                       \
  ThrowAssert(m_sync_state != MODIFIABLE);                              \
  ThrowAssert(is_valid(entity));                                        \
  ThrowAssert(bucket_ptr(entity));                                      \
  ThrowAssert(this == &stk::mesh::BulkData::get(entity));               \
  const MeshIndex &mesh_idx = mesh_index(entity);                       \
  const Bucket &b = *mesh_idx.bucket;                                   \
  unsigned bucket_ord = mesh_idx.bucket_ordinal;                        \
  return b.end_##rank_name##s(bucket_ord);                              \
}

//
// Methods defined here!
//

RANK_FUNCTION_DEFS(node, Node);
RANK_FUNCTION_DEFS(edge, Edge);
RANK_FUNCTION_DEFS(face, Face);
RANK_FUNCTION_DEFS(element, Element);

FAST_RANK_FUNCTION_DEFS(node, Node);
FAST_RANK_FUNCTION_DEFS(edge, Edge);
FAST_RANK_FUNCTION_DEFS(face, Face);
FAST_RANK_FUNCTION_DEFS(element, Element);

FUNCTION_DEF(begin, entities, Entity);
FUNCTION_DEF(end,   entities, Entity);
FUNCTION_DEF(begin, ordinals, ConnectivityOrdinal);
FUNCTION_DEF(end,   ordinals, ConnectivityOrdinal);
FUNCTION_DEF(begin, permutations, Permutation);
FUNCTION_DEF(end,   permutations, Permutation);

#undef FUNCTION_DEFS
#undef RANK_FUNCTION_DEFS
#undef FAST_RANK_FUNCTION_DEFS
#undef BEGIN_END_PAIR

inline
unsigned BulkData::num_connectivity(Entity entity, EntityRank rank) const
{
  ThrowAssert(bucket_ptr(entity));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->num_connectivity(mesh_idx.bucket_ordinal, rank);
}

inline
unsigned BulkData::find_ordinal(Entity entity, EntityRank rank, ConnectivityOrdinal ordinal) const
{
  ThrowAssert(bucket_ptr(entity));
  const MeshIndex &mesh_idx = mesh_index(entity);
  unsigned num_rels = mesh_idx.bucket->num_connectivity(mesh_idx.bucket_ordinal, rank);
  ConnectivityOrdinal const *ords = mesh_idx.bucket->begin_ordinals(mesh_idx.bucket_ordinal, rank);

  unsigned i = 0;
  for (; i < num_rels; ++i)
  {
    if (ords[i] == ordinal)
      break;
  }
  return i;
}

inline
bool BulkData::has_permutation(Entity entity, EntityRank rank) const
{
  ThrowAssert(bucket_ptr(entity));
  return bucket(entity).has_permutation(rank);
}

/** \} */

inline
void BulkData::internal_basic_part_check(const Part* part,
                                         const unsigned ent_rank,
                                         const unsigned undef_rank,
                                         bool& intersection_ok,
                                         bool& rel_target_ok,
                                         bool& rank_ok) const
{
  // const unsigned part_rank = part->primary_entity_rank();

  intersection_ok = true;
  rel_target_ok   = true;

  // Do we allow arbitrary part changes to entities regardless of part rank? For the sake of the migration, we will for now.
#ifdef SIERRA_MIGRATION
  rank_ok = true;
#else
  const unsigned part_rank = part->primary_entity_rank();
  rank_ok         = ( ent_rank == part_rank ||
                      undef_rank  == part_rank );
#endif
}

inline bool BulkData::internal_quick_verify_change_part(const Part* part,
                                                        const unsigned ent_rank,
                                                        const unsigned undef_rank) const
{
  bool intersection_ok, rel_target_ok, rank_ok;
  internal_basic_part_check(part, ent_rank, undef_rank, intersection_ok, rel_target_ok, rank_ok);
  return intersection_ok && rel_target_ok && rank_ok;
}

inline
int BulkData::entity_comm_owner(const EntityKey & key) const
{
  const int owner_rank = m_entity_comm_map.owner_rank(key);
  ThrowAssertMsg(owner_rank == InvalidProcessRank || owner_rank == parallel_owner_rank(get_entity(key)),
                 "Expect entity " << key.id() << " to have owner " <<
                 parallel_owner_rank(get_entity(key)) << " but in comm map, found " << owner_rank);
  return owner_rank;
}

inline
bool BulkData::in_receive_ghost( EntityKey key ) const
{
  // Ghost communication with owner.
  const int owner_rank = entity_comm_owner(key);
  PairIterEntityComm ec = entity_comm(key);
  return !ec.empty() && ec.front().ghost_id != 0 &&
         ec.front().proc == owner_rank;
}

inline
bool BulkData::in_receive_ghost( const Ghosting & ghost , EntityKey key ) const
{
  const int owner_rank = entity_comm_owner(key);
  return in_ghost( ghost , key , owner_rank );
}

inline
bool BulkData::in_send_ghost( EntityKey key) const
{
  // Ghost communication with non-owner.
  const int owner_rank = entity_comm_owner(key);
  PairIterEntityComm ec = entity_comm(key);
  return ! ec.empty() && ec.back().ghost_id != 0 &&
    ec.back().proc != owner_rank;
}

class LessRelation
{
public:

  LessRelation(const BulkData &mesh) : m_mesh(mesh) { }

  inline bool operator() ( const Relation & lhs , const Relation & rhs ) const
  {
    bool result = false;

    // In Sierra, relations are sorted by RelationType in addition to Rank, Identifier, and target entity key.
  #ifdef SIERRA_MIGRATION
    if (lhs.entity_rank() != rhs.entity_rank()) {
      result = lhs.entity_rank() < rhs.entity_rank();
    }
    else if (lhs.getRelationType() != rhs.getRelationType()) {
      result = lhs.getRelationType() < rhs.getRelationType();
    }
    else if (lhs.relation_ordinal() != rhs.relation_ordinal()) {
      result = lhs.relation_ordinal() < rhs.relation_ordinal();
    }
  #else
    if ( lhs.m_raw_relation.value != rhs.m_raw_relation.value ) {
      result = lhs.m_raw_relation.value < rhs.m_raw_relation.value ;
    }
  #endif
    else {
      Entity lhs_entity = lhs.entity();
      const size_t lhs_offset = m_mesh.is_valid(lhs_entity) ? lhs_entity.local_offset() : Entity::MaxEntity;
      Entity rhs_entity = rhs.entity();
      const size_t rhs_offset = m_mesh.is_valid(rhs_entity) ? rhs_entity.local_offset() : Entity::MaxEntity;
      result = lhs_offset < rhs_offset;
    }
    return result ;
  }

  bool operator() ( const Relation & lhs , Relation::raw_relation_id_type rhs ) const
    { return lhs.raw_relation_id() < rhs ; }

private:

  const BulkData &m_mesh;

  LessRelation();
};

inline
Relation::Relation(const BulkData &mesh,  Entity ent , RelationIdentifier id )
  : m_raw_relation( Relation::raw_relation_id( mesh.entity_rank(ent) , id ) ),
    m_target_entity(ent)
{
#ifdef SIERRA_MIGRATION
  setRelationType(RelationType::INVALID);
#endif
}

inline
bool Bucket::other_entities_have_single_rank(unsigned bucket_ordinal, EntityRank rank) const
{
  Entity const * const other_rels = m_dynamic_other_connectivity.begin_entities(bucket_ordinal);
  Entity const * const other_rels_end = m_dynamic_other_connectivity.end_entities(bucket_ordinal);;

  if ((other_rels == other_rels_end) && (rank != InvalidEntityRank))
    return false;

  return (m_mesh.entity_rank(*other_rels) == rank) && (m_mesh.entity_rank(*(other_rels_end - 1)) == rank);
}

} // namespace mesh
} // namespace stk

#include <stk_mesh/baseImpl/Partition.hpp>
#include <stk_mesh/base/BulkDataPartOperations.tcc>

#endif //  stk_mesh_BulkData_hpp
