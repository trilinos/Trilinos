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
#include <stk_util/parallel/DistributedIndex.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Ghosting.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Trace.hpp>
#include <stk_mesh/base/EntityComm.hpp>

#include <stk_mesh/baseImpl/EntityRepository.hpp>
#include <stk_mesh/baseImpl/BucketRepository.hpp>

#include <algorithm>
#include <map>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

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
public:

  void optimize_buckets_at_modification_end(bool b) { m_optimize_buckets = b; }

  inline static BulkData & get( const Bucket & bucket);
  inline static BulkData & get( const Entity & entity);
  inline static BulkData & get( const Ghosting & ghost);

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
      bool use_memory_pool = true );

  //------------------------------------
  /** \brief  The meta data manager for this bulk data manager. */
  const MetaData & mesh_meta_data() const { return m_mesh_meta_data ; }

  /** \brief  The parallel machine */
  ParallelMachine parallel() const { return m_parallel_machine ; }

  /** \brief  Size of the parallel machine */
  unsigned parallel_size()   const { return m_parallel_size ; }

  /** \brief  Rank of the parallel machine's local processor */
  unsigned parallel_rank()   const { return m_parallel_rank ; }

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
  void update_field_data_states() const { m_bucket_repository.update_field_data_states(); }

  /** \brief  Copy field data from src entity to Dest entity
   *           - Fields that exist on the src that don't exist on the dest will
   *             be ignored
   *           - Fields that exist on the dest that don't exist on the src will
   *             be zeroed or initialized with the Field-specified initial-value.
   */
  void copy_entity_fields( const Entity & src, Entity & dest) {
    //TODO fix const correctness for src
    Entity & non_const_src = const_cast<Entity &>(src);
    m_bucket_repository.copy_fields(
        dest.bucket(),
        dest.bucket_ordinal(),
        non_const_src.bucket(),
        non_const_src.bucket_ordinal()
        );
  }

  //------------------------------------
  /** \brief  Query all buckets of a given entity rank */
  const std::vector<Bucket*> & buckets( EntityRank rank ) const
  { return m_bucket_repository.buckets(rank); }

#ifndef SWIG //SRK
  AllBucketsRange get_bucket_range() const
  {
    return m_bucket_repository.get_bucket_range();
  }

  AllBucketsRange get_bucket_range( EntityRank rank ) const
  {
    return m_bucket_repository.get_bucket_range(rank);
  }
#endif

  /** \brief  Get entity with a given key */
  Entity * get_entity( EntityRank entity_rank , EntityId entity_id ) const {
    require_good_rank_and_id(entity_rank, entity_id);
    return m_entity_repo.get_entity( EntityKey(entity_rank, entity_id));
  }

  /** \brief  Get entity with a given key */
  Entity * get_entity( const EntityKey key ) const  {
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
  Entity & declare_entity( EntityRank ent_rank ,
      EntityId ent_id , const PartVector& parts);

  void change_entity_id( EntityId id, Entity & entity);

  /** \brief  Change the parallel-locally-owned entity's
   *          part membership by adding and/or removing parts
   *
   *  A parallel-local operation.
   *
   *  If the locally owned entity is shared or ghosted then
   *  the change will be propogated to the sharing or ghosting
   *  processes by modification_end.
   */
  void change_entity_parts( Entity & entity,
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
  void change_entity_parts( Entity & entity,
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
  bool destroy_entity( Entity * & entity );

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
      std::vector<Entity *>& requested_entities);

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
   *     part.primary_entity_rank() == e_from.entity_rank()
   *     then 'e_to' has induced membership in 'part'.
   *  2) If there exists a part relation 'part_rel' such that
   *     'e_from' is a member of part_rel.m_root and
   *     the entity relation conforms to the part relation
   *     then 'e_to' has induced membership in part_rel.m_target.
   *
   * Note that relation-declarations must be symmetric across all
   * sharers of the involved entities within a modification cycle.
   */
  void declare_relation( Entity & e_from ,
      Entity & e_to ,
      const RelationIdentifier local_id );

  /** \brief  Declare a collection of relations by simply iterating
   *          the input and calling declare_relation on each entry.
   */
  void declare_relation( Entity & entity, const std::vector<Relation> & rel);

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
  bool destroy_relation( Entity & e_from ,
                         Entity & e_to,
                         const RelationIdentifier local_id );

  //------------------------------------
  //------------------------------------
  /** \brief  All entities with communication information. */
  const std::vector<Entity*> & entity_comm() const
    { return m_entity_comm ; }

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
                        const std::vector<Entity*> & remove_receive );

  /** \brief  Empty every single Ghosting.
   *          Same result, but more efficient than, calling
   *          change_ghosting to remove every single ghosted entity.
   */
  void destroy_all_ghosting();

  /** \brief  Query all ghostings */
  const std::vector<Ghosting*> & ghostings() const { return m_ghosting ; }

  /** \brief  Entity Comm functions that are now moved to BulkData
   * These functions are only here for backwards compatibility.
   * We plan to remove all comm accessors from the Entity in the future.
   */
  PairIterEntityComm entity_comm(const EntityKey & key) const { return m_entity_comm_map.comm(key); }
  PairIterEntityComm entity_comm_sharing(const EntityKey & key) const { return m_entity_comm_map.sharing(key); }
  PairIterEntityComm entity_comm(const EntityKey & key, const Ghosting & sub ) const { return m_entity_comm_map.comm(key,sub); }
  bool entity_comm_insert( const EntityKey & key, const EntityCommInfo & val) { return m_entity_comm_map.insert(key,val); }
  bool entity_comm_erase(  const EntityKey & key, const EntityCommInfo & val) { return m_entity_comm_map.erase(key,val); }
  bool entity_comm_erase(  const EntityKey & key, const Ghosting & ghost) { return m_entity_comm_map.erase(key,ghost); }
  void entity_comm_clear_ghosting(const EntityKey & key ) { m_entity_comm_map.comm_clear_ghosting(key); }
  void entity_comm_clear(const EntityKey & key) { m_entity_comm_map.comm_clear(key); }
  void entity_comm_swap(const EntityKey & key1, const EntityKey & key2) { m_entity_comm_map.comm_swap(key1, key2); }

private:

  /** \brief  The meta data manager for this bulk data manager. */
  MetaData & meta_data() const { return m_mesh_meta_data ; }

#ifndef DOXYGEN_COMPILE

  BulkData();
  BulkData( const BulkData & );
  BulkData & operator = ( const BulkData & );

  /** \brief  Parallel index for entity keys */
  parallel::DistributedIndex          m_entities_index ;
  impl::EntityRepository              m_entity_repo ;
  impl::BucketRepository              m_bucket_repository ;
  std::vector<Entity*>                m_entity_comm ;
  std::vector<Ghosting*>              m_ghosting ; /**< Aura is [1] */

  // Other information:
  MetaData &         m_mesh_meta_data ;
  ParallelMachine    m_parallel_machine ;
  unsigned           m_parallel_size ;
  unsigned           m_parallel_rank ;
  size_t             m_sync_count ;
  BulkDataSyncState  m_sync_state ;
  bool               m_meta_data_verified ;
  bool               m_optimize_buckets;
  bool               m_mesh_finalized;
  EntityComm         m_entity_comm_map;

  /**
   * For all processors sharing an entity, find one to be the new
   * owner.
   */
  unsigned determine_new_owner( Entity & ) const ;

  /*  Entity modification consequences:
   *  1) Change entity relation => update via part relation => change parts
   *  2) Change parts => update forward relations via part relation
   *                  => update via field relation
   */
  void internal_change_entity_parts( Entity & ,
                                     const PartVector & add_parts ,
                                     const PartVector & remove_parts );

//Optional parameter 'always_propagate_internal_changes' is always true except when this function
//is being called from the sierra-framework. The fmwk redundantly does its own propagation of the
//internal part changes (mostly induced-part stuff), so it's a performance optimization to avoid
//the propagation that stk-mesh does.
  void internal_change_entity_parts( Entity & ,
                                     const OrdinalVector & add_parts ,
                                     const OrdinalVector & remove_parts,
                                     bool always_propagate_internal_changes=true);

  void internal_propagate_part_changes( Entity & entity, const PartVector & removed );
  void internal_propagate_part_changes( Entity & entity, const OrdinalVector & removed );

  void internal_change_ghosting( Ghosting & ghosts,
                                 const std::vector<EntityProc> & add_send ,
                                 const std::vector<Entity*> & remove_receive );

  bool internal_modification_end( bool regenerate_aura );
    void internal_resolve_shared_modify_delete();
    void internal_resolve_shared_modify_delete_second_pass();
    void internal_resolve_ghosted_modify_delete();
    void internal_resolve_parallel_create();
    void internal_resolve_shared_membership();

  void internal_update_distributed_index( std::vector<Entity*> & shared_new );

  /** \brief  Regenerate the shared-entity aura,
   *          adding and removing ghosted entities as necessary.
   *
   *  - a collective parallel operation.
   */
  void internal_regenerate_shared_aura();

  void internal_basic_part_check(const Part* part,
                                 const unsigned entity_rank,
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
                                         const unsigned entity_rank,
                                         const unsigned undef_rank) const;

  void internal_verify_change_parts( const MetaData   & meta ,
                                     const Entity     & entity ,
                                     const PartVector & parts ) const;

  void internal_verify_change_parts( const MetaData   & meta ,
                                     const Entity     & entity ,
                                     const OrdinalVector & parts ) const;

  //------------------------------------

  /** \name  Invariants/preconditions for MetaData.
   * \{
   */

  /** \brief  All non-const methods assert this */
  void require_ok_to_modify() const ;

  void require_entity_owner( const Entity & entity, unsigned owner) const ;

  void require_metadata_committed() const;

  void require_good_rank_and_id(EntityRank ent_rank, EntityId ent_id) const;

  void require_valid_relation( const char action[] ,
                               const BulkData & mesh ,
                               const Entity   & e_from ,
                               const Entity   & e_to );

  /** \} */

  //------------------------------------

  // FIXME: Remove this friend once unit-testing has been refactored
  friend class UnitTestModificationEndWrapper;
  friend class ::stk::mesh::MetaData;
#endif /* DOXYGEN_COMPILE */
};

BulkData & BulkData::get( const Bucket & bucket) {
  return bucket.bulk_data();
}

BulkData & BulkData::get( const Entity & entity) {
  return BulkData::get(entity.bucket());
}

BulkData & BulkData::get( const Ghosting & ghost) {
  return ghost.bulk_data();
}

/** \} */

/** \brief Free function for setting relations on entities.
 * \relates BulkData
 * \relates BucketRepository
 */
// TODO - Does this need to be in the public API? It's only used internally.
// In what context would a client ever call this?
void set_field_relations( Entity & e_from ,
                          Entity & e_to ,
                          const unsigned ident );

inline
void BulkData::internal_basic_part_check(const Part* part,
                                         const unsigned entity_rank,
                                         const unsigned undef_rank,
                                         bool& intersection_ok,
                                         bool& rel_target_ok,
                                         bool& rank_ok) const
{
  // const unsigned part_rank = part->primary_entity_rank();

  intersection_ok = part->intersection_of().empty();
  rel_target_ok   = ( part->relations().empty() ||
                      part != part->relations().begin()->m_target );
  // Do we allow arbitrary part changes to entities regardless of part rank? For the sake of the migration, we will for now.
#ifdef SIERRA_MIGRATION
  rank_ok = true;
#else
  const unsigned part_rank = part->primary_entity_rank();
  rank_ok         = ( entity_rank == part_rank ||
                      undef_rank  == part_rank );
#endif
}

inline bool BulkData::internal_quick_verify_change_part(const Part* part,
                                                        const unsigned entity_rank,
                                                        const unsigned undef_rank) const
{
  bool intersection_ok, rel_target_ok, rank_ok;
  internal_basic_part_check(part, entity_rank, undef_rank, intersection_ok, rel_target_ok, rank_ok);
  return intersection_ok && rel_target_ok && rank_ok;
}

//Optional parameter 'always_propagate_internal_changes' is always true except when this function
//is being called from the sierra-framework. The fmwk redundantly does its own propagation of the
//internal part changes (mostly induced-part stuff), so it's a performance optimization to avoid
//the propagation that stk-mesh does.
template<typename AddIterator, typename RemoveIterator>
void BulkData::change_entity_parts( Entity & entity,
                                    AddIterator begin_add_parts, AddIterator end_add_parts,
                                    RemoveIterator begin_remove_parts, RemoveIterator end_remove_parts,
                                    bool always_propagate_internal_changes)
{
  TraceIfWatching("stk::mesh::BulkData::change_entity_parts", LOG_ENTITY, entity.key());
  DiagIfWatching(LOG_ENTITY, entity.key(), "entity state: " << entity);
  //DiagIfWatching(LOG_ENTITY, entity.key(), "add_parts: " << add_parts);
  //DiagIfWatching(LOG_ENTITY, entity.key(), "remove_parts: " << remove_parts);

  require_ok_to_modify();

  require_entity_owner( entity , m_parallel_rank );

  const EntityRank entity_rank = entity.entity_rank();
  const EntityRank undef_rank  = InvalidEntityRank;

  // Transitive addition and removal:
  // 1) Include supersets of add_parts
  // 2) Do not include a remove_part if it appears in the add_parts
  // 3) Include subsets of remove_parts

  // most parts will at least have universal and topology part as supersets
  const unsigned expected_min_num_supersets = 2;

  OrdinalVector a_parts;
  a_parts.reserve( std::distance(begin_add_parts, end_add_parts) * (expected_min_num_supersets + 1) );
  for(AddIterator add_iter=begin_add_parts; add_iter!=end_add_parts; ++add_iter) {
    a_parts.push_back((*add_iter)->mesh_meta_data_ordinal());
  }
  bool quick_verify_check = true;

  for ( AddIterator ia = begin_add_parts; ia != end_add_parts ; ++ia ) {
    quick_verify_check = quick_verify_check &&
      internal_quick_verify_change_part(*ia, entity_rank, undef_rank);
    const PartVector& supersets = (*ia)->supersets();
    for(PartVector::const_iterator s_iter=supersets.begin(), s_end=supersets.end();
        s_iter!=s_end; ++s_iter) {
      a_parts.push_back((*s_iter)->mesh_meta_data_ordinal());
    }
  }

  order(a_parts);

  OrdinalVector::const_iterator a_parts_begin = a_parts.begin(),
                                a_parts_end   = a_parts.end();
  OrdinalVector r_parts ;

  for ( RemoveIterator ir = begin_remove_parts; ir != end_remove_parts ; ++ir ) {

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
      internal_quick_verify_change_part(*ir, entity_rank, undef_rank);

    if ( ! contains_ordinal( a_parts_begin, a_parts_end , (*ir)->mesh_meta_data_ordinal() ) ) {
      r_parts.push_back( (*ir)->mesh_meta_data_ordinal() );
      for ( PartVector::const_iterator  cur_part = (*ir)->subsets().begin() ;
            cur_part != (*ir)->subsets().end() ;
            ++cur_part )
        if ( entity.bucket().member ( **cur_part ) )
          r_parts.push_back ( (*cur_part)->mesh_meta_data_ordinal() );
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

  return ;
}

} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif //  stk_mesh_BulkData_hpp
