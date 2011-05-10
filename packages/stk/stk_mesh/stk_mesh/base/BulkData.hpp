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
#include <stk_mesh/baseImpl/EntityRepository.hpp>
#include <stk_mesh/base/Ghosting.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/baseImpl/BucketRepository.hpp>

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
      unsigned bucket_max_size = 1000 );

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
  /** \brief  Query the upper bound on the number of mesh entities
   *         that may be associated with a single bucket.
   */
  unsigned bucket_capacity() const { return m_bucket_repository.bucket_capacity(); }

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
   *             be zeroed
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
      const PartVector & remove_parts = PartVector() );

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
   */
  void destroy_relation( Entity & e_from ,
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
  MetaData &   m_mesh_meta_data ;
  ParallelMachine    m_parallel_machine ;
  unsigned           m_parallel_size ;
  unsigned           m_parallel_rank ;
  size_t             m_sync_count ;
  BulkDataSyncState  m_sync_state ;
  bool               m_meta_data_verified ;

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

  void internal_propagate_part_changes( Entity & entity, const PartVector & removed );

  void internal_change_ghosting( Ghosting & ghosts,
                                 const std::vector<EntityProc> & add_send ,
                                 const std::vector<Entity*> & remove_receive );

  bool internal_modification_end( bool regenerate_aura );
    void internal_resolve_shared_modify_delete();
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

} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif //  stk_mesh_BulkData_hpp
