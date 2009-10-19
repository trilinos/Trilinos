
#ifndef stk_mesh_BulkData_hpp
#define stk_mesh_BulkData_hpp

//----------------------------------------------------------------------

#include <map>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Ghosting.hpp>
#include <stk_mesh/base/Selector.hpp>

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

  ~BulkData();

  /** \brief  Construct mesh bulk data manager conformal to the given
   *          \ref stk::mesh::MetaData "meta data manager" and will
   *          distribute bulk data over the given parallel machine.
   *
   *  - The meta data manager must be committed.
   *  - The maximum number of entities per bucket may be supplied.
   *  - The bulk data is in the "ok to modify" state.
   */
  BulkData( const MetaData & mesh_meta_data ,
            ParallelMachine parallel ,
            unsigned bucket_max_size = 1000);

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
  unsigned bucket_capacity() const { return m_bucket_capacity ; }

  //------------------------------------
  /** \brief  Bulk data has two states:
   *          guaranteed to be parallel synchronized or
   *          modification in progress and may be parallel inconsistent.
   */
  bool synchronized_state() const { return m_sync_state ; }

  /** \brief  Count of the number of times that the bulk data has been
   *          parallel synchronized.  This count gets updated with
   *          each call to 'modification_end'.
   */
  size_t synchronized_count() const { return m_sync_count ; }

  /** \brief  Begin a modification phase during which the mesh bulk data
   *          could become parallel inconsistent.  This is a parallel
   *          synchronous call.
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
  void change_entity_owner( const std::vector<EntityProc> & );

  //------------------------------------
  
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
  void update_field_data_states() const ;

  //------------------------------------
  /** \brief  Query all buckets of a given entity type */
  const std::vector<Bucket*> & buckets( unsigned type ) const ;

  /** \brief  Get entity with a given key */
  /// \todo REFACTOR remove required_by argument
  Entity * get_entity( EntityType ent_type , EntityId ent_id ,
                       const char * required_by = NULL ) const ;

  //------------------------------------
  /** \brief  Create or retrieve a locally owned entity of a
   *          given type and id.
   *
   *  A parallel-local operation.
   *
   *  The entity is created as locally owned and a member of the input
   *  mesh parts.  The entity a member of the meta data's locally owned
   *  mesh part and the entity's owner_rank() == parallel_rank().
   *
   *  If two or more processes create an entity of the same type
   *  and identifier then the sharing and ownership of these entities
   *  will be resolved by the call to 'modification_end'.
   */
  Entity & declare_entity( EntityType ent_type ,
                           EntityId ent_id , const std::vector<Part*> & );

  /** \brief  Change the parallel-locally-owned entity's
   *          part membership by adding and/or removing parts
   *
   *  A parallel-local operation.
   *
   *  If the locally owned entity is shared or ghosted then
   *  the change will be propogated to the sharing or ghosting
   *  processes by modification_end.
   */
  void change_entity_parts( Entity & ,
                            const std::vector<Part*> & add_parts ,
                            const std::vector<Part*> & remove_parts = 
                                  std::vector<Part*>() );

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
   * \paragraph destroy_locally_used  Destroy Locally Used
   *
   *  Entities in the 'locally_used_part' but not in the
   *  'locally_owned_part' are shared.  These entities are deleted
   *  on the local process and removed from the sharing lists on
   *  other processes during modication_end.
   *
   * \paragraph destroy_ghosted  Destroy Ghosted
   *
   *  Entities not in the 'locally_used_part' are ghosted.
   *  These entities are removed from all ghosting lists
   *  during 'modification_end'.
   *  processes by modification_end.
   *
   *  \return  True if the request for destruction is accepted; i.e.,
   *           if the entity is not the 'to' member of a relation.
   */
  bool destroy_entity( Entity * & );

  //------------------------------------
  /** \brief  Declare a relation and its converse between
   *          entities in the same mesh.
   *
   *  A parallel-local mesh modificaton operation.
   *
   *  This mapping ( e_from , local_id , kind ) -> e_to  must be unique.
   *
   *  Relations between entities induces part membership as follows.
   *  1) If 'e_from' is a member of 'part' and
   *     part.primary_entity_type() == e_from.entity_type()
   *     then 'e_to' has induced membership in 'part'.
   *  2) If there exists a part relation 'part_rel' such that
   *     'e_from' is a member of part_rel.m_root and 
   *     the entity relation conforms to the part relation
   *     then 'e_to' has induced membership in part_rel.m_target.
   */
  void declare_relation( Entity & e_from ,
                         Entity & e_to ,
                         const unsigned local_id ,
                         const unsigned kind = 0 );

  /** \brief  Declare a collection of relations by simply iterating
   *          the input and calling declare_relation on each entry.
   */
  void declare_relation( Entity & , const std::vector<Relation> & );

  /** \brief  Remove all relations of a given kind between two entities.
   *
   *  If the relation induced a part membership for 'e_to' and 'e_to'
   *  is not shared with another processor then that part membership
   *  is removed if and only if there does not exist another relation
   *  inducing the same part membership.
   *  If 'e_to' is shared then the check for removing the induced
   *  relatinship does not occur for that entity until the call to
   *  'modification_end'.
   */
  void destroy_relation( Entity & e_from , Entity & e_to , unsigned kind = 0 );

  //------------------------------------
  //------------------------------------
  /** \brief Symmetric parallel relations for shared mesh entities.  */
  const std::vector<EntityProc> & shared_entities() const
    { return m_shares_all ; }

  //------------------------------------
  /** \brief  Query the shared-entity aura.
   *          Is likely to be stale if ownership or sharing has changed
   *          and the 'modification_end' has not been called.
   */
  Ghosting & shared_aura() const { return * m_ghosting[0] ; }

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
  void change_ghosting( Ghosting & ,
                        const std::vector<EntityProc> & add_send ,
                        const std::vector<Entity*> & remove_receive );

  /** \brief  Empty every single Ghosting.
   *          Same result, but more efficient than, calling
   *          change_ghosting to remove every single ghosted entity.
   */
  void destroy_all_ghosting();

  /** \brief  Query all ghostings */
  const std::vector<Ghosting*> & ghostings() const { return m_ghosting ; }

  //------------------------------------
  /** \brief  All non-const methods assert this */
  void assert_ok_to_modify( const char * ) const ;

  void assert_entity_owner( const char * , const Entity & , unsigned ) const ;

  void assert_good_key( const char * , const EntityKey & ) const ;

  unsigned parallel_index_rank( const EntityKey & ) const ;

  //------------------------------------
private:
  void verify_type_and_id(const char* calling_method,
                          EntityType ent_type, EntityId ent_id) const;

#ifndef DOXYGEN_COMPILE

  typedef std::map<EntityKey,Entity*> EntitySet ;

  BulkData();
  BulkData( const BulkData & );
  BulkData & operator = ( const BulkData & );

  // Containers:
  std::vector< std::vector<Bucket*> > m_buckets ;
  EntitySet                           m_entities ;
  std::vector<EntityProc>             m_shares_all ;
  std::vector<Ghosting*>              m_ghosting ; /**< Aura is [0] */

  std::vector< std::pair< EntityKey , unsigned > > m_entities_owner_index ;
  std::vector<Entity*>                m_new_entities ;
  std::vector<Entity*>                m_del_entities ;
  Bucket *                            m_bucket_nil ;

  // Other information:
  const MetaData & m_mesh_meta_data ;
  ParallelMachine  m_parallel_machine ;
  unsigned         m_parallel_size ;
  unsigned         m_parallel_rank ;
  unsigned         m_bucket_capacity ;
  size_t           m_sync_count ;
  bool             m_sync_state ;

  /*  Entity modification consequences:
   *  1) Change entity relation => update via part relation => change parts
   *  2) Change parts => update forward relations via part relation
   *                  => update via field relation
   */
  void remove_entity( Bucket * , unsigned );

  std::pair<Entity*,bool>
   internal_create_entity( const EntityKey & key , const unsigned owner );

  void internal_destroy_entity( Entity * );

  void internal_change_entity_parts( Entity & ,
                                     const PartVector & add_parts ,
                                     const PartVector & remove_parts );

  void internal_propagate_part_changes( Entity & , const PartVector & removed );

  void internal_propagate_relocation( Entity & );

  void internal_change_ghosting( Ghosting & ,
                                 const std::vector<EntityProc> & add_send ,
                                 const std::vector<Entity*> & remove_receive );

  void internal_regenerate_sharing();

  void internal_set_shared_entities();
 
  void internal_resolve_destroy_ghosted(
    const std::vector<Entity*>    & del_local ,
    const std::vector<EntityProc> & del_remote );

  void internal_resolve_destroy_shared(
    const std::vector<Entity*>    & del_local ,
    const std::vector<EntityProc> & del_remote ,
          std::vector<EntityProc> & owner_change );

  void internal_update_parallel_index(
    const std::vector<Entity*>    & new_entities ,
    const std::vector<EntityProc> & owner_change ,
          std::vector<EntityProc> & shared_new );

  void internal_resolve_created_shared(
    const std::vector<EntityProc> & new_shared );

  void internal_resolve_parallel_create_delete();

  void internal_resolve_shared_membership();

  void internal_sort_bucket_entities();

  /** \brief  Put owned entity in send list for off-process
   *          parallel index, shared, and ghosted.
   */
  void owner_send_to_all( Entity * , std::vector<EntityProc> & ) const ;

  /** \brief  Regenerate the shared-entity aura,
   *          adding and removing ghosted entities as necessary.
   *
   *  - a collective parallel operation.
   */
  void internal_regenerate_shared_aura();

  void destroy_shared_aura();

  friend class UnitTestBulkData ;
#endif /* DOXYGEN_COMPILE */

};

/** \} */

} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif //  stk_mesh_BulkData_hpp
