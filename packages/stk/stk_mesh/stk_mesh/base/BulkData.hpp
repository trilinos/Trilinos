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

#ifndef stk_mesh_BulkData_hpp
#define stk_mesh_BulkData_hpp

//----------------------------------------------------------------------

#include <stddef.h>                     // for size_t, NULL
#include <stdint.h>                     // for uint16_t
#include <algorithm>                    // for max
#include <functional>                   // for less, equal_to
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <list>                         // for list
#include <map>                          // for map, map<>::value_compare
#include <set>
#include <stk_mesh/base/Entity.hpp>     // for Entity, etc
#include <stk_mesh/base/EntityCommDatabase.hpp>  // for EntityCommDatabase
#include <stk_mesh/base/Ghosting.hpp>   // for Ghosting
#include <stk_mesh/base/Selector.hpp>   // for Selector
#include <stk_mesh/base/Trace.hpp>      // for TraceIfWatching, etc
#include <stk_mesh/base/Types.hpp>      // for MeshIndex, EntityRank, etc
#include <stk_mesh/baseImpl/BucketRepository.hpp>  // for BucketRepository
#include <stk_mesh/baseImpl/EntityRepository.hpp>  // for EntityRepository, etc
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <stk_util/util/TrackingAllocator.hpp>  // for tracking_allocator, etc
#include <string>                       // for char_traits, string
#include <utility>                      // for pair
#include <vector>                       // for vector
#include "boost/functional/hash/extensions.hpp"  // for hash
#include "boost/tuple/detail/tuple_basic.hpp"  // for get
#include "stk_mesh/base/Bucket.hpp"     // for Bucket, Bucket::size_type, etc
#include "stk_mesh/base/BucketConnectivity.hpp"  // for BucketConnectivity
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey, hash_value
#include "stk_mesh/base/FieldDataManager.hpp"
#include "stk_mesh/base/Relation.hpp"   // for Relation, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/environment/ReportHandler.hpp"  // for ThrowAssert, etc

namespace stk { namespace mesh { class FieldBase; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { namespace impl { class Partition; } } }
namespace stk { namespace mesh { struct ConnectivityMap; } }
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { namespace impl { class EntityRepository; } } }

#include "EntityCommListInfo.hpp"
#include "EntityLess.hpp"
#include "StkDebuggingMacros.hpp"
#include "SharedEntityType.hpp"

namespace stk {
namespace mesh {

class BulkData {

public:
  enum GHOSTING_ID { SHARED = 0, AURA = 1 };
  enum BulkDataSyncState { MODIFIABLE = 1 , SYNCHRONIZED = 2 };
  enum modification_optimization {MOD_END_COMPRESS_AND_SORT , MOD_END_SORT };
  enum entitySharing { NOT_MARKED=0, POSSIBLY_SHARED=1, IS_SHARED=2 };
  enum ModificationTypes {
      CHANGE_ENTITY_OWNER,
      CHANGE_GHOSTING,
      CREATE_GHOSTING,
      DECLARE_RELATION,
      DESTROY_GHOSTING,
      DESTROY_ALL_GHOSTING,
      DESTROY_RELATION,
      NumModificationTypes
  };

  enum EntityModificationTypes {
      CHANGE_ENTITY_ID,
      CHANGE_ENTITY_PARTS,
      DECLARE_ENTITY,
      DESTROY_ENTITY,
      NumEntityModificationTypes
  };
  enum PublicOrInternalMethod {
      PUBLIC = 0,
      INTERNAL,
      NumMethodTypes
  };

  /** \brief  Construct mesh bulk data manager conformal to the given
   *          \ref stk::mesh::MetaData "meta data manager" and will
   *          distribute bulk data over the given parallel machine.
   *
   *  - The maximum number of entities per bucket may be supplied.
   *  - The bulk data is in the synchronized or "locked" state.
   */
  BulkData(   MetaData & mesh_meta_data
            , ParallelMachine parallel
#ifdef SIERRA_MIGRATION
            , bool add_fmwk_data = false
#endif
            , ConnectivityMap const* arg_connectivity_map = NULL
            , FieldDataManager *field_dataManager = NULL
            , unsigned bucket_capacity = impl::BucketRepository::default_bucket_capacity
            );

  virtual ~BulkData();

  inline static BulkData & get( const Bucket & bucket);
  inline static BulkData & get( const Ghosting & ghost);

  //------------------------------------
  /** \brief  The meta data manager for this bulk data manager. */
  const MetaData & mesh_meta_data() const { return m_mesh_meta_data ; }
        MetaData & mesh_meta_data()       { return m_mesh_meta_data ; }

  /** \brief  The parallel machine */
  ParallelMachine parallel() const { return m_parallel.parallel() ; }

  /** \brief  Size of the parallel machine */
  int parallel_size()   const { return m_parallel.parallel_size() ; }

  /** \brief  Rank of the parallel machine's local processor */
  int parallel_rank()   const { return m_parallel.parallel_rank() ; }

  const ConnectivityMap & connectivity_map() const { return m_bucket_repository.connectivity_map(); }

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
   *  Input string is some description of why or where modification was
   *  started from that can be later queried.
   *
   *  \return  True if transitioned out of the guaranteed
   *           parallel consistent state to the "ok to modify" state.
   *           False if already in this state.
   */
  bool modification_begin(const std::string description = std::string("UNSPECIFIED"));

  /** \brief  Return the description given for the modification_begin
   *          call that began the most recent mesh modification cycle
   */
  const std::string& get_modification_begin_description() const;


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
  bool modification_end( modification_optimization opt = MOD_END_SORT );
  bool modification_end_for_percept( modification_optimization opt = MOD_END_SORT);
  bool modification_end_for_entity_creation( EntityRank entity_rank, modification_optimization opt = MOD_END_SORT);

  /** \brief  Give away ownership of entities to other parallel processes.
   *
   *  A parallel-synchronous operation while the mesh is not being modified
   *
   *  Each owning process inputs a list of entities and the
   *  new owning process.  Upon completion of the call the owning
   *  processes have the newly owned entities as well as the closure
   *  of those entities (without changing the ownership of the closure
   *  entities).  If a previous owner no longer needs a
   *  changed-owner entity to support the closure of a still-owned
   *  entity then the changed-owner entity is deleted from that process.
   *  All ghosts of all entities affected by the changed ownerships
   *  are deleted.  This must be called outside a mesh modification cycle
   *  as it is a self-contained atomic modification operation.  This uses
   *  enough communication that it will be most efficient to batch up all
   *  desired changes so that it can be called only once.
   */
  void change_entity_owner( const EntityProcVec & arg_change,
                            bool regenerate_aura = true,
                            modification_optimization mod_optimization = MOD_END_SORT );

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
  inline void copy_entity_fields( Entity src, Entity dst);      // CLEANUP: only used in usecases and perf tests - REMOVE from class

  //------------------------------------
  /** \brief  Query all buckets of a given entity rank
   *  Don't call inside BucketRepository member functions!
   */
  const BucketVector & buckets( EntityRank rank ) const { return m_bucket_repository.buckets(rank); }

  //iterator that traverses entities of the specified rank, in order of ascending global identifier
  const_entity_iterator begin_entities(EntityRank ent_rank) const { return m_entity_repo.begin_rank(ent_rank);}

  //end-iterator for entities of the specified rank
  const_entity_iterator end_entities(EntityRank ent_rank) const { return m_entity_repo.end_rank(ent_rank); }

  /** \brief  Get entity with a given rank and id */
  inline Entity get_entity( EntityRank ent_rank , EntityId entity_id ) const;

  /** \brief  Get entity with a given key */
  Entity get_entity( const EntityKey key ) const  { return m_entity_repo.get_entity(key); }

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

  Entity declare_entity( EntityRank ent_rank , EntityId ent_id , Part& part);

  /** This overloading of declare_entity that doesn't take a part
   * creates the new entity in the 'universal' part.
   */
  Entity declare_entity( EntityRank ent_rank , EntityId ent_id);

  /** \brief Add sharing information about a newly-created node
   *
   * All nodes that are shared between processors must have the same
   * ID across all processors, and they must have this sharing information
   * explicitly provided in the same modification cycle as it was created,
   * on all sharing processors. Otherwise, the nodes will remain unique
   * on different processors.
   */
  void add_node_sharing( Entity node, int sharing_proc );

  /** \brief Change an entity's EntityId, for Framework implementation.
   *         Only intended for use in supporting Framework implementation.
   */
  void change_entity_id( EntityId id, Entity entity);   // Deprecated?

  /** \brief  Change the parallel-locally-owned entity's
   *          part membership by adding and/or removing parts
   *
   *  If the locally owned entity is shared or ghosted then
   *  the change will be propagated to the sharing or ghosting
   *  processes by modification_end.
   */
  void change_entity_parts( Entity entity,
      const PartVector & add_parts ,
      const PartVector & remove_parts = PartVector(),
      bool always_propagate_internal_changes=true);

  /** \brief Change part-membership of the specified entities by adding
   * and/or removing parts for each entity.
   *
   * Each entity in the vector must be locally-owned.
   *
   * A parallel-synchronous operation, containing its own calls to
   * modification_begin and modification_end. The mesh must not already
   * be in the "ok-to-modify" state when this method is called.
   *
   * The supplied arguments are different on each processor (and may be
   * empty on some processors) but all processors must call this method
   * before it will complete since it includes communication to propagate
   * the part membership of any entities that are on processor boundaries.
   */
  void batch_change_entity_parts( const stk::mesh::EntityVector& entities,
                            const std::vector<PartVector>& add_parts,
                            const std::vector<PartVector>& remove_parts,
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
   *  Destruction of entities in the 'locally_owned_part' schedules
   *  all ghost copies of that entity for destruction during
   *  modification_end.  If the entity is shared with
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
  bool destroy_entity( Entity entity, bool was_ghost = false );

  //------------------------------------

  void generate_new_ids(stk::topology::rank_t rank, size_t numIdsNeeded, std::vector<stk::mesh::EntityId>& requestedIds);

  /** \brief Generate a set of entites with globally unique id's
   *
   *  Each processor fills a request vector asking for a number of new
   *  entities of the given ranks.
   *
   *  ex. request = { 0, 4,  8}
   *  request 0 entites of rank 0, 4 entites of rank 1, and 8 entites
   *  of rank 2
   */
  virtual void generate_new_entities(const std::vector<size_t>& requests,
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

  void declare_relation( Entity e_from ,
      Entity e_to ,
      const RelationIdentifier local_id,
      Permutation permutation,
      OrdinalVector& ordinal_scratch,
      PartVector& part_scratch);

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
  inline bool relation_exist( const Entity entity, EntityRank subcell_rank, RelationIdentifier subcell_id );

  /** \brief  Determine the polarity of the local side,
   *          more efficient if the local_side_id is known.
   */
  inline bool element_side_polarity( const Entity elem ,
      const Entity side , unsigned local_side_id ) const;

  /** \brief  All entities with communication information. */
  const EntityCommListInfoVector & comm_list() const { return m_entity_comm_list; }     // CLEANUP: outside of unit tests only used in percept for consistency checks (copied from BulkData) and printing

  inline VolatileFastSharedCommMapOneRank const& volatile_fast_shared_comm_map(EntityRank rank) const;  // CLEANUP: only used by FieldParallel.cpp

  /** \brief  Query the shared-entity aura.
   *          Is likely to be stale if ownership or sharing has changed
   *          and the 'modification_end' has not been called.
   */
  Ghosting & aura_ghosting() const { return *m_ghosting[AURA] ; }
  Ghosting & shared_ghosting() const { return *m_ghosting[SHARED] ; }

  /** Return the part corresponding to the specified ghosting.
   */
  inline Part& ghosting_part(const Ghosting& ghosting) const;

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
   */
  void change_ghosting( Ghosting & ghosts,
                        const std::vector<EntityProc> & add_send ,
                        const std::vector<EntityKey> & remove_receive = std::vector<EntityKey>());

  // Clear all ghosts for a particular ghosting.
  void destroy_ghosting( Ghosting& ghost_layer );

  /** \brief  Empty every single Ghosting.
   *          Same result, but more efficient than, calling
   *          change_ghosting to remove every single ghosted entity.
   */
  void destroy_all_ghosting();

  /** \brief  Vector of all ghostings */
  const std::vector<Ghosting*> & ghostings() const { return m_ghosting ; }

  /** \brief  Entity Comm functions that are now moved to BulkData
   */
  PairIterEntityComm entity_comm_map(const EntityKey & key) const { return m_entity_comm_map.comm(key); } // CLEANUP: could be replaced by comm_shared_procs outside testing (percept prints all ghostings)
  PairIterEntityComm entity_comm_map_shared(const EntityKey & key) const { return m_entity_comm_map.shared_comm_info(key); }    // CLEANUP: switch apps to use comm_shared_procs
  PairIterEntityComm entity_comm_map(const EntityKey & key, const Ghosting & sub ) const { return m_entity_comm_map.comm(key,sub); }   //CLEANUP: can replace app usage with comm_procs

  int entity_comm_map_owner(const EntityKey & key) const;       // CLEANUP: can make protected

  // Comm-related convenience methods

  bool in_shared(EntityKey key) const { return !entity_comm_map_shared(key).empty(); }         // CLEANUP: only used for testing
  bool in_shared(EntityKey key, int proc) const;         // CLEANUP: only used for testing
  bool in_receive_ghost( EntityKey key ) const;         // CLEANUP: only used for testing
  bool in_receive_ghost( const Ghosting & ghost , EntityKey entity ) const;
  bool in_send_ghost( EntityKey key) const;         // CLEANUP: only used for testing
  bool in_send_ghost( EntityKey key , int proc ) const;         // CLEANUP: only used for testing
  bool is_aura_ghosted_onto_another_proc( EntityKey key ) const;     // CLEANUP: used only by modification_end_for_entity_creation
  bool is_aura_ghosted_onto_proc( EntityKey key, int otherProc ) const;     // CLEANUP: used only by modification_end_for_entity_creation
  bool in_ghost( const Ghosting & ghost , EntityKey key , int proc ) const;     // CLEANUP: can be moved protected
  void comm_procs( EntityKey key, std::vector<int> & procs ) const; //shared and ghosted entities
  void comm_shared_procs( EntityKey key, std::vector<int> & procs ) const; // shared entities
  void shared_procs_intersection( std::vector<EntityKey> & keys, std::vector<int> & procs ) const; // CLEANUP: only used by aero
  void comm_procs( const Ghosting & ghost ,
                   EntityKey key, std::vector<int> & procs ) const;

  inline bool in_index_range(Entity entity) const;
  inline bool is_valid(Entity entity) const;
  size_t count_relations(Entity entity) const;
  bool has_no_relations(Entity entity) const;

  inline const MeshIndex& mesh_index(Entity entity) const;
  inline MeshIndex& mesh_index(Entity entity);                          //CLEANUP: move to protected for Partition unit test
  inline EntityId identifier(Entity entity) const;
  inline EntityRank entity_rank(Entity entity) const;
  inline EntityKey entity_key(Entity entity) const;
  inline EntityState state(Entity entity) const;
  inline void mark_entity(Entity entity, entitySharing sharedType);     //CLEANUP: move this and calling functions to private
  inline entitySharing is_entity_marked(Entity entity) const;           //CLEANUP: move this and calling functions to private
  inline bool add_node_sharing_called() const;                          //CLEANUP: move this and calling functions to private
  inline Bucket & bucket(Entity entity) const;
  inline Bucket * bucket_ptr(Entity entity) const;
  inline Bucket::size_type bucket_ordinal(Entity entity) const;
  inline int parallel_owner_rank(Entity entity) const;
  inline unsigned local_id(Entity entity) const;

  inline void set_local_id(Entity entity, unsigned id);

#ifdef SIERRA_MIGRATION
  typedef int FmwkId; //must be a signed type -- fmwk uses negative values sometimes
  inline FmwkId global_id(stk::mesh::Entity entity) const;
  inline const RelationVector& aux_relations(Entity entity) const;
  inline RelationVector& aux_relations(Entity entity);
  inline void set_global_id(stk::mesh::Entity entity, int id);
  void reserve_relation(stk::mesh::Entity entity, const unsigned num);
  void erase_and_clear_if_empty(stk::mesh::Entity entity, RelationIterator rel_itr);
  void internal_verify_initialization_invariant(stk::mesh::Entity entity);

  inline RelationIterator internal_begin_relation(Entity entity, const Relation::RelationType relation_type) const;
  inline RelationIterator internal_end_relation(Entity entity, const Relation::RelationType relation_type) const;
  inline void compress_relation_capacity(Entity entity);
  bool add_fmwk_data() const { return m_add_fmwk_data; }
#endif

  //
  // Connectivity getter methods. For each entity, you can get connected entities
  // of any rank (i.e. node, edge, face, element, constraint).
  // (Some of those connectivities are empty, for example an element may not
  // have any connected edges, depending whether edges exist in the mesh or not.)
  //
  // For each rank of connected entities, you can get the number of connected entities,
  // pointers to the connected entities themselves,
  // ordinals of those connected entities,
  // and permutations of the connected entities.
  //

  // Nodal connectivities are usually contiguous within a bucket. For example, an element-bucket
  // has all connected nodes for those elements allocated in contiguous memory. This should not be
  // assumed for other kinds of connectivity.

  //These connectivity getter methods are implemented by macros which are located further
  //down in this header. (Search for BEGIN_END_PAIR.)

  inline Entity const* begin(Entity entity, EntityRank rank) const;
  inline Entity const* begin_nodes(Entity entity) const;
  inline Entity const* begin_edges(Entity entity) const;
  inline Entity const* begin_faces(Entity entity) const;
  inline Entity const* begin_elements(Entity entity) const;

  // The ordinal of a connected entity is that entity's local index on the entity it
  // is connected to, as defined by standard exodus conventions. For example, the
  // connected nodes of a hex-8 element will have ordinals in the range 0 .. 7.

  // Connected entities are stored in order of ascending ordinal.
  // For cases like element-node connectivity, where an element always has all of its nodes,
  // the array of connected nodes can be indexed by ordinal.
  // For cases like element-face connectivity, where an element may not have all faces defined,
  // the array of connected faces can not always be indexed by ordinal even
  // though those faces are sorted by ordinal.
  // e.g., a hex-8 element may only have two of its six possible faces.

  inline ConnectivityOrdinal const* begin_ordinals(Entity entity, EntityRank rank) const;
  inline ConnectivityOrdinal const* begin_node_ordinals(Entity entity) const;
  inline ConnectivityOrdinal const* begin_edge_ordinals(Entity entity) const;
  inline ConnectivityOrdinal const* begin_face_ordinals(Entity entity) const;
  inline ConnectivityOrdinal const* begin_element_ordinals(Entity entity) const;

  // The permutation of a connected entity is an integer type which is used
  // to store the polarity and orientation of the entity.

  inline Permutation const* begin_permutations(Entity entity, EntityRank rank) const;
  inline Permutation const* begin_node_permutations(Entity entity) const;
  inline Permutation const* begin_edge_permutations(Entity entity) const;
  inline Permutation const* begin_face_permutations(Entity entity) const;
  inline Permutation const* begin_element_permutations(Entity entity) const;

  unsigned num_connectivity(Entity entity, EntityRank rank) const;

  inline unsigned num_nodes(Entity entity) const;
  inline unsigned num_edges(Entity entity) const;
  inline unsigned num_faces(Entity entity) const;
  inline unsigned num_elements(Entity entity) const;

  unsigned count_valid_connectivity(Entity entity, EntityRank rank) const;
  unsigned count_valid_connectivity(Entity entity) const;

  inline Entity const* end(Entity entity, EntityRank rank) const;
  inline Entity const* end_nodes(Entity entity) const;
  inline Entity const* end_edges(Entity entity) const;
  inline Entity const* end_faces(Entity entity) const;
  inline Entity const* end_elements(Entity entity) const;
  Entity const* end_constraints(Entity entity) const
  { return end(entity, stk::topology::CONSTRAINT_RANK); }
  inline ConnectivityOrdinal const* end_ordinals(Entity entity, EntityRank rank) const;
  inline ConnectivityOrdinal const* end_node_ordinals(Entity entity) const;
  inline ConnectivityOrdinal const* end_edge_ordinals(Entity entity) const;
  inline ConnectivityOrdinal const* end_face_ordinals(Entity entity) const;
  inline ConnectivityOrdinal const* end_element_ordinals(Entity entity) const;
  inline Permutation const* end_permutations(Entity entity, EntityRank rank) const;
  inline Permutation const* end_node_permutations(Entity entity) const;
  inline Permutation const* end_edge_permutations(Entity entity) const;
  inline Permutation const* end_face_permutations(Entity entity) const;
  inline Permutation const* end_element_permutations(Entity entity) const;


  // Return index (offset) of query ordinal if found, num_connectivity otherwise.
  unsigned find_ordinal(Entity entity, EntityRank rank, ConnectivityOrdinal ordinal) const;
  bool has_permutation(Entity entity, EntityRank rank) const;
  bool owned_closure(Entity entity) const { return m_closure_count[entity.local_offset()] > static_cast<uint16_t>(0); }
  void get_selected_nodes(stk::mesh::Selector selector, stk::mesh::EntityVector& nodes);
  size_t total_field_data_footprint(const FieldBase &f, EntityRank rank) const { return m_bucket_repository.total_field_data_footprint(f, rank); }
  size_t total_field_data_footprint(EntityRank rank) const;

  // Print all mesh info
  void dump_all_mesh_info(std::ostream& out = std::cout) const;

  // memoized version
  BucketVector const& get_buckets(EntityRank rank, Selector const& selector) const;

  // non-memoized version.
  void get_buckets(EntityRank rank, Selector const& selector, BucketVector & output_buckets) const;

  //
  //  Get entities of the specified rank that satisify the input selector.
  //  Note entities are returned in bucket order, though no particular order should be relied on
  //
  void get_entities(EntityRank rank, Selector const& selector, EntityVector& output_entities) const;

  bool use_entity_ids_for_resolving_sharing() const { return m_use_identifiers_for_resolving_sharing; }
  void set_use_entity_ids_for_resolving_sharing(bool input) { m_use_identifiers_for_resolving_sharing = input; }


  //Power users only.
  //Call this right after construction, before any field-data has been allocated.
  //If you call this method too late (after any field-data has been allocated, it will have no effect.
  //It turns off field-data updating so that movement of entities between buckets etc., as is done during
  //mesh-setup, will not cause corresponding churn of field-data.
  //Once the mesh is initialized with entities and relations, turn on field-data by calling the
  //method 'allocate_field_data'.
  void deactivate_field_updating();
  bool is_field_updating_active() const { return m_keep_fields_updated; }
  /** If field-data was set to not stay in sync with buckets as the mesh was populated,
   * (by calling 'deactivate_field_updating' right after construction) this call
   * causes field-data to be allocated and field-data updating is re-activated.
   * If field-data was already allocated and staying in sync, then this call is a no-op.
   */
  void allocate_field_data();

protected: //functions
  inline bool internal_set_parallel_owner_rank_but_not_comm_lists(Entity entity, int in_owner_rank);

  impl::EntityRepository &get_entity_repository() { return m_entity_repo; }

  inline void set_state(Entity entity, EntityState entity_state);
  void update_deleted_entities_container();
  std::pair<Entity, bool> internal_create_entity(EntityKey key, size_t preferred_offset = 0);

  /** \brief  Declare a collection of relations by simply iterating
   *          the input and calling declare_relation on each entry.
   */
  void internal_declare_relation( Entity entity, const std::vector<Relation> & rel);

  void internal_declare_relation( Entity entity, const std::vector<Relation> & rel,
                         OrdinalVector& ordinal_scratch, PartVector& part_scratch);

  bool internal_declare_relation(Entity e_from, Entity e_to,
                                 RelationIdentifier local_id,
                                 unsigned sync_count, bool is_back_relation,
                                 Permutation permut);
  /** methods for managing arrays of entity member-data */

  inline void log_created_parallel_copy(Entity entity);

  void internal_change_entity_owner( const std::vector<EntityProc> & arg_change,
                                     bool regenerate_aura = true,
                                     modification_optimization mod_optimization = MOD_END_SORT );

  /*  Entity modification consequences:
   *  1) Change entity relation => update via part relation => change parts
   *  2) Change parts => update forward relations via part relation
   *                  => update via field relation
   */
  void internal_change_entity_parts( Entity ,
                                     const std::vector<Part*> & add_parts ,
                                     const std::vector<Part*> & remove_parts,
                                     bool always_propagate_internal_changes=true);
  virtual bool internal_destroy_entity( Entity entity, bool was_ghost = false );

  void internal_change_ghosting( Ghosting & ghosts,
                                 const std::vector<EntityProc> & add_send ,
                                 const std::vector<EntityKey> & remove_receive,
                                 bool is_full_regen = false);

  //Optional parameter 'always_propagate_internal_changes' is always true except when this function
  //is being called from the sierra-framework. The fmwk redundantly does its own propagation of the
  //internal part changes (mostly induced-part stuff), so it's a performance optimization to avoid
  //the propagation that stk-mesh does.
  void internal_verify_and_change_entity_parts( Entity entity,
                                                const PartVector & add_parts ,
                                                const PartVector & remove_parts,
                                                bool always_propagate_internal_changes=true );

  bool internal_modification_end_for_change_entity_owner( bool regenerate_aura, modification_optimization opt );
  bool internal_modification_end_for_change_parts();

  void mark_entity_and_upward_related_entities_as_modified(Entity entity);

  void resolveUniqueIdForSharedEntityAndCreateCommMapInfoForSharingProcs(std::vector<shared_entity_type> & shared_entity_map);
  void update_shared_entities_global_ids(std::vector<shared_entity_type> & shared_entity_map);
  void resolve_entity_sharing(stk::mesh::EntityRank entityRank, std::vector<EntityKey> &entity_keys);

  void internal_resolve_shared_modify_delete();
  void update_comm_list_based_on_changes_in_comm_map();

  void internal_resolve_ghosted_modify_delete();
  void internal_resolve_shared_membership();
  void internal_update_sharing_comm_map_and_fill_list_modified_shared_entities_of_rank(stk::mesh::EntityRank entityRank, std::vector<stk::mesh::Entity> & shared_new );
  virtual void internal_update_sharing_comm_map_and_fill_list_modified_shared_entities(std::vector<stk::mesh::Entity> & shared_new );
  virtual void internal_resolve_send_ghost_membership();
  void resolve_ownership_of_modified_entities(const std::vector<stk::mesh::Entity> &shared_new);
  void move_entities_to_proper_part_ownership( const std::vector<stk::mesh::Entity> &shared_modified );

  void add_comm_list_entries_for_entities(const std::vector<stk::mesh::Entity>& shared_modified);
  bool entity_comm_map_insert(Entity entity, const EntityCommInfo & val) { return m_entity_comm_map.insert(entity_key(entity), val, parallel_owner_rank(entity)); }
  bool entity_comm_map_erase(  const EntityKey & key, const EntityCommInfo & val) { return m_entity_comm_map.erase(key,val); }
  bool entity_comm_map_erase(  const EntityKey & key, const Ghosting & ghost) { return m_entity_comm_map.erase(key,ghost); }
  void entity_comm_map_clear_ghosting(const EntityKey & key ) { m_entity_comm_map.comm_clear_ghosting(key); }
  void entity_comm_map_clear(const EntityKey & key) { m_entity_comm_map.comm_clear(key); }
  void fixup_ghosted_to_shared_nodes();
  void find_ghosted_nodes_that_need_to_be_shared(stk::mesh::EntityVector& ghosted_nodes_that_are_now_shared);

  /** \brief  Regenerate the shared-entity aura,
   *          adding and removing ghosted entities as necessary.
   *
   *  - a collective parallel operation.
   */
  void internal_regenerate_aura();
  void internal_calculate_sharing(const std::vector<EntityProc> & local_change,
                           const std::vector<EntityProc> & shared_change,
                           const std::vector<EntityProc> & ghosted_change,
                           EntityToDependentProcessorsMap & owned_node_sharing_map);

  void require_ok_to_modify() const ;
  void internal_update_fast_comm_maps();

  void internal_apply_node_sharing(const EntityToDependentProcessorsMap & owned_node_sharing_map);
  void internal_add_node_sharing( Entity node, int sharing_proc );
  void internal_remove_node_sharing( EntityKey node, int sharing_proc );
  void internal_compute_proposed_owned_closure_count(const std::vector<EntityProc> & local_change,
                                                     const std::vector<EntityProc> & shared_change,
                                                     const std::vector<EntityProc> & ghosted_change,
                                                     std::vector<uint16_t> & new_closure_count) const;
  void internal_create_new_owner_map(const std::vector<EntityProc> & local_change,
                                     const std::vector<EntityProc> & shared_change,
                                     const std::vector<EntityProc> & ghosted_change,
                                     NewOwnerMap & new_owner_map);
  void update_dependent_processor_map_from_closure_count(const std::vector<uint16_t> & proposed_closure_count,
                                                         EntityToDependentProcessorsMap & entity_to_dependent_processors_map);
  void internal_print_comm_map( std::string title);
  void internal_communicate_entity_to_dependent_processors_map(
          EntityToDependentProcessorsMap & entity_to_dependent_processors_map);

  impl::BucketRepository& bucket_repository() { return m_bucket_repository; }
  bool internal_modification_end( bool regenerate_aura, modification_optimization opt );

  bool is_entity_in_sharing_comm_map(stk::mesh::Entity entity);
  void erase_sharing_info_using_key(stk::mesh::EntityKey key, stk::mesh::BulkData::GHOSTING_ID ghostingId);
  void add_sharing_info(stk::mesh::Entity entity, stk::mesh::BulkData::GHOSTING_ID ghostingId, int sharingProc);
  void update_sharing_after_change_entity_owner();
  void get_entities_that_have_sharing(std::vector<stk::mesh::Entity> &entitiesThatHaveSharingInfo,
          stk::mesh::EntityToDependentProcessorsMap &entityKeySharing);
  void get_locally_modified_shared_entities(stk::mesh::EntityToDependentProcessorsMap &entityKeySharing, std::vector<std::pair<stk::mesh::EntityKey, int> >& sharedEntities);
  void fill_modified_entities_and_add_sharing_comm_map_info_for_shared_entities(const std::vector<std::pair<stk::mesh::EntityKey, int> > &sharedEntities,
          const std::vector<stk::mesh::Entity>& entitiesThatUsedToHaveSharingInfoBeforeCEO, std::vector<stk::mesh::Entity>& modifiedEntitiesForWhichCommMapsNeedUpdating);
  void erase_all_sharing_for_invalid_entities_on_comm_map();
  void fill_entities_that_have_lost_sharing_info(const std::vector<std::pair<stk::mesh::EntityKey, int> > &sharedEntities,
          const std::vector<stk::mesh::Entity>& entitiesThatUsedToHaveSharingInfoBeforeCEO, std::vector<stk::mesh::Entity>& modifiedEntitiesForWhichCommMapsNeedUpdating);
  void resolve_entity_ownership_and_part_membership_and_comm_list(std::vector<stk::mesh::Entity>& modifiedEntities);

  //reserves space for a new entity, or reclaims space from a previously-deleted entity
  virtual size_t generate_next_local_offset(size_t preferred_offset = 0);

  void entity_setter_debug_check(Entity entity) const
  {
    // The 0-th local_offset is special, it represents the invalid, 0-initialized entity.
    // Client should never try to set properties on this entity even though it's in the index range.
    ThrowAssert(entity.local_offset() > 0);
  }

  void entity_getter_debug_check(Entity entity) const
  {
    ThrowAssertMsg(in_index_range(entity) , "Entity has out-of-bounds offset: " << entity.local_offset() << ", maximum offset is: " << m_entity_states.size() - 1);
  }

  void fix_up_ownership(stk::mesh::Entity entity, int new_owner)
  {
      stk::mesh::EntityKey key = entity_key(entity);
      this->internal_change_owner_in_comm_data(key, new_owner);
      this->internal_set_parallel_owner_rank_but_not_comm_lists(entity, new_owner);
  }

  void require_good_rank_and_id(EntityRank ent_rank, EntityId ent_id) const;
  inline void set_synchronized_count(Entity entity, size_t sync_count);
  void remove_entity_callback(EntityRank rank, unsigned bucket_id, Bucket::size_type bucket_ord);

  bool internal_destroy_relation(Entity e_from ,
                                 Entity e_to,
                                 const RelationIdentifier local_id );
private: //functions


  inline void set_mesh_index(Entity entity, Bucket * in_bucket, Bucket::size_type ordinal );
  inline void set_entity_key(Entity entity, EntityKey key);
  void generate_send_list(const int p_rank, std::vector<EntityProc> & send_list);

  void internal_change_owner_in_comm_data(const EntityKey& key, int new_owner);
  void internal_sync_comm_list_owners();

  void internal_change_entity_key(EntityKey old_key, EntityKey new_key, Entity entity);

  void addMeshEntities(stk::topology::rank_t rank, const std::vector<stk::mesh::EntityId> new_ids,
         const std::vector<Part*> &rem, const std::vector<Part*> &add, std::vector<Entity>& requested_entities);

  // Forbidden
  BulkData();
  BulkData( const BulkData & );
  BulkData & operator = ( const BulkData & );

#ifdef GATHER_GET_BUCKETS_METRICS
  void gather_and_print_get_buckets_metrics() const;
#endif

  void gather_and_print_mesh_partitioning() const;

  // Field callbacks

  void new_bucket_callback(EntityRank rank, const PartVector& superset_parts, size_t capacity, Bucket* new_bucket);
  void new_bucket_caching(EntityRank rank, Bucket* new_bucket);

  //
  //  "fields" is an optional argument, if present copy only the listed fields.
  //
  void copy_entity_fields_callback(EntityRank dst_rank, unsigned dst_bucket_id, Bucket::size_type dst_bucket_ord,
                                   unsigned src_bucket_id, Bucket::size_type src_bucket_ord,
                                   const std::vector<FieldBase*>* fields =NULL);


  void destroy_bucket_callback(EntityRank rank, Bucket const& dying_bucket, unsigned capacity);

  // id_map, indexed by new id, maps to old id
  void reorder_buckets_callback(EntityRank rank, const std::vector<unsigned>& id_map);

  void remove_entity_field_data_callback(EntityRank rank, unsigned bucket_id, Bucket::size_type bucket_ord);
  void add_entity_callback(EntityRank rank, unsigned bucket_id, Bucket::size_type bucket_ord);

  void initialize_arrays();

  void internal_declare_relation(Entity e_from ,
                                 Entity e_to ,
                                 const RelationIdentifier local_id ,
                                 Permutation permut,
                                 OrdinalVector& ordinal_scratch,
                                 PartVector& part_scratch);

  MetaData & meta_data() const { return m_mesh_meta_data ; }

  int determine_new_owner( Entity ) const ;

  Entity internal_declare_entity( EntityRank ent_rank , EntityId ent_id ,
                                   const PartVector & parts );

  void internal_propagate_part_changes( Entity entity, const std::vector<Part*> & removed );

  Ghosting & internal_create_ghosting( const std::string & name );
  void internal_verify_inputs_and_change_ghosting(
                                    Ghosting & ghosts ,
                                    const std::vector<EntityProc> & add_send ,
                                    const std::vector<EntityKey> & remove_receive );

  void ghost_entities_and_fields(Ghosting & ghosting, const std::set<EntityProc , EntityLess>& new_send);

  bool internal_modification_end_for_entity_creation( EntityRank entity_rank, bool regenerate_aura, modification_optimization opt );
  void internal_establish_new_owner(stk::mesh::Entity entity);
  void internal_update_parts_for_shared_entity(stk::mesh::Entity entity, const bool is_entity_shared, const bool did_i_just_become_owner);
  void internal_resolve_shared_modify_delete_second_pass();

  void internal_resolve_parallel_create();

  void internal_basic_part_check(const Part* part,
                                 const unsigned ent_rank,
                                 const unsigned undef_rank,
                                 bool& intersection_ok,
                                 bool& rel_target_ok,
                                 bool& rank_ok) const;

  inline void internal_check_unpopulated_relations(Entity entity, EntityRank rank) const;

  // Returns false if there is a problem. It is expected that
  // verify_change_parts will be called if quick_verify_change_part detects
  // a problem, therefore we leave the generation of an exception to
  // verify_change_parts. We want this function to be as fast as
  // possible.
  bool internal_quick_verify_change_part(const Part* part,
                                         const unsigned ent_rank,
                                         const unsigned undef_rank) const;

  void internal_verify_change_parts( const MetaData   & meta ,
                                     const Entity entity ,
                                     const std::vector<Part*> & parts ) const;

  void require_entity_owner( const Entity entity, int owner) const ;

  void require_metadata_committed();

  bool is_good_rank_and_id(EntityRank ent_rank, EntityId ent_id) const;

  inline bool is_valid_connectivity(Entity entity, EntityRank rank) const;

  void require_valid_relation( const char action[] ,
                               const BulkData & mesh ,
                               const Entity e_from ,
                               const Entity e_to );

  // FIXME: Remove this friend once unit-testing has been refactored
  friend class UnitTestModificationEndWrapper;
  friend class ::stk::mesh::MetaData;
  friend class ::stk::mesh::impl::Partition;
  friend class ::stk::mesh::impl::EntityRepository;
  friend class ::stk::mesh::impl::BucketRepository;
  friend class stk::mesh::Bucket; // for field callback

  void reset_modification_counters();
  std::string create_modification_counts_filename() const;
  void write_modification_counts();
  void write_modification_counts_to_stream_for_method_type(std::ostream& out, enum PublicOrInternalMethod methodType);
  void write_modification_counts_to_stream(std::ostream& out);
  void write_entity_modification_entry(std::ostream& out,
                                         enum PublicOrInternalMethod methodType,
                                         EntityModificationTypes entityModification);
  void write_modification_labels_to_stream_for_method_type(std::ostream& out, enum PublicOrInternalMethod methodType);
  void write_modification_labels_to_stream(std::ostream& out);
  void write_modification_entry_label(std::ostream& out, const std::string& label, enum PublicOrInternalMethod methodType);
  void write_entity_modification_entry_label(std::ostream& out, const std::string& label, enum PublicOrInternalMethod methodType);
  std::string convert_label_for_method_type(const std::string &label, enum PublicOrInternalMethod methodType);

  struct MarkAsModified
  {
      MarkAsModified(BulkData & mesh_in) : mesh(mesh_in) {}
      void operator()(Entity entity);
      BulkData & mesh;
  };

public: // data
  mutable bool m_check_invalid_rels; // TODO REMOVE

protected: //data
  EntityCommDatabase m_entity_comm_map;
  std::vector<Ghosting*> m_ghosting;
  MetaData &m_mesh_meta_data;
  size_t m_sync_count;
  BulkDataSyncState m_sync_state;
  std::vector<int> m_mark_entity;
  bool m_add_node_sharing_called;
  std::vector<uint16_t> m_closure_count;
  std::vector<MeshIndex> m_mesh_indexes;
  impl::EntityRepository m_entity_repo;
  EntityCommListInfoVector m_entity_comm_list;
  std::list<size_t, tracking_allocator<size_t, DeletedEntityTag> > m_deleted_entities_current_modification_cycle;
  GhostReuseMap m_ghost_reuse_map;
  std::vector<EntityKey> m_entity_keys;
  std::vector<uint16_t> m_entity_states;

#ifdef SIERRA_MIGRATION
  bool m_add_fmwk_data; // flag that will add extra data to buckets to support fmwk
  std::vector<FmwkId> m_fmwk_global_ids;
  mutable std::vector<RelationVector* > m_fmwk_aux_relations;   // Relations that can't be managed by STK such as PARENT/CHILD
#endif

private: // data
  Parallel m_parallel;
  VolatileFastSharedCommMap m_volatile_fast_shared_comm_map;
  std::vector<Part*> m_ghost_parts;
  std::list<size_t, tracking_allocator<size_t, DeletedEntityTag> > m_deleted_entities;
  std::string m_modification_begin_description;
  int m_num_fields;
  bool m_keep_fields_updated;
  std::vector<size_t> m_entity_sync_counts;
  std::vector<unsigned> m_local_ids;

  //  ContiguousFieldDataManager m_default_field_data_manager;
  DefaultFieldDataManager m_default_field_data_manager;
  FieldDataManager *m_field_data_manager;
  mutable SelectorBucketMap m_selector_to_buckets_map;
  impl::BucketRepository m_bucket_repository; // needs to be destructed first!
  bool m_use_identifiers_for_resolving_sharing;

#ifdef STK_MESH_MODIFICATION_COUNTERS
  static unsigned m_num_bulk_data_counter;
  unsigned m_modification_counters[NumMethodTypes][NumModificationTypes];
  unsigned m_entity_modification_counters[NumMethodTypes][stk::topology::NUM_RANKS][NumEntityModificationTypes];
#endif

#ifdef GATHER_GET_BUCKETS_METRICS
  mutable SelectorCountMap m_selector_to_count_map;
  mutable size_t m_num_memoized_get_buckets_calls;
  mutable size_t m_num_non_memoized_get_buckets_calls;
  size_t m_num_buckets_inserted_in_cache;
  size_t m_num_buckets_removed_from_cache;
  size_t m_num_modifications;
#endif

};

} // namespace mesh
} // namespace stk

#include "BulkDataInlinedMethods.hpp"

#endif //  stk_mesh_BulkData_hpp
