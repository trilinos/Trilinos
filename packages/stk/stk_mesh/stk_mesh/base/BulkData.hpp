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
#include <stk_mesh/base/Types.hpp>      // for MeshIndex, EntityRank, etc
#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/baseImpl/BucketRepository.hpp>  // for BucketRepository
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <string>                       // for char_traits, string
#include <utility>                      // for pair
#include <vector>                       // for vector
#include <unordered_map>
#include "stk_mesh/base/Bucket.hpp"     // for Bucket, Bucket::size_type, etc
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey, hash_value
#include "stk_mesh/base/FieldDataManager.hpp"
#include "stk_mesh/base/Relation.hpp"   // for Relation, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssert, etc
#include "stk_mesh/base/ModificationSummary.hpp"
#include <stk_mesh/base/ModificationNotifier.hpp>
#include <stk_mesh/base/EntityProcMapping.hpp>
#include "stk_mesh/baseImpl/MeshModification.hpp"
#include "stk_mesh/baseImpl/elementGraph/GraphTypes.hpp"
#include <stk_mesh/baseImpl/elementGraph/MeshDiagnosticObserver.hpp>
#include "stk_mesh/baseImpl/SoloSideIdGenerator.hpp"
#include "stk_mesh/baseImpl/SideSetImpl.hpp"

namespace stk { namespace mesh { class FieldBase; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class EntityProcMapping; } }
namespace stk { namespace mesh { namespace impl { class EntityRepository; } } }
namespace stk { namespace mesh { class FaceCreator; } }
namespace stk { namespace mesh { class ElemElemGraph; } }
namespace stk { namespace mesh { class ElemElemGraphUpdater; } }
namespace stk { class CommSparse; }
namespace stk { namespace mesh { class ModificationObserver; } }
namespace stk { namespace io { class StkMeshIoBroker; } }
namespace stk { namespace mesh { namespace impl { struct RelationEntityToNode; } } }

#include "EntityCommListInfo.hpp"
#include "EntityLess.hpp"
#include "SharedEntityType.hpp"
#include "CommListUpdater.hpp"


namespace stk {
namespace mesh {
class NgpMeshManager;
struct PartStorage;
struct SideSharingData;
enum class FaceCreationBehavior;
}}

namespace stk {
namespace mesh {

void communicate_field_data(const Ghosting & ghosts, const std::vector<const FieldBase *> & fields);
void communicate_field_data(const BulkData & mesh, const std::vector<const FieldBase *> & fields);
void parallel_sum_including_ghosts(const BulkData & mesh, const std::vector<const FieldBase *> & fields);
void skin_mesh( BulkData & mesh, Selector const& element_selector, PartVector const& skin_parts, const Selector * secondary_selector);
void create_edges( BulkData & mesh, const Selector & element_selector, Part * part_to_insert_new_edges );
void internal_create_faces( BulkData & mesh, const Selector & element_selector, bool connect_faces_to_edges, FaceCreationBehavior faceCreationBehavior);
bool process_killed_elements(stk::mesh::BulkData& bulkData,
                             ElemElemGraph& elementGraph,
                             const stk::mesh::EntityVector& killedElements,
                             stk::mesh::Part& active,
                             stk::mesh::impl::ParallelSelectedInfo &remoteActiveSelector,
                             const stk::mesh::PartVector& parts_for_creating_side,
                             const stk::mesh::PartVector& boundary_mesh_parts);
stk::mesh::Entity clone_element_side(stk::mesh::BulkData &bulk,
                                     stk::mesh::EntityId id,
                                     stk::mesh::Entity elem,
                                     stk::mesh::ConnectivityOrdinal ord,
                                     const stk::mesh::PartVector &parts);

namespace impl {
stk::mesh::Entity connect_side_to_element(stk::mesh::BulkData& bulkData, stk::mesh::Entity element,
                                          stk::mesh::EntityId side_global_id, stk::mesh::ConnectivityOrdinal side_ordinal,
                                          stk::mesh::Permutation side_permutation, const stk::mesh::PartVector& parts);
}

typedef std::unordered_map<EntityKey, Entity::entity_value_type, stk::mesh::HashValueForEntityKey> GhostReuseMap;

struct sharing_info
{
    stk::mesh::Entity m_entity;
    int m_sharing_proc;
    int m_owner;
    sharing_info(stk::mesh::Entity entity, int sharing_proc, int owner) :
        m_entity(entity), m_sharing_proc(sharing_proc), m_owner(owner) {}
};

struct SharingInfoLess
{
    bool operator()(const stk::mesh::sharing_info &a, const stk::mesh::sharing_info &b)
    {
        if(a.m_entity == b.m_entity)
            return a.m_owner < b.m_owner;
        else
            return a.m_entity < b.m_entity;
    }
};


class BulkData {

public:
  enum GhostingId { SHARED = 0, AURA = 1 };
  enum EntitySharing : char { NOT_MARKED=0, POSSIBLY_SHARED=1, IS_SHARED=2, NOT_SHARED };
  enum AutomaticAuraOption { NO_AUTO_AURA, AUTO_AURA };
  using ModEndOptimization = impl::MeshModification::modification_optimization;

  /** \brief  Construct mesh bulk data manager conformal to the given
   *          \ref stk::mesh::MetaData "meta data manager" and will
   *          distribute bulk data over the given parallel machine.
   *
   *  - The maximum number of entities per bucket may be supplied.
   *  - The bulk data is in the synchronized or "locked" state.
   */
  BulkData(   MetaData & mesh_meta_data
            , ParallelMachine parallel
            , enum AutomaticAuraOption auto_aura_option = AUTO_AURA
#ifdef SIERRA_MIGRATION
            , bool add_fmwk_data = false
#endif
            , FieldDataManager *field_dataManager = NULL
            , unsigned bucket_capacity = impl::BucketRepository::default_bucket_capacity
            );

  virtual ~BulkData();

  //------------------------------------
  /** \brief  The meta data manager for this bulk data manager. */
  const MetaData & mesh_meta_data() const { return m_mesh_meta_data ; }
        MetaData & mesh_meta_data()       { return m_mesh_meta_data ; }

  /** \brief  Acquire a reference to an NgpMesh that is targeted to build configuration.
   *          Calling this method will automatically update the NgpMesh for you.
   *          Do not store a persistent pointer or reference unless you are willing
   *          to manually call update_ngp_mesh() when appropriate.
   */
  NgpMesh & get_updated_ngp_mesh() const;

  /** \brief  Perform manual update of a persistent NgpMesh instance.
   */
  void update_ngp_mesh() const;

  /** \brief  The parallel machine */
  ParallelMachine parallel() const { return m_parallel.parallel() ; }

  /** \brief  Size of the parallel machine */
  int parallel_size()   const { return m_parallel.parallel_size() ; }

  /** \brief  Rank of the parallel machine's local processor */
  int parallel_rank()   const { return m_parallel.parallel_rank() ; }

  //------------------------------------
  /** \brief  Bulk data has two states:
   *          guaranteed to be parallel synchronized or
   *          modification in progress and may be parallel inconsistent.
   */

  bool in_modifiable_state() const { return m_meshModification.in_modifiable_state(); }
  bool in_synchronized_state() const { return m_meshModification.in_synchronized_state(); }

  bool is_automatic_aura_on() const { return m_autoAuraOption == AUTO_AURA; }


  /** \brief  Count of the number of times that the bulk data has been
   *          parallel synchronized.  This count gets updated with
   *          each call to 'modification_end'.
   */
  size_t synchronized_count() const { return m_meshModification.synchronized_count() ; }

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
  bool modification_begin(const std::string description = std::string("UNSPECIFIED"))
  {
      notifier.notify_modification_begin();
      m_lastModificationDescription = description;
      return m_meshModification.modification_begin(description);
  }

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

  bool modification_end(ModEndOptimization modEndOpt = ModEndOptimization::MOD_END_SORT)
  {
      notifier.notify_started_modification_end();
      return m_meshModification.modification_end(modEndOpt);
  }

  void sort_entities(const stk::mesh::EntitySorterBase& sorter);

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
  virtual void change_entity_owner( const EntityProcVec & arg_change)
  {
      notifier.notify_elements_about_to_move_procs(arg_change);
      m_meshModification.change_entity_owner(arg_change);
      notifier.notify_elements_moved_procs(arg_change);
  }

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
  void update_field_data_states(FieldBase *field);

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
  const_entity_iterator begin_entities(EntityRank ent_rank) const;

  //end-iterator for entities of the specified rank
  const_entity_iterator end_entities(EntityRank ent_rank) const;

  /** \brief  Get entity with a given rank and id */
  Entity get_entity( EntityRank ent_rank , EntityId entity_id ) const;

  /** \brief  Get entity with a given key */
  Entity get_entity( const EntityKey key ) const;

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
  template<typename PARTVECTOR>
  Entity declare_entity( EntityRank ent_rank , EntityId ent_id , const PARTVECTOR& parts); // Mod Mark
  template<typename PARTVECTOR>
  Entity declare_node(EntityId id, const PARTVECTOR& parts);
  template<typename PARTVECTOR>
  Entity declare_edge(EntityId id, const PARTVECTOR& parts);
  template<typename PARTVECTOR>
  Entity declare_element(EntityId id, const PARTVECTOR& parts);
  template<typename PARTVECTOR>
  Entity declare_constraint(EntityId id, const PARTVECTOR& parts);
  Entity declare_entity( EntityRank ent_rank , EntityId ent_id , Part& part); // Mod Mark

  template<typename IDVECTOR>
  void declare_entities(stk::topology::rank_t rank, const IDVECTOR& new_ids,
                        const PartVector &parts, EntityVector &requested_entities);

  /** This overloading of declare_entity that doesn't take a part
   * creates the new entity in the 'universal' part.
   */

  Entity declare_solo_side(const PartVector& parts);

  Entity declare_solo_side(EntityId ent_id, const PartVector& parts);

  template<typename PARTVECTOR>
  Entity declare_element_side(Entity elem, const unsigned side_ordinal, const PARTVECTOR& add_parts = {});

  Entity declare_node(EntityId ent_id);
  Entity declare_element(EntityId ent_id);
  Entity declare_edge(EntityId ent_id);
  Entity declare_constraint(EntityId ent_id);

  EntityId get_solo_side_id();

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
  void change_entity_id( EntityId id, Entity entity);   // Deprecated? // Mod Mark move to StkTransitionBulkData?

  /** \brief  Change the parallel-locally-owned entity's
   *          part membership by adding and/or removing parts
   *
   *  If the locally owned entity is shared or ghosted then
   *  the change will be propagated to the sharing or ghosting
   *  processes by modification_end.
   */
  template<typename PARTVECTOR>
  void change_entity_parts( Entity entity, // Mod Mark
      const PARTVECTOR & add_parts ,
      const PARTVECTOR & remove_parts = PARTVECTOR());

  template<typename PARTVECTOR>
  void change_entity_parts( const EntityVector& entities,
      const PARTVECTOR & add_parts ,
      const PARTVECTOR & remove_parts = PARTVECTOR());

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
  void batch_change_entity_parts( const stk::mesh::EntityVector& entities, // Mod Mark
                            const std::vector<PartVector>& add_parts,
                            const std::vector<PartVector>& remove_parts);

  void batch_change_entity_parts( const stk::mesh::EntityVector& entities, // Mod Mark
                            const PartVector& add_parts,
                            const PartVector& remove_parts);

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
  bool destroy_entity( Entity entity, bool was_ghost = false ); // Mod Mark

  //------------------------------------

  void generate_new_ids_given_reserved_ids(stk::topology::rank_t rank, size_t numIdsNeeded, const std::vector<stk::mesh::EntityId>& reserved_ids, std::vector<stk::mesh::EntityId>& requestedIds) const;
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
  virtual void generate_new_entities(const std::vector<size_t>& requests, // Mod Mark
      std::vector<Entity>& requested_entities);

  Permutation find_permutation( const stk::topology &hr_entity_topo,
                                Entity const *higher_rank_entity_nodes,
                                const stk::topology &side_topo,
                                Entity const *side_nodes,
                                unsigned side_ordinal) const;

  bool check_permutation(Entity entity, Entity rel_entity, unsigned rel_ordinal, Permutation expected) const;

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
  void declare_relation( Entity e_from , // Mod Mark
      Entity e_to ,
      const RelationIdentifier local_id,
      Permutation permutation = stk::mesh::Permutation::INVALID_PERMUTATION);

  void declare_relation( Entity e_from ,
                         const EntityVector& to_entities);
 
  //it's ugly to have 3 scratch-space vectors in the API, but for now
  //it is a big performance improvement. TODO: improve the internals to remove
  //the need for these vectors.
  void declare_relation( Entity e_from ,
                         Entity e_to ,
                         const RelationIdentifier local_id ,
                         Permutation permut,
                         OrdinalVector& scratch1,
                         OrdinalVector& scratch2,
                         OrdinalVector& scratch3);

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
  bool destroy_relation( Entity e_from , // Mod Mark
                         Entity e_to,
                         const RelationIdentifier local_id );


  // Check if entity has a specific relation to an entity of subcell_rank
  inline bool relation_exist( const Entity entity, EntityRank subcell_rank, RelationIdentifier subcell_id );

  /** \brief  Determine the polarity of the local side,
   *          more efficient if the local_side_id is known.
   */
  inline bool element_side_polarity( const Entity elem ,
      const Entity side , unsigned local_side_id ) const;

  inline VolatileFastSharedCommMapOneRank const& volatile_fast_shared_comm_map(EntityRank rank) const;  // CLEANUP: only used by FieldParallel.cpp
  inline const std::vector<int>& all_sharing_procs(stk::mesh::EntityRank rank) const { return m_all_sharing_procs[rank]; }

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
  Ghosting & create_ghosting( const std::string & name ); // Mod Mark

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
                        const std::vector<EntityKey> & remove_receive = std::vector<EntityKey>()); // Mod Mark

  void batch_add_to_ghosting(Ghosting &ghosting, const EntityProcVec &entitiesAndDestinationProcs); // Mod Mark

  // Clear all ghosts for a particular ghosting.
  void destroy_ghosting( Ghosting& ghost_layer ); // Mod Mark

  /** \brief  Empty every single Ghosting.
   *          Same result, but more efficient than, calling
   *          change_ghosting to remove every single ghosted entity.
   */
  void destroy_all_ghosting(); // Mod Mark

  /** \brief  Vector of all ghostings */
  const std::vector<Ghosting*> & ghostings() const { return m_ghosting ; }

  size_t get_num_communicated_entities() const { return m_entity_comm_list.size(); }

  bool in_shared(EntityKey key) const { return !internal_entity_comm_map_shared(key).empty(); }
  bool in_shared(Entity entity) const;
  bool in_shared(EntityKey key, int proc) const;
  bool in_shared(Entity entity, int proc) const;
  bool in_receive_ghost( EntityKey key ) const;
  bool in_receive_ghost( Entity entity ) const;
  bool in_receive_ghost( const Ghosting & ghost , EntityKey key ) const;
  bool in_receive_ghost( const Ghosting & ghost , Entity entity ) const;
  bool in_receive_custom_ghost( EntityKey key ) const;
  bool in_send_ghost( EntityKey key) const;
  bool in_send_ghost( Entity entity) const;
  bool in_send_ghost( EntityKey key , int proc ) const;
  bool in_send_ghost( const Ghosting & ghosting, EntityKey key, int proc) const;
  bool is_aura_ghosted_onto_another_proc( EntityKey key ) const;
  bool in_ghost( const Ghosting & ghost , EntityKey key , int proc ) const;
  bool in_ghost( const Ghosting & ghost , Entity entity , int proc ) const;
  bool in_ghost( const Ghosting & ghost , Entity entity ) const;
  void shared_procs_intersection(const std::vector<EntityKey> & keys, std::vector<int> & procs ) const;
  void shared_procs_intersection(const EntityVector & entities, std::vector<int> & procs ) const;

  // Comm-related convenience methods

  bool is_communicated_with_proc(Entity entity, int proc) const;
  void comm_procs( EntityKey key, std::vector<int> & procs ) const; //shared and ghosted entities
  void comm_procs( const Ghosting & ghost , EntityKey key, std::vector<int> & procs ) const;
  void comm_shared_procs( EntityKey key, std::vector<int> & procs ) const;
  void comm_shared_procs( Entity entity, std::vector<int> & procs ) const;

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
  inline Bucket & bucket(Entity entity) const;
  inline Bucket * bucket_ptr(Entity entity) const;
  inline Bucket::size_type bucket_ordinal(Entity entity) const;
  inline int parallel_owner_rank(Entity entity) const;
  inline unsigned local_id(Entity entity) const;

  inline void set_local_id(Entity entity, unsigned id);

#ifdef SIERRA_MIGRATION
  typedef int64_t FmwkId; //must be a signed type -- fmwk uses negative values sometimes
  inline FmwkId global_id(stk::mesh::Entity entity) const;
  inline const RelationVector& aux_relations(Entity entity) const;
  inline RelationVector& aux_relations(Entity entity); // Mod Mark
  inline void set_global_id(stk::mesh::Entity entity, FmwkId id);
  void reserve_relation(stk::mesh::Entity entity, const unsigned num); // Mod Mark
  void erase_and_clear_if_empty(stk::mesh::Entity entity, RelationIterator rel_itr); // Mod Mark
  void internal_verify_initialization_invariant(stk::mesh::Entity entity);

  inline RelationIterator internal_begin_relation(Entity entity, const Relation::RelationType relation_type) const;
  inline RelationIterator internal_end_relation(Entity entity, const Relation::RelationType relation_type) const;
  inline void compress_relation_capacity(Entity entity);
  bool add_fmwk_data() const { return m_add_fmwk_data; }
#endif

  //
  //  Get the maximum ID that is allowed by the host code e.g. (64 vs. 32 bit limits)
  //
  uint64_t get_max_allowed_id() const;

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
  unsigned num_sides(Entity entity) const;

  unsigned count_valid_connectivity(Entity entity, EntityRank rank) const;

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
  bool owned_closure(Entity entity) const
  {
      return m_closure_count[entity.local_offset()] > static_cast<uint16_t>(0);
  }

  size_t total_field_data_footprint(const FieldBase &f, EntityRank rank) const { return m_bucket_repository.total_field_data_footprint(f, rank); }
  size_t total_field_data_footprint(EntityRank rank) const;

  // Print all mesh info, consider using:
  // std::ostringstream oss;
  // oss << "output." << parallel_rank();
  // std::ofstream out(oss.str(), std::ios_base::app);
  // dump_all_mesh_info(out);
  // out.close();
  void dump_all_mesh_info(std::ostream& out) const;
  void dump_mesh_per_proc(const std::string& fileNamePrefix) const;
  void dump_mesh_bucket_info(std::ostream& out, Bucket* bucket) const;

  // memoized version
  BucketVector const& get_buckets(EntityRank rank, Selector const& selector) const;

  // non-memoized version.
  void get_buckets(EntityRank rank, Selector const& selector, BucketVector & output_buckets) const;

  //
  //  Get entities of the specified rank that satisfy the input selector.
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

  void reallocate_field_data(stk::mesh::FieldBase & field);

  const std::string & get_last_modification_description() const { return m_lastModificationDescription; }

  void register_observer(std::shared_ptr<stk::mesh::ModificationObserver> observer) const;
  void unregister_observer(std::shared_ptr<ModificationObserver> observer) const;
  template<typename ObserverType>
  bool has_observer_type() const { return notifier.has_observer_type<ObserverType>(); }
  template<typename ObserverType>
  std::vector<std::shared_ptr<ObserverType>> get_observer_type() const { return notifier.get_observer_type<ObserverType>(); }

  void initialize_face_adjacent_element_graph();
  void delete_face_adjacent_element_graph();
  stk::mesh::ElemElemGraph& get_face_adjacent_element_graph();
  const stk::mesh::ElemElemGraph& get_face_adjacent_element_graph() const;
  bool has_face_adjacent_element_graph() const;

  void enable_mesh_diagnostic_rule(stk::mesh::MeshDiagnosticFlag flag);
  unsigned get_mesh_diagnostic_error_count() const ;
  void throw_on_mesh_diagnostic_error();

  size_t get_size_of_entity_index_space() const { return m_entity_keys.size(); }

  void destroy_elements_of_topology(stk::topology topologyToDelete);

  bool does_sideset_exist(const stk::mesh::Part &part) const;
  SideSet& create_sideset(const stk::mesh::Part &part, bool fromInput = false);
  const SideSet& get_sideset(const stk::mesh::Part &part) const;
  SideSet& get_sideset(const stk::mesh::Part &part);
  size_t get_number_of_sidesets() const;
  bool was_mesh_modified_since_sideset_creation();
  void clear_sidesets();
  void clear_sideset(const stk::mesh::Part &part);
  std::vector<SideSet *> get_sidesets();
  std::vector<const SideSet *> get_sidesets() const;
  void synchronize_sideset_sync_count();

  void clone_solo_side_id_generator(const stk::mesh::BulkData &oldBulk);
  void create_side_entities(const SideSet &sideSet, const stk::mesh::PartVector& parts);

  bool supports_large_ids() const {return m_supportsLargeIds; }
  void set_large_ids_flag(bool largeIds) { m_supportsLargeIds = largeIds; }

protected: //functions
  Entity declare_entity( EntityRank ent_rank , EntityId ent_id);// Mod Mark


  template<typename PARTVECTOR>
  Entity declare_element_side_with_id(const stk::mesh::EntityId globalSideId,
                              Entity elem,
                              const unsigned sideOrd,
                              const PARTVECTOR& parts);

  bool resolve_node_sharing()
  {
      return m_meshModification.resolve_node_sharing();
  }

  bool modification_end_after_node_sharing_resolution()
  {
      notifier.notify_started_modification_end();
      return m_meshModification.modification_end_after_node_sharing_resolution();
  }

  void make_mesh_parallel_consistent_after_element_death(const std::vector<sharing_info>& shared_modified,
                                                         const stk::mesh::EntityVector& deletedSides,
                                                         stk::mesh::ElemElemGraph &elementGraph,
                                                         const stk::mesh::EntityVector &killedElements,
                                                         stk::mesh::Part &activePart);

  void make_mesh_parallel_consistent_after_skinning(const std::vector<sharing_info>& sharedModified);

  bool modification_end_for_entity_creation( const std::vector<EntityRank> & entity_rank_vector,
                                             stk::mesh::impl::MeshModification::modification_optimization opt = stk::mesh::impl::MeshModification::MOD_END_SORT); // Mod Mark

  bool internal_modification_end_for_skin_mesh( EntityRank entity_rank,
                                                stk::mesh::impl::MeshModification::modification_optimization opt,
                                                stk::mesh::Selector selectedToSkin,
                                                const stk::mesh::Selector * only_consider_second_element_from_this_selector); // Mod Mark

  bool inputs_ok_and_need_ghosting(Ghosting & ghosts ,
                               const std::vector<EntityProc> & add_send ,
                               const std::vector<Entity> & remove_receive,
                               std::vector<EntityProc> &filtered_add_send);

  void internal_batch_add_to_ghosting(Ghosting &ghosting, const EntityProcVec &entitiesAndDestinationProcs); // Mod Mark

  void ghost_entities_and_fields(Ghosting & ghosting,
                                 const std::set<EntityProc , EntityLess>& new_send,
                                 bool isFullRegen = false);

  void conditionally_add_entity_to_ghosting_set(const stk::mesh::Ghosting &ghosting,
                                                stk::mesh::Entity entity,
                                                int toProc,
                                                std::set <stk::mesh::EntityProc,
                                                stk::mesh::EntityLess > &entitiesWithClosure);

  void add_closure_entities(const stk::mesh::Ghosting& ghosting,
                            const stk::mesh::EntityProcVec& entities,
                            std::set <stk::mesh::EntityProc,
                            stk::mesh::EntityLess > &entitiesWithClosure);

  const EntityCommListInfoVector & internal_comm_list() const { return m_entity_comm_list; }
  PairIterEntityComm internal_entity_comm_map(const EntityKey & key) const { return m_entity_comm_map.comm(key); }

  PairIterEntityComm internal_entity_comm_map(Entity entity) const
  {
    if (m_entitycomm[entity.local_offset()] != nullptr) {
      const EntityCommInfoVector& vec = m_entitycomm[entity.local_offset()]->comm_map;
      return PairIterEntityComm(vec.begin(), vec.end());
    }
    return PairIterEntityComm();
  }

  PairIterEntityComm internal_entity_comm_map(const EntityKey & key, const Ghosting & sub ) const { return m_entity_comm_map.comm(key,sub); }

  PairIterEntityComm internal_entity_comm_map(Entity entity, const Ghosting & sub ) const
  {
    if (m_entitycomm[entity.local_offset()] != nullptr) {
      const EntityCommInfoVector& vec = m_entitycomm[entity.local_offset()]->comm_map;
      const EntityCommInfo s_begin( sub.ordinal() ,     0 );
      const EntityCommInfo s_end(   sub.ordinal() + 1 , 0 );
    
      EntityCommInfoVector::const_iterator i = vec.begin();
      EntityCommInfoVector::const_iterator e = vec.end();
    
      i = std::lower_bound( i , e , s_begin );
      e = std::lower_bound( i , e , s_end );
    
      return PairIterEntityComm( i , e );
    }
    return PairIterEntityComm();
  }

  inline EntitySharing internal_is_entity_marked(Entity entity) const;
  PairIterEntityComm internal_entity_comm_map_shared(const EntityKey & key) const { return m_entity_comm_map.shared_comm_info(key); }

  PairIterEntityComm internal_entity_comm_map_shared(Entity entity) const
  {
    const EntityComm* entityComm = m_entitycomm[entity.local_offset()];
    if (entityComm != nullptr) {
      if (entityComm->isShared) {
        const EntityCommInfoVector& vec = entityComm->comm_map;
        return shared_comm_info_range(vec);
      }
    }
    return PairIterEntityComm();
  }

  void markEntitiesForResolvingSharingInfoUsingNodes(stk::mesh::EntityRank entityRank, bool onlyConsiderSoloSides, std::vector<shared_entity_type>& shared_entities);
  virtual void sortNodesIfNeeded(std::vector<stk::mesh::EntityKey>& nodes);

  void gather_shared_nodes(std::vector<Entity> & shared_nodes);

  inline void set_state(Entity entity, EntityState entity_state);
  inline void set_entity_key(Entity entity, EntityKey key);
  void update_deleted_entities_container();
  std::pair<Entity, bool> internal_create_entity(EntityKey key, size_t preferred_offset = 0); // Mod Mark
  std::pair<Entity, bool> internal_get_or_create_entity_with_notification(EntityKey key, size_t preferred_offset = 0);


  /** \brief  Declare a collection of relations by simply iterating
   *          the input and calling declare_relation on each entry.
   */
  void internal_declare_relation( Entity entity, const std::vector<Relation> & rel); // Mod Mark

  void internal_declare_relation( Entity entity, const std::vector<Relation> & rel,
                         OrdinalVector& scratch1); // Mod Mark

  bool internal_declare_relation(Entity e_from, Entity e_to,
                                 RelationIdentifier local_id, Permutation permut); // Mod Mark

  /** methods for managing arrays of entity member-data */

  inline void log_created_parallel_copy(Entity entity);

  void internal_change_entity_owner( const std::vector<EntityProc> & arg_change,
                                     stk::mesh::impl::MeshModification::modification_optimization mod_optimization = stk::mesh::impl::MeshModification::MOD_END_SORT );  // Mod Mark

  void internal_change_entity_parts_without_propogating_to_downward_connected_entities(Entity entity, const OrdinalVector& add_parts, const OrdinalVector& remove_parts, OrdinalVector& parts_removed, OrdinalVector& newBucketPartList, OrdinalVector& scratchSpace);
  void internal_determine_inducible_parts(Entity entity, const OrdinalVector& add_parts, const OrdinalVector& parts_removed, OrdinalVector& inducible_parts_added, OrdinalVector& inducible_parts_removed);
  void internal_determine_inducible_parts_and_propagate_to_downward_connected_entities(Entity entity, const OrdinalVector& add_parts, const OrdinalVector& parts_removed);

  /*  Entity modification consequences:
   *  1) Change entity relation => update via part relation => change parts
   *  2) Change parts => update forward relations via part relation
   *                  => update via field relation
   */
  void internal_change_entity_parts( Entity ,
                                     const OrdinalVector& add_parts ,
                                     const OrdinalVector& remove_parts,
                                     OrdinalVector& scratchOrdinalVec, OrdinalVector& scratchSpace); // Mod Mark

  bool internal_destroy_entity_with_notification(Entity entity, bool wasGhost = false); // Mod Mark
  virtual bool internal_destroy_entity(Entity entity, bool wasGhost = false);

  void internal_change_ghosting( Ghosting & ghosts,
                                 const std::vector<EntityProc> & add_send ,
                                 const std::vector<Entity> & remove_receive,
                                 bool is_full_regen = false); // Mod Mark

  void internal_change_ghosting( Ghosting & ghosts,
                                 EntityProcMapping& entityProcMapping);

  void internal_add_to_ghosting( Ghosting &ghosting, const std::vector<EntityProc> &add_send); // Mod Mark

  //Optional parameter 'always_propagate_internal_changes' is always true except when this function
  //is being called from the sierra-framework. The fmwk redundantly does its own propagation of the
  //internal part changes (mostly induced-part stuff), so it's a performance optimization to avoid
  //the propagation that stk-mesh does.
  template<typename PARTVECTOR>
  void internal_verify_and_change_entity_parts( Entity entity,
                                                const PARTVECTOR & add_parts ,
                                                const PARTVECTOR & remove_parts); // Mod Mark

  template<typename PARTVECTOR>
  void internal_verify_and_change_entity_parts( const EntityVector& entities,
                                                const PARTVECTOR & add_parts ,
                                                const PARTVECTOR & remove_parts);

  void internal_insert_all_parts_induced_from_higher_rank_entities_to_vector(stk::mesh::Entity entity,
                                                                               stk::mesh::Entity e_to,
                                                                               OrdinalVector &to_add);

  bool internal_modification_end_for_change_entity_owner( stk::mesh::impl::MeshModification::modification_optimization opt ); // Mod Mark
  bool internal_modification_end_for_change_parts(); // Mod Mark
  void internal_modification_end_for_change_ghosting(); // Mod Mark

  void mark_entity_and_upward_related_entities_as_modified(Entity entity);

  void set_common_entity_key_and_fix_ordering_of_nodes_and_update_comm_map(std::vector<shared_entity_type> & shared_entity_map);
  void find_and_delete_internal_faces(stk::mesh::EntityRank entityRank,
                                      const stk::mesh::Selector *only_consider_second_element_from_this_selector); // Mod Mark

  void filter_upward_ghost_relations(const Entity entity, std::function<void(Entity)> filter);
  EntityVector get_upward_send_ghost_relations(const Entity entity);
  EntityVector get_upward_recv_ghost_relations(const Entity entity);
  void add_entity_to_same_ghosting(Entity entity, Entity connectedGhost);
  void update_comm_list_based_on_changes_in_comm_map();

  void internal_resolve_formerly_shared_entities(const stk::mesh::EntityVector& entitiesNoLongerShared);
  void internal_resolve_ghosted_modify_delete(const stk::mesh::EntityVector& entitiesNoLongerShared);
  void internal_resolve_shared_part_membership_for_element_death(); // Mod Mark

  void remove_unneeded_induced_parts(stk::mesh::Entity entity, const EntityCommInfoVector& entity_comm_info,
          PartStorage& part_storage, stk::CommSparse& comm);

  void internal_resolve_shared_membership(const stk::mesh::EntityVector & entitiesNoLongerShared); // Mod Mark
  virtual void internal_resolve_parallel_create();

  virtual void internal_resolve_parallel_create(const std::vector<stk::mesh::EntityRank>& ranks); // Mod Mark
  void internal_resolve_parallel_create_nodes();
  void internal_resolve_parallel_create_edges_and_faces();

  void internal_update_sharing_comm_map_and_fill_list_modified_shared_entities_of_rank(stk::mesh::EntityRank entityRank, std::vector<stk::mesh::Entity> & shared_new ); // Mod Mark
  void internal_send_part_memberships_from_owner(const std::vector<EntityProc> &send_list);

  virtual void internal_update_sharing_comm_map_and_fill_list_modified_shared_entities(std::vector<stk::mesh::Entity> & shared_new );
  void extract_entity_from_shared_entity_type(const std::vector<shared_entity_type>& shared_entities, std::vector<Entity>& shared_new);
  virtual void fill_shared_entities_of_rank_while_updating_sharing_info(stk::mesh::EntityRank rank, std::vector<Entity> &shared_new);
  virtual void internal_update_sharing_comm_map_and_fill_list_modified_shared_entities_of_node_rank(stk::mesh::EntityVector& shared_new);

  virtual void internal_resolve_send_ghost_membership();
  virtual bool should_sort_buckets_by_first_entity_identifier() const { return false; }
  void resolve_ownership_of_modified_entities(const std::vector<stk::mesh::Entity> &shared_new); // Mod Mark
  void move_entities_to_proper_part_ownership( const std::vector<stk::mesh::Entity> &shared_modified ); // Mod Mark

  void add_comm_list_entries_for_entities(const std::vector<stk::mesh::Entity>& shared_modified);

  std::pair<EntityComm*,bool> entity_comm_map_insert(Entity entity, const EntityCommInfo &val)
  {
      EntityKey key = entity_key(entity);
      std::pair<EntityComm*,bool> result = m_entity_comm_map.insert(key, val, parallel_owner_rank(entity));
      if (val.ghost_id == 0) {
        result.first->isShared = true;
      }
      else {
        result.first->isGhost = true;
      }
      if(result.second)
      {
          m_entitycomm[entity.local_offset()] = result.first;
          m_modSummary.track_comm_map_insert(entity, val);
          notifier.notify_local_entity_comm_info_changed(key.rank());
      }
      return result;
  }
  void remove_entity_comm(Entity entity)
  {
    m_entitycomm[entity.local_offset()] = nullptr;
  }

  bool entity_comm_map_erase(const EntityKey &key, const EntityCommInfo &val)
  {
      bool didErase = m_entity_comm_map.erase(key, val);
      if(didErase)
      {
          m_modSummary.track_comm_map_erase(key, val);
          notifier.notify_local_entity_comm_info_changed(key.rank());
      }
      return didErase;
  }
  bool entity_comm_map_erase(const EntityKey &key, const Ghosting &ghost)
  {
      bool didErase = m_entity_comm_map.erase(key, ghost);
      if(didErase)
      {
          m_modSummary.track_comm_map_erase(key, ghost);
          notifier.notify_local_entity_comm_info_changed(key.rank());
      }
      return didErase;
  }
  void entity_comm_map_clear_ghosting(const EntityKey & key)
  {
      bool didClear = m_entity_comm_map.comm_clear_ghosting(key);
      if (didClear) {
          m_modSummary.track_comm_map_clear_ghosting(key);
          notifier.notify_local_entity_comm_info_changed(key.rank());
      }
  }
  void entity_comm_map_clear(const EntityKey & key)
  {
      bool didClear = m_entity_comm_map.comm_clear(key);
      if (didClear){
          m_modSummary.track_comm_map_clear(key);
          notifier.notify_local_entity_comm_info_changed(key.rank());
      }
  }

  /** \brief  Regenerate the shared-entity aura,
   *          adding and removing ghosted entities as necessary.
   *
   *  - a collective parallel operation.
   */
  void internal_regenerate_aura();
  void fill_list_of_entities_to_send_for_aura_ghosting(EntityProcMapping& send);

  void require_ok_to_modify() const ;
  void internal_update_fast_comm_maps();

  impl::BucketRepository& bucket_repository() { return m_bucket_repository; }

  bool is_entity_in_sharing_comm_map(stk::mesh::Entity entity);
  void erase_sharing_info_using_key(stk::mesh::EntityKey key, stk::mesh::BulkData::GhostingId ghostingId);
  void add_sharing_info(stk::mesh::Entity entity, stk::mesh::BulkData::GhostingId ghostingId, int sharingProc);
  void update_sharing_after_change_entity_owner(); // Mod Mark
  void get_entities_that_have_sharing(std::vector<stk::mesh::Entity> &entitiesThatHaveSharingInfo,
          stk::mesh::EntityToDependentProcessorsMap &entityKeySharing);
  void get_locally_modified_shared_entities(stk::mesh::EntityToDependentProcessorsMap &entityKeySharing, std::vector<std::pair<stk::mesh::EntityKey, int> >& sharedEntities);
  void fill_modified_entities_and_add_sharing_comm_map_info_for_shared_entities(const std::vector<std::pair<stk::mesh::EntityKey, int> > &sharedEntities,
          const std::vector<stk::mesh::Entity>& entitiesThatUsedToHaveSharingInfoBeforeCEO, std::vector<stk::mesh::Entity>& modifiedEntitiesForWhichCommMapsNeedUpdating);
  void erase_all_sharing_for_invalid_entities_on_comm_map();
  void fill_entities_that_have_lost_sharing_info(const std::vector<std::pair<stk::mesh::EntityKey, int> > &sharedEntities,
          const std::vector<stk::mesh::Entity>& entitiesThatUsedToHaveSharingInfoBeforeCEO, std::vector<stk::mesh::Entity>& modifiedEntitiesForWhichCommMapsNeedUpdating);
  void resolve_entity_ownership_and_part_membership_and_comm_list(std::vector<stk::mesh::Entity>& modifiedEntities); // Mod Mark

  //reserves space for a new entity, or reclaims space from a previously-deleted entity
  virtual Entity generate_new_entity(unsigned preferred_offset = 0);

  void entity_setter_debug_check(Entity entity) const
  {
    // The 0-th local_offset is special, it represents the invalid, 0-initialized entity.
    // Client should never try to set properties on this entity even though it's in the index range.
    ThrowAssert(entity.local_offset() > 0);
  }

  void entity_getter_debug_check(Entity entity) const
  {
    ThrowAssertMsg(in_index_range(entity) , "Entity has out-of-bounds offset: " << entity.local_offset() << ", maximum offset is: " << m_entity_keys.size() - 1);
  }

  void require_good_rank_and_id(EntityRank ent_rank, EntityId ent_id) const;
  void remove_entity_callback(EntityRank rank, unsigned bucket_id, Bucket::size_type bucket_ord);

  bool internal_destroy_relation(Entity e_from ,
                                 Entity e_to,
                                 const RelationIdentifier local_id ); // Mod Mark

  void check_mesh_consistency();
  bool comm_mesh_verify_parallel_consistency(std::ostream & error_log);
  void delete_shared_entities_which_are_no_longer_in_owned_closure(EntityProcVec& entitiesToRemoveFromSharing); // Mod Mark
  virtual void remove_entities_from_sharing(const EntityProcVec& entitiesToRemoveFromSharing, stk::mesh::EntityVector & entitiesNoLongerShared);
  virtual void check_if_entity_from_other_proc_exists_on_this_proc_and_update_info_if_shared(std::vector<shared_entity_type>& shared_entity_map, int proc_id, const shared_entity_type &sentity);
  void update_owner_global_key_and_sharing_proc(stk::mesh::EntityKey global_key_other_proc,  shared_entity_type& shared_entity_this_proc, int proc_id) const;
  void update_shared_entity_this_proc(EntityKey global_key_other_proc, shared_entity_type& shared_entity_this_proc, int proc_id);
  void mark_shared_sides_and_fill_list_of_sides_not_on_boundary(std::vector<shared_entity_type>& shared_entity_map,
          int proc_id, shared_entity_type &sentity, std::vector<stk::mesh::EntityKeyProc> &entities_to_send_data,
          const stk::mesh::Selector *only_consider_second_element_from_this_selector);

  void fillVectorOfSharedEntitiesByProcessor(std::vector<shared_entity_type> & potentially_shared_sides, std::vector<std::vector<shared_entity_type> > &shared_entities_by_proc);

  void unpack_shared_entities(stk::CommSparse &comm, std::vector< std::pair<int, shared_entity_type> > &shared_entities_and_proc);

  void unpackEntityFromOtherProcAndUpdateInfoIfSharedLocally(stk::CommSparse &comm, std::vector<shared_entity_type> & shared_entity_map);

  virtual void change_entity_key_and_update_sharing_info(std::vector<shared_entity_type> & shared_entity_map);
  void change_entity_key_to_match_owner(const std::vector<shared_entity_type> & potentially_shared_sides);
  void insert_sharing_info_into_comm_map(const std::vector<shared_entity_type> & potentially_shared_sides);

  inline void internal_mark_entity(Entity entity, EntitySharing sharedType);

  void internal_change_entity_key(EntityKey old_key, EntityKey new_key, Entity entity); // Mod Mark

  void resolve_incremental_ghosting_for_entity_creation_or_skin_mesh(EntityRank entity_rank, stk::mesh::Selector selectedToSkin, bool connectFacesToPreexistingGhosts);

  void internal_finish_modification_end(impl::MeshModification::modification_optimization opt); // Mod Mark

  void internal_set_owner(Entity entity, int new_owner);

  std::vector<uint64_t> internal_get_ids_in_use(stk::topology::rank_t rank, const std::vector<stk::mesh::EntityId>& reserved_ids = std::vector<stk::mesh::EntityId>()) const;

  virtual void de_induce_parts_from_nodes(const stk::mesh::EntityVector & deactivatedElements, stk::mesh::Part & activePart);

  virtual void remove_boundary_faces_from_part(stk::mesh::ElemElemGraph &graph,
                                         const stk::mesh::EntityVector & deactivatedElements,
                                         const stk::mesh::Part & activePart);

  virtual void internal_adjust_entity_and_downward_connectivity_closure_count(stk::mesh::Entity entity,
                                                                      stk::mesh::Bucket *bucket_old,
                                                                      int closureCountAdjustment); // Mod Mark

  inline void set_mesh_index(Entity entity, Bucket * in_bucket, Bucket::size_type ordinal );

  stk::mesh::impl::BucketRepository& get_bucket_repository() { return m_bucket_repository; }

  void set_modification_summary_proc_id(int proc_id) { m_modSummary.set_proc_id(proc_id); }
  virtual void notify_finished_mod_end();

  void use_elem_elem_graph_to_determine_shared_entities(std::vector<stk::mesh::Entity>& shared_entities);
  void use_nodes_to_resolve_sharing(stk::mesh::EntityRank rank, std::vector<Entity>& shared_new, bool onlyConsiderSoloSides = false);
  void change_connectivity_for_edge_or_face(stk::mesh::Entity side, const std::vector<stk::mesh::EntityKey>& node_keys);
  void update_side_elem_permutations(Entity side);
  void resolve_parallel_side_connections(std::vector<SideSharingData>& sideSharingDataToSend,
                                         std::vector<SideSharingData>& sideSharingDataReceived);
  void add_comm_map_for_sharing(const std::vector<SideSharingData>& sidesSharingData, stk::mesh::EntityVector& shared_entities);

private:
  void register_device_mesh() const;
  void unregister_device_mesh() const;

  void create_ngp_mesh_manager() const;
  void record_entity_deletion(Entity entity);
  void break_boundary_relations_and_delete_buckets(const std::vector<impl::RelationEntityToNode> & relationsToDestroy, const stk::mesh::BucketVector & bucketsToDelete);
  void delete_buckets(const stk::mesh::BucketVector & buckets);
  void mark_entities_as_deleted(stk::mesh::Bucket * bucket);

  void generate_ghosting_receive_list(const stk::mesh::Ghosting &ghosting, const std::vector <Entity> &remove_receive,
    std::vector<Entity> &recvGhosts, std::vector<bool>& ghostStatus);

  void verify_and_filter_add_send(Ghosting & ghosts, const std::vector<EntityProc> & add_send, bool &need_to_change_ghosting,
                                  bool &add_send_is_owned, std::vector <EntityProc> &filtered_add_send );

  void verify_remove_receive(Ghosting & ghosts, const std::vector<Entity> & remove_receive, bool &need_to_change_ghosting, bool &remove_receive_are_part_of_this_ghosting);

  bool check_errors_and_determine_if_ghosting_needed_in_parallel(const stk::mesh::Ghosting &ghosts,
                                        bool add_send_is_owned,
                                        bool remove_receive_are_part_of_this_ghosting,
                                        bool need_to_change_ghosting,
                                        const std::vector<EntityProc> & add_send,
                                        const std::vector<Entity> & remove_receive);

  void delete_unneeded_entries_from_the_comm_list();

  void internal_resolve_sharing_and_ghosting_for_sides(bool connectFacesToPreexistingGhosts);

#ifdef __CUDACC__
public:
#endif
  struct EntityParallelState {
    int                 from_proc;
    EntityState         state;
    EntityCommListInfo  comm_info;
    bool                remote_owned_closure;
    const BulkData* mesh;

    bool operator<(const EntityParallelState& rhs) const
    { return EntityLess(*mesh)(comm_info.entity, rhs.comm_info.entity); }
  };
#ifdef __CUDACC__
private:
#endif

  void communicate_entity_modification( const bool shared , std::vector<EntityParallelState > & data ); // Mod Mark
  bool pack_entity_modification( const bool packShared , stk::CommSparse & comm );

  virtual bool does_entity_need_orphan_protection(stk::mesh::Entity entity) const
  {
      const bool isNode = (stk::topology::NODE_RANK == entity_rank(entity));
      const bool isNotConnected = (1u == m_closure_count[entity.local_offset()]);
      const bool isOwned = bucket(entity).owned();
      const bool isCreatedState = (stk::mesh::Created == state(entity));
      return isNode && isNotConnected && isCreatedState && isOwned;
  }

  virtual bool does_entity_have_orphan_protection(stk::mesh::Entity entity) const
  {
      bool hasOrphanProtection = false;
      if (entity_rank(entity) == stk::topology::NODE_RANK && m_closure_count[entity.local_offset()] >= BulkData::orphaned_node_marking)
      {
          hasOrphanProtection = true;
      }
      return hasOrphanProtection;
  }

  // Only to be called from add_node_sharing
  void protect_orphaned_node(Entity entity)
  {
      if (does_entity_need_orphan_protection(entity))
      {
          internal_force_protect_orphaned_node(entity);
      }
  }

  void unprotect_orphaned_node(Entity entity)
  {
      if (does_entity_have_orphan_protection(entity))
      {
          internal_force_unprotect_orphaned_node(entity);
      }
  }

  void delete_sides_on_all_procs(const stk::mesh::EntityVector & deletedSides);
  void set_shared_owned_parts_and_ownership_on_comm_data(const std::vector<sharing_info>& shared_modified);
  void sync_parts_on_downward_related_entities_for_created_sides(const stk::mesh::EntityVector& locally_created_sides);
  void fill_entity_procs_for_owned_modified_or_created(std::vector<EntityProc> & send_list) const;
  stk::mesh::EntityVector get_lower_ranked_shared_entities(const stk::mesh::EntityVector& created_sides) const;

  stk::mesh::EntityVector get_nodes_to_deactivate(const stk::mesh::EntityVector & deactivatedElements, const stk::mesh::Part & activePart) const;

  inline bool internal_add_node_sharing_called() const;

  void internal_sync_comm_list_owners();

  // Forbidden
  BulkData();
  BulkData( const BulkData & );
  BulkData & operator = ( const BulkData & );

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
                                 OrdinalVector& scratch1, OrdinalVector& scratch2, OrdinalVector& scratch3); // Mod Mark

  int determine_new_owner( Entity ) const;  // Mod Mark

  template<typename PARTVECTOR>
  Entity internal_declare_entity( EntityRank ent_rank , EntityId ent_id ,
                                 const PARTVECTOR & parts ); // Mod Mark

  void internal_fill_parts_to_actually_remove(const OrdinalVector & removed,
                                              OrdinalVector &partsThatShouldStillBeInduced,
                                              OrdinalVector &delParts); // Mod Mark

  void internal_propagate_induced_part_changes_to_downward_connected_entities( Entity entity,
                                                                               const OrdinalVector & added,
                                                                               const OrdinalVector & removed ); // Mod Mark

  Ghosting & internal_create_ghosting( const std::string & name );
  void internal_verify_inputs_and_change_ghosting(
                                    Ghosting & ghosts ,
                                    const std::vector<EntityProc> & add_send ,
                                    const std::vector<EntityKey> & remove_receive );


  bool internal_modification_end_for_entity_creation( const std::vector<EntityRank> & entity_rank_vector,
                                                      impl::MeshModification::modification_optimization opt ); // Mod Mark


  void internal_establish_new_owner(stk::mesh::Entity entity);
  void internal_update_parts_for_shared_entity(stk::mesh::Entity entity, const bool is_entity_shared, const bool did_i_just_become_owner);

  inline void internal_check_unpopulated_relations(Entity entity, EntityRank rank) const;

  void internal_adjust_closure_count(Entity entity,
                                     const OrdinalVector & add_parts,
                                     const OrdinalVector & remove_parts); // Mod Mark

  void internal_fill_new_part_list_and_removed_part_list(stk::mesh::Entity entity,
                                                           const OrdinalVector & add_parts,
                                                           const OrdinalVector & remove_parts,
                                                           OrdinalVector &newBucketPartList,
                                                           OrdinalVector &parts_removed);
  void internal_move_entity_to_new_bucket(stk::mesh::Entity entity, const OrdinalVector &newBucketPartList,
                                          OrdinalVector& scratchSpace); // Mod Mark

  void internal_verify_change_parts( const MetaData   & meta ,
                                     const Entity entity ,
                                     const PartVector & parts ) const; // Mod Mark

  void require_entity_owner( const Entity entity, int owner) const ;

  bool is_good_rank_and_id(EntityRank ent_rank, EntityId ent_id) const;

  inline bool is_valid_connectivity(Entity entity, EntityRank rank) const;

  void require_valid_relation( const char action[] ,
                               const BulkData & mesh ,
                               const Entity e_from ,
                               const Entity e_to );

  friend class ::stk::mesh::MetaData;
  friend class ::stk::mesh::impl::BucketRepository;
  friend class stk::mesh::Bucket; // for field callback
  friend class Ghosting; // friend until Ghosting is refactored to be like Entity
  friend class ::stk::mesh::impl::MeshModification;
  friend class ::stk::mesh::FaceCreator;
  friend class ::stk::mesh::EntityLess;
  friend class ::stk::io::StkMeshIoBroker;
  friend class stk::mesh::DeviceMesh;

  // friends until it is decided what we're doing with Fields and Parallel and BulkData
  friend void communicate_field_data(const Ghosting & ghosts, const std::vector<const FieldBase *> & fields);
  friend void communicate_field_data(const BulkData & mesh, const std::vector<const FieldBase *> & fields);
  template <Operation Op> friend void parallel_op_including_ghosts_impl(const BulkData & mesh, const std::vector<const FieldBase *> & fields);
  friend void skin_mesh( BulkData & mesh, Selector const& element_selector, PartVector const& skin_parts, const Selector * secondary_selector);
  friend void create_edges( BulkData & mesh, const Selector & element_selector, Part * part_to_insert_new_edges );
  friend void internal_create_faces( BulkData & mesh, const Selector & element_selector, bool connect_faces_to_edges, FaceCreationBehavior faceCreationBehavior);
  friend bool process_killed_elements(stk::mesh::BulkData& bulkData,
                                        ElemElemGraph& elementGraph,
                                        const stk::mesh::EntityVector& killedElements,
                                        stk::mesh::Part& active,
                                        stk::mesh::impl::ParallelSelectedInfo &remoteActiveSelector,
                                        const stk::mesh::PartVector& parts_for_creating_side,
                                        const stk::mesh::PartVector* boundary_mesh_parts);
  friend stk::mesh::Entity clone_element_side(stk::mesh::BulkData &bulk,
                                                stk::mesh::EntityId id,
                                                stk::mesh::Entity elem,
                                                stk::mesh::ConnectivityOrdinal ord,
                                                const stk::mesh::PartVector &parts);

  friend stk::mesh::Entity impl::connect_side_to_element(stk::mesh::BulkData& bulkData, stk::mesh::Entity element,
                                                   stk::mesh::EntityId side_global_id, stk::mesh::ConnectivityOrdinal side_ordinal,
                                                   stk::mesh::Permutation side_permutation, const stk::mesh::PartVector& parts);

  void pack_owned_verify(CommSparse & commSparse);
  bool unpack_not_owned_verify(CommSparse & commSparse, std::ostream & error_log);
  bool verify_parallel_attributes( std::ostream & error_log );
  bool verify_parallel_attributes_for_bucket( Bucket const& bucket, std::ostream & error_log);
  void unpack_not_owned_verify_report_errors(Entity entity,
                                             bool bad_key,
                                             bool bad_own,
                                             bool bad_part,
                                             bool bad_rel,
                                             bool bad_comm,
                                             EntityKey recv_entity_key,
                                             int recv_owner_rank,
                                             PartVector const& recv_parts,
                                             std::vector<Relation> const& recv_relations,
                                             std::vector<int> const& recv_comm,
                                             std::ostream & error_log);
  void unpack_not_owned_verify_compare_comm_info( CommBuffer&            buf,
                                                  Entity                 entity,
                                                  EntityKey &            recv_entity_key,
                                                  int       &            recv_owner_rank,
                                                  unsigned  &            recv_comm_count,
                                                  PartVector &    recv_parts,
                                                  std::vector<Relation>& recv_relations,
                                                  std::vector<int>    &  recv_comm,
                                                  bool&                  bad_comm);


  void fillEntityCommInfoForEntity(stk::mesh::Ghosting &ghost_id,
                                   stk::mesh::BulkData &mesh,
                                   std::vector<stk::mesh::EntityKey> nodes,
                                   EntityCommInfoVector &sharing_processors);

  void determineEntitiesThatNeedGhosting(stk::mesh::BulkData &stkMeshBulkData,
                                         stk::mesh::Entity edge,
                                         std::vector<stk::mesh::Entity>& entitiesConnectedToNodes,
                                         const stk::mesh::Entity* nodes,
                                         std::set<EntityProc, EntityLess> &addGhostedEntities);

  void find_upward_connected_entities_to_ghost_onto_other_processors(stk::mesh::BulkData &mesh,
                                                                     std::set<EntityProc, EntityLess> &entitiesToGhostOntoOtherProcessors,
                                                                     EntityRank entity_rank,
                                                                     stk::mesh::Selector selected,
                                                                     bool connectFacesToPreexistingGhosts);

  void reset_add_node_sharing() { m_add_node_sharing_called = false; }

  void destroy_dependent_ghosts( Entity entity, EntityProcVec& entitiesToRemoveFromSharing );

  template<typename PARTVECTOR>
  Entity create_and_connect_side(const stk::mesh::EntityId globalSideId,
                                 Entity elem,
                                 const unsigned localSideId,
                                 const PARTVECTOR& parts);

  stk::mesh::EntityId select_side_id(Entity elem, unsigned sideOrdinal);
  void entity_comm_list_insert(Entity node);

public: // data
  mutable bool m_check_invalid_rels; // TODO REMOVE

protected: //data
  mutable ModificationNotifier notifier;
  static const uint16_t orphaned_node_marking;
  EntityCommDatabase m_entity_comm_map;
  std::vector<Ghosting*> m_ghosting;
  MetaData &m_mesh_meta_data;
  std::vector<EntitySharing> m_mark_entity; //indexed by Entity
  bool m_add_node_sharing_called;
  std::vector<uint16_t> m_closure_count; //indexed by Entity
  std::vector<MeshIndex> m_mesh_indexes; //indexed by Entity
  impl::EntityRepository* m_entity_repo;
  EntityCommListInfoVector m_entity_comm_list;
  std::vector<EntityComm*> m_entitycomm;
  std::vector<int> m_owner;
  CommListUpdater m_comm_list_updater;
  std::list<Entity::entity_value_type> m_deleted_entities_current_modification_cycle;
  GhostReuseMap m_ghost_reuse_map;
  std::vector<EntityKey> m_entity_keys; //indexed by Entity

#ifdef SIERRA_MIGRATION
  bool m_add_fmwk_data; // flag that will add extra data to buckets to support fmwk
  std::vector<FmwkId> m_fmwk_global_ids;
  mutable std::vector<RelationVector* > m_fmwk_aux_relations;   // Relations that can't be managed by STK such as PARENT/CHILD
  inline bool should_sort_faces_by_node_ids() const { return m_shouldSortFacesByNodeIds; }
  bool m_shouldSortFacesByNodeIds;
#endif
  enum AutomaticAuraOption m_autoAuraOption;
  stk::mesh::impl::MeshModification m_meshModification;

  void internal_force_protect_orphaned_node(stk::mesh::Entity entity)
  {
      m_closure_count[entity.local_offset()] += BulkData::orphaned_node_marking;
  }

  void internal_force_unprotect_orphaned_node(stk::mesh::Entity entity)
  {
      m_closure_count[entity.local_offset()] -= BulkData::orphaned_node_marking;
  }
  void set_mod_summary_proc_id(int proc_id) { m_modSummary.set_proc_id(proc_id); }

  Parallel m_parallel;

private: // data
  VolatileFastSharedCommMap m_volatile_fast_shared_comm_map;
  std::vector<std::vector<int> > m_all_sharing_procs;
  PartVector m_ghost_parts;
  std::list<Entity::entity_value_type> m_deleted_entities;
  int m_num_fields;
  bool m_keep_fields_updated;
  std::vector<unsigned> m_local_ids; //indexed by Entity

  //  ContiguousFieldDataManager m_default_field_data_manager;
  DefaultFieldDataManager m_default_field_data_manager;
  FieldDataManager *m_field_data_manager;
  mutable SelectorBucketMap m_selector_to_buckets_map;
  impl::BucketRepository m_bucket_repository; // needs to be destructed first!
  bool m_use_identifiers_for_resolving_sharing;
  std::string m_lastModificationDescription;
  stk::EmptyModificationSummary m_modSummary;
  // If needing debug info for modifications, comment out above line and uncomment line below
  //stk::ModificationSummary m_modSummary;
  std::shared_ptr<stk::mesh::MeshDiagnosticObserver> m_meshDiagnosticObserver;
  stk::mesh::ElemElemGraph* m_elemElemGraph = nullptr;
  std::shared_ptr<stk::mesh::ElemElemGraphUpdater> m_elemElemGraphUpdater;
  stk::mesh::impl::SideSetImpl<unsigned> m_sideSetData;
  mutable stk::mesh::NgpMeshManager* m_ngpMeshManager;
  mutable bool m_isDeviceMeshRegistered;

protected:
  stk::mesh::impl::SoloSideIdGenerator m_soloSideIdGenerator;
  bool m_supportsLargeIds = false;
};

void dump_mesh_info(const stk::mesh::BulkData& mesh, std::ostream&out, EntityVector ev);

void sync_to_host_and_mark_modified(MetaData& meta);

} // namespace mesh
} // namespace stk

#include "BulkDataInlinedMethods.hpp"


#endif //  stk_mesh_BulkData_hpp
