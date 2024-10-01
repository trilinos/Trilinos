// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
//
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
//
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
#include <stddef.h>                     // for size_t
#include <stdint.h>                     // for uint16_t
#include <algorithm>                    // for max
#include <stk_mesh/base/EntityIterator.hpp>
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
#include <functional>
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey, hash_value
#include "stk_mesh/base/FieldDataManager.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssert, etc
#include "stk_mesh/base/ModificationSummary.hpp"
#include <stk_mesh/base/ModificationNotifier.hpp>
#include "stk_mesh/baseImpl/MeshModification.hpp"
#include "stk_mesh/baseImpl/elementGraph/GraphTypes.hpp"
#include <stk_mesh/baseImpl/elementGraph/MeshDiagnosticObserver.hpp>
#include "NgpProfilingBlock.hpp"
#include "stk_mesh/baseImpl/SoloSideIdGenerator.hpp"
#include "stk_mesh/baseImpl/SideSetImpl.hpp"
#include "stk_mesh/baseImpl/DeviceMeshHostData.hpp"

namespace stk { namespace mesh { class EntityLess; } }
namespace stk { namespace mesh { class FieldBase; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { namespace impl { class AuraGhosting; } } }
namespace stk { namespace mesh { namespace impl { class EntityKeyMapping; } } }
namespace stk { namespace mesh { class FaceCreator; } }
namespace stk { namespace mesh { class ElemElemGraph; } }
namespace stk { namespace mesh { class ElemElemGraphUpdater; } }
namespace stk { namespace mesh { class MeshBuilder; } }
namespace stk { namespace mesh { class NgpMeshBase; } }
namespace stk { class CommSparse; }
namespace stk { namespace mesh { class ModificationObserver; } }
namespace stk { namespace io { class StkMeshIoBroker; } }
namespace stk { namespace mesh { namespace impl { struct RelationEntityToNode; } } }
namespace stk { namespace mesh { namespace impl { NgpMeshBase* get_ngp_mesh(const BulkData & bulk); } } }
namespace stk { namespace mesh { namespace impl { void set_ngp_mesh(const BulkData & bulk, NgpMeshBase * ngpMesh); } } }
namespace stk { namespace mesh { namespace impl { DeviceMeshHostData* get_ngp_mesh_host_data(const BulkData & bulk); } } }

namespace stk {
namespace tools {

template<typename T>
void replace_bulk_data(const stk::mesh::BulkData & inMesh, T & outMesh, std::function<void(T& outMesh)> op);

} }

#include "EntityCommListInfo.hpp"
#include "SharedEntityType.hpp"
#include "CommListUpdater.hpp"

namespace stk {
namespace mesh {
struct PartStorage;
struct SideSharingData;
enum class FaceCreationBehavior;
}}

namespace stk {
namespace mesh {

using ModEndOptimizationFlag = impl::MeshModification::modification_optimization;


void communicate_field_data(const Ghosting & ghosts, const std::vector<const FieldBase *> & fields);
void communicate_field_data(const BulkData & mesh, const std::vector<const FieldBase *> & fields);
void parallel_sum_including_ghosts(const BulkData & mesh, const std::vector<const FieldBase *> & fields);
void skin_mesh( BulkData & mesh, Selector const& element_selector, PartVector const& skin_parts, const Selector * secondary_selector);
void create_edges( BulkData & mesh, const Selector & element_selector, Part * part_to_insert_new_edges );
void internal_create_faces( BulkData & mesh, const Selector & element_selector, bool connect_faces_to_edges, FaceCreationBehavior faceCreationBehavior);
bool process_killed_elements(stk::mesh::BulkData& bulkData,
                             const stk::mesh::EntityVector& killedElements,
                             stk::mesh::Part& active,
                             stk::mesh::impl::ParallelSelectedInfo &remoteActiveSelector,
                             const stk::mesh::PartVector& parts_for_creating_side,
                             const stk::mesh::PartVector& boundary_mesh_parts,
                             ModEndOptimizationFlag deathModEndOpt);
stk::mesh::Entity clone_element_side(stk::mesh::BulkData &bulk,
                                     stk::mesh::EntityId id,
                                     stk::mesh::Entity elem,
                                     stk::mesh::ConnectivityOrdinal ord,
                                     const stk::mesh::PartVector &parts);
void remove_ghosts_from_remote_procs(stk::mesh::BulkData& bulk, EntityVector& recvGhostsToRemove);
void destroy_elements(stk::mesh::BulkData &bulk, stk::mesh::EntityVector &elementsToDestroy, const stk::mesh::Selector& orphansToDelete);

namespace impl {
stk::mesh::Entity connect_side_to_element(stk::mesh::BulkData& bulkData, stk::mesh::Entity element,
                                          stk::mesh::EntityId side_global_id, stk::mesh::ConnectivityOrdinal side_ordinal,
                                          stk::mesh::Permutation side_permutation, const stk::mesh::PartVector& parts);
}


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

  virtual ~BulkData();

  //------------------------------------
  /** \brief  The meta data manager for this bulk data manager. */
  const MetaData & mesh_meta_data() const { return *m_meta_data ; }
        MetaData & mesh_meta_data()       { return *m_meta_data ; }

  std::shared_ptr<MetaData> mesh_meta_data_ptr() {return m_meta_data; }
  const std::shared_ptr<MetaData> mesh_meta_data_ptr() const { return m_meta_data; }

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

  /** \brief set option for automatic maintenance of aura ghosting.
   *  If applyImmediately, then also create or destroy aura right now.
   *  Otherwise aura will be created or destroyed at the next modification_end.
   *  Synchronous/collective: must be called on all MPI ranks.
   */
  void set_automatic_aura_option(AutomaticAuraOption auraOption, bool applyImmediately=false);


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
  virtual
  bool modification_begin(const std::string description = std::string("UNSPECIFIED"))
  {
      ProfilingBlock block("mod begin:" + description);
      if(m_meshModification.in_modifiable_state()) {
        return false;
      }
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

  virtual
  bool modification_end(ModEndOptimizationFlag modEndOpt = ModEndOptimizationFlag::MOD_END_SORT)
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
  virtual bool change_entity_owner( const EntityProcVec & arg_change)
  {
      const bool anyEntitiesMoved = m_meshModification.change_entity_owner(arg_change);
      return anyEntitiesMoved;
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
  void update_field_data_states(bool rotateNgpFieldViews = false);
  void update_field_data_states(FieldBase *field, bool rotateNgpFieldViews = false);

  /** \brief  Copy field data from src entity to Dest entity
   *           - Fields that exist on the src that don't exist on the dest will
   *             be ignored
   *           - Fields that exist on the dest that don't exist on the src will
   *             be zeroed or initialized with the Field-specified initial-value.
   */
  void copy_entity_fields( Entity src, Entity dst);

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
  void add_node_sharing(Entity node, int sharing_proc );
  void add_node_sharing(const EntityProcVec& nodesAndProcs);

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

  void change_entity_parts(const Selector& selector,
                           EntityRank rank,
                           const PartVector& add_parts,
                           const PartVector& remove_parts = PartVector());

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
                            const std::vector<PartVector>& remove_parts,
                            ModEndOptimizationFlag opt = ModEndOptimizationFlag::MOD_END_SORT);

  void batch_change_entity_parts( const stk::mesh::EntityVector& entities, // Mod Mark
                            const PartVector& add_parts,
                            const PartVector& remove_parts,
                            ModEndOptimizationFlag opt = ModEndOptimizationFlag::MOD_END_SORT);

  void batch_change_entity_parts( const Selector& selector, // Mod Mark
                            EntityRank rank,
                            const PartVector& add_parts,
                            const PartVector& remove_parts,
                            ModEndOptimizationFlag opt = ModEndOptimizationFlag::MOD_END_SORT);

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

#ifndef STK_HIDE_DEPRECATED_CODE // Delete after July 31 2024
  STK_DEPRECATED_MSG("Use function in FindPermutation.hpp") Permutation find_permutation( const stk::topology &hr_entity_topo,
                                Entity const *higher_rank_entity_nodes,
                                const stk::topology &side_topo,
                                Entity const *side_nodes,
                                unsigned side_ordinal) const;

  STK_DEPRECATED_MSG("Use function in FindPermutation.hpp") bool check_permutation(Entity entity, Entity rel_entity, unsigned rel_ordinal, Permutation expected) const;
#endif

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

#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Oct 2024
  // Check if entity has a specific relation to an entity of subcell_rank
  STK_DEPRECATED bool relation_exist( const Entity entity, EntityRank subcell_rank, RelationIdentifier subcell_id );
#endif

#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Sept 2024
  STK_DEPRECATED inline VolatileFastSharedCommMapOneRank const& volatile_fast_shared_comm_map(EntityRank rank) const;
#endif
  inline HostCommMapIndices volatile_fast_shared_comm_map(EntityRank rank, int proc) const;

  const std::vector<int>& all_sharing_procs(stk::mesh::EntityRank rank) const;

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

  /** \brief send ghost entities to the specified entity-proc locations
   * return true if entities were ghosted. (false if specified ghosts were already ghosted)
   */
  virtual bool batch_add_to_ghosting(Ghosting &ghosting, const EntityProcVec &entitiesAndDestinationProcs);

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
  bool in_send_ghost( const Ghosting & ghosting, Entity entity, int proc) const;
  bool is_aura_ghosted_onto_another_proc( EntityKey key ) const;
  bool in_ghost( const Ghosting & ghost , EntityKey key , int proc ) const;
  bool in_ghost( const Ghosting & ghost , Entity entity , int proc ) const;
  bool in_ghost( const Ghosting & ghost , Entity entity ) const;
  void shared_procs_intersection(const std::vector<EntityKey> & keys, std::vector<int> & procs ) const;
  void shared_procs_intersection(const EntityVector & entities, std::vector<int> & procs ) const;

  // Comm-related convenience methods

  bool is_communicated_with_proc(Entity entity, int proc) const;
  void comm_procs(Entity entity, std::vector<int> & procs ) const;
  void comm_procs( const Ghosting & ghost , EntityKey key, std::vector<int> & procs ) const;
  void comm_shared_procs( EntityKey key, std::vector<int> & procs ) const;
  void comm_shared_procs( Entity entity, std::vector<int> & procs ) const;

  inline bool in_index_range(Entity entity) const;
  inline bool is_valid(Entity entity) const;
  size_t count_relations(Entity entity, bool onlyDownwardRelations = false) const;
  bool has_no_relations(Entity entity) const;

  inline const MeshIndex& mesh_index(Entity entity) const;
  inline MeshIndex& mesh_index(Entity entity);
  inline EntityId identifier(Entity entity) const;
  inline EntityRank entity_rank(Entity entity) const;
  inline EntityKey entity_key(Entity entity) const;
  inline EntityState state(Entity entity) const;
  inline Bucket & bucket(Entity entity) const;
  inline Bucket * bucket_ptr(Entity entity) const;
  inline unsigned bucket_ordinal(Entity entity) const;
  inline int parallel_owner_rank(Entity entity) const;
  inline unsigned local_id(Entity entity) const;

  inline void set_local_id(Entity entity, unsigned id);

#ifdef SIERRA_MIGRATION
  typedef int64_t FmwkId; //must be a signed type -- fmwk uses negative values sometimes
  inline FmwkId global_id(stk::mesh::Entity entity) const;
  inline void set_global_id(stk::mesh::Entity entity, FmwkId id);
  void initialize_global_ids();

  const RelationVector& aux_relations(Entity entity) const;
  RelationVector& aux_relations(Entity entity); // Mod Mark
  void reserve_relation(stk::mesh::Entity entity, const unsigned num); // Mod Mark
  void erase_and_clear_if_empty(stk::mesh::Entity entity, RelationIterator rel_itr); // Mod Mark
  void initialize_aux_relations();

  RelationIterator internal_begin_relation(Entity entity, const RelationType relation_type) const;
  RelationIterator internal_end_relation(Entity entity, const RelationType relation_type) const;
  void compress_relation_capacity(Entity entity);
  bool add_fmwk_data() const { return m_add_fmwk_data; }
#else
  bool add_fmwk_data() const { return false; }
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

  inline ConnectedEntities get_connected_entities(Entity entity, EntityRank rank) const;
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

  // memoized version
  BucketVector const& get_buckets(EntityRank rank, Selector const& selector) const;

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

  bool initialize_face_adjacent_element_graph();
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
  void create_side_entities(const SideSet &sideSet, const stk::mesh::PartVector& parts);

  bool supports_large_ids() const {return m_supportsLargeIds; }
  void set_large_ids_flag(bool largeIds) { m_supportsLargeIds = largeIds; }

  unsigned get_initial_bucket_capacity() const { return m_bucket_repository.get_initial_bucket_capacity(); }
  unsigned get_maximum_bucket_capacity() const { return m_bucket_repository.get_maximum_bucket_capacity(); }

  virtual bool does_entity_have_orphan_protection(stk::mesh::Entity entity) const
  {
      bool hasOrphanProtection = false;
      if (entity_rank(entity) == stk::topology::NODE_RANK && m_closure_count[entity.local_offset()] >= BulkData::orphaned_node_marking)
      {
          hasOrphanProtection = true;
      }
      return hasOrphanProtection;
  }

  bool is_mesh_consistency_check_on() const { return m_runConsistencyCheck; }

protected: //functions
  BulkData(std::shared_ptr<MetaData> mesh_meta_data,
           ParallelMachine parallel,
           enum AutomaticAuraOption auto_aura_option = AUTO_AURA,
#ifdef SIERRA_MIGRATION
           bool add_fmwk_data = false,
#endif
           std::unique_ptr<FieldDataManager> field_dataManager = std::unique_ptr<FieldDataManager>(),
           unsigned initialBucketCapacity = get_default_initial_bucket_capacity(),
           unsigned maximumBucketCapacity = get_default_maximum_bucket_capacity(),
           std::shared_ptr<impl::AuraGhosting> auraGhosting = std::shared_ptr<impl::AuraGhosting>(),
           bool createUpwardConnectivity = true);

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
                                                         stk::mesh::Part &activePart,
                                                         ModEndOptimizationFlag opt);

  void make_mesh_parallel_consistent_after_skinning(const std::vector<sharing_info>& sharedModified);

  bool modification_end_for_entity_creation( const std::vector<EntityRank> & entity_rank_vector,
                                             ModEndOptimizationFlag opt = ModEndOptimizationFlag::MOD_END_SORT); // Mod Mark

  bool internal_modification_end_for_skin_mesh( EntityRank entity_rank,
                                                ModEndOptimizationFlag opt,
                                                const stk::mesh::Selector& selectedToSkin,
                                                const stk::mesh::Selector * only_consider_second_element_from_this_selector); // Mod Mark

  bool inputs_ok_and_need_ghosting(Ghosting & ghosts ,
                               const std::vector<EntityProc> & add_send ,
                               const std::vector<Entity> & remove_receive,
                               std::vector<EntityProc> &filtered_add_send,
                               std::vector<Entity> &filtered_remove_receive);

  void internal_batch_add_to_ghosting(Ghosting &ghosting, const EntityProcVec &entitiesAndDestinationProcs); // Mod Mark

  void ghost_entities_and_fields(Ghosting & ghosting,
                                 EntityProcVec&& new_send,
                                 bool isFullRegen = false,
                                 const std::vector<EntityProc>& removedSendGhosts = std::vector<EntityProc>());

  void conditionally_add_entity_to_ghosting_set(const Ghosting &ghosting,
                                                Entity entity,
                                                int toProc,
                                                EntityProcVec& entitiesWithClosure);

  void add_closure_entities(const Ghosting& ghosting,
                            const EntityProcVec& entities,
                            EntityProcVec& entitiesWithClosure);

  const EntityCommListInfoVector & internal_comm_list() const { return m_entity_comm_list; }
  const EntityCommDatabase& internal_comm_db() const { return m_entity_comm_map; }

  PairIterEntityComm internal_entity_comm_map(const EntityKey & key) const { return m_entity_comm_map.comm(key); }

  PairIterEntityComm internal_entity_comm_map(Entity entity) const
  {
    if (m_entitycomm[entity.local_offset()] != -1) {
      return m_entity_comm_map.comm(m_entitycomm[entity.local_offset()]);
    }
    return PairIterEntityComm();
  }

  PairIterEntityComm internal_entity_comm_map(const EntityKey & key, const Ghosting & sub ) const { return m_entity_comm_map.comm(key,sub); }

  PairIterEntityComm internal_entity_comm_map(Entity entity, const Ghosting & sub ) const
  {
    if (m_entitycomm[entity.local_offset()] != -1) {
      return ghost_info_range(m_entity_comm_map.comm(m_entitycomm[entity.local_offset()]), sub.ordinal());
    }
    return PairIterEntityComm();
  }

  inline EntitySharing internal_is_entity_marked(Entity entity) const;
  PairIterEntityComm internal_entity_comm_map_shared(const EntityKey & key) const { return m_entity_comm_map.shared_comm_info(key); }

  PairIterEntityComm internal_entity_comm_map_shared(Entity entity) const
  {
    const int entityCommIndex = m_entitycomm[entity.local_offset()];
    if (entityCommIndex != -1) {
      return shared_comm_info_range(m_entity_comm_map.comm(entityCommIndex));
    }
    return PairIterEntityComm();
  }

  void markEntitiesForResolvingSharingInfoUsingNodes(stk::mesh::EntityRank entityRank, bool onlyConsiderSoloSides, std::vector<shared_entity_type>& shared_entities);
  virtual void sortNodesIfNeeded(std::vector<stk::mesh::EntityKey>& nodes);

  void gather_shared_nodes(std::vector<Entity> & shared_nodes);

  inline void set_state(Entity entity, EntityState entity_state);
  inline void set_entity_key(Entity entity, EntityKey key);
  std::pair<Entity, bool> internal_create_entity(EntityKey key, size_t preferred_offset = 0); // Mod Mark
  std::pair<Entity, bool> internal_get_or_create_entity_with_notification(EntityKey key, size_t preferred_offset = 0);


  /** \brief  Declare a collection of relations by simply iterating
   *          the input and calling declare_relation on each entry.
   */
  void internal_declare_relation( Entity entity, const RelationVector & rel,
                         OrdinalVector& scratch1); // Mod Mark

  bool internal_declare_relation(Entity e_from, Entity e_to,
                                 RelationIdentifier local_id, Permutation permut); // Mod Mark

  /** methods for managing arrays of entity member-data */

  inline void log_created_parallel_copy(Entity entity);

  void internal_change_entity_parts_without_propagating_to_downward_connected_entities(Entity entity, const OrdinalVector& add_parts, const OrdinalVector& remove_parts, OrdinalVector& parts_removed, OrdinalVector& newBucketPartList, OrdinalVector& scratchSpace);
  void internal_change_bucket_parts_without_propagating_to_downward_connected_entities(Bucket* bucket, EntityRank rank, const OrdinalVector& add_parts, const OrdinalVector& remove_parts, OrdinalVector& ranked_parts_removed, OrdinalVector& newBucketPartList);
  void internal_change_entity_parts_without_propagating_to_downward_connected_entities_with_notification(Entity entity, const OrdinalVector& add_parts, const OrdinalVector& remove_parts, OrdinalVector& parts_removed, OrdinalVector& newBucketPartList, OrdinalVector& scratchSpace);
  void internal_determine_inducible_parts(EntityRank e_rank, const OrdinalVector& add_parts, const OrdinalVector& parts_removed, OrdinalVector& inducible_parts_added, OrdinalVector& inducible_parts_removed);
  void internal_determine_inducible_parts_and_propagate_to_downward_connected_entities(
                                                        Entity entity,
                                                        const OrdinalVector& add_parts,
                                                        const OrdinalVector& parts_removed,
                                                        OrdinalVector & scratchOrdinalVec,
                                                        OrdinalVector & scratchSpace );

  void internal_determine_inducible_parts_and_propagate_to_downward_connected_entities(
                                                        Bucket* bucket,
                                                        const OrdinalVector& add_parts,
                                                        const OrdinalVector& parts_removed,
                                                        OrdinalVector & scratchOrdinalVec,
                                                        OrdinalVector & scratchSpace );

  /*  Entity modification consequences:
   *  1) Change entity relation => update via part relation => change parts
   *  2) Change parts => update forward relations via part relation
   *                  => update via field relation
   */
  void internal_change_entity_parts( Entity ,
                                     const OrdinalVector& add_parts ,
                                     const OrdinalVector& remove_parts,
                                     OrdinalVector& scratchOrdinalVec, OrdinalVector& scratchSpace); // Mod Mark

  void internal_change_entity_parts( const Selector& selector,
                                     const EntityRank rank,
                                     const OrdinalVector& add_parts ,
                                     const OrdinalVector& remove_parts,
                                     OrdinalVector& scratchOrdinalVec, OrdinalVector& scratchSpace); // Mod Mark

  bool internal_destroy_entity_with_notification(Entity entity, bool wasGhost = false); // Mod Mark
  virtual bool internal_destroy_entity(Entity entity, bool wasGhost = false);

  void internal_change_ghosting( Ghosting & ghosts,
                                 const std::vector<EntityProc> & add_send ,
                                 const std::vector<Entity> & remove_receive,
                                 bool add_send_is_globally_empty = false);

  void internal_add_to_ghosting( Ghosting &ghosting, const std::vector<EntityProc> &add_send); // Mod Mark

  template<typename PARTVECTOR>
  void internal_verify_and_change_entity_parts( Entity entity,
                                                const PARTVECTOR & add_parts ,
                                                const PARTVECTOR & remove_parts,
                                                OrdinalVector& scratchOrdinalVec,
                                                OrdinalVector& scratchSpace);

  template<typename PARTVECTOR>
  void internal_verify_and_change_entity_parts( const EntityVector& entities,
                                                const PARTVECTOR & add_parts ,
                                                const PARTVECTOR & remove_parts);

  template<typename PARTVECTOR>
  void internal_verify_and_change_entity_parts( const Selector& selector,
                                                const EntityRank rank,
                                                const PARTVECTOR & add_parts ,
                                                const PARTVECTOR & remove_parts);


  void internal_insert_all_parts_induced_from_higher_rank_entities_to_vector(stk::mesh::Entity entity,
                                                                               stk::mesh::Entity e_to,
                                                                               OrdinalVector &to_add);

  bool internal_modification_end_for_change_entity_owner(ModEndOptimizationFlag opt );
  bool internal_modification_end_for_change_parts(ModEndOptimizationFlag opt = ModEndOptimizationFlag::MOD_END_SORT);
  void internal_modification_end_for_change_ghosting();

  void mark_entity_and_upward_related_entities_as_modified(Entity entity, bool markUpDownClosureIfShared = false);

  void set_common_entity_key_and_fix_ordering_of_nodes_and_update_comm_map(std::vector<shared_entity_type> & shared_entity_map);
  void find_and_delete_internal_faces(stk::mesh::EntityRank entityRank,
                                      const stk::mesh::Selector *only_consider_second_element_from_this_selector); // Mod Mark

  void filter_upward_ghost_relations(const Entity entity, std::function<void(Entity)> filter);
  EntityVector get_upward_send_ghost_relations(const Entity entity);
  EntityVector get_upward_recv_ghost_relations(const Entity entity);
  void update_comm_list_based_on_changes_in_comm_map();

  void internal_resolve_shared_part_membership_for_element_death(); // Mod Mark

  void remove_unneeded_induced_parts(stk::mesh::Entity entity, PairIterEntityComm entity_comm_info,
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
  void internal_add_comm_list_entries(EntityCommListInfoVector& newCommListEntries);

  std::pair<int,bool> entity_comm_map_insert(Entity entity, const EntityCommInfo &val)
  {
      EntityKey key = entity_key(entity);
      std::pair<int,bool> result = m_entity_comm_map.insert(key, val, parallel_owner_rank(entity));
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
    m_entitycomm[entity.local_offset()] = -1;
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
  void internal_remove_aura();

  void require_ok_to_modify() const ;
#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Sept 2024
  void internal_update_fast_comm_maps() const;
#endif
  void internal_update_ngp_fast_comm_maps() const;
  void internal_update_all_sharing_procs() const;

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

  void entity_getter_debug_check(Entity entity) const
  {
    STK_ThrowAssertMsg(in_index_range(entity) , "Entity has out-of-bounds offset: " << entity.local_offset() << ", maximum offset is: " << m_entity_keys.size() - 1);
  }

  void require_good_rank_and_id(EntityRank ent_rank, EntityId ent_id) const;
  void remove_entity_callback(EntityRank rank, unsigned bucket_id, unsigned bucket_ord);

  bool internal_destroy_relation(Entity e_from ,
                                 Entity e_to,
                                 const RelationIdentifier local_id ); // Mod Mark

  void check_mesh_consistency();
  virtual void remove_entities_from_sharing(const EntityProcVec& entitiesToRemoveFromSharing, EntityVector & entitiesNoLongerShared);
  virtual void check_if_entity_from_other_proc_exists_on_this_proc_and_update_info_if_shared(std::vector<shared_entity_type>& shared_entity_map, int proc_id, const shared_entity_type &sentity);
  void update_owner_global_key_and_sharing_proc(stk::mesh::EntityKey global_key_other_proc,  shared_entity_type& shared_entity_this_proc, int proc_id) const;
  void update_shared_entity_this_proc(EntityKey global_key_other_proc, shared_entity_type& shared_entity_this_proc, int proc_id);
  void mark_shared_sides_and_fill_list_of_sides_not_on_boundary(std::vector<shared_entity_type>& shared_entity_map,
          int proc_id, shared_entity_type &sentity, std::vector<stk::mesh::EntityKeyProc> &entities_to_send_data,
          const stk::mesh::Selector *only_consider_second_element_from_this_selector);

  void fillVectorOfSharedEntitiesByProcessor(std::vector<shared_entity_type> & potentially_shared_sides, std::vector<std::vector<shared_entity_type> > &shared_entities_by_proc);

  void unpackEntityFromOtherProcAndUpdateInfoIfSharedLocally(stk::CommSparse &comm, std::vector<shared_entity_type> & shared_entity_map);

  virtual void change_entity_key_and_update_sharing_info(std::vector<shared_entity_type> & shared_entity_map);
  void change_entity_key_to_match_owner(const std::vector<shared_entity_type> & potentially_shared_sides);
  void insert_sharing_info_into_comm_map(const std::vector<shared_entity_type> & potentially_shared_sides);

  inline void internal_mark_entity(Entity entity, EntitySharing sharedType);

  void internal_change_entity_key(EntityKey old_key, EntityKey new_key, Entity entity); // Mod Mark

  void resolve_incremental_ghosting_for_entity_creation_or_skin_mesh(EntityRank entity_rank, const stk::mesh::Selector& selectedToSkin, bool connectFacesToPreexistingGhosts);

  void internal_finish_modification_end(ModEndOptimizationFlag opt); // Mod Mark

  void internal_set_owner(Entity entity, int new_owner);

  std::vector<uint64_t> internal_get_ids_in_use(stk::topology::rank_t rank, const std::vector<stk::mesh::EntityId>& reserved_ids = std::vector<stk::mesh::EntityId>()) const;

  virtual void de_induce_parts_from_nodes(const stk::mesh::EntityVector & deactivatedElements, stk::mesh::Part & activePart);

  virtual void remove_boundary_faces_from_part(stk::mesh::ElemElemGraph &graph,
                                         const stk::mesh::EntityVector & deactivatedElements,
                                         const stk::mesh::Part & activePart);

  virtual void internal_adjust_entity_and_downward_connectivity_closure_count(stk::mesh::Entity entity,
                                                                      stk::mesh::Bucket *bucket_old,
                                                                      int closureCountAdjustment); // Mod Mark

  inline void set_mesh_index(Entity entity, Bucket * in_bucket, unsigned ordinal );

  virtual void notify_finished_mod_end();

  void use_elem_elem_graph_to_determine_shared_entities(std::vector<stk::mesh::Entity>& shared_entities);
  void use_nodes_to_resolve_sharing(stk::mesh::EntityRank rank, std::vector<Entity>& shared_new, bool onlyConsiderSoloSides = false);
  void change_connectivity_for_edge_or_face(stk::mesh::Entity side, const std::vector<stk::mesh::EntityKey>& node_keys);
  void resolve_parallel_side_connections(std::vector<SideSharingData>& sideSharingDataToSend,
                                         std::vector<SideSharingData>& sideSharingDataReceived);
  void add_comm_map_for_sharing(const std::vector<SideSharingData>& sidesSharingData, stk::mesh::EntityVector& shared_entities);

private:
  void register_device_mesh() const;
  void unregister_device_mesh() const;

  void set_ngp_mesh(NgpMeshBase * ngpMesh) const { m_ngpMeshBase = ngpMesh; }
  NgpMeshBase * get_ngp_mesh() const { return m_ngpMeshBase; }

  void record_entity_deletion(Entity entity, bool isGhost);
  void break_boundary_relations_and_delete_buckets(const std::vector<impl::RelationEntityToNode> & relationsToDestroy, const stk::mesh::BucketVector & bucketsToDelete);
  void delete_buckets(const stk::mesh::BucketVector & buckets);
  void mark_entities_as_deleted(stk::mesh::Bucket * bucket);

  void filter_ghosting_remove_receives(const stk::mesh::Ghosting &ghosting, const std::vector <Entity> &remove_receive,
    std::vector<Entity> &recvGhosts, std::vector<bool>& ghostStatus);

  void verify_and_filter_add_send(Ghosting & ghosts, const std::vector<EntityProc> & add_send, bool &need_to_change_ghosting,
                                  bool &add_send_is_owned, std::vector <EntityProc> &filtered_add_send );

  void verify_and_filter_remove_receive(Ghosting & ghosts, const std::vector<Entity> & remove_receive, bool &need_to_change_ghosting, std::vector<Entity> & filtered_remove_receive);

  bool check_errors_and_determine_if_ghosting_needed_in_parallel(const stk::mesh::Ghosting &ghosts,
                                        bool add_send_is_owned,
                                        bool need_to_change_ghosting,
                                        const std::vector<EntityProc> & add_send);

  void delete_unneeded_entries_from_the_comm_list();

  void internal_resolve_sharing_and_ghosting_for_sides(bool connectFacesToPreexistingGhosts);

  virtual bool does_entity_need_orphan_protection(stk::mesh::Entity entity) const
  {
      const bool isNode = (stk::topology::NODE_RANK == entity_rank(entity));
      const bool isNotConnected = (1u == m_closure_count[entity.local_offset()]);
      const bool isOwned = bucket(entity).owned();
      const bool isCreatedState = (stk::mesh::Created == state(entity));
      return isNode && isNotConnected && isCreatedState && isOwned;
  }

  // Only to be called from add_node_sharing
  void protect_orphaned_node(Entity entity)
  {
      if (does_entity_need_orphan_protection(entity)) {
          internal_force_protect_orphaned_node(entity);
      }
  }

  void unprotect_orphaned_node(Entity entity)
  {
      if (does_entity_have_orphan_protection(entity)) {
          internal_force_unprotect_orphaned_node(entity);
      }
  }

  void delete_sides_on_all_procs(const stk::mesh::EntityVector & deletedSides);
  void set_shared_owned_parts_and_ownership_on_comm_data(const std::vector<sharing_info>& shared_modified);
  void sync_parts_on_downward_related_entities_for_created_sides(const stk::mesh::EntityVector& locally_created_sides);
  void fill_entity_procs_for_owned_modified_or_created(std::vector<EntityProc> & send_list) const;
  stk::mesh::EntityVector get_lower_ranked_shared_entities(const stk::mesh::EntityVector& created_sides) const;

  inline bool internal_add_node_sharing_called() const;

  // Forbidden
  BulkData();

  BulkData( const BulkData & );
  BulkData & operator = ( const BulkData & );

  void gather_and_print_mesh_partitioning() const;

  // Field callbacks

  void new_bucket_callback(EntityRank rank, const PartVector& superset_parts, size_t capacity, Bucket* new_bucket);
  void new_bucket_caching(EntityRank rank, Bucket* new_bucket);
  void remove_bucket_caching(EntityRank rank, const Selector selector);

  //
  //  "fields" is an optional argument, if present copy only the listed fields.
  //
  void copy_entity_fields_callback(EntityRank dst_rank, unsigned dst_bucket_id, unsigned dst_bucket_ord,
                                   unsigned src_bucket_id, unsigned src_bucket_ord,
                                   const std::vector<FieldBase*>* fields = nullptr);


  void destroy_bucket_callback(EntityRank rank, Bucket const& dying_bucket, unsigned capacity);

  // id_map, indexed by new id, maps to old id
  void reorder_buckets_callback(EntityRank rank, const std::vector<unsigned>& id_map);

  void remove_entity_field_data_callback(EntityRank rank, unsigned bucket_id, unsigned bucket_ord);
  void add_entity_callback(EntityRank rank, unsigned bucketId, unsigned bucketCapacity, unsigned indexInBucket);
  void reset_empty_field_data_callback(EntityRank rank, unsigned bucketId, unsigned bucketSize,
                                       unsigned bucketCapacity, const FieldVector & fields);

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

  void internal_propagate_induced_part_changes_to_downward_connected_entities(
                                                        Entity entity,
                                                        const OrdinalVector & added,
                                                        const OrdinalVector & removed,
                                                        OrdinalVector & scratchOrdinalVec,
                                                        OrdinalVector & scratchSpace );

  void internal_propagate_induced_part_changes_to_downward_connected_entities(
                                                        Bucket* bucket,
                                                        const OrdinalVector & addParts,
                                                        const OrdinalVector & removeParts,
                                                        OrdinalVector& scratchOrdinalPartsRemoved,
                                                        OrdinalVector& scratchOrdinalVec);

  Ghosting & internal_create_ghosting( const std::string & name );
  void internal_verify_inputs_and_change_ghosting(
                                    Ghosting & ghosts ,
                                    const std::vector<EntityProc> & add_send ,
                                    const std::vector<EntityKey> & remove_receive );


  bool internal_modification_end_for_entity_creation( const std::vector<EntityRank> & entity_rank_vector,
                                                      ModEndOptimizationFlag opt ); // Mod Mark


  inline void internal_check_unpopulated_relations(Entity entity, EntityRank rank) const;

  void internal_adjust_closure_count(Entity entity,
                                     const OrdinalVector & add_parts,
                                     const OrdinalVector & remove_parts); // Mod Mark

  void internal_adjust_closure_count(Bucket* bucket,
                                     const OrdinalVector & add_parts,
                                     const OrdinalVector & remove_parts); // Mod Mark

  void internal_fill_new_part_list_and_removed_part_list(stk::mesh::Bucket* bucket_old,
                                                           const OrdinalVector & add_parts,
                                                           const OrdinalVector & remove_parts,
                                                           OrdinalVector &newBucketPartList,
                                                           OrdinalVector &parts_removed);

  void internal_modify_bucket_parts(Bucket* bucket, const OrdinalVector& addParts, const OrdinalVector& removeParts);
  void internal_move_entity_to_new_bucket(stk::mesh::Entity entity, const OrdinalVector &newBucketPartList,
                                          OrdinalVector& scratchSpace); // Mod Mark

  void require_entity_owner( const Entity entity, int owner) const ;

  bool is_valid_connectivity(Entity entity, EntityRank rank) const;

  friend class MeshBuilder;
  friend class ::stk::mesh::MetaData;
  friend class ::stk::mesh::impl::BucketRepository;
  friend class stk::mesh::Bucket; // for field callback
  friend class Ghosting; // friend until Ghosting is refactored to be like Entity
  friend class ::stk::mesh::impl::MeshModification;
  friend class ::stk::mesh::impl::AuraGhosting;
  friend class ::stk::mesh::FaceCreator;
  friend class ::stk::mesh::EntityLess;
  friend class ::stk::io::StkMeshIoBroker;
  friend class stk::mesh::DeviceMesh;
  friend class stk::mesh::StkFieldSyncDebugger;
  template <typename T, template <typename> class NgpDebugger> friend class stk::mesh::DeviceField;

  // friends until it is decided what we're doing with Fields and Parallel and BulkData
  friend void communicate_field_data(const Ghosting & ghosts, const std::vector<const FieldBase *> & fields);
  friend void communicate_field_data(const BulkData & mesh, const std::vector<const FieldBase *> & fields);
  template <Operation Op> friend void parallel_op_including_ghosts_impl(const BulkData & mesh, const std::vector<const FieldBase *> & fields);
  friend void skin_mesh( BulkData & mesh, Selector const& element_selector, PartVector const& skin_parts, const Selector * secondary_selector);
  friend void create_edges( BulkData & mesh, const Selector & element_selector, Part * part_to_insert_new_edges );
  friend void internal_create_faces( BulkData & mesh, const Selector & element_selector, bool connect_faces_to_edges, FaceCreationBehavior faceCreationBehavior);
  friend bool process_killed_elements(stk::mesh::BulkData& bulkData,
                                        const stk::mesh::EntityVector& killedElements,
                                        stk::mesh::Part& active,
                                        stk::mesh::impl::ParallelSelectedInfo &remoteActiveSelector,
                                        const stk::mesh::PartVector& parts_for_creating_side,
                                        const stk::mesh::PartVector* boundary_mesh_parts,
                                        ModEndOptimizationFlag deathModEndOpt);
  friend stk::mesh::Entity clone_element_side(stk::mesh::BulkData &bulk,
                                                stk::mesh::EntityId id,
                                                stk::mesh::Entity elem,
                                                stk::mesh::ConnectivityOrdinal ord,
                                                const stk::mesh::PartVector &parts);

  friend void remove_ghosts_from_remote_procs(stk::mesh::BulkData& bulk, EntityVector& recvGhostsToRemove);
  friend void destroy_elements(stk::mesh::BulkData &bulk, stk::mesh::EntityVector &elementsToDestroy, const stk::mesh::Selector& orphansToDelete);

  friend stk::mesh::Entity impl::connect_side_to_element(stk::mesh::BulkData& bulkData, stk::mesh::Entity element,
                                                   stk::mesh::EntityId side_global_id, stk::mesh::ConnectivityOrdinal side_ordinal,
                                                   stk::mesh::Permutation side_permutation, const stk::mesh::PartVector& parts);
  friend NgpMeshBase * impl::get_ngp_mesh(const BulkData & bulk);
  friend void impl::set_ngp_mesh(const BulkData & bulk, NgpMeshBase * ngpMesh);
  friend impl::DeviceMeshHostData* impl::get_ngp_mesh_host_data(const BulkData & bulk);

  template<typename T>
  friend void stk::tools::replace_bulk_data(const stk::mesh::BulkData & in_mesh, T & out_mesh, std::function<void(T& outMesh)> op);

  void determineEntitiesThatNeedGhosting(stk::mesh::Entity edge,
                                         std::vector<stk::mesh::Entity>& entitiesConnectedToNodes,
                                         const stk::mesh::Entity* nodes,
                                         EntityProcVec& addGhostedEntities);

  void find_upward_connected_entities_to_ghost_onto_other_processors(EntityProcVec& entitiesToGhostOntoOtherProcessors,
                                                                     EntityRank entity_rank,
                                                                     const stk::mesh::Selector& selected,
                                                                     bool connectFacesToPreexistingGhosts);

  void reset_add_node_sharing() { m_add_node_sharing_called = false; }

  template<typename PARTVECTOR>
  Entity create_and_connect_side(const stk::mesh::EntityId globalSideId,
                                 Entity elem,
                                 const unsigned localSideId,
                                 const PARTVECTOR& parts);

  stk::mesh::EntityId select_side_id(Entity elem, unsigned sideOrdinal);
  void entity_comm_list_insert(Entity node);
  void init_mesh_consistency_check_mode();

public: // data
  mutable bool m_check_invalid_rels; // TODO REMOVE

protected: //data
  bool m_createUpwardConnectivity = true;
  std::shared_ptr<impl::AuraGhosting> m_auraGhosting;
  mutable ModificationNotifier notifier;
  static const uint16_t orphaned_node_marking;
  EntityCommDatabase m_entity_comm_map;
  std::vector<Ghosting*> m_ghosting;
  std::shared_ptr<MetaData> m_meta_data;
  std::vector<EntitySharing> m_mark_entity; //indexed by Entity
  bool m_add_node_sharing_called;
  std::vector<uint16_t> m_closure_count; //indexed by Entity
  std::vector<MeshIndex> m_mesh_indexes; //indexed by Entity
  std::unique_ptr<impl::EntityKeyMapping> m_entityKeyMapping;
  EntityCommListInfoVector m_entity_comm_list;
  std::vector<int> m_entitycomm;
  std::vector<int> m_owner;
  std::vector<std::pair<EntityKey,EntityCommInfo>> m_removedGhosts;
  CommListUpdater m_comm_list_updater;
  std::vector<EntityKey> m_entity_keys; //indexed by Entity
  std::unique_ptr<FieldDataManager> m_field_data_manager;
  impl::BucketRepository m_bucket_repository;

#ifdef SIERRA_MIGRATION
  bool m_add_fmwk_data; // flag that will add extra data to buckets to support fmwk
  std::vector<FmwkId> m_fmwk_global_ids;
  mutable std::vector<RelationVector* > m_fmwk_aux_relations;   // Relations that can't be managed by STK such as PARENT/CHILD
  inline bool should_sort_faces_by_node_ids() const { return m_shouldSortFacesByNodeIds; }
  bool m_shouldSortFacesByNodeIds;
#else
  inline bool should_sort_faces_by_node_ids() const { return false; }
#endif
  enum AutomaticAuraOption m_autoAuraOption;
  bool m_turningOffAutoAura;
  stk::mesh::impl::MeshModification m_meshModification;

  void internal_force_protect_orphaned_node(stk::mesh::Entity entity)
  {
      m_closure_count[entity.local_offset()] += BulkData::orphaned_node_marking;
  }

  void internal_force_unprotect_orphaned_node(stk::mesh::Entity entity)
  {
      m_closure_count[entity.local_offset()] -= BulkData::orphaned_node_marking;
  }

  Parallel m_parallel;

private: // data
  mutable VolatileFastSharedCommMap m_volatile_fast_shared_comm_map;
  mutable unsigned m_volatile_fast_shared_comm_map_sync_count;
  mutable std::shared_ptr<impl::DeviceMeshHostData> m_ngpMeshHostData;
  mutable std::vector<std::vector<int> > m_all_sharing_procs;
  mutable unsigned m_all_sharing_procs_sync_count;
  PartVector m_ghost_parts;
  int m_num_fields;
  bool m_keep_fields_updated;
  std::vector<unsigned> m_local_ids; //indexed by Entity

  mutable std::vector<SelectorBucketMap> m_selector_to_buckets_maps;
  bool m_use_identifiers_for_resolving_sharing;
  stk::EmptyModificationSummary m_modSummary;
  // If needing debug info for modifications, comment out above line and uncomment line below
  //stk::ModificationSummary m_modSummary;
  std::shared_ptr<stk::mesh::MeshDiagnosticObserver> m_meshDiagnosticObserver;
  stk::mesh::ElemElemGraph* m_elemElemGraph = nullptr;
  std::shared_ptr<stk::mesh::ElemElemGraphUpdater> m_elemElemGraphUpdater;
  stk::mesh::impl::SideSetImpl<unsigned> m_sideSetData;
  mutable stk::mesh::NgpMeshBase* m_ngpMeshBase;
  mutable bool m_isDeviceMeshRegistered;
  mutable size_t m_ngpFieldSyncBufferModCount;
  bool m_runConsistencyCheck;

protected:
  std::string m_lastModificationDescription;
  stk::mesh::impl::SoloSideIdGenerator m_soloSideIdGenerator;
  bool m_supportsLargeIds = false;
};


inline unsigned
BulkData::num_connectivity(Entity entity, EntityRank rank) const
{
  STK_ThrowAssert(bucket_ptr(entity));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->num_connectivity(mesh_idx.bucket_ordinal, rank);
}

inline unsigned
BulkData::find_ordinal(Entity entity, EntityRank rank, ConnectivityOrdinal ordinal) const
{
  STK_ThrowAssert(bucket_ptr(entity));
  const MeshIndex &mesh_idx = mesh_index(entity);
  unsigned num_rels = mesh_idx.bucket->num_connectivity(mesh_idx.bucket_ordinal, rank);
  ConnectivityOrdinal const *ords = mesh_idx.bucket->begin_ordinals(mesh_idx.bucket_ordinal, rank);
  STK_ThrowAssert(ords);

  unsigned i = 0;
  for (; i < num_rels; ++i)
  {
    if (ords[i] == ordinal)
      break;
  }
  return i;
}

inline ConnectedEntities
BulkData::get_connected_entities(Entity entity, EntityRank rank) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, rank));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->get_connected_entities(mesh_idx.bucket_ordinal, rank);
}

inline Entity const*
BulkData::begin(Entity entity, EntityRank rank) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, rank));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin(mesh_idx.bucket_ordinal, rank);
}

inline Entity const*
BulkData::begin_nodes(Entity entity) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, stk::topology::NODE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_nodes(mesh_idx.bucket_ordinal);
}

inline Entity const*
BulkData::begin_edges(Entity entity) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, stk::topology::EDGE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_edges(mesh_idx.bucket_ordinal);
}

inline Entity const*
BulkData::begin_faces(Entity entity) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, stk::topology::FACE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_faces(mesh_idx.bucket_ordinal);
}

inline Entity const*
BulkData::begin_elements(Entity entity) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, stk::topology::ELEMENT_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_elements(mesh_idx.bucket_ordinal);
}

inline ConnectivityOrdinal const*
BulkData::begin_ordinals(Entity entity, EntityRank rank) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, rank));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_ordinals(mesh_idx.bucket_ordinal, rank);
}

inline ConnectivityOrdinal const*
BulkData::begin_node_ordinals(Entity entity) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, stk::topology::NODE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_node_ordinals(mesh_idx.bucket_ordinal);
}

inline ConnectivityOrdinal const*
BulkData::begin_edge_ordinals(Entity entity) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, stk::topology::EDGE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_edge_ordinals(mesh_idx.bucket_ordinal);
}

inline ConnectivityOrdinal const*
BulkData::begin_face_ordinals(Entity entity) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, stk::topology::FACE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_face_ordinals(mesh_idx.bucket_ordinal);
}

inline ConnectivityOrdinal const*
BulkData::begin_element_ordinals(Entity entity) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, stk::topology::ELEMENT_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_element_ordinals(mesh_idx.bucket_ordinal);
}

inline Permutation const*
BulkData::begin_permutations(Entity entity, EntityRank rank) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, rank));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_permutations(mesh_idx.bucket_ordinal, rank);
}

inline Permutation const*
BulkData::begin_node_permutations(Entity entity) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, stk::topology::NODE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_node_permutations(mesh_idx.bucket_ordinal);
}

inline Permutation const*
BulkData::begin_edge_permutations(Entity entity) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, stk::topology::EDGE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_edge_permutations(mesh_idx.bucket_ordinal);
}

inline Permutation const*
BulkData::begin_face_permutations(Entity entity) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, stk::topology::FACE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_face_permutations(mesh_idx.bucket_ordinal);
}

inline Permutation const*
BulkData::begin_element_permutations(Entity entity) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, stk::topology::ELEMENT_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_element_permutations(mesh_idx.bucket_ordinal);
}

inline unsigned
BulkData::num_nodes(Entity entity) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, stk::topology::NODE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->num_nodes(mesh_idx.bucket_ordinal);
}

inline unsigned
BulkData::num_edges(Entity entity) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, stk::topology::EDGE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->num_edges(mesh_idx.bucket_ordinal);
}

inline unsigned
BulkData::num_faces(Entity entity) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, stk::topology::FACE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->num_faces(mesh_idx.bucket_ordinal);
}

inline unsigned
BulkData::num_elements(Entity entity) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, stk::topology::ELEMENT_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->num_elements(mesh_idx.bucket_ordinal);
}

inline Entity const*
BulkData::end(Entity entity, EntityRank rank) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, rank));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end(mesh_idx.bucket_ordinal, rank);
}

inline Entity const*
BulkData::end_nodes(Entity entity) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, stk::topology::NODE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_nodes(mesh_idx.bucket_ordinal);
}

inline Entity const*
BulkData::end_edges(Entity entity) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, stk::topology::EDGE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_edges(mesh_idx.bucket_ordinal);
}

inline Entity const*
BulkData::end_faces(Entity entity) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, stk::topology::FACE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_faces(mesh_idx.bucket_ordinal);
}

inline Entity const*
BulkData::end_elements(Entity entity) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, stk::topology::ELEMENT_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_elements(mesh_idx.bucket_ordinal);
}

inline ConnectivityOrdinal const*
BulkData::end_ordinals(Entity entity, EntityRank rank) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, rank));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_ordinals(mesh_idx.bucket_ordinal, rank);
}

inline ConnectivityOrdinal const*
BulkData::end_node_ordinals(Entity entity) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, stk::topology::NODE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_node_ordinals(mesh_idx.bucket_ordinal);
}

inline ConnectivityOrdinal const*
BulkData::end_edge_ordinals(Entity entity) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, stk::topology::EDGE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_edge_ordinals(mesh_idx.bucket_ordinal);
}

inline ConnectivityOrdinal const*
BulkData::end_face_ordinals(Entity entity) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, stk::topology::FACE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_face_ordinals(mesh_idx.bucket_ordinal);
}

inline ConnectivityOrdinal const*
BulkData::end_element_ordinals(Entity entity) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, stk::topology::ELEMENT_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_element_ordinals(mesh_idx.bucket_ordinal);
}

inline Permutation const*
BulkData::end_permutations(Entity entity, EntityRank rank) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, rank));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_permutations(mesh_idx.bucket_ordinal, rank);
}

inline Permutation const*
BulkData::end_node_permutations(Entity entity) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, stk::topology::NODE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_node_permutations(mesh_idx.bucket_ordinal);
}

inline Permutation const*
BulkData::end_edge_permutations(Entity entity) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, stk::topology::EDGE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_edge_permutations(mesh_idx.bucket_ordinal);
}

inline Permutation const*
BulkData::end_face_permutations(Entity entity) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, stk::topology::FACE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_face_permutations(mesh_idx.bucket_ordinal);
}

inline Permutation const*
BulkData::end_element_permutations(Entity entity) const
{
  STK_ThrowAssert(is_valid_connectivity(entity, stk::topology::ELEMENT_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_element_permutations(mesh_idx.bucket_ordinal);
}

inline bool
BulkData::has_permutation(Entity entity, EntityRank rank) const
{
  STK_ThrowAssert(bucket_ptr(entity));
  return bucket(entity).has_permutation(rank);
}

inline bool
BulkData::in_shared(Entity entity) const
{
  if (m_entitycomm[entity.local_offset()] != -1) {
    PairIterEntityComm commInfo = m_entity_comm_map.comm(m_entitycomm[entity.local_offset()]);
    return commInfo.front().ghost_id == SHARED;
  }
  return false;
}

inline bool
BulkData::in_receive_ghost( EntityKey key ) const
{
  const std::vector<Ghosting*> & ghosts= ghostings();
  for (size_t i=ghosts.size()-1;i>=AURA;--i)
  {
      if ( in_receive_ghost(*ghosts[i], key) )
          return true;
  }
  return false;
}

inline bool
BulkData::in_receive_ghost( Entity entity ) const
{
  const std::vector<Ghosting*> & ghosts= ghostings();
  for (size_t i=ghosts.size()-1;i>=AURA;--i)
  {
      if ( in_receive_ghost(*ghosts[i], entity) )
          return true;
  }
  return false;
}

inline bool
BulkData::in_receive_custom_ghost( EntityKey key ) const
{
  const std::vector<Ghosting*> & ghosts= ghostings();
  for (size_t i=ghosts.size()-1;i>AURA;--i)
  {
      if ( in_receive_ghost(*ghosts[i], key) )
          return true;
  }
  return false;
}

inline bool
BulkData::in_receive_ghost( const Ghosting & ghost , EntityKey key ) const
{
  const int owner_rank = parallel_owner_rank(get_entity(key));
  return in_ghost( ghost , key , owner_rank );
}

inline bool
BulkData::in_receive_ghost( const Ghosting & ghost , Entity entity ) const
{
  if (m_entitycomm[entity.local_offset()] == -1) {
    return false;
  }

  const int owner_rank = parallel_owner_rank(entity);
  if (owner_rank == parallel_rank()) {
    return false;
  }

  PairIterEntityComm entityCommInfo = m_entity_comm_map.comm(m_entitycomm[entity.local_offset()]);
  for(; !entityCommInfo.empty(); ++entityCommInfo) {
    if (entityCommInfo->ghost_id == ghost.ordinal()) {
      return true;
    }
  }

  return false;
}

inline bool
BulkData::in_send_ghost( EntityKey key) const
{
  return in_send_ghost(get_entity(key));
}

inline bool
BulkData::in_send_ghost( Entity entity) const
{
  const int owner_rank = parallel_owner_rank(entity);
  if (owner_rank == parallel_rank()) {
    for ( PairIterEntityComm ec = internal_entity_comm_map(entity); ! ec.empty() ; ++ec ) {
      if (ec->ghost_id != 0 && ec->proc != owner_rank) {
        return true;
      }
    }
  }

  return false;
}

#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Sept 2024
STK_DEPRECATED inline VolatileFastSharedCommMapOneRank const&
BulkData::volatile_fast_shared_comm_map(EntityRank rank) const
{
  STK_ThrowAssert(this->in_synchronized_state());
  STK_ThrowAssertMsg(rank < stk::topology::ELEMENT_RANK, "Cannot share entities of rank: " << rank);
  if (m_volatile_fast_shared_comm_map_sync_count < synchronized_count()) {
    internal_update_fast_comm_maps();
  }
  return m_volatile_fast_shared_comm_map[rank];
}
#endif

inline HostCommMapIndices
BulkData::volatile_fast_shared_comm_map(EntityRank rank, int proc) const
{
  STK_ThrowAssert(this->in_synchronized_state());
  STK_ThrowAssertMsg(rank < stk::topology::ELEMENT_RANK, "Cannot share entities of rank: " << rank);
  if (m_ngpMeshHostData == nullptr ||
      m_ngpMeshHostData->volatileFastSharedCommMapSyncCount < synchronized_count()) {
    internal_update_ngp_fast_comm_maps();
  }

  if (parallel_size() > 1)
  {
    const size_t dataBegin = m_ngpMeshHostData->hostVolatileFastSharedCommMapOffset[rank](proc);
    const size_t dataEnd = m_ngpMeshHostData->hostVolatileFastSharedCommMapOffset[rank](proc+1);
    return Kokkos::subview(m_ngpMeshHostData->hostVolatileFastSharedCommMap[rank], Kokkos::pair<size_t,size_t>(dataBegin, dataEnd));
  } else
  {
    return HostCommMapIndices("empty comm map indices", 0);
  }
}

inline Part&
BulkData::ghosting_part(const Ghosting& ghosting) const
{
  STK_ThrowRequireMsg(ghosting.ordinal() < m_ghost_parts.size(), "BulkData::ghosting_part ERROR, no part corresponds to ghosting with name="<<ghosting.name()<<" and ordinal="<<ghosting.ordinal());
  return *m_ghost_parts[ghosting.ordinal()];
}

inline bool
BulkData::in_index_range(Entity entity) const
{
  return entity.local_offset() < m_entity_keys.size();
}

inline bool
BulkData::is_valid(Entity entity) const
{
  STK_ThrowAssertMsg(in_index_range(entity),
                 "Error in stk::mesh::BulkData::is_valid, entity not in index range. "
                 " entity.local_offset()="<<entity.local_offset()<<", valid range is < "
                 << get_size_of_entity_index_space());
  return !m_meshModification.is_entity_deleted(entity.local_offset());
}

inline const MeshIndex&
BulkData::mesh_index(Entity entity) const
{
#ifndef NDEBUG
  entity_getter_debug_check(entity);
#endif

  return m_mesh_indexes[entity.local_offset()];
}

inline MeshIndex&
BulkData::mesh_index(Entity entity)
{
  STK_ThrowAssert(entity.local_offset() > 0); // setter check due to non-const

  return m_mesh_indexes[entity.local_offset()];
}

inline EntityId
BulkData::identifier(Entity entity) const
{
#ifndef NDEBUG
  entity_getter_debug_check(entity);
#endif

  return m_entity_keys[entity.local_offset()].id();
}

inline EntityRank
BulkData::entity_rank(Entity entity) const
{
#ifndef NDEBUG
  entity_getter_debug_check(entity);
#endif

  return m_entity_keys[entity.local_offset()].rank();
}

inline EntityKey
BulkData::entity_key(Entity entity) const
{
#ifndef NDEBUG
  entity_getter_debug_check(entity);
#endif

  return m_entity_keys[entity.local_offset()];
}

inline EntityState
BulkData::state(Entity entity) const
{
#ifndef NDEBUG
  entity_getter_debug_check(entity);
#endif
  return m_meshModification.get_entity_state(entity.local_offset());
}

inline void
BulkData::internal_mark_entity(Entity entity, EntitySharing sharedType)
{
    m_mark_entity[entity.local_offset()] = sharedType;
}

inline BulkData::EntitySharing
BulkData::internal_is_entity_marked(Entity entity) const
{
    return m_mark_entity[entity.local_offset()];
}

inline bool
BulkData::internal_add_node_sharing_called() const
{
  return m_add_node_sharing_called;
}

inline Bucket &
BulkData::bucket(Entity entity) const
{
  Bucket* bptr = bucket_ptr(entity);

  STK_ThrowAssertMsg(bptr != nullptr, "BulkData::bucket error, invalid bucket (nullptr) for entity with local_offset()="<<entity.local_offset());

  return *bptr;
}

inline Bucket *
BulkData::bucket_ptr(Entity entity) const
{
#ifndef NDEBUG
  entity_getter_debug_check(entity);
#endif

  return mesh_index(entity).bucket;
}

inline unsigned
BulkData::bucket_ordinal(Entity entity) const
{
#ifndef NDEBUG
  entity_getter_debug_check(entity);
#endif

  return mesh_index(entity).bucket_ordinal;
}

inline int
BulkData::parallel_owner_rank(Entity entity) const
{
#ifndef NDEBUG
  entity_getter_debug_check(entity);
#endif

  return m_owner[entity.local_offset()];
}

inline unsigned
BulkData::local_id(Entity entity) const
{
#ifndef NDEBUG
  entity_getter_debug_check(entity);
#endif

  return m_local_ids[entity.local_offset()];
}

#ifdef SIERRA_MIGRATION
inline BulkData::FmwkId
BulkData::global_id(stk::mesh::Entity entity) const
{
#ifndef NDEBUG
  entity_getter_debug_check(entity);
#endif

  return m_fmwk_global_ids[entity.local_offset()];
}

inline void
BulkData::set_global_id(stk::mesh::Entity entity, BulkData::FmwkId id)
{
  STK_ThrowAssert(entity.local_offset() > 0);

  m_modSummary.track_set_global_id(entity, id);

  m_fmwk_global_ids[entity.local_offset()] = id;
}

#endif

inline void
BulkData::set_mesh_index(Entity entity, Bucket * in_bucket, unsigned ordinal )
{
  STK_ThrowAssert(entity.local_offset() > 0);

  if (in_bucket != NULL) {
    STK_ThrowAssertMsg(in_bucket->size() >= ordinal, "Detected bad bucket/ordinal.");
  }
  MeshIndex &mesh_idx = mesh_index(entity);
  mesh_idx.bucket = in_bucket;
  mesh_idx.bucket_ordinal = ordinal;
}

inline void
BulkData::set_entity_key(Entity entity, EntityKey key)
{
  STK_ThrowAssert(entity.local_offset() > 0);

  m_entity_keys[entity.local_offset()] = key;
}

inline void
BulkData::set_state(Entity entity, EntityState entity_state)
{
  STK_ThrowAssert(entity.local_offset() > 0);
  m_meshModification.set_entity_state(entity.local_offset(), entity_state);
  m_mark_entity[entity.local_offset()] = NOT_MARKED;
}

inline void
BulkData::set_local_id(Entity entity, unsigned id)
{
  STK_ThrowAssert(entity.local_offset() > 0);

  m_local_ids[entity.local_offset()] = id;
}

namespace impl {
inline NgpMeshBase * get_ngp_mesh(const BulkData & bulk) {
  return bulk.get_ngp_mesh();
}

inline void set_ngp_mesh(const BulkData & bulk, NgpMeshBase * ngpMesh) {
  bulk.set_ngp_mesh(ngpMesh);
}

inline impl::DeviceMeshHostData* get_ngp_mesh_host_data(const BulkData& bulk)
{
  if (bulk.m_ngpMeshHostData == nullptr) {
    bulk.m_ngpMeshHostData = std::make_shared<impl::DeviceMeshHostData>();
  }
  return bulk.m_ngpMeshHostData.get();
}

}

} // namespace mesh
} // namespace stk

#endif //  stk_mesh_BulkData_hpp

