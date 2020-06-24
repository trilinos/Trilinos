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

#ifndef stk_mesh_base_impl_MeshImplUtils_hpp
#define stk_mesh_base_impl_MeshImplUtils_hpp

//----------------------------------------------------------------------

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/EntityProcMapping.hpp>
#include "stk_util/parallel/DistributedIndex.hpp"  // for DistributedIndex, etc

#include <vector>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {
namespace impl {

//----------------------------------------------------------------------
//These functions are not part of the public API of stk-mesh.
//They are intended for use internally in the implementation of
//stk-mesh capabilities.
//----------------------------------------------------------------------

void find_entities_these_nodes_have_in_common(const BulkData& mesh, stk::mesh::EntityRank rank, unsigned numNodes, const Entity* nodes, std::vector<Entity>& entity_vector);

void find_entities_with_larger_ids_these_nodes_have_in_common_and_locally_owned(stk::mesh::EntityId id, const BulkData& mesh, stk::mesh::EntityRank rank, unsigned numNodes, const Entity* nodes, std::vector<Entity>& entity_vector);

bool do_these_nodes_have_any_shell_elements_in_common(BulkData& mesh, unsigned numNodes, const Entity* nodes);

void find_locally_owned_elements_these_nodes_have_in_common(const BulkData& mesh, unsigned numNodes, const Entity* nodes, std::vector<Entity>& elems);

bool find_element_edge_ordinal_and_equivalent_nodes(BulkData& mesh, Entity element, unsigned numEdgeNodes, const Entity* edgeNodes, unsigned& elemEdgeOrdinal, Entity* elemEdgeNodes);

bool shared_entities_modified_on_any_proc(const BulkData& mesh, stk::ParallelMachine comm);

void get_ghost_data( const BulkData& bulkData, Entity entity, std::vector<EntityGhostData> & dataVector );

void connectUpwardEntityToEntity(stk::mesh::BulkData& mesh, stk::mesh::Entity upward_entity,
        stk::mesh::Entity entity, const stk::mesh::Entity* nodes);

void delete_upward_relations(stk::mesh::BulkData& bulkData,
                                             const stk::mesh::Entity& entity);

void delete_entities_and_upward_relations(stk::mesh::BulkData &bulkData, const stk::mesh::EntityVector &entities);

void internal_generate_parallel_change_lists( const BulkData & mesh ,
                                              const std::vector<EntityProc> & local_change ,
                                              std::vector<EntityProc> & shared_change ,
                                              std::vector<EntityProc> & ghosted_change );

stk::mesh::EntityVector convert_keys_to_entities(stk::mesh::BulkData &bulk, const std::vector<stk::mesh::EntityKey>& node_keys);

void internal_clean_and_verify_parallel_change(
  const BulkData & mesh ,
  std::vector<EntityProc> & local_change );

int check_no_shared_elements_or_higher(const BulkData& mesh);
int check_for_connected_nodes(const BulkData& mesh);
bool check_permutations_on_all(stk::mesh::BulkData& mesh);
void find_side_nodes(BulkData& mesh, Entity element, int side_ordinal, EntityVector & permuted_face_nodes);

inline
stk::mesh::EntityId side_id_formula(stk::mesh::EntityId elemId, unsigned sideOrdinal)
{
    //this is the side-id formula used by IO. the "+1" is because IO always uses one-based side ordinals
    return 10*elemId + sideOrdinal + 1;
}

class GlobalIdEntitySorter : public EntitySorterBase
{
public:
    virtual void sort(stk::mesh::BulkData &bulk, EntityVector& entityVector) const
    {
        std::sort(entityVector.begin(), entityVector.end(), EntityLess(bulk));
    }
};

template<class DO_THIS_FOR_ENTITY_IN_CLOSURE, class DESIRED_ENTITY>
void VisitClosureGeneral(
        const BulkData & mesh,
        Entity inputEntity,
        DO_THIS_FOR_ENTITY_IN_CLOSURE & do_this,
        DESIRED_ENTITY & desired_entity)
{
    if (desired_entity(inputEntity)) {
      do_this(inputEntity);
      const EntityRank inputEntityRank = mesh.entity_rank(inputEntity);
      for (EntityRank rank = stk::topology::NODE_RANK ; rank < inputEntityRank ; ++rank) {
          unsigned num_entities_of_rank = mesh.num_connectivity(inputEntity,rank);
          if (num_entities_of_rank > 0) {
            const bool dontRecurse = rank == stk::topology::NODE_RANK ||
                                     inputEntityRank <= stk::topology::ELEM_RANK;
            const Entity * entities = mesh.begin(inputEntity,rank);

            for (unsigned i=0 ; i<num_entities_of_rank ; ++i) {
                if (dontRecurse) {
                  if (desired_entity(entities[i])) {
                    do_this(entities[i]);
                  }
                }
                else {
                    VisitClosureGeneral(mesh,entities[i],do_this,desired_entity);
                }
            }
          }
      }
    }
}

template<class DO_THIS_FOR_ENTITY_IN_CLOSURE, typename FORWARD_ITERATOR, class DESIRED_ENTITY>
void VisitClosureGeneral(
        const stk::mesh::BulkData & mesh,
        const FORWARD_ITERATOR & start,
        const FORWARD_ITERATOR & finish,
        DO_THIS_FOR_ENTITY_IN_CLOSURE & do_this,
        DESIRED_ENTITY & desired_entity)
{
    for (FORWARD_ITERATOR entity_iterator = start ; entity_iterator != finish ; ++entity_iterator)
    {
        VisitClosureGeneral<DO_THIS_FOR_ENTITY_IN_CLOSURE,DESIRED_ENTITY>(mesh,*entity_iterator,do_this,desired_entity);
    }
}

template <typename VECTOR>
struct StoreInVector {
    StoreInVector(VECTOR & vec_in) : ev(vec_in) {}
    void operator()(stk::mesh::Entity entity) {
      ev.push_back(entity);
    }
    VECTOR & ev;
};

template <typename SET>
struct StoreInSet {
    StoreInSet(SET & set_in) : es(set_in) {}
    void operator()(stk::mesh::Entity entity) {
      es.insert(entity);
    }
    SET & es;
};

struct AlwaysVisit {
    bool operator()(Entity entity) { return true; }
};

struct OnlyVisitOnce {
    OnlyVisitOnce(const BulkData& mesh_in) : mesh(mesh_in) {}
    bool operator()(Entity entity) {
        if (mesh.is_valid(entity) && already_visited.find(entity) == already_visited.end()) {
            already_visited.insert(entity);
            return true;
        }
        return false;
    }
    const BulkData& mesh;
    std::set<Entity> already_visited;
};

struct OnlyVisitUnchanged
{
    OnlyVisitUnchanged(BulkData & mesh_in) : mesh(mesh_in) {}
    bool operator()(Entity entity){
        if (mesh.state(entity) == Unchanged) {
            return true;
        }
        return false;
    }
    BulkData & mesh;
};

struct OnlyVisitLocallyOwnedOnce {
    OnlyVisitLocallyOwnedOnce(const BulkData & mesh_in) : mesh(mesh_in), ovo(mesh_in) {}
    bool operator()(Entity entity)
    {
        return ovo(entity) && mesh.bucket(entity).owned();
    }
    const BulkData& mesh;
    OnlyVisitOnce ovo;
};

struct OnlyVisitSharedOnce {
    OnlyVisitSharedOnce(const BulkData & mesh_in) : mesh(mesh_in), ovo(mesh_in) {}
    bool operator()(Entity entity)
    {
        if (ovo(entity) && !mesh.in_shared(mesh.entity_key(entity))) { return true; }
        return false;
    }
    const BulkData & mesh;
    OnlyVisitOnce ovo;
};

struct OnlyVisitGhostsOnce
{
    OnlyVisitGhostsOnce(BulkData & mesh_in) : mesh(mesh_in), ovo(mesh_in) {}
    bool operator()(Entity entity) {
        if (ovo(entity) && mesh.in_receive_ghost(entity)) { return true; }
        return false;
    }
   BulkData & mesh;
   OnlyVisitOnce ovo;
};

template<class DO_THIS_FOR_ENTITY_IN_CLOSURE>
void VisitClosure(
        const stk::mesh::BulkData & mesh,
        stk::mesh::Entity entity_of_interest,
        DO_THIS_FOR_ENTITY_IN_CLOSURE & do_this)
{
    OnlyVisitOnce ovo(mesh);
    VisitClosureGeneral(mesh,entity_of_interest,do_this,ovo);
}


template<class DO_THIS_FOR_ENTITY_IN_CLOSURE, typename FORWARD_ITERATOR>
void VisitClosure(
        const stk::mesh::BulkData & mesh,
        const FORWARD_ITERATOR & start,
        const FORWARD_ITERATOR & finish,
        DO_THIS_FOR_ENTITY_IN_CLOSURE & do_this)
{
    OnlyVisitOnce ovo(mesh);
    VisitClosureGeneral(mesh,start,finish,do_this,ovo);
}


// cyclomatic complexity = 6
template<class DO_THIS_FOR_ENTITY_IN_UPWARD_CLOSURE, class DESIRED_ENTITY>
void VisitUpwardClosureGeneral(
        const BulkData & mesh,
        Entity entity_of_interest,
        DO_THIS_FOR_ENTITY_IN_UPWARD_CLOSURE & do_this,
        DESIRED_ENTITY & desired_entity)
{
    if (desired_entity(entity_of_interest)) {
        do_this(entity_of_interest);
        if (mesh.is_valid(entity_of_interest)) {
            EntityRank entity_of_interest_rank = mesh.entity_rank(entity_of_interest);
            EntityVector entities_of_rank_up;
            for (EntityRank rank_up = EntityRank(stk::topology::END_RANK-1) ; rank_up > entity_of_interest_rank ; --rank_up) {
                size_t num_entities_of_rank_up = mesh.num_connectivity(entity_of_interest,rank_up);
                const Entity * entity_up_it = mesh.begin(entity_of_interest,rank_up);

                for (size_t j=0 ; j<num_entities_of_rank_up ; ++j, ++entity_up_it) {
                    VisitUpwardClosureGeneral(mesh,*entity_up_it,do_this,desired_entity);
                }
            }
        }
    }
}

template<class DO_THIS_FOR_ENTITY_IN_UPWARD_CLOSURE, typename FORWARD_ITERATOR, class DESIRED_ENTITY>
void VisitUpwardClosureGeneral(
        const BulkData & mesh,
        const FORWARD_ITERATOR & start,
        const FORWARD_ITERATOR & finish,
        DO_THIS_FOR_ENTITY_IN_UPWARD_CLOSURE & do_this,
        DESIRED_ENTITY & desired_entity)
{
    for (FORWARD_ITERATOR entity_iterator = start ; entity_iterator != finish ; ++entity_iterator)
    {
        VisitUpwardClosureGeneral<DO_THIS_FOR_ENTITY_IN_UPWARD_CLOSURE,DESIRED_ENTITY>(mesh,*entity_iterator,do_this,desired_entity);
    }
}

template<class DO_THIS_FOR_ENTITY_IN_UPWARD_CLOSURE>
void VisitUpwardClosure(
        const BulkData & mesh,
        Entity entity_of_interest,
        DO_THIS_FOR_ENTITY_IN_UPWARD_CLOSURE & do_this)
{
    OnlyVisitOnce ovo(mesh);
    VisitUpwardClosureGeneral(mesh,entity_of_interest,do_this,ovo);
}

template<class DO_THIS_FOR_ENTITY_IN_UPWARD_CLOSURE, typename FORWARD_ITERATOR>
void VisitUpwardClosure(
        const BulkData & mesh,
        const FORWARD_ITERATOR & start,
        const FORWARD_ITERATOR & finish,
        DO_THIS_FOR_ENTITY_IN_UPWARD_CLOSURE & do_this)
{
    OnlyVisitOnce ovo(mesh);
    VisitUpwardClosureGeneral(mesh,start,finish,do_this,ovo);
}

template<class DO_THIS_FOR_ENTITY_IN_AURA_CLOSURE, typename FORWARD_ITERATOR, class DESIRED_ENTITY>
void VisitAuraClosureGeneral(
        const BulkData & mesh,
        const FORWARD_ITERATOR & start,
        const FORWARD_ITERATOR & finish,
        DO_THIS_FOR_ENTITY_IN_AURA_CLOSURE & do_this,
        DESIRED_ENTITY & desired_entity)
{
    std::set<Entity> entity_set;
    StoreInSet<std::set<Entity> > sis(entity_set);
    VisitClosure(mesh,start,finish,sis);
    VisitUpwardClosure(mesh,entity_set.begin(),entity_set.end(),sis);
    VisitClosureGeneral(mesh,entity_set.begin(),entity_set.end(),do_this,desired_entity);
}

template<class DO_THIS_FOR_ENTITY_IN_AURA_CLOSURE, class DESIRED_ENTITY>
void VisitAuraClosureGeneral(
        const BulkData & mesh,
        Entity entity_of_interest,
        DO_THIS_FOR_ENTITY_IN_AURA_CLOSURE & do_this,
        DESIRED_ENTITY & desired_entity)
{
    Entity * start = &entity_of_interest;
    Entity * finish = start+1;
    VisitAuraClosureGeneral(mesh,start,finish,do_this,desired_entity);
}

template<class DO_THIS_FOR_ENTITY_IN_AURA_CLOSURE, typename FORWARD_ITERATOR>
void VisitAuraClosure(
        const BulkData & mesh,
        const FORWARD_ITERATOR & start,
        const FORWARD_ITERATOR & finish,
        DO_THIS_FOR_ENTITY_IN_AURA_CLOSURE & do_this)
{
    OnlyVisitOnce ovo(mesh);
    VisitAuraClosureGeneral(mesh,start,finish,do_this,ovo);
}

template<class DO_THIS_FOR_ENTITY_IN_AURA_CLOSURE>
void VisitAuraClosure(
        const BulkData & mesh,
        Entity entity_of_interest,
        DO_THIS_FOR_ENTITY_IN_AURA_CLOSURE & do_this)
{
    OnlyVisitOnce ovo(mesh);
    VisitAuraClosureGeneral(mesh,entity_of_interest,do_this,ovo);
}

stk::parallel::DistributedIndex::KeySpanVector convert_entity_keys_to_spans( const MetaData & meta );

void get_part_ordinals_to_induce_on_lower_ranks_except_for_omits(const BulkData& mesh,
                             const Entity entity_from ,
                             const OrdinalVector       & omit ,
                             EntityRank            entity_rank_to ,
                             OrdinalVector       & induced_parts);
void get_part_ordinals_to_induce_on_lower_ranks(const BulkData& mesh,
                             const Entity entity_from ,
                             EntityRank            entity_rank_to ,
                             OrdinalVector       & induced_parts);

stk::mesh::Entity get_or_create_face_at_element_side(stk::mesh::BulkData & bulk,
                                                     stk::mesh::Entity elem,
                                                     int side_ordinal,
                                                     stk::mesh::EntityId new_face_global_id,
                                                     const stk::mesh::PartVector & parts = stk::mesh::PartVector());

template<typename PARTVECTOR>
stk::mesh::Entity connect_element_to_entity(stk::mesh::BulkData & mesh, stk::mesh::Entity elem, stk::mesh::Entity entity,
        const unsigned relationOrdinal, const PARTVECTOR& parts, stk::topology entity_top);

void connect_face_to_other_elements(stk::mesh::BulkData & bulk,
                                    stk::mesh::Entity face,
                                    stk::mesh::Entity elem_with_face,
                                    int elem_with_face_side_ordinal);


enum ShellStatus {
    NO_SHELLS = 25,
    YES_SHELLS_ONE_SHELL_ONE_SOLID = 31,
    YES_SHELLS_BOTH_SHELLS_OR_BOTH_SOLIDS = 46
};

void create_shell_status(const std::vector<stk::topology> & elements_touching_surface, stk::topology original_element_topology, std::vector<ShellStatus> & element_shell_status);

template<typename ENTITY_ID>
bool should_face_be_connected_to_element_side(std::vector<ENTITY_ID> & face_nodes,
                                              std::vector<ENTITY_ID> & element_side_nodes,
                                              stk::topology element_side_topology,
                                              ShellStatus  shell_status)
{
    bool should_connect = false;
    if(face_nodes.size() == element_side_nodes.size()) 
    {
        const stk::EquivalentPermutation equiv_result = element_side_topology.is_equivalent(face_nodes.data(), element_side_nodes.data());
        const bool nodes_match = equiv_result.is_equivalent;
        if (nodes_match) {
           if (NO_SHELLS == shell_status) {
               should_connect = true;
           }
           else {
               const unsigned permutation_of_element_side = equiv_result.permutation_number;
               const bool element_side_polarity_matches_face_nodes = permutation_of_element_side < element_side_topology.num_positive_permutations();
               if (YES_SHELLS_ONE_SHELL_ONE_SOLID == shell_status) {
                   should_connect = !element_side_polarity_matches_face_nodes;
               }
               else { // YES_SHELLS_BOTH_SHELS_OR_BOTH_SOLIDS
                   should_connect = element_side_polarity_matches_face_nodes;
               }
           }
        }
    }
    return should_connect;
}

struct StoreInEntityProcMapping {
    StoreInEntityProcMapping(BulkData & mesh_in, EntityProcMapping& epm_in)
    :mesh(mesh_in)
    ,myMapping(epm_in)
    {}

    void operator()(Entity entity) {
      myMapping.addEntityProc(entity, proc);
    }

    BulkData & mesh;
    EntityProcMapping& myMapping;
    int proc;
};

struct StoreInEntityProcSet {
    StoreInEntityProcSet(
            BulkData & mesh_in,
            std::set<stk::mesh::EntityProc, stk::mesh::EntityLess> & set_in)
    :mesh(mesh_in)
    ,myset(set_in)
    ,alreadyGhostedToProc(mesh_in.get_size_of_entity_index_space(), -1) { }

    void operator()(Entity entity) {
      if (proc != alreadyGhostedToProc[entity.local_offset()]) {
        alreadyGhostedToProc[entity.local_offset()] = proc;
        myset.insert(stk::mesh::EntityProc(entity,proc));
      }
    }

    BulkData & mesh;
    std::set<stk::mesh::EntityProc , stk::mesh::EntityLess> & myset;
    int proc;
    std::vector<int> alreadyGhostedToProc;
};

struct OnlyGhosts  {
    OnlyGhosts(BulkData & mesh_in) : mesh(mesh_in) {}
    bool operator()(Entity entity) {
      if (mesh.is_valid(entity)) {
        if (proc != mesh.parallel_owner_rank(entity)) {
          const bool isSharedWithProc = mesh.in_shared(entity, proc);
          return !isSharedWithProc;
        }
      }
      return false;
    }
    BulkData & mesh;
    int proc;
};

struct OnlyGhostsEPM  {
    OnlyGhostsEPM(BulkData & mesh_in, const EntityProcMapping& epm_in)
    : mesh(mesh_in), myMapping(epm_in) {}
    bool operator()(Entity entity) {
      if (mesh.is_valid(entity) && !myMapping.find(entity, proc)) {
        if (proc != mesh.parallel_owner_rank(entity)) {
          const bool isSharedWithProc = mesh.in_shared(entity, proc);
          return !isSharedWithProc;
        }
      }
      return false;
    }
    BulkData & mesh;
    const EntityProcMapping& myMapping;
    int proc;
};

struct OnlyNewGhosts  {
    OnlyNewGhosts(const BulkData & mesh_in, const Ghosting& ghosting_in) : mesh(mesh_in), ghosting(ghosting_in) {}
    bool operator()(Entity entity) {
      if (mesh.is_valid(entity)) {
        if (proc != mesh.parallel_owner_rank(entity)) {
          if (!mesh.in_ghost(ghosting, entity, proc)) {
            const bool isSharedWithProc = mesh.in_shared(entity, proc);
            return !isSharedWithProc;
          }
        }
      }
      return false;
    }
    const BulkData& mesh;
    const Ghosting& ghosting;
    int proc;
};

struct OnlyRecvGhosts {
  OnlyRecvGhosts(const BulkData& mesh_in, const Ghosting& ghost, const std::vector<bool>& status)
  : mesh(mesh_in), ghosting(ghost), ghostStatus(status) {}
  bool operator()(Entity entity) {
    return mesh.is_valid(entity) && mesh.in_receive_ghost(ghosting, entity) && !ghostStatus[entity.local_offset()];
  }
  const BulkData& mesh;
  const Ghosting& ghosting;
  const std::vector<bool>& ghostStatus;
};

struct VecPushBack {
  VecPushBack(std::vector<Entity>& rcvGhosts, std::vector<bool>& status)
  : recvGhosts(rcvGhosts), ghostStatus(status) {}
  void operator()(Entity entity) {
    recvGhosts.push_back(entity);
    ghostStatus[entity.local_offset()] = true;
  }
  std::vector<Entity>& recvGhosts;
  std::vector<bool>& ghostStatus;
};

void send_entity_keys_to_owners(
  BulkData & mesh ,
  const std::vector<Entity> & recvGhosts ,
        std::set< EntityProc , EntityLess > & sendGhosts);

void comm_sync_send_recv(
  BulkData & mesh ,
  std::set< EntityProc , EntityLess > & new_send ,
  std::set< EntityKey > & new_recv );

void comm_sync_send_recv(
  BulkData & mesh ,
  std::set< EntityProc , EntityLess > & new_send ,
  std::vector<Entity> & new_recv,
  std::vector<bool>& ghostStatus );

void comm_sync_aura_send_recv(
  BulkData & mesh ,
  std::vector<EntityProc> & sendGhosts,
  EntityProcMapping& entityProcMapping,
  std::vector<bool>& ghostStatus );

void insert_upward_relations(const BulkData& bulk_data, Entity rel_entity,
                             const EntityRank rank_of_orig_entity,
                             const int share_proc,
                             std::vector<EntityProc>& send);

void insert_upward_relations(const BulkData& bulk_data, Entity rel_entity,
                             const EntityRank rank_of_orig_entity,
                             const int share_proc,
                             EntityProcMapping& send);

void move_unowned_entities_for_owner_to_ghost(
  stk::mesh::BulkData & mesh ,
  std::set< stk::mesh::EntityProc , stk::mesh::EntityLess > & entitiesToGhostOntoOtherProcessors);

struct HashValueForEntityVector
{
    size_t operator()(const EntityVector &vector) const
    {
        size_t hashValue = 0;
        for(size_t i=0; i<vector.size(); i++)
        {
            hashValue += vector[i].local_offset();
        }
        return hashValue;
    }
};

void convert_part_ordinals_to_parts(const stk::mesh::MetaData& meta,
                                    const OrdinalVector& input_ordinals,
                                    stk::mesh::PartVector& output_parts);

void filter_out(OrdinalVector& vec,
                const OrdinalVector& parts,
                OrdinalVector& removed);

void merge_in(OrdinalVector& vec, const OrdinalVector& parts);

stk::mesh::ConnectivityOrdinal get_ordinal_from_side_entity(const std::vector<stk::mesh::Entity> &sides,
                                                            stk::mesh::ConnectivityOrdinal const * ordinals,
                                                            stk::mesh::Entity side);
stk::mesh::ConnectivityOrdinal get_ordinal_for_element_side_pair(const stk::mesh::BulkData &bulkData,
                                                                 stk::mesh::Entity element,
                                                                 stk::mesh::Entity side);

void fill_inducible_parts_from_list(const MetaData& meta,
                                    const OrdinalVector & partList,
                                    EntityRank rank,
                                    OrdinalVector &induciblePartsFromList);

void fill_part_list_differences(const BulkData &mesh,
                                Entity entity,
                                const PartVector &recv_parts,
                                std::set<std::string> &thisProcExtraParts,
                                std::set<std::string> &otherProcExtraParts);

void check_size_of_types();

void check_declare_element_side_inputs(const BulkData & mesh,
                                       const Entity elem,
                                       const unsigned localSideId);
} // namespace impl
} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif // stk_mesh_base_impl_MeshImplUtils_hpp

