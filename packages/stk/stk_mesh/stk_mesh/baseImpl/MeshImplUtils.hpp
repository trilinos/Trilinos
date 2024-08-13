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
#include <stk_mesh/base/EntityLess.hpp>
#include <stk_mesh/base/EntityProcMapping.hpp>
#include <stk_mesh/base/EntitySorterBase.hpp>
#include "stk_util/parallel/DistributedIndex.hpp"  // for DistributedIndex, etc

#include <vector>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {
namespace impl {

struct EntityGhostData;

//----------------------------------------------------------------------
//These functions are not part of the public API of stk-mesh.
//They are intended for use internally in the implementation of
//stk-mesh capabilities.
//----------------------------------------------------------------------

inline
bool is_in_list(Entity entity, const Entity* begin, const Entity* end) 
{
  return std::find(begin, end, entity) != end; 
}

template<typename VecType>
void
remove_entities_not_in_list(const Entity* beginList,
                            const Entity* endList,
                            VecType& entities)
{
  int numFound=0;
  for(int j=0, initialSize=entities.size(); j<initialSize; ++j) {
    if (is_in_list(entities[j], beginList, endList)) {
      if (j > numFound) {
        entities[numFound] = entities[j];
      }
      ++numFound;
    }
  }    
  entities.resize(numFound);
}

template<typename VecType>
void
remove_entities_not_connected_to_other_nodes(const BulkData& mesh,
                                             stk::mesh::EntityRank rank,
                                             unsigned numNodes,
                                             const Entity* nodes,
                                             VecType& elementsInCommon)
{
  for(unsigned i = 1; i < numNodes; ++i) {
    const ConnectedEntities conn = mesh.get_connected_entities(nodes[i], rank);
    remove_entities_not_in_list(conn.data(), conn.data()+conn.size(), elementsInCommon);
  }    
}

template<typename VecType>
void
find_entities_these_nodes_have_in_common(const BulkData& mesh,
                                         stk::mesh::EntityRank rank,
                                         unsigned numNodes,
                                         const Entity* nodes,
                                         VecType& entitiesInCommon)
{
  entitiesInCommon.clear();
  if(numNodes > 0) {
    const ConnectedEntities conn = mesh.get_connected_entities(nodes[0], rank);
    entitiesInCommon.assign(conn.data(), conn.data()+conn.size());
    remove_entities_not_connected_to_other_nodes(mesh, rank, numNodes, nodes, entitiesInCommon);
  }
}

void find_entities_with_larger_ids_these_nodes_have_in_common_and_locally_owned(stk::mesh::EntityId id, const BulkData& mesh, stk::mesh::EntityRank rank, unsigned numNodes, const Entity* nodes, std::vector<Entity>& entity_vector);

template<class Pred>
void
find_entities_these_nodes_have_in_common_and(const BulkData& mesh, EntityRank rank,
                                             unsigned numNodes, const Entity* nodes,
                                             std::vector<Entity>& entitiesInCommon,
                                             const Pred& pred)
{
  entitiesInCommon.clear();
  if(numNodes > 0) {
    const ConnectedEntities conn = mesh.get_connected_entities(nodes[0], rank);
    entitiesInCommon.reserve(conn.size());
    for(unsigned i=0; i<conn.size(); ++i) {
      if (pred(conn[i])) {
        entitiesInCommon.push_back(conn[i]);
      }
    }

    remove_entities_not_connected_to_other_nodes(mesh, rank, numNodes, nodes, entitiesInCommon);
  }
}

const EntityCommListInfo& find_entity(const BulkData& mesh,
                                      const EntityCommListInfoVector& entities,
                                      const EntityKey& key);

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

bool internal_clean_and_verify_parallel_change(
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

stk::parallel::DistributedIndex::KeySpanVector convert_entity_keys_to_spans( const MetaData & meta );

void get_part_ordinals_to_induce_on_lower_ranks_except_for_omits(const BulkData& mesh,
                             const Bucket& bucket_from ,
                             const OrdinalVector       & omit ,
                             EntityRank            entity_rank_to ,
                             OrdinalVector       & induced_parts);

void get_part_ordinals_to_induce_on_lower_ranks(const BulkData& mesh,
                                                const Bucket& bucket_from ,
                                                EntityRank entity_rank_to ,
                                                OrdinalVector& induced_parts);

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

void send_entity_keys_to_owners(
  BulkData & mesh ,
  const std::vector<Entity> & recvGhosts ,
        std::set< EntityProc , EntityLess > & sendGhosts);

void comm_sync_send_recv(const BulkData& mesh,
                         EntityProcVec& sendGhosts,
                         std::set<EntityKey>& recvGhosts);

void comm_sync_send_recv(const BulkData & mesh ,
                         const std::vector<Entity>& removeRecvGhosts,
                         EntityProcVec& newSendGhosts,
                         std::set< EntityKeyProc> & removeSendGhosts);

void comm_sync_aura_send_recv(
  BulkData & mesh ,
  std::vector<EntityProc> & sendGhosts,
  EntityProcMapping& entityProcMapping,
  std::vector<bool>& ghostStatus );

void comm_sync_nonowned_sends(
  const BulkData & mesh ,
  std::vector<EntityProc> & nonOwnedSendGhosts,
  EntityProcMapping& entityProcMapping);

void insert_upward_relations_for_owned(const BulkData& bulk_data,
                             const Entity entity,
                             const EntityRank entityRank,
                             const EntityRank maxRank,
                             const std::vector<int>& share_proc,
                             EntityProcMapping& send);

void move_unowned_entities_for_owner_to_ghost(BulkData & mesh, EntityProcVec& sendGhosts);

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

bool are_any_parts_ranked(const stk::mesh::MetaData& meta,
                          const OrdinalVector& partOrdinals);

void filter_out(OrdinalVector& vec,
                const OrdinalVector& parts,
                OrdinalVector& removed,
                bool trackRemoved = true);

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

void require_valid_relation(const char action[],
                            const BulkData& mesh,
                            const Entity e_from,
                            const Entity e_to);

bool is_valid_relation(const BulkData& mesh,
                       Entity e_from,
                       Entity e_to,
                       EntityRank e_to_rank,
                       ConnectivityOrdinal ord);

bool is_good_rank_and_id(const MetaData& meta,
                         EntityRank rank,
                         EntityId id);

EntityId get_global_max_id_in_use(const BulkData& mesh,
                                  EntityRank rank,
                                  const std::vector<Entity::entity_value_type>& deletedEntitiesCurModCycle,
                                  const std::vector<EntityId>& reservedIds = std::vector<EntityId>());

void check_declare_element_side_inputs(const BulkData & mesh,
                                       const Entity elem,
                                       const unsigned localSideId);

bool connect_edge_to_elements(stk::mesh::BulkData& bulk, stk::mesh::Entity edge);
void connect_face_to_elements(stk::mesh::BulkData& bulk, stk::mesh::Entity face);

bool has_upward_recv_ghost_connectivity(const stk::mesh::BulkData &bulk,
                                        const stk::mesh::Ghosting& ghosting,
                                        stk::mesh::Entity entity);
bool has_upward_send_ghost_connectivity(const stk::mesh::BulkData &bulk,
                                        const stk::mesh::Ghosting& ghosting,
                                        int proc,
                                        stk::mesh::Entity entity);
bool has_upward_connectivity(const stk::mesh::BulkData &bulk, stk::mesh::Entity entity);

bool can_destroy_entity(const stk::mesh::BulkData &bulk, stk::mesh::Entity entity);

void destroy_upward_connected_aura_entities(stk::mesh::BulkData &bulk,
                                            stk::mesh::Entity connectedEntity,
                                            EntityVector& scratchSpace);

void print_upward_connected_entities(stk::mesh::BulkData& bulk,
                                     stk::mesh::Entity entity,
                                     std::ostream& os);
} // namespace impl
} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif // stk_mesh_base_impl_MeshImplUtils_hpp

