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

#ifndef ELEMENT_DISCONNECT_TOOL_HPP
#define ELEMENT_DISCONNECT_TOOL_HPP

#include <array>
#include <iterator>
#include <vector>

#include "DisjointSet.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_topology/topology.hpp"
#include "stk_util/util/SortAndUnique.hpp"

namespace stk::experimental
{

using EntityContainer = std::unordered_set<stk::mesh::Entity>;

// Returns a list of target-rank entities that are adjacent to the source entities via a given edge-rank such that each
// target entity satisfies a user-supplied predicate function.
template <topology::rank_t EdgeR, topology::rank_t TargetR, typename Pred>
auto get_adjacent_entities(const mesh::BulkData& bulk, const EntityContainer& srcEntities, Pred&& predicate)
{
  EntityContainer adjEntitySet;
  for (const auto& srcEntity : srcEntities) {
    const auto& sourceEdges = bulk.get_connected_entities(srcEntity, EdgeR);
    for (const auto& sEdge : sourceEdges) {
      const auto& targetEntities = bulk.get_connected_entities(sEdge, TargetR);
      for (const auto& tEntity : targetEntities) {
        if (std::forward<Pred>(predicate)(tEntity)) adjEntitySet.insert(tEntity);
      }
    }
  }
  return adjEntitySet;
}

struct FaceInfo {
  stk::mesh::Entity face{stk::mesh::Entity::InvalidEntity};
  stk::mesh::Entity elem{stk::mesh::Entity::InvalidEntity};
  unsigned side{};
};

struct FacePair {
  FaceInfo faceOrig{};
  FaceInfo faceNew{};
};

struct RelationTriplet
{
  RelationTriplet(const stk::mesh::Entity& entity_,
      const stk::mesh::Entity& origNode_,
      const stk::mesh::Entity& newNode_,
      const stk::mesh::ConnectivityOrdinal& ordinal_)
      : entity(entity_), origNode(origNode_), newNode(newNode_), ordinal(ordinal_)
  {
  }

  stk::mesh::Entity entity{};
  stk::mesh::Entity origNode{};
  stk::mesh::Entity newNode{};
  stk::mesh::ConnectivityOrdinal ordinal{};
};

class EntityDisconnectTool
{
 public:
  EntityDisconnectTool() { check_serial_execution(); }

  template <typename ContainerT>
  EntityDisconnectTool(stk::mesh::BulkData& mesh, const ContainerT& dcFaces = {})
      : m_disconnectFaces(dcFaces.begin(), dcFaces.end()), m_mesh(&mesh)
  {
    check_serial_execution();
    initialize();
  }

  virtual ~EntityDisconnectTool() = default;

  const auto& mesh() const { return *m_mesh; }
  auto& mesh() { return *m_mesh; }

  void set_disconnect_faces(const stk::mesh::EntityVector& dcFaces);
  void set_mesh(stk::mesh::BulkData& mesh) { m_mesh = &mesh; }

  const auto& get_disjoint_set() const { return m_disjointSet; }
  auto& get_disjoint_set() { return m_disjointSet; }
  const EntityContainer& get_retained_faces() const { return m_retainedFaces; }
  const EntityContainer& get_adjacent_elements() const { return m_adjacentElements; }
  const std::vector<FacePair>& get_elem_side_pairs() const { return m_disconnectFacePairs; }

  stk::mesh::EntityIdVector determine_new_nodes();
  void modify_mesh();

 private:
  static void check_serial_execution();

  void initialize();
  void filter_exterior_faces();
  void identify_adjacent_retained_faces();
  void identify_all_affected_faces();
  void identify_adjacent_elements();
  void add_elements_adjacent_to(const EntityContainer& faces, EntityContainer& adjacentElements);
  void identify_elem_side_pairs();

  void merge_connected_faces(const EntityContainer& adjacentFaces);

  stk::mesh::EntityIdVector compute_new_node_ids(const std::vector<NodeElemKey>& newNodeKeys);
  stk::mesh::EntityIdVector get_local_node_ids() const;
  std::vector<NodeElemKey> get_new_node_keys() const;

  void assign_new_ids_to_nodes(
      const std::vector<NodeElemKey>& newNodeKeys, const stk::mesh::EntityIdVector& newNodeIds);

  static bool is_interior(const stk::mesh::BulkData& mesh, const stk::mesh::Entity& face)
  {
    return mesh.num_connectivity(face, stk::topology::ELEMENT_RANK) > 1U;
  }

  bool is_disconnect_face(const stk::mesh::Entity& face) const;
  bool is_element_connected(const stk::mesh::Entity& face) const;

  void declare_new_nodes();
  std::vector<RelationTriplet> get_relation_triplets(const EntityContainer& adjacentEntities);
  void update_entity_nodal_relations(const EntityContainer& adjacentEntities);

  void update_face_relations();
  void disconnect_face_from_element(const FaceInfo& info);
  stk::mesh::Entity declare_new_face(FaceInfo& info);
  std::vector<stk::mesh::Entity> get_side_nodes(const FaceInfo& info) const;

  stk::mesh::Entity get_first_element(const stk::mesh::Entity& entity) const;

  DisjointSet m_disjointSet;
  EntityContainer m_disconnectFaces;
  EntityContainer m_retainedFaces;
  EntityContainer m_adjacentElements;
  EntityContainer m_allAffectedFaces;
  std::vector<FacePair> m_disconnectFacePairs;
  stk::mesh::BulkData* m_mesh;
};

}  // namespace stk::experimental

#endif
