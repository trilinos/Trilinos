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

#include "stk_tools/mesh_tools/EntityDisconnectTool.hpp"
#include <algorithm>
#include <iterator>

#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/DestroyRelations.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_tools/mesh_tools/DisjointSet.hpp"
#include "stk_topology/topology.hpp"
#include "stk_util/parallel/GenerateParallelConsistentIDs.hpp"
#include "stk_util/parallel/Parallel.hpp"
#include "stk_util/parallel/ParallelReduce.hpp"
#include "stk_util/util/ReportHandler.hpp"

namespace stk::experimental
{

void EntityDisconnectTool::check_serial_execution()
{
  STK_ThrowRequireMsg(stk::parallel_machine_size(stk::parallel_machine_world()) == 1,
      "The Entity Disconnect Tool does not yet support parallel execution.");
}

void EntityDisconnectTool::initialize()
{
  filter_exterior_faces();
  identify_adjacent_retained_faces();
  identify_all_affected_faces();
  identify_adjacent_elements();
  identify_elem_side_pairs();
}

void EntityDisconnectTool::filter_exterior_faces()
{
  const auto& bulk = mesh();
#if __cplusplus >= 202002L
  std::erase_if(m_disconnectFaces, [&bulk](const stk::mesh::Entity& face) { return !is_interior(bulk, face); });
#else
  for(auto it = m_disconnectFaces.begin(); it != m_disconnectFaces.end();) {
    if (!is_interior(bulk, *it)) {
      it = m_disconnectFaces.erase(it);
    }
    else {
      ++it;
    }
  }
#endif
}

void EntityDisconnectTool::identify_adjacent_retained_faces()
{
  static constexpr auto NodeR = topology::rank_t::NODE_RANK;
  static constexpr auto FaceR = topology::rank_t::FACE_RANK;
  m_retainedFaces =
      get_adjacent_entities<NodeR, FaceR>(mesh(), m_disconnectFaces, [this](const stk::mesh::Entity& connFace) {
        return !is_disconnect_face(connFace) && is_element_connected(connFace);
      });
}

void EntityDisconnectTool::identify_all_affected_faces() {
  m_allAffectedFaces.clear();
  m_allAffectedFaces.reserve(m_disconnectFaces.size() + m_retainedFaces.size());
  m_allAffectedFaces.insert(m_disconnectFaces.begin(), m_disconnectFaces.end());
  m_allAffectedFaces.insert(m_retainedFaces.begin(), m_retainedFaces.end());
}

bool EntityDisconnectTool::is_disconnect_face(const stk::mesh::Entity& face) const
{
  return m_disconnectFaces.find(face) != m_disconnectFaces.end();
}

bool EntityDisconnectTool::is_element_connected(const stk::mesh::Entity& face) const
{
  return mesh().num_connectivity(face, stk::topology::ELEMENT_RANK) >= 1;
}

void EntityDisconnectTool::identify_adjacent_elements()
{
  m_adjacentElements.clear();
  m_adjacentElements.reserve(2U * m_allAffectedFaces.size());
  add_elements_adjacent_to(m_allAffectedFaces, m_adjacentElements);
}

void EntityDisconnectTool::add_elements_adjacent_to(
    const EntityContainer& faces, EntityContainer& adjacentElements)
{
  for (const auto& adjFace : faces) {
    const auto dcElems = mesh().get_connected_entities(adjFace, stk::topology::ELEMENT_RANK);
    for (auto e = 0U; e < dcElems.size(); ++e) {
      m_adjacentElements.insert(dcElems[e]);
    }
  }
}

void EntityDisconnectTool::identify_elem_side_pairs()
{
  m_disconnectFacePairs.clear();
  m_disconnectFacePairs.reserve(m_disconnectFaces.size());
  for (const auto& face : m_disconnectFaces) {
    if (is_interior(mesh(), face)) {
      const auto faceElems = mesh().get_connected_entities(face, stk::topology::ELEMENT_RANK);
      const auto* faceOrds = mesh().begin_ordinals(face, stk::topology::ELEMENT_RANK);
      FaceInfo side0{face, faceElems[0], faceOrds[0]};
      FaceInfo side1{face, faceElems[1], faceOrds[1]};
      m_disconnectFacePairs.emplace_back(FacePair{side0, side1});
    }
  }
}

stk::mesh::EntityIdVector EntityDisconnectTool::determine_new_nodes()
{
  m_disjointSet.fill_set(mesh(), m_adjacentElements);
  merge_connected_faces(m_retainedFaces);

  auto newNodeKeys = get_new_node_keys();
  const auto newNodeIds = compute_new_node_ids(newNodeKeys);
  assign_new_ids_to_nodes(newNodeKeys, newNodeIds);

  return newNodeIds;
}

void EntityDisconnectTool::merge_connected_faces(const EntityContainer& adjacentFaces)
{
  for (const auto& connFace : adjacentFaces) {
    if (is_interior(mesh(), connFace)) {
      const auto& faceElems = mesh().get_connected_entities(connFace, stk::topology::ELEMENT_RANK);
      const auto& faceNodes = mesh().get_connected_entities(connFace, stk::topology::NODE_RANK);
      for (auto n = 0U; n < faceNodes.size(); ++n) {
        const NodeElemKey key0(faceNodes[n], faceElems[0]);
        const NodeElemKey key1(faceNodes[n], faceElems[1]);
        m_disjointSet.merge_nodes(key0, key1);
      }
    }
  }
}

stk::mesh::EntityIdVector EntityDisconnectTool::compute_new_node_ids(const std::vector<NodeElemKey>& newNodeKeys)
{
  const auto localNodeIds = get_local_node_ids();
  const auto globalNewNodes = stk::get_global_sum(stk::parallel_machine_world(), newNodeKeys.size());
  const auto globalMaxNode =
      stk::get_global_max(stk::parallel_machine_world(), *std::max_element(localNodeIds.begin(), localNodeIds.end()));
  const auto maxEntityId = globalMaxNode + globalNewNodes + 1;
  return stk::generate_parallel_consistent_ids(maxEntityId, localNodeIds, newNodeKeys, stk::parallel_machine_world());
}

stk::mesh::EntityIdVector EntityDisconnectTool::get_local_node_ids() const
{
  const auto& bulk = mesh();
  const auto& meta = bulk.mesh_meta_data();
  stk::mesh::Selector owned = meta.universal_part() & meta.locally_owned_part();
  const auto nodeBuckets = bulk.get_buckets(stk::topology::NODE_RANK, owned);
  stk::mesh::EntityVector localNodes;
  for (const auto* bucket : nodeBuckets) {
    localNodes.insert(localNodes.end(), bucket->begin(), bucket->end());
  }
  stk::mesh::EntityIdVector localNodeIds;
  localNodeIds.reserve(localNodes.size());
  std::transform(localNodes.begin(), localNodes.end(), std::back_inserter(localNodeIds),
      [&bulk](const stk::mesh::Entity& node) { return bulk.identifier(node); });
  return localNodeIds;
}

std::vector<NodeElemKey> EntityDisconnectTool::get_new_node_keys() const
{
  std::vector<NodeElemKey> newNodeKeys;
  for (const auto& [key, djNode] : m_disjointSet) {
    if (djNode.parent == nullptr && djNode.nodeId == stk::mesh::InvalidEntityId) {
      newNodeKeys.push_back(key);
    }
  }
  return newNodeKeys;
}

void EntityDisconnectTool::assign_new_ids_to_nodes(
    const std::vector<NodeElemKey>& newNodeKeys, const std::vector<stk::mesh::EntityId>& newNodeIds)
{
  STK_ThrowRequireMsg(newNodeKeys.size() == newNodeIds.size(),
      "The number of new node IDs is not equal to the number of new nodes requested!");
  for (auto n = 0U; n < newNodeKeys.size(); ++n) {
    auto& djNode = m_disjointSet[newNodeKeys[n]];
    djNode.nodeId = newNodeIds[n];
    djNode.isNew = true;
  }
}

void EntityDisconnectTool::modify_mesh()
{
  mesh().modification_begin();
  declare_new_nodes();
  update_entity_nodal_relations(m_adjacentElements);
  update_face_relations();
  mesh().modification_end();
}

void EntityDisconnectTool::declare_new_nodes()
{
  auto& djSet = get_disjoint_set();
  for (auto& [djKey, djNode] : djSet) {
    if (djNode.isNew) {
      const auto parentEntity = djNode.origEntity;
      const auto& parentBucket = mesh().bucket(parentEntity);
      const auto newNode = mesh().declare_entity(stk::topology::NODE_RANK, djNode.nodeId, parentBucket.supersets());
      djNode.entity = newNode;
      mesh().copy_entity_fields(parentEntity, newNode);

      /* bulk.add_node_sharing(Entity node, int sharing_proc) */
    }
  }
}

std::vector<stk::mesh::Entity> EntityDisconnectTool::get_side_nodes(const FaceInfo& info) const
{
  const auto* elemNodes = mesh().begin_nodes(info.elem);
  const auto elemTopo = mesh().bucket(info.elem).topology();
  std::vector<stk::mesh::Entity> sideNodes(elemTopo.side_topology().num_nodes());
  elemTopo.side_nodes(elemNodes, info.side, sideNodes.data());
  return sideNodes;
}

stk::mesh::Entity EntityDisconnectTool::get_first_element(const stk::mesh::Entity& entity) const
{
  if (mesh().entity_rank(entity) == stk::topology::ELEMENT_RANK) {
    return entity;
  } else {
    if (mesh().num_connectivity(entity, stk::topology::ELEMENT_RANK) < 1U) {
      std::stringstream errorMsg;
      errorMsg << "Entity " << mesh().identifier(entity) << "(" << mesh().entity_rank(entity)
               << ") must be connected to at least one element.";
      STK_ThrowAssertMsg(false, errorMsg.str());
    }
    return mesh().begin(entity, stk::topology::ELEMENT_RANK)[0];
  }
}

std::vector<RelationTriplet> EntityDisconnectTool::get_relation_triplets(const EntityContainer& adjacentEntities)
{
  std::vector<RelationTriplet> relationTriplets;
  // TODO: reserve vector with sensible size
  //   Should be djTrees | filter(new trees) | tree_size | accumulate
  const auto& disjointSet = get_disjoint_set();
  for (const auto& entity : adjacentEntities) {
    const auto elem = get_first_element(entity);
    const auto numNodes = mesh().num_connectivity(entity, stk::topology::NODE_RANK);
    const auto* entityNodes = mesh().begin(entity, stk::topology::NODE_RANK);
    const auto* entityOrdinals = mesh().begin_ordinals(entity, stk::topology::NODE_RANK);
    for (unsigned n = 0; n < numNodes; ++n) {
      const auto& djNode = disjointSet.find_root(NodeElemKey(entityNodes[n], elem));
      if (djNode.isNew) {
        STK_ThrowAssert(entityNodes[n] == djNode.origEntity);
        const auto localOrd = entityOrdinals[n];
        relationTriplets.emplace_back(entity, djNode.origEntity, djNode.entity, localOrd);
      }
    }
  }
  return relationTriplets;
}

void EntityDisconnectTool::update_entity_nodal_relations(const EntityContainer& adjacentEntities)
{
  const auto relationTriplets = get_relation_triplets(adjacentEntities);
  for (auto&& [entity, oldEnt, newEnt, localOrd] : relationTriplets) {
    mesh().destroy_relation(entity, oldEnt, localOrd);
  }
  for (auto&& [entity, oldEnt, newEnt, localOrd] : relationTriplets) {
    mesh().declare_relation(entity, newEnt, localOrd);
  }
}

void EntityDisconnectTool::update_face_relations()
{
  for (auto& [sideOrig, sideNew] : m_disconnectFacePairs) {
    // Somewhat confusingly, face entity of "New" is a copy of the original face at this stage.  Using info of "New" as
    // it corresponds to the opposite element of "Orig".
    disconnect_face_from_element(sideNew);
  }
  for (auto& [sideOrig, sideNew] : m_disconnectFacePairs) {
    auto newFace = declare_new_face(sideNew);
    sideNew.face = newFace;
  }
  update_entity_nodal_relations(m_allAffectedFaces);
}

void EntityDisconnectTool::disconnect_face_from_element(const FaceInfo& info)
{
  const auto& jElem = info.elem;
  const auto& jSide = info.side;
  mesh().destroy_relation(jElem, info.face, jSide);
}

stk::mesh::Entity EntityDisconnectTool::declare_new_face(FaceInfo& info)
{
  const auto& parts = mesh().bucket(info.face).supersets();
  const auto faceId = stk::mesh::impl::side_id_formula(mesh().identifier(info.elem), info.side);
  const auto newFace = mesh().declare_entity(stk::topology::FACE_RANK, faceId, parts);
  const auto sideNodes = get_side_nodes(info);
  for (std::size_t nIndex = 0; nIndex < sideNodes.size(); ++nIndex) {
    mesh().declare_relation(newFace, sideNodes[nIndex], nIndex);
  }
  mesh().declare_relation(info.elem, newFace, info.side);
  mesh().copy_entity_fields(info.face, newFace);

  return newFace;
}

}  // namespace stk::experimental
