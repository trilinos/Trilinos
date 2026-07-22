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

#include "stk_tools/mesh_tools/DisconnectGroup.hpp"
#include "stk_tools/mesh_tools/DisconnectUtils.hpp"
#include "stk_io/IossBridge.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_util/parallel/CommSparse.hpp"
#include <algorithm>
#include <map>

namespace stk {
namespace tools {
namespace impl {

DisconnectGroup::DisconnectGroup(const stk::mesh::BulkData& bulk)
  : m_bulk(bulk),
    m_parts(stk::mesh::ConstPartVector()),
    m_entities(stk::mesh::EntityVector()),
    m_node(stk::mesh::Entity()),
    m_active(true) {}

DisconnectGroup::DisconnectGroup(const stk::mesh::BulkData& bulk, const stk::mesh::Part* part, stk::mesh::Entity node)
  : m_bulk(bulk),
    m_parts(stk::mesh::ConstPartVector()),
    m_entities(stk::mesh::EntityVector()),
    m_node(node),
    m_active(true)
{
  STK_ThrowRequire(m_bulk.is_valid(m_node));

  if(part != nullptr) {
    m_parts.push_back(part);
  }

  m_entities = get_group_elements();
  store_node_sharing_info();
  update_id();
}

DisconnectGroup::DisconnectGroup(const stk::mesh::BulkData& bulk, const BlockPair& blockPair, stk::mesh::Entity node)
  : m_bulk(bulk),
    m_parts(stk::mesh::ConstPartVector()),
    m_entities(stk::mesh::EntityVector()),
    m_node(node),
    m_active(true),
    m_blockPair(blockPair),
    m_hasBlockPair(true)
{
  STK_ThrowRequire(m_bulk.is_valid(m_node));

  if(m_blockPair.second != nullptr) {
    m_parts.push_back(m_blockPair.second);
  }

  m_entities = get_group_elements();
  store_node_sharing_info();
  update_id();
}

DisconnectGroup::DisconnectGroup(const stk::mesh::BulkData& bulk, const stk::mesh::Part* part)
  : m_bulk(bulk),
    m_parts(stk::mesh::ConstPartVector()),
    m_entities(stk::mesh::EntityVector()),
    m_node(stk::mesh::Entity()),
    m_active(true)
{
  if(part != nullptr) {
    m_parts.push_back(part);
  }
}

DisconnectGroup::DisconnectGroup(const stk::mesh::BulkData& bulk, const stk::mesh::ConstPartVector& parts, stk::mesh::Entity node)
  : m_bulk(bulk),
    m_parts(parts),
    m_entities(stk::mesh::EntityVector()),
    m_node(node),
    m_active(true)
{
  STK_ThrowRequire(m_bulk.is_valid(m_node));

  m_entities = get_group_elements();
  store_node_sharing_info();
  update_id();
  std::sort(m_parts.begin(), m_parts.end(), stk::mesh::PartLess());
}

DisconnectGroup::DisconnectGroup(const DisconnectGroup& group)
  : m_bulk(group.m_bulk),
    m_parts(group.m_parts),
    m_entities(group.m_entities),
    m_node(group.m_node),
    m_active(group.m_active),
    m_id(group.m_id),
    m_entityOwnerProcVec(group.m_entityOwnerProcVec),
    m_blockPair(group.m_blockPair),
    m_hasBlockPair(group.m_hasBlockPair)
{
}

bool DisconnectGroup::operator<(const DisconnectGroup& group) const {
  if(m_parts.empty() || group.m_parts.empty()) { return false; }

  if(m_hasBlockPair) {
    return m_blockPair < group.m_blockPair;
  }

  return (m_parts[0]->mesh_meta_data_ordinal() < group.m_parts[0]->mesh_meta_data_ordinal());
}

bool DisconnectGroup::operator==(const DisconnectGroup& group) const {
  if(m_id != -1 && group.m_id != -1) {
    return (m_id == group.m_id) && (m_node == group.m_node);
  }

  if(m_hasBlockPair) {
    return (m_entities == group.m_entities);
  }

  return (m_parts == group.m_parts) && (m_entities == group.m_entities);
}

void DisconnectGroup::update_id() {
  if(m_entities.empty())   {
    m_id =  -1;
  } else {
    m_id =  m_bulk.identifier(m_entities[0]);
  }
}

std::vector<int> DisconnectGroup::get_sharing_procs(const stk::mesh::Part& part) const
{
  for(const EntityOwnerProc& entityOwnerProc : m_entityOwnerProcVec) {
    if(entityOwnerProc.part == &part) {
      return entityOwnerProc.ownerProcVec;
    }
  }
  return std::vector<int>();
}

void DisconnectGroup::insert_owner_info(stk::mesh::Part* part, int ownerProc)
{
  EntityOwnerProc* entry = nullptr;
  for(EntityOwnerProc& entityOwnerProc : m_entityOwnerProcVec) {
    if(entityOwnerProc.part == part) {
      entry = &entityOwnerProc;
      break;
    }
  }

  if(entry == nullptr) {
    EntityOwnerProc newEntry;
    newEntry.part = part;
    newEntry.ownerProcVec.push_back(ownerProc);
    m_entityOwnerProcVec.push_back(newEntry);
  } else {
    stk::util::insert_keep_sorted_and_unique(ownerProc, entry->ownerProcVec);
  }
}

void DisconnectGroup::store_node_sharing_info()
{
  unsigned numElems = m_bulk.num_elements(m_node);
  const stk::mesh::Entity* nodeElements = m_bulk.begin_elements(m_node);

  const stk::mesh::Bucket* prevBucketPtr = nullptr;
  stk::mesh::Part* part = nullptr;
  for(unsigned i = 0; i < numElems; i++) {
    const stk::mesh::Bucket* bucketPtr = m_bulk.bucket_ptr(nodeElements[i]);
    if (bucketPtr != prevBucketPtr) {
      part = get_block_part_for_bucket(m_bulk, *bucketPtr);
      prevBucketPtr = bucketPtr;
    }
    int elemOwner = m_bulk.parallel_owner_rank(nodeElements[i]);
    insert_owner_info(part, elemOwner);
  }
}

bool DisconnectGroup::needs_to_communicate() const
{
  const stk::mesh::Entity * elems = m_bulk.begin_elements(m_node);
  const unsigned numElems = m_bulk.num_elements(m_node);

  for (unsigned i = 0; i < numElems; ++i) {
    const stk::mesh::Bucket & elemBucket = m_bulk.bucket(elems[i]);
    for(const stk::mesh::Part* part : m_parts) {
      if (elemBucket.member(*part) && elemBucket.owned()) {
        return true;
      }
    }
  }
  return false;
}

void DisconnectGroup::update_info(stk::mesh::Entity node)
{
  m_node = node;
  STK_ThrowRequire(m_bulk.is_valid(node));
  m_entities = get_group_elements();
  store_node_sharing_info();
  update_id();
}

stk::mesh::EntityVector get_elements_for_node_in_parts(const stk::mesh::BulkData& bulk, const stk::mesh::Entity node, const stk::mesh::ConstPartVector& parts)
{
  stk::mesh::EntityVector elements;
  stk::mesh::EntityLess compare(bulk);
  unsigned numElems = bulk.num_elements(node);
  const stk::mesh::Entity* nodeElements = bulk.begin_elements(node);
  stk::mesh::Selector selector = stk::mesh::selectUnion(parts);

  for(unsigned i = 0; i < numElems; i++) {
    const stk::mesh::Bucket& elemBucket = bulk.bucket(nodeElements[i]);
    if(selector(elemBucket)) {
      elements.push_back(nodeElements[i]);
    }
  }

  std::sort(elements.begin(), elements.end(), compare);
  return elements;
}

stk::mesh::EntityVector DisconnectGroup::get_group_elements() const
{
  return get_elements_for_node_in_parts(m_bulk, m_node, m_parts);
}

stk::mesh::EntityIdVector DisconnectGroup::get_group_element_ids() const
{
  stk::mesh::EntityVector elements = get_group_elements();
  stk::mesh::EntityIdVector elementIds;

  for(stk::mesh::Entity element : elements) {
    elementIds.push_back(m_bulk.identifier(element));
  }
  return elementIds;
}
void DisconnectGroup::pack_group_info(stk::CommBuffer& procBuffer, stk::mesh::EntityId newNodeId, int /*proc*/) const {
  STK_ThrowRequire(!m_parts.empty());
  stk::mesh::EntityId parentNodeId = m_bulk.identifier(m_node);

  procBuffer.pack<stk::mesh::EntityId>(parentNodeId);
  procBuffer.pack<unsigned>(m_parts.size());

  for(const stk::mesh::Part* part : m_parts) {
    procBuffer.pack<stk::mesh::PartOrdinal>(part->mesh_meta_data_ordinal());
  }
  if(m_parts.size() > 1) {
    stk::mesh::EntityIdVector elements = get_group_element_ids();
    procBuffer.pack<unsigned>(elements.size());
    procBuffer.pack<stk::mesh::EntityId>(elements.data(), elements.size());
  }
  procBuffer.pack<stk::mesh::EntityId>(newNodeId);

  procBuffer.pack<bool>(m_hasBlockPair);
  if(m_hasBlockPair) {
    procBuffer.pack<stk::mesh::PartOrdinal>(m_blockPair.first->mesh_meta_data_ordinal());
    procBuffer.pack<stk::mesh::PartOrdinal>(m_blockPair.second->mesh_meta_data_ordinal());
  }
}

void DisconnectGroup::unpack_group_info(stk::CommBuffer& procBuffer, stk::mesh::EntityId& newNodeId, int /*proc*/) {
  unsigned numParts = 0;
  stk::mesh::PartOrdinal partOrdinal;
  stk::mesh::EntityId parentNodeId;

  procBuffer.unpack<stk::mesh::EntityId>(parentNodeId);

  m_node = m_bulk.get_entity(stk::topology::NODE_RANK, parentNodeId);
  STK_ThrowRequire(m_bulk.is_valid(m_node));

  procBuffer.unpack<unsigned>(numParts);
  STK_ThrowRequire(numParts != 0u);

  for(unsigned i = 0; i < numParts; i++) {
    procBuffer.unpack<stk::mesh::PartOrdinal>(partOrdinal);
    const stk::mesh::Part* part = &m_bulk.mesh_meta_data().get_part(partOrdinal);
    m_parts.push_back(part);
  }

  if(m_parts.size() > 1) {
    unsigned numElems = 0;
    procBuffer.unpack<unsigned>(numElems);
    stk::mesh::EntityIdVector elementIds(numElems);
    procBuffer.unpack<stk::mesh::EntityId>(elementIds.data(), numElems);

    for(stk::mesh::EntityId elementId : elementIds) {
      stk::mesh::Entity element = m_bulk.get_entity(stk::topology::ELEMENT_RANK, elementId);
      STK_ThrowRequire(m_bulk.is_valid(element));
      m_entities.push_back(element);
    }
  } else {
    m_entities = get_group_elements();
  }

  procBuffer.unpack<stk::mesh::EntityId>(newNodeId);

  procBuffer.unpack<bool>(m_hasBlockPair);
  if(m_hasBlockPair) {
    procBuffer.unpack<stk::mesh::PartOrdinal>(partOrdinal);
    stk::mesh::Part* part = &m_bulk.mesh_meta_data().get_part(partOrdinal);
    m_blockPair.first = part;
    procBuffer.unpack<stk::mesh::PartOrdinal>(partOrdinal);
    part = &m_bulk.mesh_meta_data().get_part(partOrdinal);
    m_blockPair.second = part;
  }
}

std::vector<int> get_elements_owners(const stk::mesh::BulkData& bulk, const stk::mesh::EntityVector& entities)
{
  std::vector<int> entityOwners;

  for(stk::mesh::Entity entity : entities) {
    stk::util::insert_keep_sorted_and_unique(bulk.parallel_owner_rank(entity), entityOwners);
  }
  return entityOwners;
}

std::vector<int> get_elements_owners_for_node_in_parts(const stk::mesh::BulkData& bulk, const stk::mesh::Entity node,
                                                       const stk::mesh::ConstPartVector& parts)
{
  stk::mesh::EntityVector entities = get_elements_for_node_in_parts(bulk, node, parts);
  return get_elements_owners(bulk, entities);
}

std::vector<int> DisconnectGroup::get_entity_owners() const
{
  std::vector<int> entityOwners = get_elements_owners(m_bulk, m_entities);
  return entityOwners;
}

}}}
