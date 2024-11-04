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
#ifndef _DisconnectGroup_hpp_
#define _DisconnectGroup_hpp_

#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/Part.hpp"
#include "stk_mesh/base/SideSetEntry.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_tools/mesh_tools/DisconnectTypes.hpp"
#include "stk_util/parallel/ParallelComm.hpp"
#include <utility>
#include <vector>
#include <ostream>

namespace stk {
namespace mesh { class BulkData; }
namespace tools {
namespace impl {

struct EntityOwnerProc {
  stk::mesh::Part* part;
  std::vector<int> ownerProcVec;

  bool operator<(const EntityOwnerProc& rhs) const {
    return part->mesh_meta_data_ordinal() < rhs.part->mesh_meta_data_ordinal();
  }
};

class DisconnectGroup
{
public:
  DisconnectGroup(const stk::mesh::BulkData& bulk);

  DisconnectGroup(const stk::mesh::BulkData& bulk, const stk::mesh::Part* part, stk::mesh::Entity node);

  DisconnectGroup(const stk::mesh::BulkData& bulk, const BlockPair& blockPair, stk::mesh::Entity node);

  DisconnectGroup(const stk::mesh::BulkData& bulk, const stk::mesh::Part* part);

  DisconnectGroup(const stk::mesh::BulkData& bulk, const stk::mesh::ConstPartVector& parts, stk::mesh::Entity node);

  DisconnectGroup(const DisconnectGroup& group);

  bool operator<(const DisconnectGroup& group) const ;

  bool operator==(const DisconnectGroup& group) const;

  int id() const {return m_id;}

  bool needs_to_communicate() const;

  void update_info(stk::mesh::Entity node);

  const stk::mesh::ConstPartVector& get_parts() const { return m_parts; }
  const stk::mesh::EntityVector& get_entities() const { return m_entities; }

  stk::mesh::Entity get_node() const { return m_node; }
  stk::mesh::EntityId get_node_id() const { return m_bulk.identifier(m_node); }
  const stk::mesh::BulkData& get_bulk() const { return m_bulk; }

  std::vector<int> get_entity_owners() const;

  stk::mesh::EntityVector get_group_elements() const;

  stk::mesh::EntityIdVector get_group_element_ids() const;

  void pack_group_info(stk::CommBuffer& procBuffer, stk::mesh::EntityId newNodeId, int proc) const;

  void unpack_group_info(stk::CommBuffer& procBuffer, stk::mesh::EntityId& newNodeId, int proc);

  void set_active(bool flag) const { m_active = flag; }

  bool is_active() const { return m_active; }

  bool has_block_pair() const { return m_hasBlockPair; }

  BlockPair get_block_pair() const { return m_blockPair; }

  const std::vector<EntityOwnerProc> get_part_owner_proc() const { return m_entityOwnerProcVec; }
  std::vector<int> get_sharing_procs(const stk::mesh::Part& part) const;

private:
  void update_id();
  void insert_part_owner_proc(EntityOwnerProc proc) { m_entityOwnerProcVec.push_back(proc); }
  void store_node_sharing_info();
  void insert_owner_info(stk::mesh::Part* part, int ownerProc);

  const stk::mesh::BulkData& m_bulk;
  stk::mesh::ConstPartVector m_parts;
  stk::mesh::EntityVector m_entities;
  stk::mesh::Entity m_node;
  mutable bool m_active = true;
  int m_id = -1;
  std::vector<EntityOwnerProc> m_entityOwnerProcVec;
  BlockPair m_blockPair;
  bool m_hasBlockPair = false;
};

stk::mesh::EntityVector get_elements_for_node_in_parts(const stk::mesh::BulkData& bulk, const stk::mesh::Entity node, const stk::mesh::ConstPartVector& parts);
std::vector<int> get_elements_owners_for_node_in_parts(const stk::mesh::BulkData& bulk, const stk::mesh::Entity node,
                                                       const stk::mesh::ConstPartVector& parts);
std::vector<int> get_elements_owners(const stk::mesh::BulkData& bulk, const stk::mesh::EntityVector& entities);
}}}

#endif

