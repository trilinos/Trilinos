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
#ifndef _DisconnectBlocksImpl_hpp_
#define _DisconnectBlocksImpl_hpp_

#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/SideSetEntry.hpp"
#include "stk_mesh/base/Entity.hpp"
#include <utility>
#include <vector>

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Part; } }

namespace stk {
namespace tools {
namespace impl {

using NodeMapType = std::map<stk::mesh::Entity, stk::mesh::Entity>;
using BlockPairType = std::pair<stk::mesh::Part*, stk::mesh::Part*>;
using SideSetType = std::vector<stk::mesh::SideSetEntry>;

bool is_block(const stk::mesh::BulkData & bulk, stk::mesh::Part & part);

int64_t get_block_id_for_element(const stk::mesh::BulkData & bulk, stk::mesh::Entity element);

void get_nodes_for_element_side(const stk::mesh::BulkData & bulk,
                                stk::mesh::Entity element,
                                stk::mesh::ConnectivityOrdinal sideOrdinal,
                                std::vector<stk::mesh::Entity> & sideNodes);

void get_node_ordinals_for_element_side(const stk::mesh::BulkData & bulk,
                                        stk::mesh::Entity element,
                                        stk::mesh::ConnectivityOrdinal sideOrdinal,
                                        std::vector<stk::mesh::ConnectivityOrdinal> & sideNodeOrdinals);

void disconnect_side_nodes(stk::mesh::BulkData & bulk,
                           const NodeMapType nodeMap,
                           stk::mesh::Entity elem,
                           const std::vector<stk::mesh::Entity> & sideNodes,
                           const std::vector<stk::mesh::ConnectivityOrdinal> & sideNodeOrdinals);

NodeMapType get_nodes_to_disconnect(const stk::mesh::BulkData & bulk,
                                    const BlockPairType & blockPair,
                                    const std::vector<stk::mesh::SideSetEntry> & sidesToDisconnect);

void create_new_duplicate_nodes(stk::mesh::BulkData & bulk, NodeMapType & nodeMap);

void get_all_blocks_in_mesh(const stk::mesh::BulkData & bulk, stk::mesh::PartVector & blocksInMesh);

std::vector<BlockPairType> get_block_pairs_to_disconnect(const stk::mesh::BulkData & bulk);

std::vector<SideSetType> get_sidesets_to_disconnect(stk::mesh::BulkData & bulk,
                                                    const std::vector<BlockPairType> & blockPairsToDisconnect);

void disconnect_block_pair_boundary(stk::mesh::BulkData & bulk,
                                    const BlockPairType & blockPair,
                                    const SideSetType & sideSet,
                                    const NodeMapType & nodeMap);

void disconnect_face_adjacent_elements(stk::mesh::BulkData & bulk,
                                       const BlockPairType & blockPair,
                                       const NodeMapType & nodeMap);

NodeMapType get_non_face_adjacent_nodes(stk::mesh::BulkData & bulk,
                                        const std::vector<BlockPairType> & blockPairsToDisconnect);

void disconnect_non_face_adjacent_elements(stk::mesh::BulkData & bulk,
                                           const BlockPairType & blockPair,
                                           const NodeMapType & nodesToDisconnect);

} } }

#endif
