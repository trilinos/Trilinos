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
#ifndef STK_UNIT_TESTS_STK_IO_IOMeshFixture_hpp
#define STK_UNIT_TESTS_STK_IO_IOMeshFixture_hpp

#include <gtest/gtest.h>  // for AssertHelper, EXPECT_EQ, etc

#include <algorithm>
#include <stk_io/FillMesh.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>  // for StkMeshIoBroker
#include <stk_mesh/base/ExodusTranslator.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_topology/topology.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <string>
#include <vector>

namespace stk
{
namespace io
{
namespace unit_test
{

class IOMeshFixture : public stk::unit_test_util::MeshFixture
{
protected:
  stk::mesh::Part& create_io_part(const std::string& partName,
                                  int id = -1,
                                  stk::topology topo = stk::topology::HEX_8)
  {
    stk::mesh::Part& part = get_meta().declare_part_with_topology(partName, topo);
    if (id != -1) {
      get_meta().set_part_id(part, id);
    }
    stk::io::put_io_part_attribute(part);
    return part;
  }

  void move_element(const stk::mesh::EntityId elemId,
                    stk::mesh::Part& sourcePart,
                    stk::mesh::Part& destPart)
  {
    stk::mesh::Entity elem = get_bulk().get_entity(stk::topology::ELEM_RANK, elemId);
    get_bulk().batch_change_entity_parts({elem}, {&destPart}, {&sourcePart});
  }

  void move_face(const stk::mesh::EntityId faceId, stk::mesh::Part& sourcePart, stk::mesh::Part& destPart)
  {
    stk::mesh::Entity face = get_bulk().get_entity(get_meta().side_rank(), faceId);
    get_bulk().batch_change_entity_parts({face}, {&destPart}, {&sourcePart});
  }

  void add_nodes(const stk::mesh::EntityIdVector nodeIds, stk::mesh::Part& destPart)
  {
    stk::mesh::EntityVector nodes;
    nodes.reserve(nodeIds.size());

    for (stk::mesh::EntityId nodeId : nodeIds) {
      stk::mesh::Entity node = get_bulk().get_entity(stk::topology::NODE_RANK, nodeId);
      STK_ThrowRequire(get_bulk().is_valid(node));
      nodes.push_back(node);
    }

    get_bulk().batch_change_entity_parts(nodes, {&destPart}, {});
  }

  void add_node(const stk::mesh::EntityId nodeId, stk::mesh::Part& destPart) { add_nodes({nodeId}, destPart); }

  void create_side(const stk::mesh::EntityId elemId,
                   stk::mesh::ConnectivityOrdinal sideOrd,
                   stk::mesh::Part& sidePart)
  {
    stk::mesh::Entity elem = get_bulk().get_entity(stk::topology::ELEM_RANK, elemId);
    ASSERT_TRUE(get_bulk().is_valid(elem));
    get_bulk().modification_begin();
    get_bulk().declare_element_side(elem, sideOrd, stk::mesh::PartVector{&sidePart});
    get_bulk().modification_end();
  }

  stk::mesh::Selector create_block_subset_selector(const stk::mesh::PartVector& blocksToExclude)
  {
    stk::mesh::Selector meshSubsetSelector = get_meta().universal_part();
    if (!blocksToExclude.empty()) {
      stk::mesh::PartVector elemBlocks;
      stk::mesh::fill_element_block_parts(get_meta(), stk::topology::INVALID_TOPOLOGY, elemBlocks);
      for(const stk::mesh::Part* excludedBlock : blocksToExclude) {
        auto foundBlock = std::find(elemBlocks.begin(), elemBlocks.end(), excludedBlock);
        STK_ThrowRequire(foundBlock != elemBlocks.end());
        elemBlocks.erase(foundBlock);
      }
      meshSubsetSelector = stk::mesh::selectUnion(elemBlocks);
    }

    return meshSubsetSelector;
  }
};

}  // namespace unit_test
}  // namespace io
}  // namespace stk
#endif

