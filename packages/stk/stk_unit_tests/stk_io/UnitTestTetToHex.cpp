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

#include <gtest/gtest.h>
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/GetEntities.hpp>  // for get_selected_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_io/IossBridge.hpp>
#include <string>                       // for string
#include <vector>                       // for vector, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator&
#include "stk_mesh/base/Types.hpp"      // for PartVector, BucketVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine
#include "stk_mesh/base/FEMHelpers.hpp"

namespace {

void check_hex_mesh(const std::string& fileName)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(comm).create();
  stk::io::fill_mesh(fileName, *bulkPtr);
  stk::mesh::Part& hexPart = bulkPtr->mesh_meta_data().get_topology_root_part(stk::topology::HEX_8);
  stk::mesh::Part& tetPart = bulkPtr->mesh_meta_data().get_topology_root_part(stk::topology::TET_4);

  const unsigned numExpectedHexes = 1;
  const unsigned numExpectedTets = 0;

  EXPECT_EQ(numExpectedHexes, stk::mesh::count_entities(*bulkPtr, stk::topology::ELEM_RANK, hexPart));
  EXPECT_EQ(numExpectedTets, stk::mesh::count_entities(*bulkPtr, stk::topology::ELEM_RANK, tetPart));
}

TEST( StkMeshIoBroker, change_tet_to_hex_block )
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(comm).create();

  stk::io::fill_mesh("generated:1x1x1|tets", *bulkPtr);

  const unsigned numExpectedElems = 6;
  EXPECT_EQ(numExpectedElems, stk::mesh::count_entities(*bulkPtr, stk::topology::ELEM_RANK, bulkPtr->mesh_meta_data().universal_part()));

  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();

  stk::mesh::Part* block1Tet = meta.get_part("block_1");
  ASSERT_TRUE(block1Tet != nullptr);
  meta.rename(*block1Tet, "block_1_tet");
  stk::io::remove_io_part_attribute(*block1Tet);
  EXPECT_FALSE(stk::io::has_io_part_attribute(*block1Tet));
  meta.set_part_id(*block1Tet, stk::mesh::Part::INVALID_ID);

  stk::mesh::Part& block1Hex = meta.declare_part_with_topology("block_1", stk::topology::HEX_8);
  
  meta.set_part_id(block1Hex, 1);
  stk::io::put_io_part_attribute(block1Hex);

  bulkPtr->modification_begin();
  stk::mesh::PartVector hexParts = { &block1Hex };
  stk::mesh::EntityIdVector nodeIds = { 1, 2, 4, 3, 5, 6, 8, 7 };
  stk::mesh::EntityId elemId = 7;//still have 6 tet elements, ids=1..6
  stk::mesh::declare_element(*bulkPtr, hexParts, elemId, nodeIds);

  bulkPtr->destroy_elements_of_topology(stk::topology::TET_4);
  bulkPtr->modification_end();

  std::string fileName("hexmesh.g");
  stk::io::write_mesh(fileName, *bulkPtr);

  check_hex_mesh(fileName);
  unlink(fileName.c_str());
}

} // anonymous namespace

