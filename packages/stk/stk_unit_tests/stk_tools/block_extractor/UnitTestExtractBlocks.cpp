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
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/BuildMesh.hpp>
#include <stk_tools/block_extractor/ExtractBlocks.hpp>
#include "stk_mesh/base/Entity.hpp"

namespace
{
using stk::unit_test_util::build_mesh;

class MeshWithTwoBlocks : public stk::unit_test_util::MeshFixture
{
protected:
  void switch_half_mesh_to_part(stk::mesh::Part &addPart, stk::mesh::Part &removePart)
  {
    get_bulk().modification_begin();
    switch_elem_to_part(1, addPart, removePart);
    switch_elem_to_part(2, addPart, removePart);
    get_bulk().modification_end();
  }

  void switch_elem_to_part(stk::mesh::EntityId elemId, stk::mesh::Part &addPart, stk::mesh::Part &removePart)
  {
    stk::mesh::Entity elem = get_bulk().get_entity(stk::topology::ELEM_RANK, elemId);
    if(get_bulk().is_valid(elem) && get_bulk().bucket(elem).owned())
      get_bulk().change_entity_parts(elem, stk::mesh::ConstPartVector{&addPart}, stk::mesh::ConstPartVector{&removePart});
  }

  void expect_num_elems_in_part(stk::mesh::BulkData &bulk, size_t expectedNum, stk::mesh::Part& block)
  {
    stk::mesh::Selector selector(block);
    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts(bulk, counts, &selector);
    EXPECT_EQ(expectedNum, counts[stk::topology::ELEM_RANK]);
  }
};

TEST_F(MeshWithTwoBlocks, extractBlock2_onlyHaveBlock2)
{
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  stk::mesh::Part &block2 = get_meta().declare_part_with_topology("block_2", stk::topology::HEX_8);
  stk::io::fill_mesh("generated:1x1x4", get_bulk());

  stk::mesh::Part &block1 = *get_meta().get_part("block_1");
  switch_half_mesh_to_part(block2, block1);

  expect_num_elems_in_part(get_bulk(), 2, block1);
  expect_num_elems_in_part(get_bulk(), 2, block2);

  std::shared_ptr<stk::mesh::BulkData> newBulkPtr = build_mesh(get_comm());
  stk::mesh::BulkData& newBulk = *newBulkPtr;
  stk::tools::extract_blocks(get_bulk(), newBulk, {"block_2"});

  expect_num_elems_in_part(newBulk, 0, block1);
  expect_num_elems_in_part(newBulk, 2, block2);
}

TEST_F(MeshWithTwoBlocks, getNodesetPartsAndNames)
{

  setup_mesh("generated:1x1x4|nodeset:xX", stk::mesh::BulkData::AUTO_AURA);


  std::vector<std::string> ns_names = stk::tools::find_nodeset_names_from_id(get_bulk(), {2,1});

  std::vector<std::string> ns_names_expected = {"nodelist_2","nodelist_1"};

  ASSERT_EQ(ns_names.size(), ns_names_expected.size() );

  for(size_t k = 0 ; k < ns_names_expected.size(); ++k)
    EXPECT_EQ(ns_names[k], ns_names_expected[k] );

  std::vector<stk::mesh::Part*> parts;
  stk::tools::GetPartsByName( parts,get_bulk(),ns_names);
  ASSERT_EQ(static_cast<size_t>(2), (parts.size()) );
  stk::mesh::Part &ns1 = *get_meta().get_part("nodelist_1");
  EXPECT_EQ(ns1, *(parts[1]) );

  stk::mesh::Part &ns2 = *get_meta().get_part("nodelist_2");
  EXPECT_EQ(ns2, *(parts[0]) );

}

TEST_F(MeshWithTwoBlocks, getBlockPartsAndNames)
{
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  stk::mesh::Part &block2 = get_meta().declare_part_with_topology("block_2", stk::topology::HEX_8);
  stk::io::fill_mesh("generated:1x1x4|nodeset:xX", get_bulk());

  stk::mesh::Part &block1 = *get_meta().get_part("block_1");
  switch_half_mesh_to_part(block2, block1);

  std::vector<std::string> blk_names = stk::tools::GetBlockNamesFromIDs(get_bulk(), {1,-1});


  std::vector<std::string> blk_names_expected = {"block_1","block_2"};

  ASSERT_EQ(blk_names.size(), blk_names_expected.size() );

  for(size_t k = 0 ; k < blk_names_expected.size(); ++k)
    EXPECT_EQ(blk_names[k], blk_names_expected[k] );

  std::vector<stk::mesh::Part*> parts;
  stk::tools::GetPartsByName( parts,get_bulk(),blk_names);
  ASSERT_EQ(static_cast<size_t>(2), (parts.size()) );
  stk::mesh::Part &blk1 = *get_meta().get_part("block_1");
  EXPECT_EQ(blk1, *(parts[0]) );

  stk::mesh::Part &blk2 = *get_meta().get_part("block_2");
  EXPECT_EQ(blk2, *(parts[1]) );

}

TEST_F(MeshWithTwoBlocks, getOneBlockAndOneNodeset)
{
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  stk::mesh::Part &block2 = get_meta().declare_part_with_topology("block_2", stk::topology::HEX_8);
  stk::io::fill_mesh("generated:1x1x4|nodeset:xX", get_bulk());

  stk::mesh::Part &block1 = *get_meta().get_part("block_1");
  switch_half_mesh_to_part(block2, block1);

  // apparently that block two gets constructed with ID=-1...
  std::vector<std::string> blk_names = stk::tools::GetBlockNamesFromIDs(get_bulk(), {-1});
  std::vector<std::string> ns_names = stk::tools::find_nodeset_names_from_id(get_bulk(), {1});

  stk::mesh::Selector theSelector = stk::tools::GetBlockAndNodesetSelector(get_bulk(), ns_names, blk_names);

  std::vector<size_t> counts;
  stk::mesh::comm_mesh_counts(get_bulk(), counts, &theSelector);
  size_t expectedNumElem = 2; // two elements in block 2
  size_t expectedNumNode = 16; // 12 nodes associated with 2 elems in block 2, then 4 nodes from nodeset 1

  EXPECT_EQ(expectedNumElem, counts[stk::topology::ELEM_RANK]);
  EXPECT_EQ(expectedNumNode, counts[stk::topology::NODE_RANK]);

}

class MeshWithOneBlock : public stk::unit_test_util::MeshFixture
{

};

TEST_F(MeshWithOneBlock, extractBlock2_throws)
{
  setup_mesh("generated:1x1x4", stk::mesh::BulkData::AUTO_AURA);
  std::shared_ptr<stk::mesh::BulkData> newBulkPtr = build_mesh(get_comm());
  EXPECT_THROW(stk::tools::extract_blocks(get_bulk(), *newBulkPtr, {"block_2"}), std::exception);
}

}
