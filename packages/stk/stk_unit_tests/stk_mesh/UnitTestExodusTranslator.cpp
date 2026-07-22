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
#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/ExodusTranslator.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <vector>
#include <string>
#include <memory>

namespace {

std::unique_ptr<stk::mesh::BulkData> setup_2block_shell_tri_mesh(MPI_Comm comm)
{
  std::unique_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(comm).set_spatial_dimension(3).create();

  std::string meshDesc = "0,1,SHELL_TRI_3, 1,2,3, BLOCK_1\n"
                         "0,2,SHELL_TRI_3, 2,3,4, BLOCK_2";
  std::vector<double> coords = {0,0,0,  0,1,0,  1,0,0,  1,1,0};
  stk::unit_test_util::setup_text_mesh(*bulkPtr, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coords));

  return bulkPtr;
}

TEST(TestExodusTranslator, fill_element_block_parts )
{
  MPI_Comm comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::unique_ptr<stk::mesh::BulkData> bulkPtr = setup_2block_shell_tri_mesh(comm);

  stk::mesh::PartVector allBlockParts;
  stk::mesh::fill_element_block_parts(bulkPtr->mesh_meta_data(), stk::topology::INVALID_TOPOLOGY, allBlockParts);
  ASSERT_EQ(2u, allBlockParts.size());

  stk::mesh::fill_element_block_parts(bulkPtr->mesh_meta_data(), stk::topology::SHELL_TRI_3, allBlockParts);
  ASSERT_EQ(2u, allBlockParts.size());

  stk::mesh::fill_element_block_parts(bulkPtr->mesh_meta_data(), stk::topology::HEX_8, allBlockParts);
  ASSERT_EQ(0u, allBlockParts.size());
}

TEST(TestExodusTranslator, get_element_block_part)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::unique_ptr<stk::mesh::BulkData> bulkPtr = setup_2block_shell_tri_mesh(comm);

  stk::mesh::Entity node2 = bulkPtr->get_entity(stk::topology::NODE_RANK, 2);
  EXPECT_ANY_THROW(stk::mesh::get_element_block_part(*bulkPtr, node2));

  stk::mesh::Entity elem2 = bulkPtr->get_entity(stk::topology::ELEM_RANK, 2);
  stk::mesh::Part* blockPart = stk::mesh::get_element_block_part(*bulkPtr, elem2);
  ASSERT_NE(nullptr, blockPart);
  EXPECT_EQ(blockPart->name(), "BLOCK_2");
}

} // anonymous namespace

