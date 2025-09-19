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
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_tools/mesh_clone/ReplaceBulkData.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <gtest/gtest.h>

namespace
{

void verify_num_elements(const stk::mesh::BulkData& bulk, unsigned goldNumElements)
{
  EXPECT_EQ(goldNumElements, stk::mesh::count_entities(bulk, stk::topology::ELEM_RANK, bulk.mesh_meta_data().universal_part()));
}

TEST(ReplaceBulkData, generated)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD)
                                              .set_initial_bucket_capacity(1)
                                              .set_maximum_bucket_capacity(1).create();

  stk::io::fill_mesh("generated:1x1x2", *bulk);

  verify_num_elements(*bulk, 2);

  std::shared_ptr<stk::mesh::BulkData> replacementBulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD)
                                                         .set_initial_bucket_capacity(1)
                                                         .set_maximum_bucket_capacity(1).create();

  stk::io::fill_mesh("generated:1x1x1", *replacementBulk);

  stk::tools::replace_bulk_data(*replacementBulk, *bulk);

  verify_num_elements(*bulk, 1);
}

TEST(ReplaceBulkData, increment_sync_count)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD).create();

  stk::io::fill_mesh("generated:1x1x2", *bulk);

  verify_num_elements(*bulk, 2);

  std::shared_ptr<stk::mesh::BulkData> replacementBulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD).create();

  stk::io::fill_mesh("generated:1x1x1", *replacementBulk);

  unsigned syncCountBeforeReplace = bulk->synchronized_count();

  stk::tools::replace_bulk_data(*replacementBulk, *bulk);

  EXPECT_EQ(syncCountBeforeReplace+1, bulk->synchronized_count());
}

TEST(ReplaceBulkData, unchangedMetaDataParts)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD).create();
  bulk->mesh_meta_data().declare_part("newPart");

  stk::io::fill_mesh("generated:1x1x2", *bulk);

  verify_num_elements(*bulk, 2);

  std::shared_ptr<stk::mesh::BulkData> replacementBulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD).create();

  stk::io::fill_mesh("generated:1x1x1", *replacementBulk);

  stk::tools::replace_bulk_data(*replacementBulk, *bulk);

  verify_num_elements(*bulk, 1);

  stk::mesh::Part* newPart = bulk->mesh_meta_data().get_part("newPart");
  EXPECT_TRUE(newPart != nullptr);
}

}
