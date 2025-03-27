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
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_unit_test_utils/stk_mesh_fixtures/HexFixture.hpp>
#include <memory>
#include <vector>

std::shared_ptr<stk::mesh::BulkData>
setup_hex_mesh(MPI_Comm comm, size_t nx, size_t ny, size_t nz)
{
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(comm)
                                         .set_spatial_dimension(3)
                                         .set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA)
                                         .create();
  stk::mesh::fixtures::HexFixture::fill_mesh(nx,ny,nz, *bulkPtr);
  return bulkPtr;
}

TEST(GhostingSymmetry, ghostNode_fromProc0_toProc1andProc2)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  const int numProcs = stk::parallel_machine_size(comm);
  if(numProcs != 3) { GTEST_SKIP(); }
  const int myProc = stk::parallel_machine_rank(comm);

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = setup_hex_mesh(comm, 1, 1, 3);

  stk::mesh::Entity node1 = bulkPtr->get_entity(stk::topology::NODE_RANK, 1);

  //initially node1 only exists on proc 0
  if (myProc == 0) {
    EXPECT_TRUE(bulkPtr->is_valid(node1));
  }
  else {
    EXPECT_FALSE(bulkPtr->is_valid(node1));
  }

  bulkPtr->modification_begin();

  stk::mesh::Ghosting& myGhosting = bulkPtr->create_ghosting("myCustomGhosting");
  std::vector<stk::mesh::EntityProc> nodesToGhost;
  
  if (myProc == 0) {
    //proc 0 owns node1, ghost it to procs 1 and 2:
    nodesToGhost.push_back(stk::mesh::EntityProc(node1, 1));
    nodesToGhost.push_back(stk::mesh::EntityProc(node1, 2));
  }

  bulkPtr->change_ghosting(myGhosting, nodesToGhost);

  bulkPtr->modification_end();

  node1 = bulkPtr->get_entity(stk::topology::NODE_RANK, 1);
  EXPECT_TRUE(bulkPtr->is_valid(node1));

  std::vector<int> commProcs;
  bulkPtr->comm_procs(node1, commProcs);

  const bool isOwned = bulkPtr->bucket(node1).owned();
  size_t numOtherProcs = isOwned ? 2 : 1;
  //because ghost receivers know about owner, not each other

  EXPECT_EQ(numOtherProcs, commProcs.size());
}

