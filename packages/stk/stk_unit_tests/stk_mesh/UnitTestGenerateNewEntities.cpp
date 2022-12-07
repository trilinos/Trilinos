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

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <gtest/gtest.h>                             // for AssertHelper, etc
#include <stddef.h>                                  // for size_t
#include <stk_mesh/base/BulkData.hpp>                // for BulkData
#include <stk_mesh/base/MetaData.hpp>                // for MetaData
#include <stk_util/parallel/Parallel.hpp>            // for ParallelMachine
#include <vector>                                    // for vector, etc
#include "mpi.h"                                     // for MPI_COMM_WORLD, etc
#include "stk_mesh/base/Types.hpp"                   // for EntityVector, etc
#include <stk_unit_test_utils/BuildMesh.hpp>
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace {
using stk::unit_test_util::build_mesh;

TEST( UnitTestStkMeshGenerateNewEntities , testUnit )
{
  // Test BulkData's generate_new_entities method.

  stk::ParallelMachine pm(MPI_COMM_WORLD);

  const int spatial_dimension = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatial_dimension, pm);
  stk::mesh::MetaData& meta_data = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& bulk_data = *bulkPtr;

  meta_data.commit();

  const stk::mesh::PartVector no_parts;

  bulk_data.modification_begin();

  bulk_data.declare_node(bulk_data.parallel_rank() + 1, no_parts);

  bulk_data.modification_end();

  // Create a request vector for 2 new nodes on each processor
  size_t num_nodes_requested = 2;
  std::vector<size_t> requests(meta_data.entity_rank_count(), 0);
  requests[0] = num_nodes_requested;

  bulk_data.modification_begin();

  // generate_new_entities creates new blank entities of the requested ranks
  stk::mesh::EntityVector new_nodes;
  bulk_data.generate_new_entities(requests, new_nodes);
  ASSERT_EQ(new_nodes.size(), num_nodes_requested);

  // confirm that the nodes we created earlier are not in the new entities
  for (stk::mesh::EntityVector::const_iterator itr = new_nodes.begin();
       itr != new_nodes.end(); ++itr) {
    ASSERT_GT(static_cast<int>(bulk_data.identifier(*itr)), bulk_data.parallel_size());
  }

  bulk_data.modification_end();
}

}
