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

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t
#include <sstream>                      // for ostringstream, etc
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_topology/topology.hpp>    // for topology, etc
#include <string>                       // for string
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH
#include "stk_mesh/base/Types.hpp"      // for BucketVector, PartVector
#include <stk_unit_test_utils/ioUtils.hpp>  // for fill_mesh_using_stk_io
#include "stkMeshTestUtils.hpp"

namespace stk { namespace mesh { class Part; } }

namespace
{

bool entities_have_part(const stk::mesh::BulkData& bulk, stk::mesh::EntityVector entities, stk::mesh::Part* part)
{
  for(stk::mesh::Entity entity : entities) {
    stk::mesh::Bucket& bucket = bulk.bucket(entity);
    if(!bucket.member(*part)) {
      return false;
    }
  }
  return true;
}

//BEGIN
TEST(StkMeshHowTo, changeEntityPartsUsingSelector)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(communicator) > 1) { return;}

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(communicator).create();
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();

  unsigned elementCount = 10;
  stk::io::fill_mesh("generated:1x1x" + std::to_string(elementCount), *bulkPtr);

  stk::mesh::Part* block1Part = meta.get_part("block_1");
  stk::mesh::Part& newBlock2Part = meta.declare_part("block_2", stk::topology::ELEM_RANK);
  stk::mesh::Selector block1Selector(*block1Part);
  stk::mesh::EntityVector entities;

  stk::mesh::get_entities(*bulkPtr, stk::topology::ELEM_RANK,  block1Selector, entities);
  EXPECT_EQ(elementCount, entities.size());
  EXPECT_TRUE(entities_have_part(*bulkPtr, entities, block1Part));

  stk::mesh::PartVector addParts{&newBlock2Part};
  stk::mesh::PartVector removeParts{block1Part};
  bulkPtr->batch_change_entity_parts(stk::mesh::Selector(*block1Part), stk::topology::ELEM_RANK,
                                 addParts, removeParts);

  EXPECT_TRUE(entities_have_part(*bulkPtr, entities, &newBlock2Part));
  EXPECT_FALSE(entities_have_part(*bulkPtr, entities, block1Part));
}
//END
}
