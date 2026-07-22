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
#include "stk_unit_test_utils/ioUtils.hpp"

namespace stk { namespace mesh { class Part; } }

namespace
{

int get_other_proc(int myProc)
{
  return (myProc == 1) ? 0 : 1;
}

void verify_that_elem1_and_node1_are_only_valid_on_p0(const stk::mesh::BulkData &bulkData, stk::mesh::Entity elem1, stk::mesh::Entity node1)
{
  if (bulkData.parallel_rank() == 0)
  {
    EXPECT_TRUE(bulkData.is_valid(elem1));
    EXPECT_TRUE(bulkData.is_valid(node1));
  }
  else {
    EXPECT_FALSE(bulkData.is_valid(elem1)); //elem1 only valid on proc 0, initially
    EXPECT_FALSE(bulkData.is_valid(node1)); //node1 only valid on proc 0, initially
  }
}

void verify_that_elem1_and_downward_connected_entities_are_ghosted_from_p0_to_p1(const stk::mesh::BulkData &bulkData, stk::mesh::EntityId id)
{
  //now we have ghosted elem1 from proc 0 to proc 1, so it should be valid on both procs
  //when an entity is ghosted, any downward-connected entities for that entity will also
  //be ghosted. So node1 should now also be valid on both procs
  stk::mesh::Entity elem1 = bulkData.get_entity(stk::topology::ELEM_RANK, id);
  stk::mesh::Entity node1 = bulkData.get_entity(stk::topology::NODE_RANK, id);
  EXPECT_TRUE(bulkData.is_valid(elem1));
  EXPECT_TRUE(bulkData.is_valid(node1));
}

void verify_elem1_is_valid_only_on_p0(const stk::mesh::BulkData &bulk, stk::mesh::Entity elem1)
{
  if(bulk.parallel_rank() == 0)
    EXPECT_TRUE(bulk.is_valid(elem1));
  else
    EXPECT_TRUE(!bulk.is_valid(elem1));
}

void verify_elem1_is_valid_on_both_procs(const stk::mesh::BulkData &bulk, stk::mesh::EntityId id)
{
  stk::mesh::Entity elem1 = bulk.get_entity(stk::topology::ELEM_RANK, id);
  EXPECT_TRUE(bulk.is_valid(elem1));
}

//BEGIN_GHOST_ELEM
TEST(StkMeshHowTo, customGhostElem)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(communicator) == 2)
  {
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(communicator).create();
    stk::mesh::BulkData& bulkData = *bulkPtr;
    stk::io::fill_mesh("generated:1x1x4", bulkData);

    stk::mesh::EntityId id = 1;
    stk::mesh::Entity elem1 = bulkData.get_entity(stk::topology::ELEM_RANK, id);
    stk::mesh::Entity node1 = bulkData.get_entity(stk::topology::NODE_RANK, id);
    verify_that_elem1_and_node1_are_only_valid_on_p0(bulkData, elem1, node1);

    bulkData.modification_begin();
    stk::mesh::Ghosting& ghosting = bulkData.create_ghosting("custom ghost for elem 1");
    std::vector<std::pair<stk::mesh::Entity, int> > elemProcPairs;
    if (bulkData.parallel_rank() == 0)
      elemProcPairs.push_back(std::make_pair(elem1, get_other_proc(bulkData.parallel_rank())));
    bulkData.change_ghosting(ghosting, elemProcPairs);
    bulkData.modification_end();

    verify_that_elem1_and_downward_connected_entities_are_ghosted_from_p0_to_p1(bulkData, id);
  }
}

TEST(StkMeshHowTo, addElementToGhostingUsingSpecializedModificationForPerformance)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(communicator) == 2)
  {
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(communicator).create();
    stk::mesh::BulkData& bulk = *bulkPtr;
    stk::io::fill_mesh("generated:1x1x4", bulk);

    stk::mesh::EntityId elementId = 1;
    stk::mesh::Entity elem1 = bulk.get_entity(stk::topology::ELEM_RANK, elementId);
    verify_elem1_is_valid_only_on_p0(bulk, elem1);

    bulk.modification_begin();
    stk::mesh::Ghosting& ghosting = bulk.create_ghosting("my custom ghosting");
    bulk.modification_end();

    stk::mesh::EntityProcVec entityProcPairs;
    if(bulk.parallel_rank() == 0)
      entityProcPairs.push_back(stk::mesh::EntityProc(elem1, get_other_proc(bulk.parallel_rank())));

    bulk.batch_add_to_ghosting(ghosting, entityProcPairs);

    verify_elem1_is_valid_on_both_procs(bulk, elementId);
  }
}
//END_GHOST_ELEM
}
