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
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/BulkData.hpp>   // for MetaData

//BEGIN_CREATED
TEST(stkMeshHowTo, checkCreatedStateAfterMeshCreation)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(communicator);
  if (parallel_size != 1) { return; }
  const std::string fileName = "generated:1x1x1";
  stk::io::StkMeshIoBroker meshReader(communicator);
  meshReader.add_mesh_database(fileName, stk::io::READ_MESH);
  meshReader.create_input_mesh();
  meshReader.populate_bulk_data();
  stk::mesh::BulkData & bulk = meshReader.bulk_data();
  for (size_t n=1 ; n<=8 ; ++n) {
    stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, n);
    EXPECT_EQ( stk::mesh::Created, bulk.state(node) );
  }
  stk::mesh::Entity element = bulk.get_entity(stk::topology::ELEMENT_RANK,1);
  EXPECT_EQ( stk::mesh::Created, bulk.state(element) );
  bulk.modification_begin();
  bulk.modification_end();
  EXPECT_EQ( stk::mesh::Unchanged, bulk.state(element) );
}
//END_CREATED

//BEGIN_DELETED
TEST(stkMeshHowTo, checkDeletedState)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(communicator);
  if (parallel_size != 1) { return; }
  const std::string fileName = "generated:1x1x1";
  stk::io::StkMeshIoBroker meshReader(communicator);
  meshReader.add_mesh_database(fileName, stk::io::READ_MESH);
  meshReader.create_input_mesh();
  meshReader.populate_bulk_data();
  stk::mesh::BulkData & bulk = meshReader.bulk_data();
  stk::mesh::Entity element = bulk.get_entity(stk::topology::ELEMENT_RANK, 1);
  stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
  bulk.modification_begin();
  bulk.destroy_relation(element, node1, 0);
  EXPECT_TRUE(bulk.destroy_entity(node1));
  // Elements must have all their nodes defined:
  stk::mesh::Entity new_node = bulk.declare_node(25);
  bulk.declare_relation(element,new_node, 0);
  bulk.modification_end();
  EXPECT_FALSE( bulk.is_valid(node1) ); // deleted makes it not valid
  EXPECT_TRUE( bulk.in_index_range(node1));
  EXPECT_EQ( stk::mesh::Deleted, bulk.state(node1) );
  bulk.modification_begin();
  bulk.modification_end();
  EXPECT_FALSE( bulk.is_valid(node1) );
}
//END_DELETED

//BEGIN_MODIFIED
TEST(stkMeshHowTo, checkModifiedState)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(communicator);
  if (parallel_size != 1) { return; }
  const std::string fileName = "generated:1x1x1";
  stk::io::StkMeshIoBroker meshReader(communicator);
  meshReader.add_mesh_database(fileName, stk::io::READ_MESH);
  meshReader.create_input_mesh();
  meshReader.populate_bulk_data();
  stk::mesh::BulkData & bulk = meshReader.bulk_data();
  stk::mesh::Entity element = bulk.get_entity(stk::topology::ELEMENT_RANK, 1);
  stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
  bulk.modification_begin();
  bulk.destroy_relation(element, node1, 0);
  // Elements must have all their nodes defined:
  stk::mesh::Entity new_node = bulk.declare_node(25);
  bulk.declare_relation(element, new_node, 0);
  bulk.modification_end();
  EXPECT_EQ( stk::mesh::Modified, bulk.state(node1) );
  EXPECT_EQ( stk::mesh::Modified, bulk.state(element) );
}
//END_MODIFIED


//BEGIN_AURA_MODIFIED
TEST(stkMeshHowTo, checkAuraCreatedState)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(communicator);
  if (parallel_size != 2) { GTEST_SKIP(); }
  const std::string fileName = "generated:1x1x2";
  stk::io::StkMeshIoBroker meshReader(communicator);
  meshReader.add_mesh_database(fileName, stk::io::READ_MESH);
  meshReader.create_input_mesh();
  meshReader.populate_bulk_data();
  stk::mesh::BulkData & bulk = meshReader.bulk_data();
  stk::mesh::Entity element1 = bulk.get_entity(stk::topology::ELEMENT_RANK, 1);
  stk::mesh::Entity element2 = bulk.get_entity(stk::topology::ELEMENT_RANK, 2);
  if (bulk.parallel_rank() == 0) {
    EXPECT_EQ( stk::mesh::Created,  bulk.state(element1) );
    EXPECT_TRUE( bulk.in_receive_ghost(bulk.entity_key(element2)) );
    EXPECT_EQ( stk::mesh::Created, bulk.state(element2) );
  } else { // parallel_rank == 1
    EXPECT_TRUE( bulk.in_receive_ghost(bulk.entity_key(element1)) );
    EXPECT_EQ( stk::mesh::Created, bulk.state(element1) );
    EXPECT_EQ( stk::mesh::Created,  bulk.state(element2) );
  }
}
//END_AURA_MODIFIED

//BEGIN_CEO_MODIFIED
TEST(stkMeshHowTo, checkCEOModifiedState)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(communicator);
  if (parallel_size != 2) { GTEST_SKIP(); }
  const std::string fileName = "generated:1x1x2";
  stk::io::StkMeshIoBroker meshReader(communicator);
  meshReader.add_mesh_database(fileName, stk::io::READ_MESH);
  meshReader.create_input_mesh();
  meshReader.populate_bulk_data();
  stk::mesh::BulkData & bulk = meshReader.bulk_data();
  stk::mesh::Entity node5 = bulk.get_entity(stk::topology::NODE_RANK, 5);
  std::vector<stk::mesh::EntityProc> ceo_vector;
  if (bulk.parallel_rank() == 0) {
    ceo_vector.push_back(stk::mesh::EntityProc(node5,1));
  }
  bulk.change_entity_owner(ceo_vector);
  stk::mesh::Entity element1 = bulk.get_entity(stk::topology::ELEMENT_RANK, 1);
  stk::mesh::Entity element2 = bulk.get_entity(stk::topology::ELEMENT_RANK, 2);
  if (bulk.parallel_rank() == 0) {
    EXPECT_EQ( stk::mesh::Modified,  bulk.state(element1) );
    EXPECT_EQ( stk::mesh::Created, bulk.state(element2) );
  }
  else {
    EXPECT_EQ( stk::mesh::Created,  bulk.state(element1) );
    EXPECT_EQ( stk::mesh::Modified, bulk.state(element2) );
  }
}
//END_CEO_MODIFIED

//BEGIN_PARALLEL_MODIFIED
TEST(stkMeshHowTo, checkParallelConsistencyModifiedState)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(communicator);
  if (parallel_size != 2) { GTEST_SKIP(); }
  const std::string fileName = "generated:1x1x2";
  stk::io::StkMeshIoBroker meshReader(communicator);
  meshReader.add_mesh_database(fileName, stk::io::READ_MESH);
  meshReader.create_input_mesh();
  meshReader.populate_bulk_data();
  stk::mesh::BulkData & bulk = meshReader.bulk_data();
  stk::mesh::Entity element1 = bulk.get_entity(stk::topology::ELEMENT_RANK, 1);
  stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
  bulk.modification_begin();
  bulk.modification_end();
  if (bulk.parallel_rank() == 0) {
    EXPECT_EQ( stk::mesh::Unchanged, bulk.state(element1) );
  } else { // parallel_rank == 1
    EXPECT_EQ( stk::mesh::Unchanged, bulk.state(element1) );
  }
  bulk.modification_begin();
  if (bulk.parallel_rank() == 0) {
    bulk.destroy_relation(element1, node1, 0);
    stk::mesh::Entity new_node = bulk.declare_node(15);
    bulk.declare_relation(element1,new_node,0);
  }
  if (bulk.parallel_rank() == 0) {
    EXPECT_EQ( stk::mesh::Modified, bulk.state(element1) );
  } else { // parallel_rank == 1
    EXPECT_EQ( stk::mesh::Unchanged, bulk.state(element1) );
  }
  bulk.modification_end();
  if (bulk.parallel_rank() == 0) {
    EXPECT_EQ( stk::mesh::Modified, bulk.state(element1) );
  } else { // parallel_rank == 1
    EXPECT_EQ( stk::mesh::Created, bulk.state(element1) );
  }
}
//END_PARALLEL_MODIFIED

