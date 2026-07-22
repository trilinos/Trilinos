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
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MeshBuilder.hpp>   // for MeshBuilder
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, entity_rank_names
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Types.hpp"      // for EntityId
#include "stk_mesh/base/Comm.hpp"      // for comm_mesh_counts
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stkMeshTestUtils.hpp"

namespace {

void verify_global_node_count(size_t expectedTotalNumNodes, const stk::mesh::BulkData &mesh)
{
  std::vector<size_t> entity_counts;
  stk::mesh::comm_mesh_counts(mesh, entity_counts);
  EXPECT_EQ(expectedTotalNumNodes, entity_counts[stk::topology::NODE_RANK]);
}

void verify_nodes_2_and_3_are_no_longer_shared(const stk::mesh::BulkData &mesh, const stk::mesh::EntityVector &nodes)
{
  EXPECT_TRUE(mesh.is_valid(nodes[0]));
  ASSERT_TRUE(mesh.is_valid(nodes[1]));
  ASSERT_TRUE(mesh.is_valid(nodes[2]));
  EXPECT_FALSE(mesh.bucket(nodes[1]).shared());
  EXPECT_FALSE(mesh.bucket(nodes[2]).shared());
}

void verify_nodes_2_and_3_are_removed(const stk::mesh::BulkData &mesh, const stk::mesh::EntityVector &nodes)
{
  EXPECT_TRUE(mesh.is_valid(nodes[0]));
  // These nodes were deleted because the special marking for "independent"
  // nodes was removed when the nodes became connected to the element and
  // now that the element is deleted, these nodes are no longer needed on
  // proc 1.
  EXPECT_FALSE(mesh.is_valid(nodes[1]));
  EXPECT_FALSE(mesh.is_valid(nodes[2]));
}

//BEGIN
TEST(stkMeshHowTo, createSharedNodes)
{
  const unsigned spatialDimension = 2;
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  builder.set_spatial_dimension(spatialDimension);
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
  stk::mesh::MetaData& metaData = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& bulkData = *bulkPtr;
  stk::mesh::Part &triPart = metaData.declare_part_with_topology("tri_part", stk::topology::TRIANGLE_3_2D);
  metaData.commit();

  if (bulkData.parallel_size() == 2)
  {
    bulkData.modification_begin();

    const unsigned nodesPerElem = 3;
    stk::mesh::EntityIdVector elemIds = {1, 2};//one elemId for each proc
    std::vector<stk::mesh::EntityIdVector> elemNodeIds = { {1, 3, 2}, {4, 2, 3} };
    const int myproc = bulkData.parallel_rank();

    stk::mesh::Entity elem = bulkData.declare_element(elemIds[myproc], stk::mesh::ConstPartVector{&triPart});
    stk::mesh::EntityVector elemNodes(nodesPerElem);
    elemNodes[0] = bulkData.declare_node(elemNodeIds[myproc][0]);
    elemNodes[1] = bulkData.declare_node(elemNodeIds[myproc][1]);
    elemNodes[2] = bulkData.declare_node(elemNodeIds[myproc][2]);

    bulkData.declare_relation(elem, elemNodes[0], 0);
    bulkData.declare_relation(elem, elemNodes[1], 1);
    bulkData.declare_relation(elem, elemNodes[2], 2);

    int otherproc = testUtils::get_other_proc(myproc);
    bulkData.add_node_sharing(elemNodes[1], otherproc);
    bulkData.add_node_sharing(elemNodes[2], otherproc);

    bulkData.modification_end();

    const size_t expectedTotalNumNodes = 4;
    verify_global_node_count(expectedTotalNumNodes, bulkData);
  }
}
//END

//BEGININDEP
TEST(stkMeshHowTo, createIndependentSharedNodes)
{
  const unsigned spatialDimension = 2;
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  builder.set_spatial_dimension(spatialDimension);
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
  stk::mesh::BulkData& bulkData = *bulkPtr;
  stk::mesh::MetaData& metaData = bulkPtr->mesh_meta_data();
  metaData.commit();

  if (bulkData.parallel_size() == 2)
  {
    bulkData.modification_begin();

    const unsigned nodesPerProc = 3;
    std::vector<stk::mesh::EntityIdVector> nodeIds = { {1, 3, 2}, {4, 2, 3} };
    const int myproc = bulkData.parallel_rank();
    stk::mesh::EntityVector nodes(nodesPerProc);
    nodes[0] = bulkData.declare_node(nodeIds[myproc][0]);
    nodes[1] = bulkData.declare_node(nodeIds[myproc][1]);
    nodes[2] = bulkData.declare_node(nodeIds[myproc][2]);

    int otherproc = testUtils::get_other_proc(myproc);
    bulkData.add_node_sharing(nodes[1], otherproc);
    bulkData.add_node_sharing(nodes[2], otherproc);

    bulkData.modification_end();

    const size_t expectedTotalNumNodes = 4;
    verify_global_node_count(expectedTotalNumNodes, bulkData);
  }
}
//ENDINDEP

//BEGIN_INDEP_DEP
TEST(stkMeshHowTo, createIndependentSharedNodesThenAddDependence)
{
  const unsigned spatialDimension = 2;
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  builder.set_spatial_dimension(spatialDimension);
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
  stk::mesh::BulkData& bulkData = *bulkPtr;
  stk::mesh::MetaData& metaData = bulkPtr->mesh_meta_data();
  stk::mesh::Part &triPart = metaData.declare_part_with_topology("triPart", stk::topology::TRIANGLE_3_2D);
  metaData.commit();

  if(bulkData.parallel_size() == 2)
  {
    bulkData.modification_begin();

    const unsigned nodesPerProc = 3;
    std::vector<stk::mesh::EntityIdVector> nodeIds = { {1, 3, 2}, {4, 2, 3}};
    const int myproc = bulkData.parallel_rank();

    stk::mesh::EntityVector nodes(nodesPerProc);
    nodes[0] = bulkData.declare_node(nodeIds[myproc][0]);
    nodes[1] = bulkData.declare_node(nodeIds[myproc][1]);
    nodes[2] = bulkData.declare_node(nodeIds[myproc][2]);

    int otherproc = testUtils::get_other_proc(myproc);
    bulkData.add_node_sharing(nodes[1], otherproc);
    bulkData.add_node_sharing(nodes[2], otherproc);

    const size_t expectedNumNodesPriorToModEnd = 6;
    verify_global_node_count(expectedNumNodesPriorToModEnd, bulkData);

    bulkData.modification_end();

    const size_t expectedNumNodesAfterModEnd = 4; // nodes 2 and 3 are shared
    verify_global_node_count(expectedNumNodesAfterModEnd, bulkData);

    const unsigned elemsPerProc = 1;
    stk::mesh::EntityId elemIds[][elemsPerProc] = { {1}, {2}};

    bulkData.modification_begin();
    stk::mesh::Entity elem = bulkData.declare_element(elemIds[myproc][0], stk::mesh::ConstPartVector{&triPart});
    bulkData.declare_relation(elem, nodes[0], 0);
    bulkData.declare_relation(elem, nodes[1], 1);
    bulkData.declare_relation(elem, nodes[2], 2);
    EXPECT_NO_THROW(bulkData.modification_end());

    bulkData.modification_begin();
    bulkData.destroy_entity(elem);
    bulkData.modification_end();

    if(myproc == 0)
      verify_nodes_2_and_3_are_no_longer_shared(bulkData, nodes);

    else  // myproc == 1
      verify_nodes_2_and_3_are_removed(bulkData, nodes);
  }
}

//END_INDEP_DEP
}
