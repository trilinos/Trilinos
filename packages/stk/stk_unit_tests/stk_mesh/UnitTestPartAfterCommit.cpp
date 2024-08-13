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

#include <gtest/gtest.h>                // for AssertHelper, ASSERT_EQ, etc
#include <stddef.h>                     // for size_t
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/GetEntities.hpp>  // for get_selected_entities, etc
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <string>                       // for string
#include "mpi.h"                        // for MPI_COMM_WORLD, etc
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH
#include "stk_io/StkMeshIoBroker.hpp"   // for StkMeshIoBroker
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity, operator<<
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, put_field
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator-
#include "stk_mesh/base/Types.hpp"      // for EntityVector, PartVector
#include "stk_topology/topology.hpp"    // for topology, etc

void testNodesAreSelected(stk::mesh::BulkData &stkMeshBulkData,
                          const stk::mesh::EntityVector &nodes,
                          stk::mesh::Selector part1Selector)
{
  stk::mesh::EntityVector selectedNodes;
  stk::mesh::get_selected_entities(part1Selector, stkMeshBulkData.buckets(stk::topology::NODE_RANK), selectedNodes);
  ASSERT_EQ(selectedNodes.size(), nodes.size());
  for(size_t i = 0; i < nodes.size(); i++)
  {
    EXPECT_EQ(nodes[i], selectedNodes[i]);
  }
}

void testFieldOnNodes(stk::mesh::BulkData &stkMeshBulkData,
                      const stk::mesh::EntityVector &nodes,
                      stk::mesh::Part &nodePart,
                      stk::mesh::Field<double> &nodeField1)
{
  for(size_t i = 0; i < nodes.size(); ++i)
  {
    EXPECT_TRUE(stkMeshBulkData.bucket(nodes[i]).member(nodePart));
    double* fieldPtr = stk::mesh::field_data(nodeField1, nodes[i]);
    EXPECT_EQ(*static_cast<double*>(nodeField1.get_initial_value()), *fieldPtr);
  }
}

void testAddingNodesToPart(stk::mesh::BulkData &stkMeshBulkData,
                           const stk::mesh::EntityVector &nodes,
                           stk::mesh::Part &nodePart,
                           stk::mesh::Selector part1Selector,
                           stk::mesh::Field<double> &nodeField1)
{
  stkMeshBulkData.modification_begin();
  stk::mesh::PartVector addParts(1, &nodePart);
  for(size_t i = 0; i < nodes.size(); ++i)
  {
    if(stkMeshBulkData.parallel_owner_rank(nodes[i]) == stkMeshBulkData.parallel_rank())
    {
      stkMeshBulkData.change_entity_parts(nodes[i], addParts);
    }
  }
  stkMeshBulkData.modification_end();

  testFieldOnNodes(stkMeshBulkData, nodes, nodePart, nodeField1);

  testNodesAreSelected(stkMeshBulkData, nodes, part1Selector);
}

void testRemovingNodesFromPart(stk::mesh::BulkData &stkMeshBulkData,
                               const stk::mesh::EntityVector &nodes,
                               stk::mesh::Part &nodePart,
                               stk::mesh::Selector partSelector,
                               stk::mesh::Field<double> &nodeField1)
{
  stkMeshBulkData.modification_begin();
  stk::mesh::PartVector emptyAddParts;
  stk::mesh::PartVector removeParts(1, &nodePart);
  for(size_t i = 0; i < nodes.size(); ++i)
  {
    if(stkMeshBulkData.parallel_owner_rank(nodes[i]) == stkMeshBulkData.parallel_rank())
    {
      stkMeshBulkData.change_entity_parts(nodes[i], emptyAddParts, removeParts);
    }
  }
  stkMeshBulkData.modification_end();

  for(size_t i = 0; i < nodes.size(); ++i)
  {
    EXPECT_FALSE(stkMeshBulkData.bucket(nodes[i]).member(nodePart));
  }

  stk::mesh::EntityVector selectedNodes;
  stk::mesh::get_selected_entities(partSelector, stkMeshBulkData.buckets(stk::topology::NODE_RANK), selectedNodes);
  EXPECT_EQ(0u, selectedNodes.size());
}

TEST(UnitTestPartsAfterCommit, FieldsAndSelectors)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;

  stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
  const std::string generatedMeshSpecification = "generated:1x1x4";
  stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
  stkMeshIoBroker.create_input_mesh();

  stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
  stk::mesh::Part& nodePart1 = stkMeshMetaData.declare_part("nodePart1");
  stk::mesh::Field<double>& nodeField1 = stkMeshMetaData.declare_field<double>(stk::topology::NODE_RANK, "nodeField1");
  const double initialValue = 3.14;
  stk::mesh::put_field_on_mesh(nodeField1, nodePart1, &initialValue);
  stkMeshIoBroker.populate_bulk_data();

  stk::mesh::BulkData& stkMeshBulkData = stkMeshIoBroker.bulk_data();

  stk::mesh::EntityVector nodes;
  stk::mesh::get_entities(stkMeshBulkData, stk::topology::NODE_RANK, nodes);

  stk::mesh::Selector part1Selector = nodePart1;
  testAddingNodesToPart(stkMeshBulkData, nodes, nodePart1, part1Selector, nodeField1);

  EXPECT_TRUE(stkMeshMetaData.is_commit());

  stk::mesh::Part& partAfterCommit = stkMeshMetaData.declare_part("new_part");

  testAddingNodesToPart(stkMeshBulkData, nodes, partAfterCommit, part1Selector, nodeField1);

  stk::mesh::Selector newPartSelector = partAfterCommit;
  testNodesAreSelected(stkMeshBulkData, nodes, newPartSelector);

  testRemovingNodesFromPart(stkMeshBulkData, nodes, partAfterCommit, newPartSelector, nodeField1);
}

TEST(UnitTestPartsAfterCommit, PartInduction)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;

  stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
  const std::string generatedMeshSpecification = "generated:1x1x4";
  stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
  stkMeshIoBroker.create_input_mesh();

  stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
  stk::mesh::Part& firstPart = stkMeshMetaData.declare_part("firstPart", stk::topology::ELEMENT_RANK);

  stk::mesh::Field<double>& nodeField1 = stkMeshMetaData.declare_field<double>(stk::topology::NODE_RANK, "nodeField1");
  const double initialValue = 3.14;
  stk::mesh::put_field_on_mesh(nodeField1, firstPart, &initialValue);
  stkMeshIoBroker.populate_bulk_data();

  stk::mesh::BulkData& stkMeshBulkData = stkMeshIoBroker.bulk_data();

  stk::mesh::EntityVector locallyOwnedElements;
  stk::mesh::get_selected_entities(stkMeshMetaData.locally_owned_part(), stkMeshBulkData.buckets(stk::topology::ELEMENT_RANK), locallyOwnedElements);
  ASSERT_TRUE(!locallyOwnedElements.empty());

  stkMeshBulkData.modification_begin();
  stk::mesh::PartVector addParts(1, &firstPart);
  stkMeshBulkData.change_entity_parts(locallyOwnedElements[0], addParts);
  stkMeshBulkData.modification_end();

  stk::mesh::EntityVector nodesInFirstPart;
  stk::mesh::get_selected_entities(firstPart, stkMeshBulkData.buckets(stk::topology::NODE_RANK), nodesInFirstPart);

  for(size_t i=0; i<nodesInFirstPart.size(); i++)
  {
    double* fieldPtr = stk::mesh::field_data(nodeField1, nodesInFirstPart[i]);
    *fieldPtr = stkMeshBulkData.identifier(nodesInFirstPart[i]);
  }

  stk::mesh::Part& partAfterCommit = stkMeshMetaData.declare_part("partAfterCommit", stk::topology::ELEMENT_RANK);

  stkMeshBulkData.modification_begin();
  addParts[0] = &partAfterCommit;
  stkMeshBulkData.change_entity_parts(locallyOwnedElements[0], addParts);
  stkMeshBulkData.modification_end();

  stk::mesh::EntityVector nodesInPartDeclaredAfterCommit;
  stk::mesh::get_selected_entities(partAfterCommit, stkMeshBulkData.buckets(stk::topology::NODE_RANK), nodesInPartDeclaredAfterCommit);
  ASSERT_EQ(nodesInFirstPart.size(), nodesInPartDeclaredAfterCommit.size());
  for(size_t i=0; i<nodesInFirstPart.size(); i++)
  {
    EXPECT_EQ(nodesInFirstPart[i], nodesInPartDeclaredAfterCommit[i]);
  }

  for(size_t i = 0; i < nodesInFirstPart.size(); ++i)
  {
    double* fieldPtr = stk::mesh::field_data(nodeField1, nodesInFirstPart[i]);
    EXPECT_NEAR(static_cast<double>(stkMeshBulkData.identifier(nodesInFirstPart[i])), *fieldPtr, 1e-10);
  }
}

TEST(UnitTestPartsAfterCommit, SelectorOps)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  if(numProcs == 1)
  {
    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    const std::string generatedMeshSpecification = "generated:1x1x1";
    stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();

    stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
    stk::mesh::Part& partA = stkMeshMetaData.declare_part("partA", stk::topology::NODE_RANK);

    stkMeshIoBroker.populate_bulk_data();

    stk::mesh::BulkData& stkMeshBulkData = stkMeshIoBroker.bulk_data();
    stk::mesh::EntityVector nodes;
    stk::mesh::get_entities(stkMeshBulkData, stk::topology::NODE_RANK, nodes);
    ASSERT_EQ(8u, nodes.size());

    stkMeshBulkData.modification_begin();
    stk::mesh::PartVector addParts(1, &partA);
    for(size_t i = 0; i < 4; i++)
    {
      stkMeshBulkData.change_entity_parts(nodes[i], addParts);
    }
    stkMeshBulkData.modification_end();

    stk::mesh::Selector selectNodesNotInPartA = stkMeshMetaData.universal_part() - partA;
    stk::mesh::EntityVector nodesNotInPartA;
    stk::mesh::get_selected_entities(selectNodesNotInPartA,
                                     stkMeshBulkData.buckets(stk::topology::NODE_RANK),
                                     nodesNotInPartA);
    ASSERT_EQ(4u, nodesNotInPartA.size());

    stk::mesh::Part& partB = stkMeshMetaData.declare_part("partB", stk::topology::ELEMENT_RANK);

    stkMeshBulkData.modification_begin();
    addParts[0] = &partB;
    for(size_t i = 2; i < 6; i++)
    {
      stkMeshBulkData.change_entity_parts(nodes[i], addParts);
    }
    stkMeshBulkData.modification_end();

    stk::mesh::EntityVector nodesStillNotInPartA;
    stk::mesh::get_selected_entities(selectNodesNotInPartA,
                                     stkMeshBulkData.buckets(stk::topology::NODE_RANK),
                                     nodesStillNotInPartA);
    ASSERT_EQ(4u, nodesStillNotInPartA.size());
    for(size_t i = 0; i < nodesNotInPartA.size(); i++)
    {
      EXPECT_EQ(nodesNotInPartA[i], nodesStillNotInPartA[i]);
    }
  }
}
