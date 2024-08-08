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

#include <stddef.h>                     // for size_t
#include <stdlib.h>                     // for exit
#include <exception>                    // for exception
#include <iostream>                     // for ostringstream, etc
#include <iterator>                     // for distance
#include <map>                          // for _Rb_tree_const_iterator, etc
#include <stdexcept>                    // for logic_error, runtime_error
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/FieldParallel.hpp>  // for communicate_field_data, etc
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities, etc
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <gtest/gtest.h>
#include <string>                       // for string, basic_string, etc
#include <utility>                      // for pair
#include <vector>                       // for vector, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket, has_superset
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data, etc
#include "stk_mesh/base/Ghosting.hpp"   // for Ghosting
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, entity_rank_names, etc
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Relation.hpp"
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator|
#include "stk_mesh/base/Types.hpp"      // for EntityProc, EntityVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/util/PairIter.hpp"   // for PairIter
#include "stk_io/StkMeshIoBroker.hpp"
#include "stk_io/FillMesh.hpp"
#include "stk_mesh/base/MeshUtils.hpp"
#include "stk_unit_test_utils/BuildMesh.hpp"

using stk::unit_test_util::build_mesh;

TEST(UnitTestGhosting, WithChangeParts)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  if(numProcs == 3)
  {
    const int spatialDim = 3;
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, communicator);
    stk::mesh::MetaData& stkMeshMetaData = bulkPtr->mesh_meta_data();
    stk::mesh::BulkData& stkMeshBulkData = *bulkPtr;
    const std::string generatedMeshSpecification = "generated:1x1x3|sideset:xXyYzZ";
    stk::io::fill_mesh(generatedMeshSpecification, stkMeshBulkData);

    stk::mesh::EntityVector elementsOnProc;
    stk::mesh::get_selected_entities(stkMeshMetaData.locally_owned_part(),
                                     stkMeshBulkData.buckets(stk::topology::ELEMENT_RANK),
                                     elementsOnProc);
    ASSERT_EQ(1u, elementsOnProc.size());

    stk::mesh::EntityProcVec elementsToGhost;
    if(stkMeshBulkData.parallel_rank() == 0)
    {
      const int ghostToProc2 = 2;
      elementsToGhost.push_back(std::make_pair(elementsOnProc[0], ghostToProc2));
    }

    stkMeshBulkData.modification_begin();
    stk::mesh::Ghosting &ghostElemFrom0To2 = stkMeshBulkData.create_ghosting("ghostElemFrom0to2");
    stkMeshBulkData.change_ghosting(ghostElemFrom0To2, elementsToGhost);
    stkMeshBulkData.modification_end();

    std::vector<size_t> entityCounts;
    stk::mesh::count_entities(stkMeshMetaData.universal_part(), stkMeshBulkData, entityCounts);
    if(stkMeshBulkData.parallel_rank() == 0)
    {
      EXPECT_EQ(12u, entityCounts[stk::topology::NODE_RANK]);             //       __ __
      EXPECT_EQ( 0u, entityCounts[stk::topology::EDGE_RANK]);             //      |  |G |
      EXPECT_EQ( 9u, entityCounts[stk::topology::FACE_RANK]);             //      |__|__|
      EXPECT_EQ( 2u, entityCounts[stk::topology::ELEMENT_RANK]);
    }
    else if(stkMeshBulkData.parallel_rank() == 1)
    {
      EXPECT_EQ(16u, entityCounts[stk::topology::NODE_RANK]);             //       __ __ __
      EXPECT_EQ( 0u, entityCounts[stk::topology::EDGE_RANK]);             //      |G |  |G |
      EXPECT_EQ(14u, entityCounts[stk::topology::FACE_RANK]);             //      |__|__|__|
      EXPECT_EQ( 3u, entityCounts[stk::topology::ELEMENT_RANK]);
    }
    else
    {
      EXPECT_EQ(16u, entityCounts[stk::topology::NODE_RANK]);             //       __ __ __
      EXPECT_EQ( 0u, entityCounts[stk::topology::EDGE_RANK]);             //      |G |G |  |
      EXPECT_EQ(14u, entityCounts[stk::topology::FACE_RANK]);             //      |__|__|__|
      EXPECT_EQ( 3u, entityCounts[stk::topology::ELEMENT_RANK]);
    }

    stk::mesh::Part &newPart = stkMeshMetaData.declare_part("newPart");

    stkMeshBulkData.modification_begin();
    if(stkMeshBulkData.parallel_rank() == 0)
    {
      stk::mesh::PartVector partsToAdd(1, &newPart);
      stk::mesh::EntityVector sharedNodes;
      stk::mesh::get_selected_entities(stkMeshMetaData.globally_shared_part(),
                                       stkMeshBulkData.buckets(stk::topology::NODE_RANK),
                                       sharedNodes);
      for(size_t i=0; i<sharedNodes.size(); i++)
      {
        stkMeshBulkData.change_entity_parts(sharedNodes[i], partsToAdd);
      }
    }
    stkMeshBulkData.modification_end();

    stk::mesh::count_entities(stkMeshMetaData.universal_part(), stkMeshBulkData, entityCounts);
    if(stkMeshBulkData.parallel_rank() == 0)
    {
      EXPECT_EQ(12u, entityCounts[stk::topology::NODE_RANK]);             //       __ __
      EXPECT_EQ( 0u, entityCounts[stk::topology::EDGE_RANK]);             //      |  |G |
      EXPECT_EQ( 9u, entityCounts[stk::topology::FACE_RANK]);             //      |__|__|
      EXPECT_EQ( 2u, entityCounts[stk::topology::ELEMENT_RANK]);
    }
    else if(stkMeshBulkData.parallel_rank() == 1)
    {
      EXPECT_EQ(16u, entityCounts[stk::topology::NODE_RANK]);             //       __ __ __
      EXPECT_EQ( 0u, entityCounts[stk::topology::EDGE_RANK]);             //      |G |  |G |
      EXPECT_EQ(14u, entityCounts[stk::topology::FACE_RANK]);             //      |__|__|__|
      EXPECT_EQ( 3u, entityCounts[stk::topology::ELEMENT_RANK]);
    }
    else
    {
      EXPECT_EQ(16u, entityCounts[stk::topology::NODE_RANK]);             //          __ __
      EXPECT_EQ( 0u, entityCounts[stk::topology::EDGE_RANK]);             //      |  |G |  |
      EXPECT_EQ(14u, entityCounts[stk::topology::FACE_RANK]);             //      |  |__|__|
      EXPECT_EQ( 3u, entityCounts[stk::topology::ELEMENT_RANK]);          //      ^------------one face still ghosted
    }
  }
}

TEST(UnitTestGhosting, WithDeclareConstraintRelatedToRecvGhostNode)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  if(numProcs == 3)
  {
    const int spatialDim = 3;
    stk::mesh::MeshBuilder builder(communicator);
    builder.set_spatial_dimension(spatialDim);
    std::vector<std::string> rank_names(5);
    rank_names[0] = "node";
    rank_names[1] = "edge";
    rank_names[2] = "face";
    rank_names[3] = "elem";
    rank_names[4] = "constraint";
    builder.set_entity_rank_names(rank_names);
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
    stk::mesh::MetaData& stkMeshMetaData = bulkPtr->mesh_meta_data();
    stk::mesh::BulkData& stkMeshBulkData = *bulkPtr;
    const std::string generatedMeshSpecification = "generated:1x1x3|sideset:xXyYzZ";
    stk::io::fill_mesh(generatedMeshSpecification, stkMeshBulkData);

    stk::mesh::EntityVector elementsOnProc;
    stk::mesh::get_selected_entities(stkMeshMetaData.locally_owned_part(),
                                     stkMeshBulkData.buckets(stk::topology::ELEMENT_RANK),
                                     elementsOnProc);
    ASSERT_EQ(1u, elementsOnProc.size());

    stk::mesh::EntityProcVec elementsToGhost;
    if(stkMeshBulkData.parallel_rank() == 0)
    {
      const int ghostToProc2 = 2;
      elementsToGhost.push_back(std::make_pair(elementsOnProc[0], ghostToProc2));
    }

    stkMeshBulkData.modification_begin();
    stk::mesh::Ghosting &ghostElemFrom0To2 = stkMeshBulkData.create_ghosting("ghostElemFrom0to2");
    stkMeshBulkData.change_ghosting(ghostElemFrom0To2, elementsToGhost);
    stkMeshBulkData.modification_end();

    std::vector<size_t> entityCounts;
    stk::mesh::count_entities(stkMeshMetaData.universal_part(), stkMeshBulkData, entityCounts);
    if(stkMeshBulkData.parallel_rank() == 0)
    {
      EXPECT_EQ(12u, entityCounts[stk::topology::NODE_RANK]);             //       __ __
      EXPECT_EQ( 0u, entityCounts[stk::topology::EDGE_RANK]);             //      |  |G |
      EXPECT_EQ( 9u, entityCounts[stk::topology::FACE_RANK]);             //      |__|__|
      EXPECT_EQ( 2u, entityCounts[stk::topology::ELEMENT_RANK]);
    }
    else if(stkMeshBulkData.parallel_rank() == 1)
    {
      EXPECT_EQ(16u, entityCounts[stk::topology::NODE_RANK]);             //       __ __ __
      EXPECT_EQ( 0u, entityCounts[stk::topology::EDGE_RANK]);             //      |G |  |G |
      EXPECT_EQ(14u, entityCounts[stk::topology::FACE_RANK]);             //      |__|__|__|
      EXPECT_EQ( 3u, entityCounts[stk::topology::ELEMENT_RANK]);
    }
    else
    {
      EXPECT_EQ(16u, entityCounts[stk::topology::NODE_RANK]);             //       __ __ __
      EXPECT_EQ( 0u, entityCounts[stk::topology::EDGE_RANK]);             //      |G |G |  |
      EXPECT_EQ(14u, entityCounts[stk::topology::FACE_RANK]);             //      |__|__|__|
      EXPECT_EQ( 3u, entityCounts[stk::topology::ELEMENT_RANK]);
    }

    stkMeshBulkData.modification_begin();
    if(stkMeshBulkData.parallel_rank() == 2)
    {
      stk::mesh::EntityId constraintId = 1;
      stk::mesh::Entity constraint = stkMeshBulkData.declare_constraint(constraintId);
      stk::mesh::Entity node1 = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, 1);
      EXPECT_TRUE(stkMeshBulkData.bucket(node1).member(stkMeshBulkData.ghosting_part(ghostElemFrom0To2)));
      stkMeshBulkData.declare_relation(constraint, node1, 0);
    }

    fixup_ghosted_to_shared_nodes(stkMeshBulkData);
    EXPECT_NO_THROW(stkMeshBulkData.modification_end());

    stk::mesh::count_entities(stkMeshMetaData.universal_part(), stkMeshBulkData, entityCounts);
    if(stkMeshBulkData.parallel_rank() == 0)
    {
      stk::mesh::EntityId constraintId = 1;
      stk::mesh::Entity constraint = stkMeshBulkData.get_entity(stk::topology::CONSTRAINT_RANK, constraintId);
      EXPECT_TRUE(stkMeshBulkData.bucket(constraint).in_aura());
      stk::mesh::Entity node1 = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, 1);
      EXPECT_TRUE(stkMeshBulkData.bucket(node1).shared());
      EXPECT_TRUE(stkMeshBulkData.bucket(node1).owned());
      EXPECT_EQ(12u, entityCounts[stk::topology::NODE_RANK]);             //       __ __
      EXPECT_EQ( 0u, entityCounts[stk::topology::EDGE_RANK]);             //      |  |G |
      EXPECT_EQ( 9u, entityCounts[stk::topology::FACE_RANK]);             //      |__|__|
      EXPECT_EQ( 2u, entityCounts[stk::topology::ELEMENT_RANK]);
      EXPECT_EQ( 1u, entityCounts[stk::topology::CONSTRAINT_RANK]);
    }
    else if(stkMeshBulkData.parallel_rank() == 1)
    {
      EXPECT_EQ(16u, entityCounts[stk::topology::NODE_RANK]);             //       __ __ __
      EXPECT_EQ( 0u, entityCounts[stk::topology::EDGE_RANK]);             //      |G |  |G |
      EXPECT_EQ(14u, entityCounts[stk::topology::FACE_RANK]);             //      |__|__|__|
      EXPECT_EQ( 3u, entityCounts[stk::topology::ELEMENT_RANK]);
      EXPECT_EQ( 0u, entityCounts[stk::topology::CONSTRAINT_RANK]);
    }
    else
    {
      stk::mesh::EntityId constraintId = 1;
      stk::mesh::Entity constraint = stkMeshBulkData.get_entity(stk::topology::CONSTRAINT_RANK, constraintId);
      EXPECT_TRUE(stkMeshBulkData.bucket(constraint).owned());
      stk::mesh::Entity node1 = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, 1);
      EXPECT_TRUE(stkMeshBulkData.bucket(node1).shared());
      EXPECT_TRUE(!stkMeshBulkData.bucket(node1).owned());
      EXPECT_EQ(16u, entityCounts[stk::topology::NODE_RANK]);             //       __ __ __
      EXPECT_EQ( 0u, entityCounts[stk::topology::EDGE_RANK]);             //      |G |G |  |
      EXPECT_EQ(14u, entityCounts[stk::topology::FACE_RANK]);             //      |__|__|__|
      EXPECT_EQ( 3u, entityCounts[stk::topology::ELEMENT_RANK]);
      EXPECT_EQ( 1u, entityCounts[stk::topology::CONSTRAINT_RANK]);
    }
  }
}

