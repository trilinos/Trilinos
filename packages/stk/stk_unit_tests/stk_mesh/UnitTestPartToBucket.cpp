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
#include <stddef.h>                     // for size_t, NULL
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Selector.hpp>   // for operator&, Selector, etc
#include <stk_topology/topology.hpp>    // for topology, etc
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_size, etc
#include <string>                       // for string
#include "mpi.h"                        // for MPI_COMM_WORLD, etc
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Types.hpp"      // for BucketVector, EntityId

namespace
{
void runTwoHexParallelBucketTests(const std::string &generatedMeshSpecification, const size_t expectedNumBucketsPerSlice);

TEST(PartToBucket, hex)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(communicator) != 1)
  {
    return;
  }
  stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
  const std::string generatedMeshSpecification = "generated:1x1x1";
  stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
  stkMeshIoBroker.create_input_mesh();

  stkMeshIoBroker.populate_bulk_data();

  stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
  stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

  const stk::mesh::BucketVector &nodeBuckets = stkMeshBulkData.buckets(stk::topology::NODE_RANK);
  size_t expectedNodeBuckets = 1;
  EXPECT_EQ(expectedNodeBuckets, nodeBuckets.size());

  const stk::mesh::Part& locallyOwned = stkMeshMetaData.locally_owned_part();
  const stk::mesh::Part& globallyShared = stkMeshMetaData.globally_shared_part();
  const stk::mesh::Part& universal = stkMeshMetaData.universal_part();

  const stk::mesh::Bucket& nodeBucket = *nodeBuckets[0];

  EXPECT_TRUE(nodeBucket.member(locallyOwned));
  EXPECT_TRUE(nodeBucket.member(universal));
  EXPECT_FALSE(nodeBucket.member(globallyShared));

  const stk::mesh::BucketVector &elemBuckets = stkMeshBulkData.buckets(stk::topology::ELEMENT_RANK);
  size_t expectedElemBuckets = 1;
  EXPECT_EQ(expectedElemBuckets, elemBuckets.size());

  const stk::mesh::Bucket& elemBucket = *elemBuckets[0];

  EXPECT_TRUE(elemBucket.member(locallyOwned));
  EXPECT_TRUE(elemBucket.member(universal));
  EXPECT_FALSE(elemBucket.member(globallyShared));

  const stk::mesh::BucketVector &edgeBuckets = stkMeshBulkData.buckets(stk::topology::EDGE_RANK);
  size_t expectedEdgeBuckets = 0;
  EXPECT_EQ(expectedEdgeBuckets, edgeBuckets.size());

  const stk::mesh::BucketVector &faceBuckets = stkMeshBulkData.buckets(stk::topology::FACE_RANK);
  size_t expectedFaceBuckets = 0;
  EXPECT_EQ(expectedFaceBuckets, faceBuckets.size());
}

TEST(PartToBucket, hexWithSingleSideset)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(communicator) != 1)
  {
    return;
  }
  stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
  const std::string generatedMeshSpecification = "generated:1x1x1|sideset:X";
  stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
  stkMeshIoBroker.create_input_mesh();

  stkMeshIoBroker.populate_bulk_data();

  stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
  stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

  stk::mesh::Part &surface1Part = *stkMeshMetaData.get_part("surface_1");
  stk::mesh::Selector surface1NodesSelector(surface1Part);

  const stk::mesh::BucketVector &surface1NodeBuckets = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, surface1NodesSelector);
  size_t expectedSurface1NodeBuckets = 1;
  EXPECT_EQ(expectedSurface1NodeBuckets, surface1NodeBuckets.size());

  size_t numNodesInBoundary = 4;
  EXPECT_EQ(numNodesInBoundary, surface1NodeBuckets[0]->size());

  const stk::mesh::BucketVector &nodeBuckets = stkMeshBulkData.buckets(stk::topology::NODE_RANK);
  size_t expectedNodeBuckets = 2;
  EXPECT_EQ(expectedNodeBuckets, nodeBuckets.size());

  const stk::mesh::BucketVector &elemBuckets = stkMeshBulkData.buckets(stk::topology::ELEMENT_RANK);
  size_t expectedElemBuckets = 1;
  EXPECT_EQ(expectedElemBuckets, elemBuckets.size());

  const stk::mesh::BucketVector &edgeBuckets = stkMeshBulkData.buckets(stk::topology::EDGE_RANK);
  size_t expectedEdgeBuckets = 0;
  EXPECT_EQ(expectedEdgeBuckets, edgeBuckets.size());

  const stk::mesh::BucketVector &faceBuckets = stkMeshBulkData.buckets(stk::topology::FACE_RANK);
  size_t expectedFaceBuckets = 1;
  EXPECT_EQ(expectedFaceBuckets, faceBuckets.size());

  size_t numFacesInBucket = 1;
  EXPECT_EQ(numFacesInBucket, faceBuckets[0]->size());
}

/*
       C 7--------------8 D
         |\             |\
         | \      s2    | \
         |  \           |  \          y
         | C 3--------------4 D    z  |
         |   |          |   |       \ |
         |   |          | s1|        \|
       A 5---|----------6 B |         o------x
          \  |           \  |
           \ |            \ |
            \|             \|
           A 1--------------2 B
*/

TEST(PartToBucket, hexWithTwoSidesets)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(communicator) != 1)
  {
    return;
  }
  stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
  const std::string generatedMeshSpecification = "generated:1x1x1|sideset:XY";
  stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
  stkMeshIoBroker.create_input_mesh();

  stkMeshIoBroker.populate_bulk_data();

  stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
  stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

  stk::mesh::Part &surface1Part = *stkMeshMetaData.get_part("surface_1");
  stk::mesh::Selector surface1NodesSelector(surface1Part);

  const stk::mesh::BucketVector &surface1NodeBuckets = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, surface1NodesSelector);
  size_t expectedSurface1NodeBuckets = 2;
  EXPECT_EQ(expectedSurface1NodeBuckets, surface1NodeBuckets.size());

  const stk::mesh::BucketVector &nodeBuckets = stkMeshBulkData.buckets(stk::topology::NODE_RANK);
  size_t expectedNodeBuckets = 4; // A,B,C,D in above picture
  EXPECT_EQ(expectedNodeBuckets, nodeBuckets.size());

  const stk::mesh::BucketVector &elemBuckets = stkMeshBulkData.buckets(stk::topology::ELEMENT_RANK);
  size_t expectedElemBuckets = 1;
  EXPECT_EQ(expectedElemBuckets, elemBuckets.size());

  const stk::mesh::BucketVector &edgeBuckets = stkMeshBulkData.buckets(stk::topology::EDGE_RANK);
  size_t expectedEdgeBuckets = 0;
  EXPECT_EQ(expectedEdgeBuckets, edgeBuckets.size());

  const stk::mesh::BucketVector &faceBuckets = stkMeshBulkData.buckets(stk::topology::FACE_RANK);
  size_t expectedFaceBuckets = 2;
  EXPECT_EQ(expectedFaceBuckets, faceBuckets.size());

  size_t numFacesInBucket = 1;
  EXPECT_EQ(numFacesInBucket, faceBuckets[0]->size());
}

/*
    1-12    = nodes
    A       = bucket

                    A 11------------12 A
                      |\             |\                  Locally Owned
                      | \            | \                     x  x
                      |  \           |  \                 x        x
                      | A 7--------------8 A             x    A     x
                      |   |\         |   |\              x          x
                      |   | \        |   | \              x        x
                    A 9---|--\------10 A |  \                x  x
                       \  | A 3--------------4 A
                        \ |   |        \ |   |
                         \|   |         \|   |
                        A 5---|----------6 A |
                           \  |           \  |
                            \ |            \ |
                             \|             \|
                            A 1--------------2 A
*/
TEST(PartToBucket, twoHex)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(communicator) != 1)
  {
    return;
  }
  stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
  const std::string generatedMeshSpecification = "generated:1x1x2";
  stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
  stkMeshIoBroker.create_input_mesh();
  stkMeshIoBroker.populate_bulk_data();

  stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
  stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

  const stk::mesh::BucketVector &nodeBuckets = stkMeshBulkData.buckets(stk::topology::NODE_RANK);
  size_t expectedNumNodeBuckets = 1;
  EXPECT_EQ(expectedNumNodeBuckets, nodeBuckets.size());
  size_t expectedNumNodes = 12;
  EXPECT_EQ(expectedNumNodes, nodeBuckets[0]->size());

  stk::mesh::Selector locallyOwnedSelector(stkMeshMetaData.locally_owned_part());
  stk::mesh::Selector globallySharedSelector(stkMeshMetaData.globally_shared_part());
  stk::mesh::Selector ghostedSelector(!(stkMeshMetaData.locally_owned_part() | stkMeshMetaData.globally_shared_part()));

  const stk::mesh::BucketVector &locallyOwnedNodeBuckets = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, locallyOwnedSelector);
  const stk::mesh::BucketVector &globallySharedNodeBuckets = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, globallySharedSelector);
  const stk::mesh::BucketVector &ghostedNodeBuckets = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, ghostedSelector);
  EXPECT_EQ(1u, locallyOwnedNodeBuckets.size());
  EXPECT_EQ(0u, globallySharedNodeBuckets.size());
  EXPECT_EQ(0u, ghostedNodeBuckets.size());

  const stk::mesh::BucketVector &locallyOwnedElementBuckets = stkMeshBulkData.get_buckets(stk::topology::ELEM_RANK, locallyOwnedSelector);
  const stk::mesh::BucketVector &globallySharedElementBuckets = stkMeshBulkData.get_buckets(stk::topology::ELEM_RANK, globallySharedSelector);
  const stk::mesh::BucketVector &ghostedElementBuckets = stkMeshBulkData.get_buckets(stk::topology::ELEM_RANK, ghostedSelector);
  EXPECT_EQ(1u, locallyOwnedElementBuckets.size());
  EXPECT_EQ(0u, globallySharedElementBuckets.size());
  EXPECT_EQ(0u, ghostedElementBuckets.size());
}

/*
    1-12    = nodes
    A,B,C   = buckets

           C 11------------12 C
             |\             |\                  Locally Owned   Globally Shared      Ghosted
             | \            | \                     x  x             o  o             *  *
             |  \           |  \                 x        x       o        o       *        *
   proc 1    | B 7--------------8 B             x    C     x     o    B     o     *    A     *
             |   |`         |   |`              x          x     o          o     *          *
             |   | `        |   | `              x        x       o        o       *        *
           C 9---|---------10 C |  `                x  x             o  o             *  *
              \  |   `       \  |   `
               \ |    `       \ |    `
                \|     `       \|     `
               B 5--------------6 B    `
                  `      `       `      `                                       Globally Shared
                   `    B 7--------------8 B                                          /
                    `     |\       `     |\                         Locally Owned    /        Ghosted
                     `    | \       `    | \                                  x  x  /           *  *
                      `   |  \       `   |  \            y                 x    o  /x        *        *
             proc 0    `  | A 3--------------4 A      z  |                x    o     x      *          *
                        ` |   |        ` |   |         \ |                x A  o  B  x      *    C     *
                         `|   |         `|   |          \|                 x    o   x        *        *
                        B 5---|----------6 B |           o------x             x  x              *  *
                           \  |           \  |
                            \ |            \ |
                             \|             \|
                            A 1--------------2 A

    A different way of thinking about the following tests, consider looking down the Z-axis at 'slices':

       C 11------------12 C
         |              |
         |              |
         |              |
         |              |
         |              |
         |              |                          y
       C 9-------------10 C                        |
                                                   |
                B 7--------------8 B               |
                  |              |                 o------x
                  |              |
                  |              |
                  |              |
                  |              |
                  |              |
                B 5--------------6 B

                         A 3--------------4 A
                           |              |
                           |              |
                           |              |
                           |              |
                           |              |
                           |              |
                         A 1--------------2 A


*/
TEST(PartToBucket, np2TwoHex)
{
  const std::string generatedMeshSpecification = "generated:1x1x2";
  size_t expectedNumBucketsPerSlice = 1;
  runTwoHexParallelBucketTests(generatedMeshSpecification, expectedNumBucketsPerSlice);
}

/*
    1-12    = nodes
    A-L     = buckets

    2D XY-plane Slices:

       K 11------------12 L
         |              |
         |              |
         |              |
         |              |
         |              |
         |              |
       I 9-------------10 J                        y
                                                   |
                G 7--------------8 H               |
                  |              |                 |
                  |              |                 o------x
                  |              |
                  |              |
                  |              |
                  |              |
                E 5--------------6 F

                         C 3--------------4 D      Sidesets on the right (X) and top (Y)
                           |              |         - node 1 is not part of either sideset
                           |              |         - node 2 is part of X-sideset
                           |              |         - node 3 is part of Y-sideset
                           |              |         - node 4 is part of X-sideset and Y-sideset
                           |              |        so each node resides in it's own bucket per 'slice'
                           |              |
                         A 1--------------2 B

    12 buckets are needed for this problem (two hexes with two sidesets in parallel).
*/
TEST(PartToBucket, np2TwoHexTwoSidesets)
{
  const std::string generatedMeshSpecification = "generated:1x1x2|sideset:XY";
  size_t expectedNumBucketsPerSlice = 4;
  runTwoHexParallelBucketTests(generatedMeshSpecification, expectedNumBucketsPerSlice);
}

void runTwoHexParallelBucketTests(const std::string &generatedMeshSpecification, const size_t expectedNumBucketsPerSlice)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(communicator) != 2)
  {
    return;
  }
  stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
  stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
  stkMeshIoBroker.create_input_mesh();
  stkMeshIoBroker.populate_bulk_data();

  stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();
  stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();

  const stk::mesh::BucketVector &nodeBuckets = stkMeshBulkData.buckets(stk::topology::NODE_RANK);
  size_t numSlices = 3;
  size_t expectedNodeBuckets = numSlices * expectedNumBucketsPerSlice;
  ASSERT_EQ(expectedNodeBuckets, nodeBuckets.size());

  size_t numNodesPerSlice = 4;
  size_t expectedNumNodesInEveryBucket = numNodesPerSlice / expectedNumBucketsPerSlice;
  for(size_t i=0; i<nodeBuckets.size(); i++)
  {
    EXPECT_EQ(expectedNumNodesInEveryBucket, nodeBuckets[i]->size());
  }

  const stk::mesh::Part& locallyOwned = stkMeshMetaData.locally_owned_part();
  const stk::mesh::Part& globallyShared = stkMeshMetaData.globally_shared_part();
  stk::mesh::Selector locallyOwnedSelector(locallyOwned);
  stk::mesh::Selector globallySharedSelector(globallyShared);
  stk::mesh::Selector ghostedSelector(!(locallyOwned | globallyShared));
  const stk::mesh::BucketVector &locallyOwnedNodeBuckets = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, locallyOwnedSelector);
  const stk::mesh::BucketVector &globallySharedNodeBuckets = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, globallySharedSelector);
  const stk::mesh::BucketVector &ghostedNodeBuckets = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, ghostedSelector);

  size_t expectedNumLocallyOwnedNodeBuckets = 0;
  size_t expectedNumGloballySharedNodeBuckets = 0;
  size_t expectedNumGhostedNodeBuckets = 0;
  if(stk::parallel_machine_rank(communicator) == 0)
  {
    expectedNumLocallyOwnedNodeBuckets = 2 * expectedNumBucketsPerSlice;
    expectedNumGloballySharedNodeBuckets = expectedNumBucketsPerSlice;
    expectedNumGhostedNodeBuckets = expectedNumBucketsPerSlice;
  }
  else
  {
    expectedNumLocallyOwnedNodeBuckets = expectedNumBucketsPerSlice;
    expectedNumGloballySharedNodeBuckets = expectedNumBucketsPerSlice;
    expectedNumGhostedNodeBuckets = expectedNumBucketsPerSlice;
  }
  EXPECT_EQ(expectedNumLocallyOwnedNodeBuckets, locallyOwnedNodeBuckets.size());
  EXPECT_EQ(expectedNumGloballySharedNodeBuckets, globallySharedNodeBuckets.size());
  EXPECT_EQ(expectedNumGhostedNodeBuckets, ghostedNodeBuckets.size());


  const stk::mesh::BucketVector &locallyOwnedElementBuckets = stkMeshBulkData.get_buckets(stk::topology::ELEM_RANK, locallyOwnedSelector);
  const stk::mesh::BucketVector &globallySharedElementBuckets = stkMeshBulkData.get_buckets(stk::topology::ELEM_RANK, globallySharedSelector);
  const stk::mesh::BucketVector &ghostedElementBuckets = stkMeshBulkData.get_buckets(stk::topology::ELEM_RANK, ghostedSelector);
  EXPECT_EQ(1u, locallyOwnedElementBuckets.size());
  EXPECT_EQ(0u, globallySharedElementBuckets.size());
  EXPECT_EQ(1u, ghostedElementBuckets.size());
}

void checkNodeInSelectedBucket(stk::mesh::Selector selectNode, stk::mesh::EntityId expectedGlobalId, stk::mesh::BulkData &stkMeshBulkData)
{
  const stk::mesh::BucketVector &node1BucketVector = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, selectNode);
  size_t expectedNumBucketsForAnyNode = 1;
  EXPECT_EQ(expectedNumBucketsForAnyNode, node1BucketVector.size());
  size_t expectedNumNodesPerBucket = 1;
  EXPECT_EQ(expectedNumNodesPerBucket, node1BucketVector[0]->size());
  const stk::mesh::Bucket &node1Bucket = *node1BucketVector[0];
  EXPECT_EQ(expectedGlobalId, stkMeshBulkData.identifier(node1Bucket[0]));
}
/*
         7--------------8
         |\             |\
         | \      s2    | \
         |  \           |  \          y
         |   3--------------4      z  |
         |s1 |          |   |       \ |
         |   |          |   |        \|
         5---|----------6   |         o------x
          \  |     s3    \  |
           \ |            \ |
            \|             \|
             1--------------2
*/
TEST(PartToBucket, hexWithThreeSidesets)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(communicator) != 1)
  {
    return;
  }
  stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);

  //generated-mesh 'sideset:xYz' syntax adds face surfaces on 3 sides of the mesh,
  //(minimum 'x' side, maximum 'y' side, minimum 'z' side)
  //and the IO system will create corresponding stk::mesh::Parts named 'surface_1',
  //'surface_2' and 'surface_3', respectively, which are referenced in code below.
  const std::string generatedMeshSpecification = "generated:1x1x1|sideset:xYz";
  stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
  stkMeshIoBroker.create_input_mesh();
  stkMeshIoBroker.populate_bulk_data();

  stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
  stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

  stk::mesh::Part *surface1PartPtr = stkMeshMetaData.get_part("surface_1");
  stk::mesh::Part *surface2PartPtr = stkMeshMetaData.get_part("surface_2");
  stk::mesh::Part *surface3PartPtr = stkMeshMetaData.get_part("surface_3");

  EXPECT_TRUE(surface1PartPtr != NULL);
  EXPECT_TRUE(surface2PartPtr != NULL);
  EXPECT_TRUE(surface3PartPtr != NULL);

  stk::mesh::Part &surface1Part = *surface1PartPtr;
  stk::mesh::Part &surface2Part = *surface2PartPtr;
  stk::mesh::Part &surface3Part = *surface3PartPtr;

  stk::mesh::Selector surface1Selector(surface1Part);
  const stk::mesh::BucketVector &surface1NodeBuckets = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, surface1Selector);
  size_t expectedSurface1NodeBuckets = 4;
  EXPECT_EQ(expectedSurface1NodeBuckets, surface1NodeBuckets.size());

  stk::mesh::Selector selectNode1 = surface1Part & (!surface2Part) & surface3Part;
  stk::mesh::EntityId expectedGlobalId = 1;
  checkNodeInSelectedBucket(selectNode1, expectedGlobalId, stkMeshBulkData);

  stk::mesh::Selector selectNode2 = (!surface1Part) & (!surface2Part) & surface3Part;
  expectedGlobalId = 2;
  checkNodeInSelectedBucket(selectNode2, expectedGlobalId, stkMeshBulkData);

  stk::mesh::Selector selectNode3 = surface1Part & surface2Part & surface3Part;
  expectedGlobalId = 3;
  checkNodeInSelectedBucket(selectNode3, expectedGlobalId, stkMeshBulkData);

  stk::mesh::Selector selectNode4 = (!surface1Part) & surface2Part & surface3Part;
  expectedGlobalId = 4;
  checkNodeInSelectedBucket(selectNode4, expectedGlobalId, stkMeshBulkData);

  stk::mesh::Selector selectNode5 = surface1Part & (!surface2Part) & (!surface3Part);
  expectedGlobalId = 5;
  checkNodeInSelectedBucket(selectNode5, expectedGlobalId, stkMeshBulkData);

  stk::mesh::Part *block1PartPtr = stkMeshMetaData.get_part("block_1");
  EXPECT_TRUE(block1PartPtr != NULL);

  stk::mesh::Part &block1Part = *block1PartPtr;
  stk::mesh::Selector selectNode6 = block1Part & (!surface1Part) & (!surface2Part) & (!surface3Part);
  expectedGlobalId = 6;
  checkNodeInSelectedBucket(selectNode6, expectedGlobalId, stkMeshBulkData);

  stk::mesh::Selector selectNode7 = surface1Part & surface2Part & (!surface3Part);
  expectedGlobalId = 7;
  checkNodeInSelectedBucket(selectNode7, expectedGlobalId, stkMeshBulkData);

  stk::mesh::Selector selectNode8 = (!surface1Part) & surface2Part & (!surface3Part);
  expectedGlobalId = 8;
  checkNodeInSelectedBucket(selectNode8, expectedGlobalId, stkMeshBulkData);

  const stk::mesh::BucketVector &nodeBuckets = stkMeshBulkData.buckets(stk::topology::NODE_RANK);
  size_t expectedNodeBuckets = 8;
  EXPECT_EQ(expectedNodeBuckets, nodeBuckets.size());

  const stk::mesh::BucketVector &elemBuckets = stkMeshBulkData.buckets(stk::topology::ELEMENT_RANK);
  size_t expectedElemBuckets = 1;
  EXPECT_EQ(expectedElemBuckets, elemBuckets.size());

  const stk::mesh::BucketVector &edgeBuckets = stkMeshBulkData.buckets(stk::topology::EDGE_RANK);
  size_t expectedEdgeBuckets = 0;
  EXPECT_EQ(expectedEdgeBuckets, edgeBuckets.size());

  const stk::mesh::BucketVector &faceBuckets = stkMeshBulkData.buckets(stk::topology::FACE_RANK);
  size_t expectedFaceBuckets = 3;
  EXPECT_EQ(expectedFaceBuckets, faceBuckets.size());
}

}
