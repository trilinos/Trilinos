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
#include <stk_io/FillMesh.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Selector.hpp>   // for Selector, operator<<, etc
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_topology/topology.hpp>    // for topology, etc
#include <stk_mesh/base/Types.hpp>      // for BucketVector, PartVector
#include <stk_unit_test_utils/TextMesh.hpp>
#include <string>                       // for string
namespace stk { namespace mesh { class Part; } }

namespace
{
//-BEGIN
TEST(StkMeshHowTo, basicSelectorUsage)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(communicator) != 1) { GTEST_SKIP(); }
  std::unique_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(communicator)
                                                   .set_spatial_dimension(3).create();
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();

//create a simple shell-quad-4 mesh:
//       6
// 3*----*----*9
//  | E2 | E4 |
//  |    |    |
// 2*---5*----*8
//  | E1 | E3 |
//  |    |    |
// 1*----*----*7
//       4
//

  std::string meshDesc = "0,1,SHELL_QUAD_4, 1,4,2,5, block_1\n"
                         "0,2,SHELL_QUAD_4, 2,5,6,3, block_2\n"
                         "0,3,SHELL_QUAD_4, 4,7,8,5, block_3\n"
                         "0,4,SHELL_QUAD_4, 5,8,9,6, block_4\n";
  stk::unit_test_util::setup_text_mesh(*bulkPtr, meshDesc);

  stk::mesh::Part& block_1 = *meta.get_part("block_1");
  stk::mesh::Part& block_2 = *meta.get_part("block_2");
  stk::mesh::Part& block_3 = *meta.get_part("block_3");
  stk::mesh::Part& block_4 = *meta.get_part("block_4");
  stk::mesh::PartVector allBlocks = {&block_1, &block_2, &block_3, &block_4};
 
  stk::mesh::Selector allNodes = stk::mesh::selectUnion(allBlocks);
  stk::mesh::Selector onlyCenterNode = stk::mesh::selectIntersection(allBlocks);
  stk::mesh::Selector nodes456 = (block_1 | block_2) & (block_3 | block_4);
  stk::mesh::Selector nodes123 = (block_1 | block_2) - nodes456;

  EXPECT_EQ(9u, stk::mesh::count_entities(*bulkPtr, stk::topology::NODE_RANK, allNodes));
  EXPECT_EQ(1u, stk::mesh::count_entities(*bulkPtr, stk::topology::NODE_RANK, onlyCenterNode));
  EXPECT_EQ(3u, stk::mesh::count_entities(*bulkPtr, stk::topology::NODE_RANK, nodes456));
  EXPECT_EQ(3u, stk::mesh::count_entities(*bulkPtr, stk::topology::NODE_RANK, nodes123));
}

TEST(StkMeshHowTo, betterUnderstandSelectorConstruction)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(communicator) != 1) { GTEST_SKIP(); }
  std::unique_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(communicator).create();
  const std::string generatedCubeMeshSpecification = "generated:1x1x1";
  stk::io::fill_mesh(generatedCubeMeshSpecification, *bulkPtr);

  stk::mesh::Selector nothingSelector_byDefaultConstruction;
  size_t expectingZeroBuckets = 0;
  EXPECT_EQ(expectingZeroBuckets, bulkPtr->get_buckets(stk::topology::NODE_RANK, nothingSelector_byDefaultConstruction).size());

  std::ostringstream readableSelectorDescription;
  readableSelectorDescription << nothingSelector_byDefaultConstruction;
  EXPECT_STREQ("NOTHING", readableSelectorDescription.str().c_str());

  stk::mesh::Selector allSelector(!nothingSelector_byDefaultConstruction);
  size_t numberOfAllNodeBuckets = bulkPtr->buckets(stk::topology::NODE_RANK).size();
  EXPECT_EQ(numberOfAllNodeBuckets, bulkPtr->get_buckets(stk::topology::NODE_RANK, allSelector).size());
}

TEST(StkMeshHowTo, makeSureYouAreNotIntersectingNothingSelector)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(communicator) != 1) { return; }
  std::unique_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(communicator).create();
  // syntax creates faces for surface on the positive: 'x-side', 'y-side', and 'z-side'
  // of a 1x1x1 cube, these parts are given the names: 'surface_1', 'surface_2', and 'surface_3'
  const std::string generatedCubeMeshSpecification = "generated:1x1x1|sideset:XYZ";
  stk::io::fill_mesh(generatedCubeMeshSpecification, *bulkPtr);

  stk::mesh::MetaData &stkMeshMetaData = bulkPtr->mesh_meta_data();
  stk::mesh::Part *surface1Part = stkMeshMetaData.get_part("surface_1");
  stk::mesh::Part *surface2Part = stkMeshMetaData.get_part("surface_2");
  stk::mesh::Part *surface3Part = stkMeshMetaData.get_part("surface_3");
  stk::mesh::PartVector allSurfaces;
  allSurfaces.push_back(surface1Part);
  allSurfaces.push_back(surface2Part);
  allSurfaces.push_back(surface3Part);

  stk::mesh::Selector selectorIntersectingNothing;
  for (size_t surfaceIndex = 0; surfaceIndex < allSurfaces.size(); ++surfaceIndex) {
    stk::mesh::Part &surfacePart = *(allSurfaces[surfaceIndex]);
    stk::mesh::Selector surfaceSelector(surfacePart);
    selectorIntersectingNothing &= surfacePart;
  }

  size_t expectedNumberOfBucketsWhenIntersectingNothing = 0;
  stk::mesh::BucketVector selectedBuckets = bulkPtr->get_buckets(stk::topology::NODE_RANK, selectorIntersectingNothing);
  EXPECT_EQ(expectedNumberOfBucketsWhenIntersectingNothing, selectedBuckets.size());

  stk::mesh::Selector preferredBoundaryNodesSelector = stk::mesh::selectIntersection(allSurfaces);
  size_t expectedNumberOfNodeBucketsWhenIntersectingAllSurfaces = 1;
  selectedBuckets = bulkPtr->get_buckets(stk::topology::NODE_RANK, preferredBoundaryNodesSelector);
  EXPECT_EQ(expectedNumberOfNodeBucketsWhenIntersectingAllSurfaces, selectedBuckets.size());
}
//-END
}
