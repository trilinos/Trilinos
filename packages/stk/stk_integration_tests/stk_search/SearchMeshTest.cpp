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
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_search/LocalCoarseSearch.hpp>
#include <stk_unit_test_utils/Search_UnitTestUtils.hpp>
#include <stk_unit_test_utils/MeshUtilsForBoundingVolumes.hpp>
#include <stk_unit_test_utils/getOption.h>

namespace
{

TEST(StkSearch, NGP_coarse_search_mesh_elem_boxes_MORTON)
{
  using ExecSpace = Kokkos::DefaultExecutionSpace;
  MPI_Comm comm = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::mesh::MeshBuilder builder(comm);
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
 
  stk::io::fill_mesh("generated:1x9x19|sideset:xXyYzZ", *bulkPtr);
 
  Kokkos::View<FloatBoxIdent*, ExecSpace> elemBoxes =
      createBoundingBoxesForEntities<FloatBoxIdent>(*bulkPtr, stk::topology::ELEM_RANK);
  Kokkos::View<FloatBoxIdent*, ExecSpace> faceBoxes =
      createBoundingBoxesForEntities<FloatBoxIdent>(*bulkPtr, stk::topology::FACE_RANK);

  std::cout<<"Num elem-boxes: "<<elemBoxes.size()<<", num face-boxes: "<<faceBoxes.size()<<std::endl;

  stk::search::SearchMethod searchMethod = stk::search::MORTON_LBVH;
  Kokkos::View<IdentIntersection*, ExecSpace> searchResults;
  stk::search::local_coarse_search(elemBoxes, faceBoxes, searchMethod, searchResults, ExecSpace{});

  const size_t expectedSize = 2910;
  EXPECT_EQ(expectedSize, searchResults.size())<<"expected results size: "<<expectedSize<<", actual results size: "<<searchResults.size();
}

}
