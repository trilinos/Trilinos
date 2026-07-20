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

#include <stk_io/FillMesh.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include "stk_mesh/base/EntityKey.hpp"
#include <stk_mesh/base/MeshBuilder.hpp>
#include "stk_search_util/ObjectCoordinates.hpp"
#include "stk_unit_test_utils/TextMesh.hpp"

#include <gtest/gtest.h>
#include <stk_util/parallel/Parallel.hpp>

#include <string>

namespace {

TEST(ObjectCoordinates, distance_from_nearest_entity_node)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> meshPtr = stk::mesh::MeshBuilder(comm).set_spatial_dimension(3).create();

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8"
          "|coordinates:   0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1";
 
  stk::unit_test_util::setup_text_mesh(*meshPtr, meshDesc);

  stk::mesh::Entity elem1 = meshPtr->get_entity(stk::topology::ELEM_RANK, 1);
  ASSERT_TRUE(meshPtr->is_valid(elem1));

  const stk::mesh::MetaData& meta = meshPtr->mesh_meta_data();
  const stk::mesh::FieldBase* coordField = meta.coordinate_field();

  constexpr int spatialDim = 3;
  std::array<double,spatialDim> pointAtNode1 = {0.0, 0.0, 0.0};

  EXPECT_DOUBLE_EQ(0.0, stk::search::distance_from_nearest_entity_node(*meshPtr, elem1, coordField, pointAtNode1.data()));
}

}

