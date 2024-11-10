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
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, entity_rank_names, etc
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/ForEachEntity.hpp"
#include "stk_mesh/base/Types.hpp"      // for EntityId
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_io/IossBridge.hpp"
#include "stk_unit_test_utils/TextMesh.hpp"

namespace {

constexpr unsigned SpatialDimension = 3;

void create_two_tet_element_mesh(stk::mesh::BulkData &bulk)
{
  std::string meshSpec = "0, 1,TET_4, 1,2,3,4\n"
                         "0, 2,TET_4, 2,3,4,5";
  stk::unit_test_util::setup_text_mesh(bulk, meshSpec);
}

//BEGINUseFieldBLAS
TEST(stkMeshHowTo, useFieldBLAS)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  builder.set_spatial_dimension(SpatialDimension);
  std::unique_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
  stk::mesh::MetaData& metaData = bulkPtr->mesh_meta_data();

  typedef stk::mesh::Field<double> DoubleField;
  DoubleField& pressureField = metaData.declare_field<double>(stk::topology::ELEM_RANK, "pressure");
  DoubleField& displacementsField = metaData.declare_field<double>(stk::topology::NODE_RANK, "displacements");
  DoubleField& velocityField = metaData.declare_field<double>(stk::topology::NODE_RANK, "velocity");

  double initialPressureValue = 4.4;
  constexpr unsigned numValuesPerNode = 3;
  stk::mesh::put_field_on_entire_mesh_with_initial_value(pressureField, &initialPressureValue);
  stk::mesh::put_field_on_mesh(displacementsField, metaData.universal_part(), numValuesPerNode, nullptr);
  stk::mesh::put_field_on_mesh(velocityField, metaData.universal_part(), numValuesPerNode, nullptr);

  create_two_tet_element_mesh(*bulkPtr);

  //incompatible fields, elem-rank vs node-rank
  EXPECT_ANY_THROW(stk::mesh::field_copy(pressureField, displacementsField));

  stk::mesh::field_fill(99.0, displacementsField);
  stk::mesh::field_fill(10.0, velocityField);
  const double alpha = 5.0;
  stk::mesh::field_axpy(alpha, displacementsField, velocityField);

  const double expectedVal = 10.0 + alpha*99.0;

  auto expectEqualVal = [&](const stk::mesh::BulkData& bulk, stk::mesh::Entity node) {
    const double* velocityDataForNode = stk::mesh::field_data(velocityField, node);
    for(unsigned i=0; i<numValuesPerNode; ++i) {
      EXPECT_NEAR(expectedVal, velocityDataForNode[i], 1.e-8);
    }
  };

  stk::mesh::for_each_entity_run(*bulkPtr, stk::topology::NODE_RANK, expectEqualVal);
}
//ENDUseFieldBLAS

}
