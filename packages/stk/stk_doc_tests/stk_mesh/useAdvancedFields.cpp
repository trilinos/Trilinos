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
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field, etc
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/FieldBase.hpp"  // for field_scalars_per_entity, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_io/IossBridge.hpp"
#include "stk_unit_test_utils/TextMesh.hpp"

namespace {

//BEGINUseAdvancedFields
TEST(stkMeshHowTo, useAdvancedFields)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }

  const unsigned spatialDimension = 3;
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  builder.set_spatial_dimension(spatialDimension);
  std::unique_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
  stk::mesh::MetaData& metaData = bulkPtr->mesh_meta_data();

  typedef stk::mesh::Field<double> DoubleField;
  DoubleField& tensorField = metaData.declare_field<double>(stk::topology::ELEM_RANK, "tensor");
  DoubleField& variableSizeField = metaData.declare_field<double>(stk::topology::ELEM_RANK, "variableSizeField");

  stk::mesh::Part &tetPart = metaData.declare_part_with_topology("tetElementPart", stk::topology::TET_4);
  stk::mesh::Part &hexPart = metaData.declare_part_with_topology("hexElementPart", stk::topology::HEX_8);

  const int numTensorValues = 9;
  const int numCopies = 2;
  double initialTensorValue[] = { 1,  2,  3,  4,  5,  6,  7,  8,  9,
                                 11, 12, 13, 14, 15, 16, 17, 18, 19};
  stk::mesh::put_field_on_mesh(tensorField, metaData.universal_part(), numTensorValues, numCopies, initialTensorValue);
  stk::io::set_field_output_type(tensorField, stk::io::FieldOutputType::FULL_TENSOR_36);

  const int numVectorValues = 3;
  double initialVectorValue[] = {1, 2, 3, 11, 12, 13};
  stk::mesh::put_field_on_mesh(variableSizeField, tetPart, numVectorValues, initialVectorValue);
  stk::mesh::put_field_on_mesh(variableSizeField, hexPart, numVectorValues, numCopies, initialVectorValue);
  stk::io::set_field_output_type(variableSizeField, stk::io::FieldOutputType::VECTOR_3D);

  std::string meshSpec = "0,1,TET_4, 1,2,3,4, tetElementPart\n"
                         "0,2,HEX_8, 5,6,7,8,9,10,11,12, hexElementPart";
  stk::unit_test_util::setup_text_mesh(*bulkPtr, meshSpec);

  stk::mesh::Entity tetElem = bulkPtr->get_entity(stk::topology::ELEM_RANK, 1);
  stk::mesh::Entity hexElem = bulkPtr->get_entity(stk::topology::ELEM_RANK, 2);

  const int tensorScalarsPerTet = stk::mesh::field_scalars_per_entity(tensorField, tetElem);
  const int tensorScalarsPerHex = stk::mesh::field_scalars_per_entity(tensorField, hexElem);
  EXPECT_EQ(tensorScalarsPerTet, numTensorValues*numCopies);
  EXPECT_EQ(tensorScalarsPerHex, numTensorValues*numCopies);

  const int tensorExtent0PerTet = stk::mesh::field_extent0_per_entity(tensorField, tetElem);
  const int tensorExtent0PerHex = stk::mesh::field_extent0_per_entity(tensorField, hexElem);
  EXPECT_EQ(tensorExtent0PerTet, numTensorValues);
  EXPECT_EQ(tensorExtent0PerHex, numTensorValues);

  const int tensorExtent1PerTet = stk::mesh::field_extent1_per_entity(tensorField, tetElem);
  const int tensorExtent1PerHex = stk::mesh::field_extent1_per_entity(tensorField, hexElem);
  EXPECT_EQ(tensorExtent1PerTet, numCopies);
  EXPECT_EQ(tensorExtent1PerHex, numCopies);

  double* tensorData = stk::mesh::field_data(tensorField, hexElem);
  for (int i = 0; i < tensorScalarsPerHex; ++i) {
    EXPECT_EQ(initialTensorValue[i], tensorData[i]);
  }

  const int vectorScalarsPerTet = stk::mesh::field_scalars_per_entity(variableSizeField, tetElem);
  const int vectorScalarsPerHex = stk::mesh::field_scalars_per_entity(variableSizeField, hexElem);
  EXPECT_EQ(vectorScalarsPerTet, numVectorValues);
  EXPECT_EQ(vectorScalarsPerHex, numVectorValues*numCopies);

  const int vectorExtent0PerTet = stk::mesh::field_extent0_per_entity(variableSizeField, tetElem);
  const int vectorExtent0PerHex = stk::mesh::field_extent0_per_entity(variableSizeField, hexElem);
  EXPECT_EQ(vectorExtent0PerTet, numVectorValues);
  EXPECT_EQ(vectorExtent0PerHex, numVectorValues);

  const int vectorExtent1PerTet = stk::mesh::field_extent1_per_entity(variableSizeField, tetElem);
  const int vectorExtent1PerHex = stk::mesh::field_extent1_per_entity(variableSizeField, hexElem);
  EXPECT_EQ(vectorExtent1PerTet, 1);
  EXPECT_EQ(vectorExtent1PerHex, numCopies);

  double* vectorTetData = stk::mesh::field_data(variableSizeField, tetElem);
  for (int i = 0; i < vectorScalarsPerTet; ++i) {
    EXPECT_EQ(initialVectorValue[i], vectorTetData[i]);
  }

  double* vectorHexData = stk::mesh::field_data(variableSizeField, hexElem);
  for (int i = 0; i < vectorScalarsPerHex; ++i) {
    EXPECT_EQ(initialVectorValue[i], vectorHexData[i]);
  }
}
//ENDUseAdvancedFields

}
