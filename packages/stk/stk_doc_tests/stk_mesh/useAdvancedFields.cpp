// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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
#include <stk_mesh/base/CoordinateSystems.hpp>  // for Cartesian, etc
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field, etc
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/FieldBase.hpp"  // for field_scalars_per_entity, etc
#include "stk_mesh/base/Types.hpp"      // for EntityId
#include "stk_topology/topology.hpp"    // for topology, etc
namespace stk { namespace mesh { class Part; } }

namespace {

//BEGIN
TEST(stkMeshHowTo, useAdvancedFields)
{
    const unsigned spatialDimension = 3;
    stk::mesh::MetaData metaData(spatialDimension, stk::mesh::entity_rank_names());

    typedef stk::mesh::Field<double, stk::mesh::Cartesian> VectorField;
    typedef stk::mesh::Field<double, stk::mesh::FullTensor36> TensorField;
    TensorField& tensorField = metaData.declare_field<TensorField>(stk::topology::ELEM_RANK, "tensor");
    VectorField& variableSizeField = metaData.declare_field<VectorField>(stk::topology::ELEM_RANK, "variableSizeField");

    stk::mesh::Part &tetPart = metaData.declare_part_with_topology("tetElementPart", stk::topology::TET_4);
    stk::mesh::Part &hexPart = metaData.declare_part_with_topology("hexElementPart", stk::topology::HEX_8);

    double initialTensorValue[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    stk::mesh::put_field_on_entire_mesh_with_initial_value(tensorField, initialTensorValue);

    double initialVectorValue[] = {1, 2, 3, 4, 5, 6, 7, 8};
    const unsigned nodesPerTet = 4;
    stk::mesh::put_field(variableSizeField, tetPart, nodesPerTet, initialVectorValue);
    const unsigned nodesPerHex = 8;
    stk::mesh::put_field(variableSizeField, hexPart, nodesPerHex, initialVectorValue);

    metaData.commit();
    stk::mesh::BulkData mesh(metaData, MPI_COMM_WORLD);
    mesh.modification_begin();
    stk::mesh::EntityId tetId = 1;
    stk::mesh::EntityIdVector tetNodes {1, 2, 3, 4};
    stk::mesh::Entity tetElem=stk::mesh::declare_element(mesh, tetPart, tetId, tetNodes);
    stk::mesh::EntityId hexId = 2;
    stk::mesh::EntityIdVector hexNodes {5, 6, 7, 8, 9, 10, 11, 12};
    stk::mesh::Entity hexElem=stk::mesh::declare_element(mesh, hexPart, hexId, hexNodes);
    mesh.modification_end();

    const unsigned tensor_scalars_per_hex = stk::mesh::field_scalars_per_entity(tensorField, hexElem);
    const unsigned tensor_scalars_per_tet = stk::mesh::field_scalars_per_entity(tensorField, tetElem);

    EXPECT_EQ(tensor_scalars_per_hex, tensor_scalars_per_tet);
    const unsigned tensor_enum_size = stk::mesh::FullTensor36::Size;
    EXPECT_EQ(tensor_scalars_per_hex, tensor_enum_size);

    double* tensorData = stk::mesh::field_data(tensorField, hexElem);
    for(unsigned i=0; i<tensor_scalars_per_hex; i++)
    {
        EXPECT_EQ(initialTensorValue[i], tensorData[i]);
    }

    const unsigned scalars_per_tet = stk::mesh::field_scalars_per_entity(variableSizeField, tetElem);
    EXPECT_EQ(nodesPerTet, scalars_per_tet);

    const unsigned scalars_per_hex = stk::mesh::field_scalars_per_entity(variableSizeField, hexElem);
    EXPECT_EQ(nodesPerHex, scalars_per_hex);

    double* vectorHexData = stk::mesh::field_data(variableSizeField, hexElem);
    for(unsigned i=0; i<scalars_per_hex; i++)
    {
        EXPECT_EQ(initialVectorValue[i], vectorHexData[i]);
    }

    double* vectorTetData = stk::mesh::field_data(variableSizeField, tetElem);
    for(unsigned i=0; i<scalars_per_tet; i++)
    {
        EXPECT_EQ(initialVectorValue[i], vectorTetData[i]);
    }
}
//END

}
