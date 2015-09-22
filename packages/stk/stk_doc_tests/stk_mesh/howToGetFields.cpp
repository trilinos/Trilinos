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
#include <stddef.h>                     // for size_t
#include <stk_mesh/base/CoordinateSystems.hpp>  // for Cartesian3d, etc
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, entity_rank_names, etc
#include "stk_mesh/base/FieldBase.hpp"  // for field_data, etc
#include "stk_topology/topology.hpp"    // for topology, etc
namespace stk { namespace mesh { class Part; } }

namespace {

//BEGIN
TEST(stkMeshHowTo, getFields)
{
    const unsigned spatialDimension = 3;
    stk::mesh::MetaData metaData(spatialDimension, stk::mesh::entity_rank_names());

    typedef stk::mesh::Field<double> ScalarField;
    typedef stk::mesh::Field<double, stk::mesh::Cartesian3d> VectorField;
    const stk::mesh::EntityRank node_rank = stk::topology::NODE_RANK;
    const stk::mesh::EntityRank elem_rank = stk::topology::ELEM_RANK;

    struct MyFields {
        ScalarField * pressureField;
        VectorField * displacementsField;
    };
    MyFields myFields;

    const std::string pressureFieldName = "pressure";
    const std::string displacementsFieldName = "displacements";
    myFields.pressureField = &metaData.declare_field<ScalarField>(elem_rank, pressureFieldName);
    myFields.displacementsField = &metaData.declare_field<VectorField>(node_rank, displacementsFieldName);
    metaData.commit();

    // get_field requires a string-compare search through all fields of the given EntityRank

    ScalarField* pressureField = metaData.get_field<ScalarField>(elem_rank, pressureFieldName);
    ASSERT_TRUE( pressureField != NULL );
    EXPECT_EQ( myFields.pressureField->mesh_meta_data_ordinal(), pressureField->mesh_meta_data_ordinal() );

    stk::mesh::FieldBase* pressureFieldBase = metaData.get_field(elem_rank, pressureFieldName);
    ASSERT_TRUE( pressureFieldBase != NULL );
    EXPECT_EQ( myFields.pressureField->mesh_meta_data_ordinal(), pressureFieldBase->mesh_meta_data_ordinal() );

    VectorField* displacementsField = metaData.get_field<VectorField>(node_rank, displacementsFieldName);
    ASSERT_TRUE( displacementsField != NULL );
    EXPECT_EQ( myFields.displacementsField->mesh_meta_data_ordinal(), displacementsField->mesh_meta_data_ordinal() );

    stk::mesh::FieldBase* displacementsFieldBase = metaData.get_field(node_rank, displacementsFieldName);
    ASSERT_TRUE( displacementsFieldBase != NULL );
    EXPECT_EQ( myFields.displacementsField->mesh_meta_data_ordinal(), displacementsFieldBase->mesh_meta_data_ordinal() );

}
//END

}
