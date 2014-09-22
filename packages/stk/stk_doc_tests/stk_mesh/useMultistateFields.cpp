// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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
#include <ostream>                      // for basic_ostream::operator<<
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, entity_rank_names, etc
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include "stk_mesh/base/FieldState.hpp"  // for FieldState, etc
#include "stk_mesh/base/Types.hpp"      // for EntityId
#include "stk_topology/topology.hpp"    // for topology, etc

namespace {

//BEGIN
TEST(stkMeshHowTo, useMultistateField)
{
    const unsigned spatialDimension = 3;
    stk::mesh::MetaData metaData(spatialDimension, stk::mesh::entity_rank_names());

    typedef stk::mesh::Field<double> ScalarField;
    const unsigned numStates = 2;
    ScalarField& temperatureFieldStateNp1 = metaData.declare_field<ScalarField>(stk::topology::NODE_RANK, "temperature", numStates);

    double initialTemperatureValue = 1.0;
    stk::mesh::put_field_on_entire_mesh_with_initial_value(temperatureFieldStateNp1, &initialTemperatureValue);

    metaData.commit();
    stk::mesh::BulkData mesh(metaData, MPI_COMM_WORLD);
    mesh.modification_begin();
    stk::mesh::EntityId nodeId = 1;
    stk::mesh::Entity node = mesh.declare_entity(stk::topology::NODE_RANK, nodeId);
    mesh.modification_end();

    EXPECT_EQ(stk::mesh::StateNP1, temperatureFieldStateNp1.state());
    double* temperatureStateNp1 = stk::mesh::field_data(temperatureFieldStateNp1, node);
    EXPECT_EQ(initialTemperatureValue, *temperatureStateNp1);
    double newTemperatureValue = 2.0;
    *temperatureStateNp1 = newTemperatureValue;

    ScalarField& temperatureFieldStateN = temperatureFieldStateNp1.field_of_state(stk::mesh::StateN);
    double* temperatureStateN = stk::mesh::field_data(temperatureFieldStateN, node);
    EXPECT_EQ(initialTemperatureValue, *temperatureStateN);

    mesh.update_field_data_states();

    temperatureStateN = stk::mesh::field_data(temperatureFieldStateN, node);
    EXPECT_EQ(newTemperatureValue, *temperatureStateN);
}
//END

}
