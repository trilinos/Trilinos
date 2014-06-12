
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <ostream>                      // for basic_ostream::operator<<
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, entity_rank_names, etc
#include "mpi.h"                        // for MPI_COMM_WORLD
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
