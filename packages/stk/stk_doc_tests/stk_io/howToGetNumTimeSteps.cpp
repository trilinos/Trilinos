#include <gtest/gtest.h>
#include <unistd.h>                     // for unlink
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_io/DatabasePurpose.hpp>   // for DatabasePurpose::READ_MESH, etc
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <stk_topology/topology.hpp>    // for topology, etc
#include <stk_unit_test_utils/MeshFileFixture.hpp>
#include <stk_unit_test_utils/meshCreationHelpers.hpp>


#include <Ioss_SubSystem.h>


namespace
{

void create_mesh_with_single_time_step(const std::string & filename, MPI_Comm communicator)
{
    stk::unit_test_util::create_mesh_with__field_1__field_2__field_3(filename, communicator);
}

class ExodusFileWithVariables : public stk::unit_test_util::MeshFileFixture { };

//-BEGIN
TEST_F(ExodusFileWithVariables, queryingFileWithSingleTimeStep_NumTimeStepsEqualsOne)
{
    create_mesh_with_single_time_step(filename, get_comm());
    read_mesh(filename);
    EXPECT_EQ(1, stkIo.get_num_time_steps());
}

TEST_F(ExodusFileWithVariables, queryingFileWithoutTimeSteps_NumTimeStepsEqualsZero)
{
    stk::unit_test_util::create_mesh_without_time_steps(filename, get_comm());
    read_mesh(filename);
    EXPECT_EQ(0, stkIo.get_num_time_steps());
}

TEST_F(ExodusFileWithVariables, readDefinedInputFieldsFromInvalidTimeStep_throws)
{
    create_mesh_with_single_time_step(filename, get_comm());
    read_mesh(filename);
    EXPECT_THROW(stkIo.read_defined_input_fields(3), std::exception);
}

TEST_F(ExodusFileWithVariables, readDefinedInputFields_throws)
{
    stk::unit_test_util::create_mesh_without_time_steps(filename, get_comm());
    read_mesh(filename);
    EXPECT_THROW(stkIo.read_defined_input_fields(1), std::exception);
}
//-END

}
