#include <gtest/gtest.h>
#include <unistd.h>                     // for unlink
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_io/DatabasePurpose.hpp>   // for DatabasePurpose::READ_MESH, etc
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <stk_mesh/base/BulkData.hpp>   // for MetaData, put_field
#include <stk_topology/topology.hpp>    // for topology, etc
#include <stk_unit_test_utils/MeshFileFixture.hpp>
#include <stk_unit_test_utils/meshCreationHelpers.hpp>
#include <Ioss_SubSystem.h>

namespace
{

//-BEGIN
TEST(ExodusFileWithVariables, queryingFileWithTimeSteps)
{
  std::shared_ptr<stk::mesh::BulkData> meshPtr = stk::mesh::MeshBuilder(MPI_COMM_WORLD).create();
  stk::io::StkMeshIoBroker stkIo;
  stkIo.set_bulk_data(meshPtr);
  std::string fileName("generated:2x3x4|variables:element,2|times:1");
  stk::io::fill_mesh_with_fields(fileName, stkIo, *meshPtr);

  const int expectedNumTimeSteps = 1;
  EXPECT_EQ(expectedNumTimeSteps, stkIo.get_num_time_steps());

  int validTimeStep = 1;
  stkIo.read_defined_input_fields(validTimeStep);

  int invalidTimeStep = 3;
  EXPECT_ANY_THROW(stkIo.read_defined_input_fields(invalidTimeStep));
}
//-END

}
