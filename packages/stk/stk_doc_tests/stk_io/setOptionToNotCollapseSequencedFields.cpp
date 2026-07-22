#include <gtest/gtest.h>
#include <unistd.h>                     // for unlink
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_io/DatabasePurpose.hpp>   // for DatabasePurpose::READ_MESH, etc
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <stk_topology/topology.hpp>    // for topology, etc
#include <stk_unit_test_utils/MeshFileFixture.hpp>
#include <stk_unit_test_utils/meshCreationHelpers.hpp>

namespace
{

class MultipleNumberedFieldsWithSameBaseName : public stk::unit_test_util::MeshFileFixture { };

//-BEGIN
TEST_F(MultipleNumberedFieldsWithSameBaseName, whenReading_collapseToSingleStkField)
{
  stk::unit_test_util::create_mesh_with__field_1__field_2__field_3(filename, get_comm());
  read_mesh(filename);
  EXPECT_EQ(1u, get_meta().get_fields(stk::topology::ELEM_RANK).size());
}

TEST_F(MultipleNumberedFieldsWithSameBaseName, whenReadingWithoutCollapseOption_threeStkFieldsAreRead)
{
  stk::unit_test_util::create_mesh_with__field_1__field_2__field_3(filename, get_comm());
  stkIo.set_option_to_not_collapse_sequenced_fields();
  read_mesh(filename);
  EXPECT_EQ(3u, get_meta().get_fields(stk::topology::ELEM_RANK).size());
}
//-END

}
