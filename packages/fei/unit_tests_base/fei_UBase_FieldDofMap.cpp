
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <fei_iostream.hpp>
#include <fei_FieldDofMap.hpp>

#include <vector>
#include <cmath>

namespace {

TEUCHOS_UNIT_TEST(FieldDofMap, test1)
{
  fei::FieldDofMap<int> fdmap;

  fdmap.add_field(0, 1);
  fdmap.add_field(2, 2);
  fdmap.add_field(5, 5);

  int dof_id = fdmap.get_dof_id(0, 0);
  int dof_id_correct = fei::UNKNOWN;

  TEUCHOS_TEST_EQUALITY(dof_id, dof_id_correct, out, success);

  dof_id = fdmap.get_dof_id(2, 1);
  dof_id_correct = fei::UNKNOWN + 2;

  TEUCHOS_TEST_EQUALITY(dof_id, dof_id_correct, out, success);

  dof_id = fdmap.get_dof_id(5, 3);
  dof_id_correct = fei::UNKNOWN + 6;

  TEUCHOS_TEST_EQUALITY(dof_id, dof_id_correct, out, success);

  TEUCHOS_TEST_THROW(fdmap.get_dof_id(0, 1), std::runtime_error, out, success);
}

}//namespace <anonymous>

