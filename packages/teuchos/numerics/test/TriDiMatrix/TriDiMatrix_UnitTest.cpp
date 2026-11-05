#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_SerialTriDiMatrix.hpp"

namespace {

TEUCHOS_UNIT_TEST(SerialTriDiMatrix, shapeZeroLeavesNullPointers)
{
  typedef Teuchos::SerialTriDiMatrix<int, double> matrix_type;

  matrix_type A;
  TEST_EQUALITY_CONST(A.shape(3), 0);
  TEST_ASSERT(A.values() != 0);

  TEST_EQUALITY_CONST(A.shape(0), 0);
  TEST_EQUALITY_CONST(A.numRowsCols(), 0);
  TEST_ASSERT(A.empty());

  TEST_EQUALITY(A.values(), static_cast<double*>(0));
  TEST_EQUALITY(A.DL(), static_cast<double*>(0));
  TEST_EQUALITY(A.D(), static_cast<double*>(0));
  TEST_EQUALITY(A.DU(), static_cast<double*>(0));
  TEST_EQUALITY(A.DU2(), static_cast<double*>(0));
}

TEUCHOS_UNIT_TEST(SerialTriDiMatrix, shapeUninitializedZeroLeavesNullPointers)
{
  typedef Teuchos::SerialTriDiMatrix<int, double> matrix_type;

  matrix_type A;
  TEST_EQUALITY_CONST(A.shapeUninitialized(3), 0);
  TEST_ASSERT(A.values() != 0);

  TEST_EQUALITY_CONST(A.shapeUninitialized(0), 0);
  TEST_EQUALITY_CONST(A.numRowsCols(), 0);
  TEST_ASSERT(A.empty());

  TEST_EQUALITY(A.values(), static_cast<double*>(0));
  TEST_EQUALITY(A.DL(), static_cast<double*>(0));
  TEST_EQUALITY(A.D(), static_cast<double*>(0));
  TEST_EQUALITY(A.DU(), static_cast<double*>(0));
  TEST_EQUALITY(A.DU2(), static_cast<double*>(0));
}

TEUCHOS_UNIT_TEST(SerialTriDiMatrix, shapeSetsExpectedInternalPointers)
{
  typedef Teuchos::SerialTriDiMatrix<int, double> matrix_type;

  matrix_type A;
  TEST_EQUALITY_CONST(A.shape(3), 0);

  double* const values = A.values();
  TEST_ASSERT(values != 0);

  TEST_EQUALITY(A.DL(), values);
  TEST_EQUALITY(A.D(),  values + 2);
  TEST_EQUALITY(A.DU(), values + 5);
  TEST_EQUALITY(A.DU2(), values + 7);
}

} // namespace

int main(int argc, char* argv[])
{
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
