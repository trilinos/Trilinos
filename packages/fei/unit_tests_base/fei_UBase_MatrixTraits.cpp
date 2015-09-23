
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <fei_MatrixTraits.hpp>
#include <fei_MatrixTraits_FillableMat.hpp>

#include <fei_iostream.hpp>

#include <cmath>
#include <limits>

TEUCHOS_UNIT_TEST(MatrixTraits, FillableMat_1)
{
  fei::FillableMat fm;

  fm.putCoef(0, 0, 0.0);
  fm.putCoef(1, 1, 1.1);
  fm.putCoef(1, 2, 1.2);
  fm.putCoef(2, 2, 2.2);

  int numRows = 0;
  fei::MatrixTraits<fei::FillableMat>::getNumLocalRows(&fm, numRows);
  TEUCHOS_TEST_EQUALITY(numRows, 3, out, success);

  std::vector<int> indices1(2);
  std::vector<double> coefs1(2);

  fei::MatrixTraits<fei::FillableMat>::copyOutRow(&fm, 1, 2, &coefs1[0], &indices1[0]);

  TEUCHOS_TEST_EQUALITY(indices1[0], 1, out, success);
  TEUCHOS_TEST_EQUALITY(indices1[1], 2, out, success);

  const double eps = std::numeric_limits<double>::epsilon();

  TEUCHOS_TEST_EQUALITY(std::abs(coefs1[0] - 1.1) < eps, true, out, success);
  TEUCHOS_TEST_EQUALITY(std::abs(coefs1[1] - 1.2) < eps, true, out, success);
}

