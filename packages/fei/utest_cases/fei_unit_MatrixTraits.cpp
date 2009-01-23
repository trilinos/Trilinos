
#include <fei_MatrixTraits.hpp>
#include <fei_MatrixTraits_FillableMat.hpp>

#include <fei_iostream.hpp>

#include <fei_unit_MatrixTraits.hpp>

#include <cmath>
#include <limits>

void test_MatrixTraits_FillableMat_1()
{
  FEI_COUT << "testing fei::MatrixTraits<fei::FillableMat>...";

  fei::FillableMat fm;

  fm.putCoef(0, 0, 0.0);
  fm.putCoef(1, 1, 1.1);
  fm.putCoef(1, 2, 1.2);
  fm.putCoef(2, 2, 2.2);

  int numRows = 0;
  fei::MatrixTraits<fei::FillableMat>::getNumLocalRows(&fm, numRows);
  if (numRows != 3) {
    throw std::runtime_error("MatrixTraits<FillableMat> failed 1");
  }

  std::vector<int> indices1(2);
  std::vector<double> coefs1(2);

  fei::MatrixTraits<fei::FillableMat>::copyOutRow(&fm, 1, 2, &coefs1[0], &indices1[0]);

  if (indices1[0] != 1 || indices1[1] != 2) {
    throw std::runtime_error("MatrixTraits<FillableMat> failed 2");
  }

  const double eps = std::numeric_limits<double>::epsilon();

  if (std::abs(coefs1[0] - 1.1) > eps || std::abs(coefs1[1] - 1.2) > eps) {
    throw std::runtime_error("MatrixTraits<FillableMat> failed 3");
  }

  FEI_COUT << "ok"<<FEI_ENDL;
}

bool test_mtraits::run(MPI_Comm comm)
{
  test_MatrixTraits_FillableMat_1();

  return true;
}

