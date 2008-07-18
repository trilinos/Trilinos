
#include "fei_CSVec.hpp"
#include "fei_CSRMat.hpp"

#include "fei_iostream.hpp"
#include "fei_Exception.hpp"

#include "fei_unit_CSRMat_CSVec.hpp"

#include <cmath>

void test_multiply_CSRMat_CSRMat()
{
  FEI_COUT << "testing multiply_CSRMat_CSRMat...";

  fei::FillableMat fa, fb;

  fa.putCoef(0, 0, 1.0);
  fa.putCoef(0, 1, 1.0);
  fa.putCoef(1, 0, 1.0);
  fa.putCoef(1, 1, 1.0);
  fa.putCoef(2, 0, 1.0);
  fa.putCoef(2, 1, 1.0);
  fa.putCoef(2, 2, 1.0);

  fb.putCoef(0, 0, 1.0);
  fb.putCoef(0, 1, 1.0);
  fb.putCoef(1, 0, 1.0);
  fb.putCoef(1, 1, 1.0);
  fb.putCoef(2, 0, 1.0);
  fb.putCoef(2, 1, 1.0);
  fb.putCoef(2, 2, 1.0);

  fei::CSRMat A(fa), B(fb), C;

  fei::multiply_CSRMat_CSRMat(A, B, C);

  if (C.getNumRows() != 3) {
    throw fei::Exception("CSRMAT C=AB test 1 failed.");
  }

  std::vector<int>& cols = C.getGraph().packedColumnIndices;
  std::vector<double>& coefs = C.getPackedCoefs();

  if (cols.size() != 7) {
    throw fei::Exception("CSRMAT C=AB test 2 failed.");
  }

  if (cols[3] != 1) {
    throw fei::Exception("CSRMAT C=AB test 3 failed.");
  }

  if (std::abs(coefs[3] - 2.0) > 1.e-14) {
    throw fei::Exception("CSRMAT C=AB test 4 failed.");
  }

  if (std::abs(coefs[4] - 3.0) > 1.e-14) {
    throw fei::Exception("CSRMAT C=AB test 5 failed.");
  }

  if (std::abs(coefs[6] - 1.0) > 1.e-14) {
    throw fei::Exception("CSRMAT C=AB test 6 failed.");
  }

  fei::FillableMat fa2, fb2;

  fa2.putCoef(1, 1, 1.0);
  fa2.putCoef(1, 2, 1.0);
  fa2.putCoef(2, 1, 1.0);
  fa2.putCoef(2, 2, 1.0);
  fa2.putCoef(3, 1, 1.0);
  fa2.putCoef(3, 2, 1.0);
  fa2.putCoef(3, 3, 1.0);

  fb2.putCoef(1, 1, 1.0);
  fb2.putCoef(1, 2, 1.0);
  fb2.putCoef(2, 1, 1.0);
  fb2.putCoef(2, 2, 1.0);
  fb2.putCoef(3, 1, 1.0);
  fb2.putCoef(3, 2, 1.0);
  fb2.putCoef(3, 3, 1.0);

  fei::CSRMat A2(fa2), B2(fb2), C2;

  fei::multiply_CSRMat_CSRMat(A2, B2, C2);

  if (C.getNumRows() != 3) {
    throw fei::Exception("CSRMAT C=AB test 7 failed.");
  }

  std::vector<int>& cols2 = C2.getGraph().packedColumnIndices;
  std::vector<double>& coefs2 = C2.getPackedCoefs();

  if (cols2.size() != 7) {
    throw fei::Exception("CSRMAT C=AB test 8 failed.");
  }

  if (cols2[3] != 2) {
    throw fei::Exception("CSRMAT C=AB test 9 failed.");
  }

  if (std::abs(coefs2[3] - 2.0) > 1.e-14) {
    throw fei::Exception("CSRMAT C=AB test 10 failed.");
  }

  if (std::abs(coefs2[4] - 3.0) > 1.e-14) {
    throw fei::Exception("CSRMAT C=AB test 11 failed.");
  }

  if (std::abs(coefs2[6] - 1.0) > 1.e-14) {
    throw fei::Exception("CSRMAT C=AB test 12 failed.");
  }

  fa2.clear();

  if (fa2.getNumRows() != 0 || fa2.begin() != fa2.end()) {
    throw fei::Exception("FillableMat::clear() test failed.");
  }

  FEI_COUT << "ok"<<FEI_ENDL;
}

bool test_csvec::run(MPI_Comm comm)
{
  FEI_COUT <<"testing CSRMat,CSVec constructors...";

  fei::FillableVec fv;

  fv.putEntry(0, 0.0);
  fv.putEntry(1, 1.0);
  fv.addEntry(2, 2.0);

  fei::CSVec csv(fv);

  std::vector<int>& inds = csv.indices();
  std::vector<double>& coefs = csv.coefs();

  if (inds.size() != fv.size()) {
    throw fei::Exception("CSVec ctor test failed.");
  }

  fei::FillableVec::iterator iter = fv.begin(), iter_end = fv.end();
  unsigned i=0;
  for(; iter != iter_end; ++iter, ++i) {
    if (inds[i] != iter->first) {
      throw fei::Exception("CSVec ctor test 2 failed.");
    }
    if (coefs[i] != iter->second) {
      throw fei::Exception("CSVec ctor test 3 failed.");
    }
  }

  FEI_COUT << "ok" << FEI_ENDL;

  test_multiply_CSRMat_CSRMat();

  return true;
}

