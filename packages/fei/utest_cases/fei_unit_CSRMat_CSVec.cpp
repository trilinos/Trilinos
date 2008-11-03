
#include "fei_CSVec.hpp"
#include "fei_CSRMat.hpp"

#include "fei_iostream.hpp"

#include "fei_unit_CSRMat_CSVec.hpp"

#include <cmath>
#include <limits>

void test_FillableVec_1()
{
  FEI_COUT << "testing fei::FillableVec...";

  fei::FillableVec fv;

  if (fv.hasEntry(0)) {
    throw std::runtime_error("FillableVec failed 1");
  }

  bool test_passed = true;
  try {
    fv.getEntry(0);
    test_passed = false;
  }
  catch(...) {}

  if (test_passed == false) {
    throw std::runtime_error("FillableVec failed 2");
  }

  fv.addEntry(0, 0.0);
  fv.addEntry(1, 1.0);
  fv.putEntry(2, 2.0);
  fv.addEntry(2, 2.0);

  test_passed = true;
  try {
    double coef = fv.getEntry(2);
    const double fei_eps = std::numeric_limits<double>::epsilon();
    if (std::abs(coef - 4.0) > fei_eps) test_passed = false;
  }
  catch(...) {test_passed = false;}

  if (test_passed == false) {
    throw std::runtime_error("FillableVec failed 3");
  }

  if (!fv.hasEntry(1)) {
    throw std::runtime_error("FillableVec failed 4");
  }

  if (fv.size() != 3) {
    throw std::runtime_error("FillableVec failed 5");
  }

  FEI_COUT << "ok"<<FEI_ENDL;
}

void test_FillableMat_1()
{
  FEI_COUT << "testing fei::FillableMat...";

  fei::FillableMat fm;

  if (fm.hasRow(0)) {
    throw std::runtime_error("FillableMat failed 1");
  }

  bool test_passed = true;
  try {
    fm.getRow(0);
    test_passed = false;
  }
  catch(...) {}

  if (test_passed == false) {
    throw std::runtime_error("FillableMat failed 2");
  }

  fm.sumInCoef(0, 0, 0.0);
  fm.sumInCoef(1, 1, 1.0);
  fm.putCoef(2, 2, 2.0);
  fm.sumInCoef(2, 2, 2.0);

  test_passed = true;
  try {
    fei::FillableVec* row = fm.getRow(2);
    if (row->size() != 1) test_passed = false;
  }
  catch(...) {test_passed = false;}

  if (test_passed == false) {
    throw std::runtime_error("FillableMat failed 3");
  }

  if (!fm.hasRow(1)) {
    throw std::runtime_error("FillableMat failed 4");
  }

  if (fm.getNumRows() != 3) {
    throw std::runtime_error("FillableMat failed 5");
  }

  FEI_COUT << "ok"<<FEI_ENDL;
}

void test_multiply_CSRMat_CSVec()
{
  FEI_COUT << "testing multiply_CSRMat_CSVec...";

  fei::FillableMat fm;

  fm.putCoef(0, 0, 1.0);
  fm.putCoef(0, 1, 1.0);
  fm.putCoef(1, 0, 1.0);
  fm.putCoef(1, 1, 1.0);
  fm.putCoef(2, 0, 1.0);
  fm.putCoef(2, 1, 1.0);
  fm.putCoef(2, 2, 1.0);

  fei::FillableVec fv;

  fv.putEntry(0, 1.0);
  fv.putEntry(1, 1.0);
  fv.putEntry(2, 1.0);

  fei::CSVec x(fv), y;
  fei::CSRMat A(fm);

  fei::multiply_CSRMat_CSVec(A, x, y);

  if (y.size() != 3) {
    throw std::runtime_error("CSRMat y=Ax test 1 failed.");
  }

  std::vector<int>& y_ind = y.indices();
  std::vector<double>& y_coef = y.coefs();

  if (y_ind[1] != 1) {
    throw std::runtime_error("CSRMat y=Ax test 2 failed.");
  }

  if (std::abs(y_coef[1] - 2.0) > 1.e-15) {
    throw std::runtime_error("CSRMat y=Ax test 3 failed.");
  }

  if (std::abs(y_coef[2] - 3.0) > 1.e-15) {
    throw std::runtime_error("CSRMat y=Ax test 4 failed.");
  }

  FEI_COUT << "ok"<<FEI_ENDL;
  FEI_COUT << "testing multiply_trans_CSRMat_CSVec...";

  fei::multiply_trans_CSRMat_CSVec(A, x, y);

  if (y.size() != 3) {
    throw std::runtime_error("CSRMat y=A^Tx test 1 failed.");
  }

  std::vector<int>& y_ind2 = y.indices();
  std::vector<double>& y_coef2 = y.coefs();

  if (y_ind2[1] != 1) {
    throw std::runtime_error("CSRMat y=A^Tx test 2 failed.");
  }

  if (std::abs(y_coef2[1] - 3.0) > 1.e-15) {
    throw std::runtime_error("CSRMat y=A^Tx test 3 failed.");
  }

  if (std::abs(y_coef2[2] - 1.0) > 1.e-15) {
    throw std::runtime_error("CSRMat y=A^Tx test 4 failed.");
  }

  FEI_COUT << "ok"<<FEI_ENDL;
}

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
    throw std::runtime_error("CSRMAT C=AB test 1 failed.");
  }

  std::vector<int>& cols = C.getGraph().packedColumnIndices;
  std::vector<double>& coefs = C.getPackedCoefs();

  if (cols.size() != 7) {
    throw std::runtime_error("CSRMAT C=AB test 2 failed.");
  }

  if (cols[3] != 1) {
    throw std::runtime_error("CSRMAT C=AB test 3 failed.");
  }

  if (std::abs(coefs[3] - 2.0) > 1.e-14) {
    throw std::runtime_error("CSRMAT C=AB test 4 failed.");
  }

  if (std::abs(coefs[4] - 3.0) > 1.e-14) {
    throw std::runtime_error("CSRMAT C=AB test 5 failed.");
  }

  if (std::abs(coefs[6] - 1.0) > 1.e-14) {
    throw std::runtime_error("CSRMAT C=AB test 6 failed.");
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
    throw std::runtime_error("CSRMAT C=AB test 7 failed.");
  }

  std::vector<int>& cols2 = C2.getGraph().packedColumnIndices;
  std::vector<double>& coefs2 = C2.getPackedCoefs();

  if (cols2.size() != 7) {
    throw std::runtime_error("CSRMAT C=AB test 8 failed.");
  }

  if (cols2[3] != 2) {
    throw std::runtime_error("CSRMAT C=AB test 9 failed.");
  }

  if (std::abs(coefs2[3] - 2.0) > 1.e-14) {
    throw std::runtime_error("CSRMAT C=AB test 10 failed.");
  }

  if (std::abs(coefs2[4] - 3.0) > 1.e-14) {
    throw std::runtime_error("CSRMAT C=AB test 11 failed.");
  }

  if (std::abs(coefs2[6] - 1.0) > 1.e-14) {
    throw std::runtime_error("CSRMAT C=AB test 12 failed.");
  }

  fa2.clear();

  if (fa2.getNumRows() != 0 || fa2.begin() != fa2.end()) {
    throw std::runtime_error("FillableMat::clear() test failed.");
  }

  FEI_COUT << "ok"<<FEI_ENDL;
  FEI_COUT << "testing multiply_trans_CSRMat_CSRMat...";

  fei::multiply_trans_CSRMat_CSRMat(A, B, C);

  if (C.getNumRows() != 3) {
    throw std::runtime_error("CSRMAT C=A^TB test 1 failed.");
  }

  std::vector<int>& tcols = C.getGraph().packedColumnIndices;
  std::vector<double>& tcoefs = C.getPackedCoefs();

  if (tcols.size() != 9) {
    throw std::runtime_error("CSRMAT C=A^TB test 2 failed.");
  }

  if (tcols[2] != 2) {
    throw std::runtime_error("CSRMAT C=A^TB test 3 failed.");
  }

  if (std::abs(tcoefs[2] - 1.0) > 1.e-14) {
    throw std::runtime_error("CSRMAT C=A^TB test 4 failed.");
  }

  if (std::abs(tcoefs[3] - 3.0) > 1.e-14) {
    throw std::runtime_error("CSRMAT C=A^TB test 5 failed.");
  }

  if (std::abs(tcoefs[8] - 1.0) > 1.e-14) {
    throw std::runtime_error("CSRMAT C=A^TB test 6 failed.");
  }

  FEI_COUT << "ok"<<FEI_ENDL;
}

void test_csvec_add_entry()
{
  FEI_COUT << "testing fei::add_entry(CSVec& ...)...";

  std::vector<int> ind(3);
  std::vector<double> coef(3);

  ind[0] = 2; ind[1] = 5; ind[2] = 8;
  coef[0] = 2.0; coef[1] = 5.0; coef[2] = 8.0;

  fei::CSVec csv;

  for(int i=ind.size()-1; i>=0; --i) {
    fei::add_entry(csv, ind[i], coef[i]);
  }

  if (csv.indices() != ind) {
    throw std::runtime_error("add_entry(CSVec... failed 1.");
  }

  if (csv.coefs() != coef) {
    throw std::runtime_error("add_entry(CSVec... failed 2.");
  }

  FEI_COUT << "ok" << FEI_ENDL;
  FEI_COUT << "testing fei::put_entry(CSVec& ...)...";

  coef[1] = 7.0;
  fei::put_entry(csv, ind[1], 7.0);

  if (csv.indices() != ind) {
    throw std::runtime_error("put_entry(CSVec... failed 1.");
  }

  if (csv.coefs() != coef) {
    throw std::runtime_error("put_entry(CSVec... failed 2.");
  }

  FEI_COUT << "ok" << FEI_ENDL;
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
    throw std::runtime_error("CSVec ctor test failed.");
  }

  fei::FillableVec::iterator iter = fv.begin(), iter_end = fv.end();
  unsigned i=0;
  for(; iter != iter_end; ++iter, ++i) {
    if (inds[i] != iter->first) {
      throw std::runtime_error("CSVec ctor test 2 failed.");
    }
    if (coefs[i] != iter->second) {
      throw std::runtime_error("CSVec ctor test 3 failed.");
    }
  }


  fei::FillableMat fm;

  fm.sumInCoef(0, 0, 0.0);
  fm.sumInCoef(1, 1, 1.0);
  fm.sumInCoef(2, 2, 2.0);
  fm.sumInCoef(3, 3, 3.0);

  fei::CSRMat csrm(fm);

  if (csrm.getNumRows() != 4) {
    throw std::runtime_error("CSRMat ctor test failed.");
  }

  if (csrm.getGraph().packedColumnIndices.size() != 4) {
    throw std::runtime_error("CSRMat ctor test 2 failed.");
  }

  if (csrm.getPackedCoefs()[2] != 2.0) {
    throw std::runtime_error("CSRMat ctor test 3 failed.");
  }

  FEI_COUT << "ok" << FEI_ENDL;

  test_FillableVec_1();

  test_FillableMat_1();

  test_multiply_CSRMat_CSVec();

  test_multiply_CSRMat_CSRMat();

  test_csvec_add_entry();

  return true;
}

