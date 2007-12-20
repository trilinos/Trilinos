/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <cmath>

#include <test_utils/test_CSMat.hpp>
#include "fei_FillableVec.hpp"
#include "fei_FillableMat.hpp"
#include "fei_CSVec.hpp"
#include "fei_CSRMat.hpp"

#include "fei_Exception.hpp"

#undef fei_file
#define fei_file "test_CSMat.cpp"
#include <fei_ErrMacros.hpp>

test_CSMat::test_CSMat(MPI_Comm comm)
 : tester(comm)
{
}

test_CSMat::~test_CSMat()
{
}

int test_CSMat::runtests()
{
  CHK_ERR( test1() );
  CHK_ERR( test2() );
  CHK_ERR( test3() );
  CHK_ERR( test4() );
  return(0);
}

int test_CSMat::test1()
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


  fei::FillableMat fm;

  fm.sumInCoef(0, 0, 0.0);
  fm.sumInCoef(1, 1, 1.0);
  fm.sumInCoef(2, 2, 2.0);
  fm.sumInCoef(3, 3, 3.0);

  fei::CSRMat csrm(fm);

  if (csrm.getNumRows() != 4) {
    throw fei::Exception("CSRMat ctor test failed.");
  }

  if (csrm.getGraph().packedColumnIndices.size() != 4) {
    throw fei::Exception("CSRMat ctor test 2 failed.");
  }

  if (csrm.getPackedCoefs()[2] != 2.0) {
    throw fei::Exception("CSRMat ctor test 3 failed.");
  }

  FEI_COUT << "ok"<<FEI_ENDL;
  return(0);
}

int test_CSMat::test2()
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
    throw fei::Exception("CSRMat y=Ax test 1 failed.");
  }

  std::vector<int>& y_ind = y.indices();
  std::vector<double>& y_coef = y.coefs();

  if (y_ind[1] != 1) {
    throw fei::Exception("CSRMat y=Ax test 2 failed.");
  }

  if (std::abs(y_coef[1] - 2.0) > 1.e-15) {
    throw fei::Exception("CSRMat y=Ax test 3 failed.");
  }

  if (std::abs(y_coef[2] - 3.0) > 1.e-15) {
    throw fei::Exception("CSRMat y=Ax test 4 failed.");
  }

  FEI_COUT << "ok"<<FEI_ENDL;
  FEI_COUT << "testing multiply_trans_CSRMat_CSVec...";

  fei::multiply_trans_CSRMat_CSVec(A, x, y);

  if (y.size() != 3) {
    throw fei::Exception("CSRMat y=A^Tx test 1 failed.");
  }

  std::vector<int>& y_ind2 = y.indices();
  std::vector<double>& y_coef2 = y.coefs();

  if (y_ind2[1] != 1) {
    throw fei::Exception("CSRMat y=A^Tx test 2 failed.");
  }

  if (std::abs(y_coef2[1] - 3.0) > 1.e-15) {
    throw fei::Exception("CSRMat y=A^Tx test 3 failed.");
  }

  if (std::abs(y_coef2[2] - 1.0) > 1.e-15) {
    throw fei::Exception("CSRMat y=A^Tx test 4 failed.");
  }

  FEI_COUT << "ok"<<FEI_ENDL;
  return(0);
}

int test_CSMat::test3()
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

  FEI_COUT << "ok"<<FEI_ENDL;
  FEI_COUT << "testing multiply_trans_CSRMat_CSRMat...";

  fei::multiply_trans_CSRMat_CSRMat(A, B, C);

  if (C.getNumRows() != 3) {
    throw fei::Exception("CSRMAT C=A^TB test 1 failed.");
  }

  std::vector<int>& tcols = C.getGraph().packedColumnIndices;
  std::vector<double>& tcoefs = C.getPackedCoefs();

  if (tcols.size() != 9) {
    throw fei::Exception("CSRMAT C=A^TB test 2 failed.");
  }

  if (tcols[2] != 2) {
    throw fei::Exception("CSRMAT C=A^TB test 3 failed.");
  }

  if (std::abs(tcoefs[2] - 1.0) > 1.e-14) {
    throw fei::Exception("CSRMAT C=A^TB test 4 failed.");
  }

  if (std::abs(tcoefs[3] - 3.0) > 1.e-14) {
    throw fei::Exception("CSRMAT C=A^TB test 5 failed.");
  }

  if (std::abs(tcoefs[8] - 1.0) > 1.e-14) {
    throw fei::Exception("CSRMAT C=A^TB test 6 failed.");
  }

  FEI_COUT << "ok"<<FEI_ENDL;
  return(0);
}

int test_CSMat::test4()
{
  return(0);
}
