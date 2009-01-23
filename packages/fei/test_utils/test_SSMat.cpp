/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <test_utils/test_SSMat.hpp>
#include <feiArray.hpp>
#include <fei_defs.h>

#include <fei_ProcEqns.hpp>
#include <fei_EqnBuffer.hpp>
#include <fei_SSVec.hpp>
#include <fei_SSGraph.hpp>
#include <fei_SSMat.hpp>

#undef fei_file
#define fei_file "test_SSMat.cpp"
#include <fei_ErrMacros.hpp>

test_SSMat::test_SSMat(MPI_Comm comm)
 : tester(comm)
{
}

test_SSMat::~test_SSMat()
{
}

int test_SSMat::runtests()
{
  CHK_ERR( test3() );
  CHK_ERR( test4() );
  return(0);
}

int test_SSMat::test1()
{
  return(0);
}

int test_SSMat::test2()
{
  return(0);
}

int test_SSMat::test3()
{
  FEI_COUT << "testing SSMat putRow...";

  SSMat mat0, mat1;
  int numrows = 5;
  feiArray<int> indices(3);
  feiArray<double> coefs(3);
  int numcontributions=3;

  for(int i=0; i<numrows; ++i) {
    int row = i;
    for(int n=0; n<numcontributions; ++n) {
      for(int j=0; j<3; ++j) {
	indices[j] = n+j;
	coefs[j] = 1.0*(n+j+1);
      }

      mat0.sumInRow(row, indices.dataPtr(), coefs.dataPtr(), indices.length());

      for(int k=0; k<indices.length(); ++k) {
	mat1.sumInCoef(row, indices[k], coefs[k]);
      }
    }
  }

  if (mat0 != mat1) {
    throw std::runtime_error("sumInRow causes different results than sumInCoef");
  }

  FEI_COUT << "ok"<<FEI_ENDL;

  return(0);
}

int test_SSMat::test4()
{
  SSMat A, AT, B, C1, C2;

  FEI_COUT << "testing SSMat::matMat and SSMat::matTransMat...";

  int n = 5;

  for(int i=0; i<n; ++i) {
    for(int j=0; j<n; ++j) {
      double coef = (i-j)*1.0;
      A.putCoef(i,j,coef);
      AT.putCoef(j,i,coef);
      B.putCoef(i,j,coef);
    }
  }

  int err1 = A.matMat(B, C1);
  int err2 = AT.matTransMat(B, C2);

  if (err1 != 0 || err2 != 0 || C1 != C2) {
    FEI_CERR << "matMat/matTransMat test failed."<<FEI_ENDL;
    return(-1);
  }

  FEI_COUT << "ok"<<FEI_ENDL;

  return(0);
}
