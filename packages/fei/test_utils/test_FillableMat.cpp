/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_macros.hpp"
#include "fei_mpi.h"
#include <limits>
#include <cmath>
#include "test_utils/test_FillableMat.hpp"
#include "fei_Exception.hpp"
#include "fei_FillableVec.hpp"
#include "fei_FillableMat.hpp"

#undef fei_file
#define fei_file "test_FillableMat.cpp"

#include "fei_ErrMacros.hpp"

test_FillableMat::test_FillableMat(MPI_Comm comm)
 : tester(comm)
{
}

test_FillableMat::~test_FillableMat()
{
}

void test_FillableMat_1()
{
  FEI_COUT << "testing fei::FillableMat...";

  fei::FillableMat fm;

  if (fm.hasRow(0)) {
    throw fei::Exception("FillableMat failed 1");
  }

  bool test_passed = true;
  try {
    fm.getRow(0);
    test_passed = false;
  }
  catch(...) {}

  if (test_passed == false) {
    throw fei::Exception("FillableMat failed 2");
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
    throw fei::Exception("FillableMat failed 3");
  }

  if (!fm.hasRow(1)) {
    throw fei::Exception("FillableMat failed 4");
  }

  if (fm.getNumRows() != 3) {
    throw fei::Exception("FillableMat failed 5");
  }

  FEI_COUT << "ok"<<FEI_ENDL;
}

int test_FillableMat::runtests()
{
  if (numProcs_ > 1) return(0);

  test_FillableMat_1();

  CHK_ERR( test1() );
  CHK_ERR( test2() );
  CHK_ERR( test3() );
  CHK_ERR( test4() );
  CHK_ERR( test5() );
  CHK_ERR( test6() );
  CHK_ERR( test7() );
  CHK_ERR( test8() );

  return(0);
}

int test_FillableMat::test1()
{
  return(0);
}

int test_FillableMat::test2()
{
  return(0);
}

int test_FillableMat::test3()
{
  return(0);
}

int test_FillableMat::test4()
{

  return(0);
}

int test_FillableMat::test5()
{

  return(0);
}

int test_FillableMat::test6()
{

  return(0);
}

int test_FillableMat::test7()
{

  return(0);
}

int test_FillableMat::test8()
{

  return(0);
}
