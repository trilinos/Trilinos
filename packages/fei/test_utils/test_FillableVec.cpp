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
#include "test_utils/test_FillableVec.hpp"
#include "fei_Exception.hpp"
#include "fei_FillableVec.hpp"

#undef fei_file
#define fei_file "test_FillableVec.cpp"

#include "fei_ErrMacros.hpp"

test_FillableVec::test_FillableVec(MPI_Comm comm)
 : tester(comm)
{
}

test_FillableVec::~test_FillableVec()
{
}

void test_FillableVec_1()
{
  FEI_COUT << "testing fei::FillableVec...";

  fei::FillableVec fv;

  if (fv.hasEntry(0)) {
    throw fei::Exception("FillableVec failed 1");
  }

  bool test_passed = true;
  try {
    fv.getEntry(0);
    test_passed = false;
  }
  catch(...) {}

  if (test_passed == false) {
    throw fei::Exception("FillableVec failed 2");
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
    throw fei::Exception("FillableVec failed 3");
  }

  if (!fv.hasEntry(1)) {
    throw fei::Exception("FillableVec failed 4");
  }

  if (fv.size() != 3) {
    throw fei::Exception("FillableVec failed 5");
  }

  FEI_COUT << "ok"<<FEI_ENDL;
}

int test_FillableVec::runtests()
{
  if (numProcs_ > 1) return(0);

  test_FillableVec_1();

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

int test_FillableVec::test1()
{
  return(0);
}

int test_FillableVec::test2()
{
  return(0);
}

int test_FillableVec::test3()
{
  return(0);
}

int test_FillableVec::test4()
{

  return(0);
}

int test_FillableVec::test5()
{

  return(0);
}

int test_FillableVec::test6()
{

  return(0);
}

int test_FillableVec::test7()
{

  return(0);
}

int test_FillableVec::test8()
{

  return(0);
}
