
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <fei_LinearDecomposition.hpp>
#include <iostream>

namespace {

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(LinearDecomposition, test0, T)
{
  int myproc = 3;
  int numprocs = 4;
  T lowest_global = 1;
  T highest_global = 100;

  fei::LinearDecomposition<T> ld1(myproc, numprocs, lowest_global, highest_global);
  TEUCHOS_TEST_EQUALITY(ld1.which_proc(0), -1, out, success)
  TEUCHOS_TEST_EQUALITY(ld1.which_proc(5), 0, out, success)
  TEUCHOS_TEST_EQUALITY(ld1.which_proc(25), 0, out, success)
  TEUCHOS_TEST_EQUALITY(ld1.which_proc(26), 1, out, success)
  TEUCHOS_TEST_EQUALITY(ld1.which_proc(50), 1, out, success)
  TEUCHOS_TEST_EQUALITY(ld1.which_proc(51), 2, out, success)
  TEUCHOS_TEST_EQUALITY(ld1.which_proc(76), 3, out, success)
  TEUCHOS_TEST_EQUALITY(ld1.which_proc(100), 3, out, success)
  TEUCHOS_TEST_EQUALITY(ld1.which_proc(101), -1, out, success)

  lowest_global = 1;
  highest_global = 12;
  numprocs = 10;
  fei::LinearDecomposition<T> ld2(myproc, numprocs, lowest_global, highest_global);
  TEUCHOS_TEST_EQUALITY(ld2.which_proc(0), -1, out, success)
  TEUCHOS_TEST_EQUALITY(ld2.which_proc(1), 0, out, success)
  TEUCHOS_TEST_EQUALITY(ld2.which_proc(2), 0, out, success)
  TEUCHOS_TEST_EQUALITY(ld2.which_proc(3), 1, out, success)
  TEUCHOS_TEST_EQUALITY(ld2.which_proc(5), 2, out, success)
  TEUCHOS_TEST_EQUALITY(ld2.which_proc(11), 8, out, success)
  TEUCHOS_TEST_EQUALITY(ld2.which_proc(12), 9, out, success)
  TEUCHOS_TEST_EQUALITY(ld2.which_proc(13), -1, out, success)
}

#define UNIT_TEST_GROUP(TYPE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(LinearDecomposition,test0,TYPE)

typedef long int longint;
UNIT_TEST_GROUP(int)
UNIT_TEST_GROUP(longint)

}//namespace <anonymous>

