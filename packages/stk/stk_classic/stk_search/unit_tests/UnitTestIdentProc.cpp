/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/use_cases/UseCaseEnvironment.hpp>

#include <mpi.h>

#include <stk_search/IdentProc.hpp>
#include <stk_search/diag/IdentProc.hpp>

#include <iostream>

namespace stk_search_identproc_test {

void tests() {

  stk_classic::search::ident::IdentProc<uint64_t,unsigned> a(1,0), b;
  b = a;
  stk_classic::search::ident::IdentProc<uint64_t,unsigned> c(a), d(1,1), e(0,0);


  STKUNIT_ASSERT_EQUAL(a == b,true);
  STKUNIT_ASSERT_EQUAL(a == d,false);
  STKUNIT_ASSERT_EQUAL(a != d,true);
  STKUNIT_ASSERT_EQUAL(a != b,false);
  STKUNIT_ASSERT_EQUAL(a < d,true);
  STKUNIT_ASSERT_EQUAL(a < b,false);
  STKUNIT_ASSERT_EQUAL(a > e,true);
  STKUNIT_ASSERT_EQUAL(a > b,false);
  STKUNIT_ASSERT_EQUAL(a <= b,true);
  STKUNIT_ASSERT_EQUAL(a <= d,true);
  STKUNIT_ASSERT_EQUAL(a <= e,false);
  STKUNIT_ASSERT_EQUAL(a >= b,true);
  STKUNIT_ASSERT_EQUAL(a >= d,false);
  STKUNIT_ASSERT_EQUAL(a >= e,true);

  //use_case::dw() << "Test diag writer for IdentProc: " << a << std::endl;

}

} // namespace stk_search_identproc_test

STKUNIT_UNIT_TEST(UnitTestingOfIdentProc, testUnit)
{
  MPI_Barrier( MPI_COMM_WORLD );
  stk_search_identproc_test::tests();
}
