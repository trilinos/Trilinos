
#include <unit_tests/stk_utest_macros.hpp>

#include <stk_util/use_cases/UseCaseEnvironment.hpp>

#include <mpi.h>

#include <stk_search/IdentProc.hpp>
#include <stk_search/diag/IdentProc.hpp>

#include <iostream>

namespace stk_search_identproc_test {

void tests() {

  stk::search::ident::IdentProc<uint64_t,unsigned> a(1,0), b;
  b = a;
  stk::search::ident::IdentProc<uint64_t,unsigned> c(a), d(1,1), e(0,0);


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

  std::cout << "Test Insertion operator: " << a << std::endl;

  use_case::dw() << "Test diag writer for IdentProc: " << a << std::endl;

}

} // namespace stk_search_identproc_test

STKUNIT_UNIT_TEST(UnitTestingOfIdentProc, testUnit)
{
  MPI_Barrier( MPI_COMM_WORLD );
  stk_search_identproc_test::tests();
}
