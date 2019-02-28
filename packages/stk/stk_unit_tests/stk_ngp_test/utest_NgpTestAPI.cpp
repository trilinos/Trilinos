#include <stk_ngp_test/ngp_test.hpp>

namespace test = ::ngp_testing;

class NGP_TEST_F_macro : public ::ngp_testing::Test { };

NGP_TEST(NGP_TEST_macro, can_compile) {}
NGP_TEST_F(NGP_TEST_F_macro, can_compile) {}

NGP_TEST(ResetReportSize, set_max_failure_reports_per_test) {
  int initialNumReportsPerTest = test::get_max_failure_reports_per_test();

  const int reportsPerTest = 10 + initialNumReportsPerTest;
  test::set_max_failure_reports_per_test(reportsPerTest);

  EXPECT_EQ(reportsPerTest, test::get_max_failure_reports_per_test());
}



