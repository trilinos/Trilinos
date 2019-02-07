#include <stk_ngp_test/Reporter.hpp>
#include <stk_ngp_test/ngp_test.hpp>

using ngp_testing::Report;

namespace test = ::ngp_testing;

NGP_TEST(Report, truncates_strings_over_max_size) {
  const char longString[] = "this_is_a_very_long_string_that_should_be_truncated_to_the_max_length";
  const char shortString[] = "...very_long_string_that_should_be_truncated_to_the_max_length";

  const int maxNumChar = test::TruncatedString::maxNumChar;
  const int longLength = sizeof(longString);
  const int shortLength = sizeof(shortString);

  EXPECT_GE(longLength, maxNumChar);
  EXPECT_EQ(shortLength, maxNumChar);

  Report report(longString, longString);
  EXPECT_STREQ(shortString, report.condition);
  EXPECT_STREQ(shortString, report.location);
}



