#include "gtest/gtest.h"
#include "stk_util/util/string_case_compare.hpp"

namespace {
void test_equal_case(const std::string& lhs, const std::string& rhs, bool expected_result)
{
  EXPECT_EQ(stk::equal_case(lhs, rhs), expected_result);
  EXPECT_EQ(stk::equal_case(lhs.c_str(), rhs.c_str()), expected_result);
  stk::EqualCase pred;
  EXPECT_EQ(pred(lhs, rhs), expected_result);

  stk::EqualCaseUnary unary_pred(rhs);
  EXPECT_EQ(unary_pred(lhs), expected_result);
}

void test_not_equal_case(const std::string& lhs, const std::string& rhs, bool expected_result)
{
  EXPECT_EQ(stk::not_equal_case(lhs, rhs), expected_result);
  EXPECT_EQ(stk::not_equal_case(lhs.c_str(), rhs.c_str()), expected_result);
  stk::NotEqualCase pred;
  EXPECT_EQ(pred(lhs, rhs), expected_result);

  stk::NotEqualCaseUnary unary_pred(rhs);
  EXPECT_EQ(unary_pred(lhs), expected_result);
}

void test_less_case(const std::string& lhs, const std::string& rhs, bool expected_result)
{
  EXPECT_EQ(stk::less_case(lhs, rhs), expected_result);
  EXPECT_EQ(stk::less_case(lhs.c_str(), rhs.c_str()), expected_result);
  stk::LessCase pred;
  EXPECT_EQ(pred(lhs, rhs), expected_result);
}

void test_less_equal_case(const std::string& lhs, const std::string& rhs, bool expected_result)
{
  EXPECT_EQ(stk::less_equal_case(lhs, rhs), expected_result);
  EXPECT_EQ(stk::less_equal_case(lhs.c_str(), rhs.c_str()), expected_result);
  stk::LessEqualCase pred;
  EXPECT_EQ(pred(lhs, rhs), expected_result);
}

void test_greater_case(const std::string& lhs, const std::string& rhs, bool expected_result)
{
  EXPECT_EQ(stk::greater_case(lhs, rhs), expected_result);
  EXPECT_EQ(stk::greater_case(lhs.c_str(), rhs.c_str()), expected_result);
  stk::GreaterCase pred;
  EXPECT_EQ(pred(lhs, rhs), expected_result);
}

void test_greater_equal_case(const std::string& lhs, const std::string& rhs, bool expected_result)
{
  EXPECT_EQ(stk::greater_equal_case(lhs, rhs), expected_result);
  EXPECT_EQ(stk::greater_equal_case(lhs.c_str(), rhs.c_str()), expected_result);
  stk::GreaterEqualCase pred;
  EXPECT_EQ(pred(lhs, rhs), expected_result);
}
}


TEST(StringCaseCompare, EqualCase)
{
  test_equal_case("abc", "abc", true);
  test_equal_case("abc", "Abc", true);
  test_equal_case("Abc", "abc", true);

  test_equal_case("abcde", "abc", false);
  test_equal_case("abc", "abcde", false);
}

TEST(StringCaseCompare, NotEqualCase)
{
  test_not_equal_case("abc", "abc", false);
  test_not_equal_case("abc", "Abc", false);
  test_not_equal_case("Abc", "abc", false);

  test_not_equal_case("abcde", "abc", true);
  test_not_equal_case("abc", "abcde", true);
}

TEST(StringCaseCompare, LessCase)
{
  test_less_case("ab", "abc", true);
  test_less_case("aB", "abc", true);
  test_less_case("ab", "aBc", true);

  test_less_case("abc", "ab", false);
  test_less_case("aBc", "ab", false);
}

TEST(StringCaseCompare, LessEqualCase)
{
  test_less_equal_case("ab",  "abc", true);
  test_less_equal_case("aB",  "abc", true);
  test_less_equal_case("ab",  "aBc", true);
  test_less_equal_case("abc", "aBc", true);
  test_less_equal_case("aBc", "abc", true);

  test_less_equal_case("abc", "ab", false);
  test_less_equal_case("aBc", "ab", false);
}

TEST(StringCaseCompare, GreaterCase)
{
  test_greater_case("ab", "abc", false);
  test_greater_case("aB", "abc", false);
  test_greater_case("ab", "aBc", false);

  test_greater_case("abc", "ab", true);
  test_greater_case("aBc", "ab", true);
}

TEST(StringCaseCompare, GreaterEqualCase)
{
  test_greater_equal_case("ab", "abc", false);
  test_greater_equal_case("aB", "abc", false);
  test_greater_equal_case("ab", "aBc", false);

  test_greater_equal_case("abc", "aBc", true);
  test_greater_equal_case("aBc", "abc", true);
  test_greater_equal_case("abc", "ab", true);
  test_greater_equal_case("aBc", "ab", true);
}
