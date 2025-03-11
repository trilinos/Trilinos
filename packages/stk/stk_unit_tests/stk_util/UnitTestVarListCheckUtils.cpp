// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off

#include <stk_util/diag/ParserVarUtil.hpp>
#include "stk_util/diag/String.hpp"
#include <gtest/gtest.h>
#include <regex>
#include <string>
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################
namespace Funit {

  TEST(RemoveEntryTest, ValidEntries) {
    std::vector<sierra::String> valid{"=", "variable", "variable=", "variables", "variables="};
    EXPECT_TRUE(stk::diag::remove_entry(valid, "="));
    EXPECT_TRUE(stk::diag::remove_entry(valid, "variable"));
    EXPECT_TRUE(stk::diag::remove_entry(valid, "variable="));
    EXPECT_TRUE(stk::diag::remove_entry(valid, "variables"));
    EXPECT_TRUE(stk::diag::remove_entry(valid, "variables="));
}

TEST(RemoveEntryTest, InvalidEntries) {
    std::vector<sierra::String> valid{"=", "variable", "variable=", "variables", "variables="};
    EXPECT_FALSE(stk::diag::remove_entry(valid, "invalid"));
    EXPECT_FALSE(stk::diag::remove_entry(valid, "var"));
    EXPECT_FALSE(stk::diag::remove_entry(valid, "variables=="));
}

TEST(RemoveVarTest, EmptyInput) {
    std::vector<sierra::String> input;
    std::vector<sierra::String> expected;
    EXPECT_EQ(stk::diag::remove_var_eq(input), expected);
}

TEST(RemoveVarTest, NoRemoval) {
    std::vector<sierra::String> input{"foo", "bar", "baz"};
    std::vector<sierra::String> expected{"foo", "bar", "baz"};
    EXPECT_EQ(stk::diag::remove_var_eq(input), expected);
}

TEST(RemoveVarTest, RemoveFirstEntry) {
    std::vector<sierra::String> input{"variable", "bar", "baz"};
    std::vector<sierra::String> expected{"bar", "baz"};
    EXPECT_EQ(stk::diag::remove_var_eq(input), expected);
}

TEST(RemoveVarTest, RemoveFirstAndSecondEntry) {
    std::vector<sierra::String> input{"variable", "=", "baz"};
    std::vector<sierra::String> expected{"baz"};
    EXPECT_EQ(stk::diag::remove_var_eq(input), expected);
}

TEST(RemoveVarTest, RemoveFirstButNotSecondEntry) {
    std::vector<sierra::String> input{"variable", "bar", "baz"};
    std::vector<sierra::String> expected{"bar", "baz"};
    EXPECT_EQ(stk::diag::remove_var_eq(input), expected);
}

TEST(RemoveVarTest, CornerCaseFirstEntryNotRemoved) {
    std::vector<sierra::String> input{"foo", "=", "baz"};
    std::vector<sierra::String> expected{"foo", "=", "baz"};
    EXPECT_EQ(stk::diag::remove_var_eq(input), expected);
}

TEST(RemoveVarTest, CornerCaseEmptyString) {
    std::vector<sierra::String> input{""};
    std::vector<sierra::String> expected{""};
    EXPECT_EQ(stk::diag::remove_var_eq(input), expected);
}

TEST(GetVarListTest, EmptyInput) {
  std::vector<sierra::String> input;
  std::vector<sierra::String> expected;
  EXPECT_EQ(stk::diag::join_paren_vars(input), expected);
}

TEST(GetVarListTest, NoParentheses) {
  std::vector<sierra::String> input{"foo", "bar", "baz"};
  std::vector<sierra::String> expected{"foo", "bar", "baz"};
  EXPECT_EQ(stk::diag::join_paren_vars(input), expected);
}

TEST(GetVarListTest, SingleParentheses) {
  std::vector<sierra::String> input{"foo", "(x)", "baz"};
  std::vector<sierra::String> expected{"foo(x)", "baz"};
  EXPECT_EQ(stk::diag::join_paren_vars(input), expected);
}

TEST(GetVarListTest, NestedParentheses) {
  std::vector<sierra::String> input{"foo", "(xx", "(yy))"};
  EXPECT_ANY_THROW(stk::diag::join_paren_vars(input));
}

TEST(GetVarListTest, UnmatchedParentheses) {
  std::vector<sierra::String> input{"foo", "(xx", "yy"};
  std::vector<sierra::String> expected{"foo(xx,yy"};
  EXPECT_EQ(stk::diag::join_paren_vars(input), expected);
}

TEST(GetVarListTest, MultipleEntriesWithParentheses) {
  std::vector<sierra::String> input{"foo", "(xx", "1)", "qux", "(zz", "4)"};
  std::vector<sierra::String> expected{"foo(xx,1)", "qux(zz,4)"};
  EXPECT_EQ(stk::diag::join_paren_vars(input), expected);
}

TEST(GetVarListTest, SingleParenthesesSingleWord) {
  std::vector<sierra::String> input{"foo(zz)", "baz"};
  std::vector<sierra::String> expected{"foo(zz)", "baz"};
  EXPECT_EQ(stk::diag::join_paren_vars(input), expected);
}

TEST(GetVarListTest, SingleParenthesesMultipleComponents) {
  std::vector<sierra::String> input{"foo", "(x", "1)"};
  std::vector<sierra::String> expected{"foo(x,1)"};
  EXPECT_EQ(stk::diag::join_paren_vars(input), expected);
}

TEST(GetVarListTest, SingleParenthesesEveryEntryBreak) {
  std::vector<sierra::String> input{"foo", "(", "xz", "1", ")"};
  std::vector<sierra::String> expected{"foo(xz,1)"};
  EXPECT_EQ(stk::diag::join_paren_vars(input), expected);
}
TEST(GetVarListTest, SingleParenthesesTrailingEntry) {
  std::vector<sierra::String> input{"foo", "(xy", "2)a", "bat", "rat"};
  std::vector<sierra::String> expected{"foo(xy,2)a,bat,rat"};
  EXPECT_EQ(stk::diag::join_paren_vars(input), expected);
}

TEST(StripEntryTest, RemovesSingleValidPrefix) {
  std::vector<sierra::String> valid = {sierra::String("Hello")};
  sierra::String input("HelloWorld");
  std::string result = stk::diag::strip_entry(valid, input);
  EXPECT_EQ(result, "World");
}

TEST(StripEntryTest, RemovesMultipleValidPrefixes) {
  std::vector<sierra::String> valid = {sierra::String("Hello"), sierra::String("Hi")};
  sierra::String input("HelloWorld");
  std::string result = stk::diag::strip_entry(valid, input);
  EXPECT_EQ(result, "World");
}

TEST(StripEntryTest, NoValidPrefix) {
  std::vector<sierra::String> valid = {sierra::String("Goodbye")};
  sierra::String input("HelloWorld");
  std::string result = stk::diag::strip_entry(valid, input);
  EXPECT_EQ(result, "HelloWorld");
}

TEST(StripEntryTest, EmptyInputString) {
  std::vector<sierra::String> valid = {sierra::String("Hello")};
  sierra::String input("");
  std::string result = stk::diag::strip_entry(valid, input);
  EXPECT_EQ(result, "");
}

TEST(StripEntryTest, EmptyValidPrefixes) {
  std::vector<sierra::String> valid = {};
  sierra::String input("HelloWorld");
  std::string result = stk::diag::strip_entry(valid, input);
  EXPECT_EQ(result, "HelloWorld");
}

TEST(StripEntryTest, ValidPrefixAtEnd) {
  std::vector<sierra::String> valid = {sierra::String("World")};
  sierra::String input("HelloWorld");
  std::string result = stk::diag::strip_entry(valid, input);
  EXPECT_EQ(result, "HelloWorld");
}

TEST(StripEntryTest, MultipleRemovals) {
  std::vector<sierra::String> valid = {sierra::String("Hello"), sierra::String("He")};
  sierra::String input("HelloHelloWorld");
  std::string result = stk::diag::strip_entry(valid, input);
  EXPECT_EQ(result, "HelloWorld");
}

}