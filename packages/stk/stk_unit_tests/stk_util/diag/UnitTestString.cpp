// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include "gtest/gtest.h"
#include "stk_util/diag/String.hpp"      // for String, operator==, StringBase
#include "stk_util/diag/StringUtil.hpp"  // for case_strcmp
#include <string>                        // for string

TEST(StkString, case_insensitivity)
{
  sierra::String ABC("ABC");
  sierra::String abc("abc");
  bool strings_match_by_case_insensitivity = ABC == abc;
  EXPECT_TRUE(strings_match_by_case_insensitivity);
}

TEST(StkString, case_sensitive_comparison_with_sierra_string)
{
  sierra::String ABC("ABC");
  sierra::String abc("abc");
  bool strings_do_not_match = sierra::case_strcmp(ABC, abc);
  EXPECT_FALSE(strings_do_not_match);
}

TEST(StkString, case_sensitive_comparison_with_std_string)
{
  std::string ABC("ABC");
  std::string abc("abc");
  bool strings_do_not_match = sierra::case_strcmp(ABC, abc);
  EXPECT_FALSE(strings_do_not_match);
}

namespace
{
class TestClass
{
 public:
  TestClass() = default;

  const sierra::String str_member;
};
}  // namespace

TEST(StkString, default_constructible)
{
  sierra::String s;
  EXPECT_TRUE(s.empty());

  sierra::StringBase<sierra::char_label_traits> s_label;
  EXPECT_TRUE(s_label.empty());

  // If we change the default constructor to "StringBase() = default" this fails to compile unfortunately
  TestClass d;
  EXPECT_TRUE(d.str_member.empty());
}

TEST(StkString, char_label_traits_append)
{
  sierra::StringBase<sierra::char_label_traits> s_label;

  s_label.append("abc", 3);
  EXPECT_EQ("abc", s_label.s_str());

  s_label.append("def");
  EXPECT_EQ("abc_def", s_label.s_str());
}

TEST(StkJoinString, emptyContainer)
{
  const std::set<sierra::String> emptyVect;
  EXPECT_EQ("", stk::util::join(emptyVect, ", "));
}

TEST(StkJoinString, emptyIterRange)
{
  const std::set<sierra::String> emptyVect;
  EXPECT_EQ("", stk::util::join(emptyVect.begin(), emptyVect.end(), ", "));
}

TEST(StkJoinString, containerContents)
{
  const std::vector<sierra::String> theVect{"a", "b", "c"};
  EXPECT_EQ("a, b, c", stk::util::join(theVect, ", "));
}

TEST(StkJoinString, iterRangeContents)
{
  const std::vector<sierra::String> theVect{"a", "b", "c"};
  EXPECT_EQ("a, b, c", stk::util::join(theVect.begin(), theVect.end(), ", "));
}

TEST(StkJoinString, containerNonString)
{
  const std::vector<int> theVect{0, 1, 2};
  EXPECT_EQ("0, 1, 2", stk::util::join(theVect, ", "));
}

TEST(StkJoinString, iterRangeNonString)
{
  const std::vector<int> theVect{0, 1, 2};
  EXPECT_EQ("0, 1, 2", stk::util::join(theVect.begin(), theVect.end(), ", "));
}

TEST(StkConstructString, fromConstCharPtr)
{
  const char * const str = "test message";
  const sierra::String sierra_str{str};
  EXPECT_EQ(str, sierra_str);
}

TEST(StkConstructString, fromNullPtr)
{
  const char * const str = nullptr;
  const sierra::String sierra_str{str};
  EXPECT_TRUE(sierra_str.empty());
}

TEST(StkConstructString, fromStdString)
{
  const std::string str = "test message";
  const sierra::String sierra_str{str};
  EXPECT_EQ(str, sierra_str);
}

TEST(StkConstructString, fromStringLiteral)
{
  const std::string str = "test message";
  const sierra::String sierra_str{"test message"};
  EXPECT_EQ(str, sierra_str);
}

TEST(StkConstructString, fromTemporaryStdString)
{
  const std::string str = "test message";
  const sierra::String sierra_str{std::string("test message")};
  EXPECT_EQ(str, sierra_str);
}

TEST(StkConstructString, fromStdStringView)
{
  const std::string_view str = "test message";
  const sierra::String sierra_str{str};
  EXPECT_EQ(sierra::String(str), sierra_str);
}

TEST(StkAssignString, fromConstCharPtr)
{
  const char * const str = "test message";
  sierra::String sierra_str;
  sierra_str.assign(str);
  EXPECT_EQ(str, sierra_str);
}

TEST(StkAssignString, fromNullPtr)
{
  const char * const str = nullptr;
  sierra::String sierra_str;
  sierra_str.assign(str);
  EXPECT_TRUE(sierra_str.empty());
}

TEST(StkAssignString, fromStdString)
{
  const std::string str = "test message";
  sierra::String sierra_str;
  sierra_str.assign(str);
  EXPECT_EQ(str, sierra_str);
}

namespace {
  sierra::String func_taking_sierra_string(const sierra::String & sierra_str)
  {
    return sierra_str;
  }
}

// Note, following two require non explicit converting constructor for StringViewLike -> StringBase
// I think we should deprecate this behavior.
TEST(StkString, canCallSierraStringFunctionWithStandardString)
{
  const std::string str = "test message";

  const sierra::String sierra_str = func_taking_sierra_string(str);
  EXPECT_EQ(str, sierra_str);
}

TEST(StkString, canCallSierraStringFunctionWithStringView)
{
  const std::string_view str = "test message";

  //explicit conversion from string view
  const sierra::String sierra_str = func_taking_sierra_string(sierra::String(str));
  EXPECT_EQ(sierra::String(str), sierra_str);
}

TEST(StkAssignToString, fromStdStringView)
{
  const std::string_view str = "test message";
  sierra::String sierra_str;
  sierra_str = str;
  EXPECT_EQ(sierra::String(str), sierra_str);
}

TEST(StkAssignToString, fromStdString)
{
  const std::string str = "test message";
  sierra::String sierra_str;
  sierra_str = str;
  EXPECT_EQ(sierra::String(str), sierra_str);
}

TEST(StkAssignToString, fromConstCharPtr)
{
  const char * str = "test message";
  sierra::String sierra_str;
  sierra_str = str;
  EXPECT_EQ(sierra::String(str), sierra_str);
}

TEST(StkAssignToString, fromNullPtr)
{
  const char * str = nullptr;
  sierra::String sierra_str;
  sierra_str = str;
  EXPECT_EQ(sierra::String(), sierra_str);
}

TEST(StkAppendToString, fromStdStringView)
{
  const std::string_view str = "test message";
  sierra::String sierra_str{"first: "};
  sierra_str.append(str);
  EXPECT_EQ(sierra::String("first: test message"), sierra_str);
}

TEST(StkAppendToString, fromStdString)
{
  const std::string str = "test message";
  sierra::String sierra_str{"first: "};
  sierra_str.append(str);
  EXPECT_EQ(sierra::String("first: test message"), sierra_str);
}

TEST(StkAppendToString, fromConstCharPtr)
{
  const char * str = "test message";
  sierra::String sierra_str{"first: "};
  sierra_str.append(str);
  EXPECT_EQ(sierra::String("first: test message"), sierra_str);
}

TEST(StkOperatorPlusEqToString, fromStdStringView)
{
  const std::string_view str = "test message";
  sierra::String sierra_str{"first: "};
  sierra_str += str;
  EXPECT_EQ(sierra::String("first: test message"), sierra_str);
}

TEST(StkOperatorPlusEqToString, fromStdString)
{
  const std::string str = "test message";
  sierra::String sierra_str{"first: "};
  sierra_str += str;
  EXPECT_EQ(sierra::String("first: test message"), sierra_str);
}

TEST(StkOperatorPlusEqToString, fromConstCharPtr)
{
  const char * str = "test message";
  sierra::String sierra_str{"first: "};
  sierra_str += str;
  EXPECT_EQ(sierra::String("first: test message"), sierra_str);
}

TEST(StkStringComparison, equalStrings)
{
  sierra::String str_1{"abcdefg"};
  sierra::String str_2{"AbcdEfg"};
  EXPECT_TRUE(str_1 == str_2);
  EXPECT_TRUE(str_1 <= str_2);
  EXPECT_TRUE(str_1 >= str_2);
  EXPECT_FALSE(str_1 < str_2);
  EXPECT_FALSE(str_1 > str_2);
  EXPECT_FALSE(str_1 != str_2);
}

using SierraStd = std::pair<sierra::String, std::string>;
using StdSierra = std::pair<std::string, sierra::String>;
using CharSierra = std::pair<const char *, sierra::String>;
using SierraChar = std::pair<sierra::String, const char *>;
using StringViewSierra = std::pair<std::string_view, sierra::String>;
using SierraStringView = std::pair<sierra::String, std::string_view>;

template <typename T>
class StringTypePairFixture : public testing::Test 
{
public:
using Type = T;
};

typedef ::testing::Types<SierraStd, StdSierra, CharSierra, SierraChar, StringViewSierra, SierraStringView> MyTypes;
TYPED_TEST_SUITE(StringTypePairFixture, MyTypes,);
TYPED_TEST(StringTypePairFixture, str1LessThanStr2ByLength)
{
  typename TypeParam::first_type  str_1{"abcdef"};
  typename TypeParam::second_type str_2{"AbcdEfg"};
  EXPECT_FALSE(str_1 == str_2);
  EXPECT_TRUE(str_1 <= str_2);
  EXPECT_FALSE(str_1 >= str_2);
  EXPECT_TRUE(str_1 < str_2);
  EXPECT_FALSE(str_1 > str_2);
  EXPECT_TRUE(str_1 != str_2);
}

TYPED_TEST(StringTypePairFixture, str1LessThanStr2ByCharacter)
{
  typename TypeParam::first_type str_1{"abbdefg"};
  typename TypeParam::second_type str_2{"AbcdEfg"};
  EXPECT_FALSE(str_1 == str_2);
  EXPECT_TRUE(str_1 <= str_2);
  EXPECT_FALSE(str_1 >= str_2);
  EXPECT_TRUE(str_1 < str_2);
  EXPECT_FALSE(str_1 > str_2);
  EXPECT_TRUE(str_1 != str_2);
}


TYPED_TEST(StringTypePairFixture, str1GreaterThanStr2)
{
  typename TypeParam::first_type str_1{"AbcdEfg"};
  typename TypeParam::second_type str_2{"abcdef"};
  EXPECT_FALSE(str_1 == str_2);
  EXPECT_FALSE(str_1 <= str_2);
  EXPECT_TRUE(str_1 >= str_2);
  EXPECT_FALSE(str_1 < str_2);
  EXPECT_TRUE(str_1 > str_2);
  EXPECT_TRUE(str_1 != str_2);
}

TYPED_TEST(StringTypePairFixture, operatorPlus)
{
  typename TypeParam::first_type str_1{"oNe string, "};
  typename TypeParam::second_type str_2{"two string"};
  const sierra::String expected{"ONE STRING, TWO STRING"};
  EXPECT_EQ(expected, str_1 + str_2);
}

TEST(StkStringComparison, equalIdentifiers)
{
  sierra::Identifier str_1{"abc defg"};
  sierra::Identifier str_2{"Abc_dEfg"};
  EXPECT_TRUE(str_1 == str_2);
  EXPECT_TRUE(str_1 <= str_2);
  EXPECT_TRUE(str_1 >= str_2);
  EXPECT_FALSE(str_1 < str_2);
  EXPECT_FALSE(str_1 > str_2);
  EXPECT_FALSE(str_1 != str_2);
}

TEST(StkStringComparison, compareWithNullPtr)
{
  sierra::Identifier str_1{"abc defg"};
  const char * str_2 = nullptr;
  EXPECT_EQ(-1, str_1.compare(str_2));
}

TEST(StkStringComparison, compareEmptyWithNull)
{
  sierra::Identifier str_1{};
  const char * str_2 = nullptr;
  EXPECT_EQ(0, str_1.compare(str_2));
}

namespace {
  sierra::Identifier func_taking_sierra_identifier(const sierra::Identifier & identifier)
  {
    return identifier;
  }
}

TEST(StkStringIdentfierConversions, equalIdentifiers)
{
  sierra::String identifier_to_string_construction(sierra::Identifier("test"));
  sierra::String string_to_identifier_construction(sierra::String("test"));

  //test implicit move conversions -- future work to deprecate these
  func_taking_sierra_string(sierra::Identifier("test"));
  func_taking_sierra_identifier(sierra::String("test"));

  //test implicit copy conversions -- future work to deprecate these
  const sierra::Identifier sierra_identifier("test");
  const sierra::String sierra_string("test");
  func_taking_sierra_string(sierra_identifier);
  func_taking_sierra_identifier(sierra_string);

  //shouldbe nondeprecated versions
  func_taking_sierra_string(sierra::String(sierra_identifier));
  func_taking_sierra_identifier(sierra::Identifier(sierra_string));

  std::string_view test = std::string_view(sierra_string);
  test = "";
}
