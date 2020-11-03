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
#include <gtest/gtest.h>
#include <stk_tools/block_extractor/ParseCsv.hpp>

namespace
{

template<typename T>
void compare_result(const std::vector<T> & result, const std::vector<T> & expected)
{
  ASSERT_EQ(result.size(), expected.size());
  for (size_t i = 0; i < result.size(); ++i) {
    EXPECT_EQ(result[i], expected[i]);
  }
}

TEST(SplitString, emptyString)
{
  std::string input = "";
  std::vector<std::string> separated = stk::tools::split_string(input, ',');
  compare_result(separated, {});
}

TEST(SplitString, oneInteger)
{
  std::string input = "1";
  std::vector<std::string> separated = stk::tools::split_string(input, ',');
  compare_result(separated, {"1"});
}

TEST(SplitString, twoIntegers)
{
  std::string input = "13,2";
  std::vector<std::string> separated = stk::tools::split_string(input, ',');
  compare_result(separated, {"13", "2"});
}

TEST(SplitString, twoIntegersDifferentSeparator)
{
  std::string input = "1:4:9";
  std::vector<std::string> separated = stk::tools::split_string(input, ':');
  compare_result(separated, {"1", "4", "9"});
}

TEST(SplitString, tenIntegers)
{
  std::string input = "1,2,3,4,5,6,7,8,9,10";
  std::vector<std::string> separated = stk::tools::split_string(input, ',');
  compare_result(separated, {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"});
}

TEST(SplitString, spacesInString)
{
  std::string input = " 5 , 777 ";
  std::vector<std::string> separated = stk::tools::split_string(input, ',');
  compare_result(separated, {"5", "777"});
}

TEST(SplitString, integersAndRanges)
{
  std::string input = "1,10:20:2,25";
  std::vector<std::string> separated = stk::tools::split_string(input, ',');
  compare_result(separated, {"1", "10:20:2", "25"});
}

TEST(GetIdsFromStrings, empty)
{
  std::vector<std::string> input = {};
  std::vector<int> result = stk::tools::get_ids_from_strings(input);
  compare_result(result, {});
}

TEST(GetIdsFromStrings, singleInteger)
{
  std::vector<std::string> input = {"5"};
  std::vector<int> result = stk::tools::get_ids_from_strings(input);
  compare_result(result, {5});
}

TEST(GetIdsFromStrings, multipleIntegers)
{
  std::vector<std::string> input = {"1", "2", "3", "5", "8"};
  std::vector<int> result = stk::tools::get_ids_from_strings(input);
  compare_result(result, {1, 2, 3, 5, 8});
}

TEST(GetIdsFromStrings, nonInteger)
{
  std::vector<std::string> input = {"1", "block_2"};
  EXPECT_THROW(stk::tools::get_ids_from_strings(input), std::runtime_error);
}

TEST(GetIdsFromStrings, singleRangeNonInteger)
{
  std::vector<std::string> input = {"1:block_5"};
  EXPECT_THROW(stk::tools::get_ids_from_strings(input), std::runtime_error);
}

TEST(GetIdsFromStrings, singleRangeTooManyFields)
{
  std::vector<std::string> input = {"1:5:2:2"};
  EXPECT_THROW(stk::tools::get_ids_from_strings(input), std::logic_error);
}

TEST(GetIdsFromStrings, singleRangeOutOfOrder)
{
  std::vector<std::string> input = {"5:1"};
  EXPECT_THROW(stk::tools::get_ids_from_strings(input), std::logic_error);
}

TEST(GetIdsFromStrings, singleRangeZeroLimit)
{
  std::vector<std::string> input = {"0:5"};
  EXPECT_THROW(stk::tools::get_ids_from_strings(input), std::logic_error);
}

TEST(GetIdsFromStrings, singleRangeNegativeLimit)
{
  std::vector<std::string> input = {"1:-5"};
  EXPECT_THROW(stk::tools::get_ids_from_strings(input), std::logic_error);
}

TEST(GetIdsFromStrings, singleRangeNegativeStride)
{
  std::vector<std::string> input = {"1:5:-1"};
  EXPECT_THROW(stk::tools::get_ids_from_strings(input), std::logic_error);
}

TEST(GetIdsFromStrings, singleRangeNoStride)
{
  std::vector<std::string> input = {"1:5"};
  std::vector<int> result = stk::tools::get_ids_from_strings(input);
  compare_result(result, {1, 2, 3, 4, 5});
}

TEST(GetIdsFromStrings, singleRangeStride1)
{
  std::vector<std::string> input = {"1:5:1"};
  std::vector<int> result = stk::tools::get_ids_from_strings(input);
  compare_result(result, {1, 2, 3, 4, 5});
}

TEST(GetIdsFromStrings, singleRangeStride2)
{
  std::vector<std::string> input = {"1:5:2"};
  std::vector<int> result = stk::tools::get_ids_from_strings(input);
  compare_result(result, {1, 3, 5});
}

TEST(GetIdsFromStrings, singleRangeStride3)
{
  std::vector<std::string> input = {"1:5:3"};
  std::vector<int> result = stk::tools::get_ids_from_strings(input);
  compare_result(result, {1, 4});  // Not inclusive of range end
}

TEST(GetIdsFromStrings, singleRangeStrideTooBig)
{
  std::vector<std::string> input = {"1:5:5"};
  std::vector<int> result = stk::tools::get_ids_from_strings(input);
  compare_result(result, {1});  // Only first value
}

TEST(GetIdsFromStrings, mixedIntsAndRanges)
{
  std::vector<std::string> input = {"1", "10:14:2", "20"};
  std::vector<int> result = stk::tools::get_ids_from_strings(input);
  compare_result(result, {1, 10, 12, 14, 20});
}

TEST(GetIdsFromStrings, multipleRanges)
{
  std::vector<std::string> input = {"1:5", "10:14:2"};
  std::vector<int> result = stk::tools::get_ids_from_strings(input);
  compare_result(result, {1, 2, 3, 4, 5, 10, 12, 14});
}

}
