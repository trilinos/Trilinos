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

class CommaSeparatedValues : public ::testing::Test
{
protected:
};

TEST_F(CommaSeparatedValues, emptyString_emptyList)
{
    std::string input = "";
    std::vector<std::string> separated = stk::tools::get_csv(input);
    ASSERT_TRUE(separated.empty());
}

TEST_F(CommaSeparatedValues, oneItem_stringsAreEqual)
{
    std::string input = "1";
    std::vector<std::string> separated = stk::tools::get_csv(input);
    ASSERT_EQ(1u, separated.size());
    EXPECT_EQ(input, separated[0]);
}

TEST_F(CommaSeparatedValues, twoItems_twoSeparatedStrings)
{
    std::string input = "13,2";
    std::vector<std::string> separated = stk::tools::get_csv(input);
    ASSERT_EQ(2u, separated.size());
    EXPECT_EQ("13", separated[0]);
    EXPECT_EQ("2", separated[1]);
}

TEST_F(CommaSeparatedValues, tenItems_tenSeparatedStrings)
{
    std::string input = "1,2,3,4,5,6,7,8,9,10";
    std::vector<std::string> separated = stk::tools::get_csv(input);
    ASSERT_EQ(10u, separated.size());
    for(size_t i=0; i<10; i++)
        EXPECT_EQ(std::to_string(i+1), separated[i]);
}

TEST_F(CommaSeparatedValues, spaceInString_stripsSpaces)
{
    std::string input = "5, 777";
    std::vector<std::string> separated = stk::tools::get_csv(input);
    ASSERT_EQ(2u, separated.size());
    EXPECT_EQ("5", separated[0]);
    EXPECT_EQ("777", separated[1]);
}

TEST(VectorOfStrings, gettingBlockNames_addedBlockUnderscore)
{
    std::vector<std::string> input = {"5", "777"};
    std::vector<std::string> output = stk::tools::get_block_names_given_ids(input);
    ASSERT_EQ(2u, output.size());
    EXPECT_EQ("block_5", output[0]);
    EXPECT_EQ("block_777", output[1]);
}

}
