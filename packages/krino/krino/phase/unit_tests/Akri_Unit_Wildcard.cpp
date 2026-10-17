// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_Utility.hpp>

namespace krino {

TEST(Wildcard, HasWildcard)
{
  EXPECT_TRUE(utility::has_wildcard("*"));
  EXPECT_TRUE(utility::has_wildcard("block_*"));
  EXPECT_TRUE(utility::has_wildcard("*_dragon"));
  EXPECT_FALSE(utility::has_wildcard("block_1"));
  EXPECT_FALSE(utility::has_wildcard(""));
}

TEST(Wildcard, IsSupportedWildcardPattern)
{
  EXPECT_TRUE(utility::is_supported_wildcard_pattern("block_1"));
  EXPECT_TRUE(utility::is_supported_wildcard_pattern("*"));
  EXPECT_TRUE(utility::is_supported_wildcard_pattern("block_*_a"));
  EXPECT_FALSE(utility::is_supported_wildcard_pattern("**"));
  EXPECT_FALSE(utility::is_supported_wildcard_pattern("a*b*c"));
  EXPECT_FALSE(utility::is_supported_wildcard_pattern("block_?"));
  EXPECT_FALSE(utility::is_supported_wildcard_pattern("block_?*"));
}

TEST(Wildcard, PatternWithoutWildcard_MatchesOnlyExactName)
{
  EXPECT_TRUE(utility::is_wildcard_match("block_1", "block_1"));
  EXPECT_FALSE(utility::is_wildcard_match("block_1", "block_10"));
  EXPECT_FALSE(utility::is_wildcard_match("block_1", "block"));
}

TEST(Wildcard, StarAlone_MatchesAnyName)
{
  EXPECT_TRUE(utility::is_wildcard_match("*", "block_1"));
  EXPECT_TRUE(utility::is_wildcard_match("*", ""));
}

TEST(Wildcard, TrailingStar_MatchesPrefix)
{
  EXPECT_TRUE(utility::is_wildcard_match("layer*", "layer"));
  EXPECT_TRUE(utility::is_wildcard_match("layer*", "layer17"));
  EXPECT_FALSE(utility::is_wildcard_match("layer*", "lay"));
  EXPECT_FALSE(utility::is_wildcard_match("layer*", "block_layer1"));
}

TEST(Wildcard, LeadingStar_MatchesSuffix)
{
  EXPECT_TRUE(utility::is_wildcard_match("*_dragon", "block_1_dragon"));
  EXPECT_TRUE(utility::is_wildcard_match("*_dragon", "_dragon"));
  EXPECT_FALSE(utility::is_wildcard_match("*_dragon", "block_1_dragon_2"));
}

TEST(Wildcard, InteriorStar_MatchesPrefixAndSuffix)
{
  EXPECT_TRUE(utility::is_wildcard_match("block_*_a", "block_1_a"));
  EXPECT_TRUE(utility::is_wildcard_match("block_*_a", "block__a"));
  EXPECT_FALSE(utility::is_wildcard_match("block_*_a", "block_1_b"));
  EXPECT_FALSE(utility::is_wildcard_match("block_*_a", "xblock_1_a"));
}

TEST(Wildcard, NameShorterThanPrefixAndSuffix_DoesNotMatch)
{
  // "ab" starts with "ab" and ends with "b", but too short for both without overlapping.
  EXPECT_FALSE(utility::is_wildcard_match("ab*b", "ab"));
  EXPECT_TRUE(utility::is_wildcard_match("ab*b", "abb"));
}

} // namespace krino
