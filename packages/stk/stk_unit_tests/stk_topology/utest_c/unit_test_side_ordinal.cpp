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

#include "gtest/gtest.h"              // for Message, TestPartResult, Test, EXPECT_EQ, Assertion...
#include "stk_topology/topology.hpp"  // for topology, operator<<, create_superedge_topology
#include <cstddef>                    // for size_t
#include <array>                      // for array
#include <iostream>                   // for operator<<, ostringstream, ostream, basic_ostream
#include <string>                     // for string, char_traits
#include <vector>
#include <algorithm>

struct GoldTranslatorOrdinal {
  unsigned ordinal;
  unsigned localRankedOrdinal;
  stk::topology::rank_t rank;
};

struct GoldTranslatorOrdinals {
  stk::topology topology;
  std::vector<GoldTranslatorOrdinal> goldOrdinals;
  std::vector<stk::topology::rank_t> goldSideRanks;

  GoldTranslatorOrdinals(stk::topology topology_,
                         const std::vector<GoldTranslatorOrdinal>& goldOrdinals_,
                         const std::vector<stk::topology::rank_t>& goldSideRanks_)
  : topology(topology_), goldOrdinals(goldOrdinals_), goldSideRanks(goldSideRanks_)
  { }
};

namespace {

std::vector<GoldTranslatorOrdinals> get_gold_translator_values()
{
  std::vector<GoldTranslatorOrdinals>
  goldData{ {stk::topology::QUAD_4_2D, {{0, 0, stk::topology::EDGE_RANK},
                                        {1, 1, stk::topology::EDGE_RANK},
                                        {2, 2, stk::topology::EDGE_RANK},
                                        {3, 3, stk::topology::EDGE_RANK}}, {stk::topology::EDGE_RANK}},

            {stk::topology::SHELL_QUAD_4, {{0, 0, stk::topology::FACE_RANK},
                                           {1, 1, stk::topology::FACE_RANK},
                                           {2, 0, stk::topology::EDGE_RANK},
                                           {3, 1, stk::topology::EDGE_RANK},
                                           {4, 2, stk::topology::EDGE_RANK},
                                           {5, 3, stk::topology::EDGE_RANK}}, {stk::topology::FACE_RANK, stk::topology::EDGE_RANK}},

            {stk::topology::SHELL_QUAD_4_ALL_FACE_SIDES, {{0, 0, stk::topology::FACE_RANK},
                                                          {1, 1, stk::topology::FACE_RANK},
                                                          {2, 2, stk::topology::FACE_RANK},
                                                          {3, 3, stk::topology::FACE_RANK},
                                                          {4, 4, stk::topology::FACE_RANK},
                                                          {5, 5, stk::topology::FACE_RANK}}, {stk::topology::FACE_RANK}},

            {stk::topology::SHELL_TRI_3, {{0, 0, stk::topology::FACE_RANK},
                                          {1, 1, stk::topology::FACE_RANK},
                                          {2, 0, stk::topology::EDGE_RANK},
                                          {3, 1, stk::topology::EDGE_RANK},
                                          {4, 2, stk::topology::EDGE_RANK}}, {stk::topology::FACE_RANK, stk::topology::EDGE_RANK}},

            {stk::topology::SHELL_TRI_3_ALL_FACE_SIDES, {{0, 0, stk::topology::FACE_RANK},
                                                         {1, 1, stk::topology::FACE_RANK},
                                                         {2, 2, stk::topology::FACE_RANK},
                                                         {3, 3, stk::topology::FACE_RANK},
                                                         {4, 4, stk::topology::FACE_RANK}}, {stk::topology::FACE_RANK}},

            {stk::topology::HEX_8, {{0, 0, stk::topology::FACE_RANK},
                                    {1, 1, stk::topology::FACE_RANK},
                                    {2, 2, stk::topology::FACE_RANK},
                                    {3, 3, stk::topology::FACE_RANK},
                                    {4, 4, stk::topology::FACE_RANK},
                                    {5, 5, stk::topology::FACE_RANK}}, {stk::topology::FACE_RANK}},

            {stk::topology::BEAM_2, {{0, 0, stk::topology::EDGE_RANK}}, {stk::topology::EDGE_RANK}},

            {stk::topology::SHELL_LINE_2, {{0, 0, stk::topology::EDGE_RANK},
                                           {1, 1, stk::topology::EDGE_RANK}}, {stk::topology::EDGE_RANK}}
  };

  return goldData;
}

TEST(stk_topology, test_topology_side_ranks)
{
  auto goldData = get_gold_translator_values();

  for(const auto& entry : goldData) {
    auto topo = entry.topology;

    std::vector<stk::topology::rank_t> sideRanks(topo.num_side_ranks());
    topo.side_ranks(sideRanks.data());

    EXPECT_EQ(entry.goldSideRanks, sideRanks);
  }
}

TEST(stk_topology, test_ordinal_translator)
{
  auto goldData = get_gold_translator_values();

  for(const auto& entry : goldData) {
    auto topo = entry.topology;

    for(const auto& data : entry.goldOrdinals) {
      unsigned ordinal = data.ordinal;

      unsigned goldLocalRankedOrdinal = data.localRankedOrdinal;
      stk::topology::rank_t goldRank = data.rank;

      unsigned localRankedOrdinal;
      stk::topology::rank_t rank;

      topo.ranked_side_ordinal(ordinal, localRankedOrdinal, rank);

      EXPECT_EQ(goldLocalRankedOrdinal, localRankedOrdinal) << ": failure for topology: " << topo << " and ordinal: " << ordinal;
      EXPECT_EQ(goldRank, rank) << ": failure for topology: " << topo << " and ordinal: " << ordinal;

      EXPECT_EQ(ordinal, topo.side_ordinal(localRankedOrdinal, rank)) << ": failure for topology: " << topo;
    }
  }
}

TEST(stk_topology, test_ordinal_translator_for_mock_bulk_data_ordinals)
{
  auto goldData = get_gold_translator_values();

  for(const auto& entry : goldData) {
    auto topo = entry.topology;

    std::vector<unsigned> goldOrdinalVec;
    std::vector<unsigned> ordinalVec(entry.goldOrdinals.size());

    for(const auto& data : entry.goldOrdinals) {
      goldOrdinalVec.push_back(data.ordinal);
    }

    std::transform(entry.goldOrdinals.begin(), entry.goldOrdinals.end(), ordinalVec.begin(),
                   [&topo](const GoldTranslatorOrdinal& data) {
                     return topo.side_ordinal(data.localRankedOrdinal, data.rank);
                   });

    EXPECT_EQ(goldOrdinalVec, ordinalVec);
  }
}

TEST(stk_topology, test_ordinal_translator_for_non_elements)
{
  std::vector<stk::topology> topologies = {stk::topology::LINE_2,
                                           stk::topology::SHELL_SIDE_BEAM_2};

  for(auto topo : topologies) {
    bool hasSides = topo.num_sides() > 0;
    bool hasSideRanks = topo.num_side_ranks() > 0;

    EXPECT_EQ(hasSides, hasSideRanks) << topo << " has inconsistent side definitions";
  }
}

} // namespace
