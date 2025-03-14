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

#include "gtest/gtest.h"              // for Test, EXPECT_EQ, Message, TestPartResult, SuiteApiR...
#include "stk_topology/topology.hpp"  // for topology, topology::HEX_8

#include <climits>

#ifndef __IBMCPP__

namespace {

constexpr unsigned INVALID_ORDINAL = UINT_MAX;

void verify_side_ordinal_mapping(stk::topology topo,
                                 unsigned side_ordinal,
                                 unsigned expected_ranked_side_ordinal,
                                 stk::topology::rank_t expected_rank)
{
  unsigned ranked_side_ordinal;
  stk::topology::rank_t rank;

  // Convert side ordinal to a ranked side ordinal and associated rank
  // side_ordinal -> (expected_ranked_side_ordinal, expected_rank)
  topo.ranked_side_ordinal(side_ordinal, ranked_side_ordinal, rank);

  EXPECT_EQ(expected_ranked_side_ordinal, ranked_side_ordinal);
  EXPECT_EQ(expected_rank, rank);

  if(ranked_side_ordinal != INVALID_ORDINAL && expected_rank != stk::topology::INVALID_RANK) {
    // Reverse operation: (expected_ranked_side_ordinal, expected_rank) -> side_ordinal
    EXPECT_EQ(side_ordinal, topo.side_ordinal(ranked_side_ordinal, rank));
  }
}

TEST(side_ordinals, hex8)
{
  stk::topology hex8 = stk::topology::HEX_8;

  // Get all the possible ranks for sides. HEX_8 has only one .... FACE_RANK
  constexpr unsigned expected_num_side_ranks = 1;
  const std::vector<stk::topology::rank_t> expected_side_ranks{stk::topology::FACE_RANK};

  std::vector<stk::topology::rank_t> side_ranks(expected_num_side_ranks);
  hex8.side_ranks(side_ranks.data());

  EXPECT_EQ(expected_side_ranks, side_ranks);

  // Convert side ordinals to a ranked side ordinal and associated rank; and vice versa.
  // HEX_8 has only one side rank so the side ordinals are the same as the ranked side ordinals

  // HEX_8 has a total of 6 sides .. remember that ordinals are zero based
  unsigned side_ordinal;

  // Ordinal 0 -> (0, FACE_RANK)
  side_ordinal = 0;
  verify_side_ordinal_mapping(hex8, side_ordinal, 0, stk::topology::FACE_RANK);

  // Ordinal 1 -> (1, FACE_RANK)
  side_ordinal = 1;
  verify_side_ordinal_mapping(hex8, side_ordinal, 1, stk::topology::FACE_RANK);

  // Ordinal 2 -> (2, FACE_RANK)
  side_ordinal = 2;
  verify_side_ordinal_mapping(hex8, side_ordinal, 2, stk::topology::FACE_RANK);

  // Ordinal 3 -> (3, FACE_RANK)
  side_ordinal = 3;
  verify_side_ordinal_mapping(hex8, side_ordinal, 3, stk::topology::FACE_RANK);

  // Ordinal 4 -> (4, FACE_RANK)
  side_ordinal = 4;
  verify_side_ordinal_mapping(hex8, side_ordinal, 4, stk::topology::FACE_RANK);

  // Ordinal 5 -> (5, FACE_RANK)
  side_ordinal = 5;
  verify_side_ordinal_mapping(hex8, side_ordinal, 5, stk::topology::FACE_RANK);

  // Test an invalid ordinal
  // Ordinal 6 -> (INVALID, INVALID)
  side_ordinal = 6;
  verify_side_ordinal_mapping(hex8, side_ordinal, INVALID_ORDINAL, stk::topology::INVALID_RANK);
}

TEST(side_ordinals, shell_tri_3)
{
  stk::topology shell3 = stk::topology::SHELL_TRI_3;

  // Get all the possible ranks for sides. SHELL_TRI_3 has two ....FACE_RANK and EDGE_RANK
  constexpr unsigned expected_num_side_ranks = 2;
  const std::vector<stk::topology::rank_t> expected_side_ranks{stk::topology::FACE_RANK, stk::topology::EDGE_RANK};

  std::vector<stk::topology::rank_t> side_ranks(expected_num_side_ranks);
  shell3.side_ranks(side_ranks.data());

  EXPECT_EQ(expected_side_ranks, side_ranks);

  // Convert side ordinals to a ranked side ordinal and associated rank; and vice versa.
  // SHELL_TRI_3 has two side ranks so the side ordinals for FACE_RANK are the same as
  // the ranked side ordinals but the side ordinals for the EDGE_RANK will be greater
  // than the ranked side ordinals by the number of faces (2)

  // SHELL_TRI_3 has a total of 5 sides .. remember that ordinals are zero based
  unsigned side_ordinal;

  // Ordinal 0 -> (0, FACE_RANK)
  side_ordinal = 0;
  verify_side_ordinal_mapping(shell3, side_ordinal, 0, stk::topology::FACE_RANK);

  // Ordinal 1 -> (1, FACE_RANK)
  side_ordinal = 1;
  verify_side_ordinal_mapping(shell3, side_ordinal, 1, stk::topology::FACE_RANK);

  // Ordinal 2 -> (2, FACE_RANK)
  side_ordinal = 2;
  verify_side_ordinal_mapping(shell3, side_ordinal, 0, stk::topology::EDGE_RANK);

  // Ordinal 3 -> (3, FACE_RANK)
  side_ordinal = 3;
  verify_side_ordinal_mapping(shell3, side_ordinal, 1, stk::topology::EDGE_RANK);

  // Ordinal 4 -> (4, FACE_RANK)
  side_ordinal = 4;
  verify_side_ordinal_mapping(shell3, side_ordinal, 2, stk::topology::EDGE_RANK);

  // Test an invalid ordinal
  // Ordinal 5 -> (INVALID, INVALID)
  side_ordinal = 5;
  verify_side_ordinal_mapping(shell3, side_ordinal, INVALID_ORDINAL, stk::topology::INVALID_RANK);
}

TEST(side_ordinals, shell_tri_3_all_face_sides)
{
  stk::topology shell3 = stk::topology::SHELL_TRI_3_ALL_FACE_SIDES;

  // Get all the possible ranks for sides. SHELL_TRI_3_ALL_FACE_SIDES has only one .... FACE_RANK
  constexpr unsigned expected_num_side_ranks = 1;
  const std::vector<stk::topology::rank_t> expected_side_ranks{stk::topology::FACE_RANK};

  std::vector<stk::topology::rank_t> side_ranks(expected_num_side_ranks);
  shell3.side_ranks(side_ranks.data());

  EXPECT_EQ(expected_side_ranks, side_ranks);

  // Convert side ordinals to a ranked side ordinal and associated rank; and vice versa.
  // SHELL_TRI_3_ALL_FACE_SIDES has only one side rank so the side ordinals are the same
  // as the ranked side ordinals

  // SHELL_TRI_3_ALL_FACE_SIDES has a total of 5 sides .. remember that ordinals are zero based
  unsigned side_ordinal;

  // Ordinal 0 -> (0, FACE_RANK)
  side_ordinal = 0;
  verify_side_ordinal_mapping(shell3, side_ordinal, 0, stk::topology::FACE_RANK);

  // Ordinal 1 -> (1, FACE_RANK)
  side_ordinal = 1;
  verify_side_ordinal_mapping(shell3, side_ordinal, 1, stk::topology::FACE_RANK);

  // Ordinal 2 -> (2, FACE_RANK)
  side_ordinal = 2;
  verify_side_ordinal_mapping(shell3, side_ordinal, 2, stk::topology::FACE_RANK);

  // Ordinal 3 -> (3, FACE_RANK)
  side_ordinal = 3;
  verify_side_ordinal_mapping(shell3, side_ordinal, 3, stk::topology::FACE_RANK);

  // Ordinal 4 -> (4, FACE_RANK)
  side_ordinal = 4;
  verify_side_ordinal_mapping(shell3, side_ordinal, 4, stk::topology::FACE_RANK);

  // Test an invalid ordinal
  // Ordinal 5 -> (INVALID, INVALID)
  side_ordinal = 5;
  verify_side_ordinal_mapping(shell3, side_ordinal, INVALID_ORDINAL, stk::topology::INVALID_RANK);
}

}

#endif // __IBMCPP__
