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

#include "gtest/gtest.h"              // for AssertionResult, Test, Message, TestPartResult, EXP...
#include "stk_topology/topology.hpp"  // for topology, EquivalentPermutation, topology::HEX_8

namespace {

//EquivalentElements
TEST(stk_topology_understanding, equivalent_elements)
{
  stk::EquivalentPermutation areElementsEquivalent;

  {
    if (stk::topology::topology_type<stk::topology::HEX_8>::num_permutations > 1) {
      unsigned hex1[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
      unsigned hex2[8] = { 4, 7, 6, 5, 0, 3, 2, 1 };
      unsigned hex3[8] = { 4, 5, 6, 7, 0, 1, 2, 3 };

      stk::topology hex8 = stk::topology::HEX_8;

      areElementsEquivalent = hex8.is_equivalent((unsigned*)hex1, (unsigned*)hex2);
      EXPECT_TRUE(areElementsEquivalent.is_equivalent);
      areElementsEquivalent = hex8.is_equivalent((unsigned*)hex1, (unsigned*)hex3);
      EXPECT_FALSE(areElementsEquivalent.is_equivalent);
    }
  }

  {
    unsigned triangle_1[3] = {0, 1, 2};
    unsigned triangle_2[3] = {0, 2, 1};

    stk::topology triangular_shell = stk::topology::SHELL_TRIANGLE_3;

    areElementsEquivalent = triangular_shell.is_equivalent((unsigned*)triangle_1, (unsigned*)triangle_2);

    EXPECT_TRUE(areElementsEquivalent.is_equivalent);

    unsigned permutation_index = areElementsEquivalent.permutation_number;
    unsigned goldValue = 3;
    EXPECT_EQ(goldValue, permutation_index); // From previous unit test, this is the 4th permutation
  }

}
//EndEquivalentElements
}

