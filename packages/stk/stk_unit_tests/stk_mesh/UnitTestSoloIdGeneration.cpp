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

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t
#include <vector>                       // for vector, operator==
#include "mpi.h"                        // for MPI_COMM_WORLD, etc
#include "stk_mesh/base/Types.hpp"      // for EntityId, PartVector, etc
#include "stk_mesh/baseImpl/SoloSideIdGenerator.hpp"

class SoloSideIdGeneratorTester : public stk::mesh::impl::SoloSideIdGenerator
{
public:
  SoloSideIdGeneratorTester(int numProcsArg, int proc, uint64_t maxSideId)
    : SoloSideIdGenerator(numProcsArg, proc, maxSideId)
  {}

  stk::mesh::EntityId test_get_solo_side_id_using_formula(unsigned elementId, unsigned sideOrdinal)
  {
    return get_solo_side_id_using_formula(elementId, sideOrdinal);
  }
};

TEST( SoloSideIdGenerator, get_solo_side_id_using_formula )
{
  const int numProcs = 1;
  const int myProcId = 0;
  const uint64_t maxSideId = 100;
  SoloSideIdGeneratorTester soloSideIdGenerator(numProcs, myProcId,  maxSideId);

  const unsigned maxPseudoElement = soloSideIdGenerator.max_pseudo_element();
  EXPECT_EQ(17u, soloSideIdGenerator.test_get_solo_side_id_using_formula(1, 6));
  EXPECT_EQ(18u, soloSideIdGenerator.test_get_solo_side_id_using_formula(1, 7));
  EXPECT_EQ(19u, soloSideIdGenerator.test_get_solo_side_id_using_formula(1, 8));
  EXPECT_EQ(20u, soloSideIdGenerator.test_get_solo_side_id_using_formula(1, 9));
  EXPECT_EQ(27u, soloSideIdGenerator.test_get_solo_side_id_using_formula(2, 6));
  EXPECT_NO_THROW(soloSideIdGenerator.test_get_solo_side_id_using_formula(maxPseudoElement,   9));
  EXPECT_THROW(   soloSideIdGenerator.test_get_solo_side_id_using_formula(maxPseudoElement+1, 6), std::exception);
}

TEST( SoloSideIdGenerator, get_solo_side_id )
{
  const int numProcs = 1;
  const int myProcId = 0;
  const uint64_t maxSideId = 100;
  SoloSideIdGeneratorTester soloSideIdGenerator(numProcs, myProcId,  maxSideId);

  EXPECT_EQ(17u, soloSideIdGenerator.get_solo_side_id());
  EXPECT_EQ(18u, soloSideIdGenerator.get_solo_side_id());
  EXPECT_EQ(19u, soloSideIdGenerator.get_solo_side_id());
  EXPECT_EQ(20u, soloSideIdGenerator.get_solo_side_id());
  EXPECT_EQ(27u, soloSideIdGenerator.get_solo_side_id());
}

TEST( SoloSideIdGenerator, test_id_space_per_proc )
{
  int myProcId = -1;
  int numProcs = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &myProcId);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  const uint64_t maxSideId = 1000000;
  SoloSideIdGeneratorTester soloSideIdGenerator(numProcs, myProcId, maxSideId);
  uint64_t maxPseudoElement = soloSideIdGenerator.max_pseudo_element();

  for(size_t elementIndex = myProcId+1; elementIndex <= maxPseudoElement; elementIndex+=numProcs)
  {
    for(size_t ordinal = 6; ordinal < 10; ordinal++)
    {
      EXPECT_EQ(soloSideIdGenerator.get_solo_side_id(), soloSideIdGenerator.test_get_solo_side_id_using_formula(elementIndex, ordinal));
    }
  }

  EXPECT_THROW(soloSideIdGenerator.get_solo_side_id(), std::exception);
}
