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
#include <stk_util/diag/Resource2.h>
#include <stk_unit_test_utils/timer.hpp>

namespace {

TEST(Resource2, MatchTiming)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }

  sierra::String parentName("parentResource");
  sierra::Rsrc2::ResourceRoot resourceRoot(parentName);
  sierra::Rsrc2::ResourceList resourceList;

  stk::unit_test_util::BatchTimer batchTimer(MPI_COMM_WORLD);
  batchTimer.initialize_batch_timer();

  constexpr int N = 25000;

  for(int i=0; i<N; ++i) {
    sierra::String name("ABC_"+std::to_string(i));
    resourceList.insert(resourceRoot.create_value<int>(name, i));
  }

  EXPECT_EQ(static_cast<size_t>(N), resourceList.size());

  constexpr int NumRuns = 25000;
  constexpr int NumOuterRuns = 5;

  for(int i=0; i<NumOuterRuns; ++i) {
    batchTimer.start_batch_timer();
    for(int r=0; r<NumRuns; ++r) {
      sierra::String name("ABC_"+std::to_string(r));
      sierra::Rsrc2::ResourceList matchList;
      resourceList.match(name, matchList);
    }
    batchTimer.stop_batch_timer();
  }

  batchTimer.print_batch_timing(NumRuns);
}

}

