// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_ParallelErrorMessage.hpp>
#include <sstream>
#include <string>
#include <gtest/gtest.h>

namespace krino
{

template <typename T> void write_rank_message(T & os, const int rank)
{
  os << "Rank " << rank << " message.\n";
}

TEST(ParallelErrorMessage, ConcatenateErrors)
{
  stk::Parallel world_comm(MPI_COMM_WORLD);
  ParallelErrorMessage err(world_comm);

  write_rank_message(err, world_comm.parallel_rank());
  auto global_message = err.gather_message();
  EXPECT_TRUE(global_message.first);

  if (world_comm.parallel_rank() == 0)
  {
    std::ostringstream expected;
    for (int i = 0; i < world_comm.parallel_size(); ++i)
    {
      write_rank_message(expected, i);
      expected << std::endl;
    }
    EXPECT_EQ(expected.str(), global_message.second);
  }
}
}
