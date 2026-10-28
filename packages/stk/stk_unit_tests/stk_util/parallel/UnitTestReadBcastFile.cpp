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
#include "stk_util/parallel/Parallel.hpp"
#include "stk_util/parallel/mpi_filebuf.hpp"
#include <iostream>
#include <fstream>

#if defined ( STK_HAS_MPI )
TEST(Parallel, read_and_broadcast_file)
{
    stk::ParallelMachine comm = stk::parallel_machine_world();
    const int myProc = stk::parallel_machine_rank(comm);
    const int numProcs = stk::parallel_machine_size(comm);
    if (numProcs < 2 || numProcs > 4) { GTEST_SKIP(); }

    std::string expected = "line 1\nline 2\nline 3: more stuff\nline 4";
    std::string fileName("test.txt");

    if (myProc == 0) {
      std::ofstream out(fileName);
      out << expected;
    }

    stk::parallel_machine_barrier(comm);

    std::string fileContents = stk::read_and_broadcast_file(fileName, comm);

    EXPECT_EQ(expected, fileContents);
}

#endif

