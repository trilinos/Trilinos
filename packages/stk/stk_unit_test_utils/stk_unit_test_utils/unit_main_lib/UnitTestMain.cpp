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

#include <stk_util/stk_config.h>
#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>                // for InitGoogleTest, etc
#ifdef STK_HAVE_STKNGP_TEST
#include <stk_ngp_test/ngp_test.hpp>
#endif
#include <stk_unit_test_utils/ParallelGtestOutput.hpp>
#include <stk_util/parallel/Parallel.hpp>

int gl_argc = 0;
char** gl_argv = 0;

int main(int argc, char **argv)
{
  stk::parallel_machine_init(&argc, &argv);

  int returnVal = -1;
 
  {
#ifdef STK_HAVE_STKNGP_TEST
    ngp_testing::NgpTestEnvironment testEnv(&argc, argv);
#else
    testing::InitGoogleTest(&argc, argv);
#endif

    gl_argc = argc;
    gl_argv = argv;

#ifdef STK_HAS_MPI
    int procId = stk::parallel_machine_rank(MPI_COMM_WORLD);
    stk::unit_test_util::create_parallel_output(procId);
#endif

#ifdef STK_HAVE_STKNGP_TEST
    returnVal = testEnv.run_all_tests();
    testEnv.finalize();
#else
    returnVal = RUN_ALL_TESTS();
#endif
  }

  stk::parallel_machine_finalize();

  return returnVal;
}
