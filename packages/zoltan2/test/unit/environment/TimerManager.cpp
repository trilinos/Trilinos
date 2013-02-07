// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
//
// Testing the TimerManager class.
// TODO we only test that it doesn't crash.

#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_BasicIdentifierInput.hpp>
#include <Zoltan2_TimerManager.hpp>

#include <Teuchos_DefaultComm.hpp>
#include <unistd.h>

using Zoltan2::TimerManager;
using Zoltan2::MACRO_TIMERS;
using Zoltan2::MICRO_TIMERS;
using Zoltan2::BOTH_TIMERS;

typedef Zoltan2::BasicUserTypes<scalar_t, gno_t, lno_t, gno_t> myTypes_t;
typedef Zoltan2::BasicIdentifierAdapter<myTypes_t> inputAdapter_t;

void goToSleep(const RCP<const Zoltan2::Environment> &env)
{
  env->timerStart(MICRO_TIMERS, string("sleep for 5 seconds"));
  sleep(5);
  env->timerStop(MICRO_TIMERS, string("sleep for 5 seconds"));

  env->timerStart(MICRO_TIMERS, string("sleep for 3 seconds (twice)"));
  sleep(3);
  env->timerStop(MICRO_TIMERS, string("sleep for 3 seconds (twice)"));

  env->timerStart(MICRO_TIMERS, string("sleep for 2 seconds"));
  sleep(2);
  env->timerStop(MICRO_TIMERS, string("sleep for 2 seconds"));

  env->timerStart(MICRO_TIMERS, string("sleep for 3 seconds (twice)"));
  sleep(3);
  env->timerStop(MICRO_TIMERS, string("sleep for 3 seconds (twice)"));
}


int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();

  // Create a problem, requesting that Timing be enabled.

  Teuchos::ParameterList pl("test list");
  pl.set("timer_output_stream" , "std::cout");
  pl.set("timer_type" , "both_timers");
  std::vector<const scalar_t * >weights;
  std::vector<int> strides;
  Array<gno_t> someIds(10,1);
  inputAdapter_t ia(10, someIds.getRawPtr(), weights, strides);

#ifdef HAVE_ZOLTAN2_MPI
  Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &pl,
    MPI_COMM_WORLD);
#else
  Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &pl);
#endif

  // Use the timers through the environment.

  const RCP<const Zoltan2::Environment> &env = problem.getEnvironment();

  if (comm->getRank() == 0)
    std::cout << "Sleeping..." << std::endl;

  env->timerStart(MACRO_TIMERS, string("Do the sleep test"));
  goToSleep(env);
  env->timerStop(MACRO_TIMERS, string("Do the sleep test"));

  comm->barrier();

  // Should show an error
  env->timerStop(MACRO_TIMERS, string("unstarted timer"));

  problem.printTimers();

  if (comm->getRank() == 0)
    std::cout << "PASS" << std::endl;
}
