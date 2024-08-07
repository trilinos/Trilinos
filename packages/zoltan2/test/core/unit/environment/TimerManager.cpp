// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//
// Testing the TimerManager class.
// TODO we only test that it doesn't crash.

#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_BasicIdentifierAdapter.hpp>
#include <Zoltan2_TimerManager.hpp>

#include <Teuchos_DefaultComm.hpp>

#ifdef _MSC_VER
#define NOMINMAX
#include <windows.h>
#else
#include <unistd.h>
#endif


using Zoltan2::TimerManager;
using Zoltan2::MACRO_TIMERS;
using Zoltan2::MICRO_TIMERS;
using Zoltan2::BOTH_TIMERS;


static void sleep_wrap(unsigned int seconds)
{
#ifdef _MSC_VER
  Sleep(1000*seconds);
#else
  sleep(seconds);
#endif
}

void goToSleep(const RCP<const Zoltan2::Environment> &env)
{
  env->timerStart(MICRO_TIMERS, string("sleep for 5 seconds"));
  sleep_wrap(5);
  env->timerStop(MICRO_TIMERS, string("sleep for 5 seconds"));

  env->timerStart(MICRO_TIMERS, string("sleep for 3 seconds (twice)"));
  sleep_wrap(3);
  env->timerStop(MICRO_TIMERS, string("sleep for 3 seconds (twice)"));

  env->timerStart(MICRO_TIMERS, string("sleep for 2 seconds"));
  sleep_wrap(2);
  env->timerStop(MICRO_TIMERS, string("sleep for 2 seconds"));

  env->timerStart(MICRO_TIMERS, string("sleep for 3 seconds (twice)"));
  sleep_wrap(3);
  env->timerStop(MICRO_TIMERS, string("sleep for 3 seconds (twice)"));
}


int main(int narg, char *arg[])
{
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  // Create a problem, requesting that Timing be enabled.

  Teuchos::ParameterList pl("test list");
  pl.set("timer_output_stream" , "std::cout");
  pl.set("timer_type" , "both_timers");
  std::vector<const zscalar_t * >weights;
  std::vector<int> strides;
  Array<zgno_t> someIds(10,1);
  typedef Zoltan2::BasicUserTypes<zscalar_t, zlno_t, zgno_t> myTypes_t;
  typedef Zoltan2::BasicIdentifierAdapter<myTypes_t> inputAdapter_t;
  inputAdapter_t ia(10, someIds.getRawPtr(), weights, strides);

  Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &pl, comm);

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
