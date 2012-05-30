// @HEADER
// ***********************************************************************
//                Copyright message goes here.
// ***********************************************************************
//
// Testing the TimerManager class.
// TODO we only test that it doesn't crash.

#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_BasicIdentifierInput.hpp>

#include <Teuchos_DefaultComm.hpp>
#include <unistd.h>

using Zoltan2::TimerManager;
using Zoltan2::MACRO_TIMERS;
using Zoltan2::MICRO_TIMERS;
using Zoltan2::BOTH_TIMERS;

typedef Zoltan2::BasicUserTypes<scalar_t, gno_t, lno_t, gno_t> myTypes_t;
typedef Zoltan2::BasicIdentifierInput<myTypes_t> inputAdapter_t;

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
