#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include <Kokkos_Threads.hpp>
#include <Threads/Kokkos_Threads_TaskPolicy.hpp>

#include "Teuchos_CommandLineProcessor.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

typedef Kokkos::Threads exec_space;

#include "example_kokkos_team_barrier.hpp"
#define NUM_BARRIER_APPLIED 10000

class ReferenceBarrierFunctor {
public:
  typedef Kokkos::TeamPolicy<exec_space> policy_type;
  typedef exec_space execution_space;
  
  ReferenceBarrierFunctor() {}
  
  KOKKOS_INLINE_FUNCTION
  void operator() (typename policy_type::member_type member) const {
    for (int i=0;i<NUM_BARRIER_APPLIED;++i)
      member.team_barrier();
  }
};

class TestBarrierFunctor {
private:
  Kokkos::Experimental::SimpleCoreBarrierType *_corebarrier;

public:
  typedef Kokkos::TeamPolicy<exec_space> policy_type;
  typedef exec_space execution_space;

  TestBarrierFunctor(Kokkos::Experimental::SimpleCoreBarrierType *corebarrier)
    : _corebarrier(corebarrier) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (typename policy_type::member_type member) const {
    Kokkos::Experimental::SimpleTeamMemberType me(member.team_size(),
                                                  member.team_rank(),
                                                  member.team_size() - member.team_rank(),
                                                  &_corebarrier[member.league_rank()]);
    for (int i=0;i<NUM_BARRIER_APPLIED;++i)
      me.team_barrier();
  }
};

int main (int argc, char *argv[]) {

  Teuchos::CommandLineProcessor clp;
  clp.setDocString("This example program measure the performance of a team barrier on Kokkos::Threads execution space.\n");

  int nthreads = 0;
  clp.setOption("nthreads", &nthreads, "Number of threads");

  int numa = 0;
  clp.setOption("numa", &numa, "Number of numa node");

  int core_per_numa = 0;
  clp.setOption("core-per-numa", &core_per_numa, "Number of cores per numa node");

  int team_size = 1;
  clp.setOption("team-size", &team_size, "Team size");

  int league_size = 1;
  clp.setOption("league-size", &league_size, "League size");

  bool verbose = false;
  clp.setOption("enable-verbose", "disable-verbose", &verbose, "Flag for verbose printing");

  int niter = 100;
  clp.setOption("niter", &niter, "Number of iterations for testing");

  clp.recogniseAllOptions(true);
  clp.throwExceptions(false);

  Teuchos::CommandLineProcessor::EParseCommandLineReturn r_parse= clp.parse( argc, argv );

  if (r_parse == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return 0;
  if (r_parse != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL  ) return -1;

  int r_val = 0;
  {
    exec_space::initialize(nthreads, numa, core_per_numa);
    exec_space::print_configuration(cout, true);

    Kokkos::Impl::Timer timer;

    {
      timer.reset();
      for (int i=0;i<niter;++i) {
        Kokkos::parallel_for(Kokkos::TeamPolicy<exec_space>(league_size, team_size),
                             ReferenceBarrierFunctor());
      }
      exec_space::fence();
      double t = timer.seconds()/niter;
      cout << "ReferenceTeamBarrier:: time = " << t << endl;
    }
    {
      Kokkos::Experimental::SimpleCoreBarrier corebarrier[60];

      timer.reset();
      for (int i=0;i<niter;++i) {      
        Kokkos::parallel_for(Kokkos::TeamPolicy<exec_space>(league_size, team_size),
                             TestBarrierFunctor(&corebarrier[0]));
      }
      exec_space::fence();
      double t = timer.seconds()/niter;
      cout << "TestTeamBarrier:: time = " << t << endl;
    }

    exec_space::finalize();
  }

  return r_val;
  }
