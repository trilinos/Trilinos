#include <Kokkos_Core.hpp>

#include <Kokkos_Threads.hpp>
#include <Threads/Kokkos_Threads_TaskPolicy.hpp>

#include "Teuchos_CommandLineProcessor.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

typedef Kokkos::Threads exec_space;

#include "example_chol_direct_plain.hpp"

using namespace Tacho;

int main (int argc, char *argv[]) {

  Teuchos::CommandLineProcessor clp;
  clp.setDocString("This example program measure the performance of Chol direct factorization (plain) on Kokkos::Threads execution space.\n");

  int nthreads = 0;
  clp.setOption("nthreads", &nthreads, "Number of threads");

  int numa = 0;
  clp.setOption("numa", &numa, "Number of numa node");

  int core_per_numa = 0;
  clp.setOption("core-per-numa", &core_per_numa, "Number of cores per numa node");

  int max_concurrency = 1024;
  clp.setOption("max-concurrency", &max_concurrency, "Max number of concurrent tasks");

  int max_task_dependence = 3;
  clp.setOption("max-task-dependence", &max_task_dependence, "Max number of task dependence");

  int team_size = 1;
  clp.setOption("team-size", &team_size, "Team size");

  int league_size = 1;
  clp.setOption("league-size", &league_size, "League size");

  bool team_interface = true;
  clp.setOption("enable-team-interface", "disable-team-interface",
                &team_interface, "Flag for team interface");

  bool verbose = false;
  clp.setOption("enable-verbose", "disable-verbose", &verbose, "Flag for verbose printing");

  string file_input = "test.mtx";
  clp.setOption("file-input", &file_input, "Input file (MatrixMarket SPD matrix)");

  int prunecut = 0;
  clp.setOption("prunecut", &prunecut, "Leve to prune tree from bottom");

  int seed = 0;
  clp.setOption("seed", &seed, "Seed for random number generator in graph partition");

  int nrhs = 1;
  clp.setOption("nrhs", &nrhs, "Numer of right hand side");

  int nb = 1;
  clp.setOption("nb", &nrhs, "Column block size for multiple right hand side");

  bool serial = false;
  clp.setOption("enable-serial", "disable-serial", &serial, "Flag for serial solve");

  bool solve = false;
  clp.setOption("enable-solve", "disable-solve", &solve, "Flag for trisolve");

  bool check = false;
  clp.setOption("enable-check", "disable-check", &check, "Flag for check solution");

  clp.recogniseAllOptions(true);
  clp.throwExceptions(false);

  Teuchos::CommandLineProcessor::EParseCommandLineReturn r_parse= clp.parse( argc, argv );

  if (r_parse == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return 0;
  if (r_parse != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL  ) return -1;

  int r_val = 0;
  {
    exec_space::initialize(nthreads, numa, core_per_numa);
    exec_space::print_configuration(cout, true);

    r_val = exampleCholDirectPlain
      <value_type,ordinal_type,size_type,exec_space,void>
      (file_input,
       prunecut,
       seed,
       nrhs,
       nb,
       nthreads, max_concurrency, max_task_dependence, team_size,
       league_size,
       team_interface,
       serial,
       solve,
       check,
       verbose);

    exec_space::finalize();
  }

  return r_val;
}
