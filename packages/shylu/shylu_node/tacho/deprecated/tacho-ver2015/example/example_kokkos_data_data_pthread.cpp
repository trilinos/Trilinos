#include <Kokkos_Core.hpp>

#include <Kokkos_Threads.hpp>

#include "Teuchos_CommandLineProcessor.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

typedef Kokkos::Threads exec_space;
typedef double value_type;

#include "example_kokkos_data_data.hpp"

using namespace Tacho;

int main (int argc, char *argv[]) {

  Teuchos::CommandLineProcessor clp;
  clp.setDocString("This example program measure the performance of task data parallelism (barrier) on Kokkos::Threads execution space.\n");

  int nthreads = 0;
  clp.setOption("nthreads", &nthreads, "Number of threads");

  int numa = 0;
  clp.setOption("numa", &numa, "Number of numa node");

  int core_per_numa = 0;
  clp.setOption("core-per-numa", &core_per_numa, "Number of cores per numa node");

  int league_size = 1;
  clp.setOption("league-size", &league_size, "League size");

  int team_size = 1;
  clp.setOption("team-size", &team_size, "Team size");

  int ntasks = 100;
  clp.setOption("ntasks", &ntasks, "Number of tasks to be spawned");

  bool verbose = false;
  clp.setOption("enable-verbose", "disable-verbose", &verbose, "Flag for verbose printing");

  clp.recogniseAllOptions(true);
  clp.throwExceptions(false);

  Teuchos::CommandLineProcessor::EParseCommandLineReturn r_parse= clp.parse( argc, argv );

  if (r_parse == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return 0;
  if (r_parse != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL  ) return -1;

  int r_val = 0;
  {
    exec_space::initialize(nthreads, numa, core_per_numa);
    exec_space::print_configuration(cout, true);

    r_val = exampleKokkosDataData<exec_space,value_type>((ntasks > MAXTASKS ? MAXTASKS : ntasks), league_size, team_size, verbose);

    exec_space::finalize();
  }

  return r_val;
}
