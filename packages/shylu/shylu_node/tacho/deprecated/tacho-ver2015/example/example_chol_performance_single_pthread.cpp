#include <Kokkos_Core.hpp>

#include <Kokkos_Threads.hpp>
#include <Threads/Kokkos_Threads_TaskPolicy.hpp>

#include "Teuchos_CommandLineProcessor.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

typedef Kokkos::Threads exec_space;

#include "example_chol_performance_single.hpp"

using namespace Tacho;

int main (int argc, char *argv[]) {

  Teuchos::CommandLineProcessor clp;
  clp.setDocString("This example program measure the performance of Chol algorithms on Kokkos::Threads execution space.\n");

  int nthreads = 0;
  clp.setOption("nthreads", &nthreads, "Number of threads");

  int numa = 0;
  clp.setOption("numa", &numa, "Number of numa node");

  int core_per_numa = 0;
  clp.setOption("core-per-numa", &core_per_numa, "Number of cores per numa node");

  int max_task_dependence = 10;
  clp.setOption("max-task-dependence", &max_task_dependence, "Max number of task dependence");

  int team_size = 1;
  clp.setOption("team-size", &team_size, "Team size");

  int fill_level = 0;
  clp.setOption("fill-level", &fill_level, "Fill level");

  int league_size = 1;
  clp.setOption("league-size", &league_size, "League size");

  bool team_interface = true;
  clp.setOption("enable-team-interface", "disable-team-interface",
                &team_interface, "Flag for team interface");

  bool vtune_symbolic = false;
  clp.setOption("enable-vtune-symbolic", "disable-vtune-symbolic", &vtune_symbolic, "Flag for vtune symbolic");

  bool vtune_serial = false;
  clp.setOption("enable-vtune-serial", "disable-vtune-serial", &vtune_serial, "Flag for vtune serial");

  bool vtune_task = false;
  clp.setOption("enable-vtune-task", "disable-vtune-task", &vtune_task, "Flag for vtune task");

  bool verbose = false;
  clp.setOption("enable-verbose", "disable-verbose", &verbose, "Flag for verbose printing");

  string file_input = "test.mtx";
  clp.setOption("file-input", &file_input, "Input file (MatrixMarket SPD matrix)");

  int treecut = 15;
  clp.setOption("treecut", &treecut, "Level to cut tree from bottom");

  int minblksize = 0;
  clp.setOption("minblksize", &minblksize, "Minimum block size for internal reordering");

  int prunecut = 0;
  clp.setOption("prunecut", &minblksize, "Level to prune tree from bottom");

  int seed = 0;
  clp.setOption("seed", &seed, "Seed for random number generator in graph partition");

  clp.recogniseAllOptions(true);
  clp.throwExceptions(false);

  Teuchos::CommandLineProcessor::EParseCommandLineReturn r_parse= clp.parse( argc, argv );

  if (r_parse == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return 0;
  if (r_parse != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL  ) return -1;

  int r_val = 0;
  {
    exec_space::initialize(nthreads, numa, core_per_numa);
    exec_space::print_configuration(cout, true);

    r_val = exampleCholPerformanceSingle
      <value_type,ordinal_type,size_type,exec_space,void>
      (file_input,
       treecut,
       minblksize,
       prunecut,
       seed,
       nthreads,
       max_task_dependence, team_size,
       fill_level, league_size,
       team_interface,
       (nthreads != 1),
       vtune_symbolic,
       vtune_serial,
       vtune_task,
       verbose);

    exec_space::finalize();
  }

  return r_val;
}
