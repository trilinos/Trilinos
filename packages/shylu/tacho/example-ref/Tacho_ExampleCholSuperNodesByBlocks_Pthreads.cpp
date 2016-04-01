#include <Kokkos_Core.hpp>
#include <Kokkos_Threads.hpp>
#include <Threads/Kokkos_Threads_TaskPolicy.hpp>  

#include "Teuchos_CommandLineProcessor.hpp"

#include "ShyLUTacho_config.h"

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

typedef Kokkos::Threads exec_space;

#if (defined(HAVE_SHYLUTACHO_SCOTCH) && defined(HAVE_SHYLUTACHO_CHOLMOD))
#include "Tacho_ExampleCholSuperNodesByBlocks.hpp"
using namespace Tacho;
#endif

int main (int argc, char *argv[]) {

  Teuchos::CommandLineProcessor clp;
  clp.setDocString("Tacho::DenseMatrixBase examples on Pthreads execution space.\n");

  int nthreads = 0;
  clp.setOption("nthreads", &nthreads, "Number of threads");

  int numa = 0;
  clp.setOption("numa", &numa, "Number of numa node");

  int core_per_numa = 0;
  clp.setOption("core-per-numa", &core_per_numa, "Number of cores per numa node");

  bool check = true;
  clp.setOption("enable-check", "disable-check", &check, "Flag for check solution with serial execution");

  bool verbose = false;
  clp.setOption("enable-verbose", "disable-verbose", &verbose, "Flag for verbose printing");

  std::string file_input = "test.mtx";
  clp.setOption("file-input", &file_input, "Input file (MatrixMarket SPD matrix)");

  int treecut = 0;
  clp.setOption("treecut", &treecut, "Level to cut tree from bottom");

  int prunecut = 0;
  clp.setOption("prunecut", &prunecut, "Level to prune tree from bottom");

  int fill_level = -1;
  clp.setOption("fill-level", &fill_level, "Fill level");

  int rows_per_team = 4096;
  clp.setOption("rows-per-team", &rows_per_team, "Workset size");

  int max_concurrency = 4096;
  clp.setOption("max-concurrency", &max_concurrency, "Max number of concurrent tasks");

  int max_task_dependence = 3;
  clp.setOption("max-task-dependence", &max_task_dependence, "Max number of task dependence");

  int team_size = 1;
  clp.setOption("team-size", &team_size, "Team size");

  clp.recogniseAllOptions(true);
  clp.throwExceptions(false);

  Teuchos::CommandLineProcessor::EParseCommandLineReturn r_parse= clp.parse( argc, argv );

  if (r_parse == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return 0;
  if (r_parse != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL  ) return -1;

  int r_val = 0;
  {
    exec_space::initialize(nthreads, numa, core_per_numa);

#if (defined(HAVE_SHYLUTACHO_SCOTCH) && defined(HAVE_SHYLUTACHO_CHOLMOD))
    r_val = exampleCholSuperNodesByBlocks<exec_space>
      (file_input, 
       treecut, prunecut, fill_level, rows_per_team, 
       max_concurrency, max_task_dependence, team_size,
       check, verbose);
#else
    r_val = -1;
    std::cout << "Scotch or Cholmod is NOT configured in Trilinos" << std::endl;
#endif

    exec_space::finalize();
  }
  
  return r_val;
}
