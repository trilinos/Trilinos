#include <Kokkos_Core.hpp>

#include <Kokkos_Threads.hpp>
#include <Threads/Kokkos_Threads_TaskPolicy.hpp>

#include "Teuchos_CommandLineProcessor.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

typedef Kokkos::Threads exec_space;

#include "example_tri_solve_performance.hpp"

using namespace Example;

int main (int argc, char *argv[]) {

  Teuchos::CommandLineProcessor clp;
  clp.setDocString("This example program measure the performance of IChol algorithms on Kokkos::Threads execution space.\n");

  int nthreads = 1;
  clp.setOption("nthreads", &nthreads, "Number of threads");

  int max_task_dependences = 10;
  clp.setOption("max-task-depedences", &max_task_dependences, "Max number of task dependences");

  int team_size = 1;
  clp.setOption("team-size", &team_size, "Team size");

  bool verbose = false;
  clp.setOption("verbose", "non-verbose", &verbose, "Flag for verbose printing");

  string file_input = "test.mtx";
  clp.setOption("file-input", &file_input, "Input file (MatrixMarket SPD matrix)");

  int nrhs = 1;
  clp.setOption("nrhs", &nrhs, "Number of right hand side");

  int nb = nrhs;
  clp.setOption("nb", &nb, "Blocksize of right hand side");

  clp.recogniseAllOptions(true);
  clp.throwExceptions(false);

  Teuchos::CommandLineProcessor::EParseCommandLineReturn r_parse= clp.parse( argc, argv );

  if (r_parse == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return 0;
  if (r_parse != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL  ) return -1;
  
  int r_val = 0;
  {
    exec_space::initialize(nthreads);
    exec_space::print_configuration(cout, true);
    
    r_val = exampleTriSolvePerformance
      <value_type,ordinal_type,size_type,exec_space,void>
      (file_input, nrhs, nb, nthreads, max_task_dependences, team_size, verbose);

    exec_space::finalize();
  }

  return r_val;
}
