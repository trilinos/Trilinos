#include <Kokkos_Core.hpp>
#include <Kokkos_Threads.hpp>
#include "Teuchos_CommandLineProcessor.hpp"

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

typedef Kokkos::Threads exec_space;

#include "Tacho_ExampleGraphTools.hpp"

using namespace Tacho;

int main (int argc, char *argv[]) {

  Teuchos::CommandLineProcessor clp;
  clp.setDocString("Tacho::DenseMatrixBase examples on Pthreads execution space.\n");

  int nthreads = 0;
  clp.setOption("nthreads", &nthreads, "Number of threads");

  int numa = 0;
  clp.setOption("numa", &numa, "Number of numa node");

  int core_per_numa = 0;
  clp.setOption("core-per-numa", &core_per_numa, "Number of cores per numa node");

  bool verbose = false;
  clp.setOption("enable-verbose", "disable-verbose", &verbose, "Flag for verbose printing");

  std::string file_input = "test.mtx";
  clp.setOption("file-input", &file_input, "Input file (MatrixMarket SPD matrix)");

  int treecut = 0;
  clp.setOption("treecut", &treecut, "Level to cut tree from bottom");

  int prunecut = 0;
  clp.setOption("prunecut", &prunecut, "Level to prune tree from bottom");

  clp.recogniseAllOptions(true);
  clp.throwExceptions(false);

  Teuchos::CommandLineProcessor::EParseCommandLineReturn r_parse= clp.parse( argc, argv );

  if (r_parse == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return 0;
  if (r_parse != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL  ) return -1;

  int r_val = 0;
  {
    exec_space::initialize(nthreads, numa, core_per_numa);

    r_val = exampleGraphTools<exec_space>
      (file_input, treecut, prunecut, verbose);
    
    exec_space::finalize();
  }
  
  return r_val;
}
