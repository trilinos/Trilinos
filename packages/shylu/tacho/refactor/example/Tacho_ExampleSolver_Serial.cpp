#include <Kokkos_Core.hpp>

#include "Teuchos_CommandLineProcessor.hpp"

#include "ShyLUTacho_config.h"

typedef double value_type;
typedef int    ordinal_type;
typedef size_t size_type;

typedef Kokkos::Serial exec_space;

#if (defined(HAVE_SHYLUTACHO_SCOTCH) && (defined(HAVE_SHYLUTACHO_CHOLMOD) \
        || defined(HAVE_SHYLUTACHO_AMESOS)))
#include "Tacho_ExampleSolver.hpp"
using namespace Tacho;
#endif

int main (int argc, char *argv[]) {

#ifdef HAVE_SHYLUTACHO_VTUNE
  __itt_pause();
#endif

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

  int treecut = 4;
  clp.setOption("treecut", &treecut, "Level to cut tree from bottom");

  int prunecut = 0;
  clp.setOption("prunecut", &prunecut, "Level to prune tree from bottom");

  int max_concurrency = 250000;
  clp.setOption("max-concurrency", &max_concurrency, "Max number of concurrent tasks");

  int memory_pool_grain_size = 16;
  clp.setOption("memory-pool-grain-size", &memory_pool_grain_size, "Memorypool chunk size (12 - 16)");

  int mkl_nthreads = 1;
  clp.setOption("mkl-nthreads", &mkl_nthreads, "MKL threads for nested parallelism");

  int nrhs = 1;
  clp.setOption("nrhs", &nrhs, "# of right hand side");

  int mb = 256;
  clp.setOption("mb", &mb, "Dense nested blocks size");

  int nb = 4;
  clp.setOption("nb", &nb, "Column block size of right hand side");

  //int flat_maxsize = 32000;
  //clp.setOption("flat-maxsize", &flat_maxsize, "Cross-over size for flat algorithm");

  int hier_minsize = 4;
  clp.setOption("hier-minsize", &hier_minsize, "Cross-over size for hier algorithm");

  clp.recogniseAllOptions(true);
  clp.throwExceptions(false);

  Teuchos::CommandLineProcessor::EParseCommandLineReturn r_parse= clp.parse( argc, argv );
  
  if (r_parse == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return 0;
  if (r_parse != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL  ) return -1;

  int r_val = 0;
  {
    exec_space::initialize(); //nthreads, numa, core_per_numa);
    
#if (defined(HAVE_SHYLUTACHO_SCOTCH) && (defined(HAVE_SHYLUTACHO_CHOLMOD) || defined(HAVE_SHYLUTACHO_AMESOS)))
    r_val = exampleSolver<exec_space>
      (file_input, 
       treecut,
       prunecut, 
       max_concurrency, memory_pool_grain_size, mkl_nthreads,
       nrhs, mb, nb,
       hier_minsize,
       verbose);
#else
    r_val = -1;
    std::cout << "Scotch or Cholmod is NOT configured in Trilinos" << std::endl;
#endif

    exec_space::finalize();
  }
  
  return r_val;
}
