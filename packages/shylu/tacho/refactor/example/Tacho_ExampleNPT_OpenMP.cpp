#include <Kokkos_Core.hpp>

#include "Teuchos_CommandLineProcessor.hpp"

#include "ShyLUTacho_config.h"

typedef double ValueType;

#include "Tacho_ExampleNPT.hpp"

int main (int argc, char *argv[]) {
  
#ifdef HAVE_SHYLUTACHO_VTUNE
  __itt_pause();
#endif

  Teuchos::CommandLineProcessor clp;
  clp.setDocString("Tacho::DenseMatrixBase examples on Pthreads execution space.\n");

  bool verbose_blocks = false;
  clp.setOption("enable-verbose-blocks", "disable-verbose-blocks", &verbose_blocks, "Flag for verbose printing blocks");
  
  bool verbose = false;
  clp.setOption("enable-verbose", "disable-verbose", &verbose, "Flag for verbose printing");
  
  std::string file_input = "test.mtx";
  clp.setOption("file-input", &file_input, "Input file (MatrixMarket SPD matrix)");

  int max_concurrency = 1000000;
  clp.setOption("max-concurrency", &max_concurrency, "Max number of concurrent tasks");

  int memory_pool_grain_size = 16;
  clp.setOption("memory-pool-grain-size", &memory_pool_grain_size, "Memorypool chunk size (12 - 16)");
  
  int mkl_nthreads = 1;
  clp.setOption("mkl-nthreads", &mkl_nthreads, "MKL threads for nested parallelism");

  int nrhs = 0;
  clp.setOption("nrhs", &nrhs, "# of right hand side");

  int mb = 0;
  clp.setOption("mb", &mb, "Dense nested blocks size");

  int nb = 1;
  clp.setOption("nb", &nb, "Column block size of right hand side");

  // kokkos options are not catched here
  clp.recogniseAllOptions(false);
  clp.throwExceptions(false);
  
  Teuchos::CommandLineProcessor::EParseCommandLineReturn r_parse= clp.parse( argc, argv );

  if (r_parse == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return 0;
  if (r_parse != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL  ) return -1;
  
  int r_val = 0;
  {
    Kokkos::initialize(argc, argv);
    
    r_val = Tacho::Experimental::exampleNPT<Kokkos::OpenMP>
      (file_input, 
       max_concurrency, memory_pool_grain_size, mkl_nthreads,
       nrhs, mb, nb,
       verbose_blocks,
       verbose);
      
    Kokkos::finalize();
  }
  
  return r_val;
}
