#include <Kokkos_Core.hpp>

#include "Teuchos_CommandLineProcessor.hpp"
#include "ShyLUTacho_config.h"

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

typedef Kokkos::OpenMP exec_space;

#include "Tacho_ExampleDenseHerkByBlocks.hpp"

using namespace Tacho;

int main (int argc, char *argv[]) {

  Teuchos::CommandLineProcessor clp;
  clp.setDocString("This example program measure the performance of dense Herk on Kokkos::Threads execution space.\n");

  int nthreads = 0;
  clp.setOption("nthreads", &nthreads, "Number of threads");

  int numa = 0;
  clp.setOption("numa", &numa, "Number of numa node");

  int core_per_numa = 0;
  clp.setOption("core-per-numa", &core_per_numa, "Number of cores per numa node");

  int max_concurrency = 250000;
  clp.setOption("max-concurrency", &max_concurrency, "Max number of concurrent tasks");

  int memory_pool_grain_size = 16;
  clp.setOption("memory-pool-grain-size", &memory_pool_grain_size, "Memorypool chunk size (12 - 16)");

  int mkl_nthreads = 1;
  clp.setOption("mkl-nthreads", &mkl_nthreads, "MKL threads for nested parallelism");

  bool verbose = false;
  clp.setOption("enable-verbose", "disable-verbose", &verbose, "Flag for verbose printing");

  int mmin = 1000;
  clp.setOption("mmin", &mmin, "C(mmin,mmin)");

  int mmax = 8000;
  clp.setOption("mmax", &mmax, "C(mmax,mmax)");

  int minc = 1000;
  clp.setOption("minc", &minc, "Increment of m");

  int k = 1024;
  clp.setOption("k", &k, "A(mmax,k) or A(k,mmax) according to transpose flags");

  int mb = 256;
  clp.setOption("mb", &mb, "Blocksize");

  bool check = true;
  clp.setOption("enable-check", "disable-check", &check, "Flag for check solution");

  clp.recogniseAllOptions(true);
  clp.throwExceptions(false);

  Teuchos::CommandLineProcessor::EParseCommandLineReturn r_parse= clp.parse( argc, argv );

  if (r_parse == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return 0;
  if (r_parse != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL  ) return -1;

  int r_val = 0;
  {
    exec_space::initialize(nthreads, numa, core_per_numa);
    
    std::cout << std::endl << "DenseHerkByBlocks:: Upper, ConjTranspose, Variant::One (external)" << std::endl;
    r_val = exampleDenseHerkByBlocks
      <Uplo::Upper,Trans::ConjTranspose,Variant::One,exec_space>
      (mmin, mmax, minc, k, mb,
       max_concurrency, memory_pool_grain_size, mkl_nthreads,
       check,
       verbose);

    exec_space::finalize();
  }
  
  return r_val;
}

