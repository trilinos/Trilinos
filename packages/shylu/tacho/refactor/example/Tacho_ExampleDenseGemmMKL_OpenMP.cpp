#include <Kokkos_Core.hpp>

#include <Kokkos_OpenMP.hpp>

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_BLAS.hpp"

#include "ShyLUTacho_config.h"

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

typedef Kokkos::OpenMP exec_space;

#include "Tacho_ExampleDenseGemmMKL.hpp"

using namespace Tacho;

int main (int argc, char *argv[]) {
  
  Teuchos::CommandLineProcessor clp;
  clp.setDocString("This example program measure the performance of dense Gemm on Kokkos::Threads execution space.\n");
  
  int nthreads = 0;
  clp.setOption("nthreads", &nthreads, "Number of threads");

  int numa = 0;
  clp.setOption("numa", &numa, "Number of numa node");

  int core_per_numa = 0;
  clp.setOption("core-per-numa", &core_per_numa, "Number of cores per numa node");
  
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

  clp.recogniseAllOptions(true);
  clp.throwExceptions(false);

  Teuchos::CommandLineProcessor::EParseCommandLineReturn r_parse= clp.parse( argc, argv );

  if (r_parse == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return 0;
  if (r_parse != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL  ) return -1;

  int r_val = 0;
  {
    exec_space::initialize(nthreads, numa, core_per_numa);

#ifdef HAVE_SHYLUTACHO_MKL
    std::cout << "DenseGemmByBlocks:: NoTranspose, NoTranspose" << std::endl;
    mkl_set_num_threads(nthreads);
    r_val = exampleDenseGemmMKL<exec_space>
      (mmin, mmax, minc, k,
       verbose);
#else
    TACHO_TEST_FOR_ABORT( true, MSG_NOT_HAVE_PACKAGE("MKL") );
#endif

    exec_space::finalize();
  }

  return r_val;
}
