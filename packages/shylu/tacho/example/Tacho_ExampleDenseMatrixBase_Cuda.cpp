#include <Kokkos_Core.hpp>
#include <Kokkos_Cuda.hpp>
#include "Teuchos_CommandLineProcessor.hpp"

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

typedef Kokkos::Cuda exec_space;
typedef Kokkos::Impl::is_space<exec_space>::host_mirror_space::execution_space host_space;  

#include "Tacho_ExampleDenseMatrixBase.hpp"

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

  int mmin = 1000;
  clp.setOption("mmin", &mmin, "C(mmin,mmin)");

  int mmax = 8000;
  clp.setOption("mmax", &mmax, "C(mmax,mmax)");

  int minc = 1000;
  clp.setOption("minc", &minc, "Increment of m");

  clp.recogniseAllOptions(true);
  clp.throwExceptions(false);

  Teuchos::CommandLineProcessor::EParseCommandLineReturn r_parse= clp.parse( argc, argv );

  if (r_parse == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return 0;
  if (r_parse != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL  ) return -1;

  int r_val = 0;
  {
    exec_space::initialize();
    exec_space::print_configuration(std::cout, true);

    host_space::initialize(nthreads, numa, core_per_numa);
    host_space::print_configuration(std::cout, true);

    r_val = exampleDenseMatrixBase<exec_space>
      (mmin, mmax, minc, 
       verbose);
    
    exec_space::finalize();
    host_space::finalize();
  }
  
  return r_val;
}
