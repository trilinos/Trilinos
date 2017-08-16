#include <Kokkos_Core.hpp>

#include <Kokkos_Threads.hpp>
#include <Threads/Kokkos_Threads_TaskPolicy.hpp>

#include "Teuchos_CommandLineProcessor.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

typedef Kokkos::Threads exec_space;

#include "example_dense_trsm_by_blocks.hpp"

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

  int max_concurrency = 1024;
  clp.setOption("max-concurrency", &max_concurrency, "Max number of concurrent tasks");

  int max_task_dependence = 3;
  clp.setOption("max-task-dependence", &max_task_dependence, "Max number of task dependence");

  int team_size = 1;
  clp.setOption("team-size", &team_size, "Team size");

  int mkl_nthreads = 1;
  clp.setOption("mkl-nthreads", &mkl_nthreads, "MKL threads for nested parallelism");

  bool verbose = false;
  clp.setOption("enable-verbose", "disable-verbose", &verbose, "Flag for verbose printing");

  int mmin = 1000;
  clp.setOption("mmin", &mmin, "A(mmin,mmin)");

  int mmax = 8000;
  clp.setOption("mmax", &mmax, "A(mmax,mmax)");

  int minc = 1000;
  clp.setOption("minc", &minc, "Increment of m");

  int k = 1024;
  clp.setOption("k", &k, "B(mmax,k)");

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
    exec_space::print_configuration(cout, true);

    cout << "DenseTrsmByBlocks:: ConjTranspose" << endl;
    r_val = exampleDenseTrsmByBlocks
      <value_type,ordinal_type,size_type,exec_space,void>
      (mmin, mmax, minc, k, mb,
       max_concurrency, max_task_dependence, team_size, mkl_nthreads,
       check,
       verbose);
    
    exec_space::finalize();
  }

  return r_val;
}
