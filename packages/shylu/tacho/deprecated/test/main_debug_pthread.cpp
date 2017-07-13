#include <thread>

#include <Kokkos_Core.hpp>

#include <Kokkos_TaskPolicy.hpp>
#include <Threads/Kokkos_Threads_TaskPolicy.hpp>

#include "Teuchos_CommandLineProcessor.hpp"

#include "test_macrodef.hpp"
#include "test_suite.hpp"

#include "test_tri_solve_by_blocks_debug.hpp"

/// \file test_debug_pthread.hpp
/// \brief Test pthread execution space
/// \author Kyungjoo Kim (kyukim@sandia.gov)

using namespace std;
using namespace Tacho;

int g_funct_counter = 0;

int main(int argc, char *argv[]) {
  int r_val = 0;

  Teuchos::CommandLineProcessor clp;

  int nthreads = 1;
  clp.setOption("nthreads", &nthreads, "Number of threads");

  int numa = 0;
  clp.setOption("numa", &numa, "Number of numa node");

  int core_per_numa = 0;
  clp.setOption("core-per-numa", &core_per_numa, "Number of cores per numa node");

  int max_task_dependence = 10;
  clp.setOption("max-task-dependence", &max_task_dependence, "Max number of task dependence");

  int team_size = 1;
  clp.setOption("team-size", &team_size, "Team size");

  clp.recogniseAllOptions(true);
  clp.throwExceptions(false);

  Teuchos::CommandLineProcessor::EParseCommandLineReturn r_parse= clp.parse( argc, argv );

  if (r_parse != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL  ) {
    cout << "Testing Kokkos::Qthread:: Failed in parsing command line input" << endl;
    return -1;
  }
  if (r_parse == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
    return 0;
  }

  unsigned threads_count = 0;
  if (Kokkos::hwloc::available())  {
    const unsigned numa_count       = Kokkos::hwloc::get_available_numa_count();
    const unsigned cores_per_numa   = Kokkos::hwloc::get_available_cores_per_numa();
    const unsigned threads_per_core = Kokkos::hwloc::get_available_threads_per_core();
    
    const unsigned one = 1u;
    threads_count = max(one, numa_count)*max(one, cores_per_numa)*max(one, threads_per_core);
 
    cout << " = Kokkos::hwloc = " << endl
         << "NUMA count       = " << numa_count << endl
         << "Cores per NUMA   = " << cores_per_numa << endl
         << "Threads per core = " << threads_per_core << endl
         << "Threads count    = " << threads_count << endl;
  } else {
    threads_count = thread::hardware_concurrency();

    cout << " = std::thread::hardware_concurrency = " << endl
         << "Threads count    = " << threads_count << endl;
  }

  if (static_cast<unsigned int>(nthreads) > threads_count) {
    ++r_val;
    cout << "Testing Kokkos::Threads:: Failed that the given nthreads is greater than the number of threads counted" << endl;
  } else {
    Kokkos::Threads::initialize( nthreads, numa, core_per_numa );
    Kokkos::Threads::print_configuration( cout , true /* detailed */ );

    const int blk_cnt = 6, blks[blk_cnt] = { 1, 2, 4, 8, 12, 16 };
    const int nrhs_cnt = 6, nrhs[nrhs_cnt] = { 1, 2, 4, 8, 12, 16 };

    r_val += testTriSolveByBlocksDebug<double,int,int,Kokkos::Threads,void>
      ("mm_crs_input.mtx", team_size, max_task_dependence, blks[0], nrhs[2]);
    
    Kokkos::Threads::finalize();
  }

  string eval;
  __EVAL_STRING__(r_val, eval);
  cout << "Testing Kokkos::Threads::" << eval << endl;

  return r_val;
}
