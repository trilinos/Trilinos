#include <thread>

#include <Kokkos_Core.hpp>

#include <Kokkos_TaskPolicy.hpp>
#include <Threads/Kokkos_Threads_TaskPolicy.hpp>

#include "Teuchos_CommandLineProcessor.hpp"

#include "test_macrodef.hpp"
#include "test_suite.hpp"

/// \file test_pthread.hpp
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
    Kokkos::Threads::initialize( nthreads );
    Kokkos::Threads::print_configuration( cout , true /* detailed */ );
    
    //__TestSuiteDoUnitTests__(float,int,unsigned int,Kokkos::Serial,void);
    //__TestSuiteDoUnitTests__(float,long,unsigned long,Kokkos::Serial,void);
    
    __TestSuiteDoUnitTests__(double,int,unsigned int,Kokkos::Threads,void);
    // __TestSuiteDoUnitTests__(double,long,unsigned long,Kokkos::Serial,void);
    
    // __TestSuiteDoUnitTests__(complex<float>,int,unsigned int,Kokkos::Serial,void);
    // __TestSuiteDoUnitTests__(complex<float>,long,unsigned long,Kokkos::Serial,void);
    
    // __TestSuiteDoUnitTests__(complex<double>,int,unsigned int,Kokkos::Serial,void);
    // __TestSuiteDoUnitTests__(complex<double>,long,unsigned long,Kokkos::Serial,void);
    
    Kokkos::Threads::finalize();
  }

  string eval;
  __EVAL_STRING__(r_val, eval);
  cout << "Testing Kokkos::Threads::" << eval << endl;

  return r_val;
}
