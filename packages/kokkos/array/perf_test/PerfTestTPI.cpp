
#include <iostream>
#include <iomanip>

#include <Kokkos_DeviceTPI.hpp>


#include <Kokkos_DeviceTPI_macros.hpp>
#include <PerfTestHexGrad.hpp>
#include <PerfTestGramSchmidt.hpp>
#include <Kokkos_DeviceClear_macros.hpp>


namespace Test {
  void run_test_tpi(int exp)
  { 
    try {
      for (int i = 1; i < exp; ++i) {

      const int parallel_work_length = 1<<i;

      HexGrad< float , Kokkos::DeviceTPI >::test(parallel_work_length) ;
      
      std::cout << "PASSED : PerfTestTPI" << std::endl ;
  }
  catch( const std::exception & x ) {
    std::cout << "FAILED : PerfTestTPI : " << x.what() << std::endl ;
    }
  }
}