
#include <iostream>
#include <iomanip>

#include <Kokkos_DeviceHost.hpp>
#include <Kokkos_DeviceCuda.hpp>

#include <Kokkos_MDArrayView.hpp>
#include <Kokkos_MultiVectorView.hpp>
#include <Kokkos_ValueView.hpp>
#include <Kokkos_ParallelFor.hpp>
#include <Kokkos_ParallelReduce.hpp>

#include <Kokkos_DeviceCuda_macros.hpp>
#include <PerfTestHexGrad.hpp>
#include <PerfTestGramSchmidt.hpp>


//------------------------------------------------------------------------
namespace Test {
void run_test_cuda(int exp)
{
#if defined( KOKKOS_MACRO_DEVICE_FUNCTION )
  try {
    for (int i = 1; i < exp; ++i) {

      const int parallel_work_length = 1<<i;

      HexGrad< float , Kokkos::DeviceCuda >::test(parallel_work_length) ;
    }

    std::cout << "PASSED : PerfTestCuda" << std::endl ;
  }
  catch( const std::exception & x ) {
    std::cout << "FAILED : PerfTestCuda : " << x.what() << std::endl ;
  }
#else 
  std::cout << "PASSED : SKIPPED UnitTestCuda - NO DEVICE CUDA" << std::endl ;
#endif
}
}




