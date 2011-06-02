
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

      double seconds = HexGrad< float , Kokkos::DeviceCuda >::test(parallel_work_length) ;

      std::cout << "\"Cuda HexGrad\" , "
                << parallel_work_length
                << " , "
                << seconds
                << std::endl ;
    }
    std::cout << "PASSED : PerfTestCuda HexGrad" << std::endl ;
  }
  catch( const std::exception & x ) {
    std::cout << "FAILED : PerfTestCuda HexGrad : " << x.what() << std::endl ;
  }

  try {
    for (int i = 1; i < exp; ++i) {

      const int parallel_work_length = 1<<i;

      const double seconds = ModifiedGramSchmidt< float , Kokkos::DeviceCuda >::test(parallel_work_length, 32 ) ;

      std::cout << "\"GramSchmidt(32)\" , "
                << parallel_work_length
                << " , "
                << seconds
                << " , "
                << seconds / parallel_work_length
                << std::endl ;
   }

   std::cout << "PASSED : PerfTestCuda GramSchmidt" << std::endl ;
  }
  catch( const std::exception & x ) {
    std::cout << "FAILED : PerfTestCuda GramSchmidt : " << x.what() << std::endl ;
  }
#else 
  std::cout << "PASSED : SKIPPED PerfTestCuda - NO DEVICE CUDA" << std::endl ;
#endif
}
}




