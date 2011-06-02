
#include <iostream>
#include <iomanip>

#include <Kokkos_DeviceHost.hpp>
#include <Kokkos_DeviceTPI.hpp>
#include <Kokkos_ValueView.hpp>
#include <Kokkos_MultiVectorView.hpp>
#include <Kokkos_MDArrayView.hpp>
#include <Kokkos_ParallelFor.hpp>
#include <Kokkos_ParallelReduce.hpp>


#include <Kokkos_DeviceTPI_macros.hpp>
#include <PerfTestHexGrad.hpp>
#include <PerfTestGramSchmidt.hpp>
#include <Kokkos_DeviceClear_macros.hpp>


namespace Test {

void run_test_tpi(int exp)
{ 
  const int thread_count = 4 ;

  Kokkos::DeviceTPI::initialize( thread_count );

  std::cout << "\"TPI(" << thread_count << ")\" , \"length\" , \"time\" , \"time/length\"" << std::endl;

  try {
    for (int i = 1; i < exp; ++i) {

      const int parallel_work_length = 1<<i;

      const double seconds = HexGrad< float , Kokkos::DeviceTPI >::test(parallel_work_length) ;

      std::cout << "\"HexGrad\" , "
                << parallel_work_length
                << " , "
                << seconds
                << " , "
                << seconds / parallel_work_length
                << std::endl ;
   }

   std::cout << "PASSED : PerfTestTPI HexGrad" << std::endl ;
  }
  catch( const std::exception & x ) {
    std::cout << "FAILED : PerfTestTPI HexGrad : " << x.what() << std::endl ;
  }

  try {
    for (int i = 1; i < exp; ++i) {

      const int parallel_work_length = 1<<i;

      const double seconds = ModifiedGramSchmidt< float , Kokkos::DeviceTPI >::test(parallel_work_length, 32 ) ;

      std::cout << "\"GramSchmidt(32)\" , "
                << parallel_work_length
                << " , "
                << seconds
                << " , "
                << seconds / parallel_work_length
                << std::endl ;
   }

   std::cout << "PASSED : PerfTestTPI GramSchmidt" << std::endl ;
  }
  catch( const std::exception & x ) {
    std::cout << "FAILED : PerfTestTPI GramSchmidt : " << x.what() << std::endl ;
  }
}

}

