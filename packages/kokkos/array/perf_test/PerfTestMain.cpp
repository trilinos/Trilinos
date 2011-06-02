
#include <iostream>
#include <iomanip>

#include <Kokkos_DeviceHost.hpp>

#include <Kokkos_MDArrayView.hpp>
#include <Kokkos_MultiVectorView.hpp>
#include <Kokkos_ValueView.hpp>
#include <Kokkos_ParallelFor.hpp>
#include <Kokkos_ParallelReduce.hpp>

#include <Kokkos_DeviceHost_macros.hpp>
#include <PerfTestHexGrad.hpp>
#include <PerfTestGramSchmidt.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

//------------------------------------------------------------------------
namespace Test {

void run_test_host(int exp)
{ 
  try {
    for (int i = 1; i < exp; ++i) {

      const int parallel_work_length = 1<<i;

      double seconds = HexGrad< float , Kokkos::DeviceHost >::test(parallel_work_length) ;

      std::cout << "\"Host HexGrad\" , "
                << parallel_work_length
                << " , "
                << seconds
                << std::endl ;
    }
    std::cout << "PASSED : PerfTestHost HexGrad" << std::endl ;
  }
  catch( const std::exception & x ) {
    std::cout << "FAILED : PerfTestHost HexGrad : " << x.what() << std::endl ;
  }

  try {
    for (int i = 1; i < exp; ++i) {

      const int parallel_work_length = 1<<i;

      const double seconds = ModifiedGramSchmidt< float , Kokkos::DeviceHost >::test(parallel_work_length, 32 ) ;

      std::cout << "\"GramSchmidt(32)\" , "
                << parallel_work_length
                << " , "
                << seconds
                << " , "
                << seconds / parallel_work_length
                << std::endl ;
   }

   std::cout << "PASSED : PerfTestHost GramSchmidt" << std::endl ;
  }
  catch( const std::exception & x ) {
    std::cout << "FAILED : PerfTestHost GramSchmidt : " << x.what() << std::endl ;
  }

}

}

namespace Test {
  void run_test_tpi(int exp);
  void run_test_cuda(int exp);
}


int main( int argc , char ** argv )
{
  const int exp = 20 ;

  Test::run_test_host(exp);
  Test::run_test_tpi(exp);
  Test::run_test_cuda(exp);
  return 0 ;
}

