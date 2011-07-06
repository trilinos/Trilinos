
#include <iostream>
#include <iomanip>

#include <Kokkos_DeviceHost.hpp>

#include <Kokkos_MDArrayView.hpp>
#if defined( KOKKOS_DEVICE_HOST )
#include <Kokkos_DeviceHost_MDArrayView.hpp>
#endif
#if defined( KOKKOS_DEVICE_TPI )
#include <Kokkos_DeviceTPI_MDArrayView.hpp>
#endif
#if defined( KOKKOS_DEVICE_CUDA )
#include <Kokkos_DeviceCuda_MDArrayView.hpp>
#endif


#include <Kokkos_MultiVectorView.hpp>
#if defined( KOKKOS_DEVICE_HOST )
#include <Kokkos_DeviceHost_MultiVectorView.hpp>
#endif
#if defined( KOKKOS_DEVICE_TPI )
#include <Kokkos_DeviceTPI_MultiVectorView.hpp>
#endif
#if defined( KOKKOS_DEVICE_CUDA )
#include <Kokkos_DeviceCuda_MultiVectorView.hpp>
#endif


#include <Kokkos_ValueView.hpp>
#if defined( KOKKOS_DEVICE_HOST )
#include <Kokkos_DeviceHost_ValueView.hpp>
#endif
#if defined( KOKKOS_DEVICE_TPI )
#include <Kokkos_DeviceTPI_ValueView.hpp>
#endif
#if defined( KOKKOS_DEVICE_CUDA )
#include <Kokkos_DeviceCuda_ValueView.hpp>
#endif



#include <Kokkos_ParallelFor.hpp>
#include <Kokkos_ParallelReduce.hpp>

#include <Kokkos_DeviceHost_macros.hpp>
#include <PerfTestHexGrad.hpp>
#include <PerfTestGramSchmidt.hpp>
#include <PerfTestDriver.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

//------------------------------------------------------------------------

namespace Test {

void run_test_host_hexgrad( int beg , int end )
{ Test::run_test_hexgrad< Kokkos::DeviceHost>( beg , end ); }

void run_test_host_gramschmidt( int beg , int end )
{ Test::run_test_gramschmidt< Kokkos::DeviceHost>( beg , end ); }

void run_test_tpi_hexgrad(int,int);
void run_test_tpi_gramschmidt(int,int);

void run_test_cuda_hexgrad(int,int);
void run_test_cuda_gramschmidt(int,int);

}

int main( int argc , char ** argv )
{
  Test::run_test_host_hexgrad( 10 , 20 );
  Test::run_test_tpi_hexgrad(  10 , 24 );
  Test::run_test_cuda_hexgrad( 10 , 24 );

 Test::run_test_host_gramschmidt( 10 , 20 );
  Test::run_test_tpi_gramschmidt(  10 , 24 );
  Test::run_test_cuda_gramschmidt( 10 , 24 );

  return 0 ;
}

