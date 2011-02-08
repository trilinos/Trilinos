
#include <iostream>

#include <Kokkos_CudaMDArrayView.hpp>
#include <Kokkos_CudaDeviceFor.hpp>
#include <Kokkos_CudaDeviceReduce.hpp>

#include <Kokkos_test_hex_grad.hpp>
#include <Kokkos_test_reduce.hpp>

//------------------------------------------------------------------------

int main( int argc , char ** argv )
{
  test_hex_grad< float , Kokkos::CudaMap >();

  TestMomKe< float , Kokkos::CudaMap > test ;

  return 0 ;
}

