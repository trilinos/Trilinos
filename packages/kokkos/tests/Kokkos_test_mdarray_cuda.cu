
#include <iostream>

#include <Kokkos_CudaMDArrayView.hpp>
#include <Kokkos_CudaDeviceFor.hpp>

#include <Kokkos_test_hex_grad.hpp>

//------------------------------------------------------------------------

int main( int argc , char ** argv )
{
  test_hex_grad< float , Kokkos::CudaMap >();

  return 0 ;
}

