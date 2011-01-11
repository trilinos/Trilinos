
#include <iostream>
#include <Kokkos_HostMDArrayView.hpp>
#include <Kokkos_HostDeviceFor.hpp>

#include <Kokkos_test_hex_grad.hpp>

//------------------------------------------------------------------------

int main( int argc , char ** argv )
{
  test_hex_grad< float , Kokkos::HostMap >();

  return 0 ;
}

