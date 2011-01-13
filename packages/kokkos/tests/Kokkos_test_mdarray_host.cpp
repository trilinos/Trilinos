
#include <iostream>
#include <Kokkos_HostMDArrayView.hpp>
#include <Kokkos_HostDeviceFor.hpp>
#include <Kokkos_HostDeviceReduce.hpp>

#include <Kokkos_test_hex_grad.hpp>
#include <Kokkos_test_reduce.hpp>

//------------------------------------------------------------------------

int main( int argc , char ** argv )
{
  test_hex_grad< float , Kokkos::HostMap >();
  test_mom_ke< float , Kokkos::HostMap >();

  return 0 ;
}

