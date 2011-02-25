

#include <iostream>
#include <TPI.h>
#include <Kokkos_HostMDArrayView.hpp>
#include <Kokkos_HostMVectorView.hpp>
#include <Kokkos_HostDeviceFor.hpp>
#include <Kokkos_HostDeviceReduce.hpp>
#include <Kokkos_HostMath.hpp>

#include <Kokkos_test_gram_schmidt.hpp>

//------------------------------------------------------------------------

int main( int argc , char ** argv )
{
  enum { LENGTH = 100000 , COUNT = 10 };

  TPI_Init( 4 );

  test_modified_gram_schmidt<float,Kokkos::HostMap>( LENGTH , COUNT );

  return 0 ;
}

