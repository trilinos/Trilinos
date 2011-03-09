

#include <iostream>
#include <TPI.h>
#include <Kokkos_HostTPIMultiVectorView.hpp>
#include <Kokkos_HostTPIFor.hpp>
#include <Kokkos_HostTPIReduce.hpp>
#include <Kokkos_HostMath.hpp>

#include <Kokkos_test_gram_schmidt.hpp>

//------------------------------------------------------------------------

int main( int argc , char ** argv )
{
  enum { LENGTH = 100000 , COUNT = 10 };

  TPI_Init( 4 );

  test_modified_gram_schmidt<float,Kokkos::HostTPI>( LENGTH , COUNT );

  return 0 ;
}

