#include <iostream>
#include <iomanip>
#include <sys/time.h>

#include <Kokkos_Value.hpp>
#include <Kokkos_MultiVector.hpp>
#include <Kokkos_MDArray.hpp>

#include <Kokkos_Host.hpp>

#include <Kokkos_Host_macros.hpp>
#include <explicit_dynamics_app.hpp>
#include <Kokkos_Clear_macros.hpp>

namespace Test{

void test_Host( int beg, int end, int runs, int threads){

  if ( 0 < threads ) {
    Kokkos::Host::initialize( Kokkos::Host::SetThreadCount( threads ) );
  }
  else {
    Kokkos::Host::initialize( Kokkos::Host::DetectAndUseAllCores() );
    threads = Kokkos::Host::detect_core_count();
  }

  std::cout << "\"Host with threads = \" , " << threads << std::endl ;

  explicit_dynamics::driver<float,Kokkos::Host>("Host-float", beg, end, runs);
  explicit_dynamics::driver<double,Kokkos::Host>("Host-double", beg, end, runs);

  Kokkos::Host::finalize();
}//test_host

}// namespace


