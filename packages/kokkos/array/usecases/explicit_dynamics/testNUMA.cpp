#include <iostream>
#include <iomanip>
#include <sys/time.h>

#include <Kokkos_DeviceHost.hpp>
#include <Kokkos_DeviceHost_ValueView.hpp>
#include <Kokkos_DeviceHost_MultiVectorView.hpp>
#include <Kokkos_DeviceHost_MDArrayView.hpp>
#include <Kokkos_DeviceHost_ParallelFor.hpp>
#include <Kokkos_DeviceHost_ParallelReduce.hpp>

#include <Kokkos_ValueView.hpp>
#include <Kokkos_MultiVectorView.hpp>
#include <Kokkos_MDArrayView.hpp>
#include <Kokkos_ParallelFor.hpp>
#include <Kokkos_ParallelReduce.hpp>
#include <Kokkos_DeviceNUMA.hpp>

#include <Kokkos_DeviceNUMA_macros.hpp>
#include <explicit_dynamics_app.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

namespace test{

  void test_NUMA( int beg, int end, int runs){

    std::cout << std::endl
              << "Kokkos::DeviceNUMA::detect_core_count() = " 
              <<  Kokkos::DeviceNUMA::detect_core_count()
              << std::endl ;

    Kokkos::DeviceNUMA::initialize( Kokkos::DeviceNUMA::DETECT_AND_USE_ALL_CORES );

    explicit_dynamics::driver<float,Kokkos::DeviceNUMA>("NUMA float", beg, end, runs);
    explicit_dynamics::driver<double,Kokkos::DeviceNUMA>("NUMA double", beg, end, runs);

    Kokkos::DeviceNUMA::finalize();

  }//test_host

}// namespace


