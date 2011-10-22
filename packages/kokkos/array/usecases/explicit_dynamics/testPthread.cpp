#include <iostream>
#include <iomanip>
#include <sys/time.h>

#include <Kokkos_DeviceHost.hpp>
#include <Kokkos_DeviceHost_ValueView.hpp>
#include <Kokkos_DeviceHost_MultiVectorView.hpp>
#include <Kokkos_DeviceHost_MDArrayView.hpp>
#include <Kokkos_DeviceHost_ParallelFor.hpp>
#include <Kokkos_DeviceHost_ParallelReduce.hpp>

#include <Kokkos_DevicePthread.hpp>
#include <Kokkos_DevicePthread_ValueView.hpp>
#include <Kokkos_DevicePthread_MultiVectorView.hpp>
#include <Kokkos_DevicePthread_MDArrayView.hpp>
#include <Kokkos_DevicePthread_ParallelFor.hpp>
#include <Kokkos_DevicePthread_ParallelReduce.hpp>

#include <Kokkos_DevicePthread_macros.hpp>
#include <explicit_dynamics_app.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

namespace test{

  void test_Pthread( int beg, int end, int runs, int threads){

    Kokkos::DevicePthread::initialize( threads );

    explicit_dynamics::driver<float,Kokkos::DevicePthread>("Pthread float", beg, end, runs);
    explicit_dynamics::driver<double,Kokkos::DevicePthread>("Pthread double", beg, end, runs);

    Kokkos::DevicePthread::finalize();

  }//test_host

}// namespace


