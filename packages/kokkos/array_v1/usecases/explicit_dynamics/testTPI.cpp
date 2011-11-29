#include <iostream>
#include <iomanip>
#include <sys/time.h>

#include <Kokkos_DeviceHost.hpp>
#include <Kokkos_DeviceHost_ValueView.hpp>
#include <Kokkos_DeviceHost_MultiVectorView.hpp>
#include <Kokkos_DeviceHost_MDArrayView.hpp>
#include <Kokkos_DeviceHost_ParallelFor.hpp>
#include <Kokkos_DeviceHost_ParallelReduce.hpp>

#include <Kokkos_DeviceTPI.hpp>
#include <Kokkos_DeviceTPI_ValueView.hpp>
#include <Kokkos_DeviceTPI_MultiVectorView.hpp>
#include <Kokkos_DeviceTPI_MDArrayView.hpp>
#include <Kokkos_DeviceTPI_ParallelFor.hpp>
#include <Kokkos_DeviceTPI_ParallelReduce.hpp>

#include <Kokkos_DeviceTPI_macros.hpp>
#include <explicit_dynamics_app.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

namespace test{

  void test_TPI( int beg, int end, int runs, int threads){

    Kokkos::DeviceTPI::initialize( threads );

    explicit_dynamics::driver<float,Kokkos::DeviceTPI>("TPI float", beg, end, runs);
    explicit_dynamics::driver<double,Kokkos::DeviceTPI>("TPI double", beg, end, runs);

    Kokkos::DeviceTPI::finalize();

  }//test_host

}// namespace


