#include <iostream>
#include <iomanip>
#include <sys/time.h>

#include <Kokkos_DeviceHost.hpp>
#include <Kokkos_DeviceHost_ValueView.hpp>
#include <Kokkos_DeviceHost_MultiVectorView.hpp>
#include <Kokkos_DeviceHost_MDArrayView.hpp>
#include <Kokkos_DeviceHost_ParallelFor.hpp>
#include <Kokkos_DeviceHost_ParallelReduce.hpp>

#include <Kokkos_DeviceTBB.hpp>
#include <Kokkos_DeviceTBB_ValueView.hpp>
#include <Kokkos_DeviceTBB_MultiVectorView.hpp>
#include <Kokkos_DeviceTBB_MDArrayView.hpp>
#include <Kokkos_DeviceTBB_ParallelFor.hpp>
#include <Kokkos_DeviceTBB_ParallelReduce.hpp>

#include <Kokkos_DeviceTBB_macros.hpp>
#include <explicit_dynamics_app.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

namespace test{

  void test_TBB( int beg, int end, int runs){

    explicit_dynamics::driver<float,Kokkos::DeviceTBB>("TBB float", beg, end, runs);
    explicit_dynamics::driver<double,Kokkos::DeviceTBB>("TBB double", beg, end, runs);
	}//test_host

  }//test_host

}// namespace


