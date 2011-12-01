
#include <iostream>
#include <iomanip>

#include <Kokkos_DeviceHost.hpp>
#include <Kokkos_DeviceHost_ValueView.hpp>
#include <Kokkos_DeviceHost_MultiVectorView.hpp>
#include <Kokkos_DeviceHost_MDArrayView.hpp>

#include <Kokkos_DeviceFerry.hpp>
#include <Kokkos_DeviceFerry_ValueView.hpp>
#include <Kokkos_DeviceFerry_MultiVectorView.hpp>
#include <Kokkos_DeviceFerry_MDArrayView.hpp>
#include <Kokkos_DeviceFerry_ParallelFor.hpp>
#include <Kokkos_DeviceFerry_ParallelReduce.hpp>

#include <Kokkos_DeviceFerry_macros.hpp>
#include <PerfTestHexGrad.hpp>
#include <PerfTestGramSchmidt.hpp>
#include <PerfTestDriver.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

namespace Test {

void run_test_ferry_hexgrad( int beg , int end )
{
  Kokkos::DeviceFerry::initialize();
  Test::run_test_hexgrad< Kokkos::DeviceFerry>( beg , end );
};

void run_test_ferry_gramschmidt( int beg , int end )
{
  Kokkos::DeviceFerry::initialize();
  Test::run_test_gramschmidt< Kokkos::DeviceFerry>( beg , end );
};

}




