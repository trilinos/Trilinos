
#include <iostream>
#include <iomanip>

#include <Kokkos_DeviceHost.hpp>
#include <Kokkos_DeviceHost_ValueView.hpp>
#include <Kokkos_DeviceHost_MultiVectorView.hpp>
#include <Kokkos_DeviceHost_MDArrayView.hpp>

#include <Kokkos_DeviceTBB.hpp>
#include <Kokkos_DeviceTBB_ValueView.hpp>
#include <Kokkos_DeviceTBB_MultiVectorView.hpp>
#include <Kokkos_DeviceTBB_MDArrayView.hpp>
#include <Kokkos_DeviceTBB_ParallelFor.hpp>
#include <Kokkos_DeviceTBB_ParallelReduce.hpp>


#include <Kokkos_DeviceTBB_macros.hpp>
#include <PerfTestHexGrad.hpp>
#include <PerfTestGramSchmidt.hpp>
#include <PerfTestDriver.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

#include <tbb/tbb.h>
#include <tbb/task_scheduler_init.h>

namespace Test {

void run_test_tbb_hexgrad( int beg , int end )
{
//  Kokkos::DeviceTBB::initialize( 12 );
  Test::run_test_hexgrad< Kokkos::DeviceTBB>( beg , end );
//  Kokkos::DeviceTBB::finalize();
};

void run_test_tbb_gramschmidt( int beg , int end )
{
//  Kokkos::DeviceTBB::initialize( 12 );
  Test::run_test_gramschmidt< Kokkos::DeviceTBB>( beg , end );
//  Kokkos::DeviceTBB::finalize();
};

}


