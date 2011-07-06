
#include <iostream>
#include <iomanip>

#include <Kokkos_DeviceHost.hpp>
#include <Kokkos_DeviceHost_ValueView.hpp>
#include <Kokkos_DeviceHost_MultiVectorView.hpp>
#include <Kokkos_DeviceHost_MDArrayView.hpp>

#include <Kokkos_DeviceTPI.hpp>
#include <Kokkos_DeviceTPI_ValueView.hpp>
#include <Kokkos_DeviceTPI_MultiVectorView.hpp>
#include <Kokkos_DeviceTPI_MDArrayView.hpp>

#include <Kokkos_ParallelFor.hpp>
#include <Kokkos_ParallelReduce.hpp>


#include <Kokkos_DeviceTPI_macros.hpp>
#include <PerfTestHexGrad.hpp>
#include <PerfTestGramSchmidt.hpp>
#include <PerfTestDriver.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

namespace Test {

void run_test_tpi_hexgrad( int beg , int end )
{
  Kokkos::DeviceTPI::initialize( 12 );
  Test::run_test_hexgrad< Kokkos::DeviceTPI>( beg , end );
};

void run_test_tpi_gramschmidt( int beg , int end )
{
  Kokkos::DeviceTPI::initialize( 12 );
  Test::run_test_gramschmidt< Kokkos::DeviceTPI>( beg , end );
};

}


