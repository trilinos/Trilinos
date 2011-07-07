
#include <iostream>
#include <iomanip>

#include <Kokkos_DeviceHost.hpp>
#include <Kokkos_DeviceHost_MDArrayView.hpp>
#include <Kokkos_DeviceHost_MultiVectorView.hpp>
#include <Kokkos_DeviceHost_ValueView.hpp>

#include <Kokkos_DeviceCuda.hpp>
#include <Kokkos_DeviceCuda_MDArrayView.hpp>
#include <Kokkos_DeviceCuda_MultiVectorView.hpp>
#include <Kokkos_DeviceCuda_ValueView.hpp>
#include <Kokkos_DeviceCuda_ParallelFor.hpp>
#include <Kokkos_DeviceCuda_ParallelReduce.hpp>

#include <Kokkos_DeviceCuda_macros.hpp>
#include <PerfTestHexGrad.hpp>
#include <PerfTestGramSchmidt.hpp>
#include <PerfTestDriver.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

namespace Test {

void run_test_cuda_hexgrad( int beg , int end )
{
  Kokkos::DeviceCuda::initialize();
  Test::run_test_hexgrad< Kokkos::DeviceCuda>( beg , end );
};

void run_test_cuda_gramschmidt( int beg , int end )
{
  Kokkos::DeviceCuda::initialize();
  Test::run_test_gramschmidt< Kokkos::DeviceCuda>( beg , end );
};

}




