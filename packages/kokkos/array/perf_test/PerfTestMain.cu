
#include <iostream>
#include <iomanip>

#include <Kokkos_DeviceHost.hpp>
#include <Kokkos_DeviceTPI.hpp>

#include <Kokkos_MDArrayView.hpp>
#include <Kokkos_MultiVectorView.hpp>
#include <Kokkos_ValueView.hpp>
#include <Kokkos_ParallelFor.hpp>
#include <Kokkos_ParallelReduce.hpp>

#include <Kokkos_DeviceHost_macros.hpp>
#include <PerfTestHexGrad.hpp>
#include <PerfTestGramSchmidt.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

#include <Kokkos_DeviceCuda_macros.hpp>
#include <PerfTestHexGrad.hpp>
#include <PerfTestGramSchmidt.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

//------------------------------------------------------------------------

template < class T >
void run_test(int exp)
{
  for (int i = 1; i < exp; ++i) {

    const int parallel_work_length = 1<<i;

    T::test(parallel_work_length) ;
  }
}

int main( int argc , char ** argv )
{
  const int exp = 7 ;

  run_test< HexGrad< float , Kokkos::DeviceHost > >(exp);
  run_test< HexGrad< float , Kokkos::DeviceCuda > >(exp);

  return 0 ;
}

