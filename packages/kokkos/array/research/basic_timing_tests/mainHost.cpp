#include <Kokkos_DeviceHost.hpp>
#include <Kokkos_DeviceHost_ValueView.hpp>
#include <Kokkos_DeviceHost_MultiVectorView.hpp>
#include <Kokkos_DeviceHost_MDArrayView.hpp>
#include <Kokkos_DeviceHost_ParallelFor.hpp>
#include <Kokkos_DeviceHost_macros.hpp>
#include <Kokkos_DeviceHost_ParallelReduce.hpp>

#include "test_run.hpp"

int main(int argc, char* argv[]){
  unsigned int runs = 1;
  if(argc >1){
    runs = atoi(argv[1]);
  }
  test_run<double, Kokkos::DeviceHost >("Host", runs);
  return 0;
}
