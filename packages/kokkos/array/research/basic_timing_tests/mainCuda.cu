#include <Kokkos_DeviceCuda.hpp>
#include <Kokkos_DeviceCuda_ValueView.hpp>
#include <Kokkos_DeviceCuda_MultiVectorView.hpp>
#include <Kokkos_DeviceCuda_MDArrayView.hpp>
#include <Kokkos_DeviceCuda_ParallelFor.hpp>
#include <Kokkos_DeviceCuda_ParallelReduce.hpp>

#include <Kokkos_DeviceCuda_macros.hpp>
#include "test_run.hpp"

__global__ void dummy_kernel(){}

int main(int argc, char* argv[]){
  unsigned int runs = 1;
  if(argc >1){
    runs = atoi(argv[1]);
  }
  test_run<double, Kokkos::DeviceCuda >("Cuda", runs);
  return 0;
}
