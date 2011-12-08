#include <iostream>
#include <iomanip>
#include <sys/time.h>

#include <Kokkos_Value.hpp>
#include <Kokkos_MultiVector.hpp>
#include <Kokkos_MDArray.hpp>

#include <Kokkos_Host.hpp>
#include <Kokkos_Cuda.hpp>

#include <Kokkos_Cuda_macros.hpp>
#include <explicit_dynamics_app.hpp>
#include <Kokkos_Clear_macros.hpp>

__global__ void dummy_kernel(){}

namespace Test{

void test_Cuda(int beg, int end, int runs){

  Kokkos::Cuda::initialize();

  cudaFuncSetCacheConfig(dummy_kernel, cudaFuncCachePreferL1);
  dummy_kernel<<<1, 1>>>();

  std::cout << "Kokkos Cuda: " << std::endl;


  explicit_dynamics::driver<float,Kokkos::Cuda>("Cuda-float", beg, end, runs);
  explicit_dynamics::driver<double,Kokkos::Cuda>("Cuda-double", beg, end, runs);

  Kokkos::Cuda::finalize();
}

}// namespace

