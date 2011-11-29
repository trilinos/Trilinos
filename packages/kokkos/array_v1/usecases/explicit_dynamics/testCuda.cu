#include <iostream>
#include <iomanip>
#include <sys/time.h>

#include <Kokkos_DeviceHost.hpp>
#include <Kokkos_DeviceHost_ValueView.hpp>
#include <Kokkos_DeviceHost_MultiVectorView.hpp>
#include <Kokkos_DeviceHost_MDArrayView.hpp>
#include <Kokkos_DeviceHost_ParallelFor.hpp>
#include <Kokkos_DeviceHost_ParallelReduce.hpp>

#include <Kokkos_DeviceCuda.hpp>
#include <Kokkos_DeviceCuda_ValueView.hpp>
#include <Kokkos_DeviceCuda_MultiVectorView.hpp>
#include <Kokkos_DeviceCuda_MDArrayView.hpp>
#include <Kokkos_DeviceCuda_ParallelFor.hpp>
#include <Kokkos_DeviceCuda_ParallelReduce.hpp>

#include <Kokkos_DeviceCuda_macros.hpp>
#include <explicit_dynamics_app.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

__global__ void dummy_kernel(){}

namespace test{

	void test_Cuda(int beg, int end, int runs){

		cudaFuncSetCacheConfig(dummy_kernel, cudaFuncCachePreferL1);
		dummy_kernel<<<1, 1>>>();

		std::cout << "Kokkos Cuda: " << std::endl;

		Kokkos::DeviceCuda::initialize();

    explicit_dynamics::driver<float,Kokkos::DeviceCuda>("Cuda float", beg, end, runs);
    explicit_dynamics::driver<double,Kokkos::DeviceCuda>("Cuda double", beg, end, runs);
	}

}// namespace

