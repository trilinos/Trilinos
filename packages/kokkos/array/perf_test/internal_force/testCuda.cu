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
#include <internal_force_driver.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

__global__ void dummy_kernel(){}

namespace test{

	void test_Cuda(int beg, int end, int runs){

		cudaFuncSetCacheConfig(dummy_kernel, cudaFuncCachePreferL1);
		dummy_kernel<<<1, 1>>>();

		std::cout << "Kokkos Cuda: " << std::endl;

		Kokkos::DeviceCuda::initialize();

		for(int i = beg; i < end; i++){

			double min = 100000;

			int x = 10 + 5 * i;
			int y = 10 + 5 * i;
			int z = 10 + 5 * i;

			int n = x * y * z;

			for(int j = 0; j < runs; j++){

				double time = internal_force_test<double, Kokkos::DeviceCuda>( x, y, z );

				if ( 0 == j || time < min)
					min = time;
			}

			std::cout << 	std::setw(8) << n << ", " << 
							std::setw(8) << 1000 * min << ", " << 
							std::setw(8) << 1000 * min / n << std::endl;

		}//for
	}

}// namespace

