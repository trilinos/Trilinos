#include <iostream>
#include <iomanip>
#include <sys/time.h>

#include <Kokkos_DeviceHost.hpp>
#include <Kokkos_DeviceHost_ValueView.hpp>
#include <Kokkos_DeviceHost_MultiVectorView.hpp>
#include <Kokkos_DeviceHost_MDArrayView.hpp>
#include <Kokkos_DeviceHost_ParallelFor.hpp>
#include <Kokkos_DeviceHost_ParallelReduce.hpp>

#include <Kokkos_DeviceHost_macros.hpp>
#include <internal_force_driver.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

namespace test{

	void test_Host( int beg, int end, int runs){

		std::cout << "Kokkos Host: " << std::endl;

		for(int i = beg; i < end; i++){

			double min = 100000;

			int x = 10 + 5 * i;
			int y = 10 + 5 * i;
			int z = 10 + 5 * i;

			int n = x * y * z;

			for(int j = 0; j < runs; j++){

				double time = internal_force_test<double, Kokkos::DeviceHost>( x, y, z );

				if ( 0 == j || time < min )
					min = time;
			}

			std::cout <<	std::setw(8) << n << ", " << 
							std::setw(8) << 1000 * min << ", " << 
							std::setw(8) << 1000 * min / n << std::endl;

		}//for

	}//test_host

}// namespace


