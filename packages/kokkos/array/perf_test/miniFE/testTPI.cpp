#include <iostream>
#include <iomanip>
#include <sys/time.h>

#include <Kokkos_DeviceHost.hpp>
#include <Kokkos_DeviceHost_ValueView.hpp>
#include <Kokkos_DeviceHost_MultiVectorView.hpp>
#include <Kokkos_DeviceHost_MDArrayView.hpp>
#include <Kokkos_DeviceHost_ParallelFor.hpp>
#include <Kokkos_DeviceHost_ParallelReduce.hpp>

#include <Kokkos_DeviceTPI.hpp>
#include <Kokkos_DeviceTPI_ValueView.hpp>
#include <Kokkos_DeviceTPI_MultiVectorView.hpp>
#include <Kokkos_DeviceTPI_MDArrayView.hpp>
#include <Kokkos_DeviceTPI_ParallelFor.hpp>
#include <Kokkos_DeviceTPI_ParallelReduce.hpp>

#include <Kokkos_DeviceTPI_macros.hpp>
#include <CRSMesh.hpp>
#include <assemble.hpp>
#include <CRSMatrixGatherFill.hpp>
#include <Dirichlet.hpp>
#include <CG_Solve.hpp>
#include <driver.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

namespace Test{

	void test_TPI(int beg, int end , int runs,  int threads){

		std::cout << "Kokkos TPI - " << threads << " threads: " << std::endl;
		std::cout<<"Size , Mesh , Assemble , Gather , Dirichlet , Iteration , CG Time , Total"<<std::endl;
		Kokkos::DeviceTPI::initialize( threads );


		for(int i = beg ; i < end; i+=2)
 		{
			double times[10], mins[10];
			for(int j = 0; j < runs; j++){

				run_kernel<Kokkos::DeviceTPI>(i,i,i,times);
			
				if(j == 0)
				{
					mins[0] = times[0];
					mins[1] = times[1];
					mins[2] = times[2];
					mins[3] = times[3];
					mins[4] = times[4];
					mins[5] = times[5];
					mins[6] = times[6];
				}
				for(int k = 0 ; k < 7 ; k++)
				{
					if(times[k] < mins[k])
						mins[k] = times[k];
				}
			}
			std::cout	<<i*i*i;
			for(int m = 0 ; m < 7; m++)
			{
				std::cout<<", "<<(times[m] * 1000)<<"ms";
			}
			std::cout<<std::endl;
		}
		Kokkos::DeviceTPI::finalize();

	}//test_TPI

}// namespace


