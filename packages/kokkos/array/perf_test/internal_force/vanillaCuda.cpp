#include <cstdlib>
#include <iostream>
#include <iomanip>

double run_cuda_kernel(int n);

/************************************************************/

namespace test{

	void test_Original_Cuda( int beg, int end, int runs){

		std::cout << "Original Cuda: " << std::endl;

		for(int i = beg; i < end; i++){

			int n = 1 << i;
			double time = 100000, min = 100000;

			for(int j = 0; j < runs; j++){

				time = run_cuda_kernel(n);
				if (time < min)
					min = time;

			}

			std::cout << 	std::setw(8) << n << ", " << 
							std::setw(8) << 1000 * time << ", " << 
							std::setw(8) << time / n << std::endl;

		}//for

	}//test_host

}// namespace

