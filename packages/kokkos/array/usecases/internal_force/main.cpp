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

int main()
{
  std::cout << "Kokkos Host SM miniapp: " << std::endl;
  const int x = 100;
  const int y = 10;
  const int z = 10;
  const int num_elements = x * y * z;

  double time = internal_force_test<double, Kokkos::DeviceHost>( x, y, z );
  std::cout <<	std::setw(8) << num_elements << ", " <<
    std::setw(8) << 1000 * time << ", " <<
    std::setw(8) << 1000 * time / num_elements << std::endl;
  return 0;
}

#if 0
#include <iostream>
#include <cstdlib>

namespace test{

  void test_Host(int b, int e, int r);
  void test_TPI(int b, int e, int r, int t);
  void test_Cuda(int b, int e, int r);

}

int main( int argc , char ** argv ){

  int threads, runs;
  int beg = 22;
  int end = 23;

  if (argc > 1)
    threads = atoi(argv[1]);
  else
    threads = 4;

  if (argc > 2)
    runs = atoi(argv[2]);
  else
    runs = 1;

/************************************************************/
/*  argument 1 specifies the number of threads used by     	*/
/*  TPI and OpenMP; agrument 2 sets the number of runs    	*/
/*  to be averaged for the output timings          			*/
/*                             				 				*/
/*  Ex: to run with 8 threads and average 10 runs(LINUX):  	*/
/*  ./fe_test.exe 8 10                    					*/
/************************************************************/

  std::cout << std::endl << "\tStarting..." << std::endl;

  test::test_Host(beg, end, runs);
//  test::test_TPI (beg, end, runs, threads);
//  test::test_Cuda(beg, end, runs);

  return 0;

}
#endif

