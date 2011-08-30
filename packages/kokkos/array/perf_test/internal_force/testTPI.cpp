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
#include <internal_force_driver.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

namespace test{

  void test_TPI( int beg, int end, int runs, int threads){

    std::cout << "\"Kokkos TPI: \"" << std::endl;
    std::cout << "\"number of elements\" , \"total time milliseconds\" , \"grind time milliseconds/element\"" << std::endl ;
    Kokkos::DeviceTPI::initialize( threads );

    // Loop over a range of problem sizes

    for(int i = beg; i < end; i++){

      double min = 100000;

      // Number of elements in each direction for the 'box'

      int x = 10 + 5 * i;
      int y = 10 + 5 * i;
      int z = 10 + 5 * i;

      int n = x * y * z;

      // Run multiple times and take the minimum trying to
      // eliminate operating system noise from the timing.

      for(int j = 0; j < runs; j++){
        double time = internal_force_test<double, Kokkos::DeviceTPI>( x, y, z );

        if ( 0 == j || time < min)
          min = time;
      }

      // 'time' and 'min' in seconds, display result in milliseconds,
      //
      std::cout <<     std::setw(8) << n << ", " << 
              std::setw(8) << 1000 * min << ", " << 
              std::setw(8) << 1000 * min / n << std::endl;

    }//for

    Kokkos::DeviceTPI::finalize();

  }//test_host

}// namespace


