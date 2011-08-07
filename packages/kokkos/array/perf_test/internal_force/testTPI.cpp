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
#include <internal_force_test.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

namespace test{

  void test_TPI( int beg, int end, int runs, int threads){

    std::cout << "Kokkos TPI: " << std::endl;
    Kokkos::DeviceTPI::initialize( threads );

    for(int i = beg; i < end; i++){
    
      const size_t n = 1 << i;
      double min = 100000;

      for(int j = 0; j < runs; j++){
        double time = internal_force_test<double, Kokkos::DeviceTPI>( n );

        if ( 0 == j || time < min)
          min = time;
      }

      std::cout <<   std::setw(8) << n << ", " << 
              std::setw(8) << 1000 * min << ", " << 
              std::setw(8) << min / n << std::endl;

    }//for

    Kokkos::DeviceTPI::finalize();

  }//test_host

}// namespace


