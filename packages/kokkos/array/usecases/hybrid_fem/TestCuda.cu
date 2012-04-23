
#include <TestBoxMeshFixture.hpp>
#include <Implicit.hpp>

#include <Kokkos_Cuda.hpp>
#include <Kokkos_Host.hpp>

#include <Kokkos_Cuda_macros.hpp>
#include <ParallelDataMap_macros.hpp>
#include <TestBoxMeshFixture_macros.hpp>
#include <SparseLinearSystem_macros.hpp>
#include <SparseLinearSystemFill_macros.hpp>
#include <Implicit_macros.hpp>
#include <Kokkos_Clear_macros.hpp>

#include <SparseLinearSystem_Cuda.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void test_cuda( comm::Machine machine , std::istream & input )
{
  const size_t comm_rank = comm::rank( machine );
  const size_t comm_size = comm::size( machine );
  const size_t dev_size  = Kokkos::Cuda::detect_device_count();

  Kokkos::Cuda::SelectDevice select_device(0);

  for(;;) {
    std::string which ; input >> which ;

    if ( std::string("device=comm") == which ) {
      if ( dev_size < comm_size ) {
        std::ostringstream msg ;
        msg << "cuda-device-count " << dev_size
            << " < mpi-comm-size " << comm_size ;
        throw std::runtime_error(msg.str());
      }
      select_device = Kokkos::Cuda::SelectDevice( comm_rank );
    }
    else if ( std::string("device=zero") == which ) {
      select_device = Kokkos::Cuda::SelectDevice(0);
    }
    else if ( std::string("fixture") == which ) {
      size_t fixture[3] = { 10 , 20 , 50 };
      input >> fixture[0] >> fixture[1] >> fixture[2] ;

      Kokkos::Cuda::initialize( select_device );
      test_box_fixture<Kokkos::Cuda>( machine , fixture[0] , fixture[1] , fixture[2] );
      break ;
    }
    else if ( std::string("implicit") == which ) {
      size_t count_begin = 20 , count_end = 200 , count_run = 1 ;
      input >> count_begin >> count_end >> count_run ;

      Kokkos::Cuda::initialize( select_device );
      HybridFEM::Implicit::driver<double,Kokkos::Cuda>( "Cuda" , machine , count_begin , count_end , count_run );
      break ;
    }
    else {
      if ( 0 == comm::rank( machine ) ) {
        std::cout
          << "cuda device=(comm|zero) fixture NX NY NZ" << std::endl
          << "cuda device=(comm|zero) implicit Nbegin Nend Nrun" << std::endl
          << "  # if \"device=comm\" : select device as comm-rank" << std::endl
          << "  # if \"device=zero\" : select device as zero" << std::endl
          << "  # cuda-device-count = " << dev_size << std::endl ;
      }
      break ;
    }
  }

  Kokkos::Cuda::finalize();
}

