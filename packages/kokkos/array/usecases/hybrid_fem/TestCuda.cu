
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

void test_cuda( comm::Machine machine )
{
  // Big question: selection of Cuda defines for parallel processes?

  const size_t comm_rank = comm::rank( machine );
  const size_t comm_size = comm::size( machine );
  const size_t dev_size  = Kokkos::Cuda::detect_device_count();

  std::cout << "P" << comm_rank
            << ": test_cuda { comm_size = " << comm_size
            << " dev_size = " << dev_size 
            << " }" << std::endl ;

  if ( comm_size <= dev_size ) {
    Kokkos::Cuda::initialize( Kokkos::Cuda::SelectDevice( comm_rank ) );
    test_box_fixture<Kokkos::Cuda>( machine , 50 , 100 , 200 );
    Kokkos::Cuda::finalize();
  }

  HybridFEM::Implicit::driver<double,Kokkos::Cuda>( "Cuda" , machine , 0 , 0 , 0 );

}

