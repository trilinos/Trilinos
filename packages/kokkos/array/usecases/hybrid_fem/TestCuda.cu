
#include <TestBoxMeshFixture.hpp>
#include <Implicit.hpp>
#include <NonLinear.hpp>

#include <Kokkos_Cuda.hpp>
#include <Kokkos_Host.hpp>

#include <Kokkos_Cuda_macros.hpp>
#include <ParallelDataMap_macros.hpp>
#include <TestBoxMeshFixture_macros.hpp>
#include <SparseLinearSystem_macros.hpp>
#include <SparseLinearSystemFill_macros.hpp>
#include <Implicit_macros.hpp>
#include <NonLinear_macros.hpp>
#include <Kokkos_Clear_macros.hpp>

#include <SparseLinearSystem_Cuda.hpp>

//----------------------------------------------------------------------------

void test_cuda_query( comm::Machine machine )
{
  const size_t comm_rank = comm::rank( machine );
  const size_t dev_count = Kokkos::Cuda::detect_device_count();
  std::cout << "P" << comm_rank
            << ": Cuda device_count = "
            << Kokkos::Cuda::detect_device_count()
            << std::endl ;
}

//----------------------------------------------------------------------------

void test_cuda_fixture( comm::Machine machine ,
                        size_t nx , size_t ny , size_t nz )
{
  const size_t comm_rank = comm::rank( machine );
  const size_t comm_size = comm::size( machine );
  const size_t dev_count = Kokkos::Cuda::detect_device_count();
  const size_t dev_rank =
    dev_count && dev_count <= comm_size ? comm_rank % dev_count : 0 ;

  Kokkos::Cuda::SelectDevice select_device( dev_rank );
  Kokkos::Cuda::initialize( select_device );
  test_box_fixture<Kokkos::Cuda>( machine , nx , ny , nz );
  Kokkos::Cuda::finalize();
}

//----------------------------------------------------------------------------

void test_cuda_implicit( comm::Machine machine , 
                         size_t node_count_begin ,
                         size_t node_count_end ,
                         size_t count_run )
{
  const size_t comm_rank = comm::rank( machine );
  const size_t comm_size = comm::size( machine );
  const size_t dev_count = Kokkos::Cuda::detect_device_count();
  const size_t dev_rank =
    dev_count && dev_count <= comm_size ? comm_rank % dev_count : 0 ;

  Kokkos::Cuda::SelectDevice select_device( dev_rank );
  Kokkos::Cuda::initialize( select_device );
  HybridFEM::Implicit::driver<double,Kokkos::Cuda>( "Cuda" , machine , node_count_begin , node_count_end , count_run );
  Kokkos::Cuda::finalize();
}

//----------------------------------------------------------------------------

void test_cuda_nonlinear( comm::Machine machine , 
                          size_t node_count_begin ,
                          size_t node_count_end ,
                          size_t count_run )
{
  const size_t comm_rank = comm::rank( machine );
  const size_t comm_size = comm::size( machine );
  const size_t dev_count = Kokkos::Cuda::detect_device_count();
  const size_t dev_rank =
    dev_count && dev_count <= comm_size ? comm_rank % dev_count : 0 ;

  Kokkos::Cuda::SelectDevice select_device( dev_rank );
  Kokkos::Cuda::initialize( select_device );
  HybridFEM::NonLinear::driver<double,Kokkos::Cuda>( "Cuda" , machine , node_count_begin , node_count_end , count_run );
  Kokkos::Cuda::finalize();
}

//----------------------------------------------------------------------------

