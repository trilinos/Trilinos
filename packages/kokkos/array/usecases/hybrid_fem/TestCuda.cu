
#include <KokkosArray_Cuda.hpp>
#include <KokkosArray_Host.hpp>

#include <TestBoxMeshFixture.hpp>
#include <Implicit.hpp>
#include <Nonlinear.hpp>
#include <Explicit.hpp>
#include <SparseLinearSystem.hpp>

//----------------------------------------------------------------------------

void test_cuda_query( comm::Machine machine )
{
  const size_t comm_rank = comm::rank( machine );
  std::cout << "P" << comm_rank
            << ": Cuda device_count = "
            << KokkosArray::Cuda::detect_device_count()
            << std::endl ;
}

//----------------------------------------------------------------------------

void test_cuda_fixture( comm::Machine machine ,
                        size_t nx , size_t ny , size_t nz )
{
  const size_t comm_rank = comm::rank( machine );
  const size_t comm_size = comm::size( machine );
  const size_t dev_count = KokkosArray::Cuda::detect_device_count();
  const size_t dev_rank =
    dev_count && dev_count <= comm_size ? comm_rank % dev_count : 0 ;
  const size_t gang_count = 0 ;

  KokkosArray::Cuda::SelectDevice select_device( dev_rank );
  KokkosArray::Cuda::initialize( select_device );
  test_box_fixture<KokkosArray::Cuda>( machine , gang_count , nx , ny , nz );
  KokkosArray::Cuda::finalize();
}

//----------------------------------------------------------------------------

void test_cuda_implicit( comm::Machine machine , 
                         size_t elem_count_begin ,
                         size_t elem_count_end ,
                         size_t count_run )
{
  const size_t comm_rank = comm::rank( machine );
  const size_t comm_size = comm::size( machine );
  const size_t dev_count = KokkosArray::Cuda::detect_device_count();
  const size_t dev_rank =
    dev_count && dev_count <= comm_size ? comm_rank % dev_count : 0 ;
  const size_t gang_count = 0 ;

  KokkosArray::Cuda::SelectDevice select_device( dev_rank );
  KokkosArray::Cuda::initialize( select_device );
  HybridFEM::Implicit::driver<double,KokkosArray::Cuda>( "Cuda" , machine , gang_count , elem_count_begin , elem_count_end , count_run );
  KokkosArray::Cuda::finalize();
}

//----------------------------------------------------------------------------

void test_cuda_explicit( comm::Machine machine , 
                         size_t elem_count_begin ,
                         size_t elem_count_end ,
                         size_t count_run )
{
  const size_t comm_rank = comm::rank( machine );
  const size_t comm_size = comm::size( machine );
  const size_t dev_count = KokkosArray::Cuda::detect_device_count();
  const size_t dev_rank =
    dev_count && dev_count <= comm_size ? comm_rank % dev_count : 0 ;
  const size_t gang_count = 0 ;

  KokkosArray::Cuda::SelectDevice select_device( dev_rank );
  KokkosArray::Cuda::initialize( select_device );
  Explicit::driver<double,KokkosArray::Cuda>( "Cuda" , machine , gang_count , elem_count_begin , elem_count_end , count_run );
  KokkosArray::Cuda::finalize();
}

//----------------------------------------------------------------------------

void test_cuda_nonlinear( comm::Machine machine , 
                          size_t elem_count_begin ,
                          size_t elem_count_end ,
                          size_t count_run )
{
  const size_t comm_rank = comm::rank( machine );
  const size_t comm_size = comm::size( machine );
  const size_t dev_count = KokkosArray::Cuda::detect_device_count();
  const size_t dev_rank =
    dev_count && dev_count <= comm_size ? comm_rank % dev_count : 0 ;
  const size_t gang_count = 0 ;

  KokkosArray::Cuda::SelectDevice select_device( dev_rank );
  KokkosArray::Cuda::initialize( select_device );

  typedef KokkosArray::Cuda device ;
  typedef FixtureElementHex8 hex8 ;
  HybridFEM::Nonlinear::driver<double,device,hex8>( "Cuda" , machine , gang_count , elem_count_begin , elem_count_end , count_run );
  KokkosArray::Cuda::finalize();
}

void test_cuda_nonlinear_quadratic( comm::Machine machine , 
                                    size_t elem_count_begin ,
                                    size_t elem_count_end ,
                                    size_t count_run )
{
  const size_t comm_rank = comm::rank( machine );
  const size_t comm_size = comm::size( machine );
  const size_t dev_count = KokkosArray::Cuda::detect_device_count();
  const size_t dev_rank =
    dev_count && dev_count <= comm_size ? comm_rank % dev_count : 0 ;
  const size_t gang_count = 0 ;

  KokkosArray::Cuda::SelectDevice select_device( dev_rank );
  KokkosArray::Cuda::initialize( select_device );

  typedef KokkosArray::Cuda device ;
  typedef FixtureElementHex27 hex27 ;
  HybridFEM::Nonlinear::driver<double,device,hex27>( "Cuda" , machine , gang_count , elem_count_begin , elem_count_end , count_run );
  KokkosArray::Cuda::finalize();
}

//----------------------------------------------------------------------------

