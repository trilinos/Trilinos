
#include <iostream>
#include <stdexcept>
#include <limits>
#include <utility>

#include <KokkosArray_Host.hpp>

#include <BoxMeshFixture.hpp>
#include <TestBoxMeshFixture.hpp>
#include <Implicit.hpp>
#include <Nonlinear.hpp>
#include <Explicit.hpp>
#include <SparseLinearSystem.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void test_host_query( comm::Machine machine )
{
  const size_t gang_count = KokkosArray::Host::detect_gang_capacity();

  std::cout << "P" << comm::rank( machine )
            << ": Host gang_count = " << gang_count
            << " , gang_worker_capacity = "
            << KokkosArray::Host::detect_gang_worker_capacity();

  std::cout << std::endl ;
}

//----------------------------------------------------------------------------

void test_host_fixture( comm::Machine machine ,
                        size_t gang_count ,
                        size_t gang_worker_count ,
                        size_t nx , size_t ny , size_t nz )
{
  KokkosArray::Host::initialize( gang_count , gang_worker_count );
  test_box_fixture<KokkosArray::Host>( machine , gang_count , nx , ny , nz );
  KokkosArray::Host::finalize();
}

//----------------------------------------------------------------------------

void test_host_implicit( comm::Machine machine ,
                         size_t gang_count ,
                         size_t gang_worker_count ,
                         size_t elem_count_begin ,
                         size_t elem_count_end ,
                         size_t count_run )
{
  KokkosArray::Host::initialize( gang_count , gang_worker_count );
  HybridFEM::Implicit::driver<double,KokkosArray::Host>( "Host" , machine , gang_count , elem_count_begin , elem_count_end , count_run );
  KokkosArray::Host::finalize();
}

//----------------------------------------------------------------------------

void test_host_explicit( comm::Machine machine ,
                         size_t gang_count ,
                         size_t gang_worker_count ,
                         size_t elem_count_begin ,
                         size_t elem_count_end ,
                         size_t count_run )
{
  KokkosArray::Host::initialize( gang_count , gang_worker_count );
  Explicit::driver<double,KokkosArray::Host>( "Host" , machine , gang_count , elem_count_begin , elem_count_end , count_run );
  KokkosArray::Host::finalize();
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void test_host_nonlinear( comm::Machine machine ,
                          size_t gang_count ,
                          size_t gang_worker_count ,
                          size_t elem_count_begin ,
                          size_t elem_count_end ,
                          size_t count_run )
{
  KokkosArray::Host::initialize( gang_count , gang_worker_count );
  typedef FixtureElementHex8 hex8 ;
  typedef KokkosArray::Host             device ;
  HybridFEM::Nonlinear::driver<double,device,hex8>( "Host" , machine , gang_count , elem_count_begin , elem_count_end , count_run );
  KokkosArray::Host::finalize();
}

void test_host_nonlinear_quadratic( comm::Machine machine ,
                                    size_t gang_count ,
                                    size_t gang_worker_count ,
                                    size_t elem_count_begin ,
                                    size_t elem_count_end ,
                                    size_t count_run )
{
  KokkosArray::Host::initialize( gang_count , gang_worker_count );
  typedef FixtureElementHex27 hex27 ;
  typedef KokkosArray::Host              device ;
  HybridFEM::Nonlinear::driver<double,device,hex27>( "Host" , machine , gang_count , elem_count_begin , elem_count_end , count_run );
  KokkosArray::Host::finalize();
}

//----------------------------------------------------------------------------


