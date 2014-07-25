
// Must be included first on Intel-Phi systems due to
// redefinition of SEEK_SET in <mpi.h>.

#include <ParallelComm.hpp>

#include <iostream>
#include <stdexcept>
#include <limits>
#include <utility>

//----------------------------------------------------------------------------

#include <BoxMeshFixture.hpp>
#include <TestBoxMeshFixture.hpp>
#include <Implicit.hpp>
#include <Nonlinear.hpp>
#include <Explicit.hpp>
#include <SparseLinearSystem.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void test_host_fixture( comm::Machine machine ,
                        size_t gang_count ,
                        size_t gang_worker_count ,
                        size_t nx , size_t ny , size_t nz )
{
  HostExecSpace::initialize( gang_count * gang_worker_count );
  test_box_fixture<HostExecSpace>( machine , gang_count , nx , ny , nz );
  HostExecSpace::finalize();
}

//----------------------------------------------------------------------------

void test_host_implicit( comm::Machine machine ,
                         size_t gang_count ,
                         size_t gang_worker_count ,
                         size_t elem_count_begin ,
                         size_t elem_count_end ,
                         size_t count_run )
{
  HostExecSpace::initialize( gang_count * gang_worker_count );
  HybridFEM::Implicit::driver<double,HostExecSpace>( "Threads" , machine , gang_count , elem_count_begin , elem_count_end , count_run );
  HostExecSpace::finalize();
}

//----------------------------------------------------------------------------

void test_host_explicit( comm::Machine machine ,
                         size_t gang_count ,
                         size_t gang_worker_count ,
                         size_t elem_count_begin ,
                         size_t elem_count_end ,
                         size_t count_run )
{
  HostExecSpace::initialize( gang_count * gang_worker_count );
  Explicit::driver<double,HostExecSpace>( "Threads" , machine , gang_count , elem_count_begin , elem_count_end , count_run );
  HostExecSpace::finalize();
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
  HostExecSpace::initialize( gang_count * gang_worker_count );
  typedef FixtureElementHex8 hex8 ;
  typedef HostExecSpace             device ;
  HybridFEM::Nonlinear::driver<double,device,hex8>( "Threads" , machine , gang_count , elem_count_begin , elem_count_end , count_run );
  HostExecSpace::finalize();
}

void test_host_nonlinear_quadratic( comm::Machine machine ,
                                    size_t gang_count ,
                                    size_t gang_worker_count ,
                                    size_t elem_count_begin ,
                                    size_t elem_count_end ,
                                    size_t count_run )
{
  HostExecSpace::initialize( gang_count * gang_worker_count );
  typedef FixtureElementHex27 hex27 ;
  typedef HostExecSpace              device ;
  HybridFEM::Nonlinear::driver<double,device,hex27>( "Threads" , machine , gang_count , elem_count_begin , elem_count_end , count_run );
  HostExecSpace::finalize();
}

//----------------------------------------------------------------------------


