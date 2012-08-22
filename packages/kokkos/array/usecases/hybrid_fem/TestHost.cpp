
/* #define KOKKOSARRAY_BOUNDS_CHECK 1 */

#include <iostream>
#include <stdexcept>
#include <limits>
#include <utility>

#include <BoxMeshFixture.hpp>
#include <TestBoxMeshFixture.hpp>
#include <Implicit.hpp>
#include <NonLinear.hpp>
#include <Explicit.hpp>

#include <KokkosArray_Host.hpp>

#include <KokkosArray_Host_macros.hpp>
#include <ParallelDataMap_macros.hpp>
#include <TestBoxMeshFixture_macros.hpp>
#include <SparseLinearSystem_macros.hpp>
#include <SparseLinearSystemFill_macros.hpp>

#include <Implicit_macros.hpp>
#include <NonLinear_macros.hpp>
#include <NonlinearElement_macros.hpp>
#include <Explicit_macros.hpp>
#include <KokkosArray_Clear_macros.hpp>

#include <SparseLinearSystem_Host.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void test_host_query( comm::Machine machine )
{
  const size_t node_count = KokkosArray::Host::detect_node_count();

  std::cout << "P" << comm::rank( machine )
            << ": Host node_count = " << node_count
            << " , node_core_count = "
            << KokkosArray::Host::detect_node_core_count();

  std::cout << std::endl ;
}

//----------------------------------------------------------------------------

void test_host_fixture( comm::Machine machine ,
                        size_t numa_node_count ,
                        size_t numa_node_thread_count ,
                        size_t nx , size_t ny , size_t nz )
{
  KokkosArray::Host::initialize( numa_node_count , numa_node_thread_count );
  test_box_fixture<KokkosArray::Host>( machine , nx , ny , nz );
  KokkosArray::Host::finalize();
}

//----------------------------------------------------------------------------

void test_host_implicit( comm::Machine machine , 
                         size_t numa_node_count ,
                         size_t numa_node_thread_count ,
                         size_t elem_count_begin ,
                         size_t elem_count_end ,
                         size_t count_run )
{
  KokkosArray::Host::initialize( numa_node_count , numa_node_thread_count );
  HybridFEM::Implicit::driver<double,KokkosArray::Host>( "Host" , machine , elem_count_begin , elem_count_end , count_run );
  KokkosArray::Host::finalize();
}

//----------------------------------------------------------------------------

void test_host_explicit( comm::Machine machine , 
                         size_t numa_node_count ,
                         size_t numa_node_thread_count ,
                         size_t elem_count_begin ,
                         size_t elem_count_end ,
                         size_t count_run )
{
  KokkosArray::Host::initialize( numa_node_count , numa_node_thread_count );
  Explicit::driver<double,KokkosArray::Host>( "Host" , machine , elem_count_begin , elem_count_end , count_run );
  KokkosArray::Host::finalize();
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void test_host_nonlinear( comm::Machine machine , 
                          size_t numa_node_count ,
                          size_t numa_node_thread_count ,
                          size_t elem_count_begin ,
                          size_t elem_count_end ,
                          size_t count_run )
{
  KokkosArray::Host::initialize( numa_node_count , numa_node_thread_count );
  typedef FixtureElementHex8 hex8 ;
  typedef KokkosArray::Host             device ;
  HybridFEM::NonLinear::driver<double,device,hex8>( "Host" , machine , elem_count_begin , elem_count_end , count_run );
  KokkosArray::Host::finalize();
}

void test_host_nonlinear_quadratic( comm::Machine machine , 
                                    size_t numa_node_count ,
                                    size_t numa_node_thread_count ,
                                    size_t elem_count_begin ,
                                    size_t elem_count_end ,
                                    size_t count_run )
{
  KokkosArray::Host::initialize( numa_node_count , numa_node_thread_count );
  typedef FixtureElementHex27 hex27 ;
  typedef KokkosArray::Host              device ;
  HybridFEM::NonLinear::driver<double,device,hex27>( "Host" , machine , elem_count_begin , elem_count_end , count_run );
  KokkosArray::Host::finalize();
}

//----------------------------------------------------------------------------


