
#include <iostream>
#include <stdexcept>
#include <limits>
#include <utility>

#include <KokkosArray_Host.hpp>

#include <BoxMeshFixture.hpp>
#include <Explicit.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void test_host_explicit( comm::Machine machine ,
                         size_t numa_node_count ,
                         size_t numa_node_thread_count ,
                         size_t elem_count_begin ,
                         size_t elem_count_end ,
                         size_t uq_count_begin ,
                         size_t uq_count_end ,
                         size_t count_run )
{
  KokkosArray::Host::initialize( numa_node_count , numa_node_thread_count );
  Explicit::driver<double,KokkosArray::Host>( "Host" , machine ,
                                              elem_count_begin , elem_count_end ,
                                              uq_count_begin , uq_count_end ,
                                              count_run );
  KokkosArray::Host::finalize();
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

