
#include <sstream>
#include <KokkosArray_Host.hpp>
#include <KokkosArray_hwloc.hpp>

#include <TestSpmv.hpp>
#include <TestCG.hpp>

namespace Test {

int test_host()
{
 const std::pair<unsigned,unsigned> core_top =
    KokkosArray::hwloc::get_core_topology();

  const unsigned core_size =
    KokkosArray::hwloc::get_core_capacity();

  const unsigned gang_count        = core_top.first ;
  const unsigned gang_worker_count = core_top.second * core_size ;

  std::ostringstream msg ;

  msg << "Host[ " << gang_count << " x " << gang_worker_count << " ]" ;

  KokkosArray::Host::initialize( gang_count , gang_worker_count );

  test_spmv_driver<KokkosArray::Host>(msg.str().c_str());

  test_cgsolve_driver<KokkosArray::Host>(msg.str().c_str());

  KokkosArray::Host::finalize();

  return 0 ;
}

}

