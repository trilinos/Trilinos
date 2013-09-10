
#include <utility>
#include <iostream>

#include <KokkosCore_config.h>
#include <Kokkos_hwloc.hpp>
#include <Kokkos_Threads.hpp>

#if defined( KOKKOS_HAVE_CUDA )
#include <Kokkos_Cuda.hpp>
#endif

#include <feint_fwd.hpp>

int main()
{
#if defined( KOKKOS_HAVE_PTHREAD )
  {
    std::pair<unsigned,unsigned> use_cores =
      Kokkos::hwloc::get_core_topology();

    // Use 4 cores per NUMA region, unless fewer available
    use_cores.second = std::min( 4u , use_cores.second );

    // Use 2 cores per team and 1 thread/core:
    std::pair<unsigned,unsigned>
      team_topology( use_cores.first * use_cores.second / 2 , 2 );

    Kokkos::Threads::initialize( team_topology );

    std::cout << "feint< Threads , NotUsingAtomic >" << std::endl ;
    Kokkos::Example::feint< Kokkos::Threads , false >();

    std::cout << "feint< Threads , Usingtomic >" << std::endl ;
    Kokkos::Example::feint< Kokkos::Threads , true  >();

    Kokkos::Threads::finalize();
  }
#endif

#if defined( KOKKOS_HAVE_CUDA )
  {
    const unsigned device_count = Kokkos::Cuda::detect_device_count();

    // Use the last device:

    Kokkos::Cuda::initialize( Kokkos::Cuda::SelectDevice(device_count-1) );

    std::cout << "feint< Cuda , NotUsingAtomic >" << std::endl ;
    Kokkos::Example::feint< Kokkos::Cuda , false >();

    std::cout << "feint< Cuda , UsingAtomic >" << std::endl ;
    Kokkos::Example::feint< Kokkos::Cuda , true  >();

    Kokkos::Cuda::finalize();
  }
#endif
}

