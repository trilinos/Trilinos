
#include <utility>
#include <iostream>

#include <Kokkos_hwloc.hpp>
#include <Kokkos_Threads.hpp>
#include <Kokkos_Cuda.hpp>
#include <feint_fwd.hpp>

int main()
{
  enum { UseAtomic = false };

#if defined( KOKKOS_HAVE_PTHREAD )
  {
    std::cout << "feint< Threads >" << std::endl ;

    // One team of one thread.
    std::pair<unsigned,unsigned> team_topology( 1 , 1 );

    Kokkos::Threads::initialize( team_topology );
    Kokkos::Example::feint< Kokkos::Threads , UseAtomic >();
    Kokkos::Threads::finalize();
  }
#endif

#if defined( KOKKOS_HAVE_CUDA )
  {
    std::cout << "feint< Cuda >" << std::endl ;

    Kokkos::Cuda::initialize( Kokkos::Cuda::SelectDevice(0) );
    Kokkos::Example::feint< Kokkos::Cuda , UseAtomic >();
    Kokkos::Cuda::finalize();
  }
#endif
}

