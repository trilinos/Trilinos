#include "Tpetra_BlockMap.hpp"

#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION

// #include "Tpetra_ExplicitInstantiationHelpers.hpp"

#include "Tpetra_BlockMap_def.hpp"

#include <Kokkos_SerialNode.hpp>
#if defined(HAVE_KOKKOS_TBB)
#  include <Kokkos_TBBNode.hpp>
#endif
#if defined(HAVE_KOKKOS_THREADPOOL)
#  include <Kokkos_TPINode.hpp>
#endif
#if defined(HAVE_KOKKOS_THRUST)
#  include <Kokkos_ThrustGPUNode.hpp>
#endif

namespace Tpetra {

  TPETRA_BLOCKMAP_INSTANT(int,int,Kokkos::SerialNode)
#if defined(HAVE_KOKKOS_TBB)
  TPETRA_BLOCKMAP_INSTANT(int,int,Kokkos::TBBNode)
#endif
#if defined(HAVE_KOKKOS_THREADPOOL)
    TPETRA_BLOCKMAP_INSTANT(int,int,Kokkos::TPINode)
#endif
#if defined(HAVE_KOKKOS_THRUST)
    TPETRA_BLOCKMAP_INSTANT(int,int,Kokkos::ThrustGPUNode)
#endif

} // namespace Tpetra

#endif // HAVE_TPETRA_EXPLICIT_INSTANTIATION
