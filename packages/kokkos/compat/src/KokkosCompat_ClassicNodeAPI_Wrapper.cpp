#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

namespace Kokkos {
  namespace Compat {

#ifdef KOKKOS_HAVE_PTHREAD
    template<>
    void KokkosDeviceWrapperNode<Kokkos::Threads>::init(int NumTeams, int NumThreads, int Device) {
      Kokkos::Threads::initialize(NumTeams,NumThreads);
    }
#endif

#ifdef KOKKOS_HAVE_OPENMP
    template<>
    void KokkosDeviceWrapperNode<Kokkos::OpenMP>::init(int NumTeams, int NumThreads, int Device) {
      Kokkos::OpenMP::initialize(NumTeams,NumThreads);
    }
#endif

#ifdef KOKKOS_HAVE_CUDA
    template<>
    void KokkosDeviceWrapperNode<Kokkos::Cuda>::init(int NumTeams, int NumThreads, int Device) {
      Kokkos::Cuda::SelectDevice select_device(Device);
      Kokkos::Cuda::initialize(select_device);
    }
#endif

  }
}
