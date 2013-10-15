#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

namespace Kokkos {
  namespace Compat {

#ifdef KOKKOS_HAVE_PTHREAD
    template<>
    KokkosDeviceWrapperNode<Kokkos::Threads>::~KokkosDeviceWrapperNode<Kokkos::Threads>() {
      Threads::finalize();
    }
    template<>
    void KokkosDeviceWrapperNode<Kokkos::Threads>::init(int NumTeams, int NumThreads, int Device) {
      Kokkos::Threads::initialize(NumTeams,NumThreads);
    }
#endif

#ifdef KOKKOS_HAVE_OPENMP
    template<>
    KokkosDeviceWrapperNode<Kokkos::OpenMP>::~KokkosDeviceWrapperNode<Kokkos::OpenMP>() {
      OpenMP::finalize();
    }
    template<>
    void KokkosDeviceWrapperNode<Kokkos::OpenMP>::init(int NumTeams, int NumThreads, int Device) {
      Kokkos::OpenMP::initialize(NumTeams,NumThreads);
    }
#endif

#ifdef KOKKOS_HAVE_CUDA
    template<>
    KokkosDeviceWrapperNode<Kokkos::Cuda>::~KokkosDeviceWrapperNode<Kokkos::Cuda>() {
      Cuda::host_mirror_device_type::finalize();
      Cuda::finalize();
    }
    template<>
    void KokkosDeviceWrapperNode<Kokkos::Cuda>::init(int NumTeams, int NumThreads, int Device) {
      Kokkos::Cuda::SelectDevice select_device(Device);
      Kokkos::Cuda::initialize(select_device);
      Kokkos::Cuda::host_mirror_device_type::initialize(NumTeams,NumThreads);
      //Kokkos::OpenMP::initialize(NumTeams,NumThreads);
      std::cout<<"Init Devices "<<NumTeams<<" "<<NumThreads<<"\n";
    }
#endif

  }
}
