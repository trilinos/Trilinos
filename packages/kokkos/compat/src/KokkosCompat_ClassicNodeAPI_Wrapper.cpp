#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>
#include <Kokkos_HostSpace.hpp>

//static size_t KokkosHostWrapperNodeDUMMYINT = Kokkos::NEVEREVERUSEMEIWILLFINDYOU::host_space_singleton_wrapper().size();

#ifdef KOKKOS_HAVE_CUDA
  //static size_t KokkosCudaWrapperNodeDUMMYINT = Kokkos::Impl::cuda_internal_maximum_warp_count();
  template<> int Kokkos::Compat::KokkosCudaWrapperNode::count = 0;// KokkosCudaWrapperNodeDUMMYINT;
#endif

#ifdef KOKKOS_HAVE_OPENMP
  template<> int Kokkos::Compat::KokkosOpenMPWrapperNode::count = 0;
#endif

#ifdef KOKKOS_HAVE_PTHREAD
  template<> int Kokkos::Compat::KokkosThreadsWrapperNode::count = 0;
#endif


namespace Kokkos {
  namespace Compat {

#ifdef KOKKOS_HAVE_PTHREAD
    template<>
    KokkosDeviceWrapperNode<Kokkos::Threads>::~KokkosDeviceWrapperNode<Kokkos::Threads>() {
      count--;
      if((count==0) && Threads::is_initialized())
        Threads::finalize();
    }
    template<>
    void KokkosDeviceWrapperNode<Kokkos::Threads>::init(int NumTeams, int NumThreads, int Device) {
      if(!Kokkos::Threads::is_initialized())
        Kokkos::Threads::initialize(NumTeams*NumThreads);
    }
#endif

#ifdef KOKKOS_HAVE_OPENMP
    template<>
    KokkosDeviceWrapperNode<Kokkos::OpenMP>::~KokkosDeviceWrapperNode<Kokkos::OpenMP>() {
      count--;
      if((count==0) && OpenMP::is_initialized())
        OpenMP::finalize();
    }
    template<>
    void KokkosDeviceWrapperNode<Kokkos::OpenMP>::init(int NumTeams, int NumThreads, int Device) {
      if(!Kokkos::OpenMP::is_initialized())
        Kokkos::OpenMP::initialize(NumTeams*NumThreads);
    }
#endif

#ifdef KOKKOS_HAVE_CUDA
    template<>
    KokkosDeviceWrapperNode<Kokkos::Cuda>::~KokkosDeviceWrapperNode<Kokkos::Cuda>() {
      count--;
      if(count==0) {
        if(Cuda::host_mirror_device_type::is_initialized())
          Cuda::host_mirror_device_type::finalize();
        if(Cuda::is_initialized())
          Cuda::finalize();
      }
    }
    template<>
    void KokkosDeviceWrapperNode<Kokkos::Cuda>::init(int NumTeams, int NumThreads, int Device) {
      if(!Kokkos::Cuda::host_mirror_device_type::is_initialized())
        Kokkos::Cuda::host_mirror_device_type::initialize(NumTeams*NumThreads);
      Kokkos::Cuda::SelectDevice select_device(Device);
      Kokkos::Cuda::initialize(select_device);
    }
#endif
  }
}

//  Make sure HostSpace is always intialized before Devices.
//  Some TPetra Unit tests and possibly some codes have static Nodes
//  Create a static node here in order to make sure that HostSpace is initialized first

/*
#ifdef KOKKOS_HAVE_CUDA
  static Kokkos::Compat::KokkosCudaWrapperNode KOKKOSCUDANODE_;
#endif

#ifdef KOKKOS_HAVE_OPENMP
  static Kokkos::Compat::KokkosOpenMPWrapperNode KOKKOSOPENMPNODE_;
#endif

#ifdef KOKKOS_HAVE_PTHREAD
  static Kokkos::Compat::KokkosThreadsWrapperNode KOKKOSTHREADSNODE_;
#endif
*/
