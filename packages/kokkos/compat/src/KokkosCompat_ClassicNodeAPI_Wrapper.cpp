#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>
#include <Kokkos_HostSpace.hpp>

//static size_t KokkosHostWrapperNodeDUMMYINT = Kokkos::NEVEREVERUSEMEIWILLFINDYOU::host_space_singleton_wrapper().size();

namespace Kokkos {
  namespace Compat {

    // mfh 01 Jan 2014: These definitions of the class variable count
    // need to be inside the namespace.  Declaring them as "template<>
    // int Kokkos::Compat::KokkosCudaWrapperNode::count = 0" in the
    // global namespace is a C++11 extension and results in compiler
    // warnings with Clang 3.2 on MacOS X.
#ifdef KOKKOS_HAVE_CUDA
    template<> int KokkosCudaWrapperNode::count = 0;
#endif
#ifdef KOKKOS_HAVE_OPENMP
    template<> int KokkosOpenMPWrapperNode::count = 0;
#endif
#ifdef KOKKOS_HAVE_PTHREAD
    template<> int KokkosThreadsWrapperNode::count = 0;
#endif

#ifdef KOKKOS_HAVE_PTHREAD
    template<>
    KokkosDeviceWrapperNode<Kokkos::Threads>::~KokkosDeviceWrapperNode<Kokkos::Threads>() {
      count--;
      if((count==0) && Threads::is_initialized()) {
#ifdef KOKKOS_HAVE_CUDA
        if(!Impl::is_same<Kokkos::Threads,Cuda::host_mirror_device_type>::value ||
            KokkosDeviceWrapperNode<Kokkos::Cuda>::count==0)
#endif
        Threads::finalize();
      }
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
      if((count==0) && OpenMP::is_initialized()) {
#ifdef KOKKOS_HAVE_CUDA
        if(!Impl::is_same<Kokkos::OpenMP,Cuda::host_mirror_device_type>::value ||
            KokkosDeviceWrapperNode<Kokkos::Cuda>::count==0)
#endif
        OpenMP::finalize();
      }
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
        if(Cuda::host_mirror_device_type::is_initialized()) {
          // make sure that no Actual DeviceWrapper node of the mirror_device_type is in use
          if(KokkosDeviceWrapperNode<Cuda::host_mirror_device_type>::count==0)
            Cuda::host_mirror_device_type::finalize();
        }
        if(Cuda::is_initialized())
          Cuda::finalize();
      }
    }
    template<>
    void KokkosDeviceWrapperNode<Kokkos::Cuda>::init(int NumTeams, int NumThreads, int Device) {
      if(!Kokkos::Cuda::host_mirror_device_type::is_initialized())
        Kokkos::Cuda::host_mirror_device_type::initialize(NumTeams*NumThreads);
      Kokkos::Cuda::SelectDevice select_device(Device);
      if(!Kokkos::Cuda::is_initialized())
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
