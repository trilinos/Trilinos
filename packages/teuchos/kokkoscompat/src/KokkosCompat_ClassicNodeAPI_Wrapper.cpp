#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>
#include <Kokkos_Core.hpp>

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
#ifdef KOKKOS_HAVE_SERIAL
    template<> int KokkosSerialWrapperNode::count = 0;
#endif // KOKKOS_HAVE_SERIAL

#ifdef KOKKOS_HAVE_PTHREAD
    template<>
    KokkosDeviceWrapperNode<Kokkos::Threads>::
    ~KokkosDeviceWrapperNode<Kokkos::Threads> ()
    {
      count--;
      if (count == 0 && Threads::is_initialized ()) {
#ifdef KOKKOS_HAVE_CUDA
        if (! Impl::is_same<Kokkos::Threads,HostSpace::execution_space>::value ||
            KokkosDeviceWrapperNode<Kokkos::Cuda>::count == 0)
#endif
          Threads::finalize ();
      }
    }

    template<>
    void KokkosDeviceWrapperNode<Kokkos::Threads>::
    init (int NumThreads, int NumNUMA, int NumCoresPerNUMA, int Device) {
      if (! Kokkos::Threads::is_initialized ()) {
        if(NumNUMA>0 && NumCoresPerNUMA>0)
          Kokkos::Threads::initialize ( NumThreads, NumNUMA, NumCoresPerNUMA );
        else if (NumNUMA > 0)
          Kokkos::Threads::initialize ( NumThreads, NumNUMA );
        else
          Kokkos::Threads::initialize ( NumThreads );
      }
    }

    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::Threads>::name () {
      return "Threads/Wrapper";
    }
#endif

#ifdef KOKKOS_HAVE_OPENMP
    template<>
    KokkosDeviceWrapperNode<Kokkos::OpenMP>::~KokkosDeviceWrapperNode<Kokkos::OpenMP>() {
      count--;
      if (count == 0 && OpenMP::is_initialized ()) {
#ifdef KOKKOS_HAVE_CUDA
        if (! Impl::is_same<Kokkos::OpenMP, HostSpace::execution_space>::value ||
            KokkosDeviceWrapperNode<Kokkos::Cuda>::count == 0)
#endif
        OpenMP::finalize ();
      }
    }

    template<>
    void KokkosDeviceWrapperNode<Kokkos::OpenMP>::
    init (int NumThreads, int NumNUMA, int NumCoresPerNUMA, int Device) {
      if (! Kokkos::OpenMP::is_initialized ()) {
        if(NumNUMA>0 && NumCoresPerNUMA>0)
          Kokkos::OpenMP::initialize ( NumThreads, NumNUMA, NumCoresPerNUMA );
        else if (NumNUMA > 0)
          Kokkos::OpenMP::initialize ( NumThreads, NumNUMA );
        else
          Kokkos::OpenMP::initialize ( NumThreads );
      }
    }

    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::OpenMP>::name () {
      return "OpenMP/Wrapper";
    }
#endif

#ifdef KOKKOS_HAVE_SERIAL
    template<>
    KokkosDeviceWrapperNode<Kokkos::Serial>::~KokkosDeviceWrapperNode<Kokkos::Serial>() {
      count--;
      if (count == 0 && Serial::is_initialized ()) {
        // FIXME (mfh 02 Nov 2014) This doesn't look right to me
        // somehow.  What if CUDA is NOT enabled?  Wouldn't that
        // disable the first branch of the OR, which could trigger
        // even if CUDA is not enabled?
#ifdef KOKKOS_HAVE_CUDA
        if (! Impl::is_same<Kokkos::Serial, HostSpace::execution_space>::value ||
            KokkosDeviceWrapperNode<Kokkos::Cuda>::count == 0)
#endif // KOKKOS_HAVE_CUDA
          Serial::finalize ();
      }
    }

    template<>
    void KokkosDeviceWrapperNode<Kokkos::Serial>::
    init (int NumThreads, int NumNUMA, int NumCoresPerNUMA, int Device) {
      if (! Kokkos::Serial::is_initialized ()) {
          Kokkos::Serial::initialize ();
      }
    }

    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::Serial>::name () {
      return "Serial/Wrapper";
    }
#endif // KOKKOS_HAVE_SERIAL

#ifdef KOKKOS_HAVE_CUDA
    template<>
    KokkosDeviceWrapperNode<Kokkos::Cuda>::~KokkosDeviceWrapperNode<Kokkos::Cuda>() {
      count--;
      if(count==0) {
        if(HostSpace::execution_space::is_initialized()) {
          // make sure that no Actual DeviceWrapper node of the mirror_execution_space is in use
          if(KokkosDeviceWrapperNode<HostSpace::execution_space>::count==0) {
            HostSpace::execution_space::finalize();
          }
        }
        if(Cuda::is_initialized())
          Cuda::finalize();
      }
    }
    template<>
    void KokkosDeviceWrapperNode<Kokkos::Cuda>::
    init(int NumThreads, int NumNUMA, int NumCoresPerNUMA, int Device) {

      // Setting (currently) necessary environment variables for NVIDIA UVM
      #ifdef KOKKOS_USE_CUDA_UVM
        putenv("CUDA_LAUNCH_BLOCKING=1");
      #else
        throw std::runtime_error("Using CudaWrapperNode without UVM is not allowed.");
      #endif

      if(!Kokkos::HostSpace::execution_space::is_initialized()) {
        if(NumNUMA>0 && NumCoresPerNUMA>0)
          Kokkos::HostSpace::execution_space::initialize ( NumThreads, NumNUMA, NumCoresPerNUMA );
        else if (NumNUMA > 0)
          Kokkos::HostSpace::execution_space::initialize ( NumThreads, NumNUMA );
        else
          Kokkos::HostSpace::execution_space::initialize ( NumThreads );
      }
      Kokkos::Cuda::SelectDevice select_device(Device);
      if(!Kokkos::Cuda::is_initialized())
        Kokkos::Cuda::initialize(select_device);
    }
    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::Cuda>::name() {
      return std::string("Cuda/Wrapper");
    }

#endif
  }
}



