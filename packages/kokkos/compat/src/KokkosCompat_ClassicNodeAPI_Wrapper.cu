#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>
#include <Kokkos_HostSpace.hpp>
#include <cstdio>

namespace Kokkos {
  namespace Compat {

    // mfh 01 Jan 2014: These definitions of the class variable count
    // need to be inside the namespace.  Declaring them as "template<>
    // int Kokkos::Compat::KokkosCudaWrapperNode::count = 0" in the
    // global namespace is a C++11 extension and results in compiler
    // warnings with Clang 3.2 on MacOS X.
#ifdef KOKKOS_HAVE_CUDA
    template<>
    KokkosDeviceWrapperNode<Kokkos::Cuda>::~KokkosDeviceWrapperNode<Kokkos::Cuda>() {
      count--;
      if(count==0) {
        if(Cuda::host_mirror_device_type::is_initialized()) {
          // make sure that no Actual DeviceWrapper node of the mirror_device_type is in use
          if(KokkosDeviceWrapperNode<Cuda::host_mirror_device_type>::count==0) {
            //Don't try to kill me if HostSpace was already destroyed.
            //Typical reason: static global instance of node is used, which might get destroyed after
            //the static HostSpace is destroyed.
            if(Kokkos::NEVEREVERUSEMEIWILLFINDYOU::host_space_singleton_wrapper().size()>0)
              Cuda::host_mirror_device_type::finalize();
          }
        }
        if(Cuda::is_initialized())
          if(Kokkos::NEVEREVERUSEMEIWILLFINDYOU::host_space_singleton_wrapper().size()>0)
            Cuda::finalize();
      }
    }
    template<>
    void KokkosDeviceWrapperNode<Kokkos::Cuda>::init(int NumTeams, int NumThreads, int Device) {

      // Setting (currently) necessary environment variables for NVIDIA UVM
      #ifdef KOKKOS_USE_UVM
        char str[256];
        sprintf(str,"CUDA_VISIBLE_DEVICES=%i",Device);
        putenv(str);
        putenv("CUDA_LAUNCH_BLOCKING=1");
      #else
        throw std::runtime_error("Using CudaWrapperNode without UVM is not allowed.");
      #endif

      if(!Kokkos::Cuda::host_mirror_device_type::is_initialized())
        Kokkos::Cuda::host_mirror_device_type::initialize(NumTeams*NumThreads);
      Kokkos::Cuda::SelectDevice select_device(0);
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
