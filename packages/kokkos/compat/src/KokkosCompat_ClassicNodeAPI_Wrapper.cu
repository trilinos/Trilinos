#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>
#include <Kokkos_Core.hpp>
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
        if(HostSpace::execution_space::is_initialized()) {
          // make sure that no Actual DeviceWrapper node of the mirror_device_type is in use
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
