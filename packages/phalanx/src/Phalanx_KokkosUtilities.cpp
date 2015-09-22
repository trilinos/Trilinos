#include "Phalanx_KokkosUtilities.hpp"
#include "Teuchos_Assert.hpp"

namespace PHX {

  void InitializeKokkosDevice(const int& num_threads)
  {
#if defined(PHX_KOKKOS_DEVICE_TYPE_CUDA)
    Kokkos::HostSpace::execution_space::initialize();
    Kokkos::Cuda::initialize();
    //std::cout << "CUDA device has been initialized" <<std::endl;
#elif defined(PHX_KOKKOS_DEVICE_TYPE_SERIAL)
    PHX::Device::initialize();
#else
    
    if (Kokkos::hwloc::available()) {
      //std::cout <<"hwloc"<<std::endl;
      const int num_hwloc_threads = Kokkos::hwloc::get_available_numa_count()
	* Kokkos::hwloc::get_available_cores_per_numa()
	* Kokkos::hwloc::get_available_threads_per_core()
	;
      TEUCHOS_TEST_FOR_EXCEPTION(num_threads > num_hwloc_threads, std::runtime_error,
				 "Error - PHX::InitializeKokkosDevice(num_threads) requested more threads than the node supports!");
    }

    PHX::Device::initialize(num_threads);
    //std::cout << "OpenMP/Pthreads  device has been initialized" <<std::endl;
    //std::cout <<"   number of threads= " << num_threads<<std::endl;
    //std::cout << "available threads: " << omp_get_max_threads() << std::endl;
#endif
  }
  
  void FinalizeKokkosDevice()
  {
#if defined(PHX_KOKKOS_DEVICE_TYPE_CUDA)
    Kokkos::Cuda::finalize();
    Kokkos::HostSpace::execution_space::finalize();
#else
    PHX::Device::finalize();
#endif
  }

  KokkosDeviceSession::KokkosDeviceSession()
  { PHX::InitializeKokkosDevice(); }

  KokkosDeviceSession::~KokkosDeviceSession()
  { 
    PHX::FinalizeKokkosDevice();
  }
}


