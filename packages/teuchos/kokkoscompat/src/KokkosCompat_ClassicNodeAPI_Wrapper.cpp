#include "KokkosCompat_ClassicNodeAPI_Wrapper.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Kokkos {
  namespace Compat {

#ifdef KOKKOS_ENABLE_THREADS
    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::Threads>::name () {
      return "Threads/Wrapper";
    }

    template<>
    Teuchos::ParameterList
    KokkosDeviceWrapperNode<Kokkos::Threads>::getDefaultParameters () {
      return Teuchos::ParameterList ();
    }
#endif

#ifdef KOKKOS_ENABLE_OPENMP
    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::OpenMP>::name () {
      return "OpenMP/Wrapper";
    }

    template<>
    Teuchos::ParameterList
    KokkosDeviceWrapperNode<Kokkos::OpenMP>::getDefaultParameters () {
      return Teuchos::ParameterList ();
    }
#endif

#ifdef KOKKOS_ENABLE_SERIAL
    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::Serial>::name () {
      return "Serial/Wrapper";
    }

    template<>
    Teuchos::ParameterList
    KokkosDeviceWrapperNode<Kokkos::Serial>::getDefaultParameters () {
      return Teuchos::ParameterList ();
    }
#endif // KOKKOS_ENABLE_SERIAL

#ifdef KOKKOS_ENABLE_CUDA
    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::Cuda>::name() {
      return std::string("Cuda/Wrapper");
    }

    template<>
    Teuchos::ParameterList
    KokkosDeviceWrapperNode<Kokkos::Cuda>::getDefaultParameters () {
      return Teuchos::ParameterList ();
    }
#endif // KOKKOS_ENABLE_CUDA

  } // namespace Compat
} // namespace Kokkos



