#include "KokkosCompat_ClassicNodeAPI_Wrapper.hpp"
#include "Teuchos_ParameterList.hpp"


namespace Kokkos {
  namespace Compat {

#ifdef KOKKOS_ENABLE_THREADS
    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::Threads>::name () {
      return "Threads/Wrapper";
    }

#ifndef TEUCHOS_HIDE_DEPRECATED_CODE
    template<>
    TEUCHOS_DEPRECATED
    Teuchos::ParameterList
    KokkosDeviceWrapperNode<Kokkos::Threads>::getDefaultParameters () {
      return Teuchos::ParameterList ();
    }
#endif // TEUCHOS_HIDE_DEPRECATED_CODE
#endif

#ifdef KOKKOS_ENABLE_OPENMP
    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::OpenMP>::name () {
      return "OpenMP/Wrapper";
    }

#ifndef TEUCHOS_HIDE_DEPRECATED_CODE
    template<>
    TEUCHOS_DEPRECATED
    Teuchos::ParameterList
    KokkosDeviceWrapperNode<Kokkos::OpenMP>::getDefaultParameters () {
      return Teuchos::ParameterList ();
    }
#endif // TEUCHOS_HIDE_DEPRECATED_CODE
#endif

#ifdef KOKKOS_ENABLE_SERIAL
    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::Serial>::name () {
      return "Serial/Wrapper";
    }

#ifndef TEUCHOS_HIDE_DEPRECATED_CODE
    template<>
    TEUCHOS_DEPRECATED
    Teuchos::ParameterList
    KokkosDeviceWrapperNode<Kokkos::Serial>::getDefaultParameters () {
      return Teuchos::ParameterList ();
    }
#endif // TEUCHOS_HIDE_DEPRECATED_CODE
#endif // KOKKOS_ENABLE_SERIAL

#ifdef KOKKOS_ENABLE_CUDA
    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::Cuda>::name() {
      return std::string("Cuda/Wrapper");
    }

#ifndef TEUCHOS_HIDE_DEPRECATED_CODE
    template<>
    TEUCHOS_DEPRECATED
    Teuchos::ParameterList
    KokkosDeviceWrapperNode<Kokkos::Cuda>::getDefaultParameters () {
      return Teuchos::ParameterList ();
    }
#endif // TEUCHOS_HIDE_DEPRECATED_CODE
#endif // KOKKOS_ENABLE_CUDA

  } // namespace Compat
} // namespace Kokkos



