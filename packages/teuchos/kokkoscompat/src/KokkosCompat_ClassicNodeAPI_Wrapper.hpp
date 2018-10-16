#ifndef KOKKOSCOMPAT_CLASSICNODEAPI_WRAPPER_HPP
#define KOKKOSCOMPAT_CLASSICNODEAPI_WRAPPER_HPP

#include "TeuchosKokkosCompat_config.h"
#include "Kokkos_Core.hpp"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
//
// Dear users: These are just forward declarations.  Please skip
// over them and go down to KokkosDeviceWrapperNode below.  Thank
// you.
//
namespace Teuchos {
  class ParameterList;
} // namespace Teuchos
#endif // DOXYGEN_SHOULD_SKIP_THIS

namespace Kokkos {
namespace Compat {

/// \brief Node that wraps a new Kokkos execution space.
///
/// \tparam ExecutionSpace The type of the Kokkos execution space to wrap.
/// \tparam MemorySpace The Kokkos memory space in which to work.
///   Defaults to the default memory space of ExecutionSpace.
template<class ExecutionSpace,
         class MemorySpace = typename ExecutionSpace::memory_space>
class KokkosDeviceWrapperNode {
public:
  //! The Node's Kokkos execution space.
  typedef ExecutionSpace execution_space;
  //! The Node's Kokkos memory space.
  typedef MemorySpace memory_space;
  /// \brief The Node's Kokkos::Device specialization.
  ///
  /// This is just an (execution space, memory space) pair.
  typedef Kokkos::Device<execution_space, memory_space> device_type;

  /// \brief This is NOT a "classic" Node type.
  ///
  /// We will deprecate the "classic" Node types with the 11.14
  /// release of Trilinos, and remove them entirely with the 12.0
  /// release.  This Node type is safe to use.
  static constexpr bool classic = false;

  /// \brief Constructor (that takes a Teuchos::ParameterList).
  ///
  /// \param [in/out] params List of Node configuration parameters.
  ///   If empty, we use defaults.
  KokkosDeviceWrapperNode (Teuchos::ParameterList& /* params */) {}

  //! Default constructor (sets default parameters).
  KokkosDeviceWrapperNode () {}

  //! Get a filled-in set of parameters for Node, with their default values.
  static Teuchos::ParameterList getDefaultParameters ();

  void sync () const {
    execution_space::fence ();
  }

  //! Human-readable name of this Node.
  static std::string name ();
};

#ifdef KOKKOS_ENABLE_CUDA
  typedef KokkosDeviceWrapperNode<Kokkos::Cuda> KokkosCudaWrapperNode;
#endif

#ifdef KOKKOS_ENABLE_OPENMP
  typedef KokkosDeviceWrapperNode<Kokkos::OpenMP> KokkosOpenMPWrapperNode;
#endif

#ifdef KOKKOS_ENABLE_THREADS
  typedef KokkosDeviceWrapperNode<Kokkos::Threads> KokkosThreadsWrapperNode;
#endif

#ifdef KOKKOS_ENABLE_SERIAL
  typedef KokkosDeviceWrapperNode<Kokkos::Serial> KokkosSerialWrapperNode;
#endif // KOKKOS_ENABLE_SERIAL

  // The above definitions / initializations of class (static)
  // variables need to precede the first use of these variables.
  // Otherwise, CUDA 7.5 with GCC 4.8.4 emits a warning ("explicit
  // specialization of member ... must precede its first use").

} // namespace Compat
} // namespace Kokkos
#endif
