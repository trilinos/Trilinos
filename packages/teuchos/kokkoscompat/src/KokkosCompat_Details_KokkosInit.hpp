#ifndef KOKKOSCOMPAT_DETAILS_KOKKOSINIT_HPP
#define KOKKOSCOMPAT_DETAILS_KOKKOSINIT_HPP

#include "TeuchosKokkosCompat_config.h"

namespace Kokkos {
  struct InitArguments; // forward declaration
}

namespace Teuchos {
  class ParameterList; // forward declaration
}

namespace KokkosCompat {
namespace Details {

//! Has Kokkos::initialize been called yet in this run?
bool
isKokkosInitialized ();

/// \brief If Kokkos::initialize has not yet been called, call it,
///   passing in the given command-line arguments.
///
/// Kokkos::initialize initializes all enabled Kokkos execution
/// spaces.
void
initializeKokkos (int& narg, char* arg[]);

/// \brief If Kokkos::initialize has not yet been called, call it.
///   Attempt to get command-line arguments from
///   Teuchos::GlobalMPISession.
///
/// Kokkos::initialize initializes all enabled Kokkos execution
/// spaces.
void
initializeKokkos ();

/// \brief If Kokkos::initialize has not yet been called, call it.
///   Use the version that takes arguments through Kokkos' struct.
///
/// Kokkos::initialize initializes all enabled Kokkos execution
/// spaces.
void
initializeKokkos (const Kokkos::InitArguments& args);

void
getNodeParameters (int& curNumThreads,
                   int& curNumNUMA,
                   int& curNumCoresPerNUMA,
                   int& curDevice,
                   bool& verbose,
                   const Teuchos::ParameterList& params);

Teuchos::ParameterList
getDefaultNodeParameters ();

void
setUpEnvironmentForCuda ();

} // namespace Details
} // namespace KokkosCompat

#endif // KOKKOSCOMPAT_DETAILS_KOKKOSINIT_HPP
