#include "Kokkos_Version.hpp"
#include "Trilinos_version.h"

namespace Kokkos {

  std::string Kokkos_Version() { 
		return("Kokkos in Trilinos " TRILINOS_VERSION_STRING);
  }

} // namespace Kokkos 
