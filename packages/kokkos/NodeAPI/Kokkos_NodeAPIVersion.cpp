#include "Kokkos_NodeAPIVersion.hpp"
#include "Trilinos_version.h"

namespace Kokkos {
  std::string NodeAPIVersion() { 
		return("Kokkos Node API in Trilinos " TRILINOS_VERSION_STRING);
	}
}
