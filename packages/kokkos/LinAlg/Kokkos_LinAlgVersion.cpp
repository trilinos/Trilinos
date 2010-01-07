#include "Kokkos_LinAlgVersion.hpp"
#include "Trilinos_version.h"

namespace Kokkos {
  std::string LinAlgVersion() { 
		return("Kokkos Linear Algebra in Trilinos " TRILINOS_VERSION_STRING);
	}
}
