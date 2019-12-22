#include "Tsqr_Test_MpiAndKokkosScope.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Kokkos_Core.hpp"
#include <iostream>

namespace TSQR {
namespace Test {

MpiAndKokkosScope::
MpiAndKokkosScope(int* argc, char*** argv) :
  blackHole_(static_cast<std::ostream*>(new Teuchos::oblackholestream)),
  mpiScope_(new Teuchos::GlobalMPISession(argc, argv, blackHole_.get())),
  kokkosScope_(new Kokkos::ScopeGuard(*argc, *argv))
{}

Teuchos::RCP<const Teuchos::Comm<int>>
MpiAndKokkosScope::getComm() const {
  return Teuchos::DefaultComm<int>::getComm();
}

std::ostream& MpiAndKokkosScope::outStream() const {
  // Only Process 0 gets to write to cout and cerr.  The other MPI
  // processes send their output to a "black hole" (something that
  // acts like /dev/null).
  return getComm()->getRank() == 0 ? std::cout :
    static_cast<std::ostream&>(*blackHole_);
}

std::ostream& MpiAndKokkosScope::errStream() const {
  return getComm()->getRank() == 0 ? std::cerr :
    static_cast<std::ostream&>(*blackHole_);
}

} // namespace Test
} // namespace TSQR
