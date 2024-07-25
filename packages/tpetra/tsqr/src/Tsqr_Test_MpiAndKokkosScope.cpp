// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tsqr_Test_MpiAndKokkosScope.hpp"
#include "Kokkos_Core.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_CommHelpers.hpp"
#ifdef HAVE_MPI
#  include "Teuchos_DefaultMpiComm.hpp"
#  include "Teuchos_Assert.hpp"
#else
#  include "Teuchos_DefaultSerialComm.hpp"
#endif // HAVE_MPI
#include <iostream>
#include <sstream>

namespace TSQR {
namespace Test {

#ifdef HAVE_MPI
MpiScope::MpiScope(int* argc, char*** argv) {
  (void) MPI_Init(argc, argv);
}
MpiScope::~MpiScope() {
  (void) MPI_Finalize();
}
#else
MpiScope::MpiScope(int*, char***) {}
MpiScope::~MpiScope() {}
#endif // HAVE_MPI

Teuchos::RCP<const Teuchos::Comm<int>>
MpiAndKokkosScope::getDefaultComm()
{
#ifdef HAVE_MPI
  int initialized = 0;
  (void) MPI_Initialized(&initialized);
  TEUCHOS_ASSERT( initialized == 1 );

  using comm_type = Teuchos::MpiComm<int>;
  const auto comm = Teuchos::rcp(new comm_type(MPI_COMM_WORLD));
#else
  using comm_type = Teuchos::SerialComm<int>;
  const auto comm = Teuchos::rcp(new comm_type);
#endif // HAVE_MPI

  return comm;
}

MpiAndKokkosScope::
MpiAndKokkosScope(int* argc, char*** argv) :
  mpiScope_(argc, argv),
  blackHole_(new Teuchos::oblackholestream),
  comm_(getDefaultComm()),
  kokkosScope_(new Kokkos::ScopeGuard(*argc, *argv))
{}

Teuchos::RCP<const Teuchos::Comm<int>>
MpiAndKokkosScope::getComm() const {
  return comm_;
}

std::ostream& MpiAndKokkosScope::outStream() const {
  // Only Process 0 gets to write to cout and cerr.  The other MPI
  // processes send their output to a "black hole" (something that
  // acts like /dev/null).
  return comm_->getRank() == 0 ? std::cout :
    static_cast<std::ostream&>(*blackHole_);
}

std::ostream& MpiAndKokkosScope::errStream() const {
  return comm_->getRank() == 0 ? std::cerr :
    static_cast<std::ostream&>(*blackHole_);
}

} // namespace Test
} // namespace TSQR
