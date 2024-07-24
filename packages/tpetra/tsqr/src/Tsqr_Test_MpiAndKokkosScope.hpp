// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TSQR_TEST_MPIANDKOKKOSSCOPE_HPP
#define TSQR_TEST_MPIANDKOKKOSSCOPE_HPP

#include "Teuchos_RCP.hpp"
#include <memory>
#include <ostream>

namespace Kokkos {
class ScopeGuard;
} // namespace Kokkos

namespace Teuchos {
template<class OrdinalType> class Comm;
} // namespace Teuchos

namespace TSQR {
namespace Test {

class MpiScope {
public:
  MpiScope(int* argc, char*** argv);
  ~MpiScope();
};

// Scope guard for TSQR's tests, that automatically initializes and
// finalizes both MPI (if building with MPI enabled) and Kokkos.
class MpiAndKokkosScope {
public:
  MpiAndKokkosScope(int* argc, char*** argv);

  Teuchos::RCP<const Teuchos::Comm<int>> getComm() const;
  std::ostream& outStream() const;
  std::ostream& errStream() const;

private:
  static Teuchos::RCP<const Teuchos::Comm<int>> getDefaultComm();

  MpiScope mpiScope_;
  std::unique_ptr<std::ostream> blackHole_;
  Teuchos::RCP<const Teuchos::Comm<int>> comm_;
  // The only reason ever to handle a scope guard by pointer is for
  // implementation hiding via the "pImpl" (pointer to implementation)
  // idiom.
  std::unique_ptr<Kokkos::ScopeGuard> kokkosScope_;
};

} // namespace Test
} // namespace TSQR

#endif // TSQR_TEST_MPIANDKOKKOSSCOPE_HPP
