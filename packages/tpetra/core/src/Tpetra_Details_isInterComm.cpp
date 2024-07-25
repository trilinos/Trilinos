// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Details_isInterComm.hpp"
#include "Teuchos_Comm.hpp"
#ifdef HAVE_TPETRACORE_MPI
#  include "Tpetra_Details_extractMpiCommFromTeuchos.hpp"
#endif // HAVE_TPETRACORE_MPI
#include <stdexcept> // std::logic_error

namespace Tpetra {
namespace Details {

#ifdef HAVE_TPETRACORE_MPI
bool
isInterComm (const Teuchos::Comm<int>& comm)
{
  MPI_Comm rawMpiComm = extractMpiCommFromTeuchos (comm);
  int flag = 0;
  // This is a "local routine" (MPI 3.1, Section 6.6.1, p. 259).
  (void) MPI_Comm_test_inter (rawMpiComm, &flag);
  return flag != 0;
}

#else // NOT HAVE_TPETRACORE_MPI

bool
isInterComm (const Teuchos::Comm<int>& /* comm */ )
{
  return false;
}
#endif // HAVE_TPETRACORE_MPI

} // namespace Details
} // namespace Tpetra

