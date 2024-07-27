// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_EXTRACTMPICOMMFROMTEUCHOS_HPP
#define TPETRA_DETAILS_EXTRACTMPICOMMFROMTEUCHOS_HPP

/// \file Tpetra_Details_extractMpiCommFromTeuchos.hpp
/// \brief Declaration of Tpetra::Details::extractMpiCommFromTeuchos
///
/// \warning This file and its contents are implementation details of
///   Tpetra.  Users must not rely on them.

#include "TpetraCore_config.h"
#ifdef HAVE_TPETRACORE_MPI
#  include <mpi.h> // MPI_Comm
#endif // HAVE_TPETRACORE_MPI

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Teuchos {
  // forward declaration of Comm
  template<class OrdinalType> class Comm;
} // namespace Teuchos
#endif // NOT DOXYGEN_SHOULD_SKIP_THIS

namespace Tpetra {
namespace Details {

#ifdef HAVE_TPETRACORE_MPI
/// \brief Get the MPI_Comm out of the given Teuchos::Comm object.
///
/// \param comm [in] Communicator, wrapped in a Teuchos::Comm wrapper.
///   Must be an instance of Teuchos::MpiComm or Teuchos::SerialComm.
///
/// \return The wrapped MPI_Comm (if comm is a Teuchos::MpiComm), or
///   MPI_COMM_SELF (if comm is a Teuchos::SerialComm).
MPI_Comm
extractMpiCommFromTeuchos (const Teuchos::Comm<int>& comm);
#endif // HAVE_TPETRACORE_MPI

//! Is the given Comm a Teuchos::MpiComm<int> instance?
bool teuchosCommIsAnMpiComm (const Teuchos::Comm<int>& comm);

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_EXTRACTMPICOMMFROMTEUCHOS_HPP
