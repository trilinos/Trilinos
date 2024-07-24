// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_EpetraRowMatrix.hpp"
#if defined(HAVE_TPETRA_EPETRA)

#ifdef HAVE_TPETRACORE_MPI
#  include "Epetra_MpiComm.h"
#  include "Tpetra_Details_extractMpiCommFromTeuchos.hpp"
#else
#  include "Epetra_SerialComm.h"
#endif // HAVE_TPETRACORE_MPI

namespace Tpetra {
namespace Details {

#ifdef HAVE_TPETRACORE_MPI
std::shared_ptr<Epetra_Comm>
makeEpetraCommFromTeuchosComm (const Teuchos::Comm<int>& teuchosComm)
{
  using Tpetra::Details::extractMpiCommFromTeuchos;
  // NOTE (mfh 11 Oct 2017) Tpetra or Teuchos may free this MPI_Comm
  // before Epetra is done with it.  To ensure that this doesn't
  // happen, best practice is not to let the Epetra_Comm outlive the
  // input teuchosComm.
  MPI_Comm mpiComm = extractMpiCommFromTeuchos (teuchosComm);
  Epetra_MpiComm* epetraComm = new Epetra_MpiComm (mpiComm);
  return std::shared_ptr<Epetra_Comm> (static_cast<Epetra_Comm*> (epetraComm));
}
#else
std::shared_ptr<Epetra_Comm>
makeEpetraCommFromTeuchosComm (const Teuchos::Comm<int>&)
{
  return std::shared_ptr<Epetra_Comm> (static_cast<Epetra_Comm*> (new Epetra_SerialComm));
}
#endif // HAVE_TPETRACORE_MPI

} // namespace Details
} // namespace Tpetra

#endif // defined(HAVE_TPETRA_EPETRA)
