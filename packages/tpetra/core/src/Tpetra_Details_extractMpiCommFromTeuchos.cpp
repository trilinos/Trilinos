// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Details_extractMpiCommFromTeuchos.hpp"

#ifdef HAVE_TPETRACORE_MPI
#  include "Teuchos_DefaultMpiComm.hpp"
#  include "Teuchos_DefaultSerialComm.hpp"
#  include <stdexcept>
#endif // HAVE_TPETRACORE_MPI

namespace Tpetra {
namespace Details {

#ifdef HAVE_TPETRACORE_MPI
MPI_Comm
extractMpiCommFromTeuchos (const Teuchos::Comm<int>& comm)
{
  using ::Teuchos::MpiComm;
  using ::Teuchos::SerialComm;

  const MpiComm<int>* mpiComm = dynamic_cast<const MpiComm<int>* > (&comm);
  if (mpiComm != NULL) { // It's an MpiComm; extract the MPI_Comm
    MPI_Comm rawComm = * (mpiComm->getRawMpiComm ());
    return rawComm;
  }
  else {
    const SerialComm<int>* serialComm =
      dynamic_cast<const SerialComm<int>* > (&comm);
    if (serialComm != NULL) {
      return MPI_COMM_SELF; // single-process comm including this process
    }
    else {
      throw std::invalid_argument ("Tpetra::Details::extractMpiCommFromTeuchos: "
                                   "Input Teuchos::Comm is "
                                   "neither a Teuchos::MpiComm, "
                                   "nor a Teuchos::SerialComm.  "
                                   "As a result, I don't know "
                                   "how to get the MPI_Comm out of it.");
    }
  }
}
#endif // HAVE_TPETRACORE_MPI

#ifdef HAVE_TPETRACORE_MPI
bool teuchosCommIsAnMpiComm (const Teuchos::Comm<int>& comm)
{
  const Teuchos::MpiComm<int>* mpiComm =
    dynamic_cast<const Teuchos::MpiComm<int>* > (&comm);
  return mpiComm != nullptr;
}
#else
bool teuchosCommIsAnMpiComm (const Teuchos::Comm<int>&)
{
  return false;
}
#endif // HAVE_TPETRACORE_MPI

} // namespace Details
} // namespace Tpetra

