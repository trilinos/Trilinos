// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
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

