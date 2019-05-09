/*
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
*/

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
