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
