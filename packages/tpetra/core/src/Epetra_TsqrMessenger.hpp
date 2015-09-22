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

///
/// \file Epetra_TsqrMessenger.hpp
/// \brief Function for wrapping Epetra_Comm in a communicator wrapper
///   that TSQR can use.
///
/// \note (mfh 27 Oct 2010) This file is in Tpetra (rather than
///   Epetra, where it would seem to belong) as a temporary fix.
///   Otherwise, Epetra would need an optional package dependency on
///   Teuchos and Kokkos, which would break third-party code linking
///   to the Epetra library.  Third-party code should use FIND_PACKAGE
///   on Trilinos to get the correct list of libraries against which
///   to link, but we make this easy temporary fix now so they have
///   time to fix their build systems later.
///

#ifndef EPETRA_TSQRMESSENGER_HPP
#define EPETRA_TSQRMESSENGER_HPP

#include <Tpetra_ConfigDefs.hpp>

#if defined(HAVE_TPETRA_EPETRA) && defined(HAVE_TPETRA_TSQR)

// Include Epetra's MPI wrappers.
#include <Epetra_Comm.h>

// Include Teuchos' MPI wrappers.
#include <Teuchos_Comm.hpp>

#include <Teuchos_RCP.hpp>
#include <Tsqr_TeuchosMessenger.hpp>


namespace TSQR {
  namespace Epetra {

    /// \brief Return the Teuchos communicator corresponding to the
    ///   given Epetra communicator.
    ///
    /// If the input Epetra_Comm is really an Epetra_MpiComm, return a
    /// Teuchos-wrapped version of its raw MPI_Comm MPI communicator
    /// object.  Otherwise, return a Teuchos::SerialComm instance.  It
    /// should be one of these two things, but if it's not, this
    /// function throws std::invalid_argument.
    Teuchos::RCP<const Teuchos::Comm<int> >
    extractTeuchosComm (const Teuchos::RCP<const Epetra_Comm>& epetraComm);

    //! Wrap the given Epetra_Comm in an object that TSQR understands.
    template<class Datum>
    Teuchos::RCP<TSQR::MessengerBase<Datum> >
    makeTsqrMessenger (const Teuchos::RCP<const Epetra_Comm>& pComm)
    {
      typedef TSQR::MessengerBase<Datum> base_mess_type;
      typedef TSQR::TeuchosMessenger<Datum> mess_type;

      Teuchos::RCP<mess_type> pMess =
        Teuchos::rcp (new mess_type (extractTeuchosComm (pComm)));
      return Teuchos::rcp_implicit_cast<base_mess_type> (pMess);
    }
  } // namespace Epetra
} // namespace TSQR

#endif // defined(HAVE_TPETRA_EPETRA) && defined(HAVE_TPETRA_TSQR)

#endif // EPETRA_TSQRMESSENGER_HPP

