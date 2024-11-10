// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

#if !defined(TPETRA_ENABLE_DEPRECATED_CODE)
#error This file is deprecated due to Epetra removal and will be removed
#endif

#if defined(TPETRA_ENABLE_DEPRECATED_CODE) && defined(HAVE_TPETRA_EPETRA) && defined(HAVE_TPETRA_TSQR)

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
    TPETRA_DEPRECATED_MSG("epetra removal")
    Teuchos::RCP<const Teuchos::Comm<int> >
    extractTeuchosComm (const Teuchos::RCP<const Epetra_Comm>& epetraComm);

    //! Wrap the given Epetra_Comm in an object that TSQR understands.
    template<class Datum>
    TPETRA_DEPRECATED_MSG("epetra removal")
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

#endif // defined(TPETRA_ENABLE_DEPRECATED_CODE) && defined(HAVE_TPETRA_EPETRA) && defined(HAVE_TPETRA_TSQR)

#endif // EPETRA_TSQRMESSENGER_HPP

