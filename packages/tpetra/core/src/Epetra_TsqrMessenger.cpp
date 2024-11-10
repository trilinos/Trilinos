// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

///
/// \file Epetra_TsqrMessenger.cpp
/// \brief Implementation of a function for wrapping Epetra_Comm in a
///   communicator wrapper that TSQR can use.
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

#include "TpetraCore_config.h"

#if defined(HAVE_TPETRA_EPETRA) && defined(HAVE_TPETRA_TSQR)

#include <Epetra_TsqrMessenger.hpp>
#include <Epetra_ConfigDefs.h> // EPETRA_MPI
#include <Teuchos_ConfigDefs.hpp> // HAVE_TEUCHOS_MPI

// Make sure that both Epetra and Teuchos agree on whether MPI exists.
#if defined(EPETRA_MPI) && ! defined(HAVE_TEUCHOS_MPI)
#  error "The Epetra package defines EPETRA_MPI, but the Teuchos package does NOT define HAVE_TEUCHOS_MPI.  This means that one package thinks MPI exists, but the other package doesn't.  This may indicate that Trilinos was not configured correctly.  We can't continue building Trilinos in this case."
#elif ! defined(EPETRA_MPI) && defined(HAVE_TEUCHOS_MPI)
#  error "The Epetra package does not define EPETRA_MPI, but the Teuchos package DOES define HAVE_TEUCHOS_MPI.  This means that one package thinks MPI exists, but the other package doesn't.  This may indicate that Trilinos was not configured correctly.  We can't continue building Trilinos in this case."
#endif

// Include Epetra's MPI wrappers.
#ifdef EPETRA_MPI
#  include <Epetra_MpiComm.h>
#endif // EPETRA_MPI
#include <Epetra_SerialComm.h>

// Include Teuchos' MPI wrappers.
#ifdef HAVE_TEUCHOS_MPI
#  include <Teuchos_DefaultMpiComm.hpp>
#endif // HAVE_TEUCHOS_MPI
#include <Teuchos_DefaultSerialComm.hpp>


namespace TSQR {
  namespace Epetra {

    Teuchos::RCP<const Teuchos::Comm<int> >
    extractTeuchosComm (const Teuchos::RCP<const Epetra_Comm>& epetraComm)
    {
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::rcp_dynamic_cast;
      using Teuchos::rcp_implicit_cast;

#ifdef EPETRA_MPI
      RCP<const Epetra_MpiComm> epetraMpiComm =
        rcp_dynamic_cast<const Epetra_MpiComm> (epetraComm, false);
      if (! epetraMpiComm.is_null ()) {
        // It's OK to extract and use the raw MPI_Comm object, since
        // an Epetra_MpiComm doesn't own (is not responsible for
        // calling MPI_Comm_free on) its underlying MPI_Comm object.
        MPI_Comm rawMpiComm = epetraMpiComm->Comm ();
        RCP<const Teuchos::MpiComm<int> > teuchosMpiComm =
          rcp (new Teuchos::MpiComm<int> (rawMpiComm));
        return rcp_implicit_cast<const Teuchos::Comm<int> > (teuchosMpiComm);
      }
#endif // EPETRA_MPI

      RCP<const Epetra_SerialComm> epetraSerialComm =
        rcp_dynamic_cast<const Epetra_SerialComm> (epetraComm, false);
      if (! epetraSerialComm.is_null ()) {
        RCP<const Teuchos::SerialComm<int> > teuchosSerialComm =
          rcp (new Teuchos::SerialComm<int> ());
        return rcp_implicit_cast<const Teuchos::Comm<int> > (teuchosSerialComm);
      }

      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument, "The input Epetra_Comm object is neither "
        "an Epetra_MpiComm nor an Epetra_SerialComm.  This means that we don't"
        " know how to convert it into a Teuchos::Comm<int> instance.  You are "
        "probably seeing this error message as a result of trying to invoke "
        "TSQR (the Tall Skinny QR factorization) on an Epetra_MultiVector, "
        "probably as a block orthogonalization method in either Anasazi or "
        "Belos.  TSQR currently only works on an Epetra_MultiVector if the "
        "latter's communicator is either an Epetra_MpiComm (wraps MPI) or an "
        "Epetra_SerialComm (\"stub\" communicator for a single process without "
        "MPI).");
    }
  } // namespace Epetra
} // namespace TSQR

#endif // defined(HAVE_TPETRA_EPETRA) && defined(HAVE_TPETRA_TSQR)

