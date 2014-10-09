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
///
/// Method for fetching TSQR::MessengerBase instance for use in the
/// Epetra_MultiVector adaptor for TSQR.
///
/// \note (mfh 27 Oct 2010) This file is in Tpetra (rather than
/// Epetra, where it would seem to belong) as a temporary fix.
/// Otherwise, Epetra would need an optional package dependency on
/// Teuchos and Kokkos, which would break third-party code linking to
/// the Epetra library.  Third-party code should use FIND_PACKAGE on
/// Trilinos to get the correct list of libraries against which to
/// link, but we make this easy temporary fix now so they have time to
/// fix their build systems later.
///

#ifndef __Epetra_TsqrMessenger_hpp
#define __Epetra_TsqrMessenger_hpp

#include <Tpetra_ConfigDefs.hpp>

#if defined(HAVE_TPETRA_EPETRA) && defined(HAVE_TPETRA_KOKKOSTSQR)

#include <Epetra_ConfigDefs.h> // EPETRA_MPI
#include <Epetra_Comm.h>
#ifdef EPETRA_MPI
#  include <Epetra_MpiComm.h>
// #include <Epetra_MpiSmpComm.h>
#  include <Tsqr_MpiMessenger.hpp>
#endif // EPETRA_MPI
#include <Teuchos_RCP.hpp>
#include <Tsqr_TrivialMessenger.hpp>

#include <algorithm>
#include <utility> // std::pair

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Epetra {

#ifdef EPETRA_MPI
    /// If the input Epetra_Comm is really an Epetra_MpiComm, return its
    /// raw MPI_Comm MPI communicator object.  Otherwise, return
    /// MPI_COMM_NULL.  (The Epetra_Comm interface doesn't define sends
    /// and receives, which TSQR needs.  That's why TSQR wants the raw
    /// MPI_COMM object.)
    ///
    /// \return (The MPI_Comm object, and whether it's valid)
    static std::pair< MPI_Comm, bool >
    extractRawMpiComm (const Teuchos::RCP< const Epetra_Comm >& pComm)
    {
      using Teuchos::RCP;
      using Teuchos::rcp_dynamic_cast;

      MPI_Comm rawMpiComm = MPI_COMM_NULL;
      bool haveMpiComm = false;

      RCP< const Epetra_MpiComm > pMpiComm =
        rcp_dynamic_cast< const Epetra_MpiComm > (pComm, false);
      if (pMpiComm.get() == NULL)
        haveMpiComm = false;
      else
        {
          rawMpiComm = pMpiComm->Comm();
          haveMpiComm = true;
        }
      return std::make_pair (rawMpiComm, haveMpiComm);
    }
#endif // EPETRA_MPI


    /// Given a pointer to an Epetra_Comm object, return a pointer to
    /// a TSQR::MessengerBase< Datum > (actually, a pointer to an
    /// appropriate subclass thereof, depending on whether the given
    /// Epetra_Comm uses MPI).
    template< class Datum >
    Teuchos::RCP< TSQR::MessengerBase< Datum > >
    makeTsqrMessenger (const Teuchos::RCP< const Epetra_Comm >& pComm)
    {
      using Teuchos::RCP;
      using Teuchos::rcp_implicit_cast;
      typedef TSQR::MessengerBase< Datum > base_mess_type;

#ifdef EPETRA_MPI
      // If the Epetra_Comm wraps MPI, then extract the raw MPI_Comm
      // object and use it to construct a TSQR messenger object that
      // wraps MPI.  Otherwise, return a TSQR messenger object that
      // wraps trivial communication.
      std::pair< MPI_Comm, bool > results = extractRawMpiComm (pComm);
      if (results.second == true)
        {
          typedef TSQR::MPI::MpiMessenger< Datum > mess_type;

          RCP< mess_type > pMess (new mess_type (results.first));
          RCP< base_mess_type > pMessBase =
            rcp_implicit_cast< base_mess_type > (pMess);
          return pMessBase;
        }
      else
#endif // EPETRA_MPI
        {
          typedef TSQR::TrivialMessenger< Datum > mess_type;

          RCP< mess_type > pMess (new mess_type);
          RCP< base_mess_type > pMessBase =
            rcp_implicit_cast< base_mess_type > (pMess);
          return pMessBase;
        }
    }

  } // namespace Epetra
} // namespace TSQR

#endif // defined(HAVE_TPETRA_EPETRA) && defined(HAVE_TPETRA_KOKKOSTSQR)

#endif // __Epetra_TsqrMessenger_hpp

