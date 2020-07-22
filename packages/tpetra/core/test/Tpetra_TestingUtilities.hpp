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

#ifndef TPETRA_TESTINGUTILITIES_HPP_
#define TPETRA_TESTINGUTILITIES_HPP_

/// \file Tpetra_TestingUtilities.hpp
/// \brief Internal utilities for testing Tpetra.
///
/// \warning This header file and its contents are implementation
///   details of Tpetra.  Users must not rely on this file existing,
///   or on any contents of this file.

#include "Teuchos_UnitTestHarness.hpp"
#include "Tpetra_Core.hpp"
#include "TpetraCore_ETIHelperMacros.h"
#include "Teuchos_DefaultSerialComm.hpp"
#include "Teuchos_CommHelpers.hpp"

#define TPETRA_GLOBAL_SUCCESS_CHECK(out,comm,success)  \
  { \
    int tgscLclSuccess = success ? 1 : 0; \
    int tgscGblSuccess = 1; \
    Teuchos::reduceAll<int, int>(*comm, Teuchos::REDUCE_MIN, tgscLclSuccess, Teuchos::outArg (tgscGblSuccess)); \
    if (tgscGblSuccess == 1) { \
      out << "Succeeded on all processes!" << std::endl; \
    } else { \
      out << "FAILED on at least one process!" << std::endl; \
    } \
    TEST_EQUALITY_CONST(tgscGblSuccess, 1);  \
    success = (bool) tgscGblSuccess; \
  }


namespace Tpetra {
  namespace TestingUtilities {

    /// \brief Whether to test MPI functionality.
    ///
    /// \warning This variable is an implementation detail of Tpetra.
    ///   Users must not rely on this variable's name or contents, set
    ///   this variable, or rely on behavior resulting from setting
    ///   this variable.
    ///
    /// This variable is true by default.  Its value affects the
    /// behavior of getDefaultComm() (see below in this header file).
    bool testMpi = true;

    /// \brief Get the default communicator for Tpetra tests.
    ///
    /// \warning This function is an implementation detail of Tpetra.
    ///   Users must not call this function or rely on its behavior.
    ///
    /// If testMpi (see above in this header file) false, this
    /// function will return a Teuchos::SerialComm.  Otherwise, this
    /// fucntion will return the default communicator from
    /// Tpetra::getDefaultComm.  If Trilinos was built with MPI
    /// support, the resulting communicator will be a Teuchos::MpiComm
    /// that wraps <tt>MPI_COMM_WORLD</tt>.  Otherwise, it will be
    /// either a Teuchos::SerialComm, or a Teuchos::MpiComm that wraps
    /// <tt>MPI_COMM_SELF</tt>.
    Teuchos::RCP<const Teuchos::Comm<int> >
    getDefaultComm ()
    {
      if (testMpi) {
        return Tpetra::getDefaultComm ();
      }
      else {
        // Always return the same Comm instance.  Create it if it
        // hasn't already been created.  A function-static RCP has an
        // initial value of Teuchos::null.
        static Teuchos::RCP<const Teuchos::SerialComm<int> > serialComm_;
        if (serialComm_.is_null ()) {
          serialComm_ = Teuchos::rcp (new Teuchos::SerialComm<int> ());
        }
        return serialComm_;
      }
    }


  } // namespace TestingUtilities
} // namespace Tpetra

#endif // TPETRA_TESTINGUTILITIES_HPP_
