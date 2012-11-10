// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include "Xpetra_Utils.hpp"
#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

  std::string toString(Xpetra::UnderlyingLib lib) {
    if (lib == Xpetra::UseTpetra) {
      return "Tpetra";
    } else if (lib == Xpetra::UseEpetra) {
      return "Epetra";
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError, "lib != UseTpetra && lib != UseEpetra");
    }
  }

#ifdef HAVE_XPETRA_TPETRA

  Xpetra::LookupStatus toXpetra(Tpetra::LookupStatus ls) {

    if (ls == Tpetra::AllIDsPresent)
      return Xpetra::AllIDsPresent;
    if (ls == Tpetra::IDNotPresent)
      return Xpetra::IDNotPresent;

    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::RuntimeError, "Unknown LookupStatus");

  }

  Tpetra::ProfileType toTpetra(Xpetra::ProfileType pt) {

    if (pt == Xpetra::StaticProfile)
      return Tpetra::StaticProfile;
    if (pt == Xpetra::DynamicProfile)
      return Tpetra::DynamicProfile;

    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::RuntimeError, "Unknown ProfileType");

  }

  Tpetra::OptimizeOption toTpetra(Xpetra::OptimizeOption os) {

    if (os == Xpetra::DoOptimizeStorage)
      return Tpetra::DoOptimizeStorage;
    if (os == Xpetra::DoNotOptimizeStorage)
      return Tpetra::DoNotOptimizeStorage;

    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::RuntimeError, "Unknown OptimizeOption");

  }

  Tpetra::CombineMode toTpetra(Xpetra::CombineMode cm) {

    if (cm == Xpetra::ADD)
      return Tpetra::ADD;

    if (cm == Xpetra::INSERT)
      return Tpetra::INSERT;

    if (cm == Xpetra::ABSMAX)
      return Tpetra::ABSMAX;

    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::RuntimeError, "Cannot convert Xpetra::CombineMode to Tpetra::CombineMode: unsupported CombineMode.");

  }

  Tpetra::LocalGlobal toTpetra(LocalGlobal lg) {

    if (lg == Xpetra::LocallyReplicated)
      return Tpetra::LocallyReplicated;
    if (lg == Xpetra::GloballyDistributed)
      return Tpetra::GloballyDistributed;

    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::RuntimeError, "Unknown LocalGlobal");

  }

#endif // HAVE_XPETRA_TPETRA

#ifdef HAVE_XPETRA_EPETRA

  Xpetra::LookupStatus toXpetra(int ls) {
    // This function is used only to convert the return value of Epetra_BlockMap::RemoteIDList() and Epetra_DirectoryBase::GetDirectoryEntries().
    // In the current implementation of Epetra (01/2012), these functions returns 0 (= success) or 1 (= a GID is not present on any processor).

    if (ls == 0)
      return Xpetra::AllIDsPresent;
    else if (ls == 1)
      return Xpetra::IDNotPresent;

    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::RuntimeError, "Epetra returned the following error code: " << ls << ". Xpetra do not know how to interpret this error code.");
  }

  bool toEpetra(Xpetra::ProfileType pt) {

    if (pt == Xpetra::StaticProfile)
      return true;
    if (pt == Xpetra::DynamicProfile)
      return false;

    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::RuntimeError, "Unknown ProfileType");
  }

  bool toEpetra(Xpetra::OptimizeOption os) {

    if (os == Xpetra::DoOptimizeStorage)
      return true;
    if (os == Xpetra::DoNotOptimizeStorage)
      return false;

    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::RuntimeError, "Unknown OptimizeOption");

  }

  Epetra_CombineMode toEpetra(Xpetra::CombineMode cm) {
    // Note: all the CombineMode are not supported.
    // According to Chris B., the behavior in Tpetra is the same as Epetra but I prefer to limit my tests for now.
    // See also the discussion of March 22 on the Tpetra developers mailing list.

    if (cm == Xpetra::ADD)
      return Add;
    if (cm == Xpetra::INSERT)
      return Insert;
    if (cm == Xpetra::ABSMAX)
      return AbsMax;

    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::RuntimeError, "Cannot convert Xpetra::CombineMode to Epetra_CombineMode: unsupported CombineMode.");

  }

#endif // HAVE_XPETRA_EPETRA

} // namespace xpetra
