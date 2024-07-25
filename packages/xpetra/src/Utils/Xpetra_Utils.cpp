// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

#endif  // HAVE_XPETRA_TPETRA

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

#endif  // HAVE_XPETRA_EPETRA

}  // namespace Xpetra
