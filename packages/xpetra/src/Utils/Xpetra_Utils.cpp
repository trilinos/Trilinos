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
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError, "lib != UseTpetra");
  }
}

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

}  // namespace Xpetra
