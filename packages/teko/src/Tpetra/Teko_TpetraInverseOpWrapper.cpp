// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_TpetraInverseOpWrapper.hpp"

using namespace Teuchos;

namespace Teko {
namespace TpetraHelpers {

void TpetraInverseOpWrapper::apply(const Tpetra::MultiVector<ST, LO, GO, NT>& X,
                                   Tpetra::MultiVector<ST, LO, GO, NT>& Y, Teuchos::ETransp mode,
                                   ST alpha, ST beta) const {
  TpetraOperatorWrapper::apply(X, Y, mode, alpha, beta);
}

void TpetraInverseOpWrapper::applyInverse(const Tpetra::MultiVector<ST, LO, GO, NT>& X,
                                          Tpetra::MultiVector<ST, LO, GO, NT>& Y,
                                          Teuchos::ETransp mode, ST alpha, ST beta) const {
  TpetraOperatorWrapper::apply(X, Y, mode, alpha, beta);
}

}  // namespace TpetraHelpers
}  // end namespace Teko
