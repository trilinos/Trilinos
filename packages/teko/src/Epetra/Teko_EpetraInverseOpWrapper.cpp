// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_EpetraInverseOpWrapper.hpp"

using namespace Teuchos;

namespace Teko {
namespace Epetra {

int EpetraInverseOpWrapper::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
  return EpetraOperatorWrapper::ApplyInverse(X, Y);
}

int EpetraInverseOpWrapper::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
  return EpetraOperatorWrapper::Apply(X, Y);
}

}  // end namespace Epetra
}  // end namespace Teko
