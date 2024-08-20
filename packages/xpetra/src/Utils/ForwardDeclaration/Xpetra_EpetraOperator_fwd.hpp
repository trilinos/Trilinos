// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_EPETRAOPERATOR_FWD_HPP
#define XPETRA_EPETRAOPERATOR_FWD_HPP

#include <Epetra_ConfigDefs.h>

namespace Xpetra {
template <class GO>
class EpetraOperatorT;
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
typedef EpetraOperatorT<int> EpetraOperator;
#endif
}  // namespace Xpetra

#ifndef XPETRA_EPETRAOPERATOR_SHORT
#define XPETRA_EPETRAOPERATOR_SHORT
#endif

#endif  // XPETRA_EPETRAOPERATOR_FWD_HPP
