// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_EPETRAIMPORT_FWD_HPP
#define XPETRA_EPETRAIMPORT_FWD_HPP

namespace Xpetra {
template <class GO, class NO>
class EpetraImportT;
typedef EpetraImportT<int, EpetraNode> EpetraImport;
}  // namespace Xpetra

#ifndef XPETRA_EPETRAIMPORT_SHORT
#define XPETRA_EPETRAIMPORT_SHORT
#endif

#endif  // XPETRA_EPETRAIMPORT_FWD_HPP
