// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_EPETRAEXPORT_FWD_HPP
#define XPETRA_EPETRAEXPORT_FWD_HPP

namespace Xpetra {
template <class GO, class NO>
class EpetraExportT;
typedef EpetraExportT<int, EpetraNode> EpetraExport;
}  // namespace Xpetra

#ifndef XPETRA_EPETRAEXPORT_SHORT
#define XPETRA_EPETRAEXPORT_SHORT
#endif

#endif  // XPETRA_EPETRAEXPORT_FWD_HPP
