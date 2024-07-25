// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_EPETRAMAP_FWD_HPP
#define XPETRA_EPETRAMAP_FWD_HPP

namespace Xpetra {
template <class GO, class NO>
class EpetraMapT;
typedef EpetraMapT<int, EpetraNode> EpetraMap;
}  // namespace Xpetra

#ifndef XPETRA_EPETRAMAP_SHORT
#define XPETRA_EPETRAMAP_SHORT
#endif

#endif  // XPETRA_EPETRAMAP_FWD_HPP
