// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_EPETRAVECTOR_FWD_HPP
#define XPETRA_EPETRAVECTOR_FWD_HPP

namespace Xpetra {
template <class GO, class NO>
class EpetraVectorT;
typedef EpetraVectorT<int, EpetraNode> EpetraVector;
}  // namespace Xpetra

#ifndef XPETRA_EPETRAVECTOR_SHORT
#define XPETRA_EPETRAVECTOR_SHORT
#endif

#endif  // XPETRA_EPETRAVECTOR_FWD_HPP
