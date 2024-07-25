// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_EPETRAMULTIVECTOR_FWD_HPP
#define XPETRA_EPETRAMULTIVECTOR_FWD_HPP

namespace Xpetra {
template <class GO, class NO>
class EpetraMultiVectorT;
typedef EpetraMultiVectorT<int, EpetraNode> EpetraMultiVector;
}  // namespace Xpetra

#ifndef XPETRA_EPETRAMULTIVECTOR_SHORT
#define XPETRA_EPETRAMULTIVECTOR_SHORT
#endif

#endif  // XPETRA_EPETRAMULTIVECTOR_FWD_HPP
