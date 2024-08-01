// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_EPETRAINTVECTOR_FWD_HPP
#define XPETRA_EPETRAINTVECTOR_FWD_HPP

namespace Xpetra {
template <class GO, class NO>
class EpetraIntVectorT;
typedef EpetraIntVectorT<int, EpetraNode> EpetraIntVector;
}  // namespace Xpetra

#ifndef XPETRA_EPETRAINTVECTOR_SHORT
#define XPETRA_EPETRAINTVECTOR_SHORT
#endif

#endif  // XPETRA_EPETRAINTVECTOR_FWD_HPP
