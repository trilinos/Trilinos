// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RBGEN_UTILS_H
#define RBGEN_UTILS_H

class Epetra_MultiVector;

namespace RBGen {
  
  //! Method for computing the angle between two subspaces represented as Epetra_MultiVector bases.
  double BasisAngle( const Epetra_MultiVector& S, const Epetra_MultiVector& T );
  
} // end namespace RBGen

#endif // RBGEN_UTILS_H
