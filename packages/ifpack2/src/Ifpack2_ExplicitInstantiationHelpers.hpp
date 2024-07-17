// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_EXPLICITINSTANTIATIONHELPERS_HPP
#define TPETRA_EXPLICITINSTANTIATIONHELPERS_HPP

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_CrsGraph.hpp>

#define IFPACK2_INST(CLASSNAME,S,LO,GO,NO)                       \
  template class CLASSNAME<Tpetra::CrsMatrix<S, LO, GO, NO> >;

#define IFPACK2_INST_GRAPH(CLASSNAME,LO,GO) \
  template class CLASSNAME<Tpetra::CrsGraph<LO, GO> >;

#define IFPACK2_CLASS_CrsMatrix_float_int_int_defaultNode_defaultOps(CLASSNAME) \
  IFPACK2_INST(CLASSNAME,float,int,int)

#define IFPACK2_CLASS_CrsMatrix_float_short_int_defaultNode_defaultOps(CLASSNAME) \
  IFPACK2_INST(CLASSNAME,float,short,int)

#define IFPACK2_CLASS_CrsMatrix_double_int_int_defaultNode_defaultOps(CLASSNAME) \
  IFPACK2_INST(CLASSNAME,double,int,int)

#define IFPACK2_INSTANT_CRSMATRIX_FLOAT_DOUBLE_DEFAULTS(CLASSNAME) \
  IFPACK2_CLASS_CrsMatrix_double_int_int_defaultNode_defaultOps(CLASSNAME)

#define IFPACK2_INSTANT_CRSMATRIX_COMPLEX_DEFAULTS(CLASSNAME) \
  IFPACK2_INST(CLASSNAME,std::complex<double>,int,int) \
  IFPACK2_INST(CLASSNAME,std::complex<float>,int,int)

#endif // TPETRA_EXPLICITINSTANTIATIONHELPERS_HPP

