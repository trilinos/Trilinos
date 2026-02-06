// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef AMESOS2_EXPLICITINSTANTIATIONHELPER_HPP
#define AMESOS2_EXPLICITINSTANTIATIONHELPER_HPP

#include <Tpetra_CrsMatrix.hpp>

#define AMESOS2_SOLVER_TPETRA_INST(SOLVERNAME,S,LO,GO) \
  template class SOLVERNAME<Tpetra::CrsMatrix<S, LO, GO>, Tpetra::MultiVector<S, LO, GO> >

#endif  // AMESOS2_EXPLICITINSTANTIATIONHELPER_HPP
