// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef AMESOS2_CONCRETEMATRIXADAPTER_DECL_HPP
#define AMESOS2_CONCRETEMATRIXADAPTER_DECL_HPP

namespace Amesos2 {

  template <class Matrix>
  class ConcreteMatrixAdapter {};

}

#include "Amesos2_TpetraCrsMatrix_MatrixAdapter_decl.hpp"
#include "Amesos2_KokkosCrsMatrix_MatrixAdapter_decl.hpp"

#ifdef HAVE_AMESOS2_EPETRA
#  include "Amesos2_EpetraCrsMatrix_MatrixAdapter_decl.hpp"
#endif

#endif  // AMESOS2_CONCRETEMATRIXADAPTER_DECL_HPP
