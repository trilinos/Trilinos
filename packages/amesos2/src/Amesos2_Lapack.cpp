// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER



#include "Amesos2_config.h"
#ifdef HAVE_AMESOS2_EXPLICIT_INSTANTIATION

#include "Amesos2_Lapack_decl.hpp"
#include "Amesos2_Lapack_def.hpp"
#include "Amesos2_ExplicitInstantiationHelpers.hpp"
#include "TpetraCore_ETIHelperMacros.h"

#ifdef HAVE_AMESOS2_EPETRA
namespace Amesos2 {
  AMESOS2_SOLVER_EPETRA_INST(Lapack);
}
#endif


#define AMESOS2_LAPACK_LOCAL_INSTANT(S,LO,GO,N) \
  template class Amesos2::Lapack<Tpetra::CrsMatrix<S, LO, GO, N>, \
                                  Tpetra::MultiVector<S, LO, GO, N> >;

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(AMESOS2_LAPACK_LOCAL_INSTANT)

#endif  // HAVE_AMESOS2_EXPLICIT_INSTANTIATION
