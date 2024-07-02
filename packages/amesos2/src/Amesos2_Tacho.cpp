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

#include "Amesos2_Tacho_decl.hpp"
#include "Amesos2_Tacho_def.hpp"
#include "Amesos2_ExplicitInstantiationHelpers.hpp"
#include "TpetraCore_ETIHelperMacros.h"

namespace Amesos2 {
#ifdef HAVE_AMESOS2_EPETRA
  AMESOS2_SOLVER_EPETRA_INST(TachoSolver);
#endif

#define AMESOS2_TACHO_LOCAL_INSTANT(S,LO,GO,N) \
template class Amesos2::TachoSolver<Tpetra::CrsMatrix<S, LO, GO, N>, \
                                Tpetra::MultiVector<S, LO, GO, N> >;

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(AMESOS2_TACHO_LOCAL_INSTANT)

  #define AMESOS2_KOKKOS_IMPL_SOLVER_NAME TachoSolver
  #include "Amesos2_Kokkos_Impl.hpp"
}

#endif  // HAVE_AMESOS2_EXPLICIT_INSTANTIATION
