// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "TpetraCore_config.h"
#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION

#include "Stokhos_Tpetra_ETI_Helpers_UQ_PCE_Cuda.hpp"

#include "Tpetra_CrsMatrix_UQ_PCE.hpp"
#include "Tpetra_CrsMatrix_UQ_PCE_def.hpp"
#include "Kokkos_ArithTraits_MP_Vector.hpp"
#include "Tpetra_Details_getDiagCopyWithoutOffsets_def.hpp"

#define TPETRA_LOCAL_INST_N_1(N) \
  INSTANTIATE_TPETRA_UQ_PCE_N(TPETRA_CRSMATRIX_UQ_PCE_SPEC, N)

#define TPETRA_LOCAL_INST_N_2(N) \
  INSTANTIATE_TPETRA_UQ_PCE_N(TPETRA_CRSMATRIX_INSTANT, N)

INSTANTIATE_TPETRA_UQ_PCE(TPETRA_LOCAL_INST_N_1)

INSTANTIATE_TPETRA_UQ_PCE(TPETRA_LOCAL_INST_N_2)

#endif // HAVE_TPETRA_EXPLICIT_INSTANTIATION
