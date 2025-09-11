// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// ETP 4/11/20:  Do not include MueLu_ConfigDefs.hpp because it has Kokkos
// includes which need to come after the Stokhos includes.
#include "MueLu_config.hpp"



#include "MueLu_ExplicitInstantiation.hpp"
#include "Stokhos_ConfigDefs.h"

#if defined(HAVE_STOKHOS_MUELU) && defined(HAVE_MUELU_EXPLICIT_INSTANTIATION) && defined(HAVE_STOKHOS_SACADO)

#include "Stokhos_Tpetra_ETI_Helpers_MP_Vector_16_Cuda.hpp"
#include "Stokhos_MueLu_MP_Vector.hpp"

#include "MueLu_ScalarDroppingDistanceLaplacian_def.hpp"

#define MUELU_INST_S_LO_GO_N(S, LO, GO, N) \
  MUELU_ETI_SLGN_SoC(MueLu::ScalarDroppingDistanceLaplacian, S, LO, GO, N)

TPETRA_ETI_MANGLING_TYPEDEFS()

INSTANTIATE_TPETRA_MP_VECTOR_WRAPPER_NODES(MUELU_INST_S_LO_GO_N)

#endif

 //
