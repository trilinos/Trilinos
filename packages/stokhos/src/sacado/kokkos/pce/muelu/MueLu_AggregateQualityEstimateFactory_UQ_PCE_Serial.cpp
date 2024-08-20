// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MueLu_ExplicitInstantiation.hpp"
#include "Stokhos_ConfigDefs.h"

#if defined(HAVE_STOKHOS_MUELU) && defined(HAVE_MUELU_EXPLICIT_INSTANTIATION) && defined(HAVE_STOKHOS_SACADO)

#include "Stokhos_Tpetra_ETI_Helpers_UQ_PCE.hpp"
#include "Stokhos_MueLu_UQ_PCE.hpp"

#include "MueLu_AggregateQualityEstimateFactory_def.hpp"

#define MUELU_INST_S_LO_GO_N(S, LO, GO, N) \
  template class MueLu::AggregateQualityEstimateFactory<S, LO, GO, N>;

#define MUELU_INST_N(N) \
  INSTANTIATE_TPETRA_UQ_PCE_N(MUELU_INST_S_LO_GO_N, N)

TPETRA_ETI_MANGLING_TYPEDEFS()

INSTANTIATE_TPETRA_UQ_PCE_SERIAL(MUELU_INST_N)

#endif


