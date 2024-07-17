// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MueLu_ExplicitInstantiation.hpp"

#include "MueLu_TwoLevelMatlabFactory_def.hpp"

#include "TpetraCore_ETIHelperMacros.h"

#ifdef HAVE_MUELU_MATLAB
#define MUELU_LOCAL_INSTANT(S, LO, GO, N) \
  template class MueLu::TwoLevelMatlabFactory<S, LO, GO, N>;

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(MUELU_LOCAL_INSTANT)

#endif
