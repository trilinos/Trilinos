// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <MueLu_ConfigDefs.hpp>
#include "MueLu_AMGX_Setup.hpp"

namespace MueLu {

void MueLu_AMGX_initialize() {
  AMGX_SAFE_CALL(AMGX_initialize());
}

void MueLu_AMGX_initialize_plugins() {
  AMGX_SAFE_CALL(AMGX_initialize_plugins());
}

void MueLu_AMGX_finalize() {
  AMGX_SAFE_CALL(AMGX_finalize());
}

void MueLu_AMGX_finalize_plugins() {
  AMGX_print_summary();
  AMGX_SAFE_CALL(AMGX_finalize_plugins());
}
}  // namespace MueLu
