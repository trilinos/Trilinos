#// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_AMGX_SETUP_HPP
#define MUELU_AMGX_SETUP_HPP

#ifdef HAVE_MUELU_AMGX
#include <amgx_c.h>

namespace MueLu {

void MueLu_AMGX_initialize();
void MueLu_AMGX_initialize_plugins();

void MueLu_AMGX_finalize();
void MueLu_AMGX_finalize_plugins();
}  // namespace MueLu

#endif  // HAVE_MUELU_AMGX
#endif  // ifndef MUELU_AMGX_SETUP_DEF_HPP
