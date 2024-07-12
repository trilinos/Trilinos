// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MueLu_Factory.hpp"

namespace MueLu {

bool Factory::timerSync_ = false;
#ifdef HAVE_MUELU_DEBUG
Factory::multipleCallCheckEnum Factory::multipleCallCheckGlobal_ = ENABLED;
#endif

}  // namespace MueLu
