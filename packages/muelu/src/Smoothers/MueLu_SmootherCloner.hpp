// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_SMOOTHERBASECLONER_HPP
#define MUELU_SMOOTHERBASECLONER_HPP

#include <Xpetra_MultiVector_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"

#include "MueLu_SmootherBase.hpp"
#if defined(HAVE_MUELU_IFPACK2)
#include "MueLu_Ifpack2Smoother.hpp"
#endif
#include "MueLu_TrilinosSmoother.hpp"

namespace MueLu {

}  // namespace MueLu

#define MUELU_SMOOTHERBASECLONER_SHORT
#endif  // ifndef MUELU_SMOOTHERBASECLONER_HPP
