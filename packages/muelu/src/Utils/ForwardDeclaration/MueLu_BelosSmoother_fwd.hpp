// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_BELOSSMOOTHER_FWD_HPP
#define MUELU_BELOSSMOOTHER_FWD_HPP

#include "MueLu_ConfigDefs.hpp"
#if defined(HAVE_MUELU_BELOS)

namespace MueLu {
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class BelosSmoother;
}

#ifndef MUELU_BELOSSMOOTHER_SHORT
#define MUELU_BELOSSMOOTHER_SHORT
#endif

#endif

#endif  // MUELU_BELOSSMOOTHER_FWD_HPP
