// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_STRATIMIKOSSMOOTHER_FWD_HPP
#define MUELU_STRATIMIKOSSMOOTHER_FWD_HPP

#include "MueLu_ConfigDefs.hpp"
#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)

namespace MueLu {
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class StratimikosSmoother;
}

#ifndef MUELU_STRATIMIKOSSMOOTHER_SHORT
#define MUELU_STRATIMIKOSSMOOTHER_SHORT
#endif

#endif

#endif  // MUELU_STRATIMIKOSSMOOTHER_FWD_HPP
