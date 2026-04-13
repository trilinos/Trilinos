// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_IFPACK2SMOOTHER_FWD_HPP
#define MUELU_IFPACK2SMOOTHER_FWD_HPP

#include "MueLu_ConfigDefs.hpp"
#if defined(HAVE_MUELU_IFPACK2)

namespace MueLu {
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class Ifpack2Smoother;
}

#ifndef MUELU_IFPACK2SMOOTHER_SHORT
#define MUELU_IFPACK2SMOOTHER_SHORT
#endif

#endif

#endif  // MUELU_IFPACK2SMOOTHER_FWD_HPP
