// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_AMESOS2SMOOTHER_FWD_HPP
#define MUELU_AMESOS2SMOOTHER_FWD_HPP

#include "MueLu_ConfigDefs.hpp"
#if defined(HAVE_MUELU_AMESOS2)

namespace MueLu {
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class Amesos2Smoother;
}

#ifndef MUELU_AMESOS2SMOOTHER_SHORT
#define MUELU_AMESOS2SMOOTHER_SHORT
#endif

#endif

#endif  // MUELU_AMESOS2SMOOTHER_FWD_HPP
