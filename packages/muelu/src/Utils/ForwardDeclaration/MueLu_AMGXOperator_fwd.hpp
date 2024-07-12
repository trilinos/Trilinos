// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_AMGXOPERATOR_FWD_HPP
#define MUELU_AMGXOPERATOR_FWD_HPP

#include "MueLu_ConfigDefs.hpp"
#if defined(HAVE_MUELU_AMGX)

namespace MueLu {
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class AMGXOperator;
}

#ifndef MUELU_AMGXOPERATOR_SHORT
#define MUELU_AMGXOPERATOR_SHORT
#endif

#endif

#endif  // MUELU_AMGXOPERATOR_FWD_HPP
