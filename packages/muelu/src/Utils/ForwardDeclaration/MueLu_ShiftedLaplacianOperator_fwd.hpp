// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_SHIFTEDLAPLACIANOPERATOR_FWD_HPP
#define MUELU_SHIFTEDLAPLACIANOPERATOR_FWD_HPP

#include "MueLu_ConfigDefs.hpp"
#if defined(HAVE_MUELU_BELOS) && defined(HAVE_MUELU_TPETRA_INST_INT_INT)

namespace MueLu {
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class ShiftedLaplacianOperator;
}

#ifndef MUELU_SHIFTEDLAPLACIANOPERATOR_SHORT
#define MUELU_SHIFTEDLAPLACIANOPERATOR_SHORT
#endif

#endif

#endif  // MUELU_SHIFTEDLAPLACIANOPERATOR_FWD_HPP
