// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_SHIFTEDLAPLACIAN_FWD_HPP
#define MUELU_SHIFTEDLAPLACIAN_FWD_HPP

#include "MueLu_ConfigDefs.hpp"
#if defined(HAVE_MUELU_IFPACK2) && defined(HAVE_MUELU_BELOS) && defined(HAVE_MUELU_TPETRA_INST_INT_INT)

namespace MueLu {
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class ShiftedLaplacian;
}

#ifndef MUELU_SHIFTEDLAPLACIAN_SHORT
#define MUELU_SHIFTEDLAPLACIAN_SHORT
#endif

#endif

#endif  // MUELU_SHIFTEDLAPLACIAN_FWD_HPP
