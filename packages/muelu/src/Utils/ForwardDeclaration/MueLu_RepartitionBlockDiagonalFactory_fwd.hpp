// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_REPARTITIONBLOCKDIAGONALFACTORY_FWD_HPP
#define MUELU_REPARTITIONBLOCKDIAGONALFACTORY_FWD_HPP

#include "MueLu_ConfigDefs.hpp"
#ifdef HAVE_MPI

namespace MueLu {
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class RepartitionBlockDiagonalFactory;
}

#ifndef MUELU_REPARTITIONBLOCKDIAGONALFACTORY_SHORT
#define MUELU_REPARTITIONBLOCKDIAGONALFACTORY_SHORT
#endif

#endif

#endif  // MUELU_REPARTITIONBLOCKDIAGONALFACTORY_FWD_HPP
