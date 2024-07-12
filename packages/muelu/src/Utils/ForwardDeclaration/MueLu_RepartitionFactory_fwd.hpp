// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_REPARTITIONFACTORY_FWD_HPP
#define MUELU_REPARTITIONFACTORY_FWD_HPP

#include "MueLu_ConfigDefs.hpp"
#ifdef HAVE_MPI

namespace MueLu {
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class RepartitionFactory;
}

#ifndef MUELU_REPARTITIONFACTORY_SHORT
#define MUELU_REPARTITIONFACTORY_SHORT
#endif

#endif

#endif  // MUELU_REPARTITIONFACTORY_FWD_HPP
