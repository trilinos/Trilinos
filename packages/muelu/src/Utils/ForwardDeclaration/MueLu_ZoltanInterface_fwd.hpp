// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_ZOLTANINTERFACE_FWD_HPP
#define MUELU_ZOLTANINTERFACE_FWD_HPP

#include "MueLu_ConfigDefs.hpp"
#if defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MPI)

namespace MueLu {
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class ZoltanInterface;
}

#ifndef MUELU_ZOLTANINTERFACE_SHORT
#define MUELU_ZOLTANINTERFACE_SHORT
#endif

#endif

#endif  // MUELU_ZOLTANINTERFACE_FWD_HPP
