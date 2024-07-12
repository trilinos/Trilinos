// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_ZOLTAN2INTERFACE_FWD_HPP
#define MUELU_ZOLTAN2INTERFACE_FWD_HPP

#include "MueLu_ConfigDefs.hpp"
#if defined(HAVE_MUELU_ZOLTAN2) && defined(HAVE_MPI)

namespace MueLu {
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class Zoltan2Interface;
}

#ifndef MUELU_ZOLTAN2INTERFACE_SHORT
#define MUELU_ZOLTAN2INTERFACE_SHORT
#endif

#endif

#endif  // MUELU_ZOLTAN2INTERFACE_FWD_HPP
