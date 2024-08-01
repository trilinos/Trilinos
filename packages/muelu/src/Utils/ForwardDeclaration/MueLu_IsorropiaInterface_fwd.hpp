// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_ISORROPIAINTERFACE_FWD_HPP
#define MUELU_ISORROPIAINTERFACE_FWD_HPP

#include "MueLu_ConfigDefs.hpp"
#if defined(HAVE_MUELU_ISORROPIA) && defined(HAVE_MPI)

namespace MueLu {
template <class LocalOrdinal, class GlobalOrdinal, class Node>
class IsorropiaInterface;
}

#ifndef MUELU_ISORROPIAINTERFACE_SHORT
#define MUELU_ISORROPIAINTERFACE_SHORT
#endif

#endif

#endif  // MUELU_ISORROPIAINTERFACE_FWD_HPP
