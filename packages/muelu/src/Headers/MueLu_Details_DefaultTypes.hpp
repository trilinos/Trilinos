// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_USEDEFAULTTYPES_HPP
#define MUELU_USEDEFAULTTYPES_HPP

#include <Tpetra_KokkosCompat_DefaultNode.hpp>
#include "MueLu_config.hpp"

#include <Tpetra_Details_DefaultTypes.hpp>

namespace MueLu {

typedef Tpetra::Details::DefaultTypes::scalar_type DefaultScalar;

typedef int DefaultLocalOrdinal;

#if defined HAVE_MUELU_DEFAULT_GO_LONG
typedef long DefaultGlobalOrdinal;
#elif defined HAVE_MUELU_DEFAULT_GO_LONGLONG
typedef long long DefaultGlobalOrdinal;
#else
typedef int DefaultGlobalOrdinal;
#endif

typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType DefaultNode;
}  // namespace MueLu

#endif
