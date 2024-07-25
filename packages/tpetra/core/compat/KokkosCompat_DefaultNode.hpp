// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRAKOKKOSCOMPAT_DEFAULTNODE_HPP
#define TPETRAKOKKOSCOMPAT_DEFAULTNODE_HPP

#include "TpetraCore_config.h"

#if defined(TPETRA_ENABLE_DEPRECATED_CODE)
#warning "The header file Trilinos/packages/tpetra/core/compat/KokkosCompat_DefaultNode.hpp is deprecated. Use Tpetra_KokkosCompat_DefaultNode.hpp"
#include "Tpetra_KokkosCompat_DefaultNode.hpp"
#else
#error "The header file Trilinos/packages/tpetra/core/compat/KokkosCompat_DefaultNode.hpp is deprecated. Use Tpetra_KokkosCompat_DefaultNode.hpp"
#endif

namespace KokkosClassic {
    using DefaultNode [[deprecated]] = Tpetra::KokkosClassic::DefaultNode;
}

#endif // TPETRAKOKKOSCOMPAT_DEFAULTNODE_HPP

