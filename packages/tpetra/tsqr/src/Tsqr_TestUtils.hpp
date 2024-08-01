// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TSQR_TESTUTILS_HPP
#define TSQR_TESTUTILS_HPP

/// \file Tsqr_TestUtils.hpp
/// \brief Utilities for testing various TSQR components.
/// \author Mark Hoemmen

#include "TpetraTSQR_config.h"

namespace Teuchos {
  // Forward declaration of Teuchos::Comm, so that we can use
  // RCP<Comm<int> > as the argument of methods defined in this header
  // file, without needing to include Teuchos_Comm.hpp.
  template<class Ordinal>
  class Comm;
}

#endif // TSQR_TESTUTILS_HPP
