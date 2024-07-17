// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Tests
#include "TestSpMM.hpp"

// Devices
#include "Kokkos_Core.hpp"

#ifdef KOKKOS_ENABLE_THREADS
template void performance_test_driver< double, int, Kokkos::Threads>(
  const int nGrid, const int nIter,  const int ensemble_min,
  const int ensemble_max, const int ensemble_step);
#endif
