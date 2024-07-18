// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Tests
#include "TestSpMv.hpp"

// Devices
#include "Kokkos_Core.hpp"

template <typename Storage>
void mainHost(int nGrid, int nIter, KokkosSparse::DeviceConfig dev_config) {
#ifdef __MIC__
  const int entry_min = 8;
  const int entry_max = 48;
  const int entry_step = 8;
#else
  const int entry_min = 4;
  const int entry_max = 32;
  const int entry_step = 4;
#endif

  performance_test_driver<Storage,entry_min,entry_max,entry_step>(
    nGrid,nIter,dev_config);
}

#ifdef KOKKOS_ENABLE_THREADS
template void mainHost< Stokhos::StaticFixedStorage<int,double,1,Kokkos::Threads> >(int nGrid, int nIter, KokkosSparse::DeviceConfig dev_config);
#endif

#ifdef KOKKOS_ENABLE_OPENMP
template void mainHost< Stokhos::StaticFixedStorage<int,double,1,Kokkos::OpenMP> >(int nGrid, int nIter, KokkosSparse::DeviceConfig dev_config);
#endif
