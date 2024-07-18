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
#include "Stokhos_Cuda_CrsMatrix.hpp"
#include "Kokkos_CrsMatrix_MP_Vector_Cuda.hpp"

// Devices
#include "Kokkos_Core.hpp"

template <typename Storage>
void mainCuda(int nGrid, int nIter, KokkosSparse::DeviceConfig dev_config) {
  const int entry_min = 16;
  const int entry_max = 64;
  const int entry_step = 16;
  performance_test_driver<Storage,entry_min,entry_max,entry_step>(
    nGrid,nIter,dev_config);
}

template void mainCuda< Stokhos::StaticFixedStorage<int,double,1,Kokkos::Cuda> >(int nGrid, int nIter, KokkosSparse::DeviceConfig dev_config);
