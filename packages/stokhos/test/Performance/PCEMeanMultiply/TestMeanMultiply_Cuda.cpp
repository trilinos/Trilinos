// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Tests
#include "TestMeanMultiply.hpp"

// Devices
#include "Kokkos_Core.hpp"

#ifdef KOKKOS_ENABLE_CUDA
using Kokkos::Cuda;
INST_PERF_DRIVER(double, int, Cuda)
#endif
