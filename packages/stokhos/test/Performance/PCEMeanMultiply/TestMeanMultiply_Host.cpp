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

#ifdef KOKKOS_ENABLE_THREADS
using Kokkos::Threads;
INST_PERF_DRIVER(double, int, Threads)
#endif

#ifdef KOKKOS_ENABLE_OPENMP
using Kokkos::OpenMP;
INST_PERF_DRIVER(double, int, OpenMP)
#endif
