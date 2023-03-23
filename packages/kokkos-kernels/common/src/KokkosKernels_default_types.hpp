//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOSKERNELS_DEFAULT_TYPES_H
#define KOKKOSKERNELS_DEFAULT_TYPES_H

#include "Kokkos_Core.hpp"         //for LayoutLeft/LayoutRight
#include <KokkosKernels_config.h>  //for all the ETI #cmakedefine macros

#if defined(KOKKOSKERNELS_INST_ORDINAL_INT)
using default_lno_t = int;
#elif defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T)
using default_lno_t     = int64_t;
#else
using default_lno_t     = int;
#endif
// Prefer int as the default offset type, because cuSPARSE doesn't support
// size_t for rowptrs.
#if defined(KOKKOSKERNELS_INST_OFFSET_INT)
using default_size_type = int;
#elif defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)
using default_size_type = size_t;
#else
using default_size_type = int;
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
using default_layout = Kokkos::LayoutLeft;
#elif defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
using default_layout    = Kokkos::LayoutRight;
#else
using default_layout    = Kokkos::LayoutLeft;
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
using default_scalar = double;
#elif defined(KOKKOSKERNELS_INST_FLOAT)
using default_scalar    = float;
#elif defined(KOKKOSKERNELS_INST_HALF)
using default_scalar    = Kokkos::Experimental::half_t;
#elif defined(KOKKOSKERNELS_INST_BHALF)
using default_scalar = Kokkos::Experimental::bhalf_t;
#else
using default_scalar = double;
#endif

#if defined(KOKKOS_ENABLE_CUDA)
using default_device = Kokkos::Cuda;
#elif defined(KOKKOS_ENABLE_HIP)
using default_device    = Kokkos::Experimental::HIP;
#elif defined(KOKKOS_ENABLE_OPENMPTARGET)
using default_device    = Kokkos::Experimental::OpenMPTarget;
#elif defined(KOKKOS_ENABLE_OPENMP)
using default_device = Kokkos::OpenMP;
#elif defined(KOKKOS_ENABLE_PTHREAD) || defined(KOKKOS_ENABLE_THREADS)
using default_device = Kokkos::Threads;
#else
using default_device = Kokkos::Serial;
#endif

#endif  // KOKKOSKERNELS_DEFAULT_TYPES_H
