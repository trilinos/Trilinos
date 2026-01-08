// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSKERNELS_DEFAULT_TYPES_H
#define KOKKOSKERNELS_DEFAULT_TYPES_H

#include "Kokkos_Core.hpp"         //for LayoutLeft/LayoutRight
#include <KokkosKernels_config.h>  //for all the ETI #cmakedefine macros

// define a deprecated symbol = type in the global namespace
// and a non-deprecated version in Kokkos Kernels
// these deprecations were done in 4.4.
// Intel 19 doesn't seem to like deprecating a type alias
#if defined(KOKKOS_COMPILER_INTEL) && (KOKKOS_COMPILER_INTEL < 2000)
#define KK_IMPL_MAKE_TYPE_ALIAS(symbol, type) \
  using symbol = type;                        \
  namespace KokkosKernels {                   \
  using symbol = type;                        \
  }
#else
#define KK_IMPL_MAKE_TYPE_ALIAS(symbol, type)                            \
  using symbol [[deprecated("use KokkosKernels::" #symbol ".")]] = type; \
  namespace KokkosKernels {                                              \
  using symbol = type;                                                   \
  }
#endif

#if defined(KOKKOSKERNELS_INST_ORDINAL_INT)
KK_IMPL_MAKE_TYPE_ALIAS(default_lno_t, int)
#elif defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T)
KK_IMPL_MAKE_TYPE_ALIAS(default_lno_t, int64_t)
#else
// Non-ETI build: default to int
KK_IMPL_MAKE_TYPE_ALIAS(default_lno_t, int)
#endif
// Prefer int as the default offset type, because cuSPARSE doesn't support
// size_t for rowptrs.
#if defined(KOKKOSKERNELS_INST_OFFSET_INT)
KK_IMPL_MAKE_TYPE_ALIAS(default_size_type, int)
#elif defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)
KK_IMPL_MAKE_TYPE_ALIAS(default_size_type, size_t)
#else
// Non-ETI build: default to int
KK_IMPL_MAKE_TYPE_ALIAS(default_size_type, int)
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
KK_IMPL_MAKE_TYPE_ALIAS(default_layout, Kokkos::LayoutLeft)
#elif defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
KK_IMPL_MAKE_TYPE_ALIAS(default_layout, Kokkos::LayoutRight)
#else
KK_IMPL_MAKE_TYPE_ALIAS(default_layout, Kokkos::LayoutLeft)
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
KK_IMPL_MAKE_TYPE_ALIAS(default_scalar, double)
#elif defined(KOKKOSKERNELS_INST_FLOAT)
KK_IMPL_MAKE_TYPE_ALIAS(default_scalar, float)
#elif defined(KOKKOSKERNELS_INST_HALF)
KK_IMPL_MAKE_TYPE_ALIAS(default_scalar, Kokkos::Experimental::half_t)
#elif defined(KOKKOSKERNELS_INST_BHALF)
KK_IMPL_MAKE_TYPE_ALIAS(default_scalar, Kokkos::Experimental::bhalf_t)
#else
KK_IMPL_MAKE_TYPE_ALIAS(default_scalar, double)
#endif

#if defined(KOKKOS_ENABLE_CUDA)
KK_IMPL_MAKE_TYPE_ALIAS(default_device, Kokkos::Cuda)
#elif defined(KOKKOS_ENABLE_HIP)
KK_IMPL_MAKE_TYPE_ALIAS(default_device, Kokkos::HIP)
#elif defined(KOKKOS_ENABLE_OPENMPTARGET)
KK_IMPL_MAKE_TYPE_ALIAS(default_device, Kokkos::Experimental::OpenMPTarget)
#elif defined(KOKKOS_ENABLE_OPENMP)
KK_IMPL_MAKE_TYPE_ALIAS(default_device, Kokkos::OpenMP)
#elif defined(KOKKOS_ENABLE_THREADS)
KK_IMPL_MAKE_TYPE_ALIAS(default_device, Kokkos::Threads)
#else
KK_IMPL_MAKE_TYPE_ALIAS(default_device, Kokkos::Serial)
#endif

#undef KK_IMPL_MAKE_TYPE_ALIAS

#endif  // KOKKOSKERNELS_DEFAULT_TYPES_H
