// @HEADER
// *****************************************************************************
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//
// Copyright 2012 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if 0

#include "KokkosCompat_View.hpp"
#include "KokkosCompat_View_def.hpp"

#include "Kokkos_Core.hpp"

namespace Kokkos {
  namespace Compat {

#if defined(KOKKOS_ENABLE_SERIAL)
#define COMPAT_INSTANT_SERIAL(T) \
    COMPAT_INSTANT(T,Kokkos::Serial)
#else
#  define COMPAT_INSTANT_SERIAL(T)
#endif // defined(KOKKOS_ENABLE_SERIAL)

#if defined(KOKKOS_ENABLE_THREADS)
#define COMPAT_INSTANT_THREADS(T) \
    COMPAT_INSTANT(T,Kokkos::Threads)
#else
#define COMPAT_INSTANT_THREADS(T)
#endif

#if defined(KOKKOS_ENABLE_OPENMP)
#define COMPAT_INSTANT_OPENMP(T) \
    COMPAT_INSTANT(T,Kokkos::OpenMP)
#else
#define COMPAT_INSTANT_OPENMP(T)
#endif

#define COMPAT_INSTANT_ALL(T) \
    COMPAT_INSTANT_SERIAL(T) \
    COMPAT_INSTANT_THREADS(T) \
    COMPAT_INSTANT_OPENMP(T) \
    COMPAT_INSTANT(T,Kokkos::HostSpace)

    COMPAT_INSTANT_ALL(float)
    COMPAT_INSTANT_ALL(double)
    COMPAT_INSTANT_ALL(int)
    COMPAT_INSTANT_ALL(long)
    COMPAT_INSTANT_ALL(unsigned)
    COMPAT_INSTANT_ALL(unsigned long)
    COMPAT_INSTANT_ALL(char)
    COMPAT_INSTANT_ALL(short)

#if defined(KOKKOS_ENABLE_OPENMP)
#define COMPAT_INSTANT_CUDA(T) \
    COMPAT_INSTANT(T,Kokkos::Cuda)
#else
    COMPAT_INSTANT_CUDA(T)
#endif

#if defined(KOKKOS_ENABLE_CUDA)
#define COMPAT_INSTANT_CUDA_UVM(T) \
    COMPAT_INSTANT(T,Kokkos::CudaUVMSpace)
#else
    COMPAT_INSTANT_CUDA_UVM(T)
#endif

    COMPAT_INSTANT_CUDA(float)
    COMPAT_INSTANT_CUDA(double)
    COMPAT_INSTANT_CUDA(int)
    COMPAT_INSTANT_CUDA(long)
    COMPAT_INSTANT_CUDA(unsigned)
    COMPAT_INSTANT_CUDA(unsigned long)
    COMPAT_INSTANT_CUDA(char)
    COMPAT_INSTANT_CUDA(short)

    COMPAT_INSTANT_CUDA_UVM(float)
    COMPAT_INSTANT_CUDA_UVM(double)
    COMPAT_INSTANT_CUDA_UVM(int)
    COMPAT_INSTANT_CUDA_UVM(long)
    COMPAT_INSTANT_CUDA_UVM(unsigned)
    COMPAT_INSTANT_CUDA_UVM(unsigned long)
    COMPAT_INSTANT_CUDA_UVM(char)
    COMPAT_INSTANT_CUDA_UVM(short)
  } // namespace Compat
} // namespace Kokkos

#endif // 0
