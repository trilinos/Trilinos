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
#include "KokkosBatched_Util.hpp"

/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace KokkosBatched {

void print_compiler_info() {
  printf("  supported pragmas:\n");
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
  printf("    #pragma unroll\n");
#endif
#if defined(KOKKOS_ENABLE_PRAGMA_IVDEP)
  printf("    #pragma ivdep\n");
#endif
#if defined KOKKOS_ENABLE_PRAGMA_VECTOR
  printf("    #pragma vector always\n");
#endif
#if !defined(KOKKOS_DEBUG)
#if defined(KOKKOS_ENABLE_PRAGMA_SIMD)
  printf("    #pragma simd\n");
#endif
#endif
  printf("\n");

  printf("  supported avx intrinsics:\n");
#if defined(__AVX__) || defined(__AVX2__)
  printf("    __m256d : ");
#if defined(__AVX__)
  printf("AVX ");
#endif
#if defined(__AVX2__)
  printf("AVX2 ");
#endif
  printf("\n");
#endif
#if defined(__AVX512F__)
  printf("    __m512d : AVX512F");
#endif
#if defined(__FMA__)
  printf("    FMA is supported\n");
#else
  printf("    FMA is not supported\n");
#endif
}
}  // namespace KokkosBatched
