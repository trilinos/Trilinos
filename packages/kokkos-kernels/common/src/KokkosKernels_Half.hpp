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

#ifndef KOKKOSKERNELS_HALF_HPP
#define KOKKOSKERNELS_HALF_HPP

#include "Kokkos_Core.hpp"

namespace KokkosKernels {
namespace Experimental {
////////////// BEGIN FP16/binary16 limits //////////////
#define KOKKOSKERNELS_IMPL_FP16_MAX 65504.0F  // Maximum normalized number
#define KOKKOSKERNELS_IMPL_FP16_MIN \
  0.000000059604645F  // Minimum normalized positive half precision number
#define KOKKOSKERNELS_IMPL_FP16_RADIX \
  2  // Value of the base of the exponent representation. TODO: on all archs?
#define KOKKOSKERNELS_IMPL_FP16_MANT_DIG \
  15  // Number of digits in the matissa that can be represented without losing
      // precision. TODO: Confirm this
#define KOKKOSKERNELS_IMPL_FP16_MIN_EXP \
  -14  // This is the smallest possible exponent value
#define KOKKOSKERNELS_IMPL_FP16_MAX_EXP \
  15  // This is the largest possible exponent value
#define KOKKOSKERNELS_IMPL_FP16_SIGNIFICAND_BITS 10
#define KOKKOSKERNELS_IMPL_FP16_EPSILON 0.0009765625F  // 1/2^10
#define KOKKOSKERNELS_IMPL_HUGE_VALH 0x7c00            // bits [10,14] set.
////////////// END FP16/binary16 limits //////////////

////////////// BEGIN BF16/float16 limits //////////////
#define KOKKOSKERNELS_IMPL_BF16_MAX 3.38953139e38  // Maximum normalized number
#define KOKKOSKERNELS_IMPL_BF16_MIN \
  1.1754494351e-38  // Minimum normalized positive bhalf number
#define KOKKOSKERNELS_IMPL_BF16_RADIX \
  2  // Value of the base of the exponent representation. TODO: on all archs?
#define KOKKOSKERNELS_IMPL_BF16_MANT_DIG_MIN 2
#define KOKKOSKERNELS_IMPL_BF16_MANT_DIG_MAX 3
#define KOKKOSKERNELS_IMPL_BF16_MANT_DIG \
  KOKKOSKERNELS_IMPL_BF16_MANT_DIG_MIN  // Number of digits in the matissa that
                                        // can be represented without losing
                                        // precision.
#define KOKKOSKERNELS_IMPL_BF16_MIN_EXP \
  -126  // This is the smallest possible exponent value
#define KOKKOSKERNELS_IMPL_BF16_MAX_EXP \
  127  // This is the largest possible exponent value
#define KOKKOSKERNELS_IMPL_BF16_EPSILON 0.0078125F  // 1/2^7
////////////// END BF16/bfloat16 limits //////////////

}  // namespace Experimental
}  // namespace KokkosKernels
#endif  // KOKKOSKERNELS_HALF_HPP
