/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

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
