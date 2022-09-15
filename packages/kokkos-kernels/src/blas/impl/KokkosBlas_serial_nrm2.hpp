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

#ifndef KOKKOSBLAS_SERIAL_NRM2_HPP_
#define KOKKOSBLAS_SERIAL_NRM2_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>

namespace KokkosBlas {
namespace Impl {

///
/// Serial Internal Impl
/// ====================
template <typename ValueType>
KOKKOS_INLINE_FUNCTION static
    typename Kokkos::Details::InnerProductSpaceTraits<ValueType>::mag_type
    serial_nrm2(const int m, const ValueType *KOKKOS_RESTRICT X,
                const int xs0) {
  using IPT       = Kokkos::Details::InnerProductSpaceTraits<ValueType>;
  using norm_type = typename IPT::mag_type;

  norm_type nrm = Kokkos::ArithTraits<norm_type>::zero();

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int i = 0; i < m; ++i)
    nrm += IPT::norm(IPT::dot(X[i * xs0], X[i * xs0]));

  return Kokkos::ArithTraits<norm_type>::sqrt(nrm);
}

template <typename ValueType>
KOKKOS_INLINE_FUNCTION static void serial_nrm2(
    const int m, const int n, const ValueType *KOKKOS_RESTRICT X, const int xs0,
    const int xs1,
    typename Kokkos::Details::InnerProductSpaceTraits<ValueType>::mag_type
        *KOKKOS_RESTRICT R,
    const int ys0) {
  for (int vecIdx = 0; vecIdx < n; ++vecIdx)
    R[vecIdx * ys0] = serial_nrm2(m, X + vecIdx * xs1, xs0);

  return;
}

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSBLAS_SERIAL_NRM2_HPP_
