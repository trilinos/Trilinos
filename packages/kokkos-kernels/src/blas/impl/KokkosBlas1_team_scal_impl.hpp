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

#ifndef KOKKOSBLAS1_TEAM_SCAL_IMPL_HPP_
#define KOKKOSBLAS1_TEAM_SCAL_IMPL_HPP_

#include <Kokkos_Core.hpp>
#include "KokkosBlas1_serial_scal_impl.hpp"

namespace KokkosBlas {
namespace Impl {

///
/// Team Internal Impl
/// ====================
struct TeamScaleInternal {
  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member,
                                           const int m, const ScalarType alpha,
                                           /* */ ValueType *KOKKOS_RESTRICT A,
                                           const int as0) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, m),
                         [&](const int &i) { A[i * as0] *= alpha; });
    // member.team_barrier();
    return 0;
  }

  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member,
                                           const int m, const int n,
                                           const ScalarType alpha,
                                           /* */ ValueType *KOKKOS_RESTRICT A,
                                           const int as0, const int as1) {
    if (m > n) {
      Kokkos::parallel_for(
          Kokkos::TeamThreadRange(member, m), [&](const int &i) {
            SerialScaleInternal::invoke(n, alpha, A + i * as0, as1);
          });
    } else {
      Kokkos::parallel_for(
          Kokkos::TeamThreadRange(member, n), [&](const int &j) {
            SerialScaleInternal::invoke(m, alpha, A + j * as1, as0);
          });
    }
    // member.team_barrier();
    return 0;
  }
};

///
/// TeamVector Internal Impl
/// ========================
struct TeamVectorScaleInternal {
  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member,
                                           const int m, const ScalarType alpha,
                                           /* */ ValueType *KOKKOS_RESTRICT A,
                                           const int as0) {
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m),
                         [&](const int &i) { A[i * as0] *= alpha; });
    // member.team_barrier();
    return 0;
  }

  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member,
                                           const int m, const int n,
                                           const ScalarType alpha,
                                           /* */ ValueType *KOKKOS_RESTRICT A,
                                           const int as0, const int as1) {
    if (as0 > as1) {
      Kokkos::parallel_for(
          Kokkos::TeamThreadRange(member, m), [&](const int &i) {
            Kokkos::parallel_for(
                Kokkos::ThreadVectorRange(member, n),
                [&](const int &j) { A[i * as0 + j * as1] *= alpha; });
          });
    } else {
      Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(member, m), [&](const int &i) {
            Kokkos::parallel_for(
                Kokkos::TeamThreadRange(member, n),
                [&](const int &j) { A[i * as0 + j * as1] *= alpha; });
          });
    }
    // member.team_barrier();
    return 0;
  }
};

}  // namespace Impl
}  // namespace KokkosBlas

#endif
