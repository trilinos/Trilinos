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

#ifndef KOKKOSBLAS2_TEAM_GEMV_HPP_
#define KOKKOSBLAS2_TEAM_GEMV_HPP_

#include <KokkosBlas2_team_gemv_spec.hpp>

namespace KokkosBlas {
namespace Experimental {

template <class AlgoTag, class TeamType, class MatrixType, class XVector, class YVector, class ScalarType>
void KOKKOS_INLINE_FUNCTION team_gemv(const TeamType& team, const char trans, const ScalarType& alpha,
                                      const MatrixType& A, const XVector& x, const ScalarType& beta, const YVector& y) {
  if (trans == 'N' || trans == 'n')
    TeamGemv<TeamType, Trans::NoTranspose, AlgoTag>::invoke(team, alpha, A, x, beta, y);
  else if (trans == 'T' || trans == 't')
    TeamGemv<TeamType, Trans::Transpose, AlgoTag>::invoke(team, alpha, A, x, beta, y);
  else if (trans == 'C' || trans == 'c')
    TeamGemv<TeamType, Trans::ConjTranspose, AlgoTag>::invoke(team, alpha, A, x, beta, y);
  else {
    Kokkos::abort("Matrix mode not supported");
  }
}

// default AlgoTag
template <class TeamType, class MatrixType, class XVector, class YVector, class ScalarType>
void KOKKOS_INLINE_FUNCTION team_gemv(const TeamType& team, const char trans, const ScalarType& alpha,
                                      const MatrixType& A, const XVector& x, const ScalarType& beta, const YVector& y) {
  team_gemv<KokkosBlas::Algo::Gemv::Default>(team, trans, alpha, A, x, beta, y);
}

template <class AlgoTag, class TeamType, class MatrixType, class XVector, class YVector, class ScalarType>
void KOKKOS_INLINE_FUNCTION teamvector_gemv(const TeamType& team, const char trans, const ScalarType& alpha,
                                            const MatrixType& A, const XVector& x, const ScalarType& beta,
                                            const YVector& y) {
  if (trans == 'N' || trans == 'n') {
    KokkosBlas::TeamVectorGemv<TeamType, Trans::NoTranspose, AlgoTag>::invoke(team, alpha, A, x, beta, y);
  } else if (trans == 'T' || trans == 't') {
    KokkosBlas::TeamVectorGemv<TeamType, Trans::Transpose, AlgoTag>::invoke(team, alpha, A, x, beta, y);
  } else if (trans == 'C' || trans == 'c') {
    KokkosBlas::TeamVectorGemv<TeamType, Trans::ConjTranspose, AlgoTag>::invoke(team, alpha, A, x, beta, y);
  } else {
    Kokkos::abort("Matrix mode not supported");
  }
}

// default AlgoTag
template <class TeamType, class MatrixType, class XVector, class YVector, class ScalarType>
void KOKKOS_INLINE_FUNCTION team_vector_gemv(const TeamType& team, const char trans, const ScalarType& alpha,
                                             const MatrixType& A, const XVector& x, const ScalarType& beta,
                                             const YVector& y) {
  teamvector_gemv<KokkosBlas::Algo::Gemv::Default>(team, trans, alpha, A, x, beta, y);
}

}  // namespace Experimental
}  // namespace KokkosBlas

#endif
