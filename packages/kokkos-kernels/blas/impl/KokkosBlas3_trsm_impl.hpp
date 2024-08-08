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

#ifndef KOKKOSBLAS3_TRSM_IMPL_HPP_
#define KOKKOSBLAS3_TRSM_IMPL_HPP_

/// \file KokkosBlas3_trsm_impl.hpp
/// \brief Implementation(s) of triangular linear system solve (with multiple
/// RHSs) \brief Sequential fall-back implementation calls the exisiting serial
/// batched TRSM. \brief Two sequential fall-back implementations for conjugate
/// transpose case are \brief also based on the exisiting serial batched TRSM.

#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "KokkosBlas1_set_impl.hpp"
#include "KokkosBatched_Trsm_Decl.hpp"
#include "KokkosBatched_Trsm_Serial_Impl.hpp"

namespace KokkosBlas {
namespace Impl {

template <typename ScalarType, typename ValueType>
int SerialTrsmInternalLeftLowerConj(const bool use_unit_diag, const int m, const int n, const ScalarType alpha,
                                    const ValueType* KOKKOS_RESTRICT A, const int as0, const int as1,
                                    /**/ ValueType* KOKKOS_RESTRICT B, const int bs0, const int bs1) {
  typedef Kokkos::ArithTraits<ValueType> AT;

  const ScalarType one(1.0), zero(0.0);

  if (alpha == zero)
    SerialSetInternal::invoke(m, n, zero, B, bs0, bs1);
  else {
    if (alpha != one) KokkosBlas::Impl::SerialScaleInternal::invoke(m, n, alpha, B, bs0, bs1);
    if (m <= 0 || n <= 0) return 0;

    for (int p = 0; p < m; ++p) {
      const int iend = m - p - 1, jend = n;

      const ValueType* KOKKOS_RESTRICT a21 = A + (p + 1) * as0 + p * as1;

      ValueType *KOKKOS_RESTRICT b1t = B + p * bs0, *KOKKOS_RESTRICT B2 = B + (p + 1) * bs0;

      if (!use_unit_diag) {
        const ValueType alpha11 = AT::conj(A[p * as0 + p * as1]);
        for (int j = 0; j < jend; ++j) b1t[j * bs1] = b1t[j * bs1] / alpha11;
      }

      for (int i = 0; i < iend; ++i)
        for (int j = 0; j < jend; ++j) B2[i * bs0 + j * bs1] -= AT::conj(a21[i * as0]) * b1t[j * bs1];
    }
  }
  return 0;
}

template <typename ScalarType, typename ValueType>
int SerialTrsmInternalLeftUpperConj(const bool use_unit_diag, const int m, const int n, const ScalarType alpha,
                                    const ValueType* KOKKOS_RESTRICT A, const int as0, const int as1,
                                    /**/ ValueType* KOKKOS_RESTRICT B, const int bs0, const int bs1) {
  typedef Kokkos::ArithTraits<ValueType> AT;

  const ScalarType one(1.0), zero(0.0);

  if (alpha == zero)
    SerialSetInternal::invoke(m, n, zero, B, bs0, bs1);
  else {
    if (alpha != one) KokkosBlas::Impl::SerialScaleInternal::invoke(m, n, alpha, B, bs0, bs1);
    if (m <= 0 || n <= 0) return 0;

    ValueType* KOKKOS_RESTRICT B0 = B;
    for (int p = (m - 1); p >= 0; --p) {
      const int iend = p, jend = n;

      const ValueType* KOKKOS_RESTRICT a01 = A + p * as1;
      ValueType* KOKKOS_RESTRICT b1t       = B + p * bs0;

      if (!use_unit_diag) {
        const ValueType alpha11 = AT::conj(A[p * as0 + p * as1]);
        for (int j = 0; j < n; ++j) b1t[j * bs1] = b1t[j * bs1] / alpha11;
      }

      if (p > 0) {  // Note: A workaround to produce correct results for
                    // complex<double> with Intel-18.2.199
        for (int i = 0; i < iend; ++i)
          for (int j = 0; j < jend; ++j) B0[i * bs0 + j * bs1] -= AT::conj(a01[i * as0]) * b1t[j * bs1];
      }
    }
  }
  return 0;
}

template <class AViewType, class BViewType>
void SerialTrsm_Invoke(const char side[], const char uplo[], const char trans[], const char diag[],
                       typename BViewType::const_value_type& alpha, const AViewType& A, const BViewType& B) {
  using KokkosBatched::Algo;
  using KokkosBatched::Diag;

  // Side::Left, Uplo::Lower, Trans::NoTranspose
  if (((side[0] == 'L') || (side[0] == 'l')) && ((uplo[0] == 'L') || (uplo[0] == 'l')) &&
      ((trans[0] == 'N') || (trans[0] == 'n')) && ((diag[0] == 'U') || (diag[0] == 'u')))
    KokkosBatched::SerialTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, B.extent(0), B.extent(1), alpha, A.data(), A.stride(0), A.stride(1), B.data(),
        B.stride(0), B.stride(1));
  if (((side[0] == 'L') || (side[0] == 'l')) && ((uplo[0] == 'L') || (uplo[0] == 'l')) &&
      ((trans[0] == 'N') || (trans[0] == 'n')) && ((diag[0] == 'N') || (diag[0] == 'n')))
    KokkosBatched::SerialTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(
        Diag::NonUnit::use_unit_diag, B.extent(0), B.extent(1), alpha, A.data(), A.stride(0), A.stride(1), B.data(),
        B.stride(0), B.stride(1));

  // Side::Left, Uplo::Lower, Trans::Transpose
  if (((side[0] == 'L') || (side[0] == 'l')) && ((uplo[0] == 'L') || (uplo[0] == 'l')) &&
      ((trans[0] == 'T') || (trans[0] == 't')) && ((diag[0] == 'U') || (diag[0] == 'u')))
    KokkosBatched::SerialTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, B.extent(0), B.extent(1), alpha, A.data(), A.stride(1), A.stride(0), B.data(),
        B.stride(0), B.stride(1));
  if (((side[0] == 'L') || (side[0] == 'l')) && ((uplo[0] == 'L') || (uplo[0] == 'l')) &&
      ((trans[0] == 'T') || (trans[0] == 't')) && ((diag[0] == 'N') || (diag[0] == 'n')))
    KokkosBatched::SerialTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(
        Diag::NonUnit::use_unit_diag, B.extent(0), B.extent(1), alpha, A.data(), A.stride(1), A.stride(0), B.data(),
        B.stride(0), B.stride(1));

  // Side::Left, Uplo::Lower, Trans::ConjTranspose
  if (((side[0] == 'L') || (side[0] == 'l')) && ((uplo[0] == 'L') || (uplo[0] == 'l')) &&
      ((trans[0] == 'C') || (trans[0] == 'c')) && ((diag[0] == 'U') || (diag[0] == 'u')))
    SerialTrsmInternalLeftUpperConj(Diag::Unit::use_unit_diag, B.extent(0), B.extent(1), alpha, A.data(), A.stride(1),
                                    A.stride(0), B.data(), B.stride(0), B.stride(1));
  if (((side[0] == 'L') || (side[0] == 'l')) && ((uplo[0] == 'L') || (uplo[0] == 'l')) &&
      ((trans[0] == 'C') || (trans[0] == 'c')) && ((diag[0] == 'N') || (diag[0] == 'n')))
    SerialTrsmInternalLeftUpperConj(Diag::NonUnit::use_unit_diag, B.extent(0), B.extent(1), alpha, A.data(),
                                    A.stride(1), A.stride(0), B.data(), B.stride(0), B.stride(1));

  // Side::Left, Uplo::Upper, Trans::NoTranspose
  if (((side[0] == 'L') || (side[0] == 'l')) && ((uplo[0] == 'U') || (uplo[0] == 'u')) &&
      ((trans[0] == 'N') || (trans[0] == 'n')) && ((diag[0] == 'U') || (diag[0] == 'u')))
    KokkosBatched::SerialTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, B.extent(0), B.extent(1), alpha, A.data(), A.stride(0), A.stride(1), B.data(),
        B.stride(0), B.stride(1));
  if (((side[0] == 'L') || (side[0] == 'l')) && ((uplo[0] == 'U') || (uplo[0] == 'u')) &&
      ((trans[0] == 'N') || (trans[0] == 'n')) && ((diag[0] == 'N') || (diag[0] == 'n')))
    KokkosBatched::SerialTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(
        Diag::NonUnit::use_unit_diag, B.extent(0), B.extent(1), alpha, A.data(), A.stride(0), A.stride(1), B.data(),
        B.stride(0), B.stride(1));

  // Side::Left, Uplo::Upper, Trans::Transpose
  if (((side[0] == 'L') || (side[0] == 'l')) && ((uplo[0] == 'U') || (uplo[0] == 'u')) &&
      ((trans[0] == 'T') || (trans[0] == 't')) && ((diag[0] == 'U') || (diag[0] == 'u')))
    KokkosBatched::SerialTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, B.extent(0), B.extent(1), alpha, A.data(), A.stride(1), A.stride(0), B.data(),
        B.stride(0), B.stride(1));
  if (((side[0] == 'L') || (side[0] == 'l')) && ((uplo[0] == 'U') || (uplo[0] == 'u')) &&
      ((trans[0] == 'T') || (trans[0] == 't')) && ((diag[0] == 'N') || (diag[0] == 'n')))
    KokkosBatched::SerialTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(
        Diag::NonUnit::use_unit_diag, B.extent(0), B.extent(1), alpha, A.data(), A.stride(1), A.stride(0), B.data(),
        B.stride(0), B.stride(1));

  // Side::Left, Uplo::Upper, Trans::ConjTranspose
  if (((side[0] == 'L') || (side[0] == 'l')) && ((uplo[0] == 'U') || (uplo[0] == 'u')) &&
      ((trans[0] == 'C') || (trans[0] == 'c')) && ((diag[0] == 'U') || (diag[0] == 'u')))
    SerialTrsmInternalLeftLowerConj(Diag::Unit::use_unit_diag, B.extent(0), B.extent(1), alpha, A.data(), A.stride(1),
                                    A.stride(0), B.data(), B.stride(0), B.stride(1));
  if (((side[0] == 'L') || (side[0] == 'l')) && ((uplo[0] == 'U') || (uplo[0] == 'u')) &&
      ((trans[0] == 'C') || (trans[0] == 'c')) && ((diag[0] == 'N') || (diag[0] == 'n')))
    SerialTrsmInternalLeftLowerConj(Diag::NonUnit::use_unit_diag, B.extent(0), B.extent(1), alpha, A.data(),
                                    A.stride(1), A.stride(0), B.data(), B.stride(0), B.stride(1));
  ////
  // Side::Right, Uplo::Lower, Trans::NoTranspose
  if (((side[0] == 'R') || (side[0] == 'r')) && ((uplo[0] == 'L') || (uplo[0] == 'l')) &&
      ((trans[0] == 'N') || (trans[0] == 'n')) && ((diag[0] == 'U') || (diag[0] == 'u')))
    KokkosBatched::SerialTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, B.extent(1), B.extent(0), alpha, A.data(), A.stride(1), A.stride(0), B.data(),
        B.stride(1), B.stride(0));
  if (((side[0] == 'R') || (side[0] == 'r')) && ((uplo[0] == 'L') || (uplo[0] == 'l')) &&
      ((trans[0] == 'N') || (trans[0] == 'n')) && ((diag[0] == 'N') || (diag[0] == 'n')))
    KokkosBatched::SerialTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(
        Diag::NonUnit::use_unit_diag, B.extent(1), B.extent(0), alpha, A.data(), A.stride(1), A.stride(0), B.data(),
        B.stride(1), B.stride(0));

  // Side::Right, Uplo::Lower, Trans::Transpose
  if (((side[0] == 'R') || (side[0] == 'r')) && ((uplo[0] == 'L') || (uplo[0] == 'l')) &&
      ((trans[0] == 'T') || (trans[0] == 't')) && ((diag[0] == 'U') || (diag[0] == 'u')))
    KokkosBatched::SerialTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, B.extent(1), B.extent(0), alpha, A.data(), A.stride(0), A.stride(1), B.data(),
        B.stride(1), B.stride(0));
  if (((side[0] == 'R') || (side[0] == 'r')) && ((uplo[0] == 'L') || (uplo[0] == 'l')) &&
      ((trans[0] == 'T') || (trans[0] == 't')) && ((diag[0] == 'N') || (diag[0] == 'n')))
    KokkosBatched::SerialTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(
        Diag::NonUnit::use_unit_diag, B.extent(1), B.extent(0), alpha, A.data(), A.stride(0), A.stride(1), B.data(),
        B.stride(1), B.stride(0));

  // Side::Right, Uplo::Lower, Trans::ConjTranspose
  if (((side[0] == 'R') || (side[0] == 'r')) && ((uplo[0] == 'L') || (uplo[0] == 'l')) &&
      ((trans[0] == 'C') || (trans[0] == 'c')) && ((diag[0] == 'U') || (diag[0] == 'u')))
    SerialTrsmInternalLeftLowerConj(Diag::Unit::use_unit_diag, B.extent(1), B.extent(0), alpha, A.data(), A.stride(0),
                                    A.stride(1), B.data(), B.stride(1), B.stride(0));
  if (((side[0] == 'R') || (side[0] == 'r')) && ((uplo[0] == 'L') || (uplo[0] == 'l')) &&
      ((trans[0] == 'C') || (trans[0] == 'c')) && ((diag[0] == 'N') || (diag[0] == 'n')))
    SerialTrsmInternalLeftLowerConj(Diag::NonUnit::use_unit_diag, B.extent(1), B.extent(0), alpha, A.data(),
                                    A.stride(0), A.stride(1), B.data(), B.stride(1), B.stride(0));

  // Side::Right, Uplo::Upper, Trans::NoTranspose
  if (((side[0] == 'R') || (side[0] == 'r')) && ((uplo[0] == 'U') || (uplo[0] == 'u')) &&
      ((trans[0] == 'N') || (trans[0] == 'n')) && ((diag[0] == 'U') || (diag[0] == 'u')))
    KokkosBatched::SerialTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, B.extent(1), B.extent(0), alpha, A.data(), A.stride(1), A.stride(0), B.data(),
        B.stride(1), B.stride(0));
  if (((side[0] == 'R') || (side[0] == 'r')) && ((uplo[0] == 'U') || (uplo[0] == 'u')) &&
      ((trans[0] == 'N') || (trans[0] == 'n')) && ((diag[0] == 'N') || (diag[0] == 'n')))
    KokkosBatched::SerialTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(
        Diag::NonUnit::use_unit_diag, B.extent(1), B.extent(0), alpha, A.data(), A.stride(1), A.stride(0), B.data(),
        B.stride(1), B.stride(0));

  // Side::Right, Uplo::Upper, Trans::Transpose
  if (((side[0] == 'R') || (side[0] == 'r')) && ((uplo[0] == 'U') || (uplo[0] == 'u')) &&
      ((trans[0] == 'T') || (trans[0] == 't')) && ((diag[0] == 'U') || (diag[0] == 'u')))
    KokkosBatched::SerialTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(
        Diag::Unit::use_unit_diag, B.extent(1), B.extent(0), alpha, A.data(), A.stride(0), A.stride(1), B.data(),
        B.stride(1), B.stride(0));
  if (((side[0] == 'R') || (side[0] == 'r')) && ((uplo[0] == 'U') || (uplo[0] == 'u')) &&
      ((trans[0] == 'T') || (trans[0] == 't')) && ((diag[0] == 'N') || (diag[0] == 'n')))
    KokkosBatched::SerialTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(
        Diag::NonUnit::use_unit_diag, B.extent(1), B.extent(0), alpha, A.data(), A.stride(0), A.stride(1), B.data(),
        B.stride(1), B.stride(0));

  // Side::Right, Uplo::Upper, Trans::ConjTranspose
  if (((side[0] == 'R') || (side[0] == 'r')) && ((uplo[0] == 'U') || (uplo[0] == 'u')) &&
      ((trans[0] == 'C') || (trans[0] == 'c')) && ((diag[0] == 'U') || (diag[0] == 'u')))
    SerialTrsmInternalLeftUpperConj(Diag::Unit::use_unit_diag, B.extent(1), B.extent(0), alpha, A.data(), A.stride(0),
                                    A.stride(1), B.data(), B.stride(1), B.stride(0));
  if (((side[0] == 'R') || (side[0] == 'r')) && ((uplo[0] == 'U') || (uplo[0] == 'u')) &&
      ((trans[0] == 'C') || (trans[0] == 'c')) && ((diag[0] == 'N') || (diag[0] == 'n')))
    SerialTrsmInternalLeftUpperConj(Diag::NonUnit::use_unit_diag, B.extent(1), B.extent(0), alpha, A.data(),
                                    A.stride(0), A.stride(1), B.data(), B.stride(1), B.stride(0));
}

}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSBLAS3_TRSM_IMPL_HPP_
