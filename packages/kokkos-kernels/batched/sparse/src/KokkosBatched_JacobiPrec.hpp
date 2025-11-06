// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_JACOBIPREC_HPP
#define KOKKOSBATCHED_JACOBIPREC_HPP

/// \author Kim Liegeois (knliege@sandia.gov)

#include "KokkosBatched_HadamardProduct.hpp"

namespace KokkosBatched {

/// \brief Batched Jacobi Preconditioner:
///
/// \tparam ValuesViewType: Input type for the values of the diagonal

template <class ValuesViewType>
class JacobiPrec {
 public:
  using ScalarType    = typename ValuesViewType::non_const_value_type;
  using MagnitudeType = typename KokkosKernels::ArithTraits<ScalarType>::mag_type;

 private:
  ValuesViewType diag_values;
  int n_operators;
  int n_rows;
  int n_colums;
  mutable bool computed_inverse = false;

 public:
  KOKKOS_INLINE_FUNCTION
  JacobiPrec(const ValuesViewType &_diag_values) : diag_values(_diag_values) {
    n_operators = _diag_values.extent(0);
    n_rows      = _diag_values.extent(1);
    n_colums    = n_rows;
  }

  KOKKOS_INLINE_FUNCTION
  ~JacobiPrec() {}

  KOKKOS_INLINE_FUNCTION void setComputedInverse() { computed_inverse = true; }

  template <typename MemberType, typename ArgMode>
  KOKKOS_INLINE_FUNCTION void computeInverse(const MemberType &member) const {
    auto one     = KokkosKernels::ArithTraits<MagnitudeType>::one();
    auto epsilon = KokkosKernels::ArithTraits<MagnitudeType>::epsilon();
    int tooSmall = 0;
    if (std::is_same<ArgMode, Mode::Serial>::value) {
      for (int i = 0; i < n_operators; ++i)
        for (int j = 0; j < n_colums; ++j) {
          if (Kokkos::abs<ScalarType>(diag_values(i, j)) <= epsilon) {
            ++tooSmall;
            diag_values(i, j) = one;
          } else
            diag_values(i, j) = one / diag_values(i, j);
        }
    } else if (std::is_same<ArgMode, Mode::Team>::value) {
      auto diag_values_array = diag_values.data();
      auto vs0               = diag_values.stride(0);
      auto vs1               = diag_values.stride(1);

      Kokkos::parallel_reduce(
          Kokkos::TeamThreadRange(member, 0, n_operators * n_rows),
          [&](const int &iTemp, int &ltooSmall) {
            int i, j;
            getIndices<int, typename ValuesViewType::array_layout>(iTemp, n_rows, n_operators, j, i);
            if (Kokkos::abs<ScalarType>(diag_values_array[i * vs0 + j * vs1]) <= epsilon) {
              ltooSmall++;
              diag_values_array[i * vs0 + j * vs1] = one;
            } else
              diag_values_array[i * vs0 + j * vs1] = one / diag_values_array[i * vs0 + j * vs1];
          },
          tooSmall);
    } else if (std::is_same<ArgMode, Mode::TeamVector>::value) {
      auto diag_values_array = diag_values.data();
      auto vs0               = diag_values.stride(0);
      auto vs1               = diag_values.stride(1);

      Kokkos::parallel_reduce(
          Kokkos::TeamVectorRange(member, 0, n_operators * n_rows),
          [&](const int &iTemp, int &ltooSmall) {
            int i, j;
            getIndices<int, typename ValuesViewType::array_layout>(iTemp, n_rows, n_operators, j, i);
            if (Kokkos::abs<ScalarType>(diag_values_array[i * vs0 + j * vs1]) <= epsilon) {
              ltooSmall++;
              diag_values_array[i * vs0 + j * vs1] = one;
            } else
              diag_values_array[i * vs0 + j * vs1] = one / diag_values_array[i * vs0 + j * vs1];
          },
          tooSmall);
    }

    if (tooSmall > 0)
      Kokkos::printf(
          "KokkosBatched::JacobiPrec: %d entrie(s) has/have a too small "
          "magnitude and have been replaced by one, \n",
          (int)tooSmall);
    computed_inverse = true;
  }

  KOKKOS_INLINE_FUNCTION void computeInverse() const {
    auto one     = KokkosKernels::ArithTraits<MagnitudeType>::one();
    auto epsilon = KokkosKernels::ArithTraits<MagnitudeType>::epsilon();
    int tooSmall = 0;

    for (int i = 0; i < n_operators; ++i)
      for (int j = 0; j < n_colums; ++j) {
        if (Kokkos::abs<ScalarType>(diag_values(i, j)) <= epsilon) {
          ++tooSmall;
          diag_values(i, j) = one;
        } else
          diag_values(i, j) = one / diag_values(i, j);
      }

    if (tooSmall > 0)
      Kokkos::printf(
          "KokkosBatched::JacobiPrec: %d entrie(s) has/have a too small "
          "magnitude and have been replaced by one, \n",
          (int)tooSmall);
    computed_inverse = true;
  }

  template <typename ArgTrans, typename ArgMode, int sameXY, typename MemberType, typename XViewType,
            typename YViewType>
  KOKKOS_INLINE_FUNCTION void apply(const MemberType &member, const XViewType &X, const YViewType &Y) const {
    if (!computed_inverse) {
      this->computeInverse<MemberType, ArgMode>(member);
      member.team_barrier();  // Finish writing to this->diag_values
    }

    KokkosBatched::HadamardProduct<MemberType, ArgMode>::template invoke<ValuesViewType, XViewType, YViewType>(
        member, diag_values, X, Y);
  }

  template <typename ArgTrans, int sameXY, typename XViewType, typename YViewType>
  KOKKOS_INLINE_FUNCTION void apply(const XViewType &X, const YViewType &Y) const {
    if (!computed_inverse) {
      this->computeInverse();
    }

    KokkosBatched::SerialHadamardProduct::template invoke<ValuesViewType, XViewType, YViewType>(diag_values, X, Y);
  }
};

}  // namespace KokkosBatched

#endif
