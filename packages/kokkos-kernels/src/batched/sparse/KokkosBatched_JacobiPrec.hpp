//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.4
//       Copyright (2021) National Technology & Engineering
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
#ifndef __KOKKOSBATCHED_JACOBIPREC_HPP__
#define __KOKKOSBATCHED_JACOBIPREC_HPP__

/// \author Kim Liegeois (knliege@sandia.gov)

#include "KokkosBatched_HadamardProduct.hpp"

namespace KokkosBatched {

/// \brief Batched Jacobi Preconditioner:
///
/// \tparam ValuesViewType: Input type for the values of the diagonal

template <class ValuesViewType>
class JacobiPrec {
 public:
  using ScalarType = typename ValuesViewType::non_const_value_type;
  using MagnitudeType =
      typename Kokkos::Details::ArithTraits<ScalarType>::mag_type;

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

  template <typename MemberType, typename ArgMode>
  KOKKOS_INLINE_FUNCTION void computeInverse(const MemberType &member) const {
    auto one     = Kokkos::Details::ArithTraits<MagnitudeType>::one();
    auto epsilon = Kokkos::Details::ArithTraits<MagnitudeType>::epsilon();
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
      auto vs0               = diag_values.stride_0();
      auto vs1               = diag_values.stride_1();

      Kokkos::parallel_reduce(
          Kokkos::TeamThreadRange(member, 0, n_operators * n_rows),
          [&](const int &iTemp, int &ltooSmall) {
            int i, j;
            getIndices<int, typename ValuesViewType::array_layout>(
                iTemp, n_rows, n_operators, j, i);
            if (Kokkos::abs<ScalarType>(diag_values_array[i * vs0 + j * vs1]) <=
                epsilon) {
              ltooSmall++;
              diag_values_array[i * vs0 + j * vs1] = one;
            } else
              diag_values_array[i * vs0 + j * vs1] =
                  one / diag_values_array[i * vs0 + j * vs1];
          },
          tooSmall);
    } else if (std::is_same<ArgMode, Mode::TeamVector>::value) {
      auto diag_values_array = diag_values.data();
      auto vs0               = diag_values.stride_0();
      auto vs1               = diag_values.stride_1();

      Kokkos::parallel_reduce(
          Kokkos::TeamVectorRange(member, 0, n_operators * n_rows),
          [&](const int &iTemp, int &ltooSmall) {
            int i, j;
            getIndices<int, typename ValuesViewType::array_layout>(
                iTemp, n_rows, n_operators, j, i);
            if (Kokkos::abs<ScalarType>(diag_values_array[i * vs0 + j * vs1]) <=
                epsilon) {
              ltooSmall++;
              diag_values_array[i * vs0 + j * vs1] = one;
            } else
              diag_values_array[i * vs0 + j * vs1] =
                  one / diag_values_array[i * vs0 + j * vs1];
          },
          tooSmall);
    }

    if (tooSmall > 0)
      KOKKOS_IMPL_DO_NOT_USE_PRINTF(
          "KokkosBatched::JacobiPrec: %d entrie(s) has/have a too small "
          "magnitude and have been replaced by one, \n",
          (int)tooSmall);
    computed_inverse = true;
  }

  template <typename MemberType, typename XViewType, typename YViewType,
            typename ArgTrans, typename ArgMode, int sameXY>
  KOKKOS_INLINE_FUNCTION void apply(const MemberType &member,
                                    const XViewType &X,
                                    const YViewType &Y) const {
    if (!computed_inverse) {
      this->computeInverse<MemberType, ArgMode>(member);
      member.team_barrier();  // Finish writing to this->diag_values
    }

    KokkosBatched::HadamardProduct<MemberType, ArgMode>::template invoke<
        ValuesViewType, XViewType, YViewType>(member, diag_values, X, Y);
  }
};

}  // namespace KokkosBatched

#endif