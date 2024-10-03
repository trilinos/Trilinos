// clang-format off
// @HEADER
// *****************************************************************************
//                            Tacho package
//
// Copyright 2022 NTESS and the Tacho contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// clang-format on
#ifndef __TACHO_LAPACK_SERIAL_HPP__
#define __TACHO_LAPACK_SERIAL_HPP__

/// \file  Tacho_Lapack_TEAM.hpp
/// \brief BLAS wrapper
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Kokkos_Core.hpp"

namespace Tacho {

template <typename T> struct LapackSerial {
#if 0
  struct Impl {
    template <typename MemberType>
    static KOKKOS_INLINE_FUNCTION void sytrf_lower(const MemberType &member, const int m, T *KOKKOS_RESTRICT A,
                                                   const int as0, const int as1, int *KOKKOS_RESTRICT ipiv, int *info) {
      *info = 0;
      if (m <= 0)
        return;

      using arith_traits = ArithTraits<T>;
      using mag_type = typename arith_traits::mag_type;

      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m), [&](const int &i) { ipiv[i] = i + 1; });

      const T zero(0);
      const int as = as0 + as1;
      const mag_type mu = 0.6404;
      for (int p = 0; p < m; ++p) {
        const int iend = m - p - 1;

        T *KOKKOS_RESTRICT alpha11 = A + (p)*as0 + (p)*as1,
          *KOKKOS_RESTRICT a21 = A + (p + 1) * as0 + (p)     * as1,
          *KOKKOS_RESTRICT A22 = A + (p + 1) * as0 + (p + 1) * as1;

        mag_type lambda1(0);
        int idx(0);
        {
          using reducer_value_type = typename Kokkos::MaxLoc<mag_type, int>::value_type;
          reducer_value_type value;
          Kokkos::MaxLoc<mag_type, int> reducer_value(value);
          Kokkos::parallel_reduce(
              Kokkos::TeamVectorRange(member, iend),
              [&](const int &i, reducer_value_type &update) {
                const mag_type val = arith_traits::abs(a21[i * as0]);
                if (val > update.val) {
                  update.val = val;
                  update.loc = i;
                }
              },
              reducer_value);
          member.team_barrier();

          lambda1 = value.val;
          idx = value.loc;
        }

        const mag_type abs_alpha = arith_traits::abs(*alpha11);
        if (abs_alpha < mu * lambda1) {
          mag_type lambda2(0);
          {
            using reducer_value_type = typename Kokkos::Max<mag_type>::value_type;
            reducer_value_type value;
            Kokkos::Max<mag_type> reducer_value(value);
            Kokkos::parallel_reduce(
                Kokkos::TeamVectorRange(member, iend),
                [&](const int &i, reducer_value_type &update) {
                  mag_type val(0);
                  if (i < idx)
                    val = arith_traits::abs(a21[idx * as0 + i * as1]);
                  else if (i > idx)
                    val = arith_traits::abs(A22[idx * as + (i - idx) * as0]);
                  if (val > update) {
                    update = val;
                  }
                },
                reducer_value);
            member.team_barrier();
            lambda2 = value;
          }
          const mag_type abs_alpha_idx = arith_traits::abs(A22[idx * as]);
          if (abs_alpha_idx * lambda2 < mu * lambda1 * lambda1) {
            /// pivot
            Kokkos::parallel_for(Kokkos::TeamVectorRange(member, iend), [&](const int &i) {
              if (i < idx)
                swap(a21[i * as0], A22[idx * as0 + i * as1]);
              else if (i > idx)
                swap(a21[i * as0], A22[i * as0 + idx * as1]);
              else {
                swap(alpha11[0], A22[idx * as]);
                ipiv[p] = p + idx + 2;
              }
            });
          }
        }
        member.team_barrier();
        Kokkos::single(Kokkos::PerThread(member), [&]() {
          if (*info == 0 && *alpha11 == zero) {
            *info = 1+p;
          }
        });
        const T alpha = *alpha11;
        Kokkos::parallel_for(Kokkos::TeamVectorRange(member, iend), [&](const int &i) { a21[i * as0] /= alpha; });
        member.team_barrier();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, iend), [&](const int &i) {
          const T aa = a21[i * as0];
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, i + 1), [&](const int &j) {
            const T bb = a21[j * as0];
            A22[i * as0 + j * as1] -= alpha * aa * bb;
          });
        });
        member.team_barrier();
      }
    }

    template <typename MemberType>
    static KOKKOS_INLINE_FUNCTION void sytrf_lower_nopiv(const MemberType &member, const int m, T *KOKKOS_RESTRICT A,
                                                         const int as0, const int as1, int *info) {
      *info = 0;
      if (m <= 0)
        return;

      // typedef ArithTraits<T> arith_traits;
      for (int p = 0; p < m; ++p) {
        const int iend = m - p - 1;

        T *KOKKOS_RESTRICT alpha11 = A + (p)*as0 + (p)*as1, *KOKKOS_RESTRICT a21 = A + (p + 1) * as0 + (p)*as1,
                        *KOKKOS_RESTRICT A22 = A + (p + 1) * as0 + (p + 1) * as1;

        const auto alpha = *alpha11; // arith_traits::real(*alpha11);
        Kokkos::parallel_for(Kokkos::TeamVectorRange(member, iend), [&](const int &i) { a21[i * as0] /= alpha; });
        member.team_barrier();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, iend), [&](const int &i) {
          const T aa = a21[i * as0];
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, i + 1), [&](const int &j) {
            const T bb = a21[j * as0];
            A22[i * as0 + j * as1] -= alpha * aa * bb;
          });
        });
        member.team_barrier();
      }
    }
  };

  template <typename MemberType>
  static KOKKOS_INLINE_FUNCTION void potrf(const MemberType &member, const char uplo, const int m,
                                           /* */ T *KOKKOS_RESTRICT A, const int lda, int *info) {
    switch (uplo) {
    case 'U':
    case 'u': {
      Impl::potrf_upper(member, m, A, 1, lda, info);
      break;
    }
    case 'L':
    case 'l': {
      Kokkos::abort("not implemented");
      break;
    }
    default:
      Kokkos::abort("Invalid uplo character");
    }
  }

  template <typename MemberType>
  static KOKKOS_INLINE_FUNCTION void sytrf(const MemberType &member, const char uplo, const int m,
                                           /* */ T *KOKKOS_RESTRICT A, const int lda,
                                           /* */ int *KOKKOS_RESTRICT P,
                                           /* */ T *KOKKOS_RESTRICT W, int *info) {
    switch (uplo) {
    case 'U':
    case 'u': {
      Kokkos::abort("not implemented");
      break;
    }
    case 'L':
    case 'l': {
      Impl::sytrf_lower(member, m, A, 1, lda, P, info);
      break;
    }
    default:
      Kokkos::abort("Invalid uplo character");
    }
  }
#endif

  inline static int potrf(const char uplo, const int m, T *A, const int lda, int *info) {

    *info = 0;
    if (m <= 0)
      return 0;

    typedef ArithTraits<T> arith_traits;
    const typename arith_traits::mag_type zero(0);

    for (int i = 0; i < m; ++i) {
      // check pivot
      T alpha = A[i + i*lda];
      if (arith_traits::real(alpha) <= zero) {
        *info = 1+i;
        return *info;
      }
      alpha = sqrt(arith_traits::real(alpha));
      alpha = arith_traits::real(alpha);
      A[i + i*lda] = alpha;

      // scale
      for (int j = i+1; j < m; j++) { A[i + j*lda] /= alpha; }

      // update
      for (int j = i+1; j < m; j++) {
        const T aa = A[i + j*lda];
        for (int l = i+1; l <= j; l++) {
          const T bb = arith_traits::conj(A[i + l*lda]);
          A[l + j*lda] -= aa * bb;
        }
      }
    }
    return 0;
  }

  inline static int getrf(const int m, const int n, T *  A,
                          const int lda, int * ipiv, int *info) {

    *info = 0;
    if (m <= 0 || n <= 0)
      return *info;


    using arith_traits = ArithTraits<T>;
    using mag_type = typename arith_traits::mag_type;
    const T zero(0);

    int mn = (m < n ? m : n);
    for (int j = 0; j < mn; ++j) {
      // look for pivot
      int idx = j;
      mag_type alpha = arith_traits::abs(A[j + j*lda]);
      for (int i = j + 1; i < m; i++) {
        mag_type val = arith_traits::abs(A[i + j*lda]);
        if (val > alpha) {
          alpha = val;
          idx = i;
        }
      }
      if (alpha == zero) {
        *info = j*1;
        return *info;
      }
      ipiv[j] = idx+1;

      //swap
      if (idx != j) {
        for (int k = 0; k < n; k++) {
          T val = A[j + k*lda];
          A[j + k*lda] = A[idx + k*lda];
          A[idx + k*lda] = val;
        }
      }

      // scale
      for (int i = j+1; i < m; i++) {
        A[i + j*lda] /= A[j + j*lda];
      }

      // update
      for (int k = j+1; k < n; k++) {
        for (int i = j+1; i < m; i++) {
          A[i + k*lda] -= A[i + j*lda] * A[j + k*lda];
        }
      }
    }
    return *info;
  }
};

} // namespace Tacho

#endif
