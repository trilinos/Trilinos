// clang-format off
/* =====================================================================================
Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
certain rights in this software.

SCR#:2790.0

This file is part of Tacho. Tacho is open source software: you can redistribute it
and/or modify it under the terms of BSD 2-Clause License
(https://opensource.org/licenses/BSD-2-Clause). A copy of the licese is also
provided under the main directory

Questions? Kyungjoo Kim at <kyukim@sandia.gov,https://github.com/kyungjoo-kim>

Sandia National Laboratories, Albuquerque, NM, USA
===================================================================================== */
// clang-format on
#ifndef __TACHO_LAPACK_TEAM_HPP__
#define __TACHO_LAPACK_TEAM_HPP__

/// \file  Tacho_Lapack_TEAM.hpp
/// \brief BLAS wrapper
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Kokkos_Core.hpp"

namespace Tacho {

template <typename T> struct LapackTeam {
  struct Impl {
    template <typename MemberType>
    static KOKKOS_INLINE_FUNCTION void potrf_upper(const MemberType &member, const int m, T *__restrict__ A,
                                                   const int as0, const int as1, int *info) {
      if (m <= 0)
        return;

      typedef ArithTraits<T> arith_traits;
      for (int p = 0; p < m; ++p) {
        const int jend = m - p - 1;

        T *__restrict__ alpha11 = A + (p)*as0 + (p)*as1, *__restrict__ a12t = A + (p)*as0 + (p + 1) * as1,
                        *__restrict__ A22 = A + (p + 1) * as0 + (p + 1) * as1;

        Kokkos::single(Kokkos::PerTeam(member), [&]() { *alpha11 = sqrt(arith_traits::real(*alpha11)); });
        member.team_barrier();
        const auto alpha = arith_traits::real(*alpha11);
        Kokkos::parallel_for(Kokkos::TeamVectorRange(member, jend), [&](const int &j) { a12t[j * as1] /= alpha; });
        member.team_barrier();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, jend), [&](const int &j) {
          const T aa = arith_traits::conj(a12t[j * as1]);
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, j + 1), [&](const int &i) {
            const T bb = a12t[i * as1];
            A22[i * as0 + j * as1] -= aa * bb;
          });
        });
        member.team_barrier();
      }
    }

    template <typename MemberType>
    static KOKKOS_INLINE_FUNCTION void sytrf_lower(const MemberType &member, const int m, T *__restrict__ A,
                                                   const int as0, const int as1, int *__restrict__ ipiv, int *info) {
      if (m <= 0)
        return;

      using arith_traits = ArithTraits<T>;
      using mag_type = typename arith_traits::mag_type;

      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m), [&](const int &i) { ipiv[i] = i + 1; });

      const int as = as0 + as1;
      const mag_type mu = 0.6404;
      for (int p = 0; p < m; ++p) {
        const int iend = m - p - 1;

        T *__restrict__ alpha11 = A + (p)*as0 + (p)*as1, *__restrict__ a21 = A + (p + 1) * as0 + (p)*as1,
                        *__restrict__ A22 = A + (p + 1) * as0 + (p + 1) * as1;

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
    static KOKKOS_INLINE_FUNCTION void sytrf_lower_nopiv(const MemberType &member, const int m, T *__restrict__ A,
                                                         const int as0, const int as1, int *info) {
      if (m <= 0)
        return;

      // typedef ArithTraits<T> arith_traits;
      for (int p = 0; p < m; ++p) {
        const int iend = m - p - 1;

        T *__restrict__ alpha11 = A + (p)*as0 + (p)*as1, *__restrict__ a21 = A + (p + 1) * as0 + (p)*as1,
                        *__restrict__ A22 = A + (p + 1) * as0 + (p + 1) * as1;

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
                                           /* */ T *__restrict__ A, const int lda, int *info) {
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
                                           /* */ T *__restrict__ A, const int lda,
                                           /* */ int *__restrict__ P,
                                           /* */ T *__restrict__ W, int *info) {
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

  template <typename MemberType>
  static KOKKOS_INLINE_FUNCTION void getrf(const MemberType &member, const int m, const int n, T *__restrict__ A,
                                           const int as1, int *__restrict__ ipiv, int *info) {
    if (m <= 0 || n <= 0)
      return;

    using arith_traits = ArithTraits<T>;
    using mag_type = typename arith_traits::mag_type;

    const int as0 = 1;
    for (int p = 0; p < m; ++p) {
      const int iend = m - p - 1, jend = n - p - 1;
      T *__restrict__ alpha11 = A + (p)*as0 + (p)*as1, *__restrict__ AB = A + (p)*as0, *__restrict__ ABR = alpha11,
                      *__restrict__ a21 = A + (p + 1) * as0 + (p)*as1, *__restrict__ a12 = A + (p)*as0 + (p + 1) * as1,
                      *__restrict__ A22 = A + (p + 1) * as0 + (p + 1) * as1;

      {
        int idx(0);
        using reducer_value_type = typename Kokkos::MaxLoc<mag_type, int>::value_type;
        reducer_value_type value;
        Kokkos::MaxLoc<mag_type, int> reducer_value(value);
        Kokkos::parallel_reduce(
            Kokkos::TeamVectorRange(member, 1 + iend),
            [&](const int &i, reducer_value_type &update) {
              const mag_type val = arith_traits::abs(ABR[i * as0]);
              if (val > update.val) {
                update.val = val;
                update.loc = i;
              }
            },
            reducer_value);
        member.team_barrier();
        idx = value.loc;

        /// pivot
        Kokkos::single(Kokkos::PerThread(member), [&]() { ipiv[p] = p + idx + 1; });
        if (idx) {
          Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n),
                               [&](const int &j) { swap(AB[j * as1], AB[idx * as0 + j * as1]); });
          member.team_barrier();
        }
      }

      member.team_barrier();
      const T alpha = *alpha11;
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, iend), [&](const int &i) { a21[i * as0] /= alpha; });
      member.team_barrier();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, jend), [&](const int &j) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, iend),
                             [&](const int &i) { A22[i * as0 + j * as1] -= a21[i * as0] * a12[j * as1]; });
      });
      member.team_barrier();
    }
  }
};

} // namespace Tacho

#endif
