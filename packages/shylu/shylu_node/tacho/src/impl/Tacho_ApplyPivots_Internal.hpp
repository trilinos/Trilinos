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
#ifndef __TACHO_APPLY_PIVOTS_INTERNAL_HPP__
#define __TACHO_APPLY_PIVOTS_INTERNAL_HPP__

/// \file  Tacho_ApplyPivots_Internal.hpp
/// \brief Apply pivots
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

/// row exchange
template <> struct ApplyPivots<PivotMode::Flame, Side::Left, Direct::Forward, Algo::Internal> {
  template <typename ViewTypeP, typename ViewTypeA>
  KOKKOS_INLINE_FUNCTION static int invoke(const ViewTypeP &P, const ViewTypeA &A) {
    typedef typename ViewTypeA::non_const_value_type value_type;

    if (A.extent(0) == P.extent(0)) {
      if (A.span() > 0) {
        const ordinal_type m = A.extent(0), n = A.extent(1);
        for (ordinal_type i = 0; i < m; ++i) {
          const ordinal_type piv = P(i);
          if (piv == 0) {
            /// no pivot
          } else {
            /// 1x1 pivot
            /// 2x2 pivots are already converted to 1x1 type
            const ordinal_type p = i + piv;
            for (ordinal_type j = 0; j < n; ++j) {
              const value_type tmp = A(i, j);
              A(i, j) = A(p, j);
              A(p, j) = tmp;
            }
          }
        }
      }
    } else {
      Kokkos::printf("Error: ApplyPivots<Algo::Internal> A is not square\n");
    }
    return 0;
  }

  template <typename MemberType, typename ViewTypeP, typename ViewTypeA>
  KOKKOS_INLINE_FUNCTION static int invoke(MemberType &member, const ViewTypeP &P, const ViewTypeA &A) {
    KOKKOS_IF_ON_DEVICE((
    typedef typename ViewTypeA::non_const_value_type value_type;

    if (A.extent(0) == P.extent(0)) {
      if (A.span() > 0) {
        const ordinal_type m = A.extent(0), n = A.extent(1);
        Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n), [&](const ordinal_type &j) {
          for (ordinal_type i = 0; i < m; ++i) {
            const ordinal_type piv = P(i);
            if (piv == 0) {
              /// no pivot
            } else {
              /// 1x1 pivot
              const ordinal_type p = i + piv;
              const value_type tmp = A(i, j);
              A(i, j) = A(p, j);
              A(p, j) = tmp;
            }
          }
        });
      }
    } else {
      Kokkos::printf("Error: ApplyPivots<Algo::Internal> A is not square\n");
    }))
    KOKKOS_IF_ON_HOST((invoke(P, A);))
    return 0;
  }
};

/// row exchange
template <> struct ApplyPivots<PivotMode::Flame, Side::Left, Direct::Backward, Algo::Internal> {
  template <typename ViewTypeP, typename ViewTypeA>
  KOKKOS_INLINE_FUNCTION static int invoke(const ViewTypeP &P, const ViewTypeA &A) {
    typedef typename ViewTypeA::non_const_value_type value_type;

    if (A.extent(0) == P.extent(0)) {
      if (A.span() > 0) {
        const ordinal_type m = A.extent(0), n = A.extent(1);
        for (ordinal_type i = (m - 1); i >= 0; --i) {
          const ordinal_type piv = P(i);
          if (piv == 0) {
            /// no pivot
          } else {
            /// 1x1 pivot
            /// 2x2 pivots are already converted to 1x1 type
            const ordinal_type p = i + piv;
            for (ordinal_type j = 0; j < n; ++j) {
              const value_type tmp = A(i, j);
              A(i, j) = A(p, j);
              A(p, j) = tmp;
            }
          }
        }
      }
    } else {
      Kokkos::printf("Error: ApplyPivots<Algo::Internal> A is not square\n");
    }
    return 0;
  }

  template <typename MemberType, typename ViewTypeP, typename ViewTypeA>
  KOKKOS_INLINE_FUNCTION static int invoke(MemberType &member, const ViewTypeP &P, const ViewTypeA &A) {
    KOKKOS_IF_ON_DEVICE((
    typedef typename ViewTypeA::non_const_value_type value_type;

    if (A.extent(0) == P.extent(0)) {
      if (A.span() > 0) {
        const ordinal_type m = A.extent(0), n = A.extent(1);
        Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n), [&](const ordinal_type &j) {
          for (ordinal_type i = (m - 1); i >= 0; --i) {
            const ordinal_type piv = P(i);
            if (piv == 0) {
              /// no pivot
            } else {
              /// 1x1 pivot
              const ordinal_type p = i + piv;
              const value_type tmp = A(i, j);
              A(i, j) = A(p, j);
              A(p, j) = tmp;
            }
          }
        });
      }
    } else {
      Kokkos::printf("Error: ApplyPivots<Algo::Internal> A is not square\n");
    }))
    KOKKOS_IF_ON_HOST((invoke(P, A);))
    return 0;
  }
};

//   template<>
//   struct ApplyPivots<PivotMode::Lapack,Side::Left,Direct::Forward,Algo::Internal> {
//     template<typename ViewTypeP,
//              typename ViewTypeA>
//     KOKKOS_INLINE_FUNCTION
//     static int
//     invoke(const ViewTypeP &P,
//            const ViewTypeA &A) {
//       typedef typename ViewTypeA::non_const_value_type value_type;

//       if (A.extent(0) == P.extent(0)) {
//         if (A.span() > 0) {
//           const ordinal_type m = A.extent(0), n = A.extent(1);
//           for (ordinal_type i=0;i<m;++i) {
//             const ordinal_type piv = P(i);
//             if (piv == i || piv == 0) {
//               /// no pivot
//             } else if (piv > 0) {
//               /// 1x1 pivot
//               const ordinal_type p = piv - 1;
//               for (ordinal_type j=0;j<n;++j) {
//                 const value_type tmp = A(i,j);
//                 A(i,j) = A(p,j);
//                 A(p,j) = tmp;
//               }
//             } else {
//               /// 2x2 pivot
//               const int p = -piv - 1;
//               for (ordinal_type j=0;j<n;++j) {
//                 const value_type tmp_a = A(i,j), tmp_b = A(i+1,j);
//                 A(i  ,j) = A(p  ,j);
//                 A(i+1,j) = A(p+1,j);
//                 A(p,  j) = tmp_a;
//                 A(p+1,j) = tmp_b;
//               }
//             }
//           }
//         }
//       } else {
//         Kokkos::printf("Error: ApplyPivots<Algo::Internal> A is not square\n");
//       }
//       return 0;
//     }

//     template<typename MemberType,
//              typename ViewTypeP,
//              typename ViewTypeA>
//     KOKKOS_INLINE_FUNCTION
//     static int
//     invoke(MemberType &member,
//            const ViewTypeP &P,
//            const ViewTypeA &A) {
//       KOKKOS_IF_ON_DEVICE((
//       typedef typename ViewTypeA::non_const_value_type value_type;

//       if (A.extent(0) == P.extent(0)) {
//         if (A.span() > 0) {
//           const ordinal_type m = A.extent(0), n = A.extent(1);
//           Kokkos::parallel_for
//             (Kokkos::TeamVectorRange(member, n),
//              [&](const ordinal_type &j) {
//               for (ordinal_type i=0;i<m;++i) {
//                 const ordinal_type piv = P(i);
//                 if (piv == i || piv == 0) {
//                   /// no pivot
//                 } else if (piv > 0) {
//                   /// 1x1 pivot
//                   const ordinal_type p = piv - 1;
//                   const value_type tmp = A(i,j);
//                   A(i,j) = A(p,j);
//                   A(p,j) = tmp;
//                 } else {
//                   /// 2x2 pivot
//                   const int p = -piv - 1;
//                   const value_type tmp_a = A(i,j), tmp_b = A(i+1,j);
//                   A(i  ,j) = A(p  ,j);
//                   A(i+1,j) = A(p+1,j);
//                   A(p,  j) = tmp_a;
//                   A(p+1,j) = tmp_b;
//                 }
//               }
//             });
//         }
//       } else {
//         Kokkos::printf("Error: ApplyPivots<Algo::Internal> A is not square\n");
//       }))
//       KOKKOS_IF_ON_HOST((invoke(P, A);))
//       return 0;
//     }
//   };

} // namespace Tacho
#endif
