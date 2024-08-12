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
#ifndef __KOKKOSBATCHED_GESV_IMPL_HPP__
#define __KOKKOSBATCHED_GESV_IMPL_HPP__

/// \author Kim Liegeois (knliege@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include <KokkosBatched_LU_Decl.hpp>
#include "KokkosBatched_Trsm_Decl.hpp"
#include "KokkosBatched_Copy_Decl.hpp"

namespace KokkosBatched {

struct SerialStaticPivoting {
  template <class MatrixType1, class MatrixType2, class VectorType1, class VectorType2>
  KOKKOS_INLINE_FUNCTION static int invoke(const MatrixType1 A, const MatrixType2 PDAD, const VectorType1 Y,
                                           const VectorType2 PDY, const VectorType2 D2, const VectorType2 tmp_v_1,
                                           const VectorType2 tmp_v_2);
};

template <typename MemberType>
struct TeamStaticPivoting {
  template <class MatrixType1, class MatrixType2, class VectorType1, class VectorType2>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const MatrixType1 A, const MatrixType2 PDAD,
                                           const VectorType1 Y, const VectorType2 PDY, const VectorType2 D2,
                                           const VectorType2 tmp_v_1, const VectorType2 tmp_v_2);
};

template <typename MemberType>
struct TeamVectorStaticPivoting {
  template <class MatrixType1, class MatrixType2, class VectorType1, class VectorType2>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const MatrixType1 A, const MatrixType2 PDAD,
                                           const VectorType1 Y, const VectorType2 PDY, const VectorType2 D2,
                                           const VectorType2 tmp_v_1, const VectorType2 tmp_v_2);
};

template <class MatrixType1, class MatrixType2, class VectorType1, class VectorType2>
KOKKOS_INLINE_FUNCTION int SerialStaticPivoting::invoke(const MatrixType1 A, const MatrixType2 PDAD,
                                                        const VectorType1 Y, const VectorType2 PDY,
                                                        const VectorType2 D2, const VectorType2 tmp_v_1,
                                                        const VectorType2 tmp_v_2) {
  using value_type = typename MatrixType1::non_const_value_type;
  const size_t n   = A.extent(0);

  // First, the algorithm loops over the rows and columns and search
  // for the maximal absolute value per row and column.
  for (size_t i = 0; i < n; ++i) {
    D2(i)      = Kokkos::ArithTraits<value_type>::zero();
    tmp_v_1(i) = 0;
    tmp_v_2(i) = 1.;
    for (size_t j = 0; j < n; ++j) {
      if (D2(i) < Kokkos::abs(A(j, i))) D2(i) = Kokkos::abs(A(j, i));
      if (tmp_v_1(i) < Kokkos::abs(A(i, j))) tmp_v_1(i) = Kokkos::abs(A(i, j));
    }
    D2(i) = 1. / D2(i);
  }

  // Then, the inverse of the maximal value per column is used to scale
  // A by the right.
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      A(i, j) *= D2(j);
    }
  }

  // Once again, the algorithm loops over the rows and store the maximal
  // absolute value per row but after the right scalling and do a left scalling
  // of A and Y.
  value_type D1_i;
  for (size_t i = 0; i < n; ++i) {
    D1_i = Kokkos::ArithTraits<value_type>::zero();
    for (size_t j = 0; j < n; ++j) {
      if (D1_i < Kokkos::abs(A(i, j))) D1_i = Kokkos::abs(A(i, j));
    }
    D1_i = 1. / D1_i;
    for (size_t j = 0; j < n; ++j) {
      A(i, j) *= D1_i;
    }
    Y(i) *= D1_i;
  }

  // Finally, the algorithm starts to loop over the rows in an order such that
  // their initial maximal absolute value decrease (it uses the tmp_v_1 to do
  // so), then for a given row, it finds the available column with the largest
  // absolute value. If this value is zero, the algorithm failed to compute a
  // good pivot, otherwise it puts the current row to the found column index and
  // it labels the row and column index as unavailable and continue the loop
  // over the rows.
  //
  for (size_t i = 0; i < n; ++i) {
    int row_index    = 0;
    int col_index    = 0;
    value_type tmp_0 = Kokkos::ArithTraits<value_type>::zero();
    value_type tmp_1 = Kokkos::ArithTraits<value_type>::zero();
    for (size_t j = 0; j < n; ++j) {
      if (tmp_0 < tmp_v_1(j)) {
        tmp_0     = tmp_v_1(j);
        row_index = j;
      }
    }
    for (size_t j = 0; j < n; ++j) {
      if (tmp_1 < Kokkos::abs(A(row_index, j) * tmp_v_2(j))) {
        tmp_1     = Kokkos::abs(A(row_index, j) * tmp_v_2(j));
        col_index = j;
      }
    }
    if (tmp_1 == Kokkos::ArithTraits<value_type>::zero()) return 1;
    tmp_v_1(row_index) = Kokkos::ArithTraits<value_type>::zero();
    tmp_v_2(col_index) = Kokkos::ArithTraits<value_type>::zero();

    for (size_t j = 0; j < n; ++j) {
      PDAD(col_index, j) = A(row_index, j);
    }
    PDY(col_index) = Y(row_index);
  }

  return 0;
}

template <typename MemberType>
template <class MatrixType1, class MatrixType2, class VectorType1, class VectorType2>
KOKKOS_INLINE_FUNCTION int TeamStaticPivoting<MemberType>::invoke(const MemberType &member, const MatrixType1 A,
                                                                  const MatrixType2 PDAD, const VectorType1 Y,
                                                                  const VectorType2 PDY, const VectorType2 D2,
                                                                  const VectorType2 tmp_v_1,
                                                                  const VectorType2 tmp_v_2) {
  using value_type         = typename MatrixType1::non_const_value_type;
  using reducer_value_type = typename Kokkos::MaxLoc<value_type, int>::value_type;
  // This implementation follows the strategy of SerialStaticPivoting but uses
  // an extra level of parallelism.

  // Made this non-const in order to WORKAROUND issue #349 (Credit to C. Trott)
  size_t n = A.extent(0);

  Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int &i) {
    D2(i)      = Kokkos::ArithTraits<value_type>::zero();
    tmp_v_1(i) = 0;
    tmp_v_2(i) = 1.;
    for (size_t j = 0; j < n; ++j) {
      if (D2(i) < Kokkos::abs(A(j, i))) D2(i) = Kokkos::abs(A(j, i));
      if (tmp_v_1(i) < Kokkos::abs(A(i, j))) tmp_v_1(i) = Kokkos::abs(A(i, j));
    }
    D2(i) = 1. / D2(i);
  });

  Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int &i) {
    for (size_t j = 0; j < n; ++j) {
      A(i, j) *= D2(j);
    }
  });

  Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int &i) {
    value_type D1_i = Kokkos::ArithTraits<value_type>::zero();
    for (size_t j = 0; j < n; ++j) {
      if (D1_i < Kokkos::abs(A(i, j))) D1_i = Kokkos::abs(A(i, j));
    }
    D1_i = 1. / D1_i;
    for (size_t j = 0; j < n; ++j) {
      A(i, j) *= D1_i;
    }
    Y(i) *= D1_i;
  });

  for (size_t i = 0; i < n; ++i) {
    int row_index, col_index;
    reducer_value_type value;
    Kokkos::MaxLoc<value_type, int> reducer_value(value);
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(member, n),
        [&](const int &j, reducer_value_type &update) {
          if (tmp_v_1(j) > update.val) {
            update.val = tmp_v_1(j);
            update.loc = j;
          }
        },
        reducer_value);
    row_index = value.loc;
    value.loc = 0;
    value.val = Kokkos::ArithTraits<value_type>::zero();
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(member, n),
        [&](const int &j, reducer_value_type &update) {
          if (Kokkos::abs(A(row_index, j) * tmp_v_2(j)) > update.val) {
            update.val = Kokkos::abs(A(row_index, j) * tmp_v_2(j));
            update.loc = j;
          }
        },
        reducer_value);
    col_index = value.loc;
    if (value.val == Kokkos::ArithTraits<value_type>::zero()) return 1;
    tmp_v_1(row_index) = Kokkos::ArithTraits<value_type>::zero();
    tmp_v_2(col_index) = Kokkos::ArithTraits<value_type>::zero();

    for (size_t j = 0; j < n; ++j) {
      PDAD(col_index, j) = A(row_index, j);
    }
    PDY(col_index) = Y(row_index);
  }
  return 0;
}

template <typename MemberType>
template <class MatrixType1, class MatrixType2, class VectorType1, class VectorType2>
KOKKOS_INLINE_FUNCTION int TeamVectorStaticPivoting<MemberType>::invoke(const MemberType &member, const MatrixType1 A,
                                                                        const MatrixType2 PDAD, const VectorType1 Y,
                                                                        const VectorType2 PDY, const VectorType2 D2,
                                                                        const VectorType2 tmp_v_1,
                                                                        const VectorType2 tmp_v_2) {
  using value_type         = typename MatrixType1::non_const_value_type;
  using reducer_value_type = typename Kokkos::MaxLoc<value_type, int>::value_type;
  // This implementation follows the strategy of SerialStaticPivoting but uses
  // two extra levels of parallelism.

  const size_t n = A.extent(0);

  Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int &i) {
    D2(i)      = Kokkos::ArithTraits<value_type>::zero();
    tmp_v_1(i) = 0;
    tmp_v_2(i) = 1.;
    reducer_value_type value;
    Kokkos::MaxLoc<value_type, int> reducer_value(value);
    Kokkos::parallel_reduce(
        Kokkos::ThreadVectorRange(member, n),
        [&](const int &j, reducer_value_type &update) {
          if (Kokkos::abs(A(j, i)) > update.val) {
            update.val = Kokkos::abs(A(j, i));
            update.loc = j;
          }
        },
        reducer_value);
    D2(i) = 1. / value.val;
    Kokkos::parallel_reduce(
        Kokkos::ThreadVectorRange(member, n),
        [&](const int &j, reducer_value_type &update) {
          if (Kokkos::abs(A(i, j)) > update.val) {
            update.val = Kokkos::abs(A(i, j));
            update.loc = j;
          }
        },
        reducer_value);
    tmp_v_1(i) = value.val;
  });

  Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int &i) {
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, n), [&](const int &j) { A(i, j) *= D2(j); });
  });

  Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int &i) {
    value_type D1_i = Kokkos::ArithTraits<value_type>::zero();
    reducer_value_type value;
    Kokkos::MaxLoc<value_type, int> reducer_value(value);
    Kokkos::parallel_reduce(
        Kokkos::ThreadVectorRange(member, n),
        [&](const int &j, reducer_value_type &update) {
          if (Kokkos::abs(A(i, j)) > update.val) {
            update.val = Kokkos::abs(A(i, j));
            update.loc = j;
          }
        },
        reducer_value);
    D1_i = 1. / value.val;
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, n), [&](const int &j) { A(i, j) *= D1_i; });
    Y(i) *= D1_i;
  });

  for (size_t i = 0; i < n; ++i) {
    int row_index, col_index;
    reducer_value_type value{};
    Kokkos::MaxLoc<value_type, int> reducer_value(value);
    Kokkos::parallel_reduce(
        Kokkos::TeamVectorRange(member, n),
        [&](const int &j, reducer_value_type &update) {
          if (tmp_v_1(j) > update.val) {
            update.val = tmp_v_1(j);
            update.loc = j;
          }
        },
        reducer_value);
    row_index = value.loc;
    value.loc = 0;
    value.val = Kokkos::ArithTraits<value_type>::zero();
    Kokkos::parallel_reduce(
        Kokkos::TeamVectorRange(member, n),
        [&](const int &j, reducer_value_type &update) {
          if (Kokkos::abs(A(row_index, j) * tmp_v_2(j)) > update.val) {
            update.val = Kokkos::abs(A(row_index, j) * tmp_v_2(j));
            update.loc = j;
          }
        },
        reducer_value);
    col_index = value.loc;
    if (value.val == Kokkos::ArithTraits<value_type>::zero()) return 1;
    tmp_v_1(row_index) = Kokkos::ArithTraits<value_type>::zero();
    tmp_v_2(col_index) = Kokkos::ArithTraits<value_type>::zero();

    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n),
                         [&](const int &j) { PDAD(col_index, j) = A(row_index, j); });
    PDY(col_index) = Y(row_index);
  }
  return 0;
}

template <class VectorType1, class VectorType2, class VectorType3>
KOKKOS_INLINE_FUNCTION void SerialHadamard1D(const VectorType1 X, const VectorType2 D, const VectorType3 DX) {
  const size_t n = X.extent(0);

  for (size_t i = 0; i < n; ++i) {
    DX(i) = D(i) * X(i);
  }
}

template <typename MemberType, class VectorType1, class VectorType2, class VectorType3>
KOKKOS_INLINE_FUNCTION void TeamHadamard1D(const MemberType &member, const VectorType1 X, const VectorType2 D,
                                           const VectorType3 DX) {
  const size_t n = X.extent(0);

  Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const size_t &i) { DX(i) = D(i) * X(i); });
}

template <typename MemberType, class VectorType1, class VectorType2, class VectorType3>
KOKKOS_INLINE_FUNCTION void TeamVectorHadamard1D(const MemberType &member, const VectorType1 X, const VectorType2 D,
                                                 const VectorType3 DX) {
  const size_t n = X.extent(0);

  Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n), [&](const size_t &i) { DX(i) = D(i) * X(i); });
}

///
/// Serial Impl
/// ===========
template <>
struct SerialGesv<Gesv::StaticPivoting> {
  template <typename MatrixType, typename XVectorType, typename YVectorType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MatrixType A, const XVectorType X, const YVectorType Y,
                                           const MatrixType tmp) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
    static_assert(Kokkos::is_view<MatrixType>::value, "KokkosBatched::gesv: MatrixType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<XVectorType>::value, "KokkosBatched::gesv: XVectorType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<YVectorType>::value, "KokkosBatched::gesv: YVectorType is not a Kokkos::View.");
    static_assert(MatrixType::rank == 2, "KokkosBatched::gesv: MatrixType must have rank 2.");
    static_assert(XVectorType::rank == 1, "KokkosBatched::gesv: XVectorType must have rank 1.");
    static_assert(YVectorType::rank == 1, "KokkosBatched::gesv: YVectorType must have rank 1.");

    // Check compatibility of dimensions at run time.

    if (A.extent(0) != tmp.extent(0) || A.extent(1) + 4 != tmp.extent(1)) {
      Kokkos::printf(
          "KokkosBatched::gesv: dimensions of A and tmp do not match: A: "
          "%d x %d, tmp (note: its second dimension should be the second "
          "dimension of A + 4): %d x %d\n",
          (int)A.extent(0), (int)A.extent(1), (int)tmp.extent(0), (int)tmp.extent(1));
      return 1;
    }

    if (A.extent(0) != X.extent(0) || A.extent(1) != X.extent(0) || A.extent(0) != Y.extent(0)) {
      Kokkos::printf(
          "KokkosBatched::gesv: dimensions of A and X and Y do not match: A: "
          "%d x %d, X: %d, Y: %d\n",
          (int)A.extent(0), (int)A.extent(1), (int)X.extent(0), (int)Y.extent(0));
      return 1;
    }
#endif

    const int n = A.extent(0);

    auto PDAD    = Kokkos::subview(tmp, Kokkos::ALL, Kokkos::make_pair(0, n));
    auto PDY     = Kokkos::subview(tmp, Kokkos::ALL, n);
    auto D2      = Kokkos::subview(tmp, Kokkos::ALL, n + 1);
    auto tmp_v_1 = Kokkos::subview(tmp, Kokkos::ALL, n + 2);
    auto tmp_v_2 = Kokkos::subview(tmp, Kokkos::ALL, n + 3);

    if (SerialStaticPivoting::invoke(A, PDAD, Y, PDY, D2, tmp_v_1, tmp_v_2) == 1) {
      Kokkos::printf(
          "KokkosBatched::gesv: the currently implemented static pivoting "
          "failed.\n");
      return 1;
    }

    int r_val = SerialLU<Algo::Level3::Unblocked>::invoke(PDAD);

    if (r_val == 0)
      r_val = SerialTrsm<Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::Unit, Algo::Level3::Unblocked>::invoke(
          1.0, PDAD, PDY);

    if (r_val == 0)
      r_val = SerialTrsm<Side::Left, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit, Algo::Level3::Unblocked>::invoke(
          1.0, PDAD, PDY);

    if (r_val == 0) SerialHadamard1D(PDY, D2, X);
    return r_val;
  }
};

template <>
struct SerialGesv<Gesv::NoPivoting> {
  template <typename MatrixType, typename XVectorType, typename YVectorType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MatrixType A, const XVectorType X, const YVectorType Y,
                                           const MatrixType /*tmp*/) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
    static_assert(Kokkos::is_view<MatrixType>::value, "KokkosBatched::gesv: MatrixType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<XVectorType>::value, "KokkosBatched::gesv: XVectorType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<YVectorType>::value, "KokkosBatched::gesv: YVectorType is not a Kokkos::View.");
    static_assert(MatrixType::rank == 2, "KokkosBatched::gesv: MatrixType must have rank 2.");
    static_assert(XVectorType::rank == 1, "KokkosBatched::gesv: XVectorType must have rank 1.");
    static_assert(YVectorType::rank == 1, "KokkosBatched::gesv: YVectorType must have rank 1.");

    // Check compatibility of dimensions at run time.

    if (A.extent(0) != X.extent(0) || A.extent(1) != X.extent(0) || A.extent(0) != Y.extent(0)) {
      Kokkos::printf(
          "KokkosBatched::gesv: dimensions of A and X and Y do not match: A: "
          "%d x %d, X: %d, Y: %d\n",
          (int)A.extent(0), (int)A.extent(1), (int)X.extent(0), (int)Y.extent(0));
      return 1;
    }
#endif

    int r_val = SerialLU<Algo::Level3::Unblocked>::invoke(A);

    if (r_val == 0) r_val = SerialCopy<Trans::NoTranspose, 1>::invoke(Y, X);

    if (r_val == 0)
      r_val = SerialTrsm<Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::Unit, Algo::Level3::Unblocked>::invoke(
          1.0, A, X);

    if (r_val == 0)
      r_val = SerialTrsm<Side::Left, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit, Algo::Level3::Unblocked>::invoke(
          1.0, A, X);

    return r_val;
  }
};

///
/// Team Impl
/// =========

template <typename MemberType>
struct TeamGesv<MemberType, Gesv::StaticPivoting> {
  template <typename MatrixType, typename VectorType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const MatrixType A, const VectorType X,
                                           const VectorType Y) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
    static_assert(Kokkos::is_view<MatrixType>::value, "KokkosBatched::gesv: MatrixType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<VectorType>::value, "KokkosBatched::gesv: VectorType is not a Kokkos::View.");
    static_assert(MatrixType::rank == 2, "KokkosBatched::gesv: MatrixType must have rank 2.");
    static_assert(VectorType::rank == 1, "KokkosBatched::gesv: VectorType must have rank 1.");

    // Check compatibility of dimensions at run time.
    if (A.extent(0) != X.extent(0) || A.extent(1) != X.extent(0) || A.extent(0) != Y.extent(0)) {
      Kokkos::printf(
          "KokkosBatched::gesv: dimensions of A and X and Y do not match: A: "
          "%d x %d, X: %d, Y: %d\n",
          (int)A.extent(0), (int)A.extent(1), (int)X.extent(0), (int)Y.extent(0));
      return 1;
    }
#endif
    using ScratchPadMatrixViewType = Kokkos::View<typename MatrixType::non_const_value_type **,
                                                  typename MatrixType::execution_space::scratch_memory_space>;

    const int n = A.extent(0);

    ScratchPadMatrixViewType tmp(member.team_scratch(0), n, n + 4);
    auto PDAD    = Kokkos::subview(tmp, Kokkos::ALL, Kokkos::make_pair(0, n));
    auto PDY     = Kokkos::subview(tmp, Kokkos::ALL, n);
    auto D2      = Kokkos::subview(tmp, Kokkos::ALL, n + 1);
    auto tmp_v_1 = Kokkos::subview(tmp, Kokkos::ALL, n + 2);
    auto tmp_v_2 = Kokkos::subview(tmp, Kokkos::ALL, n + 3);

    if (TeamStaticPivoting<MemberType>::invoke(member, A, PDAD, Y, PDY, D2, tmp_v_1, tmp_v_2) == 1) {
      Kokkos::printf(
          "KokkosBatched::gesv: the currently implemented static pivoting "
          "failed.\n");
      return 1;
    }
    member.team_barrier();

    int r_val = TeamLU<MemberType, Algo::Level3::Unblocked>::invoke(member, PDAD);
    member.team_barrier();

    if (r_val == 0) {
      r_val = TeamTrsm<MemberType, Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::Unit,
                       Algo::Level3::Unblocked>::invoke(member, 1.0, PDAD, PDY);
      member.team_barrier();
    }

    if (r_val == 0) {
      r_val = TeamTrsm<MemberType, Side::Left, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit,
                       Algo::Level3::Unblocked>::invoke(member, 1.0, PDAD, PDY);
      member.team_barrier();
    }

    if (r_val == 0) {
      TeamHadamard1D(member, PDY, D2, X);
      member.team_barrier();
    }

    return r_val;
  }
};

template <typename MemberType>
struct TeamGesv<MemberType, Gesv::NoPivoting> {
  template <typename MatrixType, typename VectorType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const MatrixType A, const VectorType X,
                                           const VectorType Y) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
    static_assert(Kokkos::is_view<MatrixType>::value, "KokkosBatched::gesv: MatrixType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<VectorType>::value, "KokkosBatched::gesv: VectorType is not a Kokkos::View.");
    static_assert(MatrixType::rank == 2, "KokkosBatched::gesv: MatrixType must have rank 2.");
    static_assert(VectorType::rank == 1, "KokkosBatched::gesv: VectorType must have rank 1.");

    // Check compatibility of dimensions at run time.
    if (A.extent(0) != X.extent(0) || A.extent(1) != X.extent(0) || A.extent(0) != Y.extent(0)) {
      Kokkos::printf(
          "KokkosBatched::gesv: dimensions of A and X and Y do not match: A: "
          "%d x %d, X: %d, Y: %d\n",
          (int)A.extent(0), (int)A.extent(1), (int)X.extent(0), (int)Y.extent(0));
      return 1;
    }
#endif

    int r_val = TeamLU<MemberType, Algo::Level3::Unblocked>::invoke(member, A);
    member.team_barrier();

    if (r_val == 0) {
      TeamCopy<MemberType, Trans::NoTranspose, 1>::invoke(member, Y, X);
      member.team_barrier();
    }

    if (r_val == 0) {
      TeamTrsm<MemberType, Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::Unit, Algo::Level3::Unblocked>::invoke(
          member, 1.0, A, X);
      member.team_barrier();
    }

    if (r_val == 0) {
      TeamTrsm<MemberType, Side::Left, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit, Algo::Level3::Unblocked>::invoke(
          member, 1.0, A, X);
      member.team_barrier();
    }

    return r_val;
  }
};

///
/// TeamVector Impl
/// =========

template <typename MemberType>
struct TeamVectorGesv<MemberType, Gesv::StaticPivoting> {
  template <typename MatrixType, typename VectorType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const MatrixType A, const VectorType X,
                                           const VectorType Y) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
    static_assert(Kokkos::is_view<MatrixType>::value, "KokkosBatched::gesv: MatrixType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<VectorType>::value, "KokkosBatched::gesv: VectorType is not a Kokkos::View.");
    static_assert(MatrixType::rank == 2, "KokkosBatched::gesv: MatrixType must have rank 2.");
    static_assert(VectorType::rank == 1, "KokkosBatched::gesv: VectorType must have rank 1.");

    // Check compatibility of dimensions at run time.
    if (A.extent(0) != X.extent(0) || A.extent(1) != X.extent(0) || A.extent(0) != Y.extent(0)) {
      Kokkos::printf(
          "KokkosBatched::gesv: dimensions of A and X and Y do not match: A: "
          "%d x %d, X: %d, Y: %d\n",
          (int)A.extent(0), (int)A.extent(1), (int)X.extent(0), (int)Y.extent(0));
      return 1;
    }
#endif
    using ScratchPadMatrixViewType = Kokkos::View<typename MatrixType::non_const_value_type **,
                                                  typename MatrixType::execution_space::scratch_memory_space>;

    const int n = A.extent(0);

    ScratchPadMatrixViewType tmp(member.team_scratch(0), n, n + 4);
    auto PDAD    = Kokkos::subview(tmp, Kokkos::ALL, Kokkos::make_pair(0, n));
    auto PDY     = Kokkos::subview(tmp, Kokkos::ALL, n);
    auto D2      = Kokkos::subview(tmp, Kokkos::ALL, n + 1);
    auto tmp_v_1 = Kokkos::subview(tmp, Kokkos::ALL, n + 2);
    auto tmp_v_2 = Kokkos::subview(tmp, Kokkos::ALL, n + 3);

    if (TeamVectorStaticPivoting<MemberType>::invoke(member, A, PDAD, Y, PDY, D2, tmp_v_1, tmp_v_2) == 1) {
      Kokkos::printf(
          "KokkosBatched::gesv: the currently implemented static pivoting "
          "failed.\n");
      return 1;
    }

    member.team_barrier();

    int r_val = TeamLU<MemberType, Algo::Level3::Unblocked>::invoke(member, PDAD);
    member.team_barrier();

    if (r_val == 0) {
      TeamVectorTrsm<MemberType, Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::Unit,
                     Algo::Level3::Unblocked>::invoke(member, 1.0, PDAD, PDY);
      member.team_barrier();
    }

    if (r_val == 0) {
      TeamVectorTrsm<MemberType, Side::Left, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit,
                     Algo::Level3::Unblocked>::invoke(member, 1.0, PDAD, PDY);
      member.team_barrier();
    }

    if (r_val == 0) {
      TeamVectorHadamard1D(member, PDY, D2, X);
      member.team_barrier();
    }

    return r_val;
  }
};

template <typename MemberType>
struct TeamVectorGesv<MemberType, Gesv::NoPivoting> {
  template <typename MatrixType, typename VectorType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const MatrixType A, const VectorType X,
                                           const VectorType Y) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
    static_assert(Kokkos::is_view<MatrixType>::value, "KokkosBatched::gesv: MatrixType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<VectorType>::value, "KokkosBatched::gesv: VectorType is not a Kokkos::View.");
    static_assert(MatrixType::rank == 2, "KokkosBatched::gesv: MatrixType must have rank 2.");
    static_assert(VectorType::rank == 1, "KokkosBatched::gesv: VectorType must have rank 1.");

    // Check compatibility of dimensions at run time.
    if (A.extent(0) != X.extent(0) || A.extent(1) != X.extent(0) || A.extent(0) != Y.extent(0)) {
      Kokkos::printf(
          "KokkosBatched::gesv: dimensions of A and X and Y do not match: A: "
          "%d x %d, X: %d, Y: %d\n",
          (int)A.extent(0), (int)A.extent(1), (int)X.extent(0), (int)Y.extent(0));
      return 1;
    }
#endif

    int r_val = TeamLU<MemberType, Algo::Level3::Unblocked>::invoke(member, A);
    member.team_barrier();

    if (r_val == 0) {
      TeamVectorCopy<MemberType, Trans::NoTranspose, 1>::invoke(member, Y, X);
      member.team_barrier();
    }

    if (r_val == 0) {
      TeamVectorTrsm<MemberType, Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::Unit,
                     Algo::Level3::Unblocked>::invoke(member, 1.0, A, X);
      member.team_barrier();
    }

    if (r_val == 0) {
      TeamVectorTrsm<MemberType, Side::Left, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit,
                     Algo::Level3::Unblocked>::invoke(member, 1.0, A, X);
      member.team_barrier();
    }

    return r_val;
  }
};

}  // namespace KokkosBatched

#endif
