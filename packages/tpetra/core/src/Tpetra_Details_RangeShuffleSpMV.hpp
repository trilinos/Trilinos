// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_RANGESHUFFLESPMV_HPP
#define TPETRA_DETAILS_RANGESHUFFLESPMV_HPP

namespace Tpetra::Details {

#if defined(KOKKOS_ENABLE_CUDA)
template <class Matrix,
          class XVector,
          class YVector,
          class Scalar>
struct SpMVRangeShuffleFunctor {
  static constexpr int group_size      = 4;
  static constexpr int nrows_per_group = 2;
  static constexpr int unroll_len      = 8;
  static constexpr int group_shift     = 2;
  static constexpr int row_shift       = 1;

  using execution_space = typename Matrix::execution_space;
  using LO              = typename Matrix::non_const_ordinal_type;
  using offset_type     = typename Matrix::non_const_size_type;

  const Scalar alpha;
  Matrix A;
  XVector X;
  const Scalar beta;
  YVector Y;

  SpMVRangeShuffleFunctor(const char op_[],
                          const Scalar& alpha_,
                          const Matrix& A_,
                          const XVector& X_,
                          const Scalar& beta_,
                          const YVector& Y_)
    : alpha(alpha_)
    , A(A_)
    , X(X_)
    , beta(beta_)
    , Y(Y_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const LO id) const {
    using residual_value_type = typename YVector::non_const_value_type;
    using KAT                 = KokkosKernels::ArithTraits<residual_value_type>;

    const int lane     = static_cast<int>(id & static_cast<LO>(group_size - 1));
    const LO base_row  = static_cast<LO>((id >> group_shift) << row_shift);
    const LO numRows   = static_cast<LO>(A.numRows());
    const auto row_map = A.graph.row_map.data();
    const auto entries = A.graph.entries.data();
    const auto values  = A.values.data();

    residual_value_type sum[nrows_per_group];
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
    for (int i = 0; i < nrows_per_group; ++i) {
      sum[i] = KAT::zero();
    }

#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
    for (int i = 0; i < nrows_per_group; ++i) {
      const LO row = base_row + static_cast<LO>(i);
      if (row < numRows) {
        const offset_type beg      = static_cast<offset_type>(row_map[row]);
        const offset_type end      = static_cast<offset_type>(row_map[row + 1]);
        constexpr offset_type step = group_size * unroll_len;
        for (offset_type nz = beg + static_cast<offset_type>(lane); nz < end; nz += step) {
          sum[i] += values[nz] * X(entries[nz]);
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
          for (int u = 1; u < unroll_len; ++u) {
            const offset_type nzu = nz + static_cast<offset_type>(u * group_size);
            if (nzu < end) {
              sum[i] += values[nzu] * X(entries[nzu]);
            }
          }
        }
      }
    }

#if defined(__CUDA_ARCH__)
#pragma unroll
    for (int i = 0; i < nrows_per_group; ++i) {
      sum[i] += Kokkos::shfl_down(sum[i], 2, 32);
      sum[i] += Kokkos::shfl_down(sum[i], 1, 32);
    }
#endif

    if (lane == 0) {
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
      for (int i = 0; i < nrows_per_group; ++i) {
        const LO row = base_row + static_cast<LO>(i);
        if (row < numRows) {
          Y(row) = beta * Y(row) + alpha * sum[i];
        }
      }
    }
  }
};

template <class Matrix,
          class XVector,
          class YVector,
          class Scalar>
static void
spmv_range_shuffle(const char op[],
                   const Scalar& alpha,
                   const Matrix& A,
                   const XVector& X,
                   const Scalar& beta,
                   const YVector& Y) {
  using execution_space = typename Matrix::execution_space;
  using LO              = typename Matrix::non_const_ordinal_type;
  using functor_type    = SpMVRangeShuffleFunctor<Matrix, XVector, YVector, Scalar>;

  constexpr int nrows_per_group = functor_type::nrows_per_group;
  constexpr int row_shift       = functor_type::row_shift;
  constexpr int group_shift     = functor_type::group_shift;
  const LO numRows              = static_cast<LO>(A.numRows());
  const LO ngroups              = static_cast<LO>((numRows + nrows_per_group - 1) >> row_shift);
  const LO nthreads             = static_cast<LO>(ngroups << group_shift);

  functor_type func(op, alpha, A, X, beta, Y);
  Kokkos::parallel_for("spmv_range_shuffle_u8_r2",
                       Kokkos::RangePolicy<execution_space, Kokkos::IndexType<LO> >(0, nthreads),
                       func);
}
#endif

}  // namespace Tpetra::Details

#endif
