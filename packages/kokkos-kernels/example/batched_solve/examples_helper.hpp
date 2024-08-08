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

/// \brief create_saddle_point_matrices:
///
///  This function creates the matrices and the rhs of a batched saddle point
///  systems where A and Y (the right hand side) are as follows:
///
///        ___________
///       |     |   T |
///       |  B  |  C  |
///  A =  |-----+-----|
///       |  C  |  0  |
///       |_____|_____|
///
///        _____
///       |     |
///       |  D  |
///  Y =  |-----|
///       |  0  |
///       |_____|
///
///  with A in R^{n \times n}, B in R^{(n-n_2) \times (n-n_2)} and
///  where B and C are computed as follows:
///
///  1. A sequence of n-n_2 points of R^{n_dim} is generated randomly:
///     x^(0), ..., x^(n-n_2-1)
///  2. Given this sequence, the entries are computed as follows:
///     B_{(i,j)} = \| x^(i) - x^(j)\|
///     C_{(0,j)} = 1
///     C_{(i,j)} = (x^(j))_{(i-1)} for i != 0
///
///  3. D is generated randomly.
///
/// This function uses a different sequence of x and a different D for every
/// systems within the batched system.
///
/// As a consequence of its definitation, the diagonal of A is 0 for every
/// entries.
///
/// \tparam MatrixViewType: type of the batched matrices
/// \tparam VectorViewType: type of the batched vectors
///
/// \param A [in/out]: a rank 3 view that has to be prealocated that will store
/// the entries of the batched matrix. \param Y [in/out]: a rank 2 view that has
/// to be prealocated that will store the entries of the right hand side. \param
/// n_dim [in]: the dimension of the physical space where the points are
/// randomly generated (default = 3).
///

template <typename MatrixViewType, typename VectorViewType>
void create_saddle_point_matrices(const MatrixViewType &A, const VectorViewType &Y, const int n_dim = 3) {
  Kokkos::Random_XorShift64_Pool<typename MatrixViewType::device_type::execution_space> random(13718);
  const int N   = A.extent(0);
  const int n   = A.extent(1);
  const int n_2 = n_dim + 1;
  const int n_1 = n - n_2;

  MatrixViewType xs("xs", N, n_1, n_dim);
  VectorViewType ys("ys", N, n_1);

  Kokkos::fill_random(xs, random, Kokkos::reduction_identity<typename MatrixViewType::value_type>::prod());
  Kokkos::fill_random(ys, random, Kokkos::reduction_identity<typename VectorViewType::value_type>::prod());

  auto xs_host = Kokkos::create_mirror_view(xs);
  auto ys_host = Kokkos::create_mirror_view(ys);
  auto A_host  = Kokkos::create_mirror_view(A);
  auto Y_host  = Kokkos::create_mirror_view(Y);

  Kokkos::deep_copy(xs_host, xs);
  Kokkos::deep_copy(ys_host, ys);

  for (int i = 0; i < n_1; ++i) {
    for (int j = 0; j < n_1; ++j) {
      for (int l = 0; l < N; ++l) {
        auto xs_i                             = Kokkos::subview(xs_host, l, i, Kokkos::ALL);
        auto xs_j                             = Kokkos::subview(xs_host, l, j, Kokkos::ALL);
        typename MatrixViewType::value_type d = 0;
        for (int k = 0; k < n_dim; ++k) d += Kokkos::pow(xs_i(k) - xs_j(k), 2);
        d               = Kokkos::sqrt(d);
        A_host(l, i, j) = Kokkos::pow(d, 5);
      }
    }
    for (int l = 0; l < N; ++l) {
      A_host(l, i, n_1) = (typename MatrixViewType::value_type)1.0;
      A_host(l, n_1, i) = (typename MatrixViewType::value_type)1.0;
      for (int k = 0; k < n_dim; ++k) {
        A_host(l, i, n_1 + k + 1) = xs_host(l, i, k);
        A_host(l, n_1 + k + 1, i) = xs_host(l, i, k);
      }
      Y_host(l, i) = ys_host(l, i);
    }
  }
  for (int i = n_1; i < n; ++i) {
    for (int l = 0; l < N; ++l) {
      Y_host(l, i) = (typename MatrixViewType::value_type)0.0;
    }
  }

  Kokkos::deep_copy(A, A_host);
  Kokkos::deep_copy(Y, Y_host);

  Kokkos::fence();
}

template <typename IntView, typename VectorViewType>
void create_tridiagonal_batched_matrices(const int nnz, const int BlkSize, const int N, const IntView &r,
                                         const IntView &c, const VectorViewType &D, const VectorViewType &X,
                                         const VectorViewType &B) {
  Kokkos::Random_XorShift64_Pool<typename VectorViewType::device_type::execution_space> random(13718);
  Kokkos::fill_random(X, random, Kokkos::reduction_identity<typename VectorViewType::value_type>::prod());
  Kokkos::fill_random(B, random, Kokkos::reduction_identity<typename VectorViewType::value_type>::prod());

  auto D_host = Kokkos::create_mirror_view(D);
  auto r_host = Kokkos::create_mirror_view(r);
  auto c_host = Kokkos::create_mirror_view(c);

  r_host(0) = 0;

  int current_col = 0;

  for (int i = 0; i < BlkSize; ++i) {
    r_host(i + 1) = r_host(i) + (i == 0 || i == (BlkSize - 1) ? 2 : 3);
  }
  for (int i = 0; i < nnz; ++i) {
    if (i % 3 == 0) {
      for (int l = 0; l < N; ++l) {
        D_host(l, i) = typename VectorViewType::value_type(2.0);
      }
      c_host(i) = current_col;
      ++current_col;
    } else {
      for (int l = 0; l < N; ++l) {
        D_host(l, i) = typename VectorViewType::value_type(-1.0);
      }
      c_host(i) = current_col;
      if (i % 3 == 1)
        --current_col;
      else
        ++current_col;
    }
  }

  Kokkos::fence();

  Kokkos::deep_copy(D, D_host);
  Kokkos::deep_copy(r, r_host);
  Kokkos::deep_copy(c, c_host);

  Kokkos::fence();
}

template <class VType, class IntType>
void getInvDiagFromCRS(const VType &V, const IntType &r, const IntType &c, const VType &diag) {
  auto diag_values_host = Kokkos::create_mirror_view(diag);
  auto values_host      = Kokkos::create_mirror_view(V);
  auto row_ptr_host     = Kokkos::create_mirror_view(r);
  auto colIndices_host  = Kokkos::create_mirror_view(c);

  Kokkos::deep_copy(values_host, V);
  Kokkos::deep_copy(row_ptr_host, r);
  Kokkos::deep_copy(colIndices_host, c);

  int current_index;
  int N       = diag.extent(0);
  int BlkSize = diag.extent(1);

  for (int i = 0; i < BlkSize; ++i) {
    for (current_index = row_ptr_host(i); current_index < row_ptr_host(i + 1); ++current_index) {
      if (colIndices_host(current_index) == i) break;
    }
    for (int j = 0; j < N; ++j) {
      diag_values_host(j, i) = 1. / values_host(j, current_index);
    }
  }

  Kokkos::deep_copy(diag, diag_values_host);
}