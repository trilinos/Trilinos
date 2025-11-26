// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef TEST_BATCHED_DENSE_HELPER_HPP
#define TEST_BATCHED_DENSE_HELPER_HPP

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {
template <typename MatrixViewType, typename VectorViewType>
void create_tridiagonal_batched_matrices(const MatrixViewType& A, const VectorViewType& B) {
  Kokkos::Random_XorShift64_Pool<typename VectorViewType::device_type::execution_space> random(13718);
  Kokkos::fill_random(B, random, 1);

  auto A_host = Kokkos::create_mirror_view(A);

  const int N       = A.extent(0);
  const int BlkSize = A.extent(1);

  for (int l = 0; l < N; ++l) {
    for (int i = 0; i < BlkSize; ++i) {
      for (int j = i; j < BlkSize; ++j) {
        if (i == j)
          A_host(l, i, j) = typename VectorViewType::value_type(2.0);
        else if (i == j - 1) {
          A_host(l, i, j) = typename VectorViewType::value_type(-1.0);
          A_host(l, j, i) = typename VectorViewType::value_type(-1.0);
        } else {
          A_host(l, i, j) = typename VectorViewType::value_type(0.0);
          A_host(l, j, i) = typename VectorViewType::value_type(0.0);
        }
      }
    }
  }

  Kokkos::fence();

  Kokkos::deep_copy(A, A_host);

  Kokkos::fence();
}

template <typename InViewType, typename OutViewType, typename UploType>
void create_banded_triangular_matrix(InViewType& in, OutViewType& out, int k = 1, bool band_storage = true) {
  auto h_in   = Kokkos::create_mirror_view(in);
  auto h_out  = Kokkos::create_mirror_view(out);
  const int N = in.extent(0), BlkSize = in.extent(1);

  Kokkos::deep_copy(h_in, in);
  if (band_storage) {
    assert(out.extent(0) == in.extent(0));
    assert(out.extent(1) == static_cast<std::size_t>(k + 1));
    assert(out.extent(2) == in.extent(2));
    if constexpr (std::is_same_v<UploType, KokkosBatched::Uplo::Upper>) {
      for (int i0 = 0; i0 < N; i0++) {
        for (int i1 = 0; i1 < k + 1; i1++) {
          for (int i2 = i1; i2 < BlkSize; i2++) {
            h_out(i0, k - i1, i2) = h_in(i0, i2 - i1, i2);
          }
        }
      }
    } else {
      for (int i0 = 0; i0 < N; i0++) {
        for (int i1 = 0; i1 < k + 1; i1++) {
          for (int i2 = 0; i2 < BlkSize - i1; i2++) {
            h_out(i0, i1, i2) = h_in(i0, i2 + i1, i2);
          }
        }
      }
    }
  } else {
    for (std::size_t i = 0; i < InViewType::rank(); i++) {
      assert(out.extent(i) == in.extent(i));
    }

    if constexpr (std::is_same_v<UploType, KokkosBatched::Uplo::Upper>) {
      for (int i0 = 0; i0 < N; i0++) {
        for (int i1 = 0; i1 < BlkSize; i1++) {
          for (int i2 = i1; i2 < Kokkos::min(i1 + k + 1, BlkSize); i2++) {
            h_out(i0, i1, i2) = h_in(i0, i1, i2);
          }
        }
      }
    } else {
      for (int i0 = 0; i0 < N; i0++) {
        for (int i1 = 0; i1 < BlkSize; i1++) {
          for (int i2 = Kokkos::max(0, i1 - k); i2 <= i1; i2++) {
            h_out(i0, i1, i2) = h_in(i0, i1, i2);
          }
        }
      }
    }
  }
  Kokkos::deep_copy(out, h_out);
}

/// \brief Create a diagonal matrix from an input vector:
/// Copies the input vector into the diagonal of the output matrix specified
/// by the parameter k. k > 0 means that the matrix is upper-diagonal and
/// k < 0 means the lower-diagonal. k = 0 means the diagonal.
///
/// \tparam InViewType: Input type for the vector, needs to be a 2D view
/// \tparam OutViewType: Output type for the matrix, needs to be a 3D view
///
/// \param in [in]: Input batched vector, a rank 2 view
/// \param out [out]: Output batched matrix, where the diagonal compnent
/// specified by k is filled with the input vector, a rank 3 view
/// \param k [in]: The diagonal offset to be filled (default is 0).
///
template <typename InViewType, typename OutViewType>
void create_diagonal_matrix(InViewType& in, OutViewType& out, int k = 0) {
  auto h_in   = Kokkos::create_mirror_view(in);
  auto h_out  = Kokkos::create_mirror_view(out);
  const int N = in.extent(0), BlkSize = in.extent(1);

  assert(out.extent(0) == in.extent(0));
  assert(out.extent(1) == in.extent(1) + abs(k));

  int i1_start = k >= 0 ? 0 : -k;
  int i2_start = k >= 0 ? k : 0;

  // Zero clear the output matrix
  using ScalarType = typename OutViewType::non_const_value_type;
  Kokkos::deep_copy(h_out, ScalarType(0.0));

  Kokkos::deep_copy(h_in, in);
  for (int i0 = 0; i0 < N; i0++) {
    for (int i1 = 0; i1 < BlkSize; i1++) {
      h_out(i0, i1 + i1_start, i1 + i2_start) = h_in(i0, i1);
    }
  }

  Kokkos::deep_copy(out, h_out);
}

/// \brief Creates a positive definite symmetric (PDS) matrix.
/// Takes a dense random matrix and converts it to a dense pds matrix.
///
/// \tparam InViewType: Input type for the matrix, needs to be a 3D view
/// \tparam OutViewType: Output type for the matrix, needs to be a 3D view
///
/// \param in [in]: Input batched dense matrix, a rank 3 view
/// \param out [out]: Output batched dense matrix, a rank 3 view
///
template <typename InViewType, typename OutViewType>
void random_to_pds(InViewType& in, OutViewType& out) {
  auto h_in   = Kokkos::create_mirror_view(in);
  auto h_out  = Kokkos::create_mirror_view(out);
  const int N = in.extent(0), BlkSize = in.extent(1);
  using value_type = typename InViewType::non_const_value_type;

  for (std::size_t i = 0; i < InViewType::rank(); i++) {
    assert(out.extent(i) == in.extent(i));
  }

  Kokkos::deep_copy(h_in, in);
  Kokkos::deep_copy(h_out, 0.0);

  for (int i0 = 0; i0 < N; i0++) {
    // Make a hermitian matrix
    for (int i1 = 0; i1 < BlkSize; i1++) {
      for (int i2 = i1; i2 < BlkSize; i2++) {
        if (i1 == i2) {
          // Diagonal elements must be real
          h_out(i0, i1, i2) = KokkosKernels::ArithTraits<value_type>::real(h_in(i0, i1, i2));
        } else {
          // Off-diagonal elements are complex and Hermitian
          h_out(i0, i1, i2) = h_in(i0, i1, i2);
          h_out(i0, i2, i1) = KokkosKernels::ArithTraits<value_type>::conj(h_in(i0, i1, i2));
        }
      }
    }

    // Make matrix diagonal dominant
    for (int i1 = 0; i1 < BlkSize; i1++) {
      value_type sum = 0;
      for (int i2 = 0; i2 < BlkSize; i2++) {
        if (i1 != i2) {
          sum += Kokkos::abs(h_out(i0, i1, i2));
        }
      }
      h_out(i0, i1, i1) = sum + 1.0;
    }
  }
  Kokkos::deep_copy(out, h_out);
}

/// \brief Creates a banded positive definite symmetric (PDS) matrix.
/// Takes a dense diagonal dominant matrix and converts it to a banded pds matrix either
/// in banded or conventional storage.
///
/// \tparam InViewType: Input type for the matrix, needs to be a 3D view
/// \tparam OutViewType: Output type for the matrix, needs to be a 3D view
/// \tparam UploType: Type indicating whether the matrix is upper or lower triangular
///
/// \param in [in]: Input batched banded matrix, a rank 3 view
/// \param out [out]: Output batched dense matrix, a rank 3 view
/// \param k [in]: Number of sub/super-diagonals for lower/upper triangular (default is 1)
/// \param band_storage [in]: Boolean flag indicating whether the output should be in banded storage format (default is
/// true)
template <typename InViewType, typename OutViewType, typename UploType>
void create_banded_pds_matrix(InViewType& in, OutViewType& out, int k = 1, bool band_storage = true) {
  auto h_in        = Kokkos::create_mirror_view(in);
  auto h_out       = Kokkos::create_mirror_view(out);
  using value_type = typename InViewType::non_const_value_type;
  const int N = in.extent(0), BlkSize = in.extent(1);

  Kokkos::deep_copy(h_in, in);

  if (band_storage) {
    assert(out.extent(0) == in.extent(0));
    assert(out.extent(1) == static_cast<std::size_t>(k + 1));
    assert(out.extent(2) == in.extent(2));
    if constexpr (std::is_same_v<UploType, KokkosBatched::Uplo::Upper>) {
      for (int i0 = 0; i0 < N; i0++) {
        for (int i1 = 0; i1 < k + 1; i1++) {
          for (int i2 = i1; i2 < BlkSize; i2++) {
            h_out(i0, k - i1, i2) = h_in(i0, i2 - i1, i2);
          }
        }
      }
    } else {
      for (int i0 = 0; i0 < N; i0++) {
        for (int i1 = 0; i1 < k + 1; i1++) {
          for (int i2 = 0; i2 < BlkSize - i1; i2++) {
            h_out(i0, i1, i2) = h_in(i0, i2 + i1, i2);
          }
        }
      }
    }
  } else {
    for (std::size_t i = 0; i < InViewType::rank(); i++) {
      assert(out.extent(i) == in.extent(i));
    }

    if constexpr (std::is_same_v<UploType, KokkosBatched::Uplo::Upper>) {
      for (int i0 = 0; i0 < N; i0++) {
        for (int i1 = 0; i1 < BlkSize; i1++) {
          for (int i2 = i1; i2 < Kokkos::min(i1 + k + 1, BlkSize); i2++) {
            h_out(i0, i1, i2) = h_in(i0, i1, i2);
            h_out(i0, i2, i1) = KokkosKernels::ArithTraits<value_type>::conj(h_in(i0, i1, i2));
          }
        }
      }
    } else {
      for (int i0 = 0; i0 < N; i0++) {
        for (int i1 = 0; i1 < BlkSize; i1++) {
          for (int i2 = Kokkos::max(0, i1 - k); i2 <= i1; i2++) {
            h_out(i0, i1, i2) = h_in(i0, i1, i2);
            h_out(i0, i2, i1) = KokkosKernels::ArithTraits<value_type>::conj(h_in(i0, i1, i2));
          }
        }
      }
    }
  }
  Kokkos::deep_copy(out, h_out);
}

/// \brief Converts a tridiagonal matrix from band storage to conventional storage.
/// Takes a upper/lower triangular banded matrix in banded storage
/// and converts it to a conventional storage
///
/// \tparam InViewType: Input type for the matrix, needs to be a 3D view
/// \tparam OutViewType: Output type for the matrix, needs to be a 3D view
/// \tparam UploType: Type indicating whether the matrix is upper or lower triangular
///
/// \param in [in]: Input batched banded matrix, a rank 3 view
/// \param out [out]: Output batched dense matrix, a rank 3 view
/// \param k [in]: Number of sub/super-diagonals for lower/upper triangular (default is 1)
///
template <typename InViewType, typename OutViewType, typename UploType>
void banded_to_dense(InViewType& in, OutViewType& out, int k = 1) {
  auto h_in   = Kokkos::create_mirror_view(in);
  auto h_out  = Kokkos::create_mirror_view(out);
  const int N = in.extent(0), BlkSize = in.extent(2);

  Kokkos::deep_copy(h_in, in);
  assert(in.extent(0) == out.extent(0));
  assert(in.extent(1) == static_cast<std::size_t>(k + 1));
  assert(in.extent(2) == out.extent(2));
  if constexpr (std::is_same_v<UploType, KokkosBatched::Uplo::Upper>) {
    for (int i0 = 0; i0 < N; i0++) {
      for (int i1 = 0; i1 < k + 1; i1++) {
        for (int i2 = i1; i2 < BlkSize; i2++) {
          h_out(i0, i2 - i1, i2) = h_in(i0, k - i1, i2);
        }
      }
    }
  } else {
    for (int i0 = 0; i0 < N; i0++) {
      for (int i1 = 0; i1 < k + 1; i1++) {
        for (int i2 = 0; i2 < BlkSize - i1; i2++) {
          h_out(i0, i2 + i1, i2) = h_in(i0, i1, i2);
        }
      }
    }
  }
  Kokkos::deep_copy(out, h_out);
}

/// \brief Converts a matrix from band storage to conventional storage.
/// Takes a banded matrix in banded storage and converts it to a conventional storage.
///
/// \tparam InViewType: Input type for the matrix, needs to be a 3D view
/// \tparam OutViewType: Output type for the matrix, needs to be a 3D view
///
/// \param in [in]: Input batched banded matrix, a rank 3 view
/// \param out [out]: Output batched dense matrix, a rank 3 view
/// \param kl [in]: Number of subdiagonals for within the band of A (default is 1)
/// \param ku [in]: Number of superdiagonals for within the band of A (default is 1)
///
template <typename InViewType, typename OutViewType>
void banded_to_dense(InViewType& in, OutViewType& out, int kl = 1, int ku = 1) {
  auto h_in    = Kokkos::create_mirror_view(in);
  auto h_out   = Kokkos::create_mirror_view(out);
  const int Nb = out.extent(0), m = out.extent(1), n = out.extent(2);

  Kokkos::deep_copy(h_in, in);
  assert(in.extent(0) == out.extent(0));
  assert(in.extent(1) == static_cast<std::size_t>(2 * kl + ku + 1));
  assert(in.extent(2) == out.extent(2));

  for (int i0 = 0; i0 < Nb; i0++) {
    for (int i1 = 0; i1 < m; i1++) {
      for (int i2 = Kokkos::max(0, i1 - ku); i2 < Kokkos::min(n, i1 + kl + ku + 1); i2++) {
        auto row_in_banded = kl + ku + i1 - i2;
        h_out(i0, i1, i2)  = h_in(i0, row_in_banded, i2);
      }
    }
  }
  Kokkos::deep_copy(out, h_out);
}

/// \brief Converts a matrix from conventional storage to band storage
/// Takes a band matrix in conventional storage and converts it to a banded storage.
///
/// \tparam InViewType: Input type for the matrix, needs to be a 3D view
/// \tparam OutViewType: Output type for the matrix, needs to be a 3D view
///
/// \param in [in]: Input batched dense matrix, a rank 3 view
/// \param out [out]: Output batched banded matrix, a rank 3 view
/// \param kl [in]: Number of subdiagonals for within the band of A (default is 1)
/// \param ku [in]: Number of superdiagonals for within the band of A (default is 1)
///
template <typename InViewType, typename OutViewType>
void dense_to_banded(InViewType& in, OutViewType& out, int kl = 1, int ku = 1) {
  auto h_in    = Kokkos::create_mirror_view(in);
  auto h_out   = Kokkos::create_mirror_view(out);
  const int Nb = in.extent(0), m = in.extent(1), n = in.extent(2);

  Kokkos::deep_copy(h_in, in);
  assert(out.extent(0) == in.extent(0));
  assert(out.extent(1) == static_cast<std::size_t>(2 * kl + ku + 1));
  assert(out.extent(2) == in.extent(2));

  for (int i0 = 0; i0 < Nb; i0++) {
    for (int i1 = 0; i1 < m; i1++) {
      for (int i2 = Kokkos::max(0, i1 - ku); i2 < Kokkos::min(n, i1 + kl + 1); i2++) {
        // Do not use the element from the dense matrix if it is outside the band
        auto row_in_banded           = kl + ku + i1 - i2;
        h_out(i0, row_in_banded, i2) = h_in(i0, i1, i2);
      }
    }
  }
  Kokkos::deep_copy(out, h_out);
}

/// \brief Create a triangular matrix from an input matrix:
/// Copies the input matrix into the upper/lower triangular of the output matrix specified
/// by the parameter k. Zero out elements below/above the k-th diagonal.
///
/// \tparam InViewType: Input type for the matrix, needs to be a 3D view
/// \tparam OutViewType: Output type for the matrix, needs to be a 3D view
/// \tparam UploType: Type indicating whether the matrix is upper or lower triangular
/// \tparam DiagType: Type indicating whether the matrix is unit or non-unit diagonal
///
/// \param in [in]: Input batched matrix, a rank 3 view
/// \param out [out]: Output batched matrix, where the upper or lower
/// triangular components are kept, a rank 3 view
/// \param k [in]: The diagonal offset to be zero out (default is 0).
///
template <typename InViewType, typename OutViewType, typename UploType, typename DiagType>
void create_triangular_matrix(InViewType& in, OutViewType& out, int k = 0) {
  auto h_in    = Kokkos::create_mirror_view(in);
  auto h_out   = Kokkos::create_mirror_view(out);
  const int Nb = in.extent(0), m = in.extent(1), n = in.extent(2);

  Kokkos::deep_copy(h_in, in);
  Kokkos::deep_copy(h_out, 0.0);
  for (int i0 = 0; i0 < Nb; i0++) {
    for (int i1 = 0; i1 < m; i1++) {
      for (int i2 = 0; i2 < n; i2++) {
        if constexpr (std::is_same_v<UploType, KokkosBatched::Uplo::Upper>) {
          // Upper
          // Zero out elements below the k-th diagonal
          h_out(i0, i1, i2) = i2 < i1 + k ? 0.0 : h_in(i0, i1, i2);
        } else {
          // Lower
          // Zero out elements above the k-th diagonal
          h_out(i0, i1, i2) = i2 > i1 + k ? 0.0 : h_in(i0, i1, i2);
        }
      }
    }
    for (int i1 = 0; i1 < Kokkos::min(m, n); i1++) {
      if constexpr (std::is_same_v<DiagType, KokkosBatched::Diag::Unit>) {
        // Unit diagonal
        h_out(i0, i1, i1) = 1.0;
      }
    }
  }
  Kokkos::deep_copy(out, h_out);
}

}  // namespace KokkosBatched

#endif  // TEST_BATCHED_DENSE_HELPER_HPP
