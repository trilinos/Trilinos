// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_GEMM_SERIAL_IMPL_HPP
#define KOKKOSBATCHED_GEMM_SERIAL_IMPL_HPP

#include "KokkosBlas_util.hpp"
#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Gemm_Serial_Internal.hpp"

namespace KokkosBatched {
namespace Impl {
template <typename ArgTransA, typename ArgTransB, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION static int checkGemmInput([[maybe_unused]] const AViewType &A,
                                                 [[maybe_unused]] const BViewType &B,
                                                 [[maybe_unused]] const CViewType &C) {
  static_assert(Kokkos::is_view_v<AViewType>, "KokkosBatched::gemm: AViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<BViewType>, "KokkosBatched::gemm: BViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<CViewType>, "KokkosBatched::gemm: CViewType is not a Kokkos::View.");

  static_assert(AViewType::rank <= 2, "KokkosBatched::gemm: AViewType must have rank 0, 1 or 2.");
  static_assert(BViewType::rank <= 2, "KokkosBatched::gemm: BViewType must have rank 0, 1 or 2.");
  static_assert(CViewType::rank <= 2, "KokkosBatched::gemm: CViewType must have rank 0, 1 or 2.");

#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  const int m = C.extent(0), n = C.extent(1);
  const int lda = A.extent(0);
  const int ldb = B.extent(0);

  const int ka = std::is_same_v<ArgTransA, Trans::NoTranspose> ? A.extent(1) : A.extent(0);
  const int kb = std::is_same_v<ArgTransB, Trans::NoTranspose> ? B.extent(0) : B.extent(1);

  if (ka != kb) {
    Kokkos::printf(
        "KokkosBatched::gemm: Dimensions of A and B do not match: A: %d x %d, "
        "B: %d x %d\n",
        A.extent(0), A.extent(1), B.extent(0), B.extent(1));
    return 1;
  }

  const int nrowa = std::is_same_v<ArgTransA, Trans::NoTranspose> ? m : ka;
  const int nrowb = std::is_same_v<ArgTransB, Trans::NoTranspose> ? kb : n;

  if (lda < Kokkos::max(1, nrowa)) {
    Kokkos::printf(
        "KokkosBatched::gemm: leading dimension of A must not be smaller than "
        "max(1, nrowa): "
        "lda = %d, nrowa = %d\n",
        lda, nrowa);
    return 1;
  }
  if (ldb < Kokkos::max(1, nrowb)) {
    Kokkos::printf(
        "KokkosBatched::gemm: leading dimension of B must not be smaller than "
        "max(1, nrowb): "
        "ldb = %d, nrowb = %d\n",
        ldb, nrowb);
    return 1;
  }

#endif

  return 0;
}
}  // namespace Impl

///
/// Serial Impl
/// ===========

///
/// NT/NT
///

#if defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL) && defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_BATCHED) && \
    defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_COMPACT_BATCHED)
template <>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int SerialGemm<Trans::NoTranspose, Trans::NoTranspose, Algo::Gemm::CompactMKL>::invoke(
    const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta, const CViewType &C) {
  typedef typename CViewType::value_type vector_type;
  // typedef typename vector_type::value_type value_type;

  const int m = C.extent(0), n = C.extent(1), k = A.extent(1);

  static_assert(is_vector<vector_type>::value, "value type is not vector type");
  static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                "AVX, AVX2 and AVX512 is supported");
  const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

  // no error check
  int r_val = 0;
  if (A.stride(0) == 1 && B.stride(0) == 1 && C.stride(0) == 1) {
    mkl_dgemm_compact(MKL_COL_MAJOR, MKL_NOTRANS, MKL_NOTRANS, m, n, k, alpha, (const double *)A.data(), A.stride(1),
                      (const double *)B.data(), B.stride(1), beta, (double *)C.data(), C.stride(1), format,
                      (MKL_INT)vector_type::vector_length);
  } else if (A.stride(1) == 1 && B.stride(1) == 1 && C.stride(1) == 1) {
    mkl_dgemm_compact(MKL_ROW_MAJOR, MKL_NOTRANS, MKL_NOTRANS, m, n, k, alpha, (const double *)A.data(), A.stride(0),
                      (const double *)B.data(), B.stride(0), beta, (double *)C.data(), C.stride(0), format,
                      (MKL_INT)vector_type::vector_length);
  } else {
    r_val = -1;
  }
  return r_val;
}
#endif

template <>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int SerialGemm<Trans::NoTranspose, Trans::NoTranspose, Algo::Gemm::Unblocked>::invoke(
    const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta, const CViewType &C) {
  // Quick return if possible
  const int m = C.extent(0), n = C.extent(1), k = A.extent(1);
  if (m == 0 || n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

  auto info = KokkosBatched::Impl::checkGemmInput<Trans::NoTranspose, Trans::NoTranspose>(A, B, C);
  if (info) return info;

  // C = beta C + alpha A B
  // C (m x n), A(m x k), B(k x n)
  return KokkosBatched::Impl::SerialGemmInternal<Algo::Gemm::Unblocked>::invoke(
      KokkosBlas::Impl::OpID(), KokkosBlas::Impl::OpID(), C.extent(0), C.extent(1), A.extent(1), alpha, A.data(),
      A.stride(0), A.stride(1), B.data(), B.stride(0), B.stride(1), beta, C.data(), C.stride(0), C.stride(1));
}

template <>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int SerialGemm<Trans::NoTranspose, Trans::NoTranspose, Algo::Gemm::Blocked>::invoke(
    const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta, const CViewType &C) {
  // Quick return if possible
  const int m = C.extent(0), n = C.extent(1), k = A.extent(1);
  if (m == 0 || n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

  auto info = KokkosBatched::Impl::checkGemmInput<Trans::NoTranspose, Trans::NoTranspose>(A, B, C);
  if (info) return info;

  // C = beta C + alpha A B
  // C (m x n), A(m x k), B(k x n)
  return KokkosBatched::Impl::SerialGemmInternal<Algo::Gemm::Blocked>::invoke(
      KokkosBlas::Impl::OpID(), KokkosBlas::Impl::OpID(), C.extent(0), C.extent(1), A.extent(1), alpha, A.data(),
      A.stride(0), A.stride(1), B.data(), B.stride(0), B.stride(1), beta, C.data(), C.stride(0), C.stride(1));
}

///
/// T/NT
///

#if defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL) && defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_BATCHED) && \
    defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_COMPACT_BATCHED)
template <>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int SerialGemm<Trans::Transpose, Trans::NoTranspose, Algo::Gemm::CompactMKL>::invoke(
    const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta, const CViewType &C) {
  typedef typename CViewType::value_type vector_type;
  // typedef typename vector_type::value_type value_type;

  const int m = C.extent(0), n = C.extent(1), k = A.extent(0);

  static_assert(is_vector<vector_type>::value, "value type is not vector type");
  static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                "AVX, AVX2 and AVX512 is supported");
  const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

  // no error check
  int r_val = 0;
  if (A.stride(0) == 1 && B.stride(0) == 1 && C.stride(0) == 1) {
    mkl_dgemm_compact(MKL_COL_MAJOR, MKL_TRANS, MKL_NOTRANS, m, n, k, alpha, (const double *)A.data(), A.stride(1),
                      (const double *)B.data(), B.stride(1), beta, (double *)C.data(), C.stride(1), format,
                      (MKL_INT)vector_type::vector_length);
  } else if (A.stride(1) == 1 && B.stride(1) == 1 && C.stride(1) == 1) {
    mkl_dgemm_compact(MKL_ROW_MAJOR, MKL_TRANS, MKL_NOTRANS, m, n, k, alpha, (const double *)A.data(), A.stride(0),
                      (const double *)B.data(), B.stride(0), beta, (double *)C.data(), C.stride(0), format,
                      (MKL_INT)vector_type::vector_length);
  } else {
    r_val = -1;
  }
  return r_val;
}
#endif

template <>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int SerialGemm<Trans::Transpose, Trans::NoTranspose, Algo::Gemm::Unblocked>::invoke(
    const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta, const CViewType &C) {
  // Quick return if possible
  const int m = C.extent(0), n = C.extent(1), k = A.extent(0);
  if (m == 0 || n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

  auto info = KokkosBatched::Impl::checkGemmInput<Trans::Transpose, Trans::NoTranspose>(A, B, C);
  if (info) return info;

  // C = beta C + alpha A^T B
  // C (m x n), A(k x m), B(k x n)
  return KokkosBatched::Impl::SerialGemmInternal<Algo::Gemm::Unblocked>::invoke(
      KokkosBlas::Impl::OpID(), KokkosBlas::Impl::OpID(), C.extent(0), C.extent(1), A.extent(0), alpha, A.data(),
      A.stride(1), A.stride(0), B.data(), B.stride(0), B.stride(1), beta, C.data(), C.stride(0), C.stride(1));
}

template <>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int SerialGemm<Trans::Transpose, Trans::NoTranspose, Algo::Gemm::Blocked>::invoke(
    const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta, const CViewType &C) {
  // Quick return if possible
  const int m = C.extent(0), n = C.extent(1), k = A.extent(0);
  if (m == 0 || n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

  auto info = KokkosBatched::Impl::checkGemmInput<Trans::Transpose, Trans::NoTranspose>(A, B, C);
  if (info) return info;

  // C = beta C + alpha A^T B
  // C (m x n), A(k x m), B(k x n)
  return KokkosBatched::Impl::SerialGemmInternal<Algo::Gemm::Blocked>::invoke(
      KokkosBlas::Impl::OpID(), KokkosBlas::Impl::OpID(), C.extent(0), C.extent(1), A.extent(0), alpha, A.data(),
      A.stride(1), A.stride(0), B.data(), B.stride(0), B.stride(1), beta, C.data(), C.stride(0), C.stride(1));
}

///
/// C/NT
///

#if defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL) && defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_BATCHED) && \
    defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_COMPACT_BATCHED)
template <>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int SerialGemm<Trans::ConjTranspose, Trans::NoTranspose, Algo::Gemm::CompactMKL>::invoke(
    const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta, const CViewType &C) {
  typedef typename CViewType::value_type vector_type;
  // typedef typename vector_type::value_type value_type;

  const int m = C.extent(0), n = C.extent(1), k = A.extent(0);

  static_assert(is_vector<vector_type>::value, "value type is not vector type");
  static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                "AVX, AVX2 and AVX512 is supported");
  const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

  // no error check
  int r_val = 0;
  if (A.stride(0) == 1 && B.stride(0) == 1 && C.stride(0) == 1) {
    mkl_dgemm_compact(MKL_COL_MAJOR, MKL_CONJTRANS, MKL_NOTRANS, m, n, k, alpha, (const double *)A.data(), A.stride(1),
                      (const double *)B.data(), B.stride(1), beta, (double *)C.data(), C.stride(1), format,
                      (MKL_INT)vector_type::vector_length);
  } else if (A.stride(1) == 1 && B.stride(1) == 1 && C.stride(1) == 1) {
    mkl_dgemm_compact(MKL_ROW_MAJOR, MKL_CONJTRANS, MKL_NOTRANS, m, n, k, alpha, (const double *)A.data(), A.stride(0),
                      (const double *)B.data(), B.stride(0), beta, (double *)C.data(), C.stride(0), format,
                      (MKL_INT)vector_type::vector_length);
  } else {
    r_val = -1;
  }
  return r_val;
}
#endif

template <>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int SerialGemm<Trans::ConjTranspose, Trans::NoTranspose, Algo::Gemm::Unblocked>::invoke(
    const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta, const CViewType &C) {
  // Quick return if possible
  const int m = C.extent(0), n = C.extent(1), k = A.extent(0);
  if (m == 0 || n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

  auto info = KokkosBatched::Impl::checkGemmInput<Trans::ConjTranspose, Trans::NoTranspose>(A, B, C);
  if (info) return info;

  // C = beta C + alpha A^H B
  // C (m x n), A(k x m), B(k x n)
  return KokkosBatched::Impl::SerialGemmInternal<Algo::Gemm::Unblocked>::invoke(
      KokkosBlas::Impl::OpConj(), KokkosBlas::Impl::OpID(), C.extent(0), C.extent(1), A.extent(0), alpha, A.data(),
      A.stride(1), A.stride(0), B.data(), B.stride(0), B.stride(1), beta, C.data(), C.stride(0), C.stride(1));
}

template <>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int SerialGemm<Trans::ConjTranspose, Trans::NoTranspose, Algo::Gemm::Blocked>::invoke(
    const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta, const CViewType &C) {
  // Quick return if possible
  const int m = C.extent(0), n = C.extent(1), k = A.extent(0);
  if (m == 0 || n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

  auto info = KokkosBatched::Impl::checkGemmInput<Trans::ConjTranspose, Trans::NoTranspose>(A, B, C);
  if (info) return info;

  // C = beta C + alpha A^H B
  // C (m x n), A(k x m), B(k x n)
  return KokkosBatched::Impl::SerialGemmInternal<Algo::Gemm::Blocked>::invoke(
      KokkosBlas::Impl::OpConj(), KokkosBlas::Impl::OpID(), C.extent(0), C.extent(1), A.extent(0), alpha, A.data(),
      A.stride(1), A.stride(0), B.data(), B.stride(0), B.stride(1), beta, C.data(), C.stride(0), C.stride(1));
}

///
/// NT/T
///

#if defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL) && defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_BATCHED) && \
    defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_COMPACT_BATCHED)
template <>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int SerialGemm<Trans::NoTranspose, Trans::Transpose, Algo::Gemm::CompactMKL>::invoke(
    const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta, const CViewType &C) {
  typedef typename CViewType::value_type vector_type;
  // typedef typename vector_type::value_type value_type;

  const int m = C.extent(0), n = C.extent(1), k = A.extent(1);

  static_assert(is_vector<vector_type>::value, "value type is not vector type");
  static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                "AVX, AVX2 and AVX512 is supported");
  const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

  // no error check
  int r_val = 0;
  if (A.stride(0) == 1 && B.stride(0) == 1 && C.stride(0) == 1) {
    mkl_dgemm_compact(MKL_COL_MAJOR, MKL_NOTRANS, MKL_TRANS, m, n, k, alpha, (const double *)A.data(), A.stride(1),
                      (const double *)B.data(), B.stride(1), beta, (double *)C.data(), C.stride(1), format,
                      (MKL_INT)vector_type::vector_length);
  } else if (A.stride(1) == 1 && B.stride(1) == 1 && C.stride(1) == 1) {
    mkl_dgemm_compact(MKL_ROW_MAJOR, MKL_NOTRANS, MKL_TRANS, m, n, k, alpha, (const double *)A.data(), A.stride(0),
                      (const double *)B.data(), B.stride(0), beta, (double *)C.data(), C.stride(0), format,
                      (MKL_INT)vector_type::vector_length);
  } else {
    r_val = -1;
  }
  return r_val;
}
#endif

template <>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int SerialGemm<Trans::NoTranspose, Trans::Transpose, Algo::Gemm::Unblocked>::invoke(
    const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta, const CViewType &C) {
  // Quick return if possible
  const int m = C.extent(0), n = C.extent(1), k = A.extent(1);
  if (m == 0 || n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

  auto info = KokkosBatched::Impl::checkGemmInput<Trans::NoTranspose, Trans::Transpose>(A, B, C);
  if (info) return info;

  // C = beta C + alpha A B^T
  // C (m x n), A(m x k), B(n x k)
  return KokkosBatched::Impl::SerialGemmInternal<Algo::Gemm::Unblocked>::invoke(
      KokkosBlas::Impl::OpID(), KokkosBlas::Impl::OpID(), C.extent(0), C.extent(1), A.extent(1), alpha, A.data(),
      A.stride(0), A.stride(1), B.data(), B.stride(1), B.stride(0), beta, C.data(), C.stride(0), C.stride(1));
}

template <>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int SerialGemm<Trans::NoTranspose, Trans::Transpose, Algo::Gemm::Blocked>::invoke(
    const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta, const CViewType &C) {
  // Quick return if possible
  const int m = C.extent(0), n = C.extent(1), k = A.extent(1);
  if (m == 0 || n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

  auto info = KokkosBatched::Impl::checkGemmInput<Trans::NoTranspose, Trans::Transpose>(A, B, C);
  if (info) return info;

  // C = beta C + alpha A B^T
  // C (m x n), A(m x k), B(n x k)
  return KokkosBatched::Impl::SerialGemmInternal<Algo::Gemm::Blocked>::invoke(
      KokkosBlas::Impl::OpID(), KokkosBlas::Impl::OpID(), C.extent(0), C.extent(1), A.extent(1), alpha, A.data(),
      A.stride(0), A.stride(1), B.data(), B.stride(1), B.stride(0), beta, C.data(), C.stride(0), C.stride(1));
}

///
/// T/T
///

#if defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL) && defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_BATCHED) && \
    defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_COMPACT_BATCHED)
template <>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int SerialGemm<Trans::Transpose, Trans::Transpose, Algo::Gemm::CompactMKL>::invoke(
    const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta, const CViewType &C) {
  typedef typename CViewType::value_type vector_type;
  // typedef typename vector_type::value_type value_type;

  const int m = C.extent(0), n = C.extent(1), k = A.extent(0);

  static_assert(is_vector<vector_type>::value, "value type is not vector type");
  static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                "AVX, AVX2 and AVX512 is supported");
  const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

  // no error check
  int r_val = 0;
  if (A.stride(0) == 1 && B.stride(0) == 1 && C.stride(0) == 1) {
    mkl_dgemm_compact(MKL_COL_MAJOR, MKL_TRANS, MKL_TRANS, m, n, k, alpha, (const double *)A.data(), A.stride(1),
                      (const double *)B.data(), B.stride(1), beta, (double *)C.data(), C.stride(1), format,
                      (MKL_INT)vector_type::vector_length);
  } else if (A.stride(1) == 1 && B.stride(1) == 1 && C.stride(1) == 1) {
    mkl_dgemm_compact(MKL_ROW_MAJOR, MKL_TRANS, MKL_TRANS, m, n, k, alpha, (const double *)A.data(), A.stride(0),
                      (const double *)B.data(), B.stride(0), beta, (double *)C.data(), C.stride(0), format,
                      (MKL_INT)vector_type::vector_length);
  } else {
    r_val = -1;
  }
  return r_val;
}
#endif

template <>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int SerialGemm<Trans::Transpose, Trans::Transpose, Algo::Gemm::Unblocked>::invoke(
    const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta, const CViewType &C) {
  // Quick return if possible
  const int m = C.extent(0), n = C.extent(1), k = A.extent(0);
  if (m == 0 || n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

  auto info = KokkosBatched::Impl::checkGemmInput<Trans::Transpose, Trans::Transpose>(A, B, C);
  if (info) return info;

  // C = beta C + alpha A^T B^T
  // C (m x n), A(k x m), B(n x k)
  return KokkosBatched::Impl::SerialGemmInternal<Algo::Gemm::Unblocked>::invoke(
      KokkosBlas::Impl::OpID(), KokkosBlas::Impl::OpID(), C.extent(0), C.extent(1), A.extent(0), alpha, A.data(),
      A.stride(1), A.stride(0), B.data(), B.stride(1), B.stride(0), beta, C.data(), C.stride(0), C.stride(1));
}

template <>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int SerialGemm<Trans::Transpose, Trans::Transpose, Algo::Gemm::Blocked>::invoke(
    const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta, const CViewType &C) {
  // Quick return if possible
  const int m = C.extent(0), n = C.extent(1), k = A.extent(0);
  if (m == 0 || n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

  auto info = KokkosBatched::Impl::checkGemmInput<Trans::Transpose, Trans::Transpose>(A, B, C);
  if (info) return info;

  // C = beta C + alpha A^T B^T
  // C (m x n), A(k x m), B(n x k)
  return KokkosBatched::Impl::SerialGemmInternal<Algo::Gemm::Blocked>::invoke(
      KokkosBlas::Impl::OpID(), KokkosBlas::Impl::OpID(), C.extent(0), C.extent(1), A.extent(0), alpha, A.data(),
      A.stride(1), A.stride(0), B.data(), B.stride(1), B.stride(0), beta, C.data(), C.stride(0), C.stride(1));
}

///
/// C/T
///

#if defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL) && defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_BATCHED) && \
    defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_COMPACT_BATCHED)
template <>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int SerialGemm<Trans::ConjTranspose, Trans::Transpose, Algo::Gemm::CompactMKL>::invoke(
    const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta, const CViewType &C) {
  typedef typename CViewType::value_type vector_type;
  // typedef typename vector_type::value_type value_type;

  const int m = C.extent(0), n = C.extent(1), k = A.extent(0);

  static_assert(is_vector<vector_type>::value, "value type is not vector type");
  static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                "AVX, AVX2 and AVX512 is supported");
  const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

  // no error check
  int r_val = 0;
  if (A.stride(0) == 1 && B.stride(0) == 1 && C.stride(0) == 1) {
    mkl_dgemm_compact(MKL_COL_MAJOR, MKL_CONJTRANS, MKL_TRANS, m, n, k, alpha, (const double *)A.data(), A.stride(1),
                      (const double *)B.data(), B.stride(1), beta, (double *)C.data(), C.stride(1), format,
                      (MKL_INT)vector_type::vector_length);
  } else if (A.stride(1) == 1 && B.stride(1) == 1 && C.stride(1) == 1) {
    mkl_dgemm_compact(MKL_ROW_MAJOR, MKL_CONJTRANS, MKL_TRANS, m, n, k, alpha, (const double *)A.data(), A.stride(0),
                      (const double *)B.data(), B.stride(0), beta, (double *)C.data(), C.stride(0), format,
                      (MKL_INT)vector_type::vector_length);
  } else {
    r_val = -1;
  }
  return r_val;
}
#endif

template <>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int SerialGemm<Trans::ConjTranspose, Trans::Transpose, Algo::Gemm::Unblocked>::invoke(
    const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta, const CViewType &C) {
  // Quick return if possible
  const int m = C.extent(0), n = C.extent(1), k = A.extent(0);
  if (m == 0 || n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

  auto info = KokkosBatched::Impl::checkGemmInput<Trans::ConjTranspose, Trans::Transpose>(A, B, C);
  if (info) return info;

  // C = beta C + alpha A^H B^T
  // C (m x n), A(k x m), B(n x k)
  return KokkosBatched::Impl::SerialGemmInternal<Algo::Gemm::Unblocked>::invoke(
      KokkosBlas::Impl::OpConj(), KokkosBlas::Impl::OpID(), C.extent(0), C.extent(1), A.extent(0), alpha, A.data(),
      A.stride(1), A.stride(0), B.data(), B.stride(1), B.stride(0), beta, C.data(), C.stride(0), C.stride(1));
}

template <>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int SerialGemm<Trans::ConjTranspose, Trans::Transpose, Algo::Gemm::Blocked>::invoke(
    const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta, const CViewType &C) {
  // Quick return if possible
  const int m = C.extent(0), n = C.extent(1), k = A.extent(0);
  if (m == 0 || n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

  auto info = KokkosBatched::Impl::checkGemmInput<Trans::ConjTranspose, Trans::Transpose>(A, B, C);
  if (info) return info;

  // C = beta C + alpha A^H B^T
  // C (m x n), A(k x m), B(n x k)
  return KokkosBatched::Impl::SerialGemmInternal<Algo::Gemm::Blocked>::invoke(
      KokkosBlas::Impl::OpConj(), KokkosBlas::Impl::OpID(), C.extent(0), C.extent(1), A.extent(0), alpha, A.data(),
      A.stride(1), A.stride(0), B.data(), B.stride(1), B.stride(0), beta, C.data(), C.stride(0), C.stride(1));
}

///
/// NT/C
///

#if defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL) && defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_BATCHED) && \
    defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_COMPACT_BATCHED)
template <>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int SerialGemm<Trans::NoTranspose, Trans::ConjTranspose, Algo::Gemm::CompactMKL>::invoke(
    const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta, const CViewType &C) {
  typedef typename CViewType::value_type vector_type;
  // typedef typename vector_type::value_type value_type;

  const int m = C.extent(0), n = C.extent(1), k = A.extent(1);

  static_assert(is_vector<vector_type>::value, "value type is not vector type");
  static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                "AVX, AVX2 and AVX512 is supported");
  const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

  // no error check
  int r_val = 0;
  if (A.stride(0) == 1 && B.stride(0) == 1 && C.stride(0) == 1) {
    mkl_dgemm_compact(MKL_COL_MAJOR, MKL_NOTRANS, MKL_CONJTRANS, m, n, k, alpha, (const double *)A.data(), A.stride(1),
                      (const double *)B.data(), B.stride(1), beta, (double *)C.data(), C.stride(1), format,
                      (MKL_INT)vector_type::vector_length);
  } else if (A.stride(1) == 1 && B.stride(1) == 1 && C.stride(1) == 1) {
    mkl_dgemm_compact(MKL_ROW_MAJOR, MKL_NOTRANS, MKL_CONJTRANS, m, n, k, alpha, (const double *)A.data(), A.stride(0),
                      (const double *)B.data(), B.stride(0), beta, (double *)C.data(), C.stride(0), format,
                      (MKL_INT)vector_type::vector_length);
  } else {
    r_val = -1;
  }
  return r_val;
}
#endif

template <>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int SerialGemm<Trans::NoTranspose, Trans::ConjTranspose, Algo::Gemm::Unblocked>::invoke(
    const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta, const CViewType &C) {
  // Quick return if possible
  const int m = C.extent(0), n = C.extent(1), k = A.extent(1);
  if (m == 0 || n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

  auto info = KokkosBatched::Impl::checkGemmInput<Trans::NoTranspose, Trans::ConjTranspose>(A, B, C);
  if (info) return info;

  // C = beta C + alpha A B^H
  // C (m x n), A(m x k), B(n x k)
  return KokkosBatched::Impl::SerialGemmInternal<Algo::Gemm::Unblocked>::invoke(
      KokkosBlas::Impl::OpID(), KokkosBlas::Impl::OpConj(), C.extent(0), C.extent(1), A.extent(1), alpha, A.data(),
      A.stride(0), A.stride(1), B.data(), B.stride(1), B.stride(0), beta, C.data(), C.stride(0), C.stride(1));
}

template <>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int SerialGemm<Trans::NoTranspose, Trans::ConjTranspose, Algo::Gemm::Blocked>::invoke(
    const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta, const CViewType &C) {
  // Quick return if possible
  const int m = C.extent(0), n = C.extent(1), k = A.extent(1);
  if (m == 0 || n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

  auto info = KokkosBatched::Impl::checkGemmInput<Trans::NoTranspose, Trans::ConjTranspose>(A, B, C);
  if (info) return info;

  // C = beta C + alpha A B^H
  // C (m x n), A(m x k), B(n x k)
  return KokkosBatched::Impl::SerialGemmInternal<Algo::Gemm::Blocked>::invoke(
      KokkosBlas::Impl::OpID(), KokkosBlas::Impl::OpConj(), C.extent(0), C.extent(1), A.extent(1), alpha, A.data(),
      A.stride(0), A.stride(1), B.data(), B.stride(1), B.stride(0), beta, C.data(), C.stride(0), C.stride(1));
}

///
/// T/C
///

#if defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL) && defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_BATCHED) && \
    defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_COMPACT_BATCHED)
template <>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int SerialGemm<Trans::Transpose, Trans::ConjTranspose, Algo::Gemm::CompactMKL>::invoke(
    const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta, const CViewType &C) {
  typedef typename CViewType::value_type vector_type;
  // typedef typename vector_type::value_type value_type;

  const int m = C.extent(0), n = C.extent(1), k = A.extent(0);

  static_assert(is_vector<vector_type>::value, "value type is not vector type");
  static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                "AVX, AVX2 and AVX512 is supported");
  const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

  // no error check
  int r_val = 0;
  if (A.stride(0) == 1 && B.stride(0) == 1 && C.stride(0) == 1) {
    mkl_dgemm_compact(MKL_COL_MAJOR, MKL_TRANS, MKL_CONJTRANS, m, n, k, alpha, (const double *)A.data(), A.stride(1),
                      (const double *)B.data(), B.stride(1), beta, (double *)C.data(), C.stride(1), format,
                      (MKL_INT)vector_type::vector_length);
  } else if (A.stride(1) == 1 && B.stride(1) == 1 && C.stride(1) == 1) {
    mkl_dgemm_compact(MKL_ROW_MAJOR, MKL_TRANS, MKL_CONJTRANS, m, n, k, alpha, (const double *)A.data(), A.stride(0),
                      (const double *)B.data(), B.stride(0), beta, (double *)C.data(), C.stride(0), format,
                      (MKL_INT)vector_type::vector_length);
  } else {
    r_val = -1;
  }
  return r_val;
}
#endif

template <>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int SerialGemm<Trans::Transpose, Trans::ConjTranspose, Algo::Gemm::Unblocked>::invoke(
    const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta, const CViewType &C) {
  // Quick return if possible
  const int m = C.extent(0), n = C.extent(1), k = A.extent(0);
  if (m == 0 || n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

  auto info = KokkosBatched::Impl::checkGemmInput<Trans::Transpose, Trans::ConjTranspose>(A, B, C);
  if (info) return info;

  // C = beta C + alpha A^T B^H
  // C (m x n), A(k x m), B(n x k)
  return KokkosBatched::Impl::SerialGemmInternal<Algo::Gemm::Unblocked>::invoke(
      KokkosBlas::Impl::OpID(), KokkosBlas::Impl::OpConj(), C.extent(0), C.extent(1), A.extent(0), alpha, A.data(),
      A.stride(1), A.stride(0), B.data(), B.stride(1), B.stride(0), beta, C.data(), C.stride(0), C.stride(1));
}

template <>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int SerialGemm<Trans::Transpose, Trans::ConjTranspose, Algo::Gemm::Blocked>::invoke(
    const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta, const CViewType &C) {
  // Quick return if possible
  const int m = C.extent(0), n = C.extent(1), k = A.extent(0);
  if (m == 0 || n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

  auto info = KokkosBatched::Impl::checkGemmInput<Trans::Transpose, Trans::Transpose>(A, B, C);
  if (info) return info;

  // C = beta C + alpha A^T B^H
  // C (m x n), A(k x m), B(n x k)
  return KokkosBatched::Impl::SerialGemmInternal<Algo::Gemm::Blocked>::invoke(
      KokkosBlas::Impl::OpID(), KokkosBlas::Impl::OpConj(), C.extent(0), C.extent(1), A.extent(0), alpha, A.data(),
      A.stride(1), A.stride(0), B.data(), B.stride(1), B.stride(0), beta, C.data(), C.stride(0), C.stride(1));
}

///
/// C/C
///

#if defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL) && defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_BATCHED) && \
    defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_COMPACT_BATCHED)
template <>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int SerialGemm<Trans::ConjTranspose, Trans::ConjTranspose, Algo::Gemm::CompactMKL>::invoke(
    const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta, const CViewType &C) {
  typedef typename CViewType::value_type vector_type;
  // typedef typename vector_type::value_type value_type;

  const int m = C.extent(0), n = C.extent(1), k = A.extent(0);

  static_assert(is_vector<vector_type>::value, "value type is not vector type");
  static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                "AVX, AVX2 and AVX512 is supported");
  const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

  // no error check
  int r_val = 0;
  if (A.stride(0) == 1 && B.stride(0) == 1 && C.stride(0) == 1) {
    mkl_dgemm_compact(MKL_COL_MAJOR, MKL_CONJTRANS, MKL_CONJTRANS, m, n, k, alpha, (const double *)A.data(),
                      A.stride(1), (const double *)B.data(), B.stride(1), beta, (double *)C.data(), C.stride(1), format,
                      (MKL_INT)vector_type::vector_length);
  } else if (A.stride(1) == 1 && B.stride(1) == 1 && C.stride(1) == 1) {
    mkl_dgemm_compact(MKL_ROW_MAJOR, MKL_CONJTRANS, MKL_CONJTRANS, m, n, k, alpha, (const double *)A.data(),
                      A.stride(0), (const double *)B.data(), B.stride(0), beta, (double *)C.data(), C.stride(0), format,
                      (MKL_INT)vector_type::vector_length);
  } else {
    r_val = -1;
  }
  return r_val;
}
#endif

template <>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int SerialGemm<Trans::ConjTranspose, Trans::ConjTranspose, Algo::Gemm::Unblocked>::invoke(
    const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta, const CViewType &C) {
  // Quick return if possible
  const int m = C.extent(0), n = C.extent(1), k = A.extent(0);
  if (m == 0 || n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

  auto info = KokkosBatched::Impl::checkGemmInput<Trans::ConjTranspose, Trans::ConjTranspose>(A, B, C);
  if (info) return info;

  // C = beta C + alpha A^H B^H
  // C (m x n), A(k x m), B(n x k)
  return KokkosBatched::Impl::SerialGemmInternal<Algo::Gemm::Unblocked>::invoke(
      KokkosBlas::Impl::OpConj(), KokkosBlas::Impl::OpConj(), C.extent(0), C.extent(1), A.extent(0), alpha, A.data(),
      A.stride(1), A.stride(0), B.data(), B.stride(1), B.stride(0), beta, C.data(), C.stride(0), C.stride(1));
}

template <>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int SerialGemm<Trans::ConjTranspose, Trans::ConjTranspose, Algo::Gemm::Blocked>::invoke(
    const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta, const CViewType &C) {
  // Quick return if possible
  const int m = C.extent(0), n = C.extent(1), k = A.extent(0);
  if (m == 0 || n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

  auto info = KokkosBatched::Impl::checkGemmInput<Trans::ConjTranspose, Trans::ConjTranspose>(A, B, C);
  if (info) return info;

  // C = beta C + alpha A^H B^H
  // C (m x n), A(k x m), B(n x k)
  return KokkosBatched::Impl::SerialGemmInternal<Algo::Gemm::Blocked>::invoke(
      KokkosBlas::Impl::OpConj(), KokkosBlas::Impl::OpConj(), C.extent(0), C.extent(1), A.extent(0), alpha, A.data(),
      A.stride(1), A.stride(0), B.data(), B.stride(1), B.stride(0), beta, C.data(), C.stride(0), C.stride(1));
}

}  // namespace KokkosBatched

#endif
