// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_TRSM_SERIAL_IMPL_HPP
#define KOKKOSBATCHED_TRSM_SERIAL_IMPL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Trsm_Serial_Internal.hpp"

namespace KokkosBatched {
namespace Impl {
template <typename ArgSide, typename AViewType, typename BViewType>
KOKKOS_INLINE_FUNCTION static int checkTrsmInput([[maybe_unused]] const AViewType &A,
                                                 [[maybe_unused]] const BViewType &B) {
  static_assert(Kokkos::is_view_v<AViewType>, "KokkosBatched::trsm: AViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<BViewType>, "KokkosBatched::trsm: BViewType is not a Kokkos::View.");
  static_assert(AViewType::rank == 2, "KokkosBatched::trsm: AViewType must have rank 2.");
  static_assert(BViewType::rank == 1 || BViewType::rank == 2, "KokkosBatched::trsm: BViewType must have rank 1 or 2.");
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  const int m = B.extent(0), n = B.extent(1);
  const int nrowa = std::is_same_v<ArgSide, Side::Left> ? m : n;
  const int lda   = A.extent(0);

  if (lda < Kokkos::max(1, nrowa)) {
    Kokkos::printf(
        "KokkosBatched::trsm: leading dimension of A must not be smaller than "
        "max(1, nrowa): "
        "lda = %d, nrowa = %d\n",
        lda, nrowa);
    return 1;
  }
#endif
  return 0;
}
}  // namespace Impl

///
/// L/L/NT
///
/// B := inv(tril(A)) (alpha*B)
/// A(m x m), B(m x n)

#if defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL) && defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_BATCHED) && \
    defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_COMPACT_BATCHED)
template <typename ArgDiag>
struct SerialTrsm<Side::Left, Uplo::Lower, Trans::NoTranspose, ArgDiag, Algo::Trsm::CompactMKL> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    typedef typename BViewType::value_type vector_type;
    // typedef typename vector_type::value_type value_type;

    const int m = B.extent(0), n = B_extent_1;

    static_assert(is_vector<vector_type>::value, "value type is not vector type");
    static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                  "AVX, AVX2 and AVX512 is supported");
    const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

    // no error check
    int r_val = 0;
    if (A.stride(0) == 1 && B.stride(0) == 1) {
      mkl_dtrsm_compact(MKL_COL_MAJOR, MKL_LEFT, MKL_LOWER, MKL_NOTRANS,
                        ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT, m, n, alpha, (const double *)A.data(),
                        A.stride(1), (double *)B.data(), B_stride_1, format, (MKL_INT)vector_type::vector_length);
    } else if (A.stride(1) == 1 && B_stride_1 == 1) {
      mkl_dtrsm_compact(MKL_ROW_MAJOR, MKL_LEFT, MKL_LOWER, MKL_NOTRANS,
                        ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT, m, n, alpha, (const double *)A.data(),
                        A.stride(0), (double *)B.data(), B.stride(0), format, (MKL_INT)vector_type::vector_length);
    } else {
      r_val = -1;
    }
    return r_val;
  }
};
#endif

template <typename ArgDiag>
struct SerialTrsm<Side::Left, Uplo::Lower, Trans::NoTranspose, ArgDiag, Algo::Trsm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    auto info = KokkosBatched::Impl::checkTrsmInput<Side::Left>(A, B);
    if (info) return info;

    return KokkosBatched::Impl::SerialTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, B.extent(0), B_extent_1, alpha, A.data(), A.stride(0), A.stride(1), B.data(),
        B.stride(0), B_stride_1);
  }
};

template <typename ArgDiag>
struct SerialTrsm<Side::Left, Uplo::Lower, Trans::NoTranspose, ArgDiag, Algo::Trsm::Blocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    auto info = KokkosBatched::Impl::checkTrsmInput<Side::Left>(A, B);
    if (info) return info;

    static_assert(AViewType::rank() == 2 && BViewType::rank() == 2);
    return KokkosBatched::Impl::SerialTrsmInternalLeftLower<Algo::Trsm::Blocked>::invoke(
        ArgDiag::use_unit_diag, false, B.extent(0), B_extent_1, alpha, A.data(), A.stride(0), A.stride(1), B.data(),
        B.stride(0), B_stride_1);
  }
};

///
/// L/U/NT
///
/// B := inv(triu(A)) (alpha*B)
/// A(m x m), B(m x n)
#if defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL) && defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_BATCHED) && \
    defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_COMPACT_BATCHED)
template <typename ArgDiag>
struct SerialTrsm<Side::Left, Uplo::Upper, Trans::NoTranspose, ArgDiag, Algo::Trsm::CompactMKL> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    typedef typename BViewType::value_type vector_type;
    // typedef typename vector_type::value_type value_type;

    const int m = B.extent(0), n = B_extent_1;

    static_assert(is_vector<vector_type>::value, "value type is not vector type");
    static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                  "AVX, AVX2 and AVX512 is supported");
    static_assert(AViewType::rank() == 2 && BViewType::rank() == 2);
    const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

    // no error check
    int r_val = 0;
    if (A.stride(0) == 1 && B.stride(0) == 1) {
      mkl_dtrsm_compact(MKL_COL_MAJOR, MKL_LEFT, MKL_UPPER, MKL_NOTRANS,
                        ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT, m, n, alpha, (const double *)A.data(),
                        A.stride(1), (double *)B.data(), B_stride_1, format, (MKL_INT)vector_type::vector_length);
    } else if (A.stride(1) == 1 && B_stride_1 == 1) {
      mkl_dtrsm_compact(MKL_ROW_MAJOR, MKL_LEFT, MKL_UPPER, MKL_NOTRANS,
                        ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT, m, n, alpha, (const double *)A.data(),
                        A.stride(0), (double *)B.data(), B.stride(0), format, (MKL_INT)vector_type::vector_length);
    } else {
      r_val = -1;
    }
    return r_val;
  }
};
#endif

template <typename ArgDiag>
struct SerialTrsm<Side::Left, Uplo::Upper, Trans::NoTranspose, ArgDiag, Algo::Trsm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    auto info = KokkosBatched::Impl::checkTrsmInput<Side::Left>(A, B);
    if (info) return info;

    return KokkosBatched::Impl::SerialTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, B.extent(0), B_extent_1, alpha, A.data(), A.stride(0), A.stride(1), B.data(),
        B.stride(0), B_stride_1);
  }
};

template <typename ArgDiag>
struct SerialTrsm<Side::Left, Uplo::Upper, Trans::NoTranspose, ArgDiag, Algo::Trsm::Blocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    auto info = KokkosBatched::Impl::checkTrsmInput<Side::Left>(A, B);
    if (info) return info;

    return KokkosBatched::Impl::SerialTrsmInternalLeftUpper<Algo::Trsm::Blocked>::invoke(
        ArgDiag::use_unit_diag, false, B.extent(0), B_extent_1, alpha, A.data(), A.stride(0), A.stride(1), B.data(),
        B.stride(0), B_stride_1);
  }
};

///
/// L/L/T
///
/// B := inv(tril(AT)) (alpha*B)
/// A(m x m), B(m x n)

#if defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL) && defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_BATCHED) && \
    defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_COMPACT_BATCHED)
template <typename ArgDiag>
struct SerialTrsm<Side::Left, Uplo::Lower, Trans::Transpose, ArgDiag, Algo::Trsm::CompactMKL> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    typedef typename BViewType::value_type vector_type;
    // typedef typename vector_type::value_type value_type;

    const int m = B.extent(0), n = B_extent_1;

    static_assert(is_vector<vector_type>::value, "value type is not vector type");
    static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                  "AVX, AVX2 and AVX512 is supported");
    const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

    // no error check
    int r_val = 0;
    if (A.stride(0) == 1 && B.stride(0) == 1) {
      mkl_dtrsm_compact(MKL_COL_MAJOR, MKL_LEFT, MKL_LOWER, MKL_TRANS, ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT,
                        m, n, alpha, (const double *)A.data(), A.stride(1), (double *)B.data(), B_stride_1, format,
                        (MKL_INT)vector_type::vector_length);
    } else if (A.stride(1) == 1 && B_stride_1 == 1) {
      mkl_dtrsm_compact(MKL_ROW_MAJOR, MKL_LEFT, MKL_LOWER, MKL_TRANS, ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT,
                        m, n, alpha, (const double *)A.data(), A.stride(0), (double *)B.data(), B.stride(0), format,
                        (MKL_INT)vector_type::vector_length);
    } else {
      r_val = -1;
    }
    return r_val;
  }
};
#endif

template <typename ArgDiag>
struct SerialTrsm<Side::Left, Uplo::Lower, Trans::Transpose, ArgDiag, Algo::Trsm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    auto info = KokkosBatched::Impl::checkTrsmInput<Side::Left>(A, B);
    if (info) return info;

    return KokkosBatched::Impl::SerialTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, B.extent(0), B_extent_1, alpha, A.data(), A.stride(1), A.stride(0), B.data(),
        B.stride(0), B_stride_1);
  }
};

template <typename ArgDiag>
struct SerialTrsm<Side::Left, Uplo::Lower, Trans::Transpose, ArgDiag, Algo::Trsm::Blocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    auto info = KokkosBatched::Impl::checkTrsmInput<Side::Left>(A, B);
    if (info) return info;

    return KokkosBatched::Impl::SerialTrsmInternalLeftUpper<Algo::Trsm::Blocked>::invoke(
        ArgDiag::use_unit_diag, false, B.extent(0), B_extent_1, alpha, A.data(), A.stride(1), A.stride(0), B.data(),
        B.stride(0), B_stride_1);
  }
};

///
/// L/U/T
///
/// B := inv(triu(AT)) (alpha*B)
/// A(m x m), B(m x n)
#if defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL) && defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_BATCHED) && \
    defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_COMPACT_BATCHED)
template <typename ArgDiag>
struct SerialTrsm<Side::Left, Uplo::Upper, Trans::Transpose, ArgDiag, Algo::Trsm::CompactMKL> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    typedef typename BViewType::value_type vector_type;
    // typedef typename vector_type::value_type value_type;

    const int m = B.extent(0), n = B_extent_1;

    static_assert(is_vector<vector_type>::value, "value type is not vector type");
    static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                  "AVX, AVX2 and AVX512 is supported");
    const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

    // no error check
    int r_val = 0;
    if (A.stride(0) == 1 && B.stride(0) == 1) {
      mkl_dtrsm_compact(MKL_COL_MAJOR, MKL_LEFT, MKL_UPPER, MKL_TRANS, ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT,
                        m, n, alpha, (const double *)A.data(), A.stride(1), (double *)B.data(), B_stride_1, format,
                        (MKL_INT)vector_type::vector_length);
    } else if (A.stride(1) == 1 && B_stride_1 == 1) {
      mkl_dtrsm_compact(MKL_ROW_MAJOR, MKL_LEFT, MKL_UPPER, MKL_TRANS, ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT,
                        m, n, alpha, (const double *)A.data(), A.stride(0), (double *)B.data(), B.stride(0), format,
                        (MKL_INT)vector_type::vector_length);
    } else {
      r_val = -1;
    }
    return r_val;
  }
};
#endif

template <typename ArgDiag>
struct SerialTrsm<Side::Left, Uplo::Upper, Trans::Transpose, ArgDiag, Algo::Trsm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    auto info = KokkosBatched::Impl::checkTrsmInput<Side::Left>(A, B);
    if (info) return info;

    return KokkosBatched::Impl::SerialTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, B.extent(0), B_extent_1, alpha, A.data(), A.stride(1), A.stride(0), B.data(),
        B.stride(0), B_stride_1);
  }
};

template <typename ArgDiag>
struct SerialTrsm<Side::Left, Uplo::Upper, Trans::Transpose, ArgDiag, Algo::Trsm::Blocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    auto info = KokkosBatched::Impl::checkTrsmInput<Side::Left>(A, B);
    if (info) return info;

    return KokkosBatched::Impl::SerialTrsmInternalLeftLower<Algo::Trsm::Blocked>::invoke(
        ArgDiag::use_unit_diag, false, B.extent(0), B_extent_1, alpha, A.data(), A.stride(1), A.stride(0), B.data(),
        B.stride(0), B_stride_1);
  }
};

///
/// L/L/C
///
/// B := inv(tril(AH)) (alpha*B)
/// A(m x m), B(m x n)

#if defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL) && defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_BATCHED) && \
    defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_COMPACT_BATCHED)
template <typename ArgDiag>
struct SerialTrsm<Side::Left, Uplo::Lower, Trans::ConjTranspose, ArgDiag, Algo::Trsm::CompactMKL> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    typedef typename BViewType::value_type vector_type;
    // typedef typename vector_type::value_type value_type;

    const int m = B.extent(0), n = B_extent_1;

    static_assert(is_vector<vector_type>::value, "value type is not vector type");
    static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                  "AVX, AVX2 and AVX512 is supported");
    const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

    // no error check
    int r_val = 0;
    if (A.stride(0) == 1 && B.stride(0) == 1) {
      mkl_dtrsm_compact(MKL_COL_MAJOR, MKL_LEFT, MKL_LOWER, MKL_CONJTRANS,
                        ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT, m, n, alpha, (const double *)A.data(),
                        A.stride(1), (double *)B.data(), B_stride_1, format, (MKL_INT)vector_type::vector_length);
    } else if (A.stride(1) == 1 && B_stride_1 == 1) {
      mkl_dtrsm_compact(MKL_ROW_MAJOR, MKL_LEFT, MKL_LOWER, MKL_CONJTRANS,
                        ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT, m, n, alpha, (const double *)A.data(),
                        A.stride(0), (double *)B.data(), B.stride(0), format, (MKL_INT)vector_type::vector_length);
    } else {
      r_val = -1;
    }
    return r_val;
  }
};
#endif

template <typename ArgDiag>
struct SerialTrsm<Side::Left, Uplo::Lower, Trans::ConjTranspose, ArgDiag, Algo::Trsm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    auto info = KokkosBatched::Impl::checkTrsmInput<Side::Left>(A, B);
    if (info) return info;

    return KokkosBatched::Impl::SerialTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, true, B.extent(0), B_extent_1, alpha, A.data(), A.stride(1), A.stride(0), B.data(),
        B.stride(0), B_stride_1);
  }
};

// [TO DO] ConjTranspose is not supported yet

///
/// L/U/C
///
/// B := inv(triu(AH)) (alpha*B)
/// A(m x m), B(m x n)
#if defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL) && defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_BATCHED) && \
    defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_COMPACT_BATCHED)
template <typename ArgDiag>
struct SerialTrsm<Side::Left, Uplo::Upper, Trans::ConjTranspose, ArgDiag, Algo::Trsm::CompactMKL> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    typedef typename BViewType::value_type vector_type;
    // typedef typename vector_type::value_type value_type;

    const int m = B.extent(0), n = B_extent_1;

    static_assert(is_vector<vector_type>::value, "value type is not vector type");
    static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                  "AVX, AVX2 and AVX512 is supported");
    const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

    // no error check
    int r_val = 0;
    if (A.stride(0) == 1 && B.stride(0) == 1) {
      mkl_dtrsm_compact(MKL_COL_MAJOR, MKL_LEFT, MKL_UPPER, MKL_CONJTRANS,
                        ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT, m, n, alpha, (const double *)A.data(),
                        A.stride(1), (double *)B.data(), B_stride_1, format, (MKL_INT)vector_type::vector_length);
    } else if (A.stride(1) == 1 && B_stride_1 == 1) {
      mkl_dtrsm_compact(MKL_ROW_MAJOR, MKL_LEFT, MKL_UPPER, MKL_CONJTRANS,
                        ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT, m, n, alpha, (const double *)A.data(),
                        A.stride(0), (double *)B.data(), B.stride(0), format, (MKL_INT)vector_type::vector_length);
    } else {
      r_val = -1;
    }
    return r_val;
  }
};
#endif

template <typename ArgDiag>
struct SerialTrsm<Side::Left, Uplo::Upper, Trans::ConjTranspose, ArgDiag, Algo::Trsm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    auto info = KokkosBatched::Impl::checkTrsmInput<Side::Left>(A, B);
    if (info) return info;

    return KokkosBatched::Impl::SerialTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, true, B.extent(0), B_extent_1, alpha, A.data(), A.stride(1), A.stride(0), B.data(),
        B.stride(0), B_stride_1);
  }
};

// [TO DO] ConjTranspose is not supported yet

///
/// R/L/NT
///
/// B := (alpha*B) inv(tril(A))
/// A(n x n), B(m x n)

#if defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL) && defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_BATCHED) && \
    defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_COMPACT_BATCHED)
template <typename ArgDiag>
struct SerialTrsm<Side::Right, Uplo::Lower, Trans::NoTranspose, ArgDiag, Algo::Trsm::CompactMKL> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2 && BViewType::rank() == 2);
    typedef typename BViewType::value_type vector_type;
    // typedef typename vector_type::value_type value_type;

    const int m = B.extent(0), n = B.extent(1);

    static_assert(is_vector<vector_type>::value, "value type is not vector type");
    static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                  "AVX, AVX2 and AVX512 is supported");
    const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

    // no error check
    int r_val = 0;
    if (A.stride(0) == 1 && B.stride(0) == 1) {
      mkl_dtrsm_compact(MKL_COL_MAJOR, MKL_RIGHT, MKL_LOWER, MKL_NOTRANS,
                        ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT, m, n, alpha, (const double *)A.data(),
                        A.stride(1), (double *)B.data(), B.stride(1), format, (MKL_INT)vector_type::vector_length);
    } else if (A.stride(1) == 1 && B.stride(1) == 1) {
      mkl_dtrsm_compact(MKL_ROW_MAJOR, MKL_RIGHT, MKL_LOWER, MKL_NOTRANS,
                        ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT, m, n, alpha, (const double *)A.data(),
                        A.stride(0), (double *)B.data(), B.stride(0), format, (MKL_INT)vector_type::vector_length);
    } else {
      r_val = -1;
    }
    return r_val;
  }
};
#endif

template <typename ArgDiag>
struct SerialTrsm<Side::Right, Uplo::Lower, Trans::NoTranspose, ArgDiag, Algo::Trsm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2 && BViewType::rank() == 2);
    // Quick return if possible
    if (B.extent(0) == 0 || B.extent(1) == 0) return 0;

    auto info = KokkosBatched::Impl::checkTrsmInput<Side::Right>(A, B);
    if (info) return info;

    return KokkosBatched::Impl::SerialTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, B.extent(1), B.extent(0), alpha, A.data(), A.stride(1), A.stride(0), B.data(),
        B.stride(1), B.stride(0));
  }
};

template <typename ArgDiag>
struct SerialTrsm<Side::Right, Uplo::Lower, Trans::NoTranspose, ArgDiag, Algo::Trsm::Blocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2 && BViewType::rank() == 2);
    // Quick return if possible
    if (B.extent(0) == 0 || B.extent(1) == 0) return 0;

    auto info = KokkosBatched::Impl::checkTrsmInput<Side::Right>(A, B);
    if (info) return info;

    return KokkosBatched::Impl::SerialTrsmInternalLeftUpper<Algo::Trsm::Blocked>::invoke(
        ArgDiag::use_unit_diag, false, B.extent(1), B.extent(0), alpha, A.data(), A.stride(1), A.stride(0), B.data(),
        B.stride(1), B.stride(0));
  }
};

///
/// R/U/NT
///
/// B := (alpha*B) inv(triu(A))
/// A(n x n), B(m x n)
#if defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL) && defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_BATCHED) && \
    defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_COMPACT_BATCHED)
template <typename ArgDiag>
struct SerialTrsm<Side::Right, Uplo::Upper, Trans::NoTranspose, ArgDiag, Algo::Trsm::CompactMKL> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2 && BViewType::rank() == 2);
    typedef typename BViewType::value_type vector_type;
    // typedef typename vector_type::value_type value_type;

    const int m = B.extent(0), n = B.extent(1);

    static_assert(is_vector<vector_type>::value, "value type is not vector type");
    static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                  "AVX, AVX2 and AVX512 is supported");
    const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

    // no error check
    int r_val = 0;
    if (A.stride(0) == 1 && B.stride(0) == 1) {
      mkl_dtrsm_compact(MKL_COL_MAJOR, MKL_RIGHT, MKL_UPPER, MKL_NOTRANS,
                        ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT, m, n, alpha, (const double *)A.data(),
                        A.stride(1), (double *)B.data(), B.stride(1), format, (MKL_INT)vector_type::vector_length);
    } else if (A.stride(1) == 1 && B.stride(1) == 1) {
      mkl_dtrsm_compact(MKL_ROW_MAJOR, MKL_RIGHT, MKL_UPPER, MKL_NOTRANS,
                        ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT, m, n, alpha, (const double *)A.data(),
                        A.stride(0), (double *)B.data(), B.stride(0), format, (MKL_INT)vector_type::vector_length);
    } else {
      r_val = -1;
    }
    return r_val;
  }
};
#endif

template <typename ArgDiag>
struct SerialTrsm<Side::Right, Uplo::Upper, Trans::NoTranspose, ArgDiag, Algo::Trsm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2 && BViewType::rank() == 2);
    // Quick return if possible
    if (B.extent(0) == 0 || B.extent(1) == 0) return 0;

    auto info = KokkosBatched::Impl::checkTrsmInput<Side::Right>(A, B);
    if (info) return info;

    return KokkosBatched::Impl::SerialTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, B.extent(1), B.extent(0), alpha, A.data(), A.stride(1), A.stride(0), B.data(),
        B.stride(1), B.stride(0));
  }
};

template <typename ArgDiag>
struct SerialTrsm<Side::Right, Uplo::Upper, Trans::NoTranspose, ArgDiag, Algo::Trsm::Blocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2 && BViewType::rank() == 2);
    // Quick return if possible
    if (B.extent(0) == 0 || B.extent(1) == 0) return 0;

    auto info = KokkosBatched::Impl::checkTrsmInput<Side::Right>(A, B);
    if (info) return info;

    return KokkosBatched::Impl::SerialTrsmInternalLeftLower<Algo::Trsm::Blocked>::invoke(
        ArgDiag::use_unit_diag, false, B.extent(1), B.extent(0), alpha, A.data(), A.stride(1), A.stride(0), B.data(),
        B.stride(1), B.stride(0));
  }
};

///
/// R/L/T
///
/// B := (alpha*B) inv(tril(AT))
/// A(n x n), B(m x n)

#if defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL) && defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_BATCHED) && \
    defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_COMPACT_BATCHED)
template <typename ArgDiag>
struct SerialTrsm<Side::Right, Uplo::Lower, Trans::Transpose, ArgDiag, Algo::Trsm::CompactMKL> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2 && BViewType::rank() == 2);
    typedef typename BViewType::value_type vector_type;
    // typedef typename vector_type::value_type value_type;

    const int m = B.extent(0), n = B.extent(1);

    static_assert(is_vector<vector_type>::value, "value type is not vector type");
    static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                  "AVX, AVX2 and AVX512 is supported");
    const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

    // no error check
    int r_val = 0;
    if (A.stride(0) == 1 && B.stride(0) == 1) {
      mkl_dtrsm_compact(MKL_COL_MAJOR, MKL_RIGHT, MKL_LOWER, MKL_TRANS, ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT,
                        m, n, alpha, (const double *)A.data(), A.stride(1), (double *)B.data(), B.stride(1), format,
                        (MKL_INT)vector_type::vector_length);
    } else if (A.stride(1) == 1 && B.stride(1) == 1) {
      mkl_dtrsm_compact(MKL_ROW_MAJOR, MKL_RIGHT, MKL_LOWER, MKL_TRANS, ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT,
                        m, n, alpha, (const double *)A.data(), A.stride(0), (double *)B.data(), B.stride(0), format,
                        (MKL_INT)vector_type::vector_length);
    } else {
      r_val = -1;
    }
    return r_val;
  }
};
#endif

template <typename ArgDiag>
struct SerialTrsm<Side::Right, Uplo::Lower, Trans::Transpose, ArgDiag, Algo::Trsm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2 && BViewType::rank() == 2);
    // Quick return if possible
    if (B.extent(0) == 0 || B.extent(1) == 0) return 0;

    auto info = KokkosBatched::Impl::checkTrsmInput<Side::Right>(A, B);
    if (info) return info;

    return KokkosBatched::Impl::SerialTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, B.extent(1), B.extent(0), alpha, A.data(), A.stride(0), A.stride(1), B.data(),
        B.stride(1), B.stride(0));
  }
};

template <typename ArgDiag>
struct SerialTrsm<Side::Right, Uplo::Lower, Trans::Transpose, ArgDiag, Algo::Trsm::Blocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2 && BViewType::rank() == 2);
    // Quick return if possible
    if (B.extent(0) == 0 || B.extent(1) == 0) return 0;

    auto info = KokkosBatched::Impl::checkTrsmInput<Side::Right>(A, B);
    if (info) return info;

    return KokkosBatched::Impl::SerialTrsmInternalLeftLower<Algo::Trsm::Blocked>::invoke(
        ArgDiag::use_unit_diag, false, B.extent(1), B.extent(0), alpha, A.data(), A.stride(0), A.stride(1), B.data(),
        B.stride(1), B.stride(0));
  }
};

///
/// R/U/T
///
/// B := (alpha*B) inv(triu(AT))
/// A(n x n), B(m x n)
#if defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL) && defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_BATCHED) && \
    defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_COMPACT_BATCHED)
template <typename ArgDiag>
struct SerialTrsm<Side::Right, Uplo::Upper, Trans::Transpose, ArgDiag, Algo::Trsm::CompactMKL> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2 && BViewType::rank() == 2);
    typedef typename BViewType::value_type vector_type;
    // typedef typename vector_type::value_type value_type;

    const int m = B.extent(0), n = B.extent(1);

    static_assert(is_vector<vector_type>::value, "value type is not vector type");
    static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                  "AVX, AVX2 and AVX512 is supported");
    const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

    // no error check
    int r_val = 0;
    if (A.stride(0) == 1 && B.stride(0) == 1) {
      mkl_dtrsm_compact(MKL_COL_MAJOR, MKL_RIGHT, MKL_UPPER, MKL_TRANS, ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT,
                        m, n, alpha, (const double *)A.data(), A.stride(1), (double *)B.data(), B.stride(1), format,
                        (MKL_INT)vector_type::vector_length);
    } else if (A.stride(1) == 1 && B.stride(1) == 1) {
      mkl_dtrsm_compact(MKL_ROW_MAJOR, MKL_RIGHT, MKL_UPPER, MKL_TRANS, ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT,
                        m, n, alpha, (const double *)A.data(), A.stride(0), (double *)B.data(), B.stride(0), format,
                        (MKL_INT)vector_type::vector_length);
    } else {
      r_val = -1;
    }
    return r_val;
  }
};
#endif

template <typename ArgDiag>
struct SerialTrsm<Side::Right, Uplo::Upper, Trans::Transpose, ArgDiag, Algo::Trsm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2 && BViewType::rank() == 2);
    // Quick return if possible
    if (B.extent(0) == 0 || B.extent(1) == 0) return 0;

    auto info = KokkosBatched::Impl::checkTrsmInput<Side::Right>(A, B);
    if (info) return info;

    return KokkosBatched::Impl::SerialTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, B.extent(1), B.extent(0), alpha, A.data(), A.stride(0), A.stride(1), B.data(),
        B.stride(1), B.stride(0));
  }
};

template <typename ArgDiag>
struct SerialTrsm<Side::Right, Uplo::Upper, Trans::Transpose, ArgDiag, Algo::Trsm::Blocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2 && BViewType::rank() == 2);
    // Quick return if possible
    if (B.extent(0) == 0 || B.extent(1) == 0) return 0;

    auto info = KokkosBatched::Impl::checkTrsmInput<Side::Right>(A, B);
    if (info) return info;

    return KokkosBatched::Impl::SerialTrsmInternalLeftUpper<Algo::Trsm::Blocked>::invoke(
        ArgDiag::use_unit_diag, false, B.extent(1), B.extent(0), alpha, A.data(), A.stride(0), A.stride(1), B.data(),
        B.stride(1), B.stride(0));
  }
};

///
/// R/L/C
///
/// B := (alpha*B) inv(tril(AH))
/// A(n x n), B(m x n)

#if defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL) && defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_BATCHED) && \
    defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_COMPACT_BATCHED)
template <typename ArgDiag>
struct SerialTrsm<Side::Right, Uplo::Lower, Trans::ConjTranspose, ArgDiag, Algo::Trsm::CompactMKL> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2 && BViewType::rank() == 2);
    typedef typename BViewType::value_type vector_type;
    // typedef typename vector_type::value_type value_type;

    const int m = B.extent(0), n = B.extent(1);

    static_assert(is_vector<vector_type>::value, "value type is not vector type");
    static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                  "AVX, AVX2 and AVX512 is supported");
    const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

    // no error check
    int r_val = 0;
    if (A.stride(0) == 1 && B.stride(0) == 1) {
      mkl_dtrsm_compact(MKL_COL_MAJOR, MKL_RIGHT, MKL_LOWER, MKL_CONJTRANS,
                        ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT, m, n, alpha, (const double *)A.data(),
                        A.stride(1), (double *)B.data(), B.stride(1), format, (MKL_INT)vector_type::vector_length);
    } else if (A.stride(1) == 1 && B.stride(1) == 1) {
      mkl_dtrsm_compact(MKL_ROW_MAJOR, MKL_RIGHT, MKL_LOWER, MKL_CONJTRANS,
                        ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT, m, n, alpha, (const double *)A.data(),
                        A.stride(0), (double *)B.data(), B.stride(0), format, (MKL_INT)vector_type::vector_length);
    } else {
      r_val = -1;
    }
    return r_val;
  }
};
#endif

template <typename ArgDiag>
struct SerialTrsm<Side::Right, Uplo::Lower, Trans::ConjTranspose, ArgDiag, Algo::Trsm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2 && BViewType::rank() == 2);
    // Quick return if possible
    if (B.extent(0) == 0 || B.extent(1) == 0) return 0;

    auto info = KokkosBatched::Impl::checkTrsmInput<Side::Right>(A, B);
    if (info) return info;

    return KokkosBatched::Impl::SerialTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, true, B.extent(1), B.extent(0), alpha, A.data(), A.stride(0), A.stride(1), B.data(),
        B.stride(1), B.stride(0));
  }
};

// [TO DO] ConjTranspose is not supported yet

///
/// R/U/C
///
/// B := (alpha*B) inv(triu(AH))
/// A(n x n), B(m x n)
#if defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL) && defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_BATCHED) && \
    defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_COMPACT_BATCHED)
template <typename ArgDiag>
struct SerialTrsm<Side::Right, Uplo::Upper, Trans::ConjTranspose, ArgDiag, Algo::Trsm::CompactMKL> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2 && BViewType::rank() == 2);
    typedef typename BViewType::value_type vector_type;
    // typedef typename vector_type::value_type value_type;

    const int m = B.extent(0), n = B.extent(1);

    static_assert(is_vector<vector_type>::value, "value type is not vector type");
    static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                  "AVX, AVX2 and AVX512 is supported");
    const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

    // no error check
    int r_val = 0;
    if (A.stride(0) == 1 && B.stride(0) == 1) {
      mkl_dtrsm_compact(MKL_COL_MAJOR, MKL_RIGHT, MKL_UPPER, MKL_CONJTRANS,
                        ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT, m, n, alpha, (const double *)A.data(),
                        A.stride(1), (double *)B.data(), B.stride(1), format, (MKL_INT)vector_type::vector_length);
    } else if (A.stride(1) == 1 && B.stride(1) == 1) {
      mkl_dtrsm_compact(MKL_ROW_MAJOR, MKL_RIGHT, MKL_UPPER, MKL_CONJTRANS,
                        ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT, m, n, alpha, (const double *)A.data(),
                        A.stride(0), (double *)B.data(), B.stride(0), format, (MKL_INT)vector_type::vector_length);
    } else {
      r_val = -1;
    }
    return r_val;
  }
};
#endif

template <typename ArgDiag>
struct SerialTrsm<Side::Right, Uplo::Upper, Trans::ConjTranspose, ArgDiag, Algo::Trsm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    static_assert(AViewType::rank() == 2 && BViewType::rank() == 2);
    // Quick return if possible
    if (B.extent(0) == 0 || B.extent(1) == 0) return 0;

    auto info = KokkosBatched::Impl::checkTrsmInput<Side::Right>(A, B);
    if (info) return info;

    return KokkosBatched::Impl::SerialTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, true, B.extent(1), B.extent(0), alpha, A.data(), A.stride(0), A.stride(1), B.data(),
        B.stride(1), B.stride(0));
  }
};

// [TO DO] ConjTranspose is not supported yet

}  // namespace KokkosBatched

#endif
