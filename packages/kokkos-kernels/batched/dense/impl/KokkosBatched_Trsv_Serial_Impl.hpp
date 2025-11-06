// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_TRSV_SERIAL_IMPL_HPP
#define KOKKOSBATCHED_TRSV_SERIAL_IMPL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

#include <Kokkos_DynRankView.hpp>
#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Trsv_Serial_Internal.hpp"

namespace KokkosBatched {
namespace Impl {
template <typename AViewType, typename bViewType>
KOKKOS_INLINE_FUNCTION static int checkTrsvInput([[maybe_unused]] const AViewType &A,
                                                 [[maybe_unused]] const bViewType &b) {
  static_assert(Kokkos::is_view_v<AViewType> || Kokkos::is_dyn_rank_view_v<AViewType>,
                "KokkosBatched::trsv: AViewType must be either a Kokkos::View or a Kokkos::DynRankView.");
  static_assert(Kokkos::is_view_v<bViewType> || Kokkos::is_dyn_rank_view_v<bViewType>,
                "KokkosBatched::trsv: bViewType must be either a Kokkos::View or a Kokkos::DynRankView.");
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  if (A.rank() != 2) {
    Kokkos::printf(
        "KokkosBatched::trsv: A must be a rank 2 View."
        "A.rank() = %d\n",
        A.rank());
    return 1;
  }

  if (b.rank() != 1) {
    Kokkos::printf(
        "KokkosBatched::trsv: b must be a rank 1 View."
        "b.rank() = %d\n",
        b.rank());
    return 1;
  }

  // FIXME : check leading dimension is suppressed for now
  //         because of the compatibility issue with Trilinos
  // const int lda = A.extent(0), n = A.extent(1);
  // if (lda < Kokkos::max(1, n)) {
  //   Kokkos::printf(
  //       "KokkosBatched::trsv: leading dimension of A must not be smaller than "
  //       "max(1, n): "
  //       "lda = %d, n = %d\n",
  //       lda, n);
  //   return 1;
  // }

#endif
  return 0;
}
}  // namespace Impl

//// Lower non-transpose ////
#if defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL) && defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_BATCHED) && \
    defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_COMPACT_BATCHED)
template <typename ArgDiag>
struct SerialTrsv<Uplo::Lower, Trans::NoTranspose, ArgDiag, Algo::Trsv::CompactMKL> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const bViewType &b) {
    // Quick return if possible
    // if (A.extent(1) == 0) return 0;

    auto info = KokkosBatched::Impl::checkTrsvInput(A, b);
    if (info) return info;

    using vector_type = typename bViewType::value_type;
    const int m = b.extent(0), n = 1;

    static_assert(is_vector<vector_type>::value, "value type is not vector type");
#if defined(KOKKOS_ARCH_AVX512XEON)
    static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                  "AVX, AVX2 and AVX512 is supported");
#else
    static_assert(vector_type::vector_length == 4, "AVX and AVX2 is supported");
#endif
    const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

    // no error check
    int r_val = 0;
    if (A.stride(0) == 1) {
      mkl_dtrsm_compact(MKL_COL_MAJOR, MKL_LEFT, MKL_LOWER, MKL_NOTRANS,
                        ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT, m, n, alpha, (const double *)A.data(),
                        A.stride(0), (double *)b.data(), b.stride(0), format, (MKL_INT)vector_type::vector_length);
    } else if (A.stride(1) == 1) {
      mkl_dtrsm_compact(MKL_ROW_MAJOR, MKL_LEFT, MKL_LOWER, MKL_NOTRANS,
                        ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT, m, n, alpha, (const double *)A.data(),
                        A.stride(0), (double *)b.data(), b.stride(0), format, (MKL_INT)vector_type::vector_length);
    } else {
      r_val = -1;
    }
    return r_val;
  }
};
#endif

template <typename ArgDiag>
struct SerialTrsv<Uplo::Lower, Trans::NoTranspose, ArgDiag, Algo::Trsv::Unblocked> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const bViewType &b) {
    // Quick return if possible
    // if (A.extent(1) == 0) return 0;

    auto info = KokkosBatched::Impl::checkTrsvInput(A, b);
    if (info) return info;
    return KokkosBatched::Impl::SerialTrsvInternalLower<Algo::Trsv::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, A.extent(0), alpha, A.data(), A.stride(0), A.stride(1), b.data(), b.stride(0));
  }
};

template <typename ArgDiag>
struct SerialTrsv<Uplo::Lower, Trans::NoTranspose, ArgDiag, Algo::Trsv::Blocked> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const bViewType &b) {
    // Quick return if possible
    // if (A.extent(1) == 0) return 0;

    auto info = KokkosBatched::Impl::checkTrsvInput(A, b);
    if (info) return info;
    return KokkosBatched::Impl::SerialTrsvInternalLower<Algo::Trsv::Blocked>::invoke(
        ArgDiag::use_unit_diag, false, A.extent(0), alpha, A.data(), A.stride(0), A.stride(1), b.data(), b.stride(0));
  }
};

//// Lower transpose ////
#if defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL) && defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_BATCHED) && \
    defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_COMPACT_BATCHED)
template <typename ArgDiag>
struct SerialTrsv<Uplo::Lower, Trans::Transpose, ArgDiag, Algo::Trsv::CompactMKL> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const bViewType &b) {
    // Quick return if possible
    // if (A.extent(1) == 0) return 0;

    auto info = KokkosBatched::Impl::checkTrsvInput(A, b);
    if (info) return info;

    using vector_type = typename bViewType::value_type;
    const int m = b.extent(0), n = 1;

    static_assert(is_vector<vector_type>::value, "value type is not vector type");
#if defined(KOKKOS_ARCH_AVX512XEON)
    static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                  "AVX, AVX2 and AVX512 is supported");
#else
    static_assert(vector_type::vector_length == 4, "AVX and AVX2 is supported");
#endif
    const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

    // no error check
    int r_val = 0;
    if (A.stride(0) == 1) {
      mkl_dtrsm_compact(MKL_COL_MAJOR, MKL_LEFT, MKL_LOWER, MKL_TRANS, ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT,
                        m, n, alpha, (const double *)A.data(), A.stride(0), (double *)b.data(), b.stride(0), format,
                        (MKL_INT)vector_type::vector_length);
    } else if (A.stride(1) == 1) {
      mkl_dtrsm_compact(MKL_ROW_MAJOR, MKL_LEFT, MKL_LOWER, MKL_TRANS, ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT,
                        m, n, alpha, (const double *)A.data(), A.stride(0), (double *)b.data(), b.stride(0), format,
                        (MKL_INT)vector_type::vector_length);
    } else {
      r_val = -1;
    }
    return r_val;
  }
};
#endif

template <typename ArgDiag>
struct SerialTrsv<Uplo::Lower, Trans::Transpose, ArgDiag, Algo::Trsv::Unblocked> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const bViewType &b) {
    // Quick return if possible
    // if (A.extent(1) == 0) return 0;

    auto info = KokkosBatched::Impl::checkTrsvInput(A, b);
    if (info) return info;
    return KokkosBatched::Impl::SerialTrsvInternalUpper<Algo::Trsv::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, A.extent(1), alpha, A.data(), A.stride(1), A.stride(0), b.data(), b.stride(0));
  }
};

template <typename ArgDiag>
struct SerialTrsv<Uplo::Lower, Trans::Transpose, ArgDiag, Algo::Trsv::Blocked> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const bViewType &b) {
    // Quick return if possible
    // if (A.extent(1) == 0) return 0;

    auto info = KokkosBatched::Impl::checkTrsvInput(A, b);
    if (info) return info;
    return KokkosBatched::Impl::SerialTrsvInternalUpper<Algo::Trsv::Blocked>::invoke(
        ArgDiag::use_unit_diag, false, A.extent(1), alpha, A.data(), A.stride(1), A.stride(0), b.data(), b.stride(0));
  }
};

//// Lower conjugate-transpose ////
#if defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL) && defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_BATCHED) && \
    defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_COMPACT_BATCHED)
template <typename ArgDiag>
struct SerialTrsv<Uplo::Lower, Trans::ConjTranspose, ArgDiag, Algo::Trsv::CompactMKL> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const bViewType &b) {
    // Quick return if possible
    // if (A.extent(1) == 0) return 0;

    auto info = KokkosBatched::Impl::checkTrsvInput(A, b);
    if (info) return info;

    using vector_type = typename bViewType::value_type;
    const int m = b.extent(0), n = 1;

    static_assert(is_vector<vector_type>::value, "value type is not vector type");
#if defined(KOKKOS_ARCH_AVX512XEON)
    static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                  "AVX, AVX2 and AVX512 is supported");
#else
    static_assert(vector_type::vector_length == 4, "AVX and AVX2 is supported");
#endif
    const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

    // no error check
    int r_val = 0;
    if (A.stride(0) == 1) {
      mkl_dtrsm_compact(MKL_COL_MAJOR, MKL_LEFT, MKL_LOWER, MKL_CONJTRANS,
                        ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT, m, n, alpha, (const double *)A.data(),
                        A.stride(0), (double *)b.data(), b.stride(0), format, (MKL_INT)vector_type::vector_length);
    } else if (A.stride(1) == 1) {
      mkl_dtrsm_compact(MKL_ROW_MAJOR, MKL_LEFT, MKL_LOWER, MKL_CONJTRANS,
                        ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT, m, n, alpha, (const double *)A.data(),
                        A.stride(0), (double *)b.data(), b.stride(0), format, (MKL_INT)vector_type::vector_length);
    } else {
      r_val = -1;
    }
    return r_val;
  }
};
#endif

template <typename ArgDiag>
struct SerialTrsv<Uplo::Lower, Trans::ConjTranspose, ArgDiag, Algo::Trsv::Unblocked> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const bViewType &b) {
    // Quick return if possible
    // if (A.extent(1) == 0) return 0;

    auto info = KokkosBatched::Impl::checkTrsvInput(A, b);
    if (info) return info;
    return KokkosBatched::Impl::SerialTrsvInternalUpper<Algo::Trsv::Unblocked>::invoke(
        ArgDiag::use_unit_diag, true, A.extent(1), alpha, A.data(), A.stride(1), A.stride(0), b.data(), b.stride(0));
  }
};

template <typename ArgDiag>
struct SerialTrsv<Uplo::Lower, Trans::ConjTranspose, ArgDiag, Algo::Trsv::Blocked> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const bViewType &b) {
    // Quick return if possible
    // if (A.extent(1) == 0) return 0;

    auto info = KokkosBatched::Impl::checkTrsvInput(A, b);
    if (info) return info;
    return KokkosBatched::Impl::SerialTrsvInternalUpper<Algo::Trsv::Blocked>::invoke(
        ArgDiag::use_unit_diag, true, A.extent(1), alpha, A.data(), A.stride(1), A.stride(0), b.data(), b.stride(0));
  }
};

//// Upper non-transpose ////
#if defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL) && defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_BATCHED) && \
    defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_COMPACT_BATCHED)
template <typename ArgDiag>
struct SerialTrsv<Uplo::Upper, Trans::NoTranspose, ArgDiag, Algo::Trsv::CompactMKL> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const bViewType &b) {
    // Quick return if possible
    // if (A.extent(1) == 0) return 0;

    auto info = KokkosBatched::Impl::checkTrsvInput(A, b);
    if (info) return info;

    using vector_type = typename bViewType::value_type;
    const int m = b.extent(0), n = 1;

    static_assert(is_vector<vector_type>::value, "value type is not vector type");
#if defined(KOKKOS_ARCH_AVX512XEON)
    static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                  "AVX, AVX2 and AVX512 is supported");
#else
    static_assert(vector_type::vector_length == 4, "AVX and AVX2 is supported");
#endif
    const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

    // no error check
    int r_val = 0;
    if (A.stride(0) == 1) {
      mkl_dtrsm_compact(MKL_COL_MAJOR, MKL_LEFT, MKL_UPPER, MKL_NOTRANS,
                        ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT, m, n, alpha, (const double *)A.data(),
                        A.stride(0), (double *)b.data(), b.stride(0), format, (MKL_INT)vector_type::vector_length);
    } else if (A.stride(1) == 1) {
      mkl_dtrsm_compact(MKL_ROW_MAJOR, MKL_LEFT, MKL_UPPER, MKL_NOTRANS,
                        ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT, m, n, alpha, (const double *)A.data(),
                        A.stride(0), (double *)b.data(), b.stride(0), format, (MKL_INT)vector_type::vector_length);
    } else {
      r_val = -1;
    }
    return r_val;
  }
};
#endif

template <typename ArgDiag>
struct SerialTrsv<Uplo::Upper, Trans::NoTranspose, ArgDiag, Algo::Trsv::Unblocked> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const bViewType &b) {
    // Quick return if possible
    // if (A.extent(1) == 0) return 0;

    auto info = KokkosBatched::Impl::checkTrsvInput(A, b);
    if (info) return info;
    return KokkosBatched::Impl::SerialTrsvInternalUpper<Algo::Trsv::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, A.extent(0), alpha, A.data(), A.stride(0), A.stride(1), b.data(), b.stride(0));
  }
};

template <typename ArgDiag>
struct SerialTrsv<Uplo::Upper, Trans::NoTranspose, ArgDiag, Algo::Trsv::Blocked> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const bViewType &b) {
    // Quick return if possible
    // if (A.extent(1) == 0) return 0;

    auto info = KokkosBatched::Impl::checkTrsvInput(A, b);
    if (info) return info;
    return KokkosBatched::Impl::SerialTrsvInternalUpper<Algo::Trsv::Blocked>::invoke(
        ArgDiag::use_unit_diag, false, A.extent(0), alpha, A.data(), A.stride(0), A.stride(1), b.data(), b.stride(0));
  }
};

//// Upper transpose ////
#if defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL) && defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_BATCHED) && \
    defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_COMPACT_BATCHED)
template <typename ArgDiag>
struct SerialTrsv<Uplo::Upper, Trans::Transpose, ArgDiag, Algo::Trsv::CompactMKL> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const bViewType &b) {
    // Quick return if possible
    // if (A.extent(1) == 0) return 0;

    auto info = KokkosBatched::Impl::checkTrsvInput(A, b);
    if (info) return info;

    using vector_type = typename bViewType::value_type;
    const int m = b.extent(0), n = 1;

    static_assert(is_vector<vector_type>::value, "value type is not vector type");
#if defined(KOKKOS_ARCH_AVX512XEON)
    static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                  "AVX, AVX2 and AVX512 is supported");
#else
    static_assert(vector_type::vector_length == 4, "AVX and AVX2 is supported");
#endif
    const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

    // no error check
    int r_val = 0;
    if (A.stride(0) == 1) {
      mkl_dtrsm_compact(MKL_COL_MAJOR, MKL_LEFT, MKL_UPPER, MKL_TRANS, ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT,
                        m, n, alpha, (const double *)A.data(), A.stride(0), (double *)b.data(), b.stride(0), format,
                        (MKL_INT)vector_type::vector_length);
    } else if (A.stride(1) == 1) {
      mkl_dtrsm_compact(MKL_ROW_MAJOR, MKL_LEFT, MKL_UPPER, MKL_TRANS, ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT,
                        m, n, alpha, (const double *)A.data(), A.stride(0), (double *)b.data(), b.stride(0), format,
                        (MKL_INT)vector_type::vector_length);
    } else {
      r_val = -1;
    }
    return r_val;
  }
};
#endif

template <typename ArgDiag>
struct SerialTrsv<Uplo::Upper, Trans::Transpose, ArgDiag, Algo::Trsv::Unblocked> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const bViewType &b) {
    // Quick return if possible
    // if (A.extent(1) == 0) return 0;

    auto info = KokkosBatched::Impl::checkTrsvInput(A, b);
    if (info) return info;
    return KokkosBatched::Impl::SerialTrsvInternalLower<Algo::Trsv::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, A.extent(1), alpha, A.data(), A.stride(1), A.stride(0), b.data(), b.stride(0));
  }
};

template <typename ArgDiag>
struct SerialTrsv<Uplo::Upper, Trans::Transpose, ArgDiag, Algo::Trsv::Blocked> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const bViewType &b) {
    // Quick return if possible
    // if (A.extent(1) == 0) return 0;

    auto info = KokkosBatched::Impl::checkTrsvInput(A, b);
    if (info) return info;
    return KokkosBatched::Impl::SerialTrsvInternalLower<Algo::Trsv::Blocked>::invoke(
        ArgDiag::use_unit_diag, false, A.extent(1), alpha, A.data(), A.stride(1), A.stride(0), b.data(), b.stride(0));
  }
};

//// Upper conjugate-transpose ////
#if defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL) && defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_BATCHED) && \
    defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_COMPACT_BATCHED)
template <typename ArgDiag>
struct SerialTrsv<Uplo::Upper, Trans::ConjTranspose, ArgDiag, Algo::Trsv::CompactMKL> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const bViewType &b) {
    // Quick return if possible
    // if (A.extent(1) == 0) return 0;

    auto info = KokkosBatched::Impl::checkTrsvInput(A, b);
    if (info) return info;

    using vector_type = typename bViewType::value_type;
    const int m = b.extent(0), n = 1;

    static_assert(is_vector<vector_type>::value, "value type is not vector type");
#if defined(KOKKOS_ARCH_AVX512XEON)
    static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                  "AVX, AVX2 and AVX512 is supported");
#else
    static_assert(vector_type::vector_length == 4, "AVX and AVX2 is supported");
#endif
    const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

    // no error check
    int r_val = 0;
    if (A.stride(0) == 1) {
      mkl_dtrsm_compact(MKL_COL_MAJOR, MKL_LEFT, MKL_UPPER, MKL_CONJTRANS,
                        ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT, m, n, alpha, (const double *)A.data(),
                        A.stride(0), (double *)b.data(), b.stride(0), format, (MKL_INT)vector_type::vector_length);
    } else if (A.stride(1) == 1) {
      mkl_dtrsm_compact(MKL_ROW_MAJOR, MKL_LEFT, MKL_UPPER, MKL_CONJTRANS,
                        ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT, m, n, alpha, (const double *)A.data(),
                        A.stride(0), (double *)b.data(), b.stride(0), format, (MKL_INT)vector_type::vector_length);
    } else {
      r_val = -1;
    }
    return r_val;
  }
};
#endif

template <typename ArgDiag>
struct SerialTrsv<Uplo::Upper, Trans::ConjTranspose, ArgDiag, Algo::Trsv::Unblocked> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const bViewType &b) {
    // Quick return if possible
    // if (A.extent(1) == 0) return 0;

    auto info = KokkosBatched::Impl::checkTrsvInput(A, b);
    if (info) return info;
    return KokkosBatched::Impl::SerialTrsvInternalLower<Algo::Trsv::Unblocked>::invoke(
        ArgDiag::use_unit_diag, true, A.extent(1), alpha, A.data(), A.stride(1), A.stride(0), b.data(), b.stride(0));
  }
};

template <typename ArgDiag>
struct SerialTrsv<Uplo::Upper, Trans::ConjTranspose, ArgDiag, Algo::Trsv::Blocked> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const bViewType &b) {
    // Quick return if possible
    // if (A.extent(1) == 0) return 0;

    auto info = KokkosBatched::Impl::checkTrsvInput(A, b);
    if (info) return info;
    return KokkosBatched::Impl::SerialTrsvInternalLower<Algo::Trsv::Blocked>::invoke(
        ArgDiag::use_unit_diag, true, A.extent(1), alpha, A.data(), A.stride(1), A.stride(0), b.data(), b.stride(0));
  }
};

}  // namespace KokkosBatched

#endif
