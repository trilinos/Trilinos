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
#ifndef __KOKKOSBATCHED_TRSV_SERIAL_IMPL_HPP__
#define __KOKKOSBATCHED_TRSV_SERIAL_IMPL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Trsv_Serial_Internal.hpp"

namespace KokkosBatched {

///
/// Serial Impl
/// ===========

///
/// Implemented:
/// L/NT, U/NT, L/T, U/T
///
/// Not yet implemented
/// L/CT, U/CT

///
/// L/NT
///

#if defined(__KOKKOSBATCHED_ENABLE_INTEL_MKL__) && defined(__KOKKOSBATCHED_ENABLE_INTEL_MKL_BATCHED__) && \
    defined(__KOKKOSBATCHED_ENABLE_INTEL_MKL_COMPACT_BATCHED__)
template <typename ArgDiag>
struct SerialTrsv<Uplo::Lower, Trans::NoTranspose, ArgDiag, Algo::Trsv::CompactMKL> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const bViewType &b) {
    typedef typename bViewType::value_type vector_type;
    // typedef typename vector_type::value_type value_type;

    const int m = b.extent(0), n = 1;

    static_assert(is_vector<vector_type>::value, "value type is not vector type");
    static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                  "AVX, AVX2 and AVX512 is supported");
    const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

    // no error check
    int r_val = 0;
    if (A.stride_0() == 1) {
      mkl_dtrsm_compact(MKL_COL_MAJOR, MKL_LEFT, MKL_LOWER, MKL_NOTRANS,
                        ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT, m, n, alpha, (const double *)A.data(),
                        A.stride_0(), (double *)b.data(), b.stride_0(), format, (MKL_INT)vector_type::vector_length);
    } else if (A.stride_1() == 1) {
      mkl_dtrsm_compact(MKL_ROW_MAJOR, MKL_LEFT, MKL_LOWER, MKL_NOTRANS,
                        ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT, m, n, alpha, (const double *)A.data(),
                        A.stride_0(), (double *)b.data(), b.stride_0(), format, (MKL_INT)vector_type::vector_length);
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
    return SerialTrsvInternalLower<Algo::Trsv::Unblocked>::invoke(ArgDiag::use_unit_diag, A.extent(0), alpha, A.data(),
                                                                  A.stride_0(), A.stride_1(), b.data(), b.stride_0());
  }
};

template <typename ArgDiag>
struct SerialTrsv<Uplo::Lower, Trans::NoTranspose, ArgDiag, Algo::Trsv::Blocked> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const bViewType &b) {
    return SerialTrsvInternalLower<Algo::Trsv::Blocked>::invoke(ArgDiag::use_unit_diag, A.extent(0), alpha, A.data(),
                                                                A.stride_0(), A.stride_1(), b.data(), b.stride_0());
  }
};

///
/// L/T
///

#if defined(__KOKKOSBATCHED_ENABLE_INTEL_MKL__) && defined(__KOKKOSBATCHED_ENABLE_INTEL_MKL_BATCHED__) && \
    defined(__KOKKOSBATCHED_ENABLE_INTEL_MKL_COMPACT_BATCHED__)
template <typename ArgDiag>
struct SerialTrsv<Uplo::Lower, Trans::Transpose, ArgDiag, Algo::Trsv::CompactMKL> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const bViewType &b) {
    typedef typename bViewType::value_type vector_type;
    // typedef typename vector_type::value_type value_type;

    const int m = b.extent(0), n = 1;

    static_assert(is_vector<vector_type>::value, "value type is not vector type");
    static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                  "AVX, AVX2 and AVX512 is supported");
    const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

    // no error check
    int r_val = 0;
    if (A.stride_0() == 1) {
      mkl_dtrsm_compact(MKL_COL_MAJOR, MKL_LEFT, MKL_LOWER, MKL_TRANS, ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT,
                        m, n, alpha, (const double *)A.data(), A.stride_0(), (double *)b.data(), b.stride_0(), format,
                        (MKL_INT)vector_type::vector_length);
    } else if (A.stride_1() == 1) {
      mkl_dtrsm_compact(MKL_ROW_MAJOR, MKL_LEFT, MKL_LOWER, MKL_TRANS, ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT,
                        m, n, alpha, (const double *)A.data(), A.stride_0(), (double *)b.data(), b.stride_0(), format,
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
    return SerialTrsvInternalUpper<Algo::Trsv::Unblocked>::invoke(ArgDiag::use_unit_diag, A.extent(1), alpha, A.data(),
                                                                  A.stride_1(), A.stride_0(), b.data(), b.stride_0());
  }
};

template <typename ArgDiag>
struct SerialTrsv<Uplo::Lower, Trans::Transpose, ArgDiag, Algo::Trsv::Blocked> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const bViewType &b) {
    return SerialTrsvInternalUpper<Algo::Trsv::Blocked>::invoke(ArgDiag::use_unit_diag, A.extent(1), alpha, A.data(),
                                                                A.stride_1(), A.stride_0(), b.data(), b.stride_0());
  }
};

///
/// U/NT
///

#if defined(__KOKKOSBATCHED_ENABLE_INTEL_MKL__) && defined(__KOKKOSBATCHED_ENABLE_INTEL_MKL_BATCHED__) && \
    defined(__KOKKOSBATCHED_ENABLE_INTEL_MKL_COMPACT_BATCHED__)
template <typename ArgDiag>
struct SerialTrsv<Uplo::Upper, Trans::NoTranspose, ArgDiag, Algo::Trsv::CompactMKL> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const bViewType &b) {
    typedef typename bViewType::value_type vector_type;
    // typedef typename vector_type::value_type value_type;

    const int m = b.extent(0), n = 1;

    static_assert(is_vector<vector_type>::value, "value type is not vector type");
    static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                  "AVX, AVX2 and AVX512 is supported");
    const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

    // no error check
    int r_val = 0;
    if (A.stride_0() == 1) {
      mkl_dtrsm_compact(MKL_COL_MAJOR, MKL_LEFT, MKL_UPPER, MKL_NOTRANS,
                        ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT, m, n, alpha, (const double *)A.data(),
                        A.stride_0(), (double *)b.data(), b.stride_0(), format, (MKL_INT)vector_type::vector_length);
    } else if (A.stride_1() == 1) {
      mkl_dtrsm_compact(MKL_ROW_MAJOR, MKL_LEFT, MKL_UPPER, MKL_NOTRANS,
                        ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT, m, n, alpha, (const double *)A.data(),
                        A.stride_0(), (double *)b.data(), b.stride_0(), format, (MKL_INT)vector_type::vector_length);
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
    return SerialTrsvInternalUpper<Algo::Trsv::Unblocked>::invoke(ArgDiag::use_unit_diag, A.extent(0), alpha, A.data(),
                                                                  A.stride_0(), A.stride_1(), b.data(), b.stride_0());
  }
};

template <typename ArgDiag>
struct SerialTrsv<Uplo::Upper, Trans::NoTranspose, ArgDiag, Algo::Trsv::Blocked> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const bViewType &b) {
    return SerialTrsvInternalUpper<Algo::Trsv::Blocked>::invoke(ArgDiag::use_unit_diag, A.extent(0), alpha, A.data(),
                                                                A.stride_0(), A.stride_1(), b.data(), b.stride_0());
  }
};

///
/// U/T
///

#if defined(__KOKKOSBATCHED_ENABLE_INTEL_MKL__) && defined(__KOKKOSBATCHED_ENABLE_INTEL_MKL_BATCHED__) && \
    defined(__KOKKOSBATCHED_ENABLE_INTEL_MKL_COMPACT_BATCHED__)
template <typename ArgDiag>
struct SerialTrsv<Uplo::Upper, Trans::Transpose, ArgDiag, Algo::Trsv::CompactMKL> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const bViewType &b) {
    typedef typename bViewType::value_type vector_type;
    // typedef typename vector_type::value_type value_type;

    const int m = b.extent(0), n = 1;

    static_assert(is_vector<vector_type>::value, "value type is not vector type");
    static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                  "AVX, AVX2 and AVX512 is supported");
    const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

    // no error check
    int r_val = 0;
    if (A.stride_0() == 1) {
      mkl_dtrsm_compact(MKL_COL_MAJOR, MKL_LEFT, MKL_UPPER, MKL_TRANS, ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT,
                        m, n, alpha, (const double *)A.data(), A.stride_0(), (double *)b.data(), b.stride_0(), format,
                        (MKL_INT)vector_type::vector_length);
    } else if (A.stride_1() == 1) {
      mkl_dtrsm_compact(MKL_ROW_MAJOR, MKL_LEFT, MKL_UPPER, MKL_TRANS, ArgDiag::use_unit_diag ? MKL_UNIT : MKL_NONUNIT,
                        m, n, alpha, (const double *)A.data(), A.stride_0(), (double *)b.data(), b.stride_0(), format,
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
    return SerialTrsvInternalLower<Algo::Trsv::Unblocked>::invoke(ArgDiag::use_unit_diag, A.extent(1), alpha, A.data(),
                                                                  A.stride_1(), A.stride_0(), b.data(), b.stride_0());
  }
};

template <typename ArgDiag>
struct SerialTrsv<Uplo::Upper, Trans::Transpose, ArgDiag, Algo::Trsv::Blocked> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const bViewType &b) {
    return SerialTrsvInternalLower<Algo::Trsv::Blocked>::invoke(ArgDiag::use_unit_diag, A.extent(1), alpha, A.data(),
                                                                A.stride_1(), A.stride_0(), b.data(), b.stride_0());
  }
};

}  // namespace KokkosBatched

#endif
