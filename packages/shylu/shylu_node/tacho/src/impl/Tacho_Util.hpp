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
#ifndef __TACHO_UTIL_HPP__
#define __TACHO_UTIL_HPP__

// standard C includes
#include <stdio.h>
#include <string.h>

// "std" includes
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include <cmath>
#include <complex>

#include <limits>

#include "Kokkos_Core.hpp"
#include "Kokkos_Timer.hpp"

#if defined(__INTEL_MKL__)
#include "mkl.h"
#endif

#if defined(KOKKOS_ENABLE_CUDA)
#include "cublas_v2.h"
#include "cusolverDn.h"
#endif

#if defined(KOKKOS_ENABLE_HIP)
#include "rocblas/rocblas.h"
#include "rocsolver/rocsolver.h"
#endif

#include "Tacho.hpp"

/// \file Tacho_Util.hpp
/// \brief Utility functions and constant integer class like an enum class.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

const char *Version();

///
/// error macros
///
// #define MSG_NOT_YET_IMPLEMENTED(what) "Not yet implemented: " #what
// #define MSG_INVALID_INPUT(what) "Invaid input argument: " #what
#define MSG_NOT_HAVE_PACKAGE(what) "Tacho does not have a package or library: " #what
#define MSG_INVALID_TEMPLATE_ARGS "Invaid template arguments"
#define MSG_INVALID_INPUT "Invaid input arguments"
#define MSG_NOT_IMPLEMENTED "Not yet implemented"
#if !defined(KOKKOS_ENABLE_CUDA) && !defined(KOKKOS_ENABLE_HIP)
KOKKOS_FUNCTION constexpr bool runsOnCudaOrHIP() { return false; }
#else
KOKKOS_FUNCTION constexpr bool runsOnCudaOrHIP() {
  KOKKOS_IF_ON_HOST(return false;)
  KOKKOS_IF_ON_DEVICE(return true;)
}
#endif
#if !defined(KOKKOS_ENABLE_OPENMP)
KOKKOS_FUNCTION constexpr bool runsWithOMP() { return false; }
#else
KOKKOS_FUNCTION constexpr bool runsWithOMP() {
  KOKKOS_IF_ON_HOST(return true;)
  KOKKOS_IF_ON_DEVICE(return false;)
}
#endif
template <class ExecutionSpace>
inline constexpr bool run_tacho_on_host_v = !std::is_same_v<ExecutionSpace, Kokkos::DefaultExecutionSpace>
                                          || std::is_same_v<Kokkos::DefaultExecutionSpace, Kokkos::DefaultHostExecutionSpace>;

#define TACHO_TEST_FOR_ABORT(ierr, msg)                                                                                \
  if ((ierr) != 0) {                                                                                                   \
    Kokkos::printf(">> Error in file %s, line %d, error %d \n   %s\n", __FILE__, __LINE__, ierr, msg);                         \
    Kokkos::abort(">> Tacho abort\n");                                                                                 \
  }

#define TACHO_TEST_FOR_EXCEPTION(ierr, x, msg)                                                                         \
  if ((ierr) != 0) {                                                                                                   \
    fprintf(stderr, ">> Error in file %s, line %d, error %d \n", __FILE__, __LINE__, ierr);                            \
    fprintf(stderr, "   %s\n", msg);                                                                                   \
    throw x(msg);                                                                                                      \
  }

#if defined(KOKKOS_ENABLE_ASM)
#if defined(__amd64) || defined(__amd64__) || defined(__x86_64) || defined(__x86_64__)
#if !defined(_WIN32) /* IS NOT Microsoft Windows */
#define TACHO_IMPL_PAUSE asm volatile("pause\n" ::: "memory");
#else
#define TACHO_IMPL_PAUSE __asm__ __volatile__("pause\n" ::: "memory");
#endif
#elif defined(__PPC64__)
#define TACHO_IMPL_PAUSE asm volatile("or 27, 27, 27" ::: "memory");
#endif
#else
#define TACHO_IMPL_PAUSE
#endif

///
/// label size used to identify object name
///
enum : int { MaxDependenceSize = 4, ThresholdSolvePhaseUsingBlas3 = 12, CudaVectorSize = 4 };

///
/// util
///
template <typename Ta, typename Tb> KOKKOS_FORCEINLINE_FUNCTION static Ta min(const Ta a, const Tb b) {
  return (a < static_cast<Ta>(b) ? a : static_cast<Ta>(b));
}

template <typename Ta, typename Tb> KOKKOS_FORCEINLINE_FUNCTION static Ta max(const Ta a, const Tb b) {
  return (a > static_cast<Ta>(b) ? a : static_cast<Ta>(b));
}

template <typename Ta, typename Tb> KOKKOS_FORCEINLINE_FUNCTION static void swap(Ta &a, Tb &b) {
  Ta c(a);
  a = static_cast<Ta>(b);
  b = static_cast<Tb>(c);
}

KOKKOS_FORCEINLINE_FUNCTION
static void clear(char *buf, size_type bufsize) {
  KOKKOS_IF_ON_HOST(( memset(buf, 0, bufsize); ))
  KOKKOS_IF_ON_DEVICE((
    for (size_type i = 0; i < bufsize; ++i)
      buf[i] = 0;
  ))
}

template <typename MemberType>
KOKKOS_FORCEINLINE_FUNCTION static void clear(MemberType &member, char *buf, size_type bufsize) {
  KOKKOS_IF_ON_HOST(( memset(buf, 0, bufsize); ))
  KOKKOS_IF_ON_DEVICE((
    const ordinal_type team_index_range = (bufsize / CudaVectorSize) + (bufsize % CudaVectorSize > 0);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, team_index_range), [&](const int &idx) {
      const int ioff = idx * CudaVectorSize;
      const int itmp = bufsize - ioff;
      const int icnt = itmp > CudaVectorSize ? CudaVectorSize : itmp;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, icnt), [&](const int &ii) {
        const int i = ioff + ii;
        buf[i] = 0;
      });
    });
  ))
}

template <typename T1, typename T2, typename CompareType>
KOKKOS_INLINE_FUNCTION static T1 *lower_bound(T1 *first, T1 *last, const T2 &val, CompareType compare) {
  T1 *it;
  ordinal_type step = 0, count = last - first;
  while (count > 0) {
    it = first;
    step = count / 2;
    it += step;
    if (compare(*it, val)) {
      first = ++it;
      count -= step + 1;
    } else {
      count = step;
    }
  }
  return first;
}

template <size_t BufSize, typename SpaceType = Kokkos::DefaultExecutionSpace> struct Flush {
  typedef double value_type;

  // flush a large host buffer
  Kokkos::View<value_type *, SpaceType> _buf;
  Flush() : _buf("Flush::buf", BufSize / sizeof(double)) { Kokkos::deep_copy(_buf, 1); }

  KOKKOS_INLINE_FUNCTION
  void init(value_type &update) { update = 0; }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type &update, const volatile value_type &input) { update += input; }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, value_type &update) const { update += _buf[i]; }

  void run() {
    double sum = 0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<SpaceType>(0, BufSize / sizeof(double)), *this, sum);
    SpaceType().fence();
    FILE *fp = fopen("/dev/null", "w");
    fprintf(fp, "%f\n", sum);
    fclose(fp);
  }
};

template <typename ValueType> struct Random;

template <> struct Random<double> {
  Random(const unsigned int seed = 0) { srand(seed); }
  double value() { return rand() / ((double)RAND_MAX + 1.0); }
};

template <> struct Random<std::complex<double>> {
  Random(const unsigned int seed = 0) { srand(seed); }
  std::complex<double> value() {
    return std::complex<double>(rand() / ((double)RAND_MAX + 1.0), rand() / ((double)RAND_MAX + 1.0));
  }
};

template <> struct Random<Kokkos::complex<double>> {
  Random(const unsigned int seed = 0) { srand(seed); }
  Kokkos::complex<double> value() {
    return Kokkos::complex<double>(rand() / ((double)RAND_MAX + 1.0), rand() / ((double)RAND_MAX + 1.0));
  }
};

///
/// Tag struct
///
struct NullTag {
  enum : int { tag = 0 };
};
struct PivotMode {
  struct Flame {};  /// 0 base relative pivot index
  struct Lapack {}; /// 1 base index
};

struct Partition {
  enum : int { Top = 101, Bottom, Left = 201, Right, TopLeft = 301, TopRight, BottomLeft, BottomRight };
};
template <int T> struct is_valid_partition_tag {
  enum : bool {
    value =
        (T == Partition::TopLeft || T == Partition::Top || T == Partition::TopRight || T == Partition::Left ||
         T == Partition::Right || T == Partition::BottomLeft || T == Partition::Bottom || T == Partition::BottomRight)
  };
};

struct Uplo {
  enum : int { tag = 400 };
  struct Upper {
    enum : int { tag = 401 };
    static constexpr char param = 'U';
#if defined(__INTEL_MKL__)
    static constexpr CBLAS_UPLO mkl_param = CblasUpper;
#endif
#if defined(CUBLAS_VERSION)
    static constexpr cublasFillMode_t cublas_param = CUBLAS_FILL_MODE_UPPER;
#endif
#if defined(ROCBLAS_VERSION_MAJOR)
    static constexpr rocblas_fill rocblas_param = rocblas_fill_upper;
#endif
  };
  struct Lower {
    enum : int { tag = 402 };
    static constexpr char param = 'L';
#if defined(__INTEL_MKL__)
    static constexpr CBLAS_UPLO mkl_param = CblasLower;
#endif
#if defined(CUBLAS_VERSION)
    static constexpr cublasFillMode_t cublas_param = CUBLAS_FILL_MODE_LOWER;
#endif
#if defined(ROCBLAS_VERSION_MAJOR)
    static constexpr rocblas_fill rocblas_param = rocblas_fill_lower;
#endif
  };
};
template <typename T> struct is_valid_uplo_tag {
  enum : bool { value = (std::is_same<T, Uplo::Upper>::value || std::is_same<T, Uplo::Lower>::value) };
};
template <typename T> struct transpose_uplo_tag;
template <> struct transpose_uplo_tag<Uplo::Lower> { typedef Uplo::Upper type; };
template <> struct transpose_uplo_tag<Uplo::Upper> { typedef Uplo::Lower type; };

struct Side {
  enum : int { tag = 500 };
  struct Left {
    enum : int { tag = 501 };
    static constexpr char param = 'L';
#if defined(__INTEL_MKL__)
    static constexpr CBLAS_SIDE mkl_param = CblasLeft;
#endif
#if defined(CUBLAS_VERSION)
    static constexpr cublasSideMode_t cublas_param = CUBLAS_SIDE_LEFT;
#endif
#if defined(ROCBLAS_VERSION_MAJOR)
    static constexpr rocblas_side rocblas_param = rocblas_side_left;
#endif
  };
  struct Right {
    enum : int { tag = 502 };
    static constexpr char param = 'R';
#if defined(__INTEL_MKL__)
    static constexpr CBLAS_SIDE mkl_param = CblasRight;
#endif
#if defined(CUBLAS_VERSION)
    static constexpr cublasSideMode_t cublas_param = CUBLAS_SIDE_RIGHT;
#endif
#if defined(ROCBLAS_VERSION_MAJOR)
    static constexpr rocblas_side rocblas_param = rocblas_side_right;
#endif
  };
};
template <typename T> struct is_valid_side_tag {
  enum : bool { value = (std::is_same<T, Side::Left>::value || std::is_same<T, Side::Right>::value) };
};
template <typename T> struct flip_side_tag;
template <> struct flip_side_tag<Side::Left> { typedef Side::Right type; };
template <> struct flip_side_tag<Side::Right> { typedef Side::Left type; };

struct Diag {
  enum : int { tag = 600 };
  struct Unit {
    enum : int { tag = 601 };
    static constexpr char param = 'U';
#if defined(__INTEL_MKL__)
    static constexpr CBLAS_DIAG mkl_param = CblasUnit;
#endif
#if defined(CUBLAS_VERSION)
    static constexpr cublasDiagType_t cublas_param = CUBLAS_DIAG_UNIT;
#endif
#if defined(ROCBLAS_VERSION_MAJOR)
    static constexpr rocblas_diagonal rocblas_param = rocblas_diagonal_unit;
#endif
  };
  struct NonUnit {
    enum : int { tag = 602 };
    static constexpr char param = 'N';
#if defined(__INTEL_MKL__)
    static constexpr CBLAS_DIAG mkl_param = CblasNonUnit;
#endif
#if defined(CUBLAS_VERSION)
    static constexpr cublasDiagType_t cublas_param = CUBLAS_DIAG_NON_UNIT;
#endif
#if defined(ROCBLAS_VERSION_MAJOR)
    static constexpr rocblas_diagonal rocblas_param = rocblas_diagonal_non_unit;
#endif
  };
};
template <typename T> struct is_valid_diag_tag {
  enum : bool { value = (std::is_same<T, Diag::Unit>::value || std::is_same<T, Diag::NonUnit>::value) };
};

struct Trans {
  enum : int { tag = 700 };
  struct Transpose {
    enum : int { tag = 701 };
    static constexpr char param = 'T';
#if defined(__INTEL_MKL__)
    static constexpr CBLAS_TRANSPOSE mkl_param = CblasTrans;
#endif
#if defined(CUBLAS_VERSION)
    static constexpr cublasOperation_t cublas_param = CUBLAS_OP_T;
#endif
#if defined(ROCBLAS_VERSION_MAJOR)
    static constexpr rocblas_operation rocblas_param = rocblas_operation_transpose;
#endif
  };
  struct ConjTranspose {
    enum : int { tag = 702 };
    static constexpr char param = 'C';
#if defined(__INTEL_MKL__)
    static constexpr CBLAS_TRANSPOSE mkl_param = CblasConjTrans;
#endif
#if defined(CUBLAS_VERSION)
    static constexpr cublasOperation_t cublas_param = CUBLAS_OP_C;
#endif
#if defined(ROCBLAS_VERSION_MAJOR)
    static constexpr rocblas_operation rocblas_param = rocblas_operation_conjugate_transpose;
#endif
  };
  struct NoTranspose {
    enum : int { tag = 703 };
    static constexpr char param = 'N';
#if defined(__INTEL_MKL__)
    static constexpr CBLAS_TRANSPOSE mkl_param = CblasNoTrans;
#endif
#if defined(CUBLAS_VERSION)
    static constexpr cublasOperation_t cublas_param = CUBLAS_OP_N;
#endif
#if defined(ROCBLAS_VERSION_MAJOR)
    static constexpr rocblas_operation rocblas_param = rocblas_operation_none;
#endif
  };
};
template <typename T> struct is_valid_trans_tag {
  enum : bool {
    value = (std::is_same<T, Trans::Transpose>::value || std::is_same<T, Trans::ConjTranspose>::value ||
             std::is_same<T, Trans::NoTranspose>::value)
  };
};
template <typename T> struct transpose_trans_tag;
template <typename T> struct conj_transpose_trans_tag;

template <> struct transpose_trans_tag<Trans::Transpose> { typedef Trans::NoTranspose type; };
template <> struct transpose_trans_tag<Trans::ConjTranspose> { typedef Trans::NoTranspose type; };
template <> struct transpose_trans_tag<Trans::NoTranspose> { typedef Trans::Transpose type; };

template <> struct conj_transpose_trans_tag<Trans::Transpose> { typedef Trans::NoTranspose type; };
template <> struct conj_transpose_trans_tag<Trans::ConjTranspose> { typedef Trans::NoTranspose type; };
template <> struct conj_transpose_trans_tag<Trans::NoTranspose> { typedef Trans::ConjTranspose type; };

struct Direct {
  enum : int { tag = 800 };
  struct Forward {
    enum : int { tag = 801 };
  };
  struct Backward {
    enum : int { tag = 802 };
  };
};

///
/// helper functions
///
struct Conjugate {
  enum : int { tag = 801 };

  KOKKOS_FORCEINLINE_FUNCTION Conjugate() {}
  KOKKOS_FORCEINLINE_FUNCTION Conjugate(const Conjugate &b) {}

  KOKKOS_FORCEINLINE_FUNCTION float operator()(const float &v) const { return v; }
  KOKKOS_FORCEINLINE_FUNCTION double operator()(const double &v) const { return v; }
  inline std::complex<float> operator()(const std::complex<float> &v) const { return std::conj(v); }
  inline std::complex<double> operator()(const std::complex<double> &v) const { return std::conj(v); }
  KOKKOS_FORCEINLINE_FUNCTION Kokkos::complex<float> operator()(const Kokkos::complex<float> &v) const {
    return Kokkos::conj(v);
  }
  KOKKOS_FORCEINLINE_FUNCTION Kokkos::complex<double> operator()(const Kokkos::complex<double> &v) const {
    return Kokkos::conj(v);
  }
};

struct NoConjugate {
  enum : int { tag = 802 };

  KOKKOS_FORCEINLINE_FUNCTION NoConjugate() {}
  KOKKOS_FORCEINLINE_FUNCTION NoConjugate(const NoConjugate &b) {}

  KOKKOS_FORCEINLINE_FUNCTION float operator()(const float &v) const { return v; }
  KOKKOS_FORCEINLINE_FUNCTION double operator()(const double &v) const { return v; }
  inline std::complex<float> operator()(const std::complex<float> &v) const { return v; }
  inline std::complex<double> operator()(const std::complex<double> &v) const { return v; }
  KOKKOS_FORCEINLINE_FUNCTION Kokkos::complex<float> operator()(const Kokkos::complex<float> &v) const { return v; }
  KOKKOS_FORCEINLINE_FUNCTION Kokkos::complex<double> operator()(const Kokkos::complex<double> &v) const { return v; }
};

struct Algo {
  struct External {
    enum : int { tag = 1001 };
  };
  struct Internal {
    enum : int { tag = 1002 };
  };
  struct ByBlocks {
    enum : int { tag = 1003 };
  };
  struct OnDevice {
    enum : int { tag = 1004 };
  };
  struct Serial {
    enum : int { tag = 1005 };
  };

  struct Workflow {
    struct Serial {
      enum : int { tag = 2001 };
    };
    struct SerialPanel {
      enum : int { tag = 2002 };
    };
  };
};


template <bool isCudaOrHIP>
struct ActiveAlgorithm {
  using type = Algo::Internal;
};
template <>
struct ActiveAlgorithm<false> {
  using type = Algo::External;
};

template <bool withOMP>
struct ActiveHostAlgorithm {
  using type = Algo::Serial;
};
template <>
struct ActiveHostAlgorithm<false> {
  using type = Algo::External;
};

template <typename MemoryTraitsType, Kokkos::MemoryTraitsFlags flag>
using MemoryTraits = Kokkos::MemoryTraits<MemoryTraitsType::is_unmanaged | MemoryTraitsType::is_random_access |
                                          MemoryTraitsType::is_atomic | flag>;

template <typename ViewType>
using UnmanagedViewType =
    Kokkos::View<typename ViewType::data_type, typename ViewType::array_layout, typename ViewType::device_type,
                 MemoryTraits<typename ViewType::memory_traits, Kokkos::Unmanaged>>;
template <typename ViewType>
using ConstViewType = Kokkos::View<typename ViewType::const_data_type, typename ViewType::array_layout,
                                   typename ViewType::device_type, typename ViewType::memory_traits>;
template <typename ViewType> using ConstUnmanagedViewType = ConstViewType<UnmanagedViewType<ViewType>>;

using do_not_initialize_tag = Kokkos::ViewAllocateWithoutInitializing;

template <typename T> struct ExecSpaceFactory {
  static void createInstance(T &exec_instance) { exec_instance = T(); }
#if defined(KOKKOS_ENABLE_CUDA)
  static void createInstance(const cudaStream_t &s, T &exec_instance) { exec_instance = T(); }
#endif
#if defined(KOKKOS_ENABLE_HIP)
  static void createInstance(const hipStream_t &s, T &exec_instance) { exec_instance = T(); }
#endif
};

#if defined(KOKKOS_ENABLE_CUDA)
template <> struct ExecSpaceFactory<Kokkos::Cuda> {
  static void createInstance(Kokkos::Cuda &exec_instance) { exec_instance = Kokkos::Cuda(); }
  static void createInstance(const cudaStream_t &s, Kokkos::Cuda &exec_instance) { exec_instance = Kokkos::Cuda(s); }
};
#endif
#if defined(KOKKOS_ENABLE_HIP)
template <> struct ExecSpaceFactory<Kokkos::HIP> {
  static void createInstance(Kokkos::HIP &exec_instance) { exec_instance = Kokkos::HIP(); }
  static void createInstance(const hipStream_t &s, Kokkos::HIP &exec_instance) {
    exec_instance = Kokkos::HIP(s);
  }
};
#endif

} // namespace Tacho

#endif
