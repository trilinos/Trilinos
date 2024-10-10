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

#ifndef KOKKOSSPARSE_SPMV_HANDLE_HPP_
#define KOKKOSSPARSE_SPMV_HANDLE_HPP_

#include <Kokkos_Core.hpp>
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_BsrMatrix.hpp"
// Use TPL utilities for safely finalizing matrix descriptors, etc.
#include "KokkosSparse_Utils_cusparse.hpp"
#include "KokkosSparse_Utils_rocsparse.hpp"
#include "KokkosSparse_Utils_mkl.hpp"

namespace KokkosSparse {

/// SPMVAlgorithm values can be used to select different algorithms/methods for
/// performing SpMV computations.
enum SPMVAlgorithm {
  SPMV_DEFAULT,            /// Default algorithm: best overall performance for repeated
                           /// applications of SpMV.
  SPMV_FAST_SETUP,         /// Best performance in the non-reuse case, where the handle
                           /// is only used once.
  SPMV_NATIVE,             /// Use the best KokkosKernels implementation, even if a TPL
                           /// implementation is available.
  SPMV_MERGE_PATH,         /// Use algorithm optimized for matrices with
                           /// imbalanced/irregular sparsity patterns (merge path or
                           /// similar). May call a TPL. For CrsMatrix only.
  SPMV_NATIVE_MERGE_PATH,  /// Use the KokkosKernels implementation of merge
                           /// path. For CrsMatrix only.
  SPMV_BSR_V41,            /// Use experimental version 4.1 algorithm (for BsrMatrix only)
  SPMV_BSR_V42,            /// Use experimental version 4.2 algorithm (for BsrMatrix only)
  SPMV_BSR_TC              /// Use experimental tensor core algorithm (for BsrMatrix only)
};

namespace Experimental {
/// Precision to use in the tensor core implementation of Bsr SpMV
enum class Bsr_TC_Precision {
  Automatic,  ///< Use Double, unless operations match mixed precision
  Double,     ///< fp64 += fp64 * fp64
  Mixed       ///< fp32 += fp16 * fp16
};
}  // namespace Experimental

/// Get the name of a SPMVAlgorithm enum constant
inline const char* get_spmv_algorithm_name(SPMVAlgorithm a) {
  switch (a) {
    case SPMV_DEFAULT: return "SPMV_DEFAULT";
    case SPMV_FAST_SETUP: return "SPMV_FAST_SETUP";
    case SPMV_NATIVE: return "SPMV_NATIVE";
    case SPMV_MERGE_PATH: return "SPMV_MERGE_PATH";
    case SPMV_NATIVE_MERGE_PATH: return "SPMV_NATIVE_MERGE_PATH";
    case SPMV_BSR_V41: return "SPMV_BSR_V41";
    case SPMV_BSR_V42: return "SPMV_BSR_V42";
    case SPMV_BSR_TC: return "SPMV_BSR_TC";
  }
  throw std::invalid_argument("SPMVHandle::get_algorithm_name: unknown algorithm");
  return "<Unknown>";
}

/// Return true if the given algorithm is always a native (KokkosKernels)
/// implementation, and false if it may be implemented by a TPL.
inline bool is_spmv_algorithm_native(SPMVAlgorithm a) {
  switch (a) {
    case SPMV_NATIVE:
    case SPMV_NATIVE_MERGE_PATH:
    case SPMV_BSR_V41:
    case SPMV_BSR_V42:
    case SPMV_BSR_TC: return true;
    // DEFAULT, FAST_SETUP and MERGE_PATH may call TPLs
    default: return false;
  }
}

namespace Impl {

template <typename ExecutionSpace>
struct TPL_SpMV_Data {
  // Disallow default construction: must provide the initial execution space
  TPL_SpMV_Data() = delete;
  TPL_SpMV_Data(const ExecutionSpace& exec_) : exec(exec_) {}
  void set_exec_space(const ExecutionSpace& new_exec) {
    // Check if new_exec is different from (old) exec.
    // If it is, fence the old exec now.
    // That way, SPMVHandle cleanup doesn't need
    // to worry about resources still being in use on the old exec.
    if (exec != new_exec) {
      exec.fence();
      exec = new_exec;
    }
  }
  virtual ~TPL_SpMV_Data() {}
  ExecutionSpace exec;
};

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#if defined(CUSPARSE_VERSION) && (10300 <= CUSPARSE_VERSION)
// Data used by cuSPARSE >=10.3 for both single-vector (SpMV) and multi-vector
// (SpMM).
// TODO: in future, this can also be used for BSR (cuSPARSE >=12.2)
struct CuSparse10_SpMV_Data : public TPL_SpMV_Data<Kokkos::Cuda> {
  CuSparse10_SpMV_Data(const Kokkos::Cuda& exec_) : TPL_SpMV_Data(exec_) {}
  ~CuSparse10_SpMV_Data() {
    // Prefer cudaFreeAsync on the stream that last executed a spmv, but
    // async memory management was introduced in 11.2
#if (CUDA_VERSION >= 11020)
    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaFreeAsync(buffer, exec.cuda_stream()));
#else
    // Fence here to ensure spmv is not still using buffer
    // (cudaFree does not do a device synchronize)
    exec.fence();
    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaFree(buffer));
#endif
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroySpMat(mat));
  }

  cusparseSpMatDescr_t mat;
  size_t bufferSize = 0;
  void* buffer      = nullptr;
};
#endif

// Data used by cuSPARSE <10.3 for CRS, and >=9 for BSR
struct CuSparse9_SpMV_Data : public TPL_SpMV_Data<Kokkos::Cuda> {
  CuSparse9_SpMV_Data(const Kokkos::Cuda& exec_) : TPL_SpMV_Data(exec_) {}
  ~CuSparse9_SpMV_Data() { KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroyMatDescr(mat)); }

  cusparseMatDescr_t mat;
};
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
struct RocSparse_CRS_SpMV_Data : public TPL_SpMV_Data<Kokkos::HIP> {
  RocSparse_CRS_SpMV_Data(const Kokkos::HIP& exec_) : TPL_SpMV_Data(exec_) {}
  ~RocSparse_CRS_SpMV_Data() {
    // note: hipFree includes an implicit device synchronize
    KOKKOS_IMPL_HIP_SAFE_CALL(hipFree(buffer));
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_destroy_spmat_descr(mat));
  }

  rocsparse_spmat_descr mat;
  size_t bufferSize = 0;
  void* buffer      = nullptr;
};

struct RocSparse_BSR_SpMV_Data : public TPL_SpMV_Data<Kokkos::HIP> {
  RocSparse_BSR_SpMV_Data(const Kokkos::HIP& exec_) : TPL_SpMV_Data(exec_) {}
  ~RocSparse_BSR_SpMV_Data() {
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_destroy_mat_descr(mat));
#if (KOKKOSSPARSE_IMPL_ROCM_VERSION >= 50400)
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_destroy_mat_info(info));
#endif
  }

  rocsparse_mat_descr mat;
#if (KOKKOSSPARSE_IMPL_ROCM_VERSION >= 50400)
  rocsparse_mat_info info;
#endif
};
#endif

// note: header defining __INTEL_MKL__ is pulled in above by Utils_mkl.hpp
#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL

#if (__INTEL_MKL__ > 2017)
// Data for classic MKL (both CRS and BSR)
template <typename ExecutionSpace>
struct MKL_SpMV_Data : public TPL_SpMV_Data<ExecutionSpace> {
  MKL_SpMV_Data(const ExecutionSpace& exec_) : TPL_SpMV_Data<ExecutionSpace>(exec_) {}
  ~MKL_SpMV_Data() {
    KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_destroy(mat));
    // descr is just a plain-old-data struct, no cleanup to do
  }

  sparse_matrix_t mat;
  matrix_descr descr;
};
#endif

#if defined(KOKKOS_ENABLE_SYCL)
struct OneMKL_SpMV_Data : public TPL_SpMV_Data<Kokkos::Experimental::SYCL> {
  OneMKL_SpMV_Data(const Kokkos::Experimental::SYCL& exec_) : TPL_SpMV_Data(exec_) {}
  ~OneMKL_SpMV_Data() {
    // Make sure no spmv is still running with this handle, if exec uses an
    // out-of-order queue (rare case)
    if (!exec.sycl_queue().is_in_order()) exec.fence();
#if INTEL_MKL_VERSION >= 20230200
    // MKL 2023.2 and up make this async release okay even though it takes a
    // pointer to mat, which is going out of scope after this destructor
    oneapi::mkl::sparse::release_matrix_handle(exec.sycl_queue(), &mat);
#else
    // But in older versions, wait on ev_release before letting mat go out of
    // scope
    auto ev_release = oneapi::mkl::sparse::release_matrix_handle(exec.sycl_queue(), &mat);
    ev_release.wait();
#endif
  }

  oneapi::mkl::sparse::matrix_handle_t mat;
};
#endif
#endif

template <class ExecutionSpace, class MemorySpace, class Scalar, class Offset, class Ordinal>
struct SPMVHandleImpl {
  using ExecutionSpaceType = ExecutionSpace;
  // This is its own ImplType
  using ImplType = SPMVHandleImpl<ExecutionSpace, MemorySpace, Scalar, Offset, Ordinal>;
  // Do not allow const qualifier on Scalar, Ordinal, Offset (otherwise this
  // type won't match the ETI'd type). Users should not use SPMVHandleImpl
  // directly and SPMVHandle explicitly removes const, so this should never
  // happen in practice.
  static_assert(!std::is_const_v<Scalar>, "SPMVHandleImpl: Scalar must not be a const type");
  static_assert(!std::is_const_v<Offset>, "SPMVHandleImpl: Offset must not be a const type");
  static_assert(!std::is_const_v<Ordinal>, "SPMVHandleImpl: Ordinal must not be a const type");
  SPMVHandleImpl(SPMVAlgorithm algo_) : algo(algo_) {}
  ~SPMVHandleImpl() {
    if (tpl_rank1) delete tpl_rank1;
    if (tpl_rank2) delete tpl_rank2;
  }

  ImplType* get_impl() { return this; }

  /// Get the SPMVAlgorithm used by this handle
  SPMVAlgorithm get_algorithm() const { return this->algo; }

  const SPMVAlgorithm algo                 = SPMV_DEFAULT;
  TPL_SpMV_Data<ExecutionSpace>* tpl_rank1 = nullptr;
  TPL_SpMV_Data<ExecutionSpace>* tpl_rank2 = nullptr;
  // Expert tuning parameters for native SpMV
  // TODO: expose a proper Experimental interface to set these. Currently they
  // can be assigned directly in the SPMVHandle as they are public members.
  int team_size               = -1;
  int vector_length           = -1;
  int64_t rows_per_thread     = -1;
  bool force_static_schedule  = false;
  bool force_dynamic_schedule = false;
  KokkosSparse::Experimental::Bsr_TC_Precision bsr_tc_precision =
      KokkosSparse::Experimental::Bsr_TC_Precision::Automatic;
};
}  // namespace Impl

// clang-format off
/// \class SPMVHandle
/// \brief Opaque handle type for KokkosSparse::spmv. It passes the choice of
///    algorithm to the spmv implementation, and also may store internal data which can be used to
///    speed up the spmv computation.
/// \tparam DeviceType A Kokkos::Device or execution space where the spmv computation will be run.
///    Does not necessarily need to match AMatrix's device type, but its execution space needs to be able
///    to access the memory spaces of AMatrix, XVector and YVector.
/// \tparam AMatrix A specialization of KokkosSparse::CrsMatrix or
/// KokkosSparse::BsrMatrix.
///
/// SPMVHandle's internal resources are lazily allocated and initialized by the first
/// spmv call.
///
/// SPMVHandle automatically cleans up all allocated resources when it is destructed.
/// No fencing by the user is required between the final spmv and cleanup.
///
/// A SPMVHandle instance can be used in any number of calls, with any execution space
/// instance and any X/Y vectors (with matching types) each call.
///
/// \warning However, all calls to spmv with a given instance of SPMVHandle must use the
/// same matrix.
// clang-format on

template <class DeviceType, class AMatrix, class XVector, class YVector>
struct SPMVHandle
    : public Impl::SPMVHandleImpl<typename DeviceType::execution_space, typename AMatrix::memory_space,
                                  typename AMatrix::non_const_value_type, typename AMatrix::non_const_size_type,
                                  typename AMatrix::non_const_ordinal_type> {
  using ImplType = Impl::SPMVHandleImpl<typename DeviceType::execution_space, typename AMatrix::memory_space,
                                        typename AMatrix::non_const_value_type, typename AMatrix::non_const_size_type,
                                        typename AMatrix::non_const_ordinal_type>;
  // Note: these typedef names cannot shadow template parameters
  using AMatrixType        = AMatrix;
  using XVectorType        = XVector;
  using YVectorType        = YVector;
  using ExecutionSpaceType = typename DeviceType::execution_space;
  // Check all template parameters for compatibility with each other
  // NOTE: we do not require that ExecutionSpace matches
  // AMatrix::execution_space. For example, if the matrix's device is <Cuda,
  // CudaHostPinnedSpace> it is allowed to run spmv on Serial.
  static_assert(is_crs_matrix_v<AMatrix> || Experimental::is_bsr_matrix_v<AMatrix>,
                "SPMVHandle: AMatrix must be a specialization of CrsMatrix or "
                "BsrMatrix.");
  static_assert(Kokkos::is_view<XVector>::value, "SPMVHandle: XVector must be a Kokkos::View.");
  static_assert(Kokkos::is_view<YVector>::value, "SPMVHandle: YVector must be a Kokkos::View.");
  static_assert(XVector::rank() == YVector::rank(), "SPMVHandle: ranks of XVector and YVector must match.");
  static_assert(XVector::rank() == size_t(1) || YVector::rank() == size_t(2),
                "SPMVHandle: XVector and YVector must be both rank-1 or both rank-2.");
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpaceType, typename AMatrix::memory_space>::accessible,
                "SPMVHandle: AMatrix must be accessible from ExecutionSpace");
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpaceType, typename XVector::memory_space>::accessible,
                "SPMVHandle: XVector must be accessible from ExecutionSpace");
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpaceType, typename YVector::memory_space>::accessible,
                "SPMVHandle: YVector must be accessible from ExecutionSpace");

  // Prevent copying (this object does not support reference counting)
  SPMVHandle(const SPMVHandle&)            = delete;
  SPMVHandle& operator=(const SPMVHandle&) = delete;

  /// \brief Create a new SPMVHandle using the given algorithm.
  SPMVHandle(SPMVAlgorithm algo_ = SPMV_DEFAULT) : ImplType(algo_) {
    // Validate the choice of algorithm based on A's type
    if constexpr (is_crs_matrix_v<AMatrixType>) {
      switch (get_algorithm()) {
        case SPMV_BSR_V41:
        case SPMV_BSR_V42:
        case SPMV_BSR_TC:
          throw std::invalid_argument(std::string("SPMVHandle: algorithm ") + get_spmv_algorithm_name(get_algorithm()) +
                                      " cannot be used if A is a CrsMatrix");
        default:;
      }
    } else {
      switch (get_algorithm()) {
        case SPMV_MERGE_PATH:
        case SPMV_NATIVE_MERGE_PATH:
          throw std::invalid_argument(std::string("SPMVHandle: algorithm ") + get_spmv_algorithm_name(get_algorithm()) +
                                      " cannot be used if A is a BsrMatrix");
        default:;
      }
    }
  }

  /// Get the SPMVAlgorithm used by this handle
  SPMVAlgorithm get_algorithm() const {
    // Note: get_algorithm is also a method of parent ImplType, but for
    // documentation purposes it should appear directly in the public interface
    // of SPMVHandle
    return this->algo;
  }

  /// Get pointer to this as the impl type
  ImplType* get_impl() { return static_cast<ImplType*>(this); }
};

namespace Impl {
template <typename>
struct is_spmv_handle : public std::false_type {};
template <typename... P>
struct is_spmv_handle<SPMVHandle<P...>> : public std::true_type {};
template <typename... P>
struct is_spmv_handle<const SPMVHandle<P...>> : public std::true_type {};

template <typename T>
inline constexpr bool is_spmv_handle_v = is_spmv_handle<T>::value;
}  // namespace Impl

}  // namespace KokkosSparse

#endif
