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

#ifndef KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_DECL_HPP
#define KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_DECL_HPP

#include "KokkosKernels_AlwaysFalse.hpp"
#include "KokkosSparse_Utils_mkl.hpp"
#include "KokkosSparse_Utils_cusparse.hpp"
#include "KokkosKernels_tpl_handles_decl.hpp"

#if defined(KOKKOSKERNELS_ENABLE_TPL_MKL) && (__INTEL_MKL__ > 2017)
#include <mkl.h>

namespace KokkosSparse {
namespace Impl {

// MKL 2018 and above: use new interface: sparse_matrix_t and mkl_sparse_?_mv()

// Note: Scalar here is the Kokkos type, not the MKL type
template <typename Scalar, typename Handle>
inline void spmv_bsr_mkl(Handle* handle, sparse_operation_t op, Scalar alpha, Scalar beta, MKL_INT m, MKL_INT n,
                         MKL_INT b, const MKL_INT* Arowptrs, const MKL_INT* Aentries, const Scalar* Avalues,
                         const Scalar* x, Scalar* y) {
  using MKLScalar = typename KokkosSparse::Impl::KokkosToMKLScalar<Scalar>::type;
  using ExecSpace = typename Handle::ExecutionSpaceType;
  using Subhandle = KokkosSparse::Impl::MKL_SpMV_Data<ExecSpace>;
  Subhandle* subhandle;
  const MKLScalar* x_mkl = reinterpret_cast<const MKLScalar*>(x);
  MKLScalar* y_mkl       = reinterpret_cast<MKLScalar*>(y);
  if (handle->tpl_rank1) {
    subhandle = dynamic_cast<Subhandle*>(handle->tpl_rank1);
    if (!subhandle) throw std::runtime_error("KokkosSparse::spmv: subhandle is not set up for MKL BSR");
    // note: classic mkl only runs on synchronous host exec spaces, so no need
    // to call set_exec_space on the subhandle here
  } else {
    // Use the default execution space instance, as classic MKL does not use
    // a specific instance.
    subhandle             = new Subhandle(ExecSpace());
    handle->tpl_rank1     = subhandle;
    subhandle->descr.type = SPARSE_MATRIX_TYPE_GENERAL;
    subhandle->descr.mode = SPARSE_FILL_MODE_FULL;
    subhandle->descr.diag = SPARSE_DIAG_NON_UNIT;
    // Note: the create_csr routine requires non-const values even though
    // they're not actually modified
    MKLScalar* Avalues_mkl = reinterpret_cast<MKLScalar*>(const_cast<Scalar*>(Avalues));
    if constexpr (std::is_same_v<Scalar, float>) {
      KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_s_create_bsr(
          &subhandle->mat, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, m, n, b, const_cast<MKL_INT*>(Arowptrs),
          const_cast<MKL_INT*>(Arowptrs + 1), const_cast<MKL_INT*>(Aentries), Avalues_mkl));
    } else if constexpr (std::is_same_v<Scalar, double>) {
      KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_d_create_bsr(
          &subhandle->mat, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, m, n, b, const_cast<MKL_INT*>(Arowptrs),
          const_cast<MKL_INT*>(Arowptrs + 1), const_cast<MKL_INT*>(Aentries), Avalues_mkl));
    } else if constexpr (std::is_same_v<Scalar, Kokkos::complex<float>>) {
      KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_c_create_bsr(
          &subhandle->mat, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, m, n, b, const_cast<MKL_INT*>(Arowptrs),
          const_cast<MKL_INT*>(Arowptrs + 1), const_cast<MKL_INT*>(Aentries), Avalues_mkl));
    } else if constexpr (std::is_same_v<Scalar, Kokkos::complex<double>>) {
      KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_z_create_bsr(
          &subhandle->mat, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, m, n, b, const_cast<MKL_INT*>(Arowptrs),
          const_cast<MKL_INT*>(Arowptrs + 1), const_cast<MKL_INT*>(Aentries), Avalues_mkl));
    }
  }
  MKLScalar alpha_mkl = KokkosSparse::Impl::KokkosToMKLScalar<Scalar>(alpha);
  MKLScalar beta_mkl  = KokkosSparse::Impl::KokkosToMKLScalar<Scalar>(beta);
  if constexpr (std::is_same_v<Scalar, float>) {
    KOKKOSKERNELS_MKL_SAFE_CALL(
        mkl_sparse_s_mv(op, alpha_mkl, subhandle->mat, subhandle->descr, x_mkl, beta_mkl, y_mkl));
  } else if constexpr (std::is_same_v<Scalar, double>) {
    KOKKOSKERNELS_MKL_SAFE_CALL(
        mkl_sparse_d_mv(op, alpha_mkl, subhandle->mat, subhandle->descr, x_mkl, beta_mkl, y_mkl));
  } else if constexpr (std::is_same_v<Scalar, Kokkos::complex<float>>) {
    KOKKOSKERNELS_MKL_SAFE_CALL(
        mkl_sparse_c_mv(op, alpha_mkl, subhandle->mat, subhandle->descr, x_mkl, beta_mkl, y_mkl));
  } else if constexpr (std::is_same_v<Scalar, Kokkos::complex<double>>) {
    KOKKOSKERNELS_MKL_SAFE_CALL(
        mkl_sparse_z_mv(op, alpha_mkl, subhandle->mat, subhandle->descr, x_mkl, beta_mkl, y_mkl));
  }
}

// Note: Scalar here is the Kokkos type, not the MKL type
template <typename Scalar, typename Handle>
inline void spmv_mv_bsr_mkl(Handle* handle, sparse_operation_t op, Scalar alpha, Scalar beta, MKL_INT m, MKL_INT n,
                            MKL_INT b, const MKL_INT* Arowptrs, const MKL_INT* Aentries, const Scalar* Avalues,
                            const Scalar* x, MKL_INT colx, MKL_INT ldx, Scalar* y, MKL_INT ldy) {
  using MKLScalar = typename KokkosSparse::Impl::KokkosToMKLScalar<Scalar>::type;
  using ExecSpace = typename Handle::ExecutionSpaceType;
  using Subhandle = KokkosSparse::Impl::MKL_SpMV_Data<ExecSpace>;
  Subhandle* subhandle;
  const MKLScalar* x_mkl = reinterpret_cast<const MKLScalar*>(x);
  MKLScalar* y_mkl       = reinterpret_cast<MKLScalar*>(y);
  if (handle->tpl_rank2) {
    subhandle = dynamic_cast<Subhandle*>(handle->tpl_rank2);
    if (!subhandle) throw std::runtime_error("KokkosSparse::spmv: subhandle is not set up for MKL BSR");
    // note: classic mkl only runs on synchronous host exec spaces, so no need
    // to call set_exec_space on the subhandle here
  } else {
    // Use the default execution space instance, as classic MKL does not use
    // a specific instance.
    subhandle             = new Subhandle(ExecSpace());
    handle->tpl_rank2     = subhandle;
    subhandle->descr.type = SPARSE_MATRIX_TYPE_GENERAL;
    subhandle->descr.mode = SPARSE_FILL_MODE_FULL;
    subhandle->descr.diag = SPARSE_DIAG_NON_UNIT;
    // Note: the create_csr routine requires non-const values even though
    // they're not actually modified
    MKLScalar* Avalues_mkl = reinterpret_cast<MKLScalar*>(const_cast<Scalar*>(Avalues));
    if constexpr (std::is_same_v<Scalar, float>) {
      KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_s_create_bsr(
          &subhandle->mat, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, m, n, b, const_cast<MKL_INT*>(Arowptrs),
          const_cast<MKL_INT*>(Arowptrs + 1), const_cast<MKL_INT*>(Aentries), Avalues_mkl));
    } else if constexpr (std::is_same_v<Scalar, double>) {
      KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_d_create_bsr(
          &subhandle->mat, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, m, n, b, const_cast<MKL_INT*>(Arowptrs),
          const_cast<MKL_INT*>(Arowptrs + 1), const_cast<MKL_INT*>(Aentries), Avalues_mkl));
    } else if constexpr (std::is_same_v<Scalar, Kokkos::complex<float>>) {
      KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_c_create_bsr(
          &subhandle->mat, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, m, n, b, const_cast<MKL_INT*>(Arowptrs),
          const_cast<MKL_INT*>(Arowptrs + 1), const_cast<MKL_INT*>(Aentries), Avalues_mkl));
    } else if constexpr (std::is_same_v<Scalar, Kokkos::complex<double>>) {
      KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_z_create_bsr(
          &subhandle->mat, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, m, n, b, const_cast<MKL_INT*>(Arowptrs),
          const_cast<MKL_INT*>(Arowptrs + 1), const_cast<MKL_INT*>(Aentries), Avalues_mkl));
    }
  }
  MKLScalar alpha_mkl = KokkosSparse::Impl::KokkosToMKLScalar<Scalar>(alpha);
  MKLScalar beta_mkl  = KokkosSparse::Impl::KokkosToMKLScalar<Scalar>(beta);
  if constexpr (std::is_same_v<Scalar, float>) {
    KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_s_mm(op, alpha_mkl, subhandle->mat, subhandle->descr,
                                                SPARSE_LAYOUT_ROW_MAJOR, x_mkl, colx, ldx, beta_mkl, y_mkl, ldy));
  } else if constexpr (std::is_same_v<Scalar, double>) {
    KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_d_mm(op, alpha_mkl, subhandle->mat, subhandle->descr,
                                                SPARSE_LAYOUT_ROW_MAJOR, x_mkl, colx, ldx, beta_mkl, y_mkl, ldy));
  } else if constexpr (std::is_same_v<Scalar, Kokkos::complex<float>>) {
    KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_c_mm(op, alpha_mkl, subhandle->mat, subhandle->descr,
                                                SPARSE_LAYOUT_ROW_MAJOR, x_mkl, colx, ldx, beta_mkl, y_mkl, ldy));
  } else if constexpr (std::is_same_v<Scalar, Kokkos::complex<double>>) {
    KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_z_mm(op, alpha_mkl, subhandle->mat, subhandle->descr,
                                                SPARSE_LAYOUT_ROW_MAJOR, x_mkl, colx, ldx, beta_mkl, y_mkl, ldy));
  }
}

#define KOKKOSSPARSE_SPMV_MKL(SCALAR, EXECSPACE)                                                                         \
  template <>                                                                                                            \
  struct SPMV_BSRMATRIX<EXECSPACE,                                                                                       \
                        KokkosSparse::Impl::SPMVHandleImpl<EXECSPACE, Kokkos::HostSpace, SCALAR, MKL_INT, MKL_INT>,      \
                        ::KokkosSparse::Experimental::BsrMatrix<                                                         \
                            SCALAR const, MKL_INT const, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,                   \
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>, MKL_INT const>,                                     \
                        Kokkos::View<SCALAR const*, Kokkos::LayoutLeft, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,    \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                    \
                        Kokkos::View<SCALAR*, Kokkos::LayoutLeft, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,          \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                           \
                        true> {                                                                                          \
    using device_type = Kokkos::Device<EXECSPACE, Kokkos::HostSpace>;                                                    \
    using Handle      = KokkosSparse::Impl::SPMVHandleImpl<EXECSPACE, Kokkos::HostSpace, SCALAR, MKL_INT, MKL_INT>;      \
    using AMatrix     = ::KokkosSparse::Experimental::BsrMatrix<SCALAR const, MKL_INT const, device_type,                \
                                                            Kokkos::MemoryTraits<Kokkos::Unmanaged>, MKL_INT const>; \
    using XVector     = Kokkos::View<SCALAR const*, Kokkos::LayoutLeft, device_type,                                     \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;                    \
    using YVector     = Kokkos::View<SCALAR*, Kokkos::LayoutLeft, device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    using coefficient_type = typename YVector::non_const_value_type;                                                     \
                                                                                                                         \
    static void spmv_bsrmatrix(const EXECSPACE&, Handle* handle, const char mode[], const coefficient_type& alpha,       \
                               const AMatrix& A, const XVector& X, const coefficient_type& beta, const YVector& Y) {     \
      std::string label = "KokkosSparse::spmv[TPL_MKL,BSRMATRIX," + Kokkos::ArithTraits<SCALAR>::name() + "]";           \
      Kokkos::Profiling::pushRegion(label);                                                                              \
      spmv_bsr_mkl(handle, mode_kk_to_mkl(mode[0]), alpha, beta, A.numRows(), A.numCols(), A.blockDim(),                 \
                   A.graph.row_map.data(), A.graph.entries.data(), A.values.data(), X.data(), Y.data());                 \
      Kokkos::Profiling::popRegion();                                                                                    \
    }                                                                                                                    \
  };

#ifdef KOKKOS_ENABLE_SERIAL
KOKKOSSPARSE_SPMV_MKL(float, Kokkos::Serial)
KOKKOSSPARSE_SPMV_MKL(double, Kokkos::Serial)
KOKKOSSPARSE_SPMV_MKL(Kokkos::complex<float>, Kokkos::Serial)
KOKKOSSPARSE_SPMV_MKL(Kokkos::complex<double>, Kokkos::Serial)
#endif

#ifdef KOKKOS_ENABLE_OPENMP
KOKKOSSPARSE_SPMV_MKL(float, Kokkos::OpenMP)
KOKKOSSPARSE_SPMV_MKL(double, Kokkos::OpenMP)
KOKKOSSPARSE_SPMV_MKL(Kokkos::complex<float>, Kokkos::OpenMP)
KOKKOSSPARSE_SPMV_MKL(Kokkos::complex<double>, Kokkos::OpenMP)
#endif

#undef KOKKOSSPARSE_SPMV_MKL

#define KOKKOSSPARSE_SPMV_MV_MKL(SCALAR, EXECSPACE)                                                                      \
  template <>                                                                                                            \
  struct SPMV_MV_BSRMATRIX<                                                                                              \
      EXECSPACE, KokkosSparse::Impl::SPMVHandleImpl<EXECSPACE, Kokkos::HostSpace, SCALAR, MKL_INT, MKL_INT>,             \
      ::KokkosSparse::Experimental::BsrMatrix<SCALAR const, MKL_INT const,                                               \
                                              Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,                              \
                                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, MKL_INT const>,                   \
      Kokkos::View<SCALAR const**, Kokkos::LayoutLeft, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,                     \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                                      \
      Kokkos::View<SCALAR**, Kokkos::LayoutLeft, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,                           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                             \
      true, true> {                                                                                                      \
    using device_type = Kokkos::Device<EXECSPACE, Kokkos::HostSpace>;                                                    \
    using Handle      = KokkosSparse::Impl::SPMVHandleImpl<EXECSPACE, Kokkos::HostSpace, SCALAR, MKL_INT, MKL_INT>;      \
    using AMatrix     = ::KokkosSparse::Experimental::BsrMatrix<SCALAR const, MKL_INT const, device_type,                \
                                                            Kokkos::MemoryTraits<Kokkos::Unmanaged>, MKL_INT const>; \
    using XVector     = Kokkos::View<SCALAR const**, Kokkos::LayoutLeft, device_type,                                    \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;                    \
    using YVector = Kokkos::View<SCALAR**, Kokkos::LayoutLeft, device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;    \
    using coefficient_type = typename YVector::non_const_value_type;                                                     \
                                                                                                                         \
    static void spmv_mv_bsrmatrix(const EXECSPACE&, Handle* handle, const char mode[], const coefficient_type& alpha,    \
                                  const AMatrix& A, const XVector& X, const coefficient_type& beta,                      \
                                  const YVector& Y) {                                                                    \
      std::string label = "KokkosSparse::spmv_mv[TPL_MKL,BSRMATRIX," + Kokkos::ArithTraits<SCALAR>::name() + "]";        \
      Kokkos::Profiling::pushRegion(label);                                                                              \
      MKL_INT colx = static_cast<MKL_INT>(X.extent(1));                                                                  \
      MKL_INT ldx  = static_cast<MKL_INT>(X.stride_1());                                                                 \
      MKL_INT ldy  = static_cast<MKL_INT>(Y.stride_1());                                                                 \
      spmv_mv_bsr_mkl(handle, mode_kk_to_mkl(mode[0]), alpha, beta, A.numRows(), A.numCols(), A.blockDim(),              \
                      A.graph.row_map.data(), A.graph.entries.data(), A.values.data(), X.data(), colx, ldx, Y.data(),    \
                      ldy);                                                                                              \
      Kokkos::Profiling::popRegion();                                                                                    \
    }                                                                                                                    \
  };

#ifdef KOKKOS_ENABLE_SERIAL
KOKKOSSPARSE_SPMV_MV_MKL(float, Kokkos::Serial)
KOKKOSSPARSE_SPMV_MV_MKL(double, Kokkos::Serial)
KOKKOSSPARSE_SPMV_MV_MKL(Kokkos::complex<float>, Kokkos::Serial)
KOKKOSSPARSE_SPMV_MV_MKL(Kokkos::complex<double>, Kokkos::Serial)
#endif

#ifdef KOKKOS_ENABLE_OPENMP
KOKKOSSPARSE_SPMV_MV_MKL(float, Kokkos::OpenMP)
KOKKOSSPARSE_SPMV_MV_MKL(double, Kokkos::OpenMP)
KOKKOSSPARSE_SPMV_MV_MKL(Kokkos::complex<float>, Kokkos::OpenMP)
KOKKOSSPARSE_SPMV_MV_MKL(Kokkos::complex<double>, Kokkos::OpenMP)
#endif

#undef KOKKOSSPARSE_SPMV_MV_MKL

}  // namespace Impl
}  // namespace KokkosSparse

#endif  // defined(KOKKOSKERNELS_ENABLE_TPL_MKL) && (__INTEL_MKL__ > 2017)

// cuSPARSE
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#include "cusparse.h"

//
// From  https://docs.nvidia.com/cuda/cusparse/index.html#bsrmv
// Several comments on bsrmv():
// - Only blockDim > 1 is supported
// - Only CUSPARSE_OPERATION_NON_TRANSPOSE is supported
// - Only CUSPARSE_MATRIX_TYPE_GENERAL is supported.
//
#if (9000 <= CUDA_VERSION)

#include "KokkosSparse_Utils_cusparse.hpp"

namespace KokkosSparse {
namespace Impl {

template <class Handle, class AMatrix, class XVector, class YVector>
void spmv_bsr_cusparse(const Kokkos::Cuda& exec, Handle* handle, const char mode[],
                       typename YVector::const_value_type& alpha, const AMatrix& A, const XVector& x,
                       typename YVector::const_value_type& beta, const YVector& y) {
  using offset_type = typename AMatrix::non_const_size_type;
  using entry_type  = typename AMatrix::non_const_ordinal_type;
  using value_type  = typename AMatrix::non_const_value_type;

  /* initialize cusparse library */
  cusparseHandle_t cusparseHandle = KokkosKernels::Impl::CusparseSingleton::singleton().cusparseHandle;
  /* Set cuSPARSE to use the given stream until this function exits */
  KokkosSparse::Impl::TemporarySetCusparseStream tscs(cusparseHandle, exec);

  /* Set the operation mode */
  cusparseOperation_t myCusparseOperation;
  switch (toupper(mode[0])) {
    case 'N': myCusparseOperation = CUSPARSE_OPERATION_NON_TRANSPOSE; break;
    default: {
      std::cerr << "Mode " << mode << " invalid for cusparse[*]bsrmv.\n";
      throw std::invalid_argument("Invalid mode");
    }
  }

  KokkosSparse::Impl::CuSparse9_SpMV_Data* subhandle;

  if (handle->tpl_rank1) {
    subhandle = dynamic_cast<KokkosSparse::Impl::CuSparse9_SpMV_Data*>(handle->tpl_rank1);
    if (!subhandle) throw std::runtime_error("KokkosSparse::spmv: subhandle is not set up for cusparse");
    subhandle->set_exec_space(exec);
  } else {
    /* create and set the subhandle and matrix descriptor */
    subhandle         = new KokkosSparse::Impl::CuSparse9_SpMV_Data(exec);
    handle->tpl_rank1 = subhandle;
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateMatDescr(&subhandle->mat));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetMatType(subhandle->mat, CUSPARSE_MATRIX_TYPE_GENERAL));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetMatIndexBase(subhandle->mat, CUSPARSE_INDEX_BASE_ZERO));
  }

  cusparseDirection_t dirA = CUSPARSE_DIRECTION_ROW;

  /* perform the actual SpMV operation */
  static_assert(std::is_same_v<int, offset_type> && std::is_same_v<int, entry_type>,
                "With cuSPARSE non-generic API, offset and entry types must both be int. "
                "Something wrong with TPL avail logic.");
  if constexpr (std::is_same_v<value_type, float>) {
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSbsrmv(
        cusparseHandle, dirA, myCusparseOperation, A.numRows(), A.numCols(), A.nnz(),
        reinterpret_cast<float const*>(&alpha), subhandle->mat, reinterpret_cast<float const*>(A.values.data()),
        A.graph.row_map.data(), A.graph.entries.data(), A.blockDim(), reinterpret_cast<float const*>(x.data()),
        reinterpret_cast<float const*>(&beta), reinterpret_cast<float*>(y.data())));
  } else if constexpr (std::is_same_v<value_type, double>) {
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseDbsrmv(
        cusparseHandle, dirA, myCusparseOperation, A.numRows(), A.numCols(), A.nnz(),
        reinterpret_cast<double const*>(&alpha), subhandle->mat, reinterpret_cast<double const*>(A.values.data()),
        A.graph.row_map.data(), A.graph.entries.data(), A.blockDim(), reinterpret_cast<double const*>(x.data()),
        reinterpret_cast<double const*>(&beta), reinterpret_cast<double*>(y.data())));
  } else if constexpr (std::is_same_v<value_type, Kokkos::complex<float>>) {
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCbsrmv(
        cusparseHandle, dirA, myCusparseOperation, A.numRows(), A.numCols(), A.nnz(),
        reinterpret_cast<cuComplex const*>(&alpha), subhandle->mat, reinterpret_cast<cuComplex const*>(A.values.data()),
        A.graph.row_map.data(), A.graph.entries.data(), A.blockDim(), reinterpret_cast<cuComplex const*>(x.data()),
        reinterpret_cast<cuComplex const*>(&beta), reinterpret_cast<cuComplex*>(y.data())));
  } else if constexpr (std::is_same_v<value_type, Kokkos::complex<double>>) {
    KOKKOS_CUSPARSE_SAFE_CALL(
        cusparseZbsrmv(cusparseHandle, dirA, myCusparseOperation, A.numRows(), A.numCols(), A.nnz(),
                       reinterpret_cast<cuDoubleComplex const*>(&alpha), subhandle->mat,
                       reinterpret_cast<cuDoubleComplex const*>(A.values.data()), A.graph.row_map.data(),
                       A.graph.entries.data(), A.blockDim(), reinterpret_cast<cuDoubleComplex const*>(x.data()),
                       reinterpret_cast<cuDoubleComplex const*>(&beta), reinterpret_cast<cuDoubleComplex*>(y.data())));
  } else {
    static_assert(KokkosKernels::Impl::always_false_v<value_type>,
                  "Trying to call cusparse[*]bsrmv with a scalar type not "
                  "float/double, nor complex of either!");
  }
}

// Reference
// https://docs.nvidia.com/cuda/cusparse/index.html#bsrmm
// Several comments on bsrmm():
// - Only blockDim > 1 is supported
// - Only CUSPARSE_OPERATION_NON_TRANSPOSE is supported
// - Only CUSPARSE_MATRIX_TYPE_GENERAL is supported.
// - Only LayoutLeft for X and Y:
//   for X,Y LayoutLeft we want cuSparse to do
//   C = A * B + C
//   and for X,Y LayoutRight we want cuSparse to do
//   trans(C) = A * trans(B) + trans(C)
//   -> t(t(C)) = t(A * t(B)) + t(t(C))
//   ->       C = t(t(B)) * t(A) + C
//   ->       C = B * t(A) + C
//   This is impossible in cuSparse without explicitly transposing A,
//   so we just do not support LayoutRight in cuSparse TPL now (this is
//   statically asserted here)
template <class Handle, class AMatrix, class XVector, class YVector>
void spmv_mv_bsr_cusparse(const Kokkos::Cuda& exec, Handle* handle, const char mode[],
                          typename YVector::const_value_type& alpha, const AMatrix& A, const XVector& x,
                          typename YVector::const_value_type& beta, const YVector& y) {
  using offset_type = typename AMatrix::non_const_size_type;
  using entry_type  = typename AMatrix::non_const_ordinal_type;
  using value_type  = typename AMatrix::non_const_value_type;

  /* initialize cusparse library */
  cusparseHandle_t cusparseHandle = KokkosKernels::Impl::CusparseSingleton::singleton().cusparseHandle;
  /* Set cuSPARSE to use the given stream until this function exits */
  KokkosSparse::Impl::TemporarySetCusparseStream tscs(cusparseHandle, exec);

  /* Set the operation mode */
  cusparseOperation_t myCusparseOperation;
  switch (toupper(mode[0])) {
    case 'N': myCusparseOperation = CUSPARSE_OPERATION_NON_TRANSPOSE; break;
    default: {
      std::cerr << "Mode " << mode << " invalid for cusparse[*]bsrmv.\n";
      throw std::invalid_argument("Invalid mode");
    }
  }

  int colx = static_cast<int>(x.extent(1));

  // ldx and ldy should be the leading dimension (stride between columns) of X,Y
  // respectively
  const int ldx = static_cast<int>(x.stride(1));
  const int ldy = static_cast<int>(y.stride(1));

  static_assert(std::is_same_v<typename XVector::array_layout, Kokkos::LayoutLeft> &&
                    std::is_same_v<typename YVector::array_layout, Kokkos::LayoutLeft>,
                "cuSPARSE requires both X and Y to be LayoutLeft.");

  KokkosSparse::Impl::CuSparse9_SpMV_Data* subhandle;

  if (handle->tpl_rank2) {
    subhandle = dynamic_cast<KokkosSparse::Impl::CuSparse9_SpMV_Data*>(handle->tpl_rank2);
    if (!subhandle) throw std::runtime_error("KokkosSparse::spmv: subhandle is not set up for cusparse");
    subhandle->set_exec_space(exec);
  } else {
    /* create and set the subhandle and matrix descriptor */
    subhandle         = new KokkosSparse::Impl::CuSparse9_SpMV_Data(exec);
    handle->tpl_rank2 = subhandle;
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateMatDescr(&subhandle->mat));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetMatType(subhandle->mat, CUSPARSE_MATRIX_TYPE_GENERAL));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetMatIndexBase(subhandle->mat, CUSPARSE_INDEX_BASE_ZERO));
  }
  cusparseDirection_t dirA = CUSPARSE_DIRECTION_ROW;

  /* perform the actual SpMV operation */
  static_assert(std::is_same_v<int, offset_type> && std::is_same_v<int, entry_type>,
                "With cuSPARSE non-generic API, offset and entry types must both be int. "
                "Something wrong with TPL avail logic.");
  if constexpr (std::is_same_v<value_type, float>) {
    KOKKOS_CUSPARSE_SAFE_CALL(
        cusparseSbsrmm(cusparseHandle, dirA, myCusparseOperation, CUSPARSE_OPERATION_NON_TRANSPOSE, A.numRows(), colx,
                       A.numCols(), A.nnz(), reinterpret_cast<float const*>(&alpha), subhandle->mat,
                       reinterpret_cast<float const*>(A.values.data()), A.graph.row_map.data(), A.graph.entries.data(),
                       A.blockDim(), reinterpret_cast<float const*>(x.data()), ldx,
                       reinterpret_cast<float const*>(&beta), reinterpret_cast<float*>(y.data()), ldy));
  } else if constexpr (std::is_same_v<value_type, double>) {
    KOKKOS_CUSPARSE_SAFE_CALL(
        cusparseDbsrmm(cusparseHandle, dirA, myCusparseOperation, CUSPARSE_OPERATION_NON_TRANSPOSE, A.numRows(), colx,
                       A.numCols(), A.nnz(), reinterpret_cast<double const*>(&alpha), subhandle->mat,
                       reinterpret_cast<double const*>(A.values.data()), A.graph.row_map.data(), A.graph.entries.data(),
                       A.blockDim(), reinterpret_cast<double const*>(x.data()), ldx,
                       reinterpret_cast<double const*>(&beta), reinterpret_cast<double*>(y.data()), ldy));
  } else if constexpr (std::is_same_v<value_type, Kokkos::complex<float>>) {
    KOKKOS_CUSPARSE_SAFE_CALL(
        cusparseCbsrmm(cusparseHandle, dirA, myCusparseOperation, CUSPARSE_OPERATION_NON_TRANSPOSE, A.numRows(), colx,
                       A.numCols(), A.nnz(), reinterpret_cast<cuComplex const*>(&alpha), subhandle->mat,
                       reinterpret_cast<cuComplex const*>(A.values.data()), A.graph.row_map.data(),
                       A.graph.entries.data(), A.blockDim(), reinterpret_cast<cuComplex const*>(x.data()), ldx,
                       reinterpret_cast<cuComplex const*>(&beta), reinterpret_cast<cuComplex*>(y.data()), ldy));
  } else if constexpr (std::is_same_v<value_type, Kokkos::complex<double>>) {
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseZbsrmm(
        cusparseHandle, dirA, myCusparseOperation, CUSPARSE_OPERATION_NON_TRANSPOSE, A.numRows(), colx, A.numCols(),
        A.nnz(), reinterpret_cast<cuDoubleComplex const*>(&alpha), subhandle->mat,
        reinterpret_cast<cuDoubleComplex const*>(A.values.data()), A.graph.row_map.data(), A.graph.entries.data(),
        A.blockDim(), reinterpret_cast<cuDoubleComplex const*>(x.data()), ldx,
        reinterpret_cast<cuDoubleComplex const*>(&beta), reinterpret_cast<cuDoubleComplex*>(y.data()), ldy));
  } else {
    static_assert(KokkosKernels::Impl::always_false_v<value_type>,
                  "Trying to call cusparse[*]bsrmm with a scalar type not "
                  "float/double, nor complex of either!");
  }
}

#define KOKKOSSPARSE_SPMV_CUSPARSE(SCALAR, ORDINAL, OFFSET, LAYOUT, SPACE)                                          \
  template <>                                                                                                       \
  struct SPMV_BSRMATRIX<                                                                                            \
      Kokkos::Cuda, KokkosSparse::Impl::SPMVHandleImpl<Kokkos::Cuda, SPACE, SCALAR, OFFSET, ORDINAL>,               \
      ::KokkosSparse::Experimental::BsrMatrix<SCALAR const, ORDINAL const, Kokkos::Device<Kokkos::Cuda, SPACE>,     \
                                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, OFFSET const>,               \
      Kokkos::View<SCALAR const*, LAYOUT, Kokkos::Device<Kokkos::Cuda, SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                                 \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<Kokkos::Cuda, SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,  \
      true> {                                                                                                       \
    using device_type       = Kokkos::Device<Kokkos::Cuda, SPACE>;                                                  \
    using memory_trait_type = Kokkos::MemoryTraits<Kokkos::Unmanaged>;                                              \
    using Handle            = KokkosSparse::Impl::SPMVHandleImpl<Kokkos::Cuda, SPACE, SCALAR, OFFSET, ORDINAL>;     \
    using AMatrix           = ::KokkosSparse::Experimental::BsrMatrix<SCALAR const, ORDINAL const, device_type,     \
                                                            memory_trait_type, OFFSET const>;             \
    using XVector           = Kokkos::View<SCALAR const*, LAYOUT, device_type,                                      \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;         \
    using YVector           = Kokkos::View<SCALAR*, LAYOUT, device_type, memory_trait_type>;                        \
                                                                                                                    \
    using coefficient_type = typename YVector::non_const_value_type;                                                \
                                                                                                                    \
    static void spmv_bsrmatrix(const Kokkos::Cuda& exec, Handle* handle, const char mode[],                         \
                               const coefficient_type& alpha, const AMatrix& A, const XVector& x,                   \
                               const coefficient_type& beta, const YVector& y) {                                    \
      std::string label = "KokkosSparse::spmv[TPL_CUSPARSE,BSRMATRIX," + Kokkos::ArithTraits<SCALAR>::name() + "]"; \
      Kokkos::Profiling::pushRegion(label);                                                                         \
      spmv_bsr_cusparse(exec, handle, mode, alpha, A, x, beta, y);                                                  \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };

KOKKOSSPARSE_SPMV_CUSPARSE(double, int, int, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_CUSPARSE(double, int, int, Kokkos::LayoutRight, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_CUSPARSE(float, int, int, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_CUSPARSE(float, int, int, Kokkos::LayoutRight, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutRight, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutRight, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_CUSPARSE(double, int, int, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_CUSPARSE(double, int, int, Kokkos::LayoutRight, Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_CUSPARSE(float, int, int, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_CUSPARSE(float, int, int, Kokkos::LayoutRight, Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutRight, Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutRight, Kokkos::CudaUVMSpace)

#undef KOKKOSSPARSE_SPMV_CUSPARSE

// cuSparse TPL does not support LayoutRight for this operation
// only specialize for LayoutLeft
#define KOKKOSSPARSE_SPMV_MV_CUSPARSE(SCALAR, ORDINAL, OFFSET, SPACE, ETI_AVAIL)                                       \
  template <>                                                                                                          \
  struct SPMV_MV_BSRMATRIX<                                                                                            \
      Kokkos::Cuda, KokkosSparse::Impl::SPMVHandleImpl<Kokkos::Cuda, SPACE, SCALAR, OFFSET, ORDINAL>,                  \
      ::KokkosSparse::Experimental::BsrMatrix<SCALAR const, ORDINAL const, Kokkos::Device<Kokkos::Cuda, SPACE>,        \
                                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, OFFSET const>,                  \
      Kokkos::View<SCALAR const**, Kokkos::LayoutLeft, Kokkos::Device<Kokkos::Cuda, SPACE>,                            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                                    \
      Kokkos::View<SCALAR**, Kokkos::LayoutLeft, Kokkos::Device<Kokkos::Cuda, SPACE>,                                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      false, true, ETI_AVAIL> {                                                                                        \
    using device_type       = Kokkos::Device<Kokkos::Cuda, SPACE>;                                                     \
    using memory_trait_type = Kokkos::MemoryTraits<Kokkos::Unmanaged>;                                                 \
    using Handle            = KokkosSparse::Impl::SPMVHandleImpl<Kokkos::Cuda, SPACE, SCALAR, OFFSET, ORDINAL>;        \
    using AMatrix           = ::KokkosSparse::Experimental::BsrMatrix<SCALAR const, ORDINAL const, device_type,        \
                                                            memory_trait_type, OFFSET const>;                \
    using XVector           = Kokkos::View<SCALAR const**, Kokkos::LayoutLeft, device_type,                            \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;            \
    using YVector           = Kokkos::View<SCALAR**, Kokkos::LayoutLeft, device_type, memory_trait_type>;              \
                                                                                                                       \
    using coefficient_type = typename YVector::non_const_value_type;                                                   \
                                                                                                                       \
    static void spmv_mv_bsrmatrix(const Kokkos::Cuda& exec, Handle* handle, const char mode[],                         \
                                  const coefficient_type& alpha, const AMatrix& A, const XVector& x,                   \
                                  const coefficient_type& beta, const YVector& y) {                                    \
      std::string label = "KokkosSparse::spmv_mv[TPL_CUSPARSE,BSRMATRIX," + Kokkos::ArithTraits<SCALAR>::name() + "]"; \
      Kokkos::Profiling::pushRegion(label);                                                                            \
      spmv_mv_bsr_cusparse(exec, handle, mode, alpha, A, x, beta, y);                                                  \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSSPARSE_SPMV_MV_CUSPARSE(double, int, int, Kokkos::CudaSpace, true)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(double, int, int, Kokkos::CudaSpace, false)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(float, int, int, Kokkos::CudaSpace, true)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(float, int, int, Kokkos::CudaSpace, false)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::CudaSpace, true)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::CudaSpace, false)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::CudaSpace, true)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::CudaSpace, false)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(double, int, int, Kokkos::CudaUVMSpace, true)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(double, int, int, Kokkos::CudaUVMSpace, false)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(float, int, int, Kokkos::CudaUVMSpace, true)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(float, int, int, Kokkos::CudaUVMSpace, false)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::CudaUVMSpace, true)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::CudaUVMSpace, false)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::CudaUVMSpace, true)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::CudaUVMSpace, false)

#undef KOKKOSSPARSE_SPMV_MV_CUSPARSE

}  // namespace Impl
}  // namespace KokkosSparse
#endif  // (9000 <= CUDA_VERSION)

#endif  // KOKKOSKERNELS_ENABLE_TPL_CUSPARSE

// --------------------
// rocSparse
// --------------------
#if defined(KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE)

#include <rocsparse/rocsparse.h>

#include "KokkosSparse_Utils_rocsparse.hpp"

namespace KokkosSparse {
namespace Impl {

template <class Handle, class AMatrix, class XVector, class YVector>
void spmv_bsr_rocsparse(const Kokkos::HIP& exec, Handle* handle, const char mode[],
                        typename YVector::const_value_type& alpha, const AMatrix& A, const XVector& x,
                        typename YVector::const_value_type& beta, const YVector& y) {
  /*
     rocm 5.4.0 rocsparse_*bsrmv reference:
     https://rocsparse.readthedocs.io/en/rocm-5.4.0/usermanual.html#rocsparse-bsrmv-ex

     only trans = rocsparse_operation_none is supported
     only descr = rocsparse_matrix_type_general is supported

  */

  using offset_type          = typename AMatrix::non_const_size_type;
  using ordinal_type         = typename AMatrix::non_const_ordinal_type;
  using value_type           = typename AMatrix::non_const_value_type;
  using rocsparse_value_type = typename KokkosSparse::Impl::kokkos_to_rocsparse_type<value_type>::type;

  // assert ordinals and offsets are the expected types
  static_assert(std::is_same_v<offset_type, rocsparse_int>, "A offset_type must be rocsparse_int");
  static_assert(std::is_same_v<ordinal_type, rocsparse_int>, "A ordinal_type must be rocsparse_int");

  // assert all operands are the same type
  using x_value_type = typename XVector::non_const_value_type;
  using y_value_type = typename YVector::non_const_value_type;
  static_assert(std::is_same_v<value_type, x_value_type>, "A and x must have same value type");
  static_assert(std::is_same_v<value_type, y_value_type>, "A and y must have same value type");

  // assert X and Y are non-stride (pass raw pointers to TPL)
  static_assert(!std::is_same_v<typename XVector::array_layout, Kokkos::LayoutStride>, "x must be contiguous");
  static_assert(!std::is_same_v<typename YVector::array_layout, Kokkos::LayoutStride>, "y must be contiguous");

  // assert BSR data is non-stride (pass raw pointers to TPL)
  static_assert(!std::is_same_v<typename AMatrix::values_type::array_layout, Kokkos::LayoutStride>,
                "A values must be contiguous");
  static_assert(!std::is_same_v<typename AMatrix::row_map_type::array_layout, Kokkos::LayoutStride>,
                "A row_map must be contiguous");
  static_assert(!std::is_same_v<typename AMatrix::index_type::array_layout, Kokkos::LayoutStride>,
                "A entries must be contiguous");

  rocsparse_handle rocsparseHandle = KokkosKernels::Impl::RocsparseSingleton::singleton().rocsparseHandle;
  // resets handle stream to NULL when out of scope
  KokkosSparse::Impl::TemporarySetRocsparseStream tsrs(rocsparseHandle, exec);

  // set the mode
  rocsparse_operation trans;
  switch (toupper(mode[0])) {
    case 'N': trans = rocsparse_operation_none; break;
    default: {
      std::stringstream ss;
      ss << "Mode " << mode << " invalid for rocsparse_[*]bsrmv\n";
      throw std::invalid_argument(ss.str());
    }
  }

  /*
  Specify the matrix direction.
  The rocsparse_direction indicates whether a dense matrix should be parsed by
  rows or by columns, assuming column-major storage. Values: enumerator
  rocsparse_direction_row Parse the matrix by rows. enumerator
  rocsparse_direction_column Parse the matrix by columns.
  */
  // KokkosSparse Bsr matrix blocks are layoutright (row-major)
  static_assert(std::is_same_v<typename AMatrix::block_layout_type, Kokkos::LayoutRight>,
                "A blocks must be stored layout-right");
  rocsparse_direction dir = rocsparse_direction_row;

  const rocsparse_int mb             = rocsparse_int(A.numRows());  // number of block rows
  const rocsparse_int nb             = rocsparse_int(A.numCols());  // number of block cols
  const rocsparse_int nnzb           = rocsparse_int(A.nnz());      // number of non-zero blocks
  const rocsparse_value_type* alpha_ = reinterpret_cast<const rocsparse_value_type*>(&alpha);

  const rocsparse_value_type* bsr_val = reinterpret_cast<const rocsparse_value_type*>(A.values.data());
  const rocsparse_int* bsr_row_ptr    = A.graph.row_map.data();
  const rocsparse_int* bsr_col_ind    = A.graph.entries.data();
  const rocsparse_int block_dim       = rocsparse_int(A.blockDim());
  const rocsparse_value_type* x_      = reinterpret_cast<const rocsparse_value_type*>(x.data());
  const rocsparse_value_type* beta_   = reinterpret_cast<const rocsparse_value_type*>(&beta);
  rocsparse_value_type* y_            = reinterpret_cast<rocsparse_value_type*>(y.data());

  KokkosSparse::Impl::RocSparse_BSR_SpMV_Data* subhandle;
  if (handle->tpl_rank1) {
    subhandle = dynamic_cast<KokkosSparse::Impl::RocSparse_BSR_SpMV_Data*>(handle->tpl_rank1);
    if (!subhandle) throw std::runtime_error("KokkosSparse::spmv: subhandle is not set up for rocsparse BSR");
    subhandle->set_exec_space(exec);
  } else {
    subhandle         = new KokkosSparse::Impl::RocSparse_BSR_SpMV_Data(exec);
    handle->tpl_rank1 = subhandle;
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_create_mat_descr(&subhandle->mat));
    // *_ex* functions deprecated in introduced in 6+
#if KOKKOSSPARSE_IMPL_ROCM_VERSION >= 60000
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_create_mat_info(&subhandle->info));
    if constexpr (std::is_same_v<value_type, float>) {
      KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_sbsrmv_analysis(rocsparseHandle, dir, trans, mb, nb, nnzb,
                                                                subhandle->mat, bsr_val, bsr_row_ptr, bsr_col_ind,
                                                                block_dim, subhandle->info));
    } else if constexpr (std::is_same_v<value_type, double>) {
      KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_dbsrmv_analysis(rocsparseHandle, dir, trans, mb, nb, nnzb,
                                                                subhandle->mat, bsr_val, bsr_row_ptr, bsr_col_ind,
                                                                block_dim, subhandle->info));
    } else if constexpr (std::is_same_v<value_type, Kokkos::complex<float>>) {
      KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_cbsrmv_analysis(rocsparseHandle, dir, trans, mb, nb, nnzb,
                                                                subhandle->mat, bsr_val, bsr_row_ptr, bsr_col_ind,
                                                                block_dim, subhandle->info));
    } else if constexpr (std::is_same_v<value_type, Kokkos::complex<double>>) {
      KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_zbsrmv_analysis(rocsparseHandle, dir, trans, mb, nb, nnzb,
                                                                subhandle->mat, bsr_val, bsr_row_ptr, bsr_col_ind,
                                                                block_dim, subhandle->info));
    } else {
      static_assert(KokkosKernels::Impl::always_false_v<value_type>, "unsupported value type for rocsparse_*bsrmv");
    }
    // *_ex* functions introduced in 5.4.0
#elif KOKKOSSPARSE_IMPL_ROCM_VERSION < 50400
    // No analysis step in the older versions
#else
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_create_mat_info(&subhandle->info));
    if constexpr (std::is_same_v<value_type, float>) {
      KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_sbsrmv_ex_analysis(rocsparseHandle, dir, trans, mb, nb, nnzb,
                                                                   subhandle->mat, bsr_val, bsr_row_ptr, bsr_col_ind,
                                                                   block_dim, subhandle->info));
    } else if constexpr (std::is_same_v<value_type, double>) {
      KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_dbsrmv_ex_analysis(rocsparseHandle, dir, trans, mb, nb, nnzb,
                                                                   subhandle->mat, bsr_val, bsr_row_ptr, bsr_col_ind,
                                                                   block_dim, subhandle->info));
    } else if constexpr (std::is_same_v<value_type, Kokkos::complex<float>>) {
      KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_cbsrmv_ex_analysis(rocsparseHandle, dir, trans, mb, nb, nnzb,
                                                                   subhandle->mat, bsr_val, bsr_row_ptr, bsr_col_ind,
                                                                   block_dim, subhandle->info));
    } else if constexpr (std::is_same_v<value_type, Kokkos::complex<double>>) {
      KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_zbsrmv_ex_analysis(rocsparseHandle, dir, trans, mb, nb, nnzb,
                                                                   subhandle->mat, bsr_val, bsr_row_ptr, bsr_col_ind,
                                                                   block_dim, subhandle->info));
    } else {
      static_assert(KokkosKernels::Impl::always_false_v<value_type>, "unsupported value type for rocsparse_*bsrmv");
    }
#endif
  }

  // *_ex* functions deprecated in introduced in 6+
#if KOKKOSSPARSE_IMPL_ROCM_VERSION >= 60000
  if constexpr (std::is_same_v<value_type, float>) {
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_sbsrmv(rocsparseHandle, dir, trans, mb, nb, nnzb, alpha_, subhandle->mat,
                                                     bsr_val, bsr_row_ptr, bsr_col_ind, block_dim, subhandle->info, x_,
                                                     beta_, y_));
  } else if constexpr (std::is_same_v<value_type, double>) {
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_dbsrmv(rocsparseHandle, dir, trans, mb, nb, nnzb, alpha_, subhandle->mat,
                                                     bsr_val, bsr_row_ptr, bsr_col_ind, block_dim, subhandle->info, x_,
                                                     beta_, y_));
  } else if constexpr (std::is_same_v<value_type, Kokkos::complex<float>>) {
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_cbsrmv(rocsparseHandle, dir, trans, mb, nb, nnzb, alpha_, subhandle->mat,
                                                     bsr_val, bsr_row_ptr, bsr_col_ind, block_dim, subhandle->info, x_,
                                                     beta_, y_));
  } else if constexpr (std::is_same_v<value_type, Kokkos::complex<double>>) {
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_zbsrmv(rocsparseHandle, dir, trans, mb, nb, nnzb, alpha_, subhandle->mat,
                                                     bsr_val, bsr_row_ptr, bsr_col_ind, block_dim, subhandle->info, x_,
                                                     beta_, y_));
  } else {
    static_assert(KokkosKernels::Impl::always_false_v<value_type>, "unsupported value type for rocsparse_*bsrmv");
  }
  // *_ex* functions introduced in 5.4.0
#elif KOKKOSSPARSE_IMPL_ROCM_VERSION < 50400
  if constexpr (std::is_same_v<value_type, float>) {
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_sbsrmv(rocsparseHandle, dir, trans, mb, nb, nnzb, alpha_, subhandle->mat,
                                                     bsr_val, bsr_row_ptr, bsr_col_ind, block_dim, x_, beta_, y_));
  } else if constexpr (std::is_same_v<value_type, double>) {
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_dbsrmv(rocsparseHandle, dir, trans, mb, nb, nnzb, alpha_, subhandle->mat,
                                                     bsr_val, bsr_row_ptr, bsr_col_ind, block_dim, x_, beta_, y_));
  } else if constexpr (std::is_same_v<value_type, Kokkos::complex<float>>) {
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_cbsrmv(rocsparseHandle, dir, trans, mb, nb, nnzb, alpha_, subhandle->mat,
                                                     bsr_val, bsr_row_ptr, bsr_col_ind, block_dim, x_, beta_, y_));
  } else if constexpr (std::is_same_v<value_type, Kokkos::complex<double>>) {
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_zbsrmv(rocsparseHandle, dir, trans, mb, nb, nnzb, alpha_, subhandle->mat,
                                                     bsr_val, bsr_row_ptr, bsr_col_ind, block_dim, x_, beta_, y_));
  } else {
    static_assert(KokkosKernels::Impl::always_false_v<value_type>, "unsupported value type for rocsparse_*bsrmv");
  }
#else
  if constexpr (std::is_same_v<value_type, float>) {
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_sbsrmv_ex(rocsparseHandle, dir, trans, mb, nb, nnzb, alpha_,
                                                        subhandle->mat, bsr_val, bsr_row_ptr, bsr_col_ind, block_dim,
                                                        subhandle->info, x_, beta_, y_));
  } else if constexpr (std::is_same_v<value_type, double>) {
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_dbsrmv_ex(rocsparseHandle, dir, trans, mb, nb, nnzb, alpha_,
                                                        subhandle->mat, bsr_val, bsr_row_ptr, bsr_col_ind, block_dim,
                                                        subhandle->info, x_, beta_, y_));
  } else if constexpr (std::is_same_v<value_type, Kokkos::complex<float>>) {
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_cbsrmv_ex(rocsparseHandle, dir, trans, mb, nb, nnzb, alpha_,
                                                        subhandle->mat, bsr_val, bsr_row_ptr, bsr_col_ind, block_dim,
                                                        subhandle->info, x_, beta_, y_));
  } else if constexpr (std::is_same_v<value_type, Kokkos::complex<double>>) {
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_zbsrmv_ex(rocsparseHandle, dir, trans, mb, nb, nnzb, alpha_,
                                                        subhandle->mat, bsr_val, bsr_row_ptr, bsr_col_ind, block_dim,
                                                        subhandle->info, x_, beta_, y_));
  } else {
    static_assert(KokkosKernels::Impl::always_false_v<value_type>, "unsupported value type for rocsparse_*bsrmv");
  }
#endif
}  // spmv_bsr_rocsparse

#define KOKKOSSPARSE_SPMV_ROCSPARSE(SCALAR, ORDINAL, OFFSET, LAYOUT, SPACE)                                          \
  template <>                                                                                                        \
  struct SPMV_BSRMATRIX<                                                                                             \
      Kokkos::HIP, KokkosSparse::Impl::SPMVHandleImpl<Kokkos::HIP, SPACE, SCALAR, OFFSET, ORDINAL>,                  \
      ::KokkosSparse::Experimental::BsrMatrix<SCALAR const, ORDINAL const, Kokkos::Device<Kokkos::HIP, SPACE>,       \
                                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, OFFSET const>,                \
      Kokkos::View<SCALAR const*, LAYOUT, Kokkos::Device<Kokkos::HIP, SPACE>,                                        \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                                  \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<Kokkos::HIP, SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,    \
      true> {                                                                                                        \
    using device_type       = Kokkos::Device<Kokkos::HIP, SPACE>;                                                    \
    using memory_trait_type = Kokkos::MemoryTraits<Kokkos::Unmanaged>;                                               \
    using Handle            = KokkosSparse::Impl::SPMVHandleImpl<Kokkos::HIP, SPACE, SCALAR, OFFSET, ORDINAL>;       \
    using AMatrix           = ::KokkosSparse::Experimental::BsrMatrix<SCALAR const, ORDINAL const, device_type,      \
                                                            memory_trait_type, OFFSET const>;              \
    using XVector           = Kokkos::View<SCALAR const*, LAYOUT, device_type,                                       \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;          \
    using YVector           = Kokkos::View<SCALAR*, LAYOUT, device_type, memory_trait_type>;                         \
                                                                                                                     \
    using coefficient_type = typename YVector::non_const_value_type;                                                 \
                                                                                                                     \
    static void spmv_bsrmatrix(const Kokkos::HIP& exec, Handle* handle, const char mode[],                           \
                               const coefficient_type& alpha, const AMatrix& A, const XVector& x,                    \
                               const coefficient_type& beta, const YVector& y) {                                     \
      std::string label = "KokkosSparse::spmv[TPL_ROCSPARSE,BSRMATRIX," + Kokkos::ArithTraits<SCALAR>::name() + "]"; \
      Kokkos::Profiling::pushRegion(label);                                                                          \
      spmv_bsr_rocsparse(exec, handle, mode, alpha, A, x, beta, y);                                                  \
      Kokkos::Profiling::popRegion();                                                                                \
    }                                                                                                                \
  };

KOKKOSSPARSE_SPMV_ROCSPARSE(float, rocsparse_int, rocsparse_int, Kokkos::LayoutLeft, Kokkos::HIPSpace);
KOKKOSSPARSE_SPMV_ROCSPARSE(float, rocsparse_int, rocsparse_int, Kokkos::LayoutRight, Kokkos::HIPSpace);
KOKKOSSPARSE_SPMV_ROCSPARSE(double, rocsparse_int, rocsparse_int, Kokkos::LayoutLeft, Kokkos::HIPSpace);
KOKKOSSPARSE_SPMV_ROCSPARSE(double, rocsparse_int, rocsparse_int, Kokkos::LayoutRight, Kokkos::HIPSpace);
KOKKOSSPARSE_SPMV_ROCSPARSE(Kokkos::complex<float>, rocsparse_int, rocsparse_int, Kokkos::LayoutLeft, Kokkos::HIPSpace);
KOKKOSSPARSE_SPMV_ROCSPARSE(Kokkos::complex<float>, rocsparse_int, rocsparse_int, Kokkos::LayoutRight,
                            Kokkos::HIPSpace);
KOKKOSSPARSE_SPMV_ROCSPARSE(Kokkos::complex<double>, rocsparse_int, rocsparse_int, Kokkos::LayoutLeft,
                            Kokkos::HIPSpace);
KOKKOSSPARSE_SPMV_ROCSPARSE(Kokkos::complex<double>, rocsparse_int, rocsparse_int, Kokkos::LayoutRight,
                            Kokkos::HIPSpace);

#undef KOKKOSSPARSE_SPMV_ROCSPARSE

}  // namespace Impl
}  // namespace KokkosSparse

#endif  // defined(KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE)

#endif  // KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_DECL_HPP
