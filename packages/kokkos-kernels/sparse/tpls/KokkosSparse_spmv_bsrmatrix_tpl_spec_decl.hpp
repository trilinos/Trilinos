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
#include "KokkosKernels_Controls.hpp"
#include "KokkosSparse_Utils_mkl.hpp"

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
#include <mkl.h>

namespace KokkosSparse {
namespace Experimental {
namespace Impl {

#if (__INTEL_MKL__ > 2017)
// MKL 2018 and above: use new interface: sparse_matrix_t and mkl_sparse_?_mv()

using KokkosSparse::Impl::mode_kk_to_mkl;

inline matrix_descr getDescription() {
  matrix_descr A_descr;
  A_descr.type = SPARSE_MATRIX_TYPE_GENERAL;
  A_descr.mode = SPARSE_FILL_MODE_FULL;
  A_descr.diag = SPARSE_DIAG_NON_UNIT;
  return A_descr;
}

inline void spmv_block_impl_mkl(sparse_operation_t op, float alpha, float beta,
                                MKL_INT m, MKL_INT n, MKL_INT b,
                                const MKL_INT* Arowptrs,
                                const MKL_INT* Aentries, const float* Avalues,
                                const float* x, float* y) {
  sparse_matrix_t A_mkl;
  KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_s_create_bsr(
      &A_mkl, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, m, n, b,
      const_cast<MKL_INT*>(Arowptrs), const_cast<MKL_INT*>(Arowptrs + 1),
      const_cast<MKL_INT*>(Aentries), const_cast<float*>(Avalues)));

  matrix_descr A_descr = getDescription();
  KOKKOSKERNELS_MKL_SAFE_CALL(
      mkl_sparse_s_mv(op, alpha, A_mkl, A_descr, x, beta, y));
}

inline void spmv_block_impl_mkl(sparse_operation_t op, double alpha,
                                double beta, MKL_INT m, MKL_INT n, MKL_INT b,
                                const MKL_INT* Arowptrs,
                                const MKL_INT* Aentries, const double* Avalues,
                                const double* x, double* y) {
  sparse_matrix_t A_mkl;
  KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_d_create_bsr(
      &A_mkl, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, m, n, b,
      const_cast<MKL_INT*>(Arowptrs), const_cast<MKL_INT*>(Arowptrs + 1),
      const_cast<MKL_INT*>(Aentries), const_cast<double*>(Avalues)));

  matrix_descr A_descr = getDescription();
  KOKKOSKERNELS_MKL_SAFE_CALL(
      mkl_sparse_d_mv(op, alpha, A_mkl, A_descr, x, beta, y));
}

inline void spmv_block_impl_mkl(sparse_operation_t op,
                                Kokkos::complex<float> alpha,
                                Kokkos::complex<float> beta, MKL_INT m,
                                MKL_INT n, MKL_INT b, const MKL_INT* Arowptrs,
                                const MKL_INT* Aentries,
                                const Kokkos::complex<float>* Avalues,
                                const Kokkos::complex<float>* x,
                                Kokkos::complex<float>* y) {
  sparse_matrix_t A_mkl;
  KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_c_create_bsr(
      &A_mkl, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, m, n, b,
      const_cast<MKL_INT*>(Arowptrs), const_cast<MKL_INT*>(Arowptrs + 1),
      const_cast<MKL_INT*>(Aentries), (MKL_Complex8*)Avalues));

  MKL_Complex8 alpha_mkl{alpha.real(), alpha.imag()};
  MKL_Complex8 beta_mkl{beta.real(), beta.imag()};
  matrix_descr A_descr = getDescription();
  KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_c_mv(
      op, alpha_mkl, A_mkl, A_descr, reinterpret_cast<const MKL_Complex8*>(x),
      beta_mkl, reinterpret_cast<MKL_Complex8*>(y)));
}

inline void spmv_block_impl_mkl(sparse_operation_t op,
                                Kokkos::complex<double> alpha,
                                Kokkos::complex<double> beta, MKL_INT m,
                                MKL_INT n, MKL_INT b, const MKL_INT* Arowptrs,
                                const MKL_INT* Aentries,
                                const Kokkos::complex<double>* Avalues,
                                const Kokkos::complex<double>* x,
                                Kokkos::complex<double>* y) {
  sparse_matrix_t A_mkl;
  KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_z_create_bsr(
      &A_mkl, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, m, n, b,
      const_cast<MKL_INT*>(Arowptrs), const_cast<MKL_INT*>(Arowptrs + 1),
      const_cast<MKL_INT*>(Aentries), (MKL_Complex16*)Avalues));

  matrix_descr A_descr = getDescription();
  MKL_Complex16 alpha_mkl{alpha.real(), alpha.imag()};
  MKL_Complex16 beta_mkl{beta.real(), beta.imag()};
  KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_z_mv(
      op, alpha_mkl, A_mkl, A_descr, reinterpret_cast<const MKL_Complex16*>(x),
      beta_mkl, reinterpret_cast<MKL_Complex16*>(y)));
}

inline void spm_mv_block_impl_mkl(sparse_operation_t op, float alpha,
                                  float beta, MKL_INT m, MKL_INT n, MKL_INT b,
                                  const MKL_INT* Arowptrs,
                                  const MKL_INT* Aentries, const float* Avalues,
                                  const float* x, MKL_INT colx, MKL_INT ldx,
                                  float* y, MKL_INT ldy) {
  sparse_matrix_t A_mkl;
  KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_s_create_bsr(
      &A_mkl, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, m, n, b,
      const_cast<MKL_INT*>(Arowptrs), const_cast<MKL_INT*>(Arowptrs + 1),
      const_cast<MKL_INT*>(Aentries), const_cast<float*>(Avalues)));

  matrix_descr A_descr = getDescription();
  KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_s_mm(op, alpha, A_mkl, A_descr,
                                              SPARSE_LAYOUT_ROW_MAJOR, x, colx,
                                              ldx, beta, y, ldy));
}

inline void spm_mv_block_impl_mkl(sparse_operation_t op, double alpha,
                                  double beta, MKL_INT m, MKL_INT n, MKL_INT b,
                                  const MKL_INT* Arowptrs,
                                  const MKL_INT* Aentries,
                                  const double* Avalues, const double* x,
                                  MKL_INT colx, MKL_INT ldx, double* y,
                                  MKL_INT ldy) {
  sparse_matrix_t A_mkl;
  KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_d_create_bsr(
      &A_mkl, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, m, n, b,
      const_cast<MKL_INT*>(Arowptrs), const_cast<MKL_INT*>(Arowptrs + 1),
      const_cast<MKL_INT*>(Aentries), const_cast<double*>(Avalues)));

  matrix_descr A_descr = getDescription();
  KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_d_mm(op, alpha, A_mkl, A_descr,
                                              SPARSE_LAYOUT_ROW_MAJOR, x, colx,
                                              ldx, beta, y, ldy));
}

inline void spm_mv_block_impl_mkl(
    sparse_operation_t op, Kokkos::complex<float> alpha,
    Kokkos::complex<float> beta, MKL_INT m, MKL_INT n, MKL_INT b,
    const MKL_INT* Arowptrs, const MKL_INT* Aentries,
    const Kokkos::complex<float>* Avalues, const Kokkos::complex<float>* x,
    MKL_INT colx, MKL_INT ldx, Kokkos::complex<float>* y, MKL_INT ldy) {
  sparse_matrix_t A_mkl;
  KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_c_create_bsr(
      &A_mkl, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, m, n, b,
      const_cast<MKL_INT*>(Arowptrs), const_cast<MKL_INT*>(Arowptrs + 1),
      const_cast<MKL_INT*>(Aentries), (MKL_Complex8*)Avalues));

  MKL_Complex8 alpha_mkl{alpha.real(), alpha.imag()};
  MKL_Complex8 beta_mkl{beta.real(), beta.imag()};
  matrix_descr A_descr = getDescription();
  KOKKOSKERNELS_MKL_SAFE_CALL(
      mkl_sparse_c_mm(op, alpha_mkl, A_mkl, A_descr, SPARSE_LAYOUT_ROW_MAJOR,
                      reinterpret_cast<const MKL_Complex8*>(x), colx, ldx,
                      beta_mkl, reinterpret_cast<MKL_Complex8*>(y), ldy));
}

inline void spm_mv_block_impl_mkl(
    sparse_operation_t op, Kokkos::complex<double> alpha,
    Kokkos::complex<double> beta, MKL_INT m, MKL_INT n, MKL_INT b,
    const MKL_INT* Arowptrs, const MKL_INT* Aentries,
    const Kokkos::complex<double>* Avalues, const Kokkos::complex<double>* x,
    MKL_INT colx, MKL_INT ldx, Kokkos::complex<double>* y, MKL_INT ldy) {
  sparse_matrix_t A_mkl;
  KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_z_create_bsr(
      &A_mkl, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, m, n, b,
      const_cast<MKL_INT*>(Arowptrs), const_cast<MKL_INT*>(Arowptrs + 1),
      const_cast<MKL_INT*>(Aentries), (MKL_Complex16*)Avalues));

  matrix_descr A_descr = getDescription();
  MKL_Complex16 alpha_mkl{alpha.real(), alpha.imag()};
  MKL_Complex16 beta_mkl{beta.real(), beta.imag()};
  KOKKOSKERNELS_MKL_SAFE_CALL(
      mkl_sparse_z_mm(op, alpha_mkl, A_mkl, A_descr, SPARSE_LAYOUT_ROW_MAJOR,
                      reinterpret_cast<const MKL_Complex16*>(x), colx, ldx,
                      beta_mkl, reinterpret_cast<MKL_Complex16*>(y), ldy));
}

#endif

#if (__INTEL_MKL__ == 2017)

inline void spmv_block_impl_mkl(char mode, float alpha, float beta, MKL_INT m,
                                MKL_INT n, MKL_INT b, const MKL_INT* Arowptrs,
                                const MKL_INT* Aentries, const float* Avalues,
                                const float* x, float* y) {
  mkl_sbsrmv(&mode, &m, &n, &b, &alpha, "G**C", Avalues, Aentries, Arowptrs,
             Arowptrs + 1, x, &beta, y);
}

inline void spmv_block_impl_mkl(char mode, double alpha, double beta, MKL_INT m,
                                MKL_INT n, MKL_INT b, const MKL_INT* Arowptrs,
                                const MKL_INT* Aentries, const double* Avalues,
                                const double* x, double* y) {
  mkl_dbsrmv(&mode, &m, &n, &b, &alpha, "G**C", Avalues, Aentries, Arowptrs,
             Arowptrs + 1, x, &beta, y);
}

inline void spmv_block_impl_mkl(char mode, Kokkos::complex<float> alpha,
                                Kokkos::complex<float> beta, MKL_INT m,
                                MKL_INT n, MKL_INT b, const MKL_INT* Arowptrs,
                                const MKL_INT* Aentries,
                                const Kokkos::complex<float>* Avalues,
                                const Kokkos::complex<float>* x,
                                Kokkos::complex<float>* y) {
  const MKL_Complex8* alpha_mkl = reinterpret_cast<const MKL_Complex8*>(&alpha);
  const MKL_Complex8* beta_mkl  = reinterpret_cast<const MKL_Complex8*>(&beta);
  const MKL_Complex8* Avalues_mkl =
      reinterpret_cast<const MKL_Complex8*>(Avalues);
  const MKL_Complex8* x_mkl = reinterpret_cast<const MKL_Complex8*>(x);
  MKL_Complex8* y_mkl       = reinterpret_cast<MKL_Complex8*>(y);
  mkl_cbsrmv(&mode, &m, &n, &b, alpha_mkl, "G**C", Avalues_mkl, Aentries,
             Arowptrs, Arowptrs + 1, x_mkl, beta_mkl, y_mkl);
}

inline void spmv_block_impl_mkl(char mode, Kokkos::complex<double> alpha,
                                Kokkos::complex<double> beta, MKL_INT m,
                                MKL_INT n, MKL_INT b, const MKL_INT* Arowptrs,
                                const MKL_INT* Aentries,
                                const Kokkos::complex<double>* Avalues,
                                const Kokkos::complex<double>* x,
                                Kokkos::complex<double>* y) {
  const MKL_Complex16* alpha_mkl =
      reinterpret_cast<const MKL_Complex16*>(&alpha);
  const MKL_Complex16* beta_mkl = reinterpret_cast<const MKL_Complex16*>(&beta);
  const MKL_Complex16* Avalues_mkl =
      reinterpret_cast<const MKL_Complex16*>(Avalues);
  const MKL_Complex16* x_mkl = reinterpret_cast<const MKL_Complex16*>(x);
  MKL_Complex16* y_mkl       = reinterpret_cast<MKL_Complex16*>(y);
  mkl_zbsrmv(&mode, &m, &n, &b, alpha_mkl, "G**C", Avalues_mkl, Aentries,
             Arowptrs, Arowptrs + 1, x_mkl, beta_mkl, y_mkl);
}

inline void spm_mv_block_impl_mkl(char mode, float alpha, float beta, MKL_INT m,
                                  MKL_INT n, MKL_INT b, const MKL_INT* Arowptrs,
                                  const MKL_INT* Aentries, const float* Avalues,
                                  const float* x, MKL_INT colx, MKL_INT ldx,
                                  float* y, MKL_INT ldy) {
  mkl_sbsrmm(&mode, &m, &n, &colx, &b, &alpha, "G**C", Avalues, Aentries,
             Arowptrs, Arowptrs + 1, x, &beta, y);
}

inline void spm_mv_block_impl_mkl(
    char mode, double alpha, double beta, MKL_INT m, MKL_INT n, MKL_INT b,
    const MKL_INT* Arowptrs, const MKL_INT* Aentries, const double* Avalues,
    const double* x, MKL_INT colx, MKL_INT ldx, double* y, MKL_INT ldy) {
  mkl_dbsrmm(&mode, &m, &n, &colx, &b, &alpha, "G**C", Avalues, Aentries,
             Arowptrs, Arowptrs + 1, x, ldx, &beta, y, ldy);
}

inline void spm_mv_block_impl_mkl(char mode, Kokkos::complex<float> alpha,
                                  Kokkos::complex<float> beta, MKL_INT m,
                                  MKL_INT n, MKL_INT b, const MKL_INT* Arowptrs,
                                  const MKL_INT* Aentries,
                                  const Kokkos::complex<float>* Avalues,
                                  const Kokkos::complex<float>* x, MKL_INT colx,
                                  MKL_INT ldx, Kokkos::complex<float>* y,
                                  MKL_INT ldy) {
  const MKL_Complex8* alpha_mkl = reinterpret_cast<const MKL_Complex8*>(&alpha);
  const MKL_Complex8* beta_mkl  = reinterpret_cast<const MKL_Complex8*>(&beta);
  const MKL_Complex8* Avalues_mkl =
      reinterpret_cast<const MKL_Complex8*>(Avalues);
  const MKL_Complex8* x_mkl = reinterpret_cast<const MKL_Complex8*>(x);
  MKL_Complex8* y_mkl       = reinterpret_cast<MKL_Complex8*>(y);
  mkl_cbsrmv(&mode, &m, &n, &colx, &b, alpha_mkl, "G**C", Avalues_mkl, Aentries,
             Arowptrs, Arowptrs + 1, x_mkl, ldx, beta_mkl, y_mkl, ldy);
}

inline void spm_mv_block_impl_mkl(char mode, Kokkos::complex<double> alpha,
                                  Kokkos::complex<double> beta, MKL_INT m,
                                  MKL_INT n, MKL_INT b, const MKL_INT* Arowptrs,
                                  const MKL_INT* Aentries,
                                  const Kokkos::complex<double>* Avalues,
                                  const Kokkos::complex<double>* x,
                                  MKL_INT colx, MKL_INT ldx,
                                  Kokkos::complex<double>* y, MKL_INT ldy) {
  const MKL_Complex16* alpha_mkl =
      reinterpret_cast<const MKL_Complex16*>(&alpha);
  const MKL_Complex16* beta_mkl = reinterpret_cast<const MKL_Complex16*>(&beta);
  const MKL_Complex16* Avalues_mkl =
      reinterpret_cast<const MKL_Complex16*>(Avalues);
  const MKL_Complex16* x_mkl = reinterpret_cast<const MKL_Complex16*>(x);
  MKL_Complex16* y_mkl       = reinterpret_cast<MKL_Complex16*>(y);
  mkl_zbsrmv(&mode, &m, &n, &colx, &b, alpha_mkl, "G**C", Avalues_mkl, Aentries,
             Arowptrs, Arowptrs + 1, x_mkl, ldx, beta_mkl, y_mkl, ldy);
}

#endif

#define KOKKOSSPARSE_SPMV_MKL(SCALAR, EXECSPACE, COMPILE_LIBRARY)              \
  template <>                                                                  \
  struct SPMV_BSRMATRIX<                                                       \
      SCALAR const, MKL_INT const,                                             \
      Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,                            \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, MKL_INT const, SCALAR const*,   \
      Kokkos::LayoutLeft, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,        \
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>, SCALAR*, \
      Kokkos::LayoutLeft, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,        \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, true, COMPILE_LIBRARY> {        \
    using device_type = Kokkos::Device<EXECSPACE, Kokkos::HostSpace>;          \
    using AMatrix =                                                            \
        BsrMatrix<SCALAR const, MKL_INT const, device_type,                    \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged>, MKL_INT const>;     \
    using XVector = Kokkos::View<                                              \
        SCALAR const*, Kokkos::LayoutLeft, device_type,                        \
        Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;       \
    using YVector = Kokkos::View<SCALAR*, Kokkos::LayoutLeft, device_type,     \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;     \
    using coefficient_type = typename YVector::non_const_value_type;           \
                                                                               \
    static void spmv_bsrmatrix(                                                \
        const KokkosKernels::Experimental::Controls& /*controls*/,             \
        const char mode[], const coefficient_type& alpha, const AMatrix& A,    \
        const XVector& X, const coefficient_type& beta, const YVector& Y) {    \
      std::string label = "KokkosSparse::spmv[TPL_MKL,BSRMATRIX" +             \
                          Kokkos::ArithTraits<SCALAR>::name() + "]";           \
      Kokkos::Profiling::pushRegion(label);                                    \
      spmv_block_impl_mkl(mode_kk_to_mkl(mode[0]), alpha, beta, A.numRows(),   \
                          A.numCols(), A.blockDim(), A.graph.row_map.data(),   \
                          A.graph.entries.data(), A.values.data(), X.data(),   \
                          Y.data());                                           \
      Kokkos::Profiling::popRegion();                                          \
    }                                                                          \
  };

#ifdef KOKKOS_ENABLE_SERIAL
KOKKOSSPARSE_SPMV_MKL(float, Kokkos::Serial, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MKL(double, Kokkos::Serial,
                      KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MKL(Kokkos::complex<float>, Kokkos::Serial,
                      KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MKL(Kokkos::complex<double>, Kokkos::Serial,
                      KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
#endif

#ifdef KOKKOS_ENABLE_OPENMP
KOKKOSSPARSE_SPMV_MKL(float, Kokkos::OpenMP, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MKL(double, Kokkos::OpenMP,
                      KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MKL(Kokkos::complex<float>, Kokkos::OpenMP,
                      KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MKL(Kokkos::complex<double>, Kokkos::OpenMP,
                      KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
#endif

#undef KOKKOSSPARSE_SPMV_MKL

#define KOKKOSSPARSE_SPMV_MV_MKL(SCALAR, EXECSPACE, COMPILE_LIBRARY)           \
  template <>                                                                  \
  struct SPMV_MV_BSRMATRIX<                                                    \
      SCALAR const, MKL_INT const,                                             \
      Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,                            \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, MKL_INT const, SCALAR const**,  \
      Kokkos::LayoutLeft, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,        \
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>,          \
      SCALAR**, Kokkos::LayoutLeft,                                            \
      Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,                            \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, true, true, COMPILE_LIBRARY> {  \
    using device_type = Kokkos::Device<EXECSPACE, Kokkos::HostSpace>;          \
    using AMatrix =                                                            \
        BsrMatrix<SCALAR const, MKL_INT const, device_type,                    \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged>, MKL_INT const>;     \
    using XVector = Kokkos::View<                                              \
        SCALAR const**, Kokkos::LayoutLeft, device_type,                       \
        Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;       \
    using YVector = Kokkos::View<SCALAR**, Kokkos::LayoutLeft, device_type,    \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;     \
    using coefficient_type = typename YVector::non_const_value_type;           \
                                                                               \
    static void spmv_mv_bsrmatrix(                                             \
        const KokkosKernels::Experimental::Controls& /*controls*/,             \
        const char mode[], const coefficient_type& alpha, const AMatrix& A,    \
        const XVector& X, const coefficient_type& beta, const YVector& Y) {    \
      std::string label = "KokkosSparse::spmv[TPL_MKL,BSRMATRIX" +             \
                          Kokkos::ArithTraits<SCALAR>::name() + "]";           \
      Kokkos::Profiling::pushRegion(label);                                    \
      MKL_INT colx = static_cast<MKL_INT>(X.extent(1));                        \
      MKL_INT ldx  = static_cast<MKL_INT>(X.stride_1());                       \
      MKL_INT ldy  = static_cast<MKL_INT>(Y.stride_1());                       \
      spm_mv_block_impl_mkl(mode_kk_to_mkl(mode[0]), alpha, beta, A.numRows(), \
                            A.numCols(), A.blockDim(), A.graph.row_map.data(), \
                            A.graph.entries.data(), A.values.data(), X.data(), \
                            colx, ldx, Y.data(), ldy);                         \
      Kokkos::Profiling::popRegion();                                          \
    }                                                                          \
  };

#ifdef KOKKOS_ENABLE_SERIAL
KOKKOSSPARSE_SPMV_MV_MKL(float, Kokkos::Serial,
                         KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MV_MKL(double, Kokkos::Serial,
                         KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MV_MKL(Kokkos::complex<float>, Kokkos::Serial,
                         KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MV_MKL(Kokkos::complex<double>, Kokkos::Serial,
                         KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
#endif

#ifdef KOKKOS_ENABLE_OPENMP
KOKKOSSPARSE_SPMV_MV_MKL(float, Kokkos::OpenMP,
                         KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MV_MKL(double, Kokkos::OpenMP,
                         KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MV_MKL(Kokkos::complex<float>, Kokkos::OpenMP,
                         KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MV_MKL(Kokkos::complex<double>, Kokkos::OpenMP,
                         KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
#endif

#undef KOKKOSSPARSE_SPMV_MV_MKL

}  // namespace Impl
}  // namespace Experimental
}  // namespace KokkosSparse

#endif  // KOKKOSKERNELS_ENABLE_TPL_MKL

// cuSPARSE
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#include "cusparse.h"
#include "KokkosSparse_Utils_cusparse.hpp"

//
// From  https://docs.nvidia.com/cuda/cusparse/index.html#bsrmv
// Several comments on bsrmv():
// - Only blockDim > 1 is supported
// - Only CUSPARSE_OPERATION_NON_TRANSPOSE is supported
// - Only CUSPARSE_MATRIX_TYPE_GENERAL is supported.
//
namespace KokkosSparse {
namespace Experimental {
namespace Impl {

template <class AMatrix, class XVector, class YVector>
void spmv_block_impl_cusparse(
    const KokkosKernels::Experimental::Controls& controls, const char mode[],
    typename YVector::non_const_value_type const& alpha, const AMatrix& A,
    const XVector& x, typename YVector::non_const_value_type const& beta,
    const YVector& y) {
  using offset_type = typename AMatrix::non_const_size_type;
  using entry_type  = typename AMatrix::non_const_ordinal_type;
  using value_type  = typename AMatrix::non_const_value_type;

  /* initialize cusparse library */
  cusparseHandle_t cusparseHandle = controls.getCusparseHandle();

  /* Set the operation mode */
  cusparseOperation_t myCusparseOperation;
  switch (toupper(mode[0])) {
    case 'N': myCusparseOperation = CUSPARSE_OPERATION_NON_TRANSPOSE; break;
    default: {
      std::cerr << "Mode " << mode << " invalid for cusparse[*]bsrmv.\n";
      throw std::invalid_argument("Invalid mode");
    }
  }

#if (9000 <= CUDA_VERSION)

  /* create and set the matrix descriptor */
  cusparseMatDescr_t descrA = 0;
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateMatDescr(&descrA));
  KOKKOS_CUSPARSE_SAFE_CALL(
      cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL));
  KOKKOS_CUSPARSE_SAFE_CALL(
      cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO));
  cusparseDirection_t dirA = CUSPARSE_DIRECTION_ROW;

  /* perform the actual SpMV operation */
  if ((std::is_same<int, offset_type>::value) &&
      (std::is_same<int, entry_type>::value)) {
    if (std::is_same<value_type, float>::value) {
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseSbsrmv(
          cusparseHandle, dirA, myCusparseOperation, A.numRows(), A.numCols(),
          A.nnz(), reinterpret_cast<float const*>(&alpha), descrA,
          reinterpret_cast<float const*>(A.values.data()),
          A.graph.row_map.data(), A.graph.entries.data(), A.blockDim(),
          reinterpret_cast<float const*>(x.data()),
          reinterpret_cast<float const*>(&beta),
          reinterpret_cast<float*>(y.data())));
    } else if (std::is_same<value_type, double>::value) {
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseDbsrmv(
          cusparseHandle, dirA, myCusparseOperation, A.numRows(), A.numCols(),
          A.nnz(), reinterpret_cast<double const*>(&alpha), descrA,
          reinterpret_cast<double const*>(A.values.data()),
          A.graph.row_map.data(), A.graph.entries.data(), A.blockDim(),
          reinterpret_cast<double const*>(x.data()),
          reinterpret_cast<double const*>(&beta),
          reinterpret_cast<double*>(y.data())));
    } else if (std::is_same<value_type, Kokkos::complex<float>>::value) {
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseCbsrmv(
          cusparseHandle, dirA, myCusparseOperation, A.numRows(), A.numCols(),
          A.nnz(), reinterpret_cast<cuComplex const*>(&alpha), descrA,
          reinterpret_cast<cuComplex const*>(A.values.data()),
          A.graph.row_map.data(), A.graph.entries.data(), A.blockDim(),
          reinterpret_cast<cuComplex const*>(x.data()),
          reinterpret_cast<cuComplex const*>(&beta),
          reinterpret_cast<cuComplex*>(y.data())));
    } else if (std::is_same<value_type, Kokkos::complex<double>>::value) {
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseZbsrmv(
          cusparseHandle, dirA, myCusparseOperation, A.numRows(), A.numCols(),
          A.nnz(), reinterpret_cast<cuDoubleComplex const*>(&alpha), descrA,
          reinterpret_cast<cuDoubleComplex const*>(A.values.data()),
          A.graph.row_map.data(), A.graph.entries.data(), A.blockDim(),
          reinterpret_cast<cuDoubleComplex const*>(x.data()),
          reinterpret_cast<cuDoubleComplex const*>(&beta),
          reinterpret_cast<cuDoubleComplex*>(y.data())));
    } else {
      throw std::logic_error(
          "Trying to call cusparse[*]bsrmv with a scalar type not "
          "float/double, "
          "nor complex of either!");
    }
  } else {
    throw std::logic_error(
        "With cuSPARSE pre-10.0, offset and entry types must be int. "
        "Something wrong with TPL avail logic.");
  }

  KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroyMatDescr(descrA));
#endif  // (9000 <= CUDA_VERSION)
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
//   so we just do not support LayoutRight in cuSparse TPL now
//
template <
    class AMatrix, class XVector, class YVector,
    std::enable_if_t<std::is_same<Kokkos::LayoutLeft,
                                  typename XVector::array_layout>::value &&
                         std::is_same<Kokkos::LayoutLeft,
                                      typename YVector::array_layout>::value,
                     bool> = true>
void spm_mv_block_impl_cusparse(
    const KokkosKernels::Experimental::Controls& controls, const char mode[],
    typename YVector::non_const_value_type const& alpha, const AMatrix& A,
    const XVector& x, typename YVector::non_const_value_type const& beta,
    const YVector& y) {
  using offset_type = typename AMatrix::non_const_size_type;
  using entry_type  = typename AMatrix::non_const_ordinal_type;
  using value_type  = typename AMatrix::non_const_value_type;

  /* initialize cusparse library */
  cusparseHandle_t cusparseHandle = controls.getCusparseHandle();

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

  // ldx and ldy should be the leading dimension of X,Y respectively
  const int ldx = static_cast<int>(x.extent(0));
  const int ldy = static_cast<int>(y.extent(0));

#if (9000 <= CUDA_VERSION)

  /* create and set the matrix descriptor */
  cusparseMatDescr_t descrA = 0;
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateMatDescr(&descrA));
  KOKKOS_CUSPARSE_SAFE_CALL(
      cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL));
  KOKKOS_CUSPARSE_SAFE_CALL(
      cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO));
  cusparseDirection_t dirA = CUSPARSE_DIRECTION_ROW;

  /* perform the actual SpMV operation */
  if ((std::is_same<int, offset_type>::value) &&
      (std::is_same<int, entry_type>::value)) {
    if (std::is_same<value_type, float>::value) {
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseSbsrmm(
          cusparseHandle, dirA, myCusparseOperation,
          CUSPARSE_OPERATION_NON_TRANSPOSE, A.numRows(), colx, A.numCols(),
          A.nnz(), reinterpret_cast<float const*>(&alpha), descrA,
          reinterpret_cast<float const*>(A.values.data()),
          A.graph.row_map.data(), A.graph.entries.data(), A.blockDim(),
          reinterpret_cast<float const*>(x.data()), ldx,
          reinterpret_cast<float const*>(&beta),
          reinterpret_cast<float*>(y.data()), ldy));
    } else if (std::is_same<value_type, double>::value) {
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseDbsrmm(
          cusparseHandle, dirA, myCusparseOperation,
          CUSPARSE_OPERATION_NON_TRANSPOSE, A.numRows(), colx, A.numCols(),
          A.nnz(), reinterpret_cast<double const*>(&alpha), descrA,
          reinterpret_cast<double const*>(A.values.data()),
          A.graph.row_map.data(), A.graph.entries.data(), A.blockDim(),
          reinterpret_cast<double const*>(x.data()), ldx,
          reinterpret_cast<double const*>(&beta),
          reinterpret_cast<double*>(y.data()), ldy));
    } else if (std::is_same<value_type, Kokkos::complex<float>>::value) {
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseCbsrmm(
          cusparseHandle, dirA, myCusparseOperation,
          CUSPARSE_OPERATION_NON_TRANSPOSE, A.numRows(), colx, A.numCols(),
          A.nnz(), reinterpret_cast<cuComplex const*>(&alpha), descrA,
          reinterpret_cast<cuComplex const*>(A.values.data()),
          A.graph.row_map.data(), A.graph.entries.data(), A.blockDim(),
          reinterpret_cast<cuComplex const*>(x.data()), ldx,
          reinterpret_cast<cuComplex const*>(&beta),
          reinterpret_cast<cuComplex*>(y.data()), ldy));
    } else if (std::is_same<value_type, Kokkos::complex<double>>::value) {
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseZbsrmm(
          cusparseHandle, dirA, myCusparseOperation,
          CUSPARSE_OPERATION_NON_TRANSPOSE, A.numRows(), colx, A.numCols(),
          A.nnz(), reinterpret_cast<cuDoubleComplex const*>(&alpha), descrA,
          reinterpret_cast<cuDoubleComplex const*>(A.values.data()),
          A.graph.row_map.data(), A.graph.entries.data(), A.blockDim(),
          reinterpret_cast<cuDoubleComplex const*>(x.data()), ldx,
          reinterpret_cast<cuDoubleComplex const*>(&beta),
          reinterpret_cast<cuDoubleComplex*>(y.data()), ldy));
    } else {
      throw std::logic_error(
          "Trying to call cusparse[*]bsrmm with a scalar type not "
          "float/double, "
          "nor complex of either!");
    }
  } else {
    throw std::logic_error(
        "With cuSPARSE pre-10.0, offset and entry types must be int. "
        "Something wrong with TPL avail logic.");
  }

  KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroyMatDescr(descrA));
#endif  // (9000 <= CUDA_VERSION)
}

#define KOKKOSSPARSE_SPMV_CUSPARSE(SCALAR, ORDINAL, OFFSET, LAYOUT, SPACE,     \
                                   COMPILE_LIBRARY)                            \
  template <>                                                                  \
  struct SPMV_BSRMATRIX<                                                       \
      SCALAR const, ORDINAL const, Kokkos::Device<Kokkos::Cuda, SPACE>,        \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, OFFSET const, SCALAR const*,    \
      LAYOUT, Kokkos::Device<Kokkos::Cuda, SPACE>,                             \
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>, SCALAR*, \
      LAYOUT, Kokkos::Device<Kokkos::Cuda, SPACE>,                             \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, true, COMPILE_LIBRARY> {        \
    using device_type       = Kokkos::Device<Kokkos::Cuda, SPACE>;             \
    using memory_trait_type = Kokkos::MemoryTraits<Kokkos::Unmanaged>;         \
    using AMatrix = BsrMatrix<SCALAR const, ORDINAL const, device_type,        \
                              memory_trait_type, OFFSET const>;                \
    using XVector = Kokkos::View<                                              \
        SCALAR const*, LAYOUT, device_type,                                    \
        Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;       \
    using YVector =                                                            \
        Kokkos::View<SCALAR*, LAYOUT, device_type, memory_trait_type>;         \
    using Controls = KokkosKernels::Experimental::Controls;                    \
                                                                               \
    using coefficient_type = typename YVector::non_const_value_type;           \
                                                                               \
    static void spmv_bsrmatrix(const Controls& controls, const char mode[],    \
                               const coefficient_type& alpha,                  \
                               const AMatrix& A, const XVector& x,             \
                               const coefficient_type& beta,                   \
                               const YVector& y) {                             \
      std::string label = "KokkosSparse::spmv[TPL_CUSPARSE,BSRMATRIX" +        \
                          Kokkos::ArithTraits<SCALAR>::name() + "]";           \
      Kokkos::Profiling::pushRegion(label);                                    \
      spmv_block_impl_cusparse(controls, mode, alpha, A, x, beta, y);          \
      Kokkos::Profiling::popRegion();                                          \
    }                                                                          \
  };

#if (9000 <= CUDA_VERSION)
KOKKOSSPARSE_SPMV_CUSPARSE(double, int, int, Kokkos::LayoutLeft,
                           Kokkos::CudaSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(double, int, int, Kokkos::LayoutRight,
                           Kokkos::CudaSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(float, int, int, Kokkos::LayoutLeft,
                           Kokkos::CudaSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(float, int, int, Kokkos::LayoutRight,
                           Kokkos::CudaSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<double>, int, int,
                           Kokkos::LayoutLeft, Kokkos::CudaSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<double>, int, int,
                           Kokkos::LayoutRight, Kokkos::CudaSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft,
                           Kokkos::CudaSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<float>, int, int,
                           Kokkos::LayoutRight, Kokkos::CudaSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(double, int, int, Kokkos::LayoutLeft,
                           Kokkos::CudaUVMSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(double, int, int, Kokkos::LayoutRight,
                           Kokkos::CudaUVMSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(float, int, int, Kokkos::LayoutLeft,
                           Kokkos::CudaUVMSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(float, int, int, Kokkos::LayoutRight,
                           Kokkos::CudaUVMSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<double>, int, int,
                           Kokkos::LayoutLeft, Kokkos::CudaUVMSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<double>, int, int,
                           Kokkos::LayoutRight, Kokkos::CudaUVMSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft,
                           Kokkos::CudaUVMSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<float>, int, int,
                           Kokkos::LayoutRight, Kokkos::CudaUVMSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
#endif  // (9000 <= CUDA_VERSION)

#undef KOKKOSSPARSE_SPMV_CUSPARSE

// cuSparse TPL does not support LayoutRight for this operation
// only specialize for LayoutLeft
#define KOKKOSSPARSE_SPMV_MV_CUSPARSE(SCALAR, ORDINAL, OFFSET, SPACE,          \
                                      ETI_AVAIL)                               \
  template <>                                                                  \
  struct SPMV_MV_BSRMATRIX<                                                    \
      SCALAR const, ORDINAL const, Kokkos::Device<Kokkos::Cuda, SPACE>,        \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, OFFSET const, SCALAR const**,   \
      Kokkos::LayoutLeft, Kokkos::Device<Kokkos::Cuda, SPACE>,                 \
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>,          \
      SCALAR**, Kokkos::LayoutLeft, Kokkos::Device<Kokkos::Cuda, SPACE>,       \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, false, true, ETI_AVAIL> {       \
    using device_type       = Kokkos::Device<Kokkos::Cuda, SPACE>;             \
    using memory_trait_type = Kokkos::MemoryTraits<Kokkos::Unmanaged>;         \
    using AMatrix = BsrMatrix<SCALAR const, ORDINAL const, device_type,        \
                              memory_trait_type, OFFSET const>;                \
    using XVector = Kokkos::View<                                              \
        SCALAR const**, Kokkos::LayoutLeft, device_type,                       \
        Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;       \
    using YVector  = Kokkos::View<SCALAR**, Kokkos::LayoutLeft, device_type,   \
                                 memory_trait_type>;                          \
    using Controls = KokkosKernels::Experimental::Controls;                    \
                                                                               \
    using coefficient_type = typename YVector::non_const_value_type;           \
                                                                               \
    static void spmv_mv_bsrmatrix(const Controls& controls, const char mode[], \
                                  const coefficient_type& alpha,               \
                                  const AMatrix& A, const XVector& x,          \
                                  const coefficient_type& beta,                \
                                  const YVector& y) {                          \
      std::string label = "KokkosSparse::spmv[TPL_CUSPARSE,BSRMATRIX" +        \
                          Kokkos::ArithTraits<SCALAR>::name() + "]";           \
      Kokkos::Profiling::pushRegion(label);                                    \
      spm_mv_block_impl_cusparse(controls, mode, alpha, A, x, beta, y);        \
      Kokkos::Profiling::popRegion();                                          \
    }                                                                          \
  };

#if (9000 <= CUDA_VERSION)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(double, int, int, Kokkos::CudaSpace, true)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(double, int, int, Kokkos::CudaSpace, false)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(float, int, int, Kokkos::CudaSpace, true)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(float, int, int, Kokkos::CudaSpace, false)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<double>, int, int,
                              Kokkos::CudaSpace, true)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<double>, int, int,
                              Kokkos::CudaSpace, false)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<float>, int, int,
                              Kokkos::CudaSpace, true)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<float>, int, int,
                              Kokkos::CudaSpace, false)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(double, int, int, Kokkos::CudaUVMSpace, true)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(double, int, int, Kokkos::CudaUVMSpace, false)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(float, int, int, Kokkos::CudaUVMSpace, true)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(float, int, int, Kokkos::CudaUVMSpace, false)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<double>, int, int,
                              Kokkos::CudaUVMSpace, true)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<double>, int, int,
                              Kokkos::CudaUVMSpace, false)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<float>, int, int,
                              Kokkos::CudaUVMSpace, true)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<float>, int, int,
                              Kokkos::CudaUVMSpace, false)

#endif  // (9000 <= CUDA_VERSION)

#undef KOKKOSSPARSE_SPMV_MV_CUSPARSE

}  // namespace Impl
}  // namespace Experimental
}  // namespace KokkosSparse

#endif  // KOKKOSKERNELS_ENABLE_TPL_CUSPARSE

// --------------------
// rocSparse
// --------------------
#if defined(KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE)

#include <rocsparse/rocsparse.h>

#include "KokkosSparse_Utils_rocsparse.hpp"

namespace KokkosSparse {
namespace Experimental {
namespace Impl {

template <class AMatrix, class XVector, class YVector>
void spmv_block_impl_rocsparse(
    const KokkosKernels::Experimental::Controls& controls, const char mode[],
    typename YVector::non_const_value_type const& alpha, const AMatrix& A,
    const XVector& x, typename YVector::non_const_value_type const& beta,
    const YVector& y) {
  /*
     rocm 5.4.0 rocsparse_*bsrmv reference:
     https://rocsparse.readthedocs.io/en/rocm-5.4.0/usermanual.html#rocsparse-bsrmv-ex

     only trans = rocsparse_operation_none is supported
     only descr = rocsparse_matrix_type_general is supported

  */

  using offset_type  = typename AMatrix::non_const_size_type;
  using ordinal_type = typename AMatrix::non_const_ordinal_type;
  using value_type   = typename AMatrix::non_const_value_type;
  using rocsparse_value_type =
      typename KokkosSparse::Impl::kokkos_to_rocsparse_type<value_type>::type;

  // assert ordinals and offsets are the expected types
  static_assert(std::is_same_v<offset_type, rocsparse_int>,
                "A offset_type must be rocsparse_int");
  static_assert(std::is_same_v<ordinal_type, rocsparse_int>,
                "A ordinal_type must be rocsparse_int");

  // assert all operands are the same type
  using x_value_type = typename XVector::non_const_value_type;
  using y_value_type = typename YVector::non_const_value_type;
  static_assert(std::is_same_v<value_type, x_value_type>,
                "A and x must have same value type");
  static_assert(std::is_same_v<value_type, y_value_type>,
                "A and y must have same value type");

  // assert X and Y are non-stride (pass raw pointers to TPL)
  static_assert(
      !std::is_same_v<typename XVector::array_layout, Kokkos::LayoutStride>,
      "x must be contiguous");
  static_assert(
      !std::is_same_v<typename YVector::array_layout, Kokkos::LayoutStride>,
      "y must be contiguous");

  // assert BSR data is non-stride (pass raw pointers to TPL)
  static_assert(!std::is_same_v<typename AMatrix::values_type::array_layout,
                                Kokkos::LayoutStride>,
                "A values must be contiguous");
  static_assert(!std::is_same_v<typename AMatrix::row_map_type::array_layout,
                                Kokkos::LayoutStride>,
                "A row_map must be contiguous");
  static_assert(!std::is_same_v<typename AMatrix::index_type::array_layout,
                                Kokkos::LayoutStride>,
                "A entries must be contiguous");

  rocsparse_handle handle = controls.getRocsparseHandle();

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
  static_assert(
      std::is_same_v<typename AMatrix::block_layout_type, Kokkos::LayoutRight>,
      "A blocks must be stored layout-right");
  rocsparse_direction dir = rocsparse_direction_row;

  const rocsparse_int mb = rocsparse_int(A.numRows());  // number of block rows
  const rocsparse_int nb = rocsparse_int(A.numCols());  // number of block cols
  const rocsparse_int nnzb =
      rocsparse_int(A.nnz());  // number of non-zero blocks
  const rocsparse_value_type* alpha_ =
      reinterpret_cast<const rocsparse_value_type*>(&alpha);

  const rocsparse_value_type* bsr_val =
      reinterpret_cast<const rocsparse_value_type*>(A.values.data());
  const rocsparse_int* bsr_row_ptr = A.graph.row_map.data();
  const rocsparse_int* bsr_col_ind = A.graph.entries.data();
  const rocsparse_int block_dim    = rocsparse_int(A.blockDim());
  const rocsparse_value_type* x_ =
      reinterpret_cast<const rocsparse_value_type*>(x.data());
  const rocsparse_value_type* beta_ =
      reinterpret_cast<const rocsparse_value_type*>(&beta);
  rocsparse_value_type* y_ = reinterpret_cast<rocsparse_value_type*>(y.data());

  rocsparse_mat_descr descr;
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_create_mat_descr(&descr));
  rocsparse_mat_info info;
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_create_mat_info(&info));

  // *_ex* functions introduced in 5.4.0
#if KOKKOSSPARSE_IMPL_ROCM_VERSION < 50400
  if constexpr (std::is_same_v<value_type, float>) {
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_sbsrmv(
        handle, dir, trans, mb, nb, nnzb, alpha_, descr, bsr_val, bsr_row_ptr,
        bsr_col_ind, block_dim, x_, beta_, y_));
  } else if constexpr (std::is_same_v<value_type, double>) {
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_dbsrmv(
        handle, dir, trans, mb, nb, nnzb, alpha_, descr, bsr_val, bsr_row_ptr,
        bsr_col_ind, block_dim, x_, beta_, y_));
  } else if constexpr (std::is_same_v<value_type, Kokkos::complex<float>>) {
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_cbsrmv(
        handle, dir, trans, mb, nb, nnzb, alpha_, descr, bsr_val, bsr_row_ptr,
        bsr_col_ind, block_dim, x_, beta_, y_));
  } else if constexpr (std::is_same_v<value_type, Kokkos::complex<double>>) {
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_zbsrmv(
        handle, dir, trans, mb, nb, nnzb, alpha_, descr, bsr_val, bsr_row_ptr,
        bsr_col_ind, block_dim, x_, beta_, y_));
  } else {
    static_assert(KokkosKernels::Impl::always_false_v<value_type>,
                  "unsupported value type for rocsparse_*bsrmv");
  }
#else
  if constexpr (std::is_same_v<value_type, float>) {
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_sbsrmv_ex_analysis(
        handle, dir, trans, mb, nb, nnzb, descr, bsr_val, bsr_row_ptr,
        bsr_col_ind, block_dim, info));
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_sbsrmv_ex(
        handle, dir, trans, mb, nb, nnzb, alpha_, descr, bsr_val, bsr_row_ptr,
        bsr_col_ind, block_dim, info, x_, beta_, y_));
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_bsrsv_clear(handle, info));
  } else if constexpr (std::is_same_v<value_type, double>) {
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_dbsrmv_ex_analysis(
        handle, dir, trans, mb, nb, nnzb, descr, bsr_val, bsr_row_ptr,
        bsr_col_ind, block_dim, info));
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_dbsrmv_ex(
        handle, dir, trans, mb, nb, nnzb, alpha_, descr, bsr_val, bsr_row_ptr,
        bsr_col_ind, block_dim, info, x_, beta_, y_));
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_bsrsv_clear(handle, info));
  } else if constexpr (std::is_same_v<value_type, Kokkos::complex<float>>) {
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_cbsrmv_ex_analysis(
        handle, dir, trans, mb, nb, nnzb, descr, bsr_val, bsr_row_ptr,
        bsr_col_ind, block_dim, info));
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_cbsrmv_ex(
        handle, dir, trans, mb, nb, nnzb, alpha_, descr, bsr_val, bsr_row_ptr,
        bsr_col_ind, block_dim, info, x_, beta_, y_));
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_bsrsv_clear(handle, info));
  } else if constexpr (std::is_same_v<value_type, Kokkos::complex<double>>) {
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_zbsrmv_ex_analysis(
        handle, dir, trans, mb, nb, nnzb, descr, bsr_val, bsr_row_ptr,
        bsr_col_ind, block_dim, info));
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_zbsrmv_ex(
        handle, dir, trans, mb, nb, nnzb, alpha_, descr, bsr_val, bsr_row_ptr,
        bsr_col_ind, block_dim, info, x_, beta_, y_));
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_bsrsv_clear(handle, info));
  } else {
    static_assert(KokkosKernels::Impl::always_false_v<value_type>,
                  "unsupported value type for rocsparse_*bsrmv");
  }
#endif
  rocsparse_destroy_mat_descr(descr);
  rocsparse_destroy_mat_info(info);

}  // spmv_block_impl_rocsparse

#define KOKKOSSPARSE_SPMV_ROCSPARSE(SCALAR, ORDINAL, OFFSET, LAYOUT, SPACE,    \
                                    COMPILE_LIBRARY)                           \
  template <>                                                                  \
  struct SPMV_BSRMATRIX<                                                       \
      SCALAR const, ORDINAL const, Kokkos::Device<Kokkos::HIP, SPACE>,         \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, OFFSET const, SCALAR const*,    \
      LAYOUT, Kokkos::Device<Kokkos::HIP, SPACE>,                              \
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>, SCALAR*, \
      LAYOUT, Kokkos::Device<Kokkos::HIP, SPACE>,                              \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, true, COMPILE_LIBRARY> {        \
    using device_type       = Kokkos::Device<Kokkos::HIP, SPACE>;              \
    using memory_trait_type = Kokkos::MemoryTraits<Kokkos::Unmanaged>;         \
    using AMatrix = BsrMatrix<SCALAR const, ORDINAL const, device_type,        \
                              memory_trait_type, OFFSET const>;                \
    using XVector = Kokkos::View<                                              \
        SCALAR const*, LAYOUT, device_type,                                    \
        Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;       \
    using YVector =                                                            \
        Kokkos::View<SCALAR*, LAYOUT, device_type, memory_trait_type>;         \
    using Controls = KokkosKernels::Experimental::Controls;                    \
                                                                               \
    using coefficient_type = typename YVector::non_const_value_type;           \
                                                                               \
    static void spmv_bsrmatrix(const Controls& controls, const char mode[],    \
                               const coefficient_type& alpha,                  \
                               const AMatrix& A, const XVector& x,             \
                               const coefficient_type& beta,                   \
                               const YVector& y) {                             \
      std::string label = "KokkosSparse::spmv[TPL_ROCSPARSE,BSRMATRIX" +       \
                          Kokkos::ArithTraits<SCALAR>::name() + "]";           \
      Kokkos::Profiling::pushRegion(label);                                    \
      spmv_block_impl_rocsparse(controls, mode, alpha, A, x, beta, y);         \
      Kokkos::Profiling::popRegion();                                          \
    }                                                                          \
  };

KOKKOSSPARSE_SPMV_ROCSPARSE(float, rocsparse_int, rocsparse_int,
                            Kokkos::LayoutLeft, Kokkos::HIPSpace,
                            KOKKOSKERNELS_IMPL_COMPILE_LIBRARY);
KOKKOSSPARSE_SPMV_ROCSPARSE(float, rocsparse_int, rocsparse_int,
                            Kokkos::LayoutRight, Kokkos::HIPSpace,
                            KOKKOSKERNELS_IMPL_COMPILE_LIBRARY);
KOKKOSSPARSE_SPMV_ROCSPARSE(double, rocsparse_int, rocsparse_int,
                            Kokkos::LayoutLeft, Kokkos::HIPSpace,
                            KOKKOSKERNELS_IMPL_COMPILE_LIBRARY);
KOKKOSSPARSE_SPMV_ROCSPARSE(double, rocsparse_int, rocsparse_int,
                            Kokkos::LayoutRight, Kokkos::HIPSpace,
                            KOKKOSKERNELS_IMPL_COMPILE_LIBRARY);
KOKKOSSPARSE_SPMV_ROCSPARSE(Kokkos::complex<float>, rocsparse_int,
                            rocsparse_int, Kokkos::LayoutLeft, Kokkos::HIPSpace,
                            KOKKOSKERNELS_IMPL_COMPILE_LIBRARY);
KOKKOSSPARSE_SPMV_ROCSPARSE(Kokkos::complex<float>, rocsparse_int,
                            rocsparse_int, Kokkos::LayoutRight,
                            Kokkos::HIPSpace,
                            KOKKOSKERNELS_IMPL_COMPILE_LIBRARY);
KOKKOSSPARSE_SPMV_ROCSPARSE(Kokkos::complex<double>, rocsparse_int,
                            rocsparse_int, Kokkos::LayoutLeft, Kokkos::HIPSpace,
                            KOKKOSKERNELS_IMPL_COMPILE_LIBRARY);
KOKKOSSPARSE_SPMV_ROCSPARSE(Kokkos::complex<double>, rocsparse_int,
                            rocsparse_int, Kokkos::LayoutRight,
                            Kokkos::HIPSpace,
                            KOKKOSKERNELS_IMPL_COMPILE_LIBRARY);

#undef KOKKOSSPARSE_SPMV_ROCSPARSE

}  // namespace Impl
}  // namespace Experimental
}  // namespace KokkosSparse

#endif  // defined(KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE)

#endif  // KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_DECL_HPP
