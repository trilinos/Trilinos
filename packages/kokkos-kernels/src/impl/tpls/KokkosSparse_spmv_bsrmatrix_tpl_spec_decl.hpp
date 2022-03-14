/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSKERNELS_KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_DECL_HPP
#define KOKKOSKERNELS_KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_DECL_HPP

#include "KokkosKernels_Controls.hpp"

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
#include <mkl.h>

namespace KokkosSparse {
namespace Experimental {
namespace Impl {

#if (__INTEL_MKL__ > 2017)
// MKL 2018 and above: use new interface: sparse_matrix_t and mkl_sparse_?_mv()

namespace BSR {
inline void mkl_safe_call(int errcode) {
  if (errcode != SPARSE_STATUS_SUCCESS)
    throw std::runtime_error("MKL returned non-success error code");
}

inline sparse_operation_t mode_kk_to_mkl(char mode_kk) {
  switch (toupper(mode_kk)) {
    case 'N': return SPARSE_OPERATION_NON_TRANSPOSE;
    case 'T': return SPARSE_OPERATION_TRANSPOSE;
    case 'H': return SPARSE_OPERATION_CONJUGATE_TRANSPOSE;
    default:;
  }
  throw std::invalid_argument(
      "Invalid mode for MKL (should be one of N, T, H)");
}
}  // namespace BSR

using BSR::mkl_safe_call;
using BSR::mode_kk_to_mkl;

inline matrix_descr getDescription() {
  matrix_descr A_descr;
  A_descr.type = SPARSE_MATRIX_TYPE_GENERAL;
  A_descr.mode = SPARSE_FILL_MODE_FULL;
  A_descr.diag = SPARSE_DIAG_NON_UNIT;
  return A_descr;
}

inline void spmv_block_impl_mkl(sparse_operation_t op, float alpha, float beta,
                                int m, int n, int b, const int* Arowptrs,
                                const int* Aentries, const float* Avalues,
                                const float* x, float* y) {
  sparse_matrix_t A_mkl;
  mkl_safe_call(mkl_sparse_s_create_bsr(
      &A_mkl, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, m, n, b,
      const_cast<int*>(Arowptrs), const_cast<int*>(Arowptrs + 1),
      const_cast<int*>(Aentries), const_cast<float*>(Avalues)));

  matrix_descr A_descr = getDescription();
  mkl_safe_call(mkl_sparse_s_mv(op, alpha, A_mkl, A_descr, x, beta, y));
}

inline void spmv_block_impl_mkl(sparse_operation_t op, double alpha,
                                double beta, int m, int n, int b,
                                const int* Arowptrs, const int* Aentries,
                                const double* Avalues, const double* x,
                                double* y) {
  sparse_matrix_t A_mkl;
  mkl_safe_call(mkl_sparse_d_create_bsr(
      &A_mkl, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, m, n, b,
      const_cast<int*>(Arowptrs), const_cast<int*>(Arowptrs + 1),
      const_cast<int*>(Aentries), const_cast<double*>(Avalues)));

  matrix_descr A_descr = getDescription();
  mkl_safe_call(mkl_sparse_d_mv(op, alpha, A_mkl, A_descr, x, beta, y));
}

inline void spmv_block_impl_mkl(sparse_operation_t op,
                                Kokkos::complex<float> alpha,
                                Kokkos::complex<float> beta, int m, int n,
                                int b, const int* Arowptrs, const int* Aentries,
                                const Kokkos::complex<float>* Avalues,
                                const Kokkos::complex<float>* x,
                                Kokkos::complex<float>* y) {
  sparse_matrix_t A_mkl;
  mkl_safe_call(mkl_sparse_c_create_bsr(
      &A_mkl, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, m, n, b,
      const_cast<int*>(Arowptrs), const_cast<int*>(Arowptrs + 1),
      const_cast<int*>(Aentries), (MKL_Complex8*)Avalues));

  MKL_Complex8& alpha_mkl = reinterpret_cast<MKL_Complex8&>(alpha);
  MKL_Complex8& beta_mkl  = reinterpret_cast<MKL_Complex8&>(beta);
  matrix_descr A_descr    = getDescription();
  mkl_safe_call(mkl_sparse_c_mv(op, alpha_mkl, A_mkl, A_descr,
                                reinterpret_cast<const MKL_Complex8*>(x),
                                beta_mkl, reinterpret_cast<MKL_Complex8*>(y)));
}

inline void spmv_block_impl_mkl(sparse_operation_t op,
                                Kokkos::complex<double> alpha,
                                Kokkos::complex<double> beta, int m, int n,
                                int b, const int* Arowptrs, const int* Aentries,
                                const Kokkos::complex<double>* Avalues,
                                const Kokkos::complex<double>* x,
                                Kokkos::complex<double>* y) {
  sparse_matrix_t A_mkl;
  mkl_safe_call(mkl_sparse_z_create_bsr(
      &A_mkl, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, m, n, b,
      const_cast<int*>(Arowptrs), const_cast<int*>(Arowptrs + 1),
      const_cast<int*>(Aentries), (MKL_Complex16*)Avalues));

  matrix_descr A_descr     = getDescription();
  MKL_Complex16& alpha_mkl = reinterpret_cast<MKL_Complex16&>(alpha);
  MKL_Complex16& beta_mkl  = reinterpret_cast<MKL_Complex16&>(beta);
  mkl_safe_call(mkl_sparse_z_mv(op, alpha_mkl, A_mkl, A_descr,
                                reinterpret_cast<const MKL_Complex16*>(x),
                                beta_mkl, reinterpret_cast<MKL_Complex16*>(y)));
}

inline void spm_mv_block_impl_mkl(sparse_operation_t op, float alpha,
                                  float beta, int m, int n, int b,
                                  const int* Arowptrs, const int* Aentries,
                                  const float* Avalues, const float* x,
                                  int colx, int ldx, float* y, int ldy) {
  sparse_matrix_t A_mkl;
  mkl_safe_call(mkl_sparse_s_create_bsr(
      &A_mkl, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, m, n, b,
      const_cast<int*>(Arowptrs), const_cast<int*>(Arowptrs + 1),
      const_cast<int*>(Aentries), const_cast<float*>(Avalues)));

  matrix_descr A_descr = getDescription();
  mkl_safe_call(mkl_sparse_s_mm(op, alpha, A_mkl, A_descr,
                                SPARSE_LAYOUT_ROW_MAJOR, x, colx, ldx, beta, y,
                                ldy));
}

inline void spm_mv_block_impl_mkl(sparse_operation_t op, double alpha,
                                  double beta, int m, int n, int b,
                                  const int* Arowptrs, const int* Aentries,
                                  const double* Avalues, const double* x,
                                  int colx, int ldx, double* y, int ldy) {
  sparse_matrix_t A_mkl;
  mkl_safe_call(mkl_sparse_d_create_bsr(
      &A_mkl, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, m, n, b,
      const_cast<int*>(Arowptrs), const_cast<int*>(Arowptrs + 1),
      const_cast<int*>(Aentries), const_cast<double*>(Avalues)));

  matrix_descr A_descr = getDescription();
  mkl_safe_call(mkl_sparse_d_mm(op, alpha, A_mkl, A_descr,
                                SPARSE_LAYOUT_ROW_MAJOR, x, colx, ldx, beta, y,
                                ldy));
}

inline void spm_mv_block_impl_mkl(sparse_operation_t op,
                                  Kokkos::complex<float> alpha,
                                  Kokkos::complex<float> beta, int m, int n,
                                  int b, const int* Arowptrs,
                                  const int* Aentries,
                                  const Kokkos::complex<float>* Avalues,
                                  const Kokkos::complex<float>* x, int colx,
                                  int ldx, Kokkos::complex<float>* y, int ldy) {
  sparse_matrix_t A_mkl;
  mkl_safe_call(mkl_sparse_c_create_bsr(
      &A_mkl, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, m, n, b,
      const_cast<int*>(Arowptrs), const_cast<int*>(Arowptrs + 1),
      const_cast<int*>(Aentries), (MKL_Complex8*)Avalues));

  MKL_Complex8& alpha_mkl = reinterpret_cast<MKL_Complex8&>(alpha);
  MKL_Complex8& beta_mkl  = reinterpret_cast<MKL_Complex8&>(beta);
  matrix_descr A_descr    = getDescription();
  mkl_safe_call(
      mkl_sparse_c_mm(op, alpha_mkl, A_mkl, A_descr, SPARSE_LAYOUT_ROW_MAJOR,
                      reinterpret_cast<const MKL_Complex8*>(x), colx, ldx,
                      beta_mkl, reinterpret_cast<MKL_Complex8*>(y), ldy));
}

inline void spm_mv_block_impl_mkl(
    sparse_operation_t op, Kokkos::complex<double> alpha,
    Kokkos::complex<double> beta, int m, int n, int b, const int* Arowptrs,
    const int* Aentries, const Kokkos::complex<double>* Avalues,
    const Kokkos::complex<double>* x, int colx, int ldx,
    Kokkos::complex<double>* y, int ldy) {
  sparse_matrix_t A_mkl;
  mkl_safe_call(mkl_sparse_z_create_bsr(
      &A_mkl, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, m, n, b,
      const_cast<int*>(Arowptrs), const_cast<int*>(Arowptrs + 1),
      const_cast<int*>(Aentries), (MKL_Complex16*)Avalues));

  matrix_descr A_descr     = getDescription();
  MKL_Complex16& alpha_mkl = reinterpret_cast<MKL_Complex16&>(alpha);
  MKL_Complex16& beta_mkl  = reinterpret_cast<MKL_Complex16&>(beta);
  mkl_safe_call(
      mkl_sparse_z_mm(op, alpha_mkl, A_mkl, A_descr, SPARSE_LAYOUT_ROW_MAJOR,
                      reinterpret_cast<const MKL_Complex16*>(x), colx, ldx,
                      beta_mkl, reinterpret_cast<MKL_Complex16*>(y), ldy));
}

#endif

#if (__INTEL_MKL__ == 2017)

inline void spmv_block_impl_mkl(char mode, float alpha, float beta, int m,
                                int n, int b, const int* Arowptrs,
                                const int* Aentries, const float* Avalues,
                                const float* x, float* y) {
  mkl_sbsrmv(&mode, &m, &n, &b, &alpha, "G**C", Avalues, Aentries, Arowptrs,
             Arowptrs + 1, x, &beta, y);
}

inline void spmv_block_impl_mkl(char mode, double alpha, double beta, int m,
                                int n, int b, const int* Arowptrs,
                                const int* Aentries, const double* Avalues,
                                const double* x, double* y) {
  mkl_dbsrmv(&mode, &m, &n, &b, &alpha, "G**C", Avalues, Aentries, Arowptrs,
             Arowptrs + 1, x, &beta, y);
}

inline void spmv_block_impl_mkl(char mode, Kokkos::complex<float> alpha,
                                Kokkos::complex<float> beta, int m, int n,
                                int b, const int* Arowptrs, const int* Aentries,
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
                                Kokkos::complex<double> beta, int m, int n,
                                int b, const int* Arowptrs, const int* Aentries,
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

inline void spm_mv_block_impl_mkl(char mode, float alpha, float beta, int m,
                                  int n, int b, const int* Arowptrs,
                                  const int* Aentries, const float* Avalues,
                                  const float* x, int colx, int ldx, float* y,
                                  int ldy) {
  mkl_sbsrmm(&mode, &m, &n, &colx, &b, &alpha, "G**C", Avalues, Aentries,
             Arowptrs, Arowptrs + 1, x, &beta, y);
}

inline void spm_mv_block_impl_mkl(char mode, double alpha, double beta, int m,
                                  int n, int b, const int* Arowptrs,
                                  const int* Aentries, const double* Avalues,
                                  const double* x, int colx, int ldx, double* y,
                                  int ldy) {
  mkl_dbsrmm(&mode, &m, &n, &colx, &b, &alpha, "G**C", Avalues, Aentries,
             Arowptrs, Arowptrs + 1, x, ldx, &beta, y, ldy);
}

inline void spm_mv_block_impl_mkl(char mode, Kokkos::complex<float> alpha,
                                  Kokkos::complex<float> beta, int m, int n,
                                  int b, const int* Arowptrs,
                                  const int* Aentries,
                                  const Kokkos::complex<float>* Avalues,
                                  const Kokkos::complex<float>* x, int colx,
                                  int ldx, Kokkos::complex<float>* y, int ldy) {
  const MKL_Complex8* alpha_mkl = reinterpret_cast<const MKL_Complex8*>(&alpha);
  const MKL_Complex8* beta_mkl  = reinterpret_cast<const MKL_Complex8*>(&beta);
  const MKL_Complex8* Avalues_mkl =
      reinterpret_cast<const MKL_Complex8*>(Avalues);
  const MKL_Complex8* x_mkl = reinterpret_cast<const MKL_Complex8*>(x);
  MKL_Complex8* y_mkl       = reinterpret_cast<MKL_Complex8*>(y);
  mkl_cbsrmv(&mode, &m, &n, &colx, &b, alpha_mkl, "G**C", Avalues_mkl, Aentries,
             Arowptrs, Arowptrs + 1, x_mkl, ldx, beta_mkl, y_mkl, ldy);
}

inline void spm_mv_block_impl_mkl(
    char mode, Kokkos::complex<double> alpha, Kokkos::complex<double> beta,
    int m, int n, int b, const int* Arowptrs, const int* Aentries,
    const Kokkos::complex<double>* Avalues, const Kokkos::complex<double>* x,
    int colx, int ldx, Kokkos::complex<double>* y, int ldy) {
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
      SCALAR const, int const, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,   \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, int const, SCALAR const*,       \
      Kokkos::LayoutLeft, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,        \
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>, SCALAR*, \
      Kokkos::LayoutLeft, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,        \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, true, COMPILE_LIBRARY> {        \
    using device_type = Kokkos::Device<EXECSPACE, Kokkos::HostSpace>;          \
    using AMatrix =                                                            \
        BsrMatrix<SCALAR const, int const, device_type,                        \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged>, int const>;         \
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
      SCALAR const, int const, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,   \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, int const, SCALAR const**,      \
      Kokkos::LayoutLeft, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,        \
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>,          \
      SCALAR**, Kokkos::LayoutLeft,                                            \
      Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,                            \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, true, true, COMPILE_LIBRARY> {  \
    using device_type = Kokkos::Device<EXECSPACE, Kokkos::HostSpace>;          \
    using AMatrix =                                                            \
        BsrMatrix<SCALAR const, int const, device_type,                        \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged>, int const>;         \
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
      int colx = static_cast<int>(X.extent(1));                                \
      int ldx  = static_cast<int>(X.stride_1());                               \
      int ldy  = static_cast<int>(Y.stride_1());                               \
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

#endif

// cuSPARSE
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#include "cusparse.h"
#include "KokkosKernels_SparseUtils_cusparse.hpp"

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
    } break;
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
#endif  // CUDA_VERSION
}

// Reference
// https://docs.nvidia.com/cuda/cusparse/index.html#bsrmm
// Several comments on bsrmm():
// - Only blockDim > 1 is supported
// - Only CUSPARSE_OPERATION_NON_TRANSPOSE is supported
// - Only CUSPARSE_MATRIX_TYPE_GENERAL is supported.
//
template <class AMatrix, class XVector, class YVector>
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
    } break;
  }

  int colx = static_cast<int>(x.extent(1));
  int ldx  = static_cast<int>(x.stride_1());
  int ldy  = static_cast<int>(y.stride_1());

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
#endif  // CUDA_VERSION
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
#endif

#undef KOKKOSSPARSE_SPMV_CUSPARSE

#define KOKKOSSPARSE_SPMV_MV_CUSPARSE(SCALAR, ORDINAL, OFFSET, LAYOUT, SPACE,  \
                                      COMPILE_LIBRARY)                         \
  template <>                                                                  \
  struct SPMV_MV_BSRMATRIX<                                                    \
      SCALAR const, ORDINAL const, Kokkos::Device<Kokkos::Cuda, SPACE>,        \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, OFFSET const, SCALAR const**,   \
      LAYOUT, Kokkos::Device<Kokkos::Cuda, SPACE>,                             \
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>,          \
      SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, SPACE>,                   \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, true, true, COMPILE_LIBRARY> {  \
    using device_type       = Kokkos::Device<Kokkos::Cuda, SPACE>;             \
    using memory_trait_type = Kokkos::MemoryTraits<Kokkos::Unmanaged>;         \
    using AMatrix = BsrMatrix<SCALAR const, ORDINAL const, device_type,        \
                              memory_trait_type, OFFSET const>;                \
    using XVector = Kokkos::View<                                              \
        SCALAR const**, LAYOUT, device_type,                                   \
        Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;       \
    using YVector =                                                            \
        Kokkos::View<SCALAR**, LAYOUT, device_type, memory_trait_type>;        \
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
KOKKOSSPARSE_SPMV_MV_CUSPARSE(double, int, int, Kokkos::LayoutLeft,
                              Kokkos::CudaSpace,
                              KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(double, int, int, Kokkos::LayoutRight,
                              Kokkos::CudaSpace,
                              KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(float, int, int, Kokkos::LayoutLeft,
                              Kokkos::CudaSpace,
                              KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(float, int, int, Kokkos::LayoutRight,
                              Kokkos::CudaSpace,
                              KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<double>, int, int,
                              Kokkos::LayoutLeft, Kokkos::CudaSpace,
                              KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<double>, int, int,
                              Kokkos::LayoutRight, Kokkos::CudaSpace,
                              KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<float>, int, int,
                              Kokkos::LayoutLeft, Kokkos::CudaSpace,
                              KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<float>, int, int,
                              Kokkos::LayoutRight, Kokkos::CudaSpace,
                              KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(double, int, int, Kokkos::LayoutLeft,
                              Kokkos::CudaUVMSpace,
                              KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(double, int, int, Kokkos::LayoutRight,
                              Kokkos::CudaUVMSpace,
                              KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(float, int, int, Kokkos::LayoutLeft,
                              Kokkos::CudaUVMSpace,
                              KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(float, int, int, Kokkos::LayoutRight,
                              Kokkos::CudaUVMSpace,
                              KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<double>, int, int,
                              Kokkos::LayoutLeft, Kokkos::CudaUVMSpace,
                              KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<double>, int, int,
                              Kokkos::LayoutRight, Kokkos::CudaUVMSpace,
                              KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<float>, int, int,
                              Kokkos::LayoutLeft, Kokkos::CudaUVMSpace,
                              KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<float>, int, int,
                              Kokkos::LayoutRight, Kokkos::CudaUVMSpace,
                              KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
#endif

#undef KOKKOSSPARSE_SPMV_MV_CUSPARSE

}  // namespace Impl
}  // namespace Experimental
}  // namespace KokkosSparse

#endif

#endif  // KOKKOSKERNELS_KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_DECL_HPP
