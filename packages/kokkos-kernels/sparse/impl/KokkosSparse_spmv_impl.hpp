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

#ifndef KOKKOSSPARSE_IMPL_SPMV_DEF_HPP_
#define KOKKOSSPARSE_IMPL_SPMV_DEF_HPP_

#include <sstream>

#include "KokkosKernels_Controls.hpp"
#include "Kokkos_InnerProductSpaceTraits.hpp"
#include "KokkosBlas1_scal.hpp"
#include "KokkosKernels_ExecSpaceUtils.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_spmv_handle.hpp"
#include "KokkosSparse_spmv_impl_omp.hpp"
#include "KokkosSparse_spmv_impl_merge.hpp"
#include "KokkosKernels_Error.hpp"

namespace KokkosSparse {
namespace Impl {

// This TransposeFunctor is functional, but not necessarily performant.
template <class execution_space, class AMatrix, class XVector, class YVector, bool conjugate>
struct SPMV_Transpose_Functor {
  typedef typename AMatrix::non_const_ordinal_type ordinal_type;
  typedef typename AMatrix::non_const_value_type value_type;
  typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type team_member;
  typedef Kokkos::ArithTraits<value_type> ATV;
  typedef typename YVector::non_const_value_type coefficient_type;
  typedef typename YVector::non_const_value_type y_value_type;

  const coefficient_type alpha;
  AMatrix m_A;
  XVector m_x;
  YVector m_y;
  ordinal_type rows_per_team = 0;

  SPMV_Transpose_Functor(const coefficient_type& alpha_, const AMatrix& m_A_, const XVector& m_x_, const YVector& m_y_)
      : alpha(alpha_), m_A(m_A_), m_x(m_x_), m_y(m_y_) {}

  KOKKOS_INLINE_FUNCTION void operator()(const ordinal_type iRow) const {
    const auto row                = m_A.rowConst(iRow);
    const ordinal_type row_length = row.length;
    for (ordinal_type iEntry = 0; iEntry < row_length; iEntry++) {
      const value_type val   = conjugate ? ATV::conj(row.value(iEntry)) : row.value(iEntry);
      const ordinal_type ind = row.colidx(iEntry);
      Kokkos::atomic_add(&m_y(ind), static_cast<y_value_type>(alpha * val * m_x(iRow)));
    }
  }

  KOKKOS_INLINE_FUNCTION void operator()(const team_member& dev) const {
    const ordinal_type teamWork = dev.league_rank() * rows_per_team;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(dev, rows_per_team), [&](ordinal_type loop) {
      // iRow represents a row of the matrix, so its correct type is
      // ordinal_type.
      const ordinal_type iRow = teamWork + loop;
      if (iRow >= m_A.numRows()) {
        return;
      }

      const auto row                = m_A.rowConst(iRow);
      const ordinal_type row_length = row.length;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(dev, row_length), [&](ordinal_type iEntry) {
        const value_type val   = conjugate ? ATV::conj(row.value(iEntry)) : row.value(iEntry);
        const ordinal_type ind = row.colidx(iEntry);
        Kokkos::atomic_add(&m_y(ind), static_cast<y_value_type>(alpha * val * m_x(iRow)));
      });
    });
  }
};

template <class execution_space, class AMatrix, class XVector, class YVector, int dobeta, bool conjugate>
struct SPMV_Functor {
  typedef typename AMatrix::non_const_ordinal_type ordinal_type;
  typedef typename AMatrix::non_const_value_type value_type;
  typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type team_member;
  typedef Kokkos::ArithTraits<value_type> ATV;

  const value_type alpha;
  AMatrix m_A;
  XVector m_x;
  const value_type beta;
  YVector m_y;

  const ordinal_type rows_per_team;

  SPMV_Functor(const value_type alpha_, const AMatrix m_A_, const XVector m_x_, const value_type beta_,
               const YVector m_y_, const int rows_per_team_)
      : alpha(alpha_), m_A(m_A_), m_x(m_x_), beta(beta_), m_y(m_y_), rows_per_team(rows_per_team_) {
    static_assert(static_cast<int>(XVector::rank) == 1, "XVector must be a rank 1 View.");
    static_assert(static_cast<int>(YVector::rank) == 1, "YVector must be a rank 1 View.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const ordinal_type iRow) const {
    using y_value_type = typename YVector::non_const_value_type;
    if (iRow >= m_A.numRows()) {
      return;
    }
    const KokkosSparse::SparseRowViewConst<AMatrix> row = m_A.rowConst(iRow);
    const ordinal_type row_length                       = static_cast<ordinal_type>(row.length);
    y_value_type sum                                    = 0;

    for (ordinal_type iEntry = 0; iEntry < row_length; iEntry++) {
      const value_type val = conjugate ? ATV::conj(row.value(iEntry)) : row.value(iEntry);
      sum += val * m_x(row.colidx(iEntry));
    }

    sum *= alpha;

    if (dobeta == 0) {
      m_y(iRow) = sum;
    } else {
      m_y(iRow) = beta * m_y(iRow) + sum;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member& dev) const {
    using y_value_type = typename YVector::non_const_value_type;

    Kokkos::parallel_for(Kokkos::TeamThreadRange(dev, 0, rows_per_team), [&](const ordinal_type& loop) {
      const ordinal_type iRow = static_cast<ordinal_type>(dev.league_rank()) * rows_per_team + loop;
      if (iRow >= m_A.numRows()) {
        return;
      }
      const KokkosSparse::SparseRowViewConst<AMatrix> row = m_A.rowConst(iRow);
      const ordinal_type row_length                       = static_cast<ordinal_type>(row.length);
      y_value_type sum                                    = 0;

      Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(dev, row_length),
          [&](const ordinal_type& iEntry, y_value_type& lsum) {
            const value_type val = conjugate ? ATV::conj(row.value(iEntry)) : row.value(iEntry);
            lsum += val * m_x(row.colidx(iEntry));
          },
          sum);

      Kokkos::single(Kokkos::PerThread(dev), [&]() {
        sum *= alpha;

        if (dobeta == 0) {
          m_y(iRow) = sum;
        } else {
          m_y(iRow) = beta * m_y(iRow) + sum;
        }
      });
    });
  }
};

template <class execution_space>
int64_t spmv_launch_parameters(int64_t numRows, int64_t nnz, int64_t rows_per_thread, int& team_size,
                               int& vector_length) {
  int64_t rows_per_team;
  int64_t nnz_per_row = nnz / numRows;

  if (nnz_per_row < 1) nnz_per_row = 1;

  int max_vector_length = 1;
#ifdef KOKKOS_ENABLE_CUDA
  if (std::is_same<execution_space, Kokkos::Cuda>::value) max_vector_length = 32;
#endif
#ifdef KOKKOS_ENABLE_HIP
  if (std::is_same<execution_space, Kokkos::HIP>::value) max_vector_length = 64;
#endif

  if (vector_length < 1) {
    vector_length = 1;
    while (vector_length < max_vector_length && vector_length * 6 < nnz_per_row) vector_length *= 2;
  }

  // Determine rows per thread
  if (rows_per_thread < 1) {
    if (KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>())
      rows_per_thread = 1;
    else {
      if (nnz_per_row < 20 && nnz > 5000000) {
        rows_per_thread = 256;
      } else
        rows_per_thread = 64;
    }
  }

  if (team_size < 1) {
    if (KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>()) {
      team_size = 256 / vector_length;
    } else {
      team_size = 1;
    }
  }

  rows_per_team = rows_per_thread * team_size;

  if (rows_per_team < 0) {
    int64_t nnz_per_team = 4096;
    int64_t conc         = execution_space().concurrency();
    while ((conc * nnz_per_team * 4 > nnz) && (nnz_per_team > 256)) nnz_per_team /= 2;
    rows_per_team = (nnz_per_team + nnz_per_row - 1) / nnz_per_row;
  }

  return rows_per_team;
}

// spmv_beta_no_transpose: version for CPU execution spaces (RangePolicy or
// trivial serial impl used)
template <class execution_space, class Handle, class AMatrix, class XVector, class YVector, int dobeta, bool conjugate,
          typename std::enable_if<!KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>()>::type* = nullptr>
static void spmv_beta_no_transpose(const execution_space& exec, Handle* handle,
                                   typename YVector::const_value_type& alpha, const AMatrix& A, const XVector& x,
                                   typename YVector::const_value_type& beta, const YVector& y) {
  typedef typename AMatrix::non_const_ordinal_type ordinal_type;

  if (A.numRows() <= static_cast<ordinal_type>(0)) {
    return;
  }
#if defined(KOKKOS_ENABLE_SERIAL)
  if (std::is_same<execution_space, Kokkos::Serial>::value) {
    /// serial impl
    typedef typename AMatrix::non_const_value_type value_type;
    typedef typename AMatrix::non_const_size_type size_type;
    typedef Kokkos::ArithTraits<value_type> ATV;

    const size_type* KOKKOS_RESTRICT row_map_ptr    = A.graph.row_map.data();
    const ordinal_type* KOKKOS_RESTRICT col_idx_ptr = A.graph.entries.data();
    const value_type* KOKKOS_RESTRICT values_ptr    = A.values.data();

    typename YVector::non_const_value_type* KOKKOS_RESTRICT y_ptr = y.data();
    typename XVector::const_value_type* KOKKOS_RESTRICT x_ptr     = x.data();

    const typename YVector::non_const_value_type zero(0);
    const ordinal_type nrow = A.numRows();
    if (alpha == zero) {
      if (dobeta == 0) {
        /// not working with kkosDev2_CUDA110_GCC92_cpp17/
        /// memset(y_ptr, 0, sizeof(typename YVector::value_type)*nrow);
        for (int i = 0; i < nrow; ++i) y_ptr[i] = zero;
      } else if (dobeta == 1) {
        /// so nothing
      } else {
        for (int i = 0; i < nrow; ++i) y_ptr[i] *= beta;
      }
    } else {
      for (int i = 0; i < nrow; ++i) {
        const int jbeg = row_map_ptr[i];
        const int jend = row_map_ptr[i + 1];
        int j          = jbeg;

        {
          const int jdist = (jend - jbeg) / 4;
          typename YVector::non_const_value_type tmp1(0), tmp2(0), tmp3(0), tmp4(0);
          for (int jj = 0; jj < jdist; ++jj) {
            const value_type value1                   = conjugate ? ATV::conj(values_ptr[j]) : values_ptr[j];
            const value_type value2                   = conjugate ? ATV::conj(values_ptr[j + 1]) : values_ptr[j + 1];
            const value_type value3                   = conjugate ? ATV::conj(values_ptr[j + 2]) : values_ptr[j + 2];
            const value_type value4                   = conjugate ? ATV::conj(values_ptr[j + 3]) : values_ptr[j + 3];
            const int col_idx1                        = col_idx_ptr[j];
            const int col_idx2                        = col_idx_ptr[j + 1];
            const int col_idx3                        = col_idx_ptr[j + 2];
            const int col_idx4                        = col_idx_ptr[j + 3];
            const typename XVector::value_type x_val1 = x_ptr[col_idx1];
            const typename XVector::value_type x_val2 = x_ptr[col_idx2];
            const typename XVector::value_type x_val3 = x_ptr[col_idx3];
            const typename XVector::value_type x_val4 = x_ptr[col_idx4];
            tmp1 += value1 * x_val1;
            tmp2 += value2 * x_val2;
            tmp3 += value3 * x_val3;
            tmp4 += value4 * x_val4;
            j += 4;
          }
          for (; j < jend; ++j) {
            const value_type value = conjugate ? ATV::conj(values_ptr[j]) : values_ptr[j];
            const int col_idx      = col_idx_ptr[j];
            tmp1 += value * x_ptr[col_idx];
          }
          if (dobeta == 0) {
            y_ptr[i] = alpha * (tmp1 + tmp2 + tmp3 + tmp4);
          } else if (dobeta == 1) {
            y_ptr[i] += alpha * (tmp1 + tmp2 + tmp3 + tmp4);
          } else {
            const auto y_val = y_ptr[i] * beta;
            y_ptr[i]         = y_val + alpha * (tmp1 + tmp2 + tmp3 + tmp4);
          }
        }
      }
    }
    return;
  }
#endif

#ifdef KOKKOS_ENABLE_OPENMP
  if ((std::is_same<execution_space, Kokkos::OpenMP>::value) &&
      (std::is_same<typename std::remove_cv<typename AMatrix::value_type>::type, double>::value) &&
      (std::is_same<typename XVector::non_const_value_type, double>::value) &&
      (std::is_same<typename YVector::non_const_value_type, double>::value) &&
      ((int)A.graph.row_block_offsets.extent(0) == (int)omp_get_max_threads() + 1) &&
      (((uintptr_t)(const void*)(x.data()) % 64) == 0) && (((uintptr_t)(const void*)(y.data()) % 64) == 0) &&
      !conjugate) {
    // Note BMK: this case is typically not called in practice even for OpenMP,
    // since it requires row_block_offsets to have been computed in the graph.
    // Also, as this is raw OpenMP the execution space instance is not used
    spmv_raw_openmp_no_transpose<AMatrix, XVector, YVector>(alpha, A, x, beta, y);
    return;
  }
#endif

  bool use_dynamic_schedule = handle->force_dynamic_schedule;
  bool use_static_schedule  = handle->force_static_schedule;
  SPMV_Functor<execution_space, AMatrix, XVector, YVector, dobeta, conjugate> func(alpha, A, x, beta, y, 1);
  if (((A.nnz() > 10000000) || use_dynamic_schedule) && !use_static_schedule)
    Kokkos::parallel_for("KokkosSparse::spmv<NoTranspose,Dynamic>",
                         Kokkos::RangePolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>>(exec, 0, A.numRows()),
                         func);
  else
    Kokkos::parallel_for("KokkosSparse::spmv<NoTranspose,Static>",
                         Kokkos::RangePolicy<execution_space, Kokkos::Schedule<Kokkos::Static>>(exec, 0, A.numRows()),
                         func);
}

// spmv_beta_no_transpose: version for GPU execution spaces (TeamPolicy used)
template <class execution_space, class Handle, class AMatrix, class XVector, class YVector, int dobeta, bool conjugate,
          typename std::enable_if<KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>()>::type* = nullptr>
static void spmv_beta_no_transpose(const execution_space& exec, Handle* handle,
                                   typename YVector::const_value_type& alpha, const AMatrix& A, const XVector& x,
                                   typename YVector::const_value_type& beta, const YVector& y) {
  typedef typename AMatrix::non_const_ordinal_type ordinal_type;

  if (A.numRows() <= static_cast<ordinal_type>(0)) {
    return;
  }

  bool use_dynamic_schedule = handle->force_dynamic_schedule;
  bool use_static_schedule  = handle->force_static_schedule;
  int team_size             = handle->team_size;
  int vector_length         = handle->vector_length;
  int64_t rows_per_thread   = handle->rows_per_thread;

  int64_t rows_per_team =
      spmv_launch_parameters<execution_space>(A.numRows(), A.nnz(), rows_per_thread, team_size, vector_length);
  int64_t worksets = (y.extent(0) + rows_per_team - 1) / rows_per_team;

  SPMV_Functor<execution_space, AMatrix, XVector, YVector, dobeta, conjugate> func(alpha, A, x, beta, y, rows_per_team);

  if (((A.nnz() > 10000000) || use_dynamic_schedule) && !use_static_schedule) {
    Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>> policy(1, 1);
    if (team_size < 0)
      policy = Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>>(exec, worksets, Kokkos::AUTO,
                                                                                      vector_length);
    else
      policy = Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>>(exec, worksets, team_size,
                                                                                      vector_length);
    Kokkos::parallel_for("KokkosSparse::spmv<NoTranspose,Dynamic>", policy, func);
  } else {
    Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>> policy(1, 1);
    if (team_size < 0)
      policy = Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>>(exec, worksets, Kokkos::AUTO,
                                                                                     vector_length);
    else
      policy = Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>>(exec, worksets, team_size,
                                                                                     vector_length);
    Kokkos::parallel_for("KokkosSparse::spmv<NoTranspose,Static>", policy, func);
  }
}

// spmv_beta_transpose: version for CPU execution spaces (RangePolicy or trivial
// serial impl used)
template <class execution_space, class AMatrix, class XVector, class YVector, int dobeta, bool conjugate,
          typename std::enable_if<!KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>()>::type* = nullptr>
static void spmv_beta_transpose(const execution_space& exec, typename YVector::const_value_type& alpha,
                                const AMatrix& A, const XVector& x, typename YVector::const_value_type& beta,
                                const YVector& y) {
  using ordinal_type  = typename AMatrix::non_const_ordinal_type;
  using size_type     = typename AMatrix::non_const_size_type;
  using y_scalar_type = typename YVector::non_const_value_type;

  if (A.numRows() <= static_cast<ordinal_type>(0)) {
    return;
  }

  // We need to scale y first ("scaling" by zero just means filling
  // with zeros), since the functor works by atomic-adding into y.
  if (0 == dobeta || y_scalar_type(0) == beta) {
    Kokkos::deep_copy(exec, y, y_scalar_type(0));
  } else if (dobeta != 1) {
    KokkosBlas::scal(exec, y, beta, y);
  }

#if defined(KOKKOS_ENABLE_SERIAL) || defined(KOKKOS_ENABLE_OPENMP) || defined(KOKKOS_ENABLE_THREADS)
  {
    if (exec.concurrency() == 1) {
      /// serial impl
      typedef typename AMatrix::non_const_value_type value_type;
      typedef Kokkos::ArithTraits<value_type> ATV;
      const size_type* KOKKOS_RESTRICT row_map_ptr    = A.graph.row_map.data();
      const ordinal_type* KOKKOS_RESTRICT col_idx_ptr = A.graph.entries.data();
      const value_type* KOKKOS_RESTRICT values_ptr    = A.values.data();

      typename YVector::value_type* KOKKOS_RESTRICT y_ptr = y.data();
      typename XVector::value_type* KOKKOS_RESTRICT x_ptr = x.data();

      const typename YVector::non_const_value_type zero(0);
      const ordinal_type nrow = A.numRows();
      if (alpha == zero) {
        /// do nothing
      } else {
        for (int i = 0; i < nrow; ++i) {
          const int jbeg  = row_map_ptr[i];
          const int jend  = row_map_ptr[i + 1];
          const int jdist = (jend - jbeg) / 4;

          const typename XVector::const_value_type x_val = alpha * x_ptr[i];
          int j                                          = jbeg;
          for (int jj = 0; jj < jdist; ++jj) {
            const value_type value1 = conjugate ? ATV::conj(values_ptr[j]) : values_ptr[j];
            const value_type value2 = conjugate ? ATV::conj(values_ptr[j + 1]) : values_ptr[j + 1];
            const value_type value3 = conjugate ? ATV::conj(values_ptr[j + 2]) : values_ptr[j + 2];
            const value_type value4 = conjugate ? ATV::conj(values_ptr[j + 3]) : values_ptr[j + 3];
            const int col_idx1      = col_idx_ptr[j];
            const int col_idx2      = col_idx_ptr[j + 1];
            const int col_idx3      = col_idx_ptr[j + 2];
            const int col_idx4      = col_idx_ptr[j + 3];
            y_ptr[col_idx1] += value1 * x_val;
            y_ptr[col_idx2] += value2 * x_val;
            y_ptr[col_idx3] += value3 * x_val;
            y_ptr[col_idx4] += value4 * x_val;
            j += 4;
          }
          for (; j < jend; ++j) {
            const value_type value = conjugate ? ATV::conj(values_ptr[j]) : values_ptr[j];
            const int col_idx      = col_idx_ptr[j];
            y_ptr[col_idx] += value * x_val;
          }
        }
      }
      return;
    }
  }
#endif

  typedef SPMV_Transpose_Functor<execution_space, AMatrix, XVector, YVector, conjugate> OpType;
  typename AMatrix::const_ordinal_type nrow = A.numRows();
  Kokkos::parallel_for("KokkosSparse::spmv<Transpose>", Kokkos::RangePolicy<execution_space>(exec, 0, nrow),
                       OpType(alpha, A, x, y));
}

// spmv_beta_transpose: version for GPU execution spaces (TeamPolicy used)
template <class execution_space, class AMatrix, class XVector, class YVector, int dobeta, bool conjugate,
          typename std::enable_if<KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>()>::type* = nullptr>
static void spmv_beta_transpose(const execution_space& exec, typename YVector::const_value_type& alpha,
                                const AMatrix& A, const XVector& x, typename YVector::const_value_type& beta,
                                const YVector& y) {
  using ordinal_type  = typename AMatrix::non_const_ordinal_type;
  using size_type     = typename AMatrix::non_const_size_type;
  using y_scalar_type = typename YVector::non_const_value_type;

  if (A.numRows() <= static_cast<ordinal_type>(0)) {
    return;
  }

  // We need to scale y first ("scaling" by zero just means filling
  // with zeros), since the functor works by atomic-adding into y.
  if (0 == dobeta || y_scalar_type(0) == beta) {
    Kokkos::deep_copy(exec, y, y_scalar_type(0));
  } else if (dobeta != 1) {
    KokkosBlas::scal(exec, y, beta, y);
  }

  // Assuming that no row contains duplicate entries, NNZPerRow
  // cannot be more than the number of columns of the matrix.  Thus,
  // the appropriate type is ordinal_type.
  const ordinal_type NNZPerRow = A.nnz() / A.numRows();

  int vector_length     = 1;
  int max_vector_length = 1;
#ifdef KOKKOS_ENABLE_CUDA
  if (std::is_same<execution_space, Kokkos::Cuda>::value) max_vector_length = 32;
#endif
#ifdef KOKKOS_ENABLE_HIP
  if (std::is_same<execution_space, Kokkos::HIP>::value) max_vector_length = 64;
#endif
  while ((vector_length * 2 * 3 <= NNZPerRow) && (vector_length < max_vector_length)) vector_length *= 2;

  typedef SPMV_Transpose_Functor<execution_space, AMatrix, XVector, YVector, conjugate> OpType;

  typename AMatrix::const_ordinal_type nrow = A.numRows();

  OpType op(alpha, A, x, y);

  const ordinal_type rows_per_thread = RowsPerThread<execution_space>(NNZPerRow);
  const ordinal_type team_size = Kokkos::TeamPolicy<execution_space>(exec, rows_per_thread, Kokkos::AUTO, vector_length)
                                     .team_size_recommended(op, Kokkos::ParallelForTag());
  const ordinal_type rows_per_team = rows_per_thread * team_size;
  op.rows_per_team                 = rows_per_team;
  const size_type nteams           = (nrow + rows_per_team - 1) / rows_per_team;
  Kokkos::parallel_for("KokkosSparse::spmv<Transpose>",
                       Kokkos::TeamPolicy<execution_space>(exec, nteams, team_size, vector_length), op);
}

template <class execution_space, class Handle, class AMatrix, class XVector, class YVector, int dobeta>
static void spmv_beta(const execution_space& exec, Handle* handle, const char mode[],
                      typename YVector::const_value_type& alpha, const AMatrix& A, const XVector& x,
                      typename YVector::const_value_type& beta, const YVector& y) {
  if (mode[0] == NoTranspose[0]) {
    if (handle->algo == SPMV_MERGE_PATH || handle->algo == SPMV_NATIVE_MERGE_PATH) {
      SpmvMergeHierarchical<execution_space, AMatrix, XVector, YVector>::spmv(exec, mode, alpha, A, x, beta, y);
    } else {
      spmv_beta_no_transpose<execution_space, Handle, AMatrix, XVector, YVector, dobeta, false>(exec, handle, alpha, A,
                                                                                                x, beta, y);
    }
  } else if (mode[0] == Conjugate[0]) {
    if (handle->algo == SPMV_MERGE_PATH || handle->algo == SPMV_NATIVE_MERGE_PATH) {
      SpmvMergeHierarchical<execution_space, AMatrix, XVector, YVector>::spmv(exec, mode, alpha, A, x, beta, y);
    } else {
      spmv_beta_no_transpose<execution_space, Handle, AMatrix, XVector, YVector, dobeta, true>(exec, handle, alpha, A,
                                                                                               x, beta, y);
    }
  } else if (mode[0] == Transpose[0]) {
    spmv_beta_transpose<execution_space, AMatrix, XVector, YVector, dobeta, false>(exec, alpha, A, x, beta, y);
  } else if (mode[0] == ConjugateTranspose[0]) {
    spmv_beta_transpose<execution_space, AMatrix, XVector, YVector, dobeta, true>(exec, alpha, A, x, beta, y);
  } else {
    std::stringstream ss;
    ss << __FILE__ << ":" << __LINE__ << " Invalid transpose mode " << mode << " for KokkosSparse::spmv()";
    KokkosKernels::Impl::throw_runtime_exception(ss.str());
  }
}

// Functor for implementing transpose and conjugate transpose sparse
// matrix-vector multiply with multivector (2-D View) input and
// output.  This functor works, but is not necessarily performant.
template <class execution_space, class AMatrix, class XVector, class YVector, int doalpha, int dobeta, bool conjugate>
struct SPMV_MV_Transpose_Functor {
  typedef typename AMatrix::non_const_ordinal_type ordinal_type;
  typedef typename AMatrix::non_const_value_type A_value_type;
  typedef typename YVector::non_const_value_type y_value_type;
  typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type team_member;
  typedef typename YVector::non_const_value_type coefficient_type;

  const coefficient_type alpha;
  AMatrix m_A;
  XVector m_x;
  const coefficient_type beta;
  YVector m_y;

  const ordinal_type n;
  ordinal_type rows_per_team = 0;

  SPMV_MV_Transpose_Functor(const coefficient_type& alpha_, const AMatrix& m_A_, const XVector& m_x_,
                            const coefficient_type& beta_, const YVector& m_y_)
      : alpha(alpha_), m_A(m_A_), m_x(m_x_), beta(beta_), m_y(m_y_), n(m_x_.extent(1)) {}

  KOKKOS_INLINE_FUNCTION void operator()(const ordinal_type iRow) const {
    const auto row                = m_A.rowConst(iRow);
    const ordinal_type row_length = row.length;

    for (ordinal_type iEntry = 0; iEntry < row_length; iEntry++) {
      const A_value_type val =
          conjugate ? Kokkos::ArithTraits<A_value_type>::conj(row.value(iEntry)) : row.value(iEntry);
      const ordinal_type ind = row.colidx(iEntry);

      if (doalpha != 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (ordinal_type k = 0; k < n; ++k) {
          Kokkos::atomic_add(&m_y(ind, k), static_cast<y_value_type>(alpha * val * m_x(iRow, k)));
        }
      } else {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (ordinal_type k = 0; k < n; ++k) {
          Kokkos::atomic_add(&m_y(ind, k), static_cast<y_value_type>(val * m_x(iRow, k)));
        }
      }
    }
  }

  KOKKOS_INLINE_FUNCTION void operator()(const team_member& dev) const {
    const ordinal_type teamWork = dev.league_rank() * rows_per_team;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(dev, rows_per_team), [&](ordinal_type loop) {
      // iRow represents a row of the matrix, so its correct type is
      // ordinal_type.
      const ordinal_type iRow = teamWork + loop;
      if (iRow >= m_A.numRows()) {
        return;
      }

      const auto row                = m_A.rowConst(iRow);
      const ordinal_type row_length = row.length;

      Kokkos::parallel_for(Kokkos::ThreadVectorRange(dev, row_length), [&](ordinal_type iEntry) {
        const A_value_type val =
            conjugate ? Kokkos::ArithTraits<A_value_type>::conj(row.value(iEntry)) : row.value(iEntry);
        const ordinal_type ind = row.colidx(iEntry);

        if (doalpha != 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (ordinal_type k = 0; k < n; ++k) {
            Kokkos::atomic_add(&m_y(ind, k), static_cast<y_value_type>(alpha * val * m_x(iRow, k)));
          }
        } else {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (ordinal_type k = 0; k < n; ++k) {
            Kokkos::atomic_add(&m_y(ind, k), static_cast<y_value_type>(val * m_x(iRow, k)));
          }
        }
      });
    });
  }
};

template <class execution_space, class AMatrix, class XVector, class YVector, int doalpha, int dobeta, bool conjugate>
struct SPMV_MV_LayoutLeft_Functor {
  typedef typename AMatrix::non_const_ordinal_type ordinal_type;
  typedef typename AMatrix::non_const_value_type A_value_type;
  typedef typename YVector::non_const_value_type y_value_type;
  typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type team_member;
  typedef typename YVector::non_const_value_type coefficient_type;

  const coefficient_type alpha;
  AMatrix m_A;
  XVector m_x;
  const coefficient_type beta;
  YVector m_y;
  //! The number of columns in the input and output MultiVectors.
  ordinal_type n;
  ordinal_type rows_per_thread;
  int vector_length;

  SPMV_MV_LayoutLeft_Functor(const coefficient_type& alpha_, const AMatrix& m_A_, const XVector& m_x_,
                             const coefficient_type& beta_, const YVector& m_y_, const ordinal_type rows_per_thread_,
                             int vector_length_)
      : alpha(alpha_),
        m_A(m_A_),
        m_x(m_x_),
        beta(beta_),
        m_y(m_y_),
        n(m_x_.extent(1)),
        rows_per_thread(rows_per_thread_),
        vector_length(vector_length_) {}

  template <int UNROLL>
  KOKKOS_INLINE_FUNCTION void strip_mine(const team_member& dev, const ordinal_type& iRow,
                                         const ordinal_type& kk) const {
    y_value_type sum[UNROLL];

#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
    for (int k = 0; k < UNROLL; ++k) {
      sum[k] = Kokkos::ArithTraits<y_value_type>::zero();
    }

    const auto row = m_A.rowConst(iRow);

    // The correct type of iEntry is ordinal_type, the type of the
    // number of columns in the (local) matrix.  This is because we
    // assume either that rows have no duplicate entries, or that rows
    // never have enough duplicate entries to overflow ordinal_type.

    Kokkos::parallel_for(Kokkos::ThreadVectorRange(dev, row.length), [&](ordinal_type iEntry) {
      const A_value_type val =
          conjugate ? Kokkos::ArithTraits<A_value_type>::conj(row.value(iEntry)) : row.value(iEntry);
      const ordinal_type ind = row.colidx(iEntry);
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        sum[k] += val * m_x(ind, kk + k);
      }
    });

    if (doalpha == -1) {
      for (int ii = 0; ii < UNROLL; ++ii) {
        y_value_type sumt;
        Kokkos::parallel_reduce(
            Kokkos::ThreadVectorRange(dev, vector_length),
            [&](ordinal_type, y_value_type& lsum) {
              // in this context, sum[ii] is a partial sum ii on one of the
              // vector lanes.
              lsum -= sum[ii];
            },
            sumt);
        sum[ii] = sumt;
        // that was an all-reduce, so sum[ii] is the same on every vector lane
      }
    } else {
      for (int ii = 0; ii < UNROLL; ++ii) {
        y_value_type sumt;
        Kokkos::parallel_reduce(
            Kokkos::ThreadVectorRange(dev, vector_length),
            [&](ordinal_type, y_value_type& lsum) {
              // in this context, sum[ii] is a partial sum ii on one of the
              // vector lanes.
              lsum += sum[ii];
            },
            sumt);
        if (doalpha == 1)
          sum[ii] = sumt;
        else
          sum[ii] = sumt * alpha;
      }
    }

    if (dobeta == 0) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(dev, UNROLL), [&](ordinal_type k) { m_y(iRow, kk + k) = sum[k]; });
    } else if (dobeta == 1) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(dev, UNROLL),
                           [&](ordinal_type k) { m_y(iRow, kk + k) = m_y(iRow, kk + k) + sum[k]; });
    } else if (dobeta == -1) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(dev, UNROLL),
                           [&](ordinal_type k) { m_y(iRow, kk + k) = -m_y(iRow, kk + k) + sum[k]; });
    } else {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(dev, UNROLL),
                           [&](ordinal_type k) { m_y(iRow, kk + k) = beta * m_y(iRow, kk + k) + sum[k]; });
    }
  }

  template <int UNROLL>
  KOKKOS_INLINE_FUNCTION void strip_mine(const ordinal_type& iRow, const ordinal_type& kk) const {
    y_value_type sum[UNROLL];

#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
    for (int k = 0; k < UNROLL; ++k) {
      sum[k] = Kokkos::ArithTraits<y_value_type>::zero();
    }

    const auto row = m_A.rowConst(iRow);

    // The correct type of iEntry is ordinal_type, the type of the
    // number of columns in the (local) matrix.  This is because we
    // assume either that rows have no duplicate entries, or that rows
    // never have enough duplicate entries to overflow ordinal_type.

    for (ordinal_type iEntry = 0; iEntry < row.length; iEntry++) {
      const A_value_type val =
          conjugate ? Kokkos::ArithTraits<A_value_type>::conj(row.value(iEntry)) : row.value(iEntry);
      const ordinal_type ind = row.colidx(iEntry);
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        if (doalpha == 1)
          sum[k] += val * m_x(ind, kk + k);
        else if (doalpha == -1)
          sum[k] -= val * m_x(ind, kk + k);
        else
          sum[k] += alpha * val * m_x(ind, kk + k);
      }
    }

    if (dobeta == 0) {
      for (ordinal_type k = 0; k < UNROLL; k++) m_y(iRow, kk + k) = sum[k];
    } else if (dobeta == 1) {
      for (ordinal_type k = 0; k < UNROLL; k++) m_y(iRow, kk + k) = m_y(iRow, kk + k) + sum[k];
    } else if (dobeta == -1) {
      for (ordinal_type k = 0; k < UNROLL; k++) m_y(iRow, kk + k) = -m_y(iRow, kk + k) + sum[k];
    } else {
      for (ordinal_type k = 0; k < UNROLL; k++) m_y(iRow, kk + k) = beta * m_y(iRow, kk + k) + sum[k];
    }
  }

  KOKKOS_INLINE_FUNCTION void strip_mine_1(const team_member& dev, const ordinal_type& iRow) const {
    const auto row = m_A.rowConst(iRow);

    // The correct type of iEntry is ordinal_type, the type of the
    // number of columns in the (local) matrix.  This is because we
    // assume either that rows have no duplicate entries, or that rows
    // never have enough duplicate entries to overflow ordinal_type.

    y_value_type sum;
    Kokkos::parallel_reduce(
        Kokkos::ThreadVectorRange(dev, row.length),
        [&](ordinal_type iEntry, y_value_type& lsum) {
          const A_value_type val =
              conjugate ? Kokkos::ArithTraits<A_value_type>::conj(row.value(iEntry)) : row.value(iEntry);
          lsum += val * m_x(row.colidx(iEntry), 0);
        },
        sum);
    Kokkos::single(Kokkos::PerThread(dev), [&]() {
      if (doalpha == -1) {
        sum = -sum;
      } else if (doalpha * doalpha != 1) {
        sum *= alpha;
      }

      if (dobeta == 0) {
        m_y(iRow, 0) = sum;
      } else if (dobeta == 1) {
        m_y(iRow, 0) += sum;
      } else if (dobeta == -1) {
        m_y(iRow, 0) = -m_y(iRow, 0) + sum;
      } else {
        m_y(iRow, 0) = beta * m_y(iRow, 0) + sum;
      }
    });
  }

  KOKKOS_INLINE_FUNCTION void strip_mine_1(const ordinal_type& iRow) const {
    const auto row = m_A.rowConst(iRow);

    // The correct type of iEntry is ordinal_type, the type of the
    // number of columns in the (local) matrix.  This is because we
    // assume either that rows have no duplicate entries, or that rows
    // never have enough duplicate entries to overflow ordinal_type.

    y_value_type sum = y_value_type();
    for (ordinal_type iEntry = 0; iEntry < row.length; iEntry++) {
      const A_value_type val =
          conjugate ? Kokkos::ArithTraits<A_value_type>::conj(row.value(iEntry)) : row.value(iEntry);
      sum += val * m_x(row.colidx(iEntry), 0);
    }
    if (doalpha == -1) {
      sum = -sum;
    } else if (doalpha != 1) {
      sum *= alpha;
    }

    if (dobeta == 0) {
      m_y(iRow, 0) = sum;
    } else if (dobeta == 1) {
      m_y(iRow, 0) += sum;
    } else if (dobeta == -1) {
      m_y(iRow, 0) = -m_y(iRow, 0) + sum;
    } else {
      m_y(iRow, 0) = beta * m_y(iRow, 0) + sum;
    }
  }

  KOKKOS_INLINE_FUNCTION void operator()(const ordinal_type& iRow) const {
    // mfh 20 Mar 2015, 07 Jun 2016: This is ordinal_type because it
    // needs to have the same type as n.
    ordinal_type kk = 0;

#ifdef KOKKOS_FAST_COMPILE
    for (; kk + 4 <= n; kk += 4) {
      strip_mine<4>(dev, iRow, kk);
    }
    for (; kk < n; ++kk) {
      strip_mine<1>(dev, iRow, kk);
    }
#else
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
    if ((n > 8) && (n % 8 == 1)) {
      strip_mine<9>(iRow, kk);
      kk += 9;
    }
    for (; kk + 8 <= n; kk += 8) strip_mine<8>(iRow, kk);
    if (kk < n) {
      switch (n - kk) {
#else   // NOT a GPU
    if ((n > 16) && (n % 16 == 1)) {
      strip_mine<17>(iRow, kk);
      kk += 17;
    }

    for (; kk + 16 <= n; kk += 16) {
      strip_mine<16>(iRow, kk);
    }

    if (kk < n) {
      switch (n - kk) {
        case 15: strip_mine<15>(iRow, kk); break;

        case 14: strip_mine<14>(iRow, kk); break;

        case 13: strip_mine<13>(iRow, kk); break;

        case 12: strip_mine<12>(iRow, kk); break;

        case 11: strip_mine<11>(iRow, kk); break;

        case 10: strip_mine<10>(iRow, kk); break;

        case 9: strip_mine<9>(iRow, kk); break;

        case 8: strip_mine<8>(iRow, kk); break;
#endif  // if/else: __CUDA_ARCH__ or __HIP_DEVICE_COMPILE__
        case 7: strip_mine<7>(iRow, kk); break;

        case 6: strip_mine<6>(iRow, kk); break;

        case 5: strip_mine<5>(iRow, kk); break;

        case 4: strip_mine<4>(iRow, kk); break;

        case 3: strip_mine<3>(iRow, kk); break;

        case 2: strip_mine<2>(iRow, kk); break;

        case 1: strip_mine_1(iRow); break;
      }
    }
#endif  // KOKKOS_FAST_COMPILE
  }

  KOKKOS_INLINE_FUNCTION void operator()(const team_member& dev) const {
    for (ordinal_type loop = 0; loop < rows_per_thread; ++loop) {
      // iRow indexes over (local) rows of the matrix, so its correct
      // type is ordinal_type.

      const ordinal_type iRow = (dev.league_rank() * dev.team_size() + dev.team_rank()) * rows_per_thread + loop;
      if (iRow >= m_A.numRows()) {
        return;
      }

      // mfh 20 Mar 2015, 07 Jun 2016: This is ordinal_type because it
      // needs to have the same type as n.
      ordinal_type kk = 0;

#ifdef KOKKOS_FAST_COMPILE
      for (; kk + 4 <= n; kk += 4) {
        strip_mine<4>(dev, iRow, kk);
      }
      for (; kk < n; ++kk) {
        strip_mine<1>(dev, iRow, kk);
      }
#else
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
      if ((n > 8) && (n % 8 == 1)) {
        strip_mine<9>(dev, iRow, kk);
        kk += 9;
      }
      for (; kk + 8 <= n; kk += 8) strip_mine<8>(dev, iRow, kk);
      if (kk < n) {
        switch (n - kk) {
#else   // NOT a GPU
      if ((n > 16) && (n % 16 == 1)) {
        strip_mine<17>(dev, iRow, kk);
        kk += 17;
      }

      for (; kk + 16 <= n; kk += 16) {
        strip_mine<16>(dev, iRow, kk);
      }

      if (kk < n) {
        switch (n - kk) {
          case 15: strip_mine<15>(dev, iRow, kk); break;

          case 14: strip_mine<14>(dev, iRow, kk); break;

          case 13: strip_mine<13>(dev, iRow, kk); break;

          case 12: strip_mine<12>(dev, iRow, kk); break;

          case 11: strip_mine<11>(dev, iRow, kk); break;

          case 10: strip_mine<10>(dev, iRow, kk); break;

          case 9: strip_mine<9>(dev, iRow, kk); break;

          case 8: strip_mine<8>(dev, iRow, kk); break;
#endif  // if/else: __CUDA_ARCH__ or __HIP_DEVICE_COMPILE__
          case 7: strip_mine<7>(dev, iRow, kk); break;

          case 6: strip_mine<6>(dev, iRow, kk); break;

          case 5: strip_mine<5>(dev, iRow, kk); break;

          case 4: strip_mine<4>(dev, iRow, kk); break;

          case 3: strip_mine<3>(dev, iRow, kk); break;

          case 2: strip_mine<2>(dev, iRow, kk); break;

          case 1: strip_mine_1(dev, iRow); break;
        }
      }
#endif  // KOKKOS_FAST_COMPILE
    }
  }
};

// spmv_alpha_beta_mv_no_transpose: version for CPU execution spaces
// (RangePolicy)
template <class execution_space, class AMatrix, class XVector, class YVector, int doalpha, int dobeta, bool conjugate,
          typename std::enable_if<!KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>()>::type* = nullptr>
static void spmv_alpha_beta_mv_no_transpose(const execution_space& exec,
                                            const typename YVector::non_const_value_type& alpha, const AMatrix& A,
                                            const XVector& x, const typename YVector::non_const_value_type& beta,
                                            const YVector& y) {
  using ordinal_type = typename AMatrix::non_const_ordinal_type;

  if (A.numRows() <= static_cast<ordinal_type>(0)) {
    return;
  }
  if (doalpha == 0) {
    if (dobeta != 1) {
      KokkosBlas::scal(exec, y, beta, y);
    }
    return;
  } else {
    // Assuming that no row contains duplicate entries, NNZPerRow
    // cannot be more than the number of columns of the matrix.  Thus,
    // the appropriate type is ordinal_type.
    const ordinal_type NNZPerRow = A.nnz() / A.numRows();
    const int vector_length      = 1;

#ifndef KOKKOS_FAST_COMPILE  // This uses templated functions on doalpha and
                             // dobeta and will produce 16 kernels

    typedef SPMV_MV_LayoutLeft_Functor<execution_space, AMatrix, XVector, YVector, doalpha, dobeta, conjugate> OpType;
    OpType op(alpha, A, x, beta, y, RowsPerThread<execution_space>(NNZPerRow), vector_length);

    typename AMatrix::const_ordinal_type nrow = A.numRows();

    Kokkos::parallel_for("KokkosSparse::spmv<MV,NoTranspose>", Kokkos::RangePolicy<execution_space>(exec, 0, nrow), op);

#else   // KOKKOS_FAST_COMPILE this will only instantiate one Kernel for
        // alpha/beta

    typedef SPMV_MV_LayoutLeft_Functor<execution_space, AMatrix, XVector, YVector, 2, 2, conjugate> OpType;

    typename AMatrix::const_ordinal_type nrow = A.numRows();

    OpType op(alpha, A, x, beta, y, RowsPerThread<execution_space>(NNZPerRow), vector_length);

    Kokkos::parallel_for("KokkosSparse::spmv<MV,NoTranspose>", Kokkos::RangePolicy<execution_space>(exec, 0, nrow), op);
#endif  // KOKKOS_FAST_COMPILE
  }
}

// spmv_alpha_beta_mv_no_transpose: version for GPU execution spaces
// (TeamPolicy)
template <class execution_space, class AMatrix, class XVector, class YVector, int doalpha, int dobeta, bool conjugate,
          typename std::enable_if<KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>()>::type* = nullptr>
static void spmv_alpha_beta_mv_no_transpose(const execution_space& exec,
                                            const typename YVector::non_const_value_type& alpha, const AMatrix& A,
                                            const XVector& x, const typename YVector::non_const_value_type& beta,
                                            const YVector& y) {
  using ordinal_type = typename AMatrix::non_const_ordinal_type;
  using size_type    = typename AMatrix::non_const_size_type;

  if (A.numRows() <= static_cast<ordinal_type>(0)) {
    return;
  }
  if (doalpha == 0) {
    if (dobeta != 1) {
      KokkosBlas::scal(exec, y, beta, y);
    }
    return;
  } else {
    // Assuming that no row contains duplicate entries, NNZPerRow
    // cannot be more than the number of columns of the matrix.  Thus,
    // the appropriate type is ordinal_type.
    const ordinal_type NNZPerRow = A.nnz() / A.numRows();

    ordinal_type vector_length = 1;
    while ((vector_length * 2 * 3 <= NNZPerRow) && (vector_length < 8)) vector_length *= 2;

#ifndef KOKKOS_FAST_COMPILE  // This uses templated functions on doalpha and
                             // dobeta and will produce 16 kernels

    typedef SPMV_MV_LayoutLeft_Functor<execution_space, AMatrix, XVector, YVector, doalpha, dobeta, conjugate> OpType;
    OpType op(alpha, A, x, beta, y, RowsPerThread<execution_space>(NNZPerRow), vector_length);

    typename AMatrix::const_ordinal_type nrow = A.numRows();

    const ordinal_type rows_per_thread = RowsPerThread<execution_space>(NNZPerRow);
    const ordinal_type team_size =
        Kokkos::TeamPolicy<execution_space>(exec, rows_per_thread, Kokkos::AUTO, vector_length)
            .team_size_recommended(op, Kokkos::ParallelForTag());
    const ordinal_type rows_per_team = rows_per_thread * team_size;
    const size_type nteams           = (nrow + rows_per_team - 1) / rows_per_team;
    Kokkos::parallel_for("KokkosSparse::spmv<MV,NoTranspose>",
                         Kokkos::TeamPolicy<execution_space>(exec, nteams, team_size, vector_length), op);

#else   // KOKKOS_FAST_COMPILE this will only instantiate one Kernel for
        // alpha/beta

    typedef SPMV_MV_LayoutLeft_Functor<execution_space, AMatrix, XVector, YVector, 2, 2, conjugate> OpType;

    typename AMatrix::const_ordinal_type nrow = A.numRows();

    OpType op(alpha, A, x, beta, y, RowsPerThread<execution_space>(NNZPerRow), vector_length);

    const ordinal_type rows_per_thread = RowsPerThread<execution_space>(NNZPerRow);
    const ordinal_type team_size =
        Kokkos::TeamPolicy<execution_space>(exec, rows_per_thread, Kokkos::AUTO, vector_length)
            .team_size_recommended(op, Kokkos::ParallelForTag());
    const ordinal_type rows_per_team = rows_per_thread * team_size;
    const size_type nteams           = (nrow + rows_per_team - 1) / rows_per_team;
    Kokkos::parallel_for("KokkosSparse::spmv<MV,NoTranspose>",
                         Kokkos::TeamPolicy<execution_space>(exec, nteams, team_size, vector_length), op);
#endif  // KOKKOS_FAST_COMPILE
  }
}

// spmv_alpha_beta_mv_transpose: version for CPU execution spaces (RangePolicy)
template <class execution_space, class AMatrix, class XVector, class YVector, int doalpha, int dobeta, bool conjugate,
          typename std::enable_if<!KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>()>::type* = nullptr>
static void spmv_alpha_beta_mv_transpose(const execution_space& exec,
                                         const typename YVector::non_const_value_type& alpha, const AMatrix& A,
                                         const XVector& x, const typename YVector::non_const_value_type& beta,
                                         const YVector& y) {
  using ordinal_type = typename AMatrix::non_const_ordinal_type;

  if (A.numRows() <= static_cast<ordinal_type>(0)) {
    return;
  }

  // We need to scale y first ("scaling" by zero just means filling
  // with zeros), since the functor works by atomic-adding into y.
  if (dobeta != 1) {
    KokkosBlas::scal(exec, y, beta, y);
  }

  if (doalpha != 0) {
#ifndef KOKKOS_FAST_COMPILE  // This uses templated functions on doalpha and
                             // dobeta and will produce 16 kernels

    typedef SPMV_MV_Transpose_Functor<execution_space, AMatrix, XVector, YVector, doalpha, dobeta, conjugate> OpType;
    OpType op(alpha, A, x, beta, y);

    const ordinal_type nrow = A.numRows();
    Kokkos::parallel_for("KokkosSparse::spmv<MV,Transpose>", Kokkos::RangePolicy<execution_space>(exec, 0, nrow), op);

#else  // KOKKOS_FAST_COMPILE this will only instantiate one Kernel for
       // alpha/beta

    typedef SPMV_MV_Transpose_Functor<execution_space, AMatrix, XVector, YVector, 2, 2, conjugate, SizeType> OpType;

    typename AMatrix::const_ordinal_type nrow = A.numRows();
    Kokkos::parallel_for("KokkosSparse::spmv<MV,Transpose>", Kokkos::RangePolicy<execution_space>(exec, 0, nrow), op);

#endif  // KOKKOS_FAST_COMPILE
  }
}

// spmv_alpha_beta_mv_transpose: version for GPU execution spaces (TeamPolicy)
template <class execution_space, class AMatrix, class XVector, class YVector, int doalpha, int dobeta, bool conjugate,
          typename std::enable_if<KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>()>::type* = nullptr>
static void spmv_alpha_beta_mv_transpose(const execution_space& exec,
                                         const typename YVector::non_const_value_type& alpha, const AMatrix& A,
                                         const XVector& x, const typename YVector::non_const_value_type& beta,
                                         const YVector& y) {
  using ordinal_type = typename AMatrix::non_const_ordinal_type;
  using size_type    = typename AMatrix::non_const_size_type;

  if (A.numRows() <= static_cast<ordinal_type>(0)) {
    return;
  }

  // We need to scale y first ("scaling" by zero just means filling
  // with zeros), since the functor works by atomic-adding into y.
  if (dobeta != 1) {
    KokkosBlas::scal(exec, y, beta, y);
  }

  if (doalpha != 0) {
    // Assuming that no row contains duplicate entries, NNZPerRow
    // cannot be more than the number of columns of the matrix.  Thus,
    // the appropriate type is ordinal_type.
    const ordinal_type NNZPerRow = static_cast<ordinal_type>(A.nnz() / A.numRows());

    ordinal_type vector_length = 1;
    // Transpose functor uses atomics which can't be vectorized on CPU
    while ((vector_length * 2 * 3 <= NNZPerRow) && (vector_length < 8)) vector_length *= 2;

#ifndef KOKKOS_FAST_COMPILE  // This uses templated functions on doalpha and
                             // dobeta and will produce 16 kernels

    typedef SPMV_MV_Transpose_Functor<execution_space, AMatrix, XVector, YVector, doalpha, dobeta, conjugate> OpType;
    OpType op(alpha, A, x, beta, y);

    const ordinal_type nrow            = A.numRows();
    const ordinal_type rows_per_thread = RowsPerThread<execution_space>(NNZPerRow);
    const ordinal_type team_size =
        Kokkos::TeamPolicy<execution_space>(exec, rows_per_thread, Kokkos::AUTO, vector_length)
            .team_size_recommended(op, Kokkos::ParallelForTag());
    const ordinal_type rows_per_team = rows_per_thread * team_size;
    op.rows_per_team                 = rows_per_team;
    const size_type nteams           = (nrow + rows_per_team - 1) / rows_per_team;
    Kokkos::parallel_for("KokkosSparse::spmv<MV,Transpose>",
                         Kokkos::TeamPolicy<execution_space>(exec, nteams, team_size, vector_length), op);

#else  // KOKKOS_FAST_COMPILE this will only instantiate one Kernel for
       // alpha/beta

    typedef SPMV_MV_Transpose_Functor<execution_space, AMatrix, XVector, YVector, 2, 2, conjugate, SizeType> OpType;

    typename AMatrix::const_ordinal_type nrow = A.numRows();
    OpType op(alpha, A, x, beta, y);

    const ordinal_type rows_per_thread = RowsPerThread<execution_space>(NNZPerRow);
    const ordinal_type team_size =
        Kokkos::TeamPolicy<execution_space>(exec, rows_per_thread, Kokkos::AUTO, vector_length)
            .team_size_recommended(op, Kokkos::ParallelForTag());
    const ordinal_type rows_per_team = rows_per_thread * team_size;
    op.rows_per_team                 = rows_per_team;
    const size_type nteams           = (nrow + rows_per_team - 1) / rows_per_team;
    Kokkos::parallel_for("KokkosSparse::spmv<MV,Transpose>",
                         Kokkos::TeamPolicy<execution_space>(exec, nteams, team_size, vector_length), op);

#endif  // KOKKOS_FAST_COMPILE
  }
}

template <class execution_space, class AMatrix, class XVector, class YVector, int doalpha, int dobeta>
static void spmv_alpha_beta_mv(const execution_space& exec, const char mode[],
                               const typename YVector::non_const_value_type& alpha, const AMatrix& A, const XVector& x,
                               const typename YVector::non_const_value_type& beta, const YVector& y) {
  if (mode[0] == NoTranspose[0]) {
    spmv_alpha_beta_mv_no_transpose<execution_space, AMatrix, XVector, YVector, doalpha, dobeta, false>(exec, alpha, A,
                                                                                                        x, beta, y);
  } else if (mode[0] == Conjugate[0]) {
    spmv_alpha_beta_mv_no_transpose<execution_space, AMatrix, XVector, YVector, doalpha, dobeta, true>(exec, alpha, A,
                                                                                                       x, beta, y);
  } else if (mode[0] == Transpose[0]) {
    spmv_alpha_beta_mv_transpose<execution_space, AMatrix, XVector, YVector, doalpha, dobeta, false>(exec, alpha, A, x,
                                                                                                     beta, y);
  } else if (mode[0] == ConjugateTranspose[0]) {
    spmv_alpha_beta_mv_transpose<execution_space, AMatrix, XVector, YVector, doalpha, dobeta, true>(exec, alpha, A, x,
                                                                                                    beta, y);
  } else {
    std::stringstream ss;
    ss << __FILE__ << ":" << __LINE__ << " Invalid transpose mode " << mode << " for KokkosSparse::spmv()";
    KokkosKernels::Impl::throw_runtime_exception(ss.str());
  }
}

template <class execution_space, class AMatrix, class XVector, class YVector, int doalpha>
void spmv_alpha_mv(const execution_space& exec, const char mode[], const typename YVector::non_const_value_type& alpha,
                   const AMatrix& A, const XVector& x, const typename YVector::non_const_value_type& beta,
                   const YVector& y) {
  typedef typename YVector::non_const_value_type coefficient_type;
  typedef Kokkos::ArithTraits<coefficient_type> KAT;

  if (beta == KAT::zero()) {
    spmv_alpha_beta_mv<execution_space, AMatrix, XVector, YVector, doalpha, 0>(exec, mode, alpha, A, x, beta, y);
  } else if (beta == KAT::one()) {
    spmv_alpha_beta_mv<execution_space, AMatrix, XVector, YVector, doalpha, 1>(exec, mode, alpha, A, x, beta, y);
  } else if (beta == -KAT::one()) {
    spmv_alpha_beta_mv<execution_space, AMatrix, XVector, YVector, doalpha, -1>(exec, mode, alpha, A, x, beta, y);
  } else {
    spmv_alpha_beta_mv<execution_space, AMatrix, XVector, YVector, doalpha, 2>(exec, mode, alpha, A, x, beta, y);
  }
}

}  // namespace Impl
}  // namespace KokkosSparse

#endif  // KOKKOSSPARSE_IMPL_SPMV_DEF_HPP_
