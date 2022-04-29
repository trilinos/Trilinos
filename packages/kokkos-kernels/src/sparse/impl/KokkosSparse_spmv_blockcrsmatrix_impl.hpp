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

#ifndef KOKKOSKERNELS_KOKKOSSPARSE_SPMV_BLOCKCRSMATRIX_IMPL_HPP
#define KOKKOSKERNELS_KOKKOSSPARSE_SPMV_BLOCKCRSMATRIX_IMPL_HPP

#include "KokkosBlas.hpp"
#include "KokkosBatched_Gemv_Serial_Internal.hpp"
#include "KokkosBatched_Gemm_Serial_Internal.hpp"
#include "KokkosKernels_ExecSpaceUtils.hpp"
#include "KokkosSparse_spmv_impl.hpp"

namespace KokkosSparse {
namespace Experimental {
namespace Impl {
namespace BCRS {

template <class AMatrix, class XVector, class YVector>
struct BCRS_GEMV_Functor {
  typedef typename AMatrix::execution_space execution_space;
  typedef typename AMatrix::non_const_value_type value_type;
  typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type team_member;
  typedef Kokkos::Details::ArithTraits<value_type> ATV;

  //! Nonconst version of the type of column indices in the sparse matrix.
  typedef typename AMatrix::non_const_ordinal_type ordinal_type;
  //! Nonconst version of the type of row offsets in the sparse matrix.
  typedef typename AMatrix::non_const_size_type size_type;

  const value_type alpha;
  AMatrix m_A;
  XVector m_x;
  YVector m_y;

  const ordinal_type block_dim;
  const ordinal_type blocks_per_team;

  bool conjugate = false;

  BCRS_GEMV_Functor(const value_type alpha_, const AMatrix m_A_,
                    const XVector m_x_, const YVector m_y_,
                    const int blocks_per_team_, bool conj_)
      : alpha(alpha_),
        m_A(m_A_),
        m_x(m_x_),
        m_y(m_y_),
        block_dim(m_A_.blockDim()),
        blocks_per_team(blocks_per_team_),
        conjugate(conj_) {
    static_assert(static_cast<int>(XVector::rank) == 1,
                  "XVector must be a rank 1 View.");
    static_assert(static_cast<int>(YVector::rank) == 1,
                  "YVector must be a rank 1 View.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const ordinal_type iBlock) const {
    const auto ystart = iBlock * block_dim;
    const auto start  = m_A.graph.row_map(iBlock);
    const ordinal_type count =
        static_cast<ordinal_type>(m_A.graph.row_map(iBlock + 1) - start);
    const auto row   = m_A.block_row_Const(iBlock);
    const auto beta1 = static_cast<value_type>(1);
    //
    if (conjugate) {
      for (ordinal_type ic = 0; ic < count; ++ic) {
        const auto Aview  = row.block(ic);
        const auto xstart = row.block_colidx(ic) * block_dim;
        for (ordinal_type ii = 0; ii < block_dim; ++ii) {
          value_type t(0);
          for (ordinal_type jj = 0; jj < block_dim; ++jj) {
            const auto aval =
                Kokkos::ArithTraits<value_type>::conj(Aview(ii, jj));
            t += aval * m_x(xstart + jj);
          }
          m_y(ystart + ii) += alpha * t;
        }
      }
    } else {
      for (ordinal_type ic = 0; ic < count; ++ic) {
        const auto Aview  = row.block(ic);
        const auto xstart = row.block_colidx(ic) * block_dim;
        KokkosBatched::SerialGemvInternal<KokkosBatched::Algo::Gemv::Blocked>::
            invoke<value_type, value_type>(
                block_dim, block_dim, alpha, Aview.data(), Aview.stride_0(),
                Aview.stride_1(), &m_x(xstart), m_x.stride_0(), beta1,
                &m_y(ystart), m_y.stride_0());
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member &dev) const {
    using y_value_type = typename YVector::non_const_value_type;
    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(dev, 0, blocks_per_team),
        [&](const ordinal_type &loop) {
          const ordinal_type iBlock =
              static_cast<ordinal_type>(dev.league_rank()) * blocks_per_team +
              loop;
          if (iBlock >= m_A.numRows()) {
            return;
          }
          const auto start = m_A.graph.row_map(iBlock);
          const ordinal_type count =
              static_cast<ordinal_type>(m_A.graph.row_map(iBlock + 1) - start);
          const auto row = m_A.block_row_Const(iBlock);
          //
          auto yview = Kokkos::subview(
              m_y, Kokkos::make_pair(iBlock * block_dim,
                                     iBlock * block_dim + block_dim));
          //
          for (ordinal_type ir = 0; ir < block_dim; ++ir) {
            y_value_type sum = 0;

            Kokkos::parallel_reduce(
                Kokkos::ThreadVectorRange(dev, count),
                [&](const ordinal_type &iEntry, y_value_type &lsum) {
                  const auto start_col = row.block_colidx(iEntry) * block_dim;
                  for (ordinal_type jr = 0; jr < block_dim; ++jr) {
                    const value_type val =
                        conjugate
                            ? ATV::conj(row.local_block_value(iEntry, ir, jr))
                            : row.local_block_value(iEntry, ir, jr);
                    lsum += val * m_x(start_col + jr);
                  }
                },
                sum);

            Kokkos::single(Kokkos::PerThread(dev), [&]() {
              sum *= alpha;
              yview(ir) += sum;
            });
          }
        });
  }
};

/* ******************* */

//
// spMatVec_no_transpose: version for CPU execution spaces
// (RangePolicy or trivial serial impl used)
//
template <class AT, class AO, class AD, class AS, class AlphaType,
          class XVector, class BetaType, class YVector,
          typename std::enable_if<!KokkosKernels::Impl::kk_is_gpu_exec_space<
              typename YVector::execution_space>()>::type * = nullptr>
void spMatVec_no_transpose(
    const KokkosKernels::Experimental::Controls &controls,
    const AlphaType &alpha,
    const KokkosSparse::Experimental::BlockCrsMatrix<
        AT, AO, AD, Kokkos::MemoryTraits<Kokkos::Unmanaged>, AS> &A,
    const XVector &x, const BetaType &beta, YVector &y, bool useConjugate) {
  // This is required to maintain semantics of KokkosKernels native SpMV:
  // if y contains NaN but beta = 0, the result y should be filled with 0.
  // For example, this is useful for passing in uninitialized y and beta=0.
  if (beta == Kokkos::ArithTraits<BetaType>::zero())
    Kokkos::deep_copy(y, Kokkos::ArithTraits<BetaType>::zero());
  else
    KokkosBlas::scal(y, beta, y);

  //
  // Treat the case y <- alpha * A * x + beta * y
  //

  typedef KokkosSparse::Experimental::BlockCrsMatrix<
      AT, AO, AD, Kokkos::MemoryTraits<Kokkos::Unmanaged>, AS>
      AMatrix_Internal;

  bool use_dynamic_schedule = false;  // Forces the use of a dynamic schedule
  bool use_static_schedule  = false;  // Forces the use of a static schedule
  if (controls.isParameter("schedule")) {
    if (controls.getParameter("schedule") == "dynamic") {
      use_dynamic_schedule = true;
    } else if (controls.getParameter("schedule") == "static") {
      use_static_schedule = true;
    }
  }

  BCRS_GEMV_Functor<AMatrix_Internal, XVector, YVector> func(alpha, A, x, y, 1,
                                                             useConjugate);
  if (((A.nnz() > 10000000) || use_dynamic_schedule) && !use_static_schedule) {
    Kokkos::parallel_for(
        "KokkosSparse::bcrs_spmv<NoTranspose,Dynamic>",
        Kokkos::RangePolicy<
            typename AMatrix_Internal::device_type::execution_space,
            Kokkos::Schedule<Kokkos::Dynamic>>(0, A.numRows()),
        func);
  } else {
    Kokkos::parallel_for(
        "KokkosSparse::bcrs_spmv<NoTranspose,Static>",
        Kokkos::RangePolicy<
            typename AMatrix_Internal::device_type::execution_space,
            Kokkos::Schedule<Kokkos::Static>>(0, A.numRows()),
        func);
  }
}

/* ******************* */

//
// spMatVec_no_transpose: version for GPU execution spaces (TeamPolicy used)
//
template <class AT, class AO, class AD, class AS, class AlphaType,
          class XVector, class BetaType, class YVector,
          typename std::enable_if<KokkosKernels::Impl::kk_is_gpu_exec_space<
              typename YVector::execution_space>()>::type * = nullptr>
void spMatVec_no_transpose(
    const KokkosKernels::Experimental::Controls &controls,
    const AlphaType &alpha,
    const KokkosSparse::Experimental::BlockCrsMatrix<
        AT, AO, AD, Kokkos::MemoryTraits<Kokkos::Unmanaged>, AS> &A,
    const XVector &x, const BetaType &beta, YVector &y, bool useConjugate) {
  if (A.numRows() <= static_cast<AO>(0)) {
    return;
  }

  // We need to scale y first ("scaling" by zero just means filling
  // with zeros), since the functor updates y (by adding alpha Op(A) x).
  KokkosBlas::scal(y, beta, y);

  typedef KokkosSparse::Experimental::BlockCrsMatrix<
      AT, AO, AD, Kokkos::MemoryTraits<Kokkos::Unmanaged>, AS>
      AMatrix_Internal;
  typedef typename AMatrix_Internal::execution_space execution_space;

  bool use_dynamic_schedule = false;  // Forces the use of a dynamic schedule
  bool use_static_schedule  = false;  // Forces the use of a static schedule
  if (controls.isParameter("schedule")) {
    if (controls.getParameter("schedule") == "dynamic") {
      use_dynamic_schedule = true;
    } else if (controls.getParameter("schedule") == "static") {
      use_static_schedule = true;
    }
  }
  int team_size             = -1;
  int vector_length         = -1;
  int64_t blocks_per_thread = -1;

  //
  // Use the controls to allow the user to pass in some tuning parameters.
  //
  if (controls.isParameter("team size")) {
    team_size = std::stoi(controls.getParameter("team size"));
  }
  if (controls.isParameter("vector length")) {
    vector_length = std::stoi(controls.getParameter("vector length"));
  }
  if (controls.isParameter("rows per thread")) {
    blocks_per_thread = std::stoll(controls.getParameter("rows per thread"));
  }

  //
  // Use the existing launch parameters routine from SPMV
  //
  int64_t blocks_per_team =
      KokkosSparse::Impl::spmv_launch_parameters<execution_space>(
          A.numRows(), A.nnz(), blocks_per_thread, team_size, vector_length);
  int64_t worksets = (A.numRows() + blocks_per_team - 1) / blocks_per_team;

  AMatrix_Internal A_internal = A;

  BCRS_GEMV_Functor<AMatrix_Internal, XVector, YVector> func(
      alpha, A_internal, x, y, blocks_per_team, useConjugate);

  if (((A.nnz() > 10000000) || use_dynamic_schedule) && !use_static_schedule) {
    Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>>
        policy(1, 1);
    if (team_size < 0)
      policy = Kokkos::TeamPolicy<execution_space,
                                  Kokkos::Schedule<Kokkos::Dynamic>>(
          worksets, Kokkos::AUTO, vector_length);
    else
      policy = Kokkos::TeamPolicy<execution_space,
                                  Kokkos::Schedule<Kokkos::Dynamic>>(
          worksets, team_size, vector_length);
    Kokkos::parallel_for("KokkosSparse::bcrs_spmv<NoTranspose,Dynamic>", policy,
                         func);
  } else {
    Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>>
        policy(1, 1);
    if (team_size < 0)
      policy =
          Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>>(
              worksets, Kokkos::AUTO, vector_length);
    else
      policy =
          Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>>(
              worksets, team_size, vector_length);
    Kokkos::parallel_for("KokkosSparse::bcrs_spmv<NoTranspose, Static>", policy,
                         func);
  }
}

/* ******************* */

template <class AMatrix, class XVector, class YVector>
struct BCRS_GEMV_Transpose_Functor {
  typedef typename AMatrix::execution_space execution_space;
  typedef typename AMatrix::non_const_value_type value_type;
  typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type team_member;
  typedef Kokkos::Details::ArithTraits<value_type> ATV;

  //! Nonconst version of the type of column indices in the sparse matrix.
  typedef typename AMatrix::non_const_ordinal_type ordinal_type;
  //! Nonconst version of the type of row offsets in the sparse matrix.
  typedef typename AMatrix::non_const_size_type size_type;

  const value_type alpha;

  AMatrix m_A;
  XVector m_x;
  YVector m_y;

  const ordinal_type block_dim;
  const ordinal_type blocks_per_team;

  bool conjugate = false;

  BCRS_GEMV_Transpose_Functor(const value_type alpha_, const AMatrix m_A_,
                              const XVector m_x_, const YVector m_y_,
                              const int blocks_per_team_, bool conj_)
      : alpha(alpha_),
        m_A(m_A_),
        m_x(m_x_),
        m_y(m_y_),
        block_dim(m_A_.blockDim()),
        blocks_per_team(blocks_per_team_),
        conjugate(conj_) {
    static_assert(static_cast<int>(XVector::rank) == 1,
                  "XVector must be a rank 1 View.");
    static_assert(static_cast<int>(YVector::rank) == 1,
                  "YVector must be a rank 1 View.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const ordinal_type iBlock) const {
    //
    // Assume that alpha is not zero
    //
    const auto xstart = iBlock * block_dim;
    const auto xview =
        Kokkos::subview(m_x, Kokkos::make_pair(xstart, xstart + block_dim));
    const auto start = m_A.graph.row_map(iBlock);
    const ordinal_type count =
        static_cast<ordinal_type>(m_A.graph.row_map(iBlock + 1) - start);
    const auto row    = m_A.block_row_Const(iBlock);
    const auto beta1  = static_cast<value_type>(1);
    const auto alpha1 = beta1;
    if (conjugate) {
      for (ordinal_type ic = 0; ic < count; ++ic) {
        const auto Aview  = row.block(ic);
        const auto ystart = row.block_colidx(ic) * block_dim;
        for (ordinal_type jj = 0; jj < block_dim; ++jj) {
          value_type t(0);
          for (ordinal_type ii = 0; ii < block_dim; ++ii) {
            const auto aval =
                Kokkos::ArithTraits<value_type>::conj(Aview(ii, jj));
            t += aval * xview(ii);
          }
          t *= alpha;
          Kokkos::atomic_add(&m_y(ystart + jj), t);
        }
      }
    } else {
      for (ordinal_type ic = 0; ic < count; ++ic) {
        const auto Aview  = row.block(ic);
        const auto ystart = row.block_colidx(ic) * block_dim;
        for (ordinal_type jj = 0; jj < block_dim; ++jj) {
          value_type t(0);
          KokkosBatched::SerialGemvInternal<
              KokkosBatched::Algo::Gemv::Blocked>::invoke<value_type,
                                                          value_type>(
              1, block_dim, alpha1, Aview.data() + jj, Aview.stride_1(),
              Aview.stride_0(), xview.data(), xview.stride_0(), beta1, &t, 1);
          t *= alpha;
          Kokkos::atomic_add(&m_y(ystart + jj), t);
        }
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member &dev) const {
    using y_value_type = typename YVector::non_const_value_type;
    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(dev, 0, blocks_per_team),
        [&](const ordinal_type &loop) {
          const ordinal_type iBlock =
              static_cast<ordinal_type>(dev.league_rank()) * blocks_per_team +
              loop;
          if (iBlock >= m_A.numRows()) {
            return;
          }
          const auto start = m_A.graph.row_map(iBlock);
          const ordinal_type count =
              static_cast<ordinal_type>(m_A.graph.row_map(iBlock + 1) - start);
          const auto row = m_A.block_row_Const(iBlock);
          //
          for (ordinal_type ir = 0; ir < block_dim; ++ir) {
            Kokkos::parallel_for(
                Kokkos::ThreadVectorRange(dev, count),
                [&](const ordinal_type &iEntry) {
                  for (ordinal_type jr = 0; jr < block_dim; ++jr) {
                    const value_type val =
                        conjugate
                            ? ATV::conj(row.local_block_value(iEntry, jr, ir))
                            : row.local_block_value(iEntry, jr, ir);
                    const ordinal_type ind = row.block_colidx(iEntry);
                    Kokkos::atomic_add(
                        &m_y(block_dim * ind + ir),
                        static_cast<y_value_type>(
                            alpha * val * m_x(block_dim * iBlock + jr)));
                  }
                });
          }
        });
  }
};

/* ******************* */

/// \brief  spMatVec_transpose: version for CPU execution spaces (RangePolicy or
/// trivial serial impl used)
template <class AT, class AO, class AD, class AS, class AlphaType,
          class XVector, class BetaType, class YVector,
          typename std::enable_if<!KokkosKernels::Impl::kk_is_gpu_exec_space<
              typename YVector::execution_space>()>::type * = nullptr>
void spMatVec_transpose(
    const KokkosKernels::Experimental::Controls &controls,
    const AlphaType &alpha,
    const KokkosSparse::Experimental::BlockCrsMatrix<
        AT, AO, AD, Kokkos::MemoryTraits<Kokkos::Unmanaged>, AS> &A,
    const XVector &x, const BetaType &beta, YVector &y, bool useConjugate) {
  // This is required to maintain semantics of KokkosKernels native SpMV:
  // if y contains NaN but beta = 0, the result y should be filled with 0.
  // For example, this is useful for passing in uninitialized y and beta=0.
  if (beta == Kokkos::ArithTraits<BetaType>::zero())
    Kokkos::deep_copy(y, Kokkos::ArithTraits<BetaType>::zero());
  else
    KokkosBlas::scal(y, beta, y);

  if (alpha == Kokkos::ArithTraits<AlphaType>::zero()) return;

  //
  // Treat the case y <- alpha * A^T * x + beta * y
  //

  typedef KokkosSparse::Experimental::BlockCrsMatrix<
      AT, AO, AD, Kokkos::MemoryTraits<Kokkos::Unmanaged>, AS>
      AMatrix_Internal;

  AMatrix_Internal A_internal = A;

  bool use_dynamic_schedule = false;  // Forces the use of a dynamic schedule
  bool use_static_schedule  = false;  // Forces the use of a static schedule
  if (controls.isParameter("schedule")) {
    if (controls.getParameter("schedule") == "dynamic") {
      use_dynamic_schedule = true;
    } else if (controls.getParameter("schedule") == "static") {
      use_static_schedule = true;
    }
  }

  BCRS_GEMV_Transpose_Functor<AMatrix_Internal, XVector, YVector> func(
      alpha, A_internal, x, y, 1, useConjugate);
  if (((A.nnz() > 10000000) || use_dynamic_schedule) && !use_static_schedule) {
    Kokkos::parallel_for(
        "KokkosSparse::blockcrs_spmv<Transpose,Dynamic>",
        Kokkos::RangePolicy<
            typename AMatrix_Internal::device_type::execution_space,
            Kokkos::Schedule<Kokkos::Dynamic>>(0, A.numRows()),
        func);
  } else {
    Kokkos::parallel_for(
        "KokkosSparse::blockcrs_spmv<Transpose,Static>",
        Kokkos::RangePolicy<
            typename AMatrix_Internal::device_type::execution_space,
            Kokkos::Schedule<Kokkos::Static>>(0, A.numRows()),
        func);
  }
}

//
// spMatVec_transpose: version for GPU execution spaces (TeamPolicy used)
//
template <class AMatrix, class AlphaType, class XVector, class BetaType,
          class YVector,
          typename std::enable_if<KokkosKernels::Impl::kk_is_gpu_exec_space<
              typename YVector::execution_space>()>::type * = nullptr>
void spMatVec_transpose(const KokkosKernels::Experimental::Controls &controls,
                        const AlphaType &alpha, const AMatrix &A,
                        const XVector &x, const BetaType &beta, YVector &y,
                        bool useConjugate) {
  if (A.numRows() <= 0) {
    return;
  }

  // We need to scale y first ("scaling" by zero just means filling
  // with zeros), since the functor works by atomic-adding into y.
  KokkosBlas::scal(y, beta, y);

  typedef typename AMatrix::execution_space execution_space;

  bool use_dynamic_schedule = false;  // Forces the use of a dynamic schedule
  bool use_static_schedule  = false;  // Forces the use of a static schedule
  if (controls.isParameter("schedule")) {
    if (controls.getParameter("schedule") == "dynamic") {
      use_dynamic_schedule = true;
    } else if (controls.getParameter("schedule") == "static") {
      use_static_schedule = true;
    }
  }
  int team_size             = -1;
  int vector_length         = -1;
  int64_t blocks_per_thread = -1;

  //
  // Use the controls to allow the user to pass in some tuning parameters.
  //
  if (controls.isParameter("team size")) {
    team_size = std::stoi(controls.getParameter("team size"));
  }
  if (controls.isParameter("vector length")) {
    vector_length = std::stoi(controls.getParameter("vector length"));
  }
  if (controls.isParameter("rows per thread")) {
    blocks_per_thread = std::stoll(controls.getParameter("rows per thread"));
  }

  //
  // Use the existing launch parameters routine from SPMV
  //
  int64_t blocks_per_team =
      KokkosSparse::Impl::spmv_launch_parameters<execution_space>(
          A.numRows(), A.nnz(), blocks_per_thread, team_size, vector_length);
  int64_t worksets = (A.numRows() + blocks_per_team - 1) / blocks_per_team;

  BCRS_GEMV_Transpose_Functor<AMatrix, XVector, YVector> func(
      alpha, A, x, y, blocks_per_team, useConjugate);

  if (((A.nnz() > 10000000) || use_dynamic_schedule) && !use_static_schedule) {
    Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>>
        policy(1, 1);
    if (team_size < 0)
      policy = Kokkos::TeamPolicy<execution_space,
                                  Kokkos::Schedule<Kokkos::Dynamic>>(
          worksets, Kokkos::AUTO, vector_length);
    else
      policy = Kokkos::TeamPolicy<execution_space,
                                  Kokkos::Schedule<Kokkos::Dynamic>>(
          worksets, team_size, vector_length);
    Kokkos::parallel_for("KokkosSparse::bspmv<Transpose,Dynamic>", policy,
                         func);
  } else {
    Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>>
        policy(1, 1);
    if (team_size < 0)
      policy =
          Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>>(
              worksets, Kokkos::AUTO, vector_length);
    else
      policy =
          Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>>(
              worksets, team_size, vector_length);
    Kokkos::parallel_for("KokkosSparse::bspmv<Transpose, Static>", policy,
                         func);
  }
}

/* ******************* */

template <class AMatrix, class XVector, class YVector>
struct BCRS_GEMM_Functor {
  typedef typename AMatrix::execution_space execution_space;
  typedef typename AMatrix::non_const_value_type value_type;
  typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type team_member;
  typedef Kokkos::Details::ArithTraits<value_type> ATV;

  //! Nonconst version of the type of column indices in the sparse matrix.
  typedef typename AMatrix::non_const_ordinal_type ordinal_type;
  //! Nonconst version of the type of row offsets in the sparse matrix.
  typedef typename AMatrix::non_const_size_type size_type;

  const value_type alpha;
  AMatrix m_A;
  XVector m_x;
  YVector m_y;
  const ordinal_type block_dim;
  const ordinal_type num_rhs;

  const ordinal_type blocks_per_team;

  bool conjugate = false;

  BCRS_GEMM_Functor(const value_type alpha_, const AMatrix m_A_,
                    const XVector m_x_, const YVector m_y_,
                    const int blocks_per_team_, bool conj_)
      : alpha(alpha_),
        m_A(m_A_),
        m_x(m_x_),
        m_y(m_y_),
        block_dim(m_A_.blockDim()),
        num_rhs(m_x_.extent(1)),
        blocks_per_team(blocks_per_team_),
        conjugate(conj_) {
    static_assert(static_cast<int>(XVector::rank) == 2,
                  "XVector must be a rank 2 View.");
    static_assert(static_cast<int>(YVector::rank) == 2,
                  "YVector must be a rank 2 View.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const ordinal_type iBlock) const {
    //
    const auto ystart = iBlock * block_dim;
    const auto start  = m_A.graph.row_map(iBlock);
    const ordinal_type count =
        static_cast<ordinal_type>(m_A.graph.row_map(iBlock + 1) - start);
    const auto row = m_A.block_row_Const(iBlock);
    //
    for (ordinal_type ic = 0; ic < count; ++ic) {
      const auto Aview  = row.block(ic);
      const auto xstart = row.block_colidx(ic) * block_dim;
      for (ordinal_type jr = 0; jr < num_rhs; ++jr) {
        for (ordinal_type ii = 0; ii < block_dim; ++ii) {
          value_type t(0);
          for (ordinal_type jj = 0; jj < block_dim; ++jj) {
            const auto aval =
                (conjugate)
                    ? Kokkos::ArithTraits<value_type>::conj(Aview(ii, jj))
                    : Aview(ii, jj);
            t += aval * m_x(xstart + jj, jr);
          }
          m_y(ystart + ii, jr) += alpha * t;
        }
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member &dev) const {
    using y_value_type = typename YVector::non_const_value_type;
    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(dev, 0, blocks_per_team),
        [&](const ordinal_type &loop) {
          const ordinal_type iBlock =
              static_cast<ordinal_type>(dev.league_rank()) * blocks_per_team +
              loop;
          if (iBlock >= m_A.numRows()) {
            return;
          }
          //
          const auto start = m_A.graph.row_map(iBlock);
          const ordinal_type count =
              static_cast<ordinal_type>(m_A.graph.row_map(iBlock + 1) - start);
          const auto row  = m_A.block_row_Const(iBlock);
          const auto nrhs = num_rhs;
          //
          for (ordinal_type ic = 0; ic < nrhs; ++ic) {
            for (ordinal_type ir = 0; ir < block_dim; ++ir) {
              y_value_type sum = 0;

              Kokkos::parallel_reduce(
                  Kokkos::ThreadVectorRange(dev, count),
                  [&](const ordinal_type &iEntry, y_value_type &lsum) {
                    const auto start_col = row.block_colidx(iEntry) * block_dim;
                    for (ordinal_type jr = 0; jr < block_dim; ++jr) {
                      const value_type val =
                          conjugate
                              ? ATV::conj(row.local_block_value(iEntry, ir, jr))
                              : row.local_block_value(iEntry, ir, jr);
                      lsum += val * m_x(start_col + jr, ic);
                    }
                  },
                  sum);

              Kokkos::single(Kokkos::PerThread(dev), [&]() {
                sum *= alpha;
                m_y(iBlock * block_dim + ir, ic) += sum;
              });
            }
          }
          //
        });
  }
};

/* ******************* */

//
// spMatMultiVec_no_transpose: version for CPU execution spaces
// (RangePolicy or trivial serial impl used)
//
template <class AT, class AO, class AD, class AS, class AlphaType,
          class XVector, class BetaType, class YVector,
          typename std::enable_if<!KokkosKernels::Impl::kk_is_gpu_exec_space<
              typename YVector::execution_space>()>::type * = nullptr>
void spMatMultiVec_no_transpose(
    const KokkosKernels::Experimental::Controls &controls,
    const AlphaType &alpha,
    const KokkosSparse::Experimental::BlockCrsMatrix<
        AT, AO, AD, Kokkos::MemoryTraits<Kokkos::Unmanaged>, AS> &A,
    const XVector &x, const BetaType &beta, YVector &y, bool useConjugate) {
  // This is required to maintain semantics of KokkosKernels native SpMV:
  // if y contains NaN but beta = 0, the result y should be filled with 0.
  // For example, this is useful for passing in uninitialized y and beta=0.
  if (beta == Kokkos::ArithTraits<BetaType>::zero())
    Kokkos::deep_copy(y, Kokkos::ArithTraits<BetaType>::zero());
  else
    KokkosBlas::scal(y, beta, y);
  //
  // Treat the case y <- alpha * A * x + beta * y
  //
  typedef KokkosSparse::Experimental::BlockCrsMatrix<
      AT, AO, AD, Kokkos::MemoryTraits<Kokkos::Unmanaged>, AS>
      AMatrix_Internal;

  bool use_dynamic_schedule = false;  // Forces the use of a dynamic schedule
  bool use_static_schedule  = false;  // Forces the use of a static schedule
  if (controls.isParameter("schedule")) {
    if (controls.getParameter("schedule") == "dynamic") {
      use_dynamic_schedule = true;
    } else if (controls.getParameter("schedule") == "static") {
      use_static_schedule = true;
    }
  }

  BCRS_GEMM_Functor<AMatrix_Internal, XVector, YVector> func(alpha, A, x, y, 1,
                                                             useConjugate);
  if (((A.nnz() > 10000000) || use_dynamic_schedule) && !use_static_schedule) {
    Kokkos::parallel_for(
        "KokkosSparse::bcrs_spm_mv<NoTranspose,Dynamic>",
        Kokkos::RangePolicy<
            typename AMatrix_Internal::device_type::execution_space,
            Kokkos::Schedule<Kokkos::Dynamic>>(0, A.numRows()),
        func);
  } else {
    Kokkos::parallel_for(
        "KokkosSparse::bcrs_spm_mv<NoTranspose,Static>",
        Kokkos::RangePolicy<
            typename AMatrix_Internal::device_type::execution_space,
            Kokkos::Schedule<Kokkos::Static>>(0, A.numRows()),
        func);
  }
}

/* ******************* */

//
// spMatMultiVec_no_transpose: version for GPU execution spaces (TeamPolicy
// used)
//
template <class AT, class AO, class AD, class AS, class AlphaType,
          class XVector, class BetaType, class YVector,
          typename std::enable_if<KokkosKernels::Impl::kk_is_gpu_exec_space<
              typename YVector::execution_space>()>::type * = nullptr>
void spMatMultiVec_no_transpose(
    const KokkosKernels::Experimental::Controls &controls,
    const AlphaType &alpha,
    const KokkosSparse::Experimental::BlockCrsMatrix<
        AT, AO, AD, Kokkos::MemoryTraits<Kokkos::Unmanaged>, AS> &A,
    const XVector &x, const BetaType &beta, YVector &y, bool useConjugate) {
  if (A.numRows() <= static_cast<AO>(0)) {
    return;
  }

  KokkosBlas::scal(y, beta, y);

  typedef KokkosSparse::Experimental::BlockCrsMatrix<
      AT, AO, AD, Kokkos::MemoryTraits<Kokkos::Unmanaged>, AS>
      AMatrix_Internal;
  typedef typename AMatrix_Internal::execution_space execution_space;

  bool use_dynamic_schedule = false;  // Forces the use of a dynamic schedule
  bool use_static_schedule  = false;  // Forces the use of a static schedule
  if (controls.isParameter("schedule")) {
    if (controls.getParameter("schedule") == "dynamic") {
      use_dynamic_schedule = true;
    } else if (controls.getParameter("schedule") == "static") {
      use_static_schedule = true;
    }
  }
  int team_size             = -1;
  int vector_length         = -1;
  int64_t blocks_per_thread = -1;

  //
  // Use the controls to allow the user to pass in some tuning parameters.
  //
  if (controls.isParameter("team size")) {
    team_size = std::stoi(controls.getParameter("team size"));
  }
  if (controls.isParameter("vector length")) {
    vector_length = std::stoi(controls.getParameter("vector length"));
  }
  if (controls.isParameter("rows per thread")) {
    blocks_per_thread = std::stoll(controls.getParameter("rows per thread"));
  }

  //
  // Use the existing launch parameters routine from SPMV
  //
  int64_t blocks_per_team =
      KokkosSparse::Impl::spmv_launch_parameters<execution_space>(
          A.numRows(), A.nnz(), blocks_per_thread, team_size, vector_length);
  int64_t worksets = (A.numRows() + blocks_per_team - 1) / blocks_per_team;

  AMatrix_Internal A_internal = A;

  BCRS_GEMM_Functor<AMatrix_Internal, XVector, YVector> func(
      alpha, A_internal, x, y, blocks_per_team, useConjugate);

  if (((A.nnz() > 10000000) || use_dynamic_schedule) && !use_static_schedule) {
    Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>>
        policy(1, 1);
    if (team_size < 0)
      policy = Kokkos::TeamPolicy<execution_space,
                                  Kokkos::Schedule<Kokkos::Dynamic>>(
          worksets, Kokkos::AUTO, vector_length);
    else
      policy = Kokkos::TeamPolicy<execution_space,
                                  Kokkos::Schedule<Kokkos::Dynamic>>(
          worksets, team_size, vector_length);
    Kokkos::parallel_for("KokkosSparse::bcrs_spm_mv<NoTranspose,Dynamic>",
                         policy, func);
  } else {
    Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>>
        policy(1, 1);
    if (team_size < 0)
      policy =
          Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>>(
              worksets, Kokkos::AUTO, vector_length);
    else
      policy =
          Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>>(
              worksets, team_size, vector_length);
    Kokkos::parallel_for("KokkosSparse::bcrs_spm_mv<NoTranspose, Static>",
                         policy, func);
  }
}

/* ******************* */

template <class AMatrix, class XVector, class YVector>
struct BCRS_GEMM_Transpose_Functor {
  typedef typename AMatrix::execution_space execution_space;
  typedef typename AMatrix::non_const_value_type value_type;
  typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type team_member;
  typedef Kokkos::Details::ArithTraits<value_type> ATV;

  //! Nonconst version of the type of column indices in the sparse matrix.
  typedef typename AMatrix::non_const_ordinal_type ordinal_type;
  //! Nonconst version of the type of row offsets in the sparse matrix.
  typedef typename AMatrix::non_const_size_type size_type;

  const value_type alpha;
  AMatrix m_A;
  XVector m_x;
  YVector m_y;
  const ordinal_type block_dim;
  const ordinal_type num_rhs;

  const ordinal_type blocks_per_team;

  bool conjugate = false;

  BCRS_GEMM_Transpose_Functor(const value_type alpha_, const AMatrix m_A_,
                              const XVector m_x_, const YVector m_y_,
                              const int blocks_per_team_, bool conj_)
      : alpha(alpha_),
        m_A(m_A_),
        m_x(m_x_),
        m_y(m_y_),
        block_dim(m_A_.blockDim()),
        num_rhs(m_x_.extent(1)),
        blocks_per_team(blocks_per_team_),
        conjugate(conj_) {
    static_assert(static_cast<int>(XVector::rank) == 2,
                  "XVector must be a rank 2 View.");
    static_assert(static_cast<int>(YVector::rank) == 2,
                  "YVector must be a rank 2 View.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const ordinal_type iBlock) const {
    //
    const auto xstart = iBlock * block_dim;
    const auto xview  = Kokkos::subview(
        m_x, Kokkos::make_pair(xstart, xstart + block_dim), Kokkos::ALL());
    const auto start = m_A.graph.row_map(iBlock);
    const ordinal_type count =
        static_cast<ordinal_type>(m_A.graph.row_map(iBlock + 1) - start);
    const auto row    = m_A.block_row_Const(iBlock);
    const auto beta1  = static_cast<value_type>(1);
    const auto alpha1 = beta1;
    const auto ldx    = m_x.stride_1();
    //
    if (conjugate) {
      for (ordinal_type ic = 0; ic < count; ++ic) {
        const auto Aview  = row.block(ic);
        const auto ystart = row.block_colidx(ic) * block_dim;
        for (ordinal_type jr = 0; jr < num_rhs; ++jr) {
          for (ordinal_type jj = 0; jj < block_dim; ++jj) {
            value_type t(0);
            for (ordinal_type ii = 0; ii < block_dim; ++ii) {
              const auto aval =
                  Kokkos::ArithTraits<value_type>::conj(Aview(ii, jj));
              t += aval * xview(ii, jr);
            }
            t *= alpha;
            Kokkos::atomic_add(&m_y(ystart + jj, jr), t);
          }
        }
      }
    } else {
      for (ordinal_type ic = 0; ic < count; ++ic) {
        const auto Aview  = row.block(ic);
        const auto ystart = row.block_colidx(ic) * block_dim;
        for (ordinal_type jr = 0; jr < num_rhs; ++jr) {
          for (ordinal_type jj = 0; jj < block_dim; ++jj) {
            value_type t(0);
            KokkosBatched::SerialGemvInternal<
                KokkosBatched::Algo::Gemv::Blocked>::invoke<value_type,
                                                            value_type>(
                1, block_dim, alpha1, Aview.data() + jj, Aview.stride_1(),
                Aview.stride_0(), xview.data() + jr * ldx, xview.stride_0(),
                beta1, &t, 1);
            t *= alpha;
            Kokkos::atomic_add(&m_y(ystart + jj, jr), t);
          }
        }
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member &dev) const {
    using y_value_type = typename YVector::non_const_value_type;
    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(dev, 0, blocks_per_team),
        [&](const ordinal_type &loop) {
          const ordinal_type iBlock =
              static_cast<ordinal_type>(dev.league_rank()) * blocks_per_team +
              loop;
          if (iBlock >= m_A.numRows()) {
            return;
          }
          //
          const auto start = m_A.graph.row_map(iBlock);
          const ordinal_type count =
              static_cast<ordinal_type>(m_A.graph.row_map(iBlock + 1) - start);
          const auto row  = m_A.block_row_Const(iBlock);
          const auto nrhs = m_x.extent(1);
          //
          for (size_t ic = 0; ic < nrhs; ++ic) {
            for (ordinal_type ir = 0; ir < block_dim; ++ir) {
              Kokkos::parallel_for(
                  Kokkos::ThreadVectorRange(dev, count),
                  [&](const ordinal_type &iEntry) {
                    for (ordinal_type jr = 0; jr < block_dim; ++jr) {
                      const value_type val =
                          conjugate
                              ? ATV::conj(row.local_block_value(iEntry, jr, ir))
                              : row.local_block_value(iEntry, jr, ir);
                      const ordinal_type ind = row.block_colidx(iEntry);
                      Kokkos::atomic_add(
                          &m_y(block_dim * ind + ir, ic),
                          static_cast<y_value_type>(
                              alpha * val * m_x(block_dim * iBlock + jr, ic)));
                    }
                  });
            }
          }
          //
        });
  }
};

/* ******************* */

/// \brief  spMatMultiVec_transpose: version for CPU execution spaces
/// (RangePolicy or trivial serial impl used)
template <class AT, class AO, class AD, class AS, class AlphaType,
          class XVector, class BetaType, class YVector,
          typename std::enable_if<!KokkosKernels::Impl::kk_is_gpu_exec_space<
              typename YVector::execution_space>()>::type * = nullptr>
void spMatMultiVec_transpose(
    const KokkosKernels::Experimental::Controls &controls,
    const AlphaType &alpha,
    const KokkosSparse::Experimental::BlockCrsMatrix<
        AT, AO, AD, Kokkos::MemoryTraits<Kokkos::Unmanaged>, AS> &A,
    const XVector &x, const BetaType &beta, YVector &y, bool useConjugate) {
  // This is required to maintain semantics of KokkosKernels native SpMV:
  // if y contains NaN but beta = 0, the result y should be filled with 0.
  // For example, this is useful for passing in uninitialized y and beta=0.
  if (beta == Kokkos::ArithTraits<BetaType>::zero())
    Kokkos::deep_copy(y, Kokkos::ArithTraits<BetaType>::zero());
  else
    KokkosBlas::scal(y, beta, y);
  //
  // Treat the case y <- alpha * A^T * x + beta * y
  //
  typedef KokkosSparse::Experimental::BlockCrsMatrix<
      AT, AO, AD, Kokkos::MemoryTraits<Kokkos::Unmanaged>, AS>
      AMatrix_Internal;
  typedef typename AMatrix_Internal::execution_space execution_space;

  AMatrix_Internal A_internal = A;

  bool use_dynamic_schedule = false;  // Forces the use of a dynamic schedule
  bool use_static_schedule  = false;  // Forces the use of a static schedule
  if (controls.isParameter("schedule")) {
    if (controls.getParameter("schedule") == "dynamic") {
      use_dynamic_schedule = true;
    } else if (controls.getParameter("schedule") == "static") {
      use_static_schedule = true;
    }
  }

  BCRS_GEMM_Transpose_Functor<AMatrix_Internal, XVector, YVector> func(
      alpha, A_internal, x, y, 1, useConjugate);
  if (((A.nnz() > 10000000) || use_dynamic_schedule) && !use_static_schedule) {
    Kokkos::parallel_for(
        "KokkosSparse::blockcrs_spm_mv<Transpose,Dynamic>",
        Kokkos::RangePolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>>(
            0, A.numRows()),
        func);
  } else {
    Kokkos::parallel_for(
        "KokkosSparse::blockcrs_spm_mv<Transpose,Static>",
        Kokkos::RangePolicy<execution_space, Kokkos::Schedule<Kokkos::Static>>(
            0, A.numRows()),
        func);
  }
}

//
// spMatMultiVec_transpose: version for GPU execution spaces (TeamPolicy used)
//
template <class AMatrix, class AlphaType, class XVector, class BetaType,
          class YVector,
          typename std::enable_if<KokkosKernels::Impl::kk_is_gpu_exec_space<
              typename YVector::execution_space>()>::type * = nullptr>
void spMatMultiVec_transpose(
    const KokkosKernels::Experimental::Controls &controls,
    const AlphaType &alpha, const AMatrix &A, const XVector &x,
    const BetaType &beta, YVector &y, bool useConjugate) {
  if (A.numRows() <= 0) {
    return;
  }

  KokkosBlas::scal(y, beta, y);

  typedef typename AMatrix::execution_space execution_space;

  bool use_dynamic_schedule = false;  // Forces the use of a dynamic schedule
  bool use_static_schedule  = false;  // Forces the use of a static schedule
  if (controls.isParameter("schedule")) {
    if (controls.getParameter("schedule") == "dynamic") {
      use_dynamic_schedule = true;
    } else if (controls.getParameter("schedule") == "static") {
      use_static_schedule = true;
    }
  }
  int team_size             = -1;
  int vector_length         = -1;
  int64_t blocks_per_thread = -1;

  //
  // Use the controls to allow the user to pass in some tuning
  // parameters.
  //
  if (controls.isParameter("team size")) {
    team_size = std::stoi(controls.getParameter("team size"));
  }
  if (controls.isParameter("vector length")) {
    vector_length = std::stoi(controls.getParameter("vector length"));
  }
  if (controls.isParameter("rows per thread")) {
    blocks_per_thread = std::stoll(controls.getParameter("rows per thread"));
  }

  //
  // Use the existing launch parameters routine from SPMV
  //
  int64_t blocks_per_team =
      KokkosSparse::Impl::spmv_launch_parameters<execution_space>(
          A.numRows(), A.nnz(), blocks_per_thread, team_size, vector_length);
  int64_t worksets = (A.numRows() + blocks_per_team - 1) / blocks_per_team;

  BCRS_GEMM_Transpose_Functor<AMatrix, XVector, YVector> func(
      alpha, A, x, y, blocks_per_team, useConjugate);

  if (((A.nnz() > 10000000) || use_dynamic_schedule) && !use_static_schedule) {
    Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>>
        policy(1, 1);
    if (team_size < 0)
      policy = Kokkos::TeamPolicy<execution_space,
                                  Kokkos::Schedule<Kokkos::Dynamic>>(
          worksets, Kokkos::AUTO, vector_length);
    else
      policy = Kokkos::TeamPolicy<execution_space,
                                  Kokkos::Schedule<Kokkos::Dynamic>>(
          worksets, team_size, vector_length);
    Kokkos::parallel_for("KokkosSparse::blockcrs_spm_mv<Transpose,Dynamic>",
                         policy, func);
  } else {
    Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>>
        policy(1, 1);
    if (team_size < 0)
      policy =
          Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>>(
              worksets, Kokkos::AUTO, vector_length);
    else
      policy =
          Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>>(
              worksets, team_size, vector_length);
    Kokkos::parallel_for("KokkosSparse::blockcrs_spm_mv<Transpose, Static>",
                         policy, func);
  }
}

/* ******************* */

}  // namespace BCRS

}  // namespace Impl
}  // namespace Experimental
}  // namespace KokkosSparse

#endif  // KOKKOSKERNELS_KOKKOSSPARSE_SPMV_BLOCKCRSMATRIX_IMPL_HPP
