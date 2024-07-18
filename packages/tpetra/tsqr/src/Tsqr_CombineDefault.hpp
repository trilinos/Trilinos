// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Tsqr_CombineDefault.hpp
/// \brief Default copy-in, copy-out implementation of TSQR::Combine.

#ifndef TSQR_COMBINEDEFAULT_HPP
#define TSQR_COMBINEDEFAULT_HPP

#include "Tsqr_Combine.hpp"
#include "Tsqr_Impl_Lapack.hpp"
#include "Tsqr_Matrix.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace TSQR {

  /// \class CombineDefault
  /// \brief Default copy-in, copy-out implementation of \c TSQR::Combine.
  ///
  /// This is a default implementation of TSQR::Combine, which
  /// TSQR::Combine may use (via a "has-a" relationship) if it doesn't
  /// have a specialized, faster implementation.  This default
  /// implementation copies the inputs into a contiguous matrix
  /// buffer, operates on them there via standard LAPACK calls, and
  /// copies out the results again.  It truncates to zero any values
  /// that should be zero because of the input's structure (e.g.,
  /// upper triangular).
  template<class Ordinal, class Scalar>
  class CombineDefault : public Combine<Ordinal, Scalar> {
  public:
    using ordinal_type = Ordinal;
    using scalar_type = Scalar;
    using const_mat_view_type = MatView<ordinal_type, const Scalar>;
    using mat_view_type = MatView<ordinal_type, Scalar>;

    ~CombineDefault () override = default;

    /// \brief Does the R factor have a nonnegative diagonal?
    ///
    /// CombineDefault implements a QR factorization (of a matrix with
    /// a special structure).  Some, but not all, QR factorizations
    /// produce an R factor whose diagonal may include negative
    /// entries.  This Boolean tells you whether CombineDefault
    /// promises to compute an R factor whose diagonal entries are all
    /// nonnegative.
    bool
    QR_produces_R_factor_with_nonnegative_diagonal () const override
    {
      // FIXME (mfh 19 Dec 2019) This _should_ depend on Impl::Lapack.
      return false;
    }

    ordinal_type
    work_size (const ordinal_type num_rows_Q,
               const ordinal_type num_cols_Q,
               const ordinal_type num_cols_C) const override
    {
      const int ncols = num_cols_Q < num_cols_C ?
        num_cols_C : num_cols_Q;
      const int nrows = num_rows_Q + ncols;
      const int lda = nrows;

      const int lwork1 =
        lapack_.compute_QR_lwork (nrows, ncols, nullptr, lda);
      TEUCHOS_ASSERT( lwork1 >= num_cols_Q );

      const int ldc = nrows;
      const int lwork2 =
        lapack_.apply_Q_factor_lwork ('L', 'N',
                                      nrows, num_cols_C, num_cols_Q,
                                      nullptr, lda, nullptr,
                                      nullptr, ldc);
      TEUCHOS_ASSERT( lwork2 >= 0 );
      return std::max (lwork1, lwork2);
    }

    void
    factor_first (const MatView<ordinal_type, Scalar>& A,
                  Scalar tau[],
                  Scalar work[],
                  const ordinal_type lwork) override
    {
      lapack_.compute_QR (A.extent (0), A.extent (1),
                          A.data (), A.stride (1),
                          tau, work, lwork);
    }

    void
    factor_first (Matrix<ordinal_type, Scalar>& A,
                  Scalar tau[],
                  Scalar work[],
                  const ordinal_type lwork)
    {
      MatView<ordinal_type, Scalar> A_view
        (A.extent (0), A.extent (1), A.data (), A.stride (1));
      this->factor_first (A_view, tau, work, lwork);
    }

    void
    apply_first (const ApplyType& applyType,
                 const MatView<ordinal_type, const Scalar>& A,
                 const Scalar tau[],
                 const MatView<ordinal_type, Scalar>& C,
                 Scalar work[],
                 const ordinal_type lwork) override
    {
      const ordinal_type nrows = A.extent(0);
      const ordinal_type ncols_C = C.extent(1);
      const ordinal_type ncols_A = A.extent(1);
      const ordinal_type lda = A.stride(1);
      const ordinal_type ldc = C.stride(1);

      // LAPACK has the nice feature that it only reads the first
      // letter of input strings that specify things like which side
      // to which to apply the operator, or whether to apply the
      // transpose.  That means we can make the strings more verbose,
      // as in "Left" here for the SIDE parameter.
      const std::string trans = applyType.toString ();
      lapack_.apply_Q_factor ('L', trans[0], nrows, ncols_C, ncols_A,
                              A.data(), lda, tau, C.data(), ldc,
                              work, static_cast<int> (lwork));
    }

    void
    factor_inner (const MatView<ordinal_type, Scalar>& R,
                  const MatView<ordinal_type, Scalar>& A,
                  Scalar tau[],
                  Scalar work[],
                  const ordinal_type lwork) override
    {
      const ordinal_type m = A.extent (0);
      const ordinal_type n = A.extent (1);
      const ordinal_type lda = A.stride (1);
      factor_inner_impl (m, n, R.data (), R.stride (1),
                         A.data (), lda, tau, work, lwork);
    }

    void
    apply_inner (const ApplyType& apply_type,
                 const MatView<ordinal_type, const Scalar>& A,
                 const Scalar tau[],
                 const MatView<ordinal_type, Scalar>& C_top,
                 const MatView<ordinal_type, Scalar>& C_bot,
                 Scalar work[],
                 const ordinal_type lwork) override
    {
      const ordinal_type m = A.extent (0);
      TEUCHOS_ASSERT( m == ordinal_type (C_bot.extent (0)) );
      const ordinal_type ncols_Q = A.extent (1);
      const ordinal_type ncols_C = C_top.extent (1);
      TEUCHOS_ASSERT( ncols_C == ordinal_type (C_bot.extent (1)) );
      const ordinal_type numRows = ncols_Q + m;

      A_buf_.reshape (numRows, ncols_Q);
      deep_copy (A_buf_, Scalar {});
      auto A_buf_top_bot = partition_2x1 (A_buf_.view (), ncols_Q);
      deep_copy (A_buf_top_bot.second, A);

      C_buf_.reshape (numRows, ncols_C);
      deep_copy (C_buf_, Scalar {});
      auto C_buf_top_bot = partition_2x1 (C_buf_.view (), ncols_Q);
      deep_copy (C_buf_top_bot.first, C_top);
      deep_copy (C_buf_top_bot.second, C_bot);

      const std::string trans = apply_type.toString ();
      lapack_.apply_Q_factor ('L', trans[0],
                              numRows, ncols_C, ncols_Q,
                              A_buf_.data (), A_buf_.stride (1), tau,
                              C_buf_.data (), C_buf_.stride (1),
                              work, lwork);
      // Copy back the results.
      deep_copy (C_top, C_buf_top_bot.first);
      deep_copy (C_bot, C_buf_top_bot.second);
    }

  private:
    void
    factor_inner_impl (const ordinal_type m,
                       const ordinal_type n,
                       Scalar R[],
                       const ordinal_type ldr,
                       Scalar A[],
                       const ordinal_type lda,
                       Scalar tau[],
                       Scalar work[],
                       const ordinal_type lwork)
    {
      const ordinal_type numRows = m + n;

      A_buf_.reshape (numRows, n);
      deep_copy (A_buf_, Scalar {});
      // R might be a view of the upper triangle of a cache block, but
      // we only want to include the upper triangle in the
      // factorization.  Thus, only copy the upper triangle of R into
      // the appropriate place in the buffer.
      MatView<ordinal_type, Scalar> R_view (n, n, R, ldr);
      MatView<ordinal_type, Scalar> A_buf_top (n, n, A_buf_.data(),
                                          A_buf_.stride(1));
      deep_copy (A_buf_top, R_view);

      MatView<ordinal_type, Scalar> A_view (m, n, A, lda);
      MatView<ordinal_type, Scalar> A_buf_bot (m, n, &A_buf_(n, 0),
                                          A_buf_.stride(1));
      deep_copy (A_buf_bot, A_view);
      lapack_.compute_QR (numRows, n, A_buf_.data (),
                          A_buf_.stride (1), tau, work, lwork);
      // Copy back the results.  R might be a view of the upper
      // triangle of a cache block, so only copy into the upper
      // triangle of R.
      copy_upper_triangle (R_view, A_buf_top);
      deep_copy (A_view, A_buf_bot);
    }

  public:
    void
    factor_pair (const MatView<ordinal_type, Scalar>& R_top,
                 const MatView<ordinal_type, Scalar>& R_bot,
                 Scalar tau[],
                 Scalar work[],
                 const ordinal_type lwork) override
    {
      const ordinal_type numRows = ordinal_type(2) * R_top.extent (1);
      const ordinal_type numCols = R_top.extent (1);

      A_buf_.reshape (numRows, numCols);
      deep_copy (A_buf_, Scalar {});
      auto A_buf_tb = partition_2x1 (A_buf_.view (), numCols);
      // Copy the inputs into the compute buffer.  Only touch the
      // upper triangles of R_top and R_bot, since they each may be
      // views of some cache block (where the strict lower triangle
      // contains things we don't want to include in the
      // factorization).
      copy_upper_triangle (A_buf_tb.first, R_top);
      copy_upper_triangle (A_buf_tb.second, R_bot);

      lapack_.compute_QR (numRows, numCols,
                          A_buf_.data(), A_buf_.stride(1),
                          tau, work, lwork);
      // Copy back the results.  Only read the upper triangles of the
      // two n by n row blocks of A_buf_ (this means we don't have to
      // zero out the strict lower triangles), and only touch the
      // upper triangles of R_top and R_bot.
      copy_upper_triangle (R_top, A_buf_tb.first);
      copy_upper_triangle (R_bot, A_buf_tb.second);
    }

    void
    apply_pair (const ApplyType& apply_type,
                const MatView<ordinal_type, const Scalar>& R_bot,
                const Scalar tau[],
                const MatView<ordinal_type, Scalar>& C_top,
                const MatView<ordinal_type, Scalar>& C_bot,
                Scalar work[],
                const ordinal_type lwork) override
    {
      const ordinal_type ncols_C = C_top.extent (1);
      const ordinal_type ncols_Q = R_bot.extent (1);
      const ordinal_type numRows = ordinal_type(2) * ncols_Q;

      A_buf_.reshape (numRows, ncols_Q);
      deep_copy (A_buf_, Scalar {});
      auto A_buf_tb = partition_2x1 (A_buf_.view (), ncols_Q);
      copy_upper_triangle (A_buf_tb.second, R_bot);

      C_buf_.reshape (numRows, ncols_C);
      auto C_buf_tb = partition_2x1 (C_buf_.view (), ncols_Q);
      deep_copy (C_buf_tb.first, C_top);
      deep_copy (C_buf_tb.second, C_bot);

      const std::string trans = apply_type.toString ();
      lapack_.apply_Q_factor ('L', trans[0], numRows, ncols_C,
                              ncols_Q, A_buf_.data (),
                              A_buf_.stride (1), tau,
                              C_buf_.data (), C_buf_.stride (1),
                              work, lwork);
      // Copy back the results.
      deep_copy (C_top, C_buf_tb.first);
      deep_copy (C_bot, C_buf_tb.second);
    }

  private:
    Impl::Lapack<Scalar> lapack_;
    Matrix<ordinal_type, Scalar> A_buf_;
    Matrix<ordinal_type, Scalar> C_buf_;
  };
} // namespace TSQR

#endif // TSQR_COMBINEDEFAULT_HPP
