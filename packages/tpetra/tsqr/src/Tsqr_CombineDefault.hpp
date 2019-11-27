//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************
//@HEADER

/// \file Tsqr_CombineDefault.hpp
/// \brief Default copy-in, copy-out implementation of \c TSQR::Combine.
///
#ifndef __TSQR_CombineDefault_hpp
#define __TSQR_CombineDefault_hpp

#include "Teuchos_ScalarTraits.hpp"
#include "Tsqr_ApplyType.hpp"
#include "Tsqr_Impl_Lapack.hpp"
#include "Tsqr_Matrix.hpp"

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
  class CombineDefault {
  public:
    typedef Ordinal ordinal_type;
    typedef Scalar scalar_type;
    typedef typename Teuchos::ScalarTraits< Scalar >::magnitudeType magnitude_type;
    typedef MatView<Ordinal, const Scalar> const_mat_view_type;
    typedef MatView<Ordinal, Scalar> mat_view_type;

    /// \brief Does the R factor have a nonnegative diagonal?
    ///
    /// CombineDefault implements a QR factorization (of a matrix with
    /// a special structure).  Some, but not all, QR factorizations
    /// produce an R factor whose diagonal may include negative
    /// entries.  This Boolean tells you whether CombineDefault
    /// promises to compute an R factor whose diagonal entries are all
    /// nonnegative.
    static bool QR_produces_R_factor_with_nonnegative_diagonal()
    {
      return false; // lapack_type::QR_produces_R_factor_with_nonnegative_diagonal();
    }

    void
    factor_first (const MatView<Ordinal, Scalar>& A,
                  Scalar tau[],
                  Scalar work[])
    {
      const int lwork = A.extent (1);
      lapack_.compute_QR (A.extent (0), A.extent (1),
                          A.data (), A.stride (1),
                          tau, work, lwork);
    }

    void
    factor_first (Matrix<Ordinal, Scalar>& A,
                  Scalar tau[],
                  Scalar work[])
    {
      MatView<Ordinal, Scalar> A_view
        (A.extent (0), A.extent (1), A.data (), A.stride (1));
      factor_first (A_view, tau, work);
    }

    void
    apply_first (const ApplyType& applyType,
                 const MatView<Ordinal, const Scalar>& A,
                 const Scalar tau[],
                 const MatView<Ordinal, Scalar>& C,
                 Scalar work[])
    {
      const Ordinal nrows = A.extent(0);
      const Ordinal ncols_C = C.extent(1);
      const Ordinal ncols_A = A.extent(1);
      const Ordinal lda = A.stride(1);
      const Ordinal ldc = C.stride(1);

      // LAPACK has the nice feature that it only reads the first
      // letter of input strings that specify things like which side
      // to which to apply the operator, or whether to apply the
      // transpose.  That means we can make the strings more verbose,
      // as in "Left" here for the SIDE parameter.
      const std::string trans = applyType.toString ();
      const int lwork = ncols_C;
      lapack_.apply_Q_factor ('L', trans[0], nrows, ncols_C, ncols_A,
                              A.data(), lda, tau, C.data(), ldc,
                              work, lwork);
    }

    void
    apply_inner (const ApplyType& apply_type,
                 const Ordinal m,
                 const Ordinal ncols_C,
                 const Ordinal ncols_Q,
                 const Scalar A[],
                 const Ordinal lda,
                 const Scalar tau[],
                 Scalar C_top[],
                 const Ordinal ldc_top,
                 Scalar C_bot[],
                 const Ordinal ldc_bot,
                 Scalar work[])
    {
      const Ordinal numRows = m + ncols_Q;

      A_buf_.reshape (numRows, ncols_Q);
      deep_copy (A_buf_, Scalar {});
      const_mat_view_type A_bot (m, ncols_Q, A, lda);
      mat_view_type A_buf_bot (m, ncols_Q, &A_buf_(ncols_Q, 0), A_buf_.stride(1));
      deep_copy (A_buf_bot, A_bot);

      C_buf_.reshape (numRows, ncols_C);
      deep_copy (C_buf_, Scalar {});
      mat_view_type C_buf_top (ncols_Q, ncols_C, &C_buf_(0, 0), C_buf_.stride(1));
      mat_view_type C_buf_bot (m, ncols_C, &C_buf_(ncols_Q, 0), C_buf_.stride(1));
      mat_view_type C_top_view (ncols_Q, ncols_C, C_top, ldc_top);
      mat_view_type C_bot_view (m, ncols_C, C_bot, ldc_bot);
      deep_copy (C_buf_top, C_top_view);
      deep_copy (C_buf_bot, C_bot_view);

      const std::string trans = apply_type.toString ();
      const int lwork = ncols_C;
      lapack_.apply_Q_factor ('L', trans[0], numRows, ncols_C, ncols_Q,
                              A_buf_.data(), A_buf_.stride(1), tau,
                              C_buf_.data(), C_buf_.stride(1),
                              work, lwork);
      // Copy back the results.
      deep_copy (C_top_view, C_buf_top);
      deep_copy (C_bot_view, C_buf_bot);
    }

    void
    factor_inner (const MatView<Ordinal, Scalar>& R,
                  const MatView<Ordinal, Scalar>& A,
                  Scalar tau[],
                  Scalar work[])
    {
      const Ordinal m = A.extent(0);
      const Ordinal n = A.extent(1);
      factor_inner_impl (m, n, R.data(), R.stride(1),
                         A.data(), A.stride(1), tau, work);
    }

  private:
    void
    factor_inner_impl (const Ordinal m,
                       const Ordinal n,
                       Scalar R[],
                       const Ordinal ldr,
                       Scalar A[],
                       const Ordinal lda,
                       Scalar tau[],
                       Scalar work[])
    {
      const Ordinal numRows = m + n;

      A_buf_.reshape (numRows, n);
      deep_copy (A_buf_, Scalar {});
      // R might be a view of the upper triangle of a cache block, but
      // we only want to include the upper triangle in the
      // factorization.  Thus, only copy the upper triangle of R into
      // the appropriate place in the buffer.
      MatView<Ordinal, Scalar> R_view (n, n, R, ldr);
      MatView<Ordinal, Scalar> A_buf_top (n, n, A_buf_.data(),
                                          A_buf_.stride(1));
      deep_copy (A_buf_top, R_view);

      MatView<Ordinal, Scalar> A_view (m, n, A, lda);
      MatView<Ordinal, Scalar> A_buf_bot (m, n, &A_buf_(n, 0),
                                          A_buf_.stride(1));
      deep_copy (A_buf_bot, A_view);

      const int lwork = n;
      lapack_.compute_QR (numRows, n, A_buf_.data(), A_buf_.stride(1),
                          tau, work, lwork);
      // Copy back the results.  R might be a view of the upper
      // triangle of a cache block, so only copy into the upper
      // triangle of R.
      copy_upper_triangle (n, n, R, ldr, A_buf_top.data(),
                           A_buf_top.stride(1));
      deep_copy (A_view, A_buf_bot);
    }

  public:
    void
    factor_pair (const MatView<Ordinal, Scalar>& R_top,
                 const MatView<Ordinal, Scalar>& R_bot,
                 Scalar tau[],
                 Scalar work[])
    {
      const Ordinal numRows = Ordinal(2) * R_top.extent (1);
      const Ordinal numCols = R_top.extent (1);

      A_buf_.reshape (numRows, numCols);
      deep_copy (A_buf_, Scalar {});
      MatView<Ordinal, Scalar> A_buf_top (numCols, numCols,
                                          &A_buf_(0, 0),
                                          A_buf_.stride(1));
      MatView<Ordinal, Scalar> A_buf_bot (numCols, numCols,
                                          &A_buf_(numCols, 0),
                                          A_buf_.stride(1));
      // Copy the inputs into the compute buffer.  Only touch the
      // upper triangles of R_top and R_bot, since they each may be
      // views of some cache block (where the strict lower triangle
      // contains things we don't want to include in the
      // factorization).
      copy_upper_triangle (numCols, numCols,
                           A_buf_top.data(), A_buf_top.stride(1),
                           R_top.data(), R_top.stride(1));
      copy_upper_triangle (numCols, numCols,
                           A_buf_bot.data(), A_buf_bot.stride(1),
                           R_bot.data(), R_bot.stride(1));

      const int lwork = static_cast<int> (numCols);
      lapack_.compute_QR (numRows, numCols,
                          A_buf_.data(), A_buf_.stride(1),
                          tau, work, lwork);
      // Copy back the results.  Only read the upper triangles of the
      // two n by n row blocks of A_buf_ (this means we don't have to
      // zero out the strict lower triangles), and only touch the
      // upper triangles of R_top and R_bot.
      copy_upper_triangle (numCols, numCols,
                           R_top.data(), R_top.stride(1),
                           A_buf_top.data(), A_buf_top.stride(1));
      copy_upper_triangle (numCols, numCols,
                           R_bot.data(), R_bot.stride(1),
                           A_buf_bot.data(), A_buf_bot.stride(1));
    }

    void
    apply_pair (const ApplyType& apply_type,
                const Ordinal ncols_C,
                const Ordinal ncols_Q,
                const Scalar R_bot[],
                const Ordinal ldr_bot,
                const Scalar tau[],
                Scalar C_top[],
                const Ordinal ldc_top,
                Scalar C_bot[],
                const Ordinal ldc_bot,
                Scalar work[])
    {
      const Ordinal numRows = Ordinal(2) * ncols_Q;

      A_buf_.reshape (numRows, ncols_Q);
      deep_copy (A_buf_, Scalar {});
      copy_upper_triangle (ncols_Q, ncols_Q,
                           &A_buf_(ncols_Q, 0), A_buf_.stride(1),
                           R_bot, ldr_bot);
      C_buf_.reshape (numRows, ncols_C);

      using view_type = MatView<Ordinal, Scalar>;
      view_type C_top_view (ncols_Q, ncols_C, C_top, ldc_top);
      view_type C_buf_top (ncols_Q, ncols_C,
                           C_buf_.data (), C_buf_.stride (1));
      deep_copy (C_buf_top, C_top_view);

      view_type C_bot_view (ncols_Q, ncols_C, C_bot, ldc_bot);
      view_type C_buf_bot (ncols_Q, ncols_C,
                           &C_buf_(ncols_Q, 0), C_buf_.stride (1));
      deep_copy (C_buf_bot, C_bot_view);

      const int lwork = ncols_Q;
      const std::string trans = apply_type.toString ();
      lapack_.apply_Q_factor ('L', trans[0], numRows, ncols_C, ncols_Q,
                              A_buf_.data(), A_buf_.stride(1), tau,
                              C_buf_.data(), C_buf_.stride(1),
                              work, lwork);
      // Copy back the results.
      deep_copy (C_top_view, C_buf_top);
      deep_copy (C_bot_view, C_buf_bot);
    }

  private:
    Impl::Lapack<Scalar> lapack_;
    Matrix<Ordinal, Scalar> A_buf_;
    Matrix<Ordinal, Scalar> C_buf_;
  };
} // namespace TSQR

#endif // __TSQR_CombineDefault_hpp
