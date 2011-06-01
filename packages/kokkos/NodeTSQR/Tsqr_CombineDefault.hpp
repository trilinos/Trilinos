//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef __TSQR_CombineDefault_hpp
#define __TSQR_CombineDefault_hpp

#include <Teuchos_ScalarTraits.hpp>

#include <Tsqr_ApplyType.hpp>
#include <Tsqr_Lapack.hpp>
#include <Tsqr_Matrix.hpp>

#include <algorithm>
#include <sstream>
#include <stdexcept>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  /// \class CombineDefault
  /// \brief Default implementation of TSQR::Combine
  ///
  /// This is a default implementation of TSQR::Combine, which
  /// TSQR::Combine may use (via a "has-a" relationship) if it doesn't
  /// have a specialized, faster implementation.  This default
  /// implementation copies the inputs into a contiguous matrix
  /// buffer, operates on them there via standard LAPACK calls, and
  /// copies out the results again.  It truncates to zero any values
  /// that should be zero because of the input's structure (e.g.,
  /// upper triangular).
  template< class Ordinal, class Scalar >
  class CombineDefault {
  private:
    typedef LAPACK< Ordinal, Scalar > lapack_type;

  public:
    typedef Scalar scalar_type;
    typedef typename Teuchos::ScalarTraits< Scalar >::magnitudeType magnitude_type;
    typedef Ordinal ordinal_type;

    CombineDefault () {}

    /// Whether or not the QR factorizations computed by methods of
    /// this class produce an R factor with all nonnegative diagonal
    /// entries.  
    static bool QR_produces_R_factor_with_nonnegative_diagonal()
    { 
      return lapack_type::QR_produces_R_factor_with_nonnegative_diagonal();
    }

    void
    factor_first (const Ordinal nrows,
		  const Ordinal ncols,
		  Scalar A[],
		  const Ordinal lda,
		  Scalar tau[],
		  Scalar work[])
    {
      // info must be an int, not a LocalOrdinal, since LAPACK
      // routines always (???) use int for the INFO output argument.
      int info = 0;
      lapack_.GEQR2 (nrows, ncols, A, lda, tau, work, &info);
      if (info != 0)
	{
	  std::ostringstream os;
	  os << "TSQR::CombineDefault::factor_first(): LAPACK\'s "
	     << "GEQR2 failed with INFO = " << info;
	  throw std::logic_error (os.str());
	}
    }

    void
    apply_first (const ApplyType& applyType,
		 const Ordinal nrows,
		 const Ordinal ncols_C,
		 const Ordinal ncols_A,
		 const Scalar A[],
		 const Ordinal lda,
		 const Scalar tau[],
		 Scalar C[],
		 const Ordinal ldc,
		 Scalar work[])
    {
      int info = 0;
      // LAPACK has the nice feature that it only reads the first
      // letter of input strings that specify things like which side
      // to which to apply the operator, or whether to apply the
      // transpose.  That means we can make the strings more verbose,
      // as in "Left" here for the SIDE parameter.
      lapack_.ORM2R ("Left", applyType.toString().c_str(), 
		     nrows, ncols_C, ncols_A,
		     A, lda, tau, 
		     C, ldc, work, &info);
      if (info != 0)
	{
	  std::ostringstream os;
	  os << "TSQR::CombineDefault::apply_first(): LAPACK\'s "
	     << "ORM2R failed with INFO = " << info;
	  throw std::logic_error (os.str());
	}
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
      typedef ConstMatView< Ordinal, Scalar > const_matview_type;
      typedef MatView< Ordinal, Scalar > matview_type;
      const Ordinal numRows = m + ncols_Q;

      A_buf_.reshape (numRows, ncols_Q);
      A_buf_.fill (Scalar(0));
      const_matview_type A_bot (m, ncols_Q, A, lda);
      matview_type A_buf_bot (m, ncols_Q, &A_buf_(ncols_Q, 0), A_buf_.lda());
      A_buf_bot.copy (A_bot);

      C_buf_.reshape (numRows, ncols_C);
      C_buf_.fill (Scalar(0));
      matview_type C_buf_top (ncols_Q, ncols_C, &C_buf_(0, 0), C_buf_.lda());
      matview_type C_buf_bot (m, ncols_C, &C_buf_(ncols_Q, 0), C_buf_.lda());
      matview_type C_top_view (ncols_Q, ncols_C, C_top, ldc_top);
      matview_type C_bot_view (m, ncols_C, C_bot, ldc_bot);
      C_buf_top.copy (C_top_view);
      C_buf_bot.copy (C_bot_view);

      int info = 0;
      lapack_.ORM2R ("Left", apply_type.toString().c_str(), 
		     numRows, ncols_C, ncols_Q,
		     A_buf_.get(), A_buf_.lda(), tau,
		     C_buf_.get(), C_buf_.lda(),
		     work, &info);
      if (info != 0)
	{
	  std::ostringstream os;
	  os << "TSQR::CombineDefault::apply_inner(): LAPACK\'s "
	     << "ORM2R failed with INFO = " << info;
	  throw std::logic_error (os.str());
	}
      // Copy back the results.
      C_top_view.copy (C_buf_top);
      C_bot_view.copy (C_buf_bot);
    }

    void
    factor_inner (const Ordinal m,
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
      A_buf_.fill (Scalar(0));
      // R might be a view of the upper triangle of a cache block, but
      // we only want to include the upper triangle in the
      // factorization.  Thus, only copy the upper triangle of R into
      // the appropriate place in the buffer.
      copy_upper_triangle (n, n, &A_buf_(0, 0), A_buf_.lda(), R, ldr);
      copy_matrix (m, n, &A_buf_(n, 0), A_buf_.lda(), A, lda);

      int info = 0;
      lapack_.GEQR2 (numRows, n, A_buf_.get(), A_buf_.lda(), tau, work, &info);
      if (info != 0)
	throw std::logic_error ("TSQR::CombineDefault: GEQR2 failed");

      // Copy back the results.  R might be a view of the upper
      // triangle of a cache block, so only copy into the upper
      // triangle of R.
      copy_upper_triangle (n, n, R, ldr, &A_buf_(0, 0), A_buf_.lda());
      copy_matrix (m, n, A, lda, &A_buf_(n, 0), A_buf_.lda());
    }

    void
    factor_pair (const Ordinal n,
		 Scalar R_top[],
		 const Ordinal ldr_top,
		 Scalar R_bot[],
		 const Ordinal ldr_bot,
		 Scalar tau[],
		 Scalar work[])
    {
      const Ordinal numRows = Ordinal(2) * n;

      A_buf_.reshape (numRows, n);
      A_buf_.fill (Scalar(0));
      // Copy the inputs into the compute buffer.  Only touch the
      // upper triangles of R_top and R_bot, since they each may be
      // views of some cache block (where the strict lower triangle
      // contains things we don't want to include in the
      // factorization).
      copy_upper_triangle (n, n, &A_buf_(0, 0), A_buf_.lda(), R_top, ldr_top);
      copy_upper_triangle (n, n, &A_buf_(n, 0), A_buf_.lda(), R_bot, ldr_bot);

      int info = 0;
      lapack_.GEQR2 (numRows, n, A_buf_.get(), A_buf_.lda(), tau, work, &info);
      if (info != 0)
	{
	  std::ostringstream os;
	  os << "TSQR::CombineDefault::factor_pair(): "
	     << "GEQR2 failed with INFO = " << info;
	  throw std::logic_error (os.str());
	}

      // Copy back the results.  Only read the upper triangles of the
      // two n by n row blocks of A_buf_ (this means we don't have to
      // zero out the strict lower triangles), and only touch the
      // upper triangles of R_top and R_bot.
      copy_upper_triangle (n, n, R_top, ldr_top, &A_buf_(0, 0), A_buf_.lda());
      copy_upper_triangle (n, n, R_bot, ldr_bot, &A_buf_(n, 0), A_buf_.lda());
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
      A_buf_.fill (Scalar(0));
      copy_upper_triangle (ncols_Q, ncols_Q, 
			   &A_buf_(ncols_Q, 0), A_buf_.lda(), 
			   R_bot, ldr_bot);
      C_buf_.reshape (numRows, ncols_C);
      copy_matrix (ncols_Q, ncols_C, &C_buf_(0, 0), C_buf_.lda(), C_top, ldc_top);
      copy_matrix (ncols_Q, ncols_C, &C_buf_(ncols_Q, 0), C_buf_.lda(), C_bot, ldc_bot);

      int info = 0;
      lapack_.ORM2R ("Left", apply_type.toString().c_str(), 
		     numRows, ncols_C, ncols_Q,
		     A_buf_.get(), A_buf_.lda(), tau,
		     C_buf_.get(), C_buf_.lda(),
		     work, &info);
      if (info != 0)
	throw std::logic_error ("TSQR::CombineDefault: ORM2R failed");

      // Copy back the results.
      copy_matrix (ncols_Q, ncols_C, C_top, ldc_top, &C_buf_(0, 0), C_buf_.lda());
      copy_matrix (ncols_Q, ncols_C, C_bot, ldc_bot, &C_buf_(ncols_Q, 0), C_buf_.lda());
    }

  private:
    LAPACK< Ordinal, Scalar > lapack_;
    Matrix< Ordinal, Scalar > A_buf_;
    Matrix< Ordinal, Scalar > C_buf_;
  };


} // namespace TSQR

#endif // __TSQR_CombineDefault_hpp
