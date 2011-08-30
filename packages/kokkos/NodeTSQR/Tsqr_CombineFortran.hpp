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

/// \file Tsqr_CombineFortran.hpp
/// \brief Interface to Fortran 9x back end of \c TSQR::Combine.
///
#ifndef __TSQR_CombineFortran_hpp
#define __TSQR_CombineFortran_hpp

#include <Tsqr_ApplyType.hpp>
#include <Tsqr_MatView.hpp>
#include <Tsqr_ScalarTraits.hpp>
#include <Tsqr_CombineDefault.hpp>


namespace TSQR {

  /// \class CombineFortran
  /// \brief Interface to Fortran 9x back end of \c TSQR::Combine.
  ///
  /// TSQR::Combine has three implementations: CombineDefault,
  /// CombineNative, and CombineFortran.  The latter, implemented in
  /// this file, is a C++ front end to a Fortran 9x implementation.
  /// CombineFortran is not templated on the Ordinal type, because the
  /// Fortran implementation uses int for that.
  ///
  template<class Scalar, bool is_complex = ScalarTraits<Scalar>::is_complex >
  class CombineFortran {
  private:
    typedef CombineDefault<int, Scalar> combine_default_type;

  public:
    typedef Scalar scalar_type;
    typedef typename ScalarTraits<Scalar>::magnitude_type magnitude_type;
    typedef int ordinal_type;

    CombineFortran () {}

    /// \brief Does the R factor have a nonnegative diagonal?
    ///
    /// CombineFortran implements a QR factorization (of a matrix with
    /// a special structure).  Some, but not all, QR factorizations
    /// produce an R factor whose diagonal may include negative
    /// entries.  This Boolean tells you whether CombineFortran
    /// promises to compute an R factor whose diagonal entries are all
    /// nonnegative.
    static bool QR_produces_R_factor_with_nonnegative_diagonal();

    void
    factor_first (const ordinal_type nrows,
		  const ordinal_type ncols,
		  Scalar A[],
		  const ordinal_type lda,
		  Scalar tau[],
		  Scalar work[]) const;
    
    void
    apply_first (const ApplyType& applyType,
		 const ordinal_type nrows,
		 const ordinal_type ncols_C,
		 const ordinal_type ncols_A,
		 const Scalar A[],
		 const ordinal_type lda,
		 const Scalar tau[],
		 Scalar C[],
		 const ordinal_type ldc,
		 Scalar work[]) const;

    void
    apply_inner (const ApplyType& apply_type,
		 const ordinal_type m,
		 const ordinal_type ncols_C,
		 const ordinal_type ncols_Q,
		 const Scalar A[],
		 const ordinal_type lda,
		 const Scalar tau[],
		 Scalar C_top[],
		 const ordinal_type ldc_top,
		 Scalar C_bot[],
		 const ordinal_type ldc_bot,
		 Scalar work[]) const;

    void
    factor_inner (const ordinal_type m,
		  const ordinal_type n,
		  Scalar R[],
		  const ordinal_type ldr,
		  Scalar A[],
		  const ordinal_type lda,
		  Scalar tau[],
		  Scalar work[]) const;

    void
    factor_pair (const ordinal_type n,
		 Scalar R_top[],
		 const ordinal_type ldr_top,
		 Scalar R_bot[],
		 const ordinal_type ldr_bot,
		 Scalar tau[],
		 Scalar work[]) const;
    
    void
    apply_pair (const ApplyType& apply_type,
		const ordinal_type ncols_C, 
		const ordinal_type ncols_Q, 
		const Scalar R_bot[], 
		const ordinal_type ldr_bot,
		const Scalar tau[], 
		Scalar C_top[], 
		const ordinal_type ldc_top, 
		Scalar C_bot[], 
		const ordinal_type ldc_bot, 
		Scalar work[]) const;

  private:
    mutable combine_default_type default_;
  };

  // "Forward declaration" for the real-arithmetic case.  The Fortran
  // back end works well here for Scalar = {float, double}.
  template< class Scalar >
  class CombineFortran< Scalar, false > {
  private:
    typedef CombineDefault< int, Scalar > combine_default_type;

  public:
    typedef Scalar scalar_type;
    typedef typename ScalarTraits< Scalar >::magnitude_type magnitude_type;
    typedef int ordinal_type;

    CombineFortran () {}

    static bool QR_produces_R_factor_with_nonnegative_diagonal() {
      typedef LAPACK< int, Scalar > lapack_type;

      return lapack_type::QR_produces_R_factor_with_nonnegative_diagonal() &&
	combine_default_type::QR_produces_R_factor_with_nonnegative_diagonal();
    }

    void
    factor_first (const ordinal_type nrows,
		  const ordinal_type ncols,
		  Scalar A[],
		  const ordinal_type lda,
		  Scalar tau[],
		  Scalar work[]) const
    {
      return default_.factor_first (nrows, ncols, A, lda, tau, work);
    }
    
    void
    apply_first (const ApplyType& applyType,
		 const ordinal_type nrows,
		 const ordinal_type ncols_C,
		 const ordinal_type ncols_A,
		 const Scalar A[],
		 const ordinal_type lda,
		 const Scalar tau[],
		 Scalar C[],
		 const ordinal_type ldc,
		 Scalar work[]) const
    {
      return default_.apply_first (applyType, nrows, ncols_C, ncols_A, 
				   A, lda, tau, 
				   C, ldc, work);
    }

    void
    apply_inner (const ApplyType& apply_type,
		 const ordinal_type m,
		 const ordinal_type ncols_C,
		 const ordinal_type ncols_Q,
		 const Scalar A[],
		 const ordinal_type lda,
		 const Scalar tau[],
		 Scalar C_top[],
		 const ordinal_type ldc_top,
		 Scalar C_bot[],
		 const ordinal_type ldc_bot,
		 Scalar work[]) const;

    void
    factor_inner (const ordinal_type m,
		  const ordinal_type n,
		  Scalar R[],
		  const ordinal_type ldr,
		  Scalar A[],
		  const ordinal_type lda,
		  Scalar tau[],
		  Scalar work[]) const;

    void
    factor_pair (const ordinal_type n,
		 Scalar R_top[],
		 const ordinal_type ldr_top,
		 Scalar R_bot[],
		 const ordinal_type ldr_bot,
		 Scalar tau[],
		 Scalar work[]) const;
    
    void
    apply_pair (const ApplyType& apply_type,
		const ordinal_type ncols_C, 
		const ordinal_type ncols_Q, 
		const Scalar R_bot[], 
		const ordinal_type ldr_bot,
		const Scalar tau[], 
		Scalar C_top[], 
		const ordinal_type ldc_top, 
		Scalar C_bot[], 
		const ordinal_type ldc_bot, 
		Scalar work[]) const;

  private:
    mutable combine_default_type default_;
  };


  // "Forward declaration" for complex-arithmetic version of
  // CombineFortran.  The Fortran code doesn't actually work for this
  // case, so we implement everything using CombineDefault.  This
  // will likely result in an ~2x slowdown for typical use cases.
  template< class Scalar >
  class CombineFortran< Scalar, true > {
  private:
    typedef CombineDefault< int, Scalar > combine_default_type;

  public:
    typedef Scalar scalar_type;
    typedef typename ScalarTraits< Scalar >::magnitude_type magnitude_type;
    typedef int ordinal_type;

    CombineFortran () {}

    static bool QR_produces_R_factor_with_nonnegative_diagonal() {
      return combine_default_type::QR_produces_R_factor_with_nonnegative_diagonal();
    }

    void
    factor_first (const ordinal_type nrows,
		  const ordinal_type ncols,
		  Scalar A[],
		  const ordinal_type lda,
		  Scalar tau[],
		  Scalar work[]) const
    {
      return default_.factor_first (nrows, ncols, A, lda, tau, work);
    }
    
    void
    apply_first (const ApplyType& applyType,
		 const ordinal_type nrows,
		 const ordinal_type ncols_C,
		 const ordinal_type ncols_A,
		 const Scalar A[],
		 const ordinal_type lda,
		 const Scalar tau[],
		 Scalar C[],
		 const ordinal_type ldc,
		 Scalar work[]) const
    {
      return default_.apply_first (applyType, nrows, ncols_C, ncols_A, 
				   A, lda, tau, 
				   C, ldc, work);
    }

    void
    apply_inner (const ApplyType& apply_type,
		 const ordinal_type m,
		 const ordinal_type ncols_C,
		 const ordinal_type ncols_Q,
		 const Scalar A[],
		 const ordinal_type lda,
		 const Scalar tau[],
		 Scalar C_top[],
		 const ordinal_type ldc_top,
		 Scalar C_bot[],
		 const ordinal_type ldc_bot,
		 Scalar work[]) const
    {
      default_.apply_inner (apply_type, m, ncols_C, ncols_Q, 
			    A, lda, tau, 
			    C_top, ldc_top, C_bot, ldc_bot, work);
    }

    void
    factor_inner (const ordinal_type m,
		  const ordinal_type n,
		  Scalar R[],
		  const ordinal_type ldr,
		  Scalar A[],
		  const ordinal_type lda,
		  Scalar tau[],
		  Scalar work[]) const
    {
      default_.factor_inner (m, n, R, ldr, A, lda, tau, work);
    }

    void
    factor_pair (const ordinal_type n,
		 Scalar R_top[],
		 const ordinal_type ldr_top,
		 Scalar R_bot[],
		 const ordinal_type ldr_bot,
		 Scalar tau[],
		 Scalar work[]) const
    {
      default_.factor_pair (n, R_top, ldr_top, R_bot, ldr_bot, tau, work);
    }
    
    void
    apply_pair (const ApplyType& apply_type,
		const ordinal_type ncols_C, 
		const ordinal_type ncols_Q, 
		const Scalar R_bot[], 
		const ordinal_type ldr_bot,
		const Scalar tau[], 
		Scalar C_top[], 
		const ordinal_type ldc_top, 
		Scalar C_bot[], 
		const ordinal_type ldc_bot, 
		Scalar work[]) const
    {
      default_.apply_pair (apply_type, ncols_C, ncols_Q, 
			   R_bot, ldr_bot, tau, 
			   C_top, ldc_top, C_bot, ldc_bot, work);
    }

  private:
    // Default implementation of TSQR::Combine copies data in and out
    // of a single matrix, which is given to LAPACK.  It's slow
    // because we expect the number of columns to be small, so copying
    // overhead is significant.  Experiments have shown a ~2x slowdown
    // due to copying overhead.
    mutable CombineDefault< ordinal_type, scalar_type > default_;
  };


} // namespace TSQR

#endif // __TSQR_CombineFortran_hpp
