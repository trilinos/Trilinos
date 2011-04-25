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

#ifndef __TSQR_CombineNative_hpp
#define __TSQR_CombineNative_hpp

#include "Teuchos_ScalarTraits.hpp"

#include "Tsqr_ApplyType.hpp"
#include "Tsqr_Blas.hpp"
#include "Tsqr_Lapack.hpp"
#include "Tsqr_CombineDefault.hpp"

// #define TSQR_COMBINE_NATIVE_DEBUG 1
#ifdef TSQR_COMBINE_NATIVE_DEBUG
#include "Tsqr_Util.hpp"
#include <iostream>
using std::cerr;
using std::endl;
#endif // TSQR_COMBINE_NATIVE_DEBUG

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  /// \class CombineNative
  /// \brief Interface to C++ back end of TSQR::Combine
  ///
  /// TSQR::Combine has three implementations: CombineDefault,
  /// CombineNative, and CombineFortran.  CombineNative, implemented
  /// in this file, is a fully C++ (therefore "native," as opposed to
  /// CombineFortran (implemented in Fortran) or CombineNative
  /// (implemented by wrappers around LAPACK calls)) implementation.
  ///
  /// \warning CombineNative has no complex-arithmetic implementation
  /// yet.  It's not hard to implement this (use LAPACK's ZGEQR2(P)
  /// and ZUNM2R as models), but it will take time that the author
  /// doesn't have at the moment.
  ///
  template< class Ordinal, class Scalar, bool isComplex = Teuchos::ScalarTraits< Scalar >::isComplex >
  class CombineNative 
  {
  public:
    typedef Scalar scalar_type;
    typedef typename Teuchos::ScalarTraits< Scalar >::magnitudeType magnitude_type;
    typedef Ordinal ordinal_type;

  private:
    typedef BLAS< ordinal_type, scalar_type > blas_type;
    typedef LAPACK< ordinal_type, scalar_type > lapack_type;
    typedef CombineDefault< ordinal_type, scalar_type > combine_default_type;

  public:

    CombineNative () {}

    /// Whether or not the QR factorizations computed by methods of
    /// this class produce an R factor with all nonnegative diagonal
    /// entries.  It depends on LAPACK because this implementation
    /// invokes one of {LARFGP, LARFP, LARFG} in order to compute
    /// Householder reflectors; only LAPACK versions >= 3.2 have one
    /// of {LARFGP, LARFP}, which is necessary to ensure that the BETA
    /// output of the function is always nonnegative.
    static bool QR_produces_R_factor_with_nonnegative_diagonal() { 
      return lapack_type::QR_produces_R_factor_with_nonnegative_diagonal() &&
	combine_default_type::QR_produces_R_factor_with_nonnegative_diagonal();
    }

    void
    factor_first (const Ordinal nrows,
		  const Ordinal ncols,
		  Scalar A[],
		  const Ordinal lda,
		  Scalar tau[],
		  Scalar work[]) const
    {
      return default_.factor_first (nrows, ncols, A, lda, tau, work);
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
		 Scalar work[]) const
    {
      return default_.apply_first (applyType, nrows, ncols_C, ncols_A, 
				   A, lda, tau, C, ldc, work);
    }

    void
    apply_inner (const ApplyType& applyType,
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
		 Scalar work[]) const;

    void
    factor_inner (const Ordinal m,
		  const Ordinal n,
		  Scalar R[],
		  const Ordinal ldr,
		  Scalar A[],
		  const Ordinal lda,
		  Scalar tau[],
		  Scalar work[]) const;

    void
    factor_pair (const Ordinal n,
		 Scalar R_top[],
		 const Ordinal ldr_top,
		 Scalar R_bot[],
		 const Ordinal ldr_bot,
		 Scalar tau[],
		 Scalar work[]) const;

    void
    apply_pair (const ApplyType& applyType,
		const Ordinal ncols_C,
		const Ordinal ncols_Q,
		const Scalar R_bot[],
		const Ordinal ldr_bot,
		const Scalar tau[],
		Scalar C_top[],
		const Ordinal ldc_top,
		Scalar C_bot[],
		const Ordinal ldc_bot,
		Scalar work[]) const;

  private:
    mutable combine_default_type default_;
  };
		 

  /// "Forward declaration" for the real-arithmetic case.
  ///
  template< class Ordinal, class Scalar >
  class CombineNative< Ordinal, Scalar, false >
  {
  public:
    typedef Scalar scalar_type;
    typedef typename Teuchos::ScalarTraits< Scalar >::magnitudeType magnitude_type;
    typedef Ordinal ordinal_type;

  private:
    typedef BLAS< ordinal_type, scalar_type > blas_type;
    typedef LAPACK< ordinal_type, scalar_type > lapack_type;
    typedef CombineDefault< ordinal_type, scalar_type > combine_default_type;

  public:
    CombineNative () {}

    static bool QR_produces_R_factor_with_nonnegative_diagonal() { 
      return lapack_type::QR_produces_R_factor_with_nonnegative_diagonal() &&
	combine_default_type::QR_produces_R_factor_with_nonnegative_diagonal();
    }

    void
    factor_first (const Ordinal nrows,
		  const Ordinal ncols,
		  Scalar A[],
		  const Ordinal lda,
		  Scalar tau[],
		  Scalar work[]) const
    {
      return default_.factor_first (nrows, ncols, A, lda, tau, work);
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
		 Scalar work[]) const
    {
      return default_.apply_first (applyType, nrows, ncols_C, ncols_A, 
				   A, lda, tau, C, ldc, work);
    }

    void
    apply_inner (const ApplyType& applyType,
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
		 Scalar work[]) const;

    void
    factor_inner (const Ordinal m,
		  const Ordinal n,
		  Scalar R[],
		  const Ordinal ldr,
		  Scalar A[],
		  const Ordinal lda,
		  Scalar tau[],
		  Scalar work[]) const;

    void
    factor_pair (const Ordinal n,
		 Scalar R_top[],
		 const Ordinal ldr_top,
		 Scalar R_bot[],
		 const Ordinal ldr_bot,
		 Scalar tau[],
		 Scalar work[]) const;

    void
    apply_pair (const ApplyType& applyType,
		const Ordinal ncols_C,
		const Ordinal ncols_Q,
		const Scalar R_bot[],
		const Ordinal ldr_bot,
		const Scalar tau[],
		Scalar C_top[],
		const Ordinal ldc_top,
		Scalar C_bot[],
		const Ordinal ldc_bot,
		Scalar work[]) const;

  private:
    mutable combine_default_type default_;
  };


  /// "Forward declaration" for the complex-arithmetic case.
  ///
  template< class Ordinal, class Scalar >
  class CombineNative< Ordinal, Scalar, true >
  {
  public:
    typedef Scalar scalar_type;
    typedef typename Teuchos::ScalarTraits< Scalar >::magnitudeType magnitude_type;
    typedef Ordinal ordinal_type;

  private:
    typedef BLAS< ordinal_type, scalar_type > blas_type;
    typedef LAPACK< ordinal_type, scalar_type > lapack_type;
    typedef CombineDefault< ordinal_type, scalar_type > combine_default_type;

  public:
    CombineNative () {}

    static bool QR_produces_R_factor_with_nonnegative_diagonal() { 
      return lapack_type::QR_produces_R_factor_with_nonnegative_diagonal() &&
	combine_default_type::QR_produces_R_factor_with_nonnegative_diagonal();
    }

    void
    factor_first (const Ordinal nrows,
		  const Ordinal ncols,
		  Scalar A[],
		  const Ordinal lda,
		  Scalar tau[],
		  Scalar work[]) const
    {
      return default_.factor_first (nrows, ncols, A, lda, tau, work);
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
		 Scalar work[]) const
    {
      return default_.apply_first (applyType, nrows, ncols_C, ncols_A, 
				   A, lda, tau, C, ldc, work);
    }

    void
    apply_inner (const ApplyType& applyType,
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
		 Scalar work[]) const
    {
      return default_.apply_inner (applyType, m, ncols_C, ncols_Q,
				   A, lda, tau,
				   C_top, ldc_top, C_bot, ldc_bot,
				   work);
    }

    void
    factor_inner (const Ordinal m,
		  const Ordinal n,
		  Scalar R[],
		  const Ordinal ldr,
		  Scalar A[],
		  const Ordinal lda,
		  Scalar tau[],
		  Scalar work[]) const
    {
      return default_.factor_inner (m, n, R, ldr, A, lda, tau, work);
    }

    void
    factor_pair (const Ordinal n,
		 Scalar R_top[],
		 const Ordinal ldr_top,
		 Scalar R_bot[],
		 const Ordinal ldr_bot,
		 Scalar tau[],
		 Scalar work[]) const
    {
      return default_.factor_pair (n, R_top, ldr_top, R_bot, ldr_bot, tau, work);
    }

    void
    apply_pair (const ApplyType& applyType,
		const Ordinal ncols_C,
		const Ordinal ncols_Q,
		const Scalar R_bot[],
		const Ordinal ldr_bot,
		const Scalar tau[],
		Scalar C_top[],
		const Ordinal ldc_top,
		Scalar C_bot[],
		const Ordinal ldc_bot,
		Scalar work[]) const
    {
      return default_.apply_pair (applyType, ncols_C, ncols_Q,
				  R_bot, ldr_bot, tau,
				  C_top, ldc_top, C_bot, ldc_bot,
				  work);
    }

  private:
    mutable combine_default_type default_;
  };
		 

  template< class Ordinal, class Scalar >
  void
  CombineNative< Ordinal, Scalar, false >::
  factor_inner (const Ordinal m,
		const Ordinal n,
		Scalar R[],
		const Ordinal ldr,
		Scalar A[],
		const Ordinal lda,
		Scalar tau[],
		Scalar work[]) const
  {
    const Scalar ZERO(0), ONE(1);
    lapack_type lapack;
    blas_type blas;

    for (Ordinal k = 0; k < n; ++k)
      work[k] = ZERO;

    for (Ordinal k = 0; k < n-1; ++k)
      {
	Scalar& R_kk = R[ k + k * ldr ];
	Scalar* const A_1k = &A[ 0 + k * lda ];
	Scalar* const A_1kp1 = &A[ 0 + (k+1) * lda ];

	lapack.LARFP (m + 1, R_kk, A_1k, 1, tau[k]);
	blas.GEMV ("T", m, n-k-1, ONE, A_1kp1, lda, A_1k, 1, ZERO, work, 1);

	for (Ordinal j = k+1; j < n; ++j)
	  {
	    Scalar& R_kj = R[ k + j*ldr ];

	    work[j-k-1] += R_kj;
	    R_kj -= tau[k] * work[j-k-1];
	  }
	blas.GER (m, n-k-1, -tau[k], A_1k, 1, work, 1, A_1kp1, lda);
      }
    Scalar& R_nn = R[ (n-1) + (n-1) * ldr ];
    Scalar* const A_1n = &A[ 0 + (n-1) * lda ];

    lapack.LARFP (m+1, R_nn, A_1n, 1, tau[n-1]);
  }


  template< class Ordinal, class Scalar >
  void
  CombineNative< Ordinal, Scalar, false >::
  apply_inner (const ApplyType& applyType,
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
	       Scalar work[]) const
  {
    const Scalar ZERO(0);
    blas_type blas;

    //Scalar* const y = work;
    for (Ordinal i = 0; i < ncols_C; ++i)
      work[i] = ZERO;
    
    Ordinal j_start, j_end, j_step;
    if (applyType == ApplyType::NoTranspose)
      {
	j_start = ncols_Q - 1;
	j_end = -1; // exclusive
	j_step = -1;
      }
    else
      {
	j_start = 0; 
	j_end = ncols_Q; // exclusive
	j_step = +1;
      }
    for (Ordinal j = j_start; j != j_end; j += j_step)
      {
	const Scalar* const A_1j = &A[ 0 + j*lda ];

	//blas.GEMV ("T", m, ncols_C, ONE, C_bot, ldc_bot, A_1j, 1, ZERO, &y[0], 1);	
	for (Ordinal i = 0; i < ncols_C; ++i)
	  {
	    work[i] = ZERO;
	    for (Ordinal k = 0; k < m; ++k)
	      work[i] += A_1j[k] * C_bot[ k + i*ldc_bot ];

	    work[i] += C_top[ j + i*ldc_top ];
	  }
	for (Ordinal k = 0; k < ncols_C; ++k)
	  C_top[ j + k*ldc_top ] -= tau[j] * work[k];

	blas.GER (m, ncols_C, -tau[j], A_1j, 1, work, 1, C_bot, ldc_bot);
      }
  }


  template< class Ordinal, class Scalar >
  void
  CombineNative< Ordinal, Scalar, false >::
  factor_pair (const Ordinal n,
	       Scalar R_top[],
	       const Ordinal ldr_top,
	       Scalar R_bot[],
	       const Ordinal ldr_bot,
	       Scalar tau[],
	       Scalar work[]) const
  {
    const Scalar ZERO(0), ONE(1);
    lapack_type lapack;
    blas_type blas;

    for (Ordinal k = 0; k < n; ++k)
      work[k] = ZERO;

    for (Ordinal k = 0; k < n-1; ++k)
      {
	Scalar& R_top_kk = R_top[ k + k * ldr_top ];
	Scalar* const R_bot_1k = &R_bot[ 0 + k * ldr_bot ];
	Scalar* const R_bot_1kp1 = &R_bot[ 0 + (k+1) * ldr_bot ];

	// k+2: 1 element in R_top (R_top(k,k)), and k+1 elements in
	// R_bot (R_bot(1:k,k), in 1-based indexing notation).
	lapack.LARFP (k+2, R_top_kk, R_bot_1k, 1, tau[k]);
	// One-based indexing, Matlab version of the GEMV call below:
	// work(1:k) := R_bot(1:k,k+1:n)' * R_bot(1:k,k) 
	blas.GEMV ("T", k+1, n-k-1, ONE, R_bot_1kp1, ldr_bot, R_bot_1k, 1, ZERO, work, 1);
	
	for (Ordinal j = k+1; j < n; ++j)
	  {
	    Scalar& R_top_kj = R_top[ k + j*ldr_top ];
	    work[j-k-1] += R_top_kj;
	    R_top_kj -= tau[k] * work[j-k-1];
	  }
	blas.GER (k+1, n-k-1, -tau[k], R_bot_1k, 1, work, 1, R_bot_1kp1, ldr_bot);
      }
    Scalar& R_top_nn = R_top[ (n-1) + (n-1)*ldr_top ];
    Scalar* const R_bot_1n = &R_bot[ 0 + (n-1)*ldr_bot ];

    // n+1: 1 element in R_top (n,n), and n elements in R_bot (the
    // whole last column).
    lapack.LARFP (n+1, R_top_nn, R_bot_1n, 1, tau[n-1]);
  }


  template< class Ordinal, class Scalar >
  void
  CombineNative< Ordinal, Scalar, false >::
  apply_pair (const ApplyType& applyType,
	      const Ordinal ncols_C,
	      const Ordinal ncols_Q,
	      const Scalar R_bot[],
	      const Ordinal ldr_bot,
	      const Scalar tau[],
	      Scalar C_top[],
	      const Ordinal ldc_top,
	      Scalar C_bot[],
	      const Ordinal ldc_bot,
	      Scalar work[]) const
  {
    const Scalar ZERO(0);
    blas_type blas;

    for (Ordinal i = 0; i < ncols_C; ++i)
      work[i] = ZERO;
    
    Ordinal j_start, j_end, j_step;
    if (applyType == ApplyType::NoTranspose)
      {
	j_start = ncols_Q - 1;
	j_end = -1; // exclusive
	j_step = -1;
      }
    else
      {
	j_start = 0; 
	j_end = ncols_Q; // exclusive
	j_step = +1;
      }
    for (Ordinal j_Q = j_start; j_Q != j_end; j_Q += j_step)
      { // Using Householder reflector stored in column j_Q of R_bot
	const Scalar* const R_bot_col = &R_bot[ 0 + j_Q*ldr_bot ];

	// In 1-based indexing notation, with k in 1, 2, ..., ncols_C
	// (inclusive): (Output is length ncols_C row vector)
	//
	// work(1:j) := R_bot(1:j,j)' * C_bot(1:j, 1:ncols_C) - C_top(j, 1:ncols_C)
	for (Ordinal j_C = 0; j_C < ncols_C; ++j_C)
	  { // For each column j_C of [C_top; C_bot], update row j_Q
	    // of C_top and rows 1:j_Q of C_bot.  (Again, this is in
	    // 1-based indexing notation.

	    Scalar work_j_C = ZERO;
	    const Scalar* const C_bot_col = &C_bot[ 0 + j_C*ldc_bot ];

	    for (Ordinal k = 0; k <= j_Q; ++k)
	      work_j_C += R_bot_col[k] * C_bot_col[k];

	    work_j_C += C_top[ j_Q + j_C*ldc_top ];
	    work[j_C] = work_j_C;
	  }
	for (Ordinal j_C = 0; j_C < ncols_C; ++j_C)
	  C_top[ j_Q + j_C*ldc_top ] -= tau[j_Q] * work[j_C];

	blas.GER (j_Q+1, ncols_C, -tau[j_Q], R_bot_col, 1, work, 1, C_bot, ldc_bot);
      }
  }
} // namespace TSQR



#endif // __TSQR_CombineNative_hpp
