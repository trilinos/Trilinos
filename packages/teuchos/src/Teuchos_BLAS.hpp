// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// ***********************************************************************
// @HEADER

// Kris
// 06.16.03 -- Start over from scratch
// 06.16.03 -- Initial templatization (Tpetra_BLAS.cpp is no longer needed)
// 06.18.03 -- Changed xxxxx_() function calls to XXXXX_F77()
//          -- Added warning messages for generic calls
// 07.08.03 -- Move into Teuchos package/namespace
// 07.24.03 -- The first iteration of BLAS generics is nearing completion. Caveats:
//             * TRSM isn't finished yet; it works for one case at the moment (left side, upper tri., no transpose, no unit diag.)
//             * Many of the generic implementations are quite inefficient, ugly, or both. I wrote these to be easy to debug, not for efficiency or legibility. The next iteration will improve both of these aspects as much as possible.
//             * Very little verification of input parameters is done, save for the character-type arguments (TRANS, etc.) which is quite robust.
//             * All of the routines that make use of both an incx and incy parameter (which includes much of the L1 BLAS) are set up to work iff incx == incy && incx > 0. Allowing for differing/negative values of incx/incy should be relatively trivial.
//             * All of the L2/L3 routines assume that the entire matrix is being used (that is, if A is mxn, lda = m); they don't work on submatrices yet. This *should* be a reasonably trivial thing to fix, as well.
//          -- Removed warning messages for generic calls
// 08.08.03 -- TRSM now works for all cases where SIDE == L and DIAG == N. DIAG == U is implemented but does not work correctly; SIDE == R is not yet implemented.
// 08.14.03 -- TRSM now works for all cases and accepts (and uses) leading-dimension information.
// 09.26.03 -- character input replaced with enumerated input to cause compiling errors and not run-time errors ( suggested by RAB ).

#ifndef _TEUCHOS_BLAS_HPP_
#define _TEUCHOS_BLAS_HPP_

/*! \file Teuchos_BLAS.hpp
    \brief Templated interface class to BLAS routines.
*/

#include "Teuchos_BLAS_wrappers.hpp"
#include "Teuchos_BLAS_types.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_OrdinalTraits.hpp"

/*! \class Teuchos::BLAS
    \brief The Templated BLAS Wrapper Class.

    The Teuchos::BLAS class provides functionality similar to the BLAS
    (Basic Linear Algebra Subprograms).  The BLAS provide portable, high-
    performance implementations of kernels such as dense vector multiplication,
    dot products, dense matrix-vector multiplication and dense matrix-matrix
    multiplication.

    The standard BLAS interface is Fortran-specific.  Unfortunately, the 
    interface between C++ and Fortran is not standard across all computer
    platforms.  The Teuchos_BLAS class provides C++ bindings for the BLAS
    kernels in order to insulate the rest of Petra from the details of 
    C++ to Fortran translation.

    In addition to giving access the standard BLAS functionality.
    Teuchos::BLAS also provide functionality for any <ScalarType> class that
    defines the +, - * and / operators.

    Teuchos::BLAS is a single memory image interface only.  This is appropriate 
    since the standard BLAS are only specified for serial execution 
    (or shared memory parallel).

    \note
    <ol>
            <li>These templates are specialized to use the Fortran BLAS routines for
            scalar types \c float and \c double.

            <li>If Teuchos is configured with \c --enable-teuchos-complex then these templates
            are specialized for scalar types \c complex<float> and \c complex<double> also.
    </ol>
*/

namespace Teuchos
{
  extern const char ESideChar[];
  extern const char ETranspChar[];
  extern const char EUploChar[];
  extern const char EDiagChar[];

  template<typename OrdinalType, typename ScalarType>
  class BLAS
  {    
  public:
    //@{ \name Constructor/Destructor.

    //! Default constructor.
    inline BLAS(void) {};

    //! Copy constructor.
    inline BLAS(const BLAS& BLAS_source) {};

    //! Destructor.
    inline virtual ~BLAS(void) {};
    //@}

    //@{ \name Level 1 BLAS Routines.

    //! Computes a Givens plane rotation.
    void ROTG(ScalarType* da, ScalarType* db, ScalarType* c, ScalarType* s) const;

    //! Scale the vector \c x by the constant \c alpha.
    void SCAL(const OrdinalType n, const ScalarType alpha, ScalarType* x, const OrdinalType incx) const;

    //! Copy the vector \c x to the vector \c y.
    void COPY(const OrdinalType n, const ScalarType* x, const OrdinalType incx, ScalarType* y, const OrdinalType incy) const;

    //! Perform the operation: \c y \c <- \c y+alpha*x.
    void AXPY(const OrdinalType n, const ScalarType alpha, const ScalarType* x, const OrdinalType incx, ScalarType* y, const OrdinalType incy) const;

    //! Sum the absolute values of the entries of \c x.
    ScalarType ASUM(const OrdinalType n, const ScalarType* x, const OrdinalType incx) const;

    //! Form the dot product of the vectors \c x and \c y.
    ScalarType DOT(const OrdinalType n, const ScalarType* x, const OrdinalType incx, const ScalarType* y, const OrdinalType incy) const;

    //! Compute the 2-norm of the vector \c x.
    ScalarType NRM2(const OrdinalType n, const ScalarType* x, const OrdinalType incx) const;

    //! Return the index of the element of \c x with the maximum magnitude.
    OrdinalType IAMAX(const OrdinalType n, const ScalarType* x, const OrdinalType incx) const;

    //@}

    //@{ \name Level 2 BLAS Routines.

    //! Performs the matrix-vector operation:  \c y \c <- \c alpha*A*x+beta*y or \c y \c <- \c alpha*A'*x+beta*y where \c A is a general \c m by \c n matrix.
    void GEMV(ETransp trans, const OrdinalType m, const OrdinalType n, const ScalarType alpha, const ScalarType* A, 
	      const OrdinalType lda, const ScalarType* x, const OrdinalType incx, const ScalarType beta, ScalarType* y, const OrdinalType incy) const;

    //! Performs the matrix-vector operation:  \c x \c <- \c A*x or \c x \c <- \c A'*x where \c A is a unit/non-unit \c n by \c n upper/lower triangular matrix.
    void TRMV(EUplo uplo, ETransp trans, EDiag diag, const OrdinalType n, const ScalarType* A, 
	      const OrdinalType lda, ScalarType* x, const OrdinalType incx) const;

    //! Performs the rank 1 operation:  \c A \c <- \c alpha*x*y'+A.
    void GER(const OrdinalType m, const OrdinalType n, const ScalarType alpha, const ScalarType* x, const OrdinalType incx, 
	     const ScalarType* y, const OrdinalType incy, ScalarType* A, const OrdinalType lda) const;
    //@}
    
    //@{ \name Level 3 BLAS Routines. 

    //! Performs the matrix-matrix operation: \c C \c <- \c alpha*op(A)*op(B)+beta*C where \c op(A) is either \c A or \c A', \c op(B) is either \c B or \c B', and C is an \c m by \c k matrix.
    void GEMM(ETransp transa, ETransp transb, const OrdinalType m, const OrdinalType n, const OrdinalType k, const ScalarType alpha, const ScalarType* A, const OrdinalType lda, const ScalarType* B, const OrdinalType ldb, const ScalarType beta, ScalarType* C, const OrdinalType ldc) const;

    //! Performs the matrix-matrix operation: \c C \c <- \c alpha*A*B+beta*C or \c C \c <- \c alpha*B*A+beta*C where \c A is an \c m by \c m or \c n by \c n symmetric matrix and \c B is a general matrix.
    void SYMM(ESide side, EUplo uplo, const OrdinalType m, const OrdinalType n, const ScalarType alpha, const ScalarType* A, const OrdinalType lda, const ScalarType* B, const OrdinalType ldb, const ScalarType beta, ScalarType* C, const OrdinalType ldc) const;

    //! Performs the matrix-matrix operation: \c C \c <- \c alpha*op(A)*B+beta*C or \c C \c <- \c alpha*B*op(A)+beta*C where \c op(A) is an unit/non-unit, upper/lower triangular matrix and \c B is a general matrix.
    void TRMM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const OrdinalType m, const OrdinalType n,
                const ScalarType alpha, const ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb) const;

    //! Solves the matrix equations:  \c op(A)*X=alpha*B or \c X*op(A)=alpha*B where \c X and \c B are \c m by \c n matrices, \c A is a unit/non-unit, upper/lower triangular matrix and \c op(A) is \c A or \c A'.  The matrix \c X is overwritten on \c B.
    void TRSM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const OrdinalType m, const OrdinalType n,
                const ScalarType alpha, const ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb) const;
    //@}
  };

//------------------------------------------------------------------------------------------
//      LEVEL 1 BLAS ROUTINES  
//------------------------------------------------------------------------------------------
    
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::ROTG(ScalarType* da, ScalarType* db, ScalarType* c, ScalarType* s) const
  {
    ScalarType roe, scale, r;
    ScalarType zero = ScalarTraits<ScalarType>::zero();
    ScalarType one = ScalarTraits<ScalarType>::one();

    if ( ScalarTraits<ScalarType>::magnitude( *da ) > ScalarTraits<ScalarType>::magnitude( *db ) ) { roe = *da; }
    scale = ScalarTraits<ScalarType>::magnitude( *da ) + ScalarTraits<ScalarType>::magnitude( *db );
    if ( scale == zero ) // There is nothing to do.
    {
      *c = one;
      *s = zero;
      *da = zero; *db = zero;
    } else { // Compute the Givens rotation.
      r = scale*ScalarTraits<ScalarType>::squareroot( ( *da/scale)*(*da/scale) + (*db/scale)*(*db/scale) );
      if ( roe < zero ) { r *= -one; }
      *c = *da / r;
      *s = *db / r;
      *db = ScalarTraits<ScalarType>::one();
      if( ScalarTraits<ScalarType>::magnitude( *da ) > ScalarTraits<ScalarType>::magnitude( *db ) ){ *db = *s; }
      if( ScalarTraits<ScalarType>::magnitude( *db ) >= ScalarTraits<ScalarType>::magnitude( *da ) &&
	   *c != ScalarTraits<ScalarType>::zero() ) { *db = one / *c; }
      *da = r;
    }
  } /* end ROTG */
      
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::SCAL(const OrdinalType n, const ScalarType alpha, ScalarType* x, const OrdinalType incx) const
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    OrdinalType i, ix = izero;
    if ( n > izero ) {
        // Set the initial index (ix).
        if (incx < izero) { ix = (-n+ione)*incx; } 
        // Scale the vector.
        for(i = izero; i < n; i++)
        {
            x[ix] *= alpha;
            ix += incx;
        }
    }
  } /* end SCAL */

  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::COPY(const OrdinalType n, const ScalarType* x, const OrdinalType incx, ScalarType* y, const OrdinalType incy) const
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    OrdinalType i, ix = izero, iy = izero;
    if ( n > izero ) {
	// Set the initial indices (ix, iy).
    	if (incx < izero) { ix = (-n+ione)*incx; }
    	if (incy < izero) { iy = (-n+ione)*incy; }

        for(i = izero; i < n; i++)
          {
	    y[iy] = x[ix];
	    ix += incx;
	    iy += incy;
          }
      }
  } /* end COPY */

  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::AXPY(const OrdinalType n, const ScalarType alpha, const ScalarType* x, const OrdinalType incx, ScalarType* y, const OrdinalType incy) const
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    OrdinalType i, ix = izero, iy = izero;
    if( n > izero && alpha != ScalarTraits<ScalarType>::zero())
      {
	// Set the initial indices (ix, iy).
    	if (incx < izero) { ix = (-n+ione)*incx; }
    	if (incy < izero) { iy = (-n+ione)*incy; }

        for(i = izero; i < n; i++)
          {
	    y[iy] += alpha * x[ix];
	    ix += incx;
	    iy += incy;
          }
      }
  } /* end AXPY */

  template<typename OrdinalType, typename ScalarType>
  ScalarType BLAS<OrdinalType, ScalarType>::ASUM(const OrdinalType n, const ScalarType* x, const OrdinalType incx) const
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    ScalarType result = ScalarTraits<ScalarType>::zero();
    OrdinalType i, ix = izero;
    if( n > izero ) {
	// Set the initial indices
	if (incx < izero) { ix = (-n+ione)*incx; }

    	for(i = izero; i < n; i++)
          {
	    result += ScalarTraits<ScalarType>::magnitude(x[ix]);
	    ix += incx;
          }
    } 
   return result;
  } /* end ASUM */
  
  template<typename OrdinalType, typename ScalarType>
  ScalarType BLAS<OrdinalType, ScalarType>::DOT(const OrdinalType n, const ScalarType* x, const OrdinalType incx, const ScalarType* y, const OrdinalType incy) const
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    ScalarType result = ScalarTraits<ScalarType>::zero();
    OrdinalType i, ix = izero, iy = izero;
    if( n > izero ) 
      {
	// Set the initial indices (ix,iy). 	    
	if (incx < izero) { ix = (-n+ione)*incx; }
	if (incy < izero) { iy = (-n+ione)*incy; }

	for(i = izero; i < n; i++)
	  {
	    result += x[ix] * y[iy];
	    ix += incx;
	    iy += incy;
	  }
      }
    return result;
  } /* end DOT */
  
  template<typename OrdinalType, typename ScalarType>
  ScalarType BLAS<OrdinalType, ScalarType>::NRM2(const OrdinalType n, const ScalarType* x, const OrdinalType incx) const
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    ScalarType result = ScalarTraits<ScalarType>::zero();
    OrdinalType i, ix = izero;
    if ( n > izero ) 
      {
	// Set the initial index.
	if (incx < izero) { ix = (-n+ione)*incx; }	
    
	for(i = izero; i < n; i++)
      	  {
	    result += x[ix] * x[ix];
	    ix += incx;
       	  }
    	result = ScalarTraits<ScalarType>::squareroot(result);
      }	
    return result;
  } /* end NRM2 */
  
  template<typename OrdinalType, typename ScalarType>
  OrdinalType BLAS<OrdinalType, ScalarType>::IAMAX(const OrdinalType n, const ScalarType* x, const OrdinalType incx) const
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    OrdinalType result = izero, ix = izero, i;
    ScalarType maxval;

    if ( n > izero ) 
      {
	if (incx < izero) { ix = (-n+ione)*incx; }
	maxval = ScalarTraits<ScalarType>::magnitude(x[ix]);
	ix += incx;
    	for(i = ione; i < n; i++)
      	  {
	    if(ScalarTraits<ScalarType>::magnitude(x[ix]) > maxval)
	      {
	    	result = i;
	        maxval = ScalarTraits<ScalarType>::magnitude(x[ix]);
	      }
	    ix += incx;
	  }
      }
    return result + 1; // the BLAS I?AMAX functions return 1-indexed (Fortran-style) values
  } /* end IAMAX */

//------------------------------------------------------------------------------------------
//      LEVEL 2 BLAS ROUTINES
//------------------------------------------------------------------------------------------

  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::GEMV(ETransp trans, const OrdinalType m, const OrdinalType n, const ScalarType alpha, const ScalarType* A, const OrdinalType lda, const ScalarType* x, const OrdinalType incx, const ScalarType beta, ScalarType* y, const OrdinalType incy) const
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    ScalarType zero = ScalarTraits<ScalarType>::zero();
    ScalarType one = ScalarTraits<ScalarType>::one();
    bool BadArgument = false;

    // Quick return if there is nothing to do!
    if( m == izero || n == izero || ( alpha == zero && beta == one ) ){ return; }
    
    // Otherwise, we need to check the argument list.
    if( m < izero ) { 
	cout << "BLAS::GEMV Error: M == " << m << endl;	    
	BadArgument = true;
    }
    if( n < izero ) { 
	cout << "BLAS::GEMV Error: N == " << n << endl;	    
	BadArgument = true;
    }
    if( lda < m ) { 
	cout << "BLAS::GEMV Error: LDA < MAX(1,M)"<< endl;	    
	BadArgument = true;
    }
    if( incx == izero ) {
	cout << "BLAS::GEMV Error: INCX == 0"<< endl;
	BadArgument = true;
    }
    if( incy == izero ) {
	cout << "BLAS::GEMV Error: INCY == 0"<< endl;
	BadArgument = true;
    }

    if(!BadArgument) {
      OrdinalType i, j, lenx, leny, ix, iy, jx, jy; 
      OrdinalType kx = izero, ky = izero;
      ScalarType temp;

      // Determine the lengths of the vectors x and y.
      if(ETranspChar[trans] == 'N') {
	lenx = n;
	leny = m;
      } else {
	lenx = m;
	leny = n;
      }

      // Set the starting pointers for the vectors x and y if incx/y < 0.
      if (incx < izero ) { kx =  (ione - lenx)*incx; }
      if (incy < izero ) { ky =  (ione - leny)*incy; }

      // Form y = beta*y
      ix = kx; iy = ky;
      if(beta != one) {
	if (incy == ione) {
	  if (beta == zero) {
	    for(i = izero; i < leny; i++) { y[i] = zero; }
	  } else {
	    for(i = izero; i < leny; i++) { y[i] *= beta; }
	  }
	} else {
	  if (beta == zero) {
	    for(i = izero; i < leny; i++) {
	      y[iy] = zero;
	      iy += incy;
	    }
	  } else {
	    for(i = izero; i < leny; i++) {
	      y[iy] *= beta;
	      iy += incy;
	    }
	  }
	}
      }
	
      // Return if we don't have to do anything more.
      if(alpha == zero) { return; }

      if( ETranspChar[trans] == 'N' ) {
	// Form y = alpha*A*y
	jx = kx;
	if (incy == ione) {
	  for(j = izero; j < n; j++) {
	    if (x[jx] != zero) {
	      temp = alpha*x[jx];
	      for(i = izero; i < m; i++) {
		y[i] += temp*A[j*lda + i];
	      }
	    }
	    jx += incx;
	  }
	} else {
	  for(j = izero; j < n; j++) {
	    if (x[jx] != zero) {
	      temp = alpha*x[jx];
	      iy = ky;
	      for(i = izero; i < m; i++) {
		y[iy] += temp*A[j*lda + i];
		iy += incy;
	      }
	    }
	    jx += incx;
	  }
	}
      } else {
	jy = ky;
	if (incx == ione) {
	  for(j = izero; j < n; j++) {
	    temp = zero;
	    for(i = izero; i < m; i++) {
	      temp += A[j*lda + i]*x[i];
	    }
	    y[jy] += alpha*temp;
	    jy += incy;
	  }
	} else {
	  for(j = izero; j < n; j++) {
	    temp = zero;
	    ix = kx;
	    for (i = izero; i < m; i++) {
	      temp += A[j*lda + i]*x[ix];
	      ix += incx;
	    }
	    y[jy] += alpha*temp;
	    jy += incy;
	  }
	}
      }
    } /* if (!BadArgument) */
  } /* end GEMV */

 template<typename OrdinalType, typename ScalarType>
 void BLAS<OrdinalType, ScalarType>::TRMV(EUplo uplo, ETransp trans, EDiag diag, const OrdinalType n, const ScalarType* A, const OrdinalType lda, ScalarType* x, const OrdinalType incx) const
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    ScalarType zero = ScalarTraits<ScalarType>::zero();
    ScalarType one = ScalarTraits<ScalarType>::one();
    bool BadArgument = false;

    // Quick return if there is nothing to do!
    if( n == izero ){ return; }
    
    // Otherwise, we need to check the argument list.
    if( n < izero ) { 
      cout << "BLAS::TRMV Error: N == " << n << endl;	    
      BadArgument = true;
    }
    if( lda < n ) { 
      cout << "BLAS::TRMV Error: LDA < MAX(1,N)"<< endl;	    
      BadArgument = true;
    }
    if( incx == izero ) {
      cout << "BLAS::TRMV Error: INCX == 0"<< endl;
      BadArgument = true;
    }

    if(!BadArgument) {
      OrdinalType i, j, ix, jx, kx = izero;
      ScalarType temp;
      bool NoUnit = (EDiagChar[diag] == 'N');

      // Set the starting pointer for the vector x if incx < 0.
      if (incx < izero) { kx = (-n+ione)*incx; }

      // Start the operations for a nontransposed triangular matrix 
      if (ETranspChar[trans] == 'N') {
	/* Compute x = A*x */
	if (EUploChar[uplo] == 'U') {
	  /* A is an upper triangular matrix */
	  if (incx == ione) {
	    for (j=izero; j<n; j++) {
	      if (x[j] != zero) {
		temp = x[j];
		for (i=izero; i<j; i++) {
		  x[i] += temp*A[j*lda + i];
		}
		if (NoUnit) 
		  x[j] *= A[j*lda + j];
	      }
	    }
	  } else {
	    jx = kx;
	    for (j=izero; j<n; j++) {
	      if (x[jx] != zero) {
		temp = x[jx];
		ix = kx;
		for (i=izero; i<j; i++) {
		  x[ix] += temp*A[j*lda + i];
		  ix += incx;
		}
		if (NoUnit)
		  x[jx] *= A[j*lda + j];
	      }
	      jx += incx;
	    }
	  } /* if (incx == ione) */
	} else { /* A is a lower triangular matrix */
	  if (incx == ione) {
	    for (j=n-ione; j>-ione; j--) {
	      if (x[j] != zero) {
		temp = x[j];
		for (i=n-ione; i>j; i--) {
		  x[i] += temp*A[j*lda + i];
		}
		if (NoUnit)
		  x[j] *= A[j*lda + j];
	      }
	    }
	  } else {
	    kx += (n-ione)*incx;
	    jx = kx;
	    for (j=n-ione; j>-ione; j--) {
	      if (x[jx] != zero) {
		temp = x[jx];
		ix = kx;
		for (i=n-ione; i>j; i--) {
		  x[ix] += temp*A[j*lda + i];
		  ix -= incx;
		}
		if (NoUnit) 
		  x[jx] *= A[j*lda + j];
	      }
	      jx -= incx;
	    }
	  }
	} /* if (EUploChar[uplo]=='U') */
      } else { /* A is transposed/conjugated */
	/* Compute x = A'*x */
	if (EUploChar[uplo]=='U') {
	  /* A is an upper triangular matrix */
	  if (incx == ione) {
	    for (j=n-ione; j>-ione; j--) {
	      temp = x[j];
	      if (NoUnit)
		temp *= A[j*lda + j];
	      for (i=j-ione; i>-ione; i--) {
		temp += A[j*lda + i]*x[i];
	      }
	      x[j] = temp;
	    }
	  } else {
	    jx = kx + (n-ione)*incx;
	    for (j=n-ione; j>-ione; j--) {
	      temp = x[jx];
	      ix = jx;
	      if (NoUnit)
		temp *= A[j*lda + j];
	      for (i=j-ione; i>-ione; i--) {
		ix -= incx;
		temp += A[j*lda + i]*x[ix];
	      }
	      x[jx] = temp;
	      jx -= incx;
	    }
	  }
	} else {
	  /* A is a lower triangular matrix */
	  if (incx == ione) {
	    for (j=izero; j<n; j++) {
	      temp = x[j];
	      if (NoUnit)
		temp *= A[j*lda + j];
	      for (i=j+ione; i<n; i++) {
		temp += A[j*lda + i]*x[i];
	      }
	      x[j] = temp;
	    }
	  } else {
	    jx = kx;
	    for (j=izero; j<n; j++) {
	      temp = x[jx];
	      ix = jx;
	      if (NoUnit) 
		temp *= A[j*lda + j];
	      for (i=j+ione; i<n; i++) {
		ix += incx;
		temp += A[j*lda + i]*x[ix];
	      }
	      x[jx] = temp;
	      jx += incx;	      
	    }
	  }
	} /* if (EUploChar[uplo]=='U') */
      } /* if (ETranspChar[trans]=='N') */
    } /* if (!BadArgument) */
  } /* end TRMV */
        
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::GER(const OrdinalType m, const OrdinalType n, const ScalarType alpha, const ScalarType* x, const OrdinalType incx, const ScalarType* y, const OrdinalType incy, ScalarType* A, const OrdinalType lda) const
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    ScalarType zero = ScalarTraits<ScalarType>::zero();
    ScalarType one = ScalarTraits<ScalarType>::one();
    bool BadArgument = false;

    // Quick return if there is nothing to do!
    if( m == izero || n == izero || alpha == zero ){ return; }
    
    // Otherwise, we need to check the argument list.
    if( m < izero ) { 
	cout << "BLAS::GER Error: M == " << m << endl;	    
	BadArgument = true;
    }
    if( n < izero ) { 
	cout << "BLAS::GER Error: N == " << n << endl;	    
	BadArgument = true;
    }
    if( lda < m ) { 
	cout << "BLAS::GER Error: LDA < MAX(1,M)"<< endl;	    
	BadArgument = true;
    }
    if( incx == 0 ) {
	cout << "BLAS::GER Error: INCX == 0"<< endl;
	BadArgument = true;
    }
    if( incy == 0 ) {
	cout << "BLAS::GER Error: INCY == 0"<< endl;
	BadArgument = true;
    }

    if(!BadArgument) {
      OrdinalType i, j, ix, jy = izero, kx = izero;
      ScalarType temp;

      // Set the starting pointers for the vectors x and y if incx/y < 0.
      if (incx < izero) { kx = (-m+ione)*incx; }
      if (incy < izero) { jy = (-n+ione)*incy; }

      // Start the operations for incx == 1
      if( incx == ione ) {
	for( j=izero; j<n; j++ ) {
	  if ( y[jy] != zero ) {
	    temp = alpha*y[jy];
	    for ( i=izero; i<m; i++ ) {
	      A[j*lda + i] += x[i]*temp;
	    }
	  }
	  jy += incy;
	}
      } 
      else { // Start the operations for incx != 1
	for( j=izero; j<n; j++ ) {
	  if ( y[jy] != zero ) {
	    temp = alpha*y[jy];
	    ix = kx;
	    for( i=izero; i<m; i++ ) {
	      A[j*lda + i] += x[ix]*temp;
	      ix += incx;
	    }
	  }
	  jy += incy;
	}
      }
    } /* if(!BadArgument) */
  } /* end GER */
  
//------------------------------------------------------------------------------------------
//      LEVEL 3 BLAS ROUTINES
//------------------------------------------------------------------------------------------
        
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::GEMM(ETransp transa, ETransp transb, const OrdinalType m, const OrdinalType n, const OrdinalType k, const ScalarType alpha, const ScalarType* A, const OrdinalType lda, const ScalarType* B, const OrdinalType ldb, const ScalarType beta, ScalarType* C, const OrdinalType ldc) const
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    ScalarType zero = ScalarTraits<ScalarType>::zero();
    ScalarType one = ScalarTraits<ScalarType>::one();
    OrdinalType i, j, p;
    OrdinalType NRowA = m, NRowB = k;
    ScalarType temp;
    bool BadArgument = false;

    // Change dimensions of matrix if either matrix is transposed
    if( !(ETranspChar[transa]=='N') ) {
      NRowA = k;
    }
    if( !(ETranspChar[transb]=='N') ) {
      NRowB = n;
    }

    // Quick return if there is nothing to do!
    if( (m==izero) || (n==izero) || (((alpha==zero)||(k==izero)) && (beta==one)) ){ return; }
    if( m < izero ) { 
      cout << "BLAS::GEMM Error: M == " << m << endl;	    
      BadArgument = true;
    }
    if( n < izero ) { 
      cout << "BLAS::GEMM Error: N == " << n << endl;	    
      BadArgument = true;
    }
    if( k < izero ) { 
      cout << "BLAS::GEMM Error: K == " << k << endl;	    
      BadArgument = true;
    }
    if( lda < NRowA ) { 
      cout << "BLAS::GEMM Error: LDA < MAX(1,M)"<< endl;	    
      BadArgument = true;
    }
    if( ldb < NRowB ) { 
      cout << "BLAS::GEMM Error: LDB < MAX(1,K)"<< endl;	    
      BadArgument = true;
    }
     if( ldc < m ) { 
      cout << "BLAS::GEMM Error: LDC < MAX(1,M)"<< endl;	    
      BadArgument = true;
    }

    if(!BadArgument) {

      // Only need to scale the resulting matrix C.
      if( alpha == zero ) {
	if( beta == zero ) {
	  for (j=izero; j<n; j++) {
	    for (i=izero; i<m; i++) {
	      C[j*ldc + i] = zero;
	    }
	  }
	} else {
	  for (j=izero; j<n; j++) {
	    for (i=izero; i<m; i++) {
	      C[j*ldc + i] *= beta;
	    }
	  }
	}
	return;
      }
      //
      // Now start the operations.
      //
      if ( ETranspChar[transb]=='N' ) {
	if ( ETranspChar[transa]=='N' ) {
	  // Compute C = alpha*A*B + beta*C
	  for (j=izero; j<n; j++) {
	    if( beta == zero ) {
	      for (i=izero; i<m; i++){
		C[j*ldc + i] = zero;
	      }
	    } else if( beta != one ) {
	      for (i=izero; i<m; i++){
		C[j*ldc + i] *= beta;
	      }
	    }
	    for (p=izero; p<k; p++){
	      if (B[j*ldb + p] != zero ){
		temp = alpha*B[j*ldb + p];
		for (i=izero; i<m; i++) {
		  C[j*ldc + i] += temp*A[p*lda + i];
		}
	      }
	    }
	  }
	} else {
	  // Compute C = alpha*A'*B + beta*C
	  for (j=izero; j<n; j++) {
	    for (i=izero; i<m; i++) {
	      temp = zero;
	      for (p=izero; p<k; p++) {
		temp += A[i*lda + p]*B[j*ldb + p];
	      }
	      if (beta == zero) {
		C[j*ldc + i] = alpha*temp;
	      } else {
		C[j*ldc + i] = alpha*temp + beta*C[j*ldc + i];
	      }
	    }
	  }
	}
      } else {
	if ( ETranspChar[transa]=='N' ) {
	  // Compute C = alpha*A*B' + beta*C
	  for (j=izero; j<n; j++) {
	    if (beta == zero) {
	      for (i=izero; i<m; i++) {
		C[j*ldc + i] = zero;
	      } 
	    } else if ( beta != one ) {
	      for (i=izero; i<m; i++) {
		C[j*ldc + i] *= beta;
	      }
	    }
	    for (p=izero; p<k; p++) {
	      if (B[p*ldb + j] != zero) {
		temp = alpha*B[p*ldb + j];
		for (i=izero; i<m; i++) {
		  C[j*ldc + i] += temp*A[p*lda + i];
		}
	      }
	    }
	  }
	} else {
	  // Compute C += alpha*A'*B' + beta*C
	  for (j=izero; j<n; j++) {
	    for (i=izero; i<m; i++) {
	      temp = zero;
	      for (p=izero; p<k; p++) {
		temp += A[i*lda + p]*B[p*ldb + j];
	      }
	      if (beta == zero) {
		C[j*ldc + i] = alpha*temp;
	      } else {
		C[j*ldc + i] = alpha*temp + beta*C[j*ldc + i];
	      }
	    }
	  }
	} // end if (ETranspChar[transa]=='N') ...
      } // end if (ETranspChar[transb]=='N') ...
    } // end if (!BadArgument) ...
  } // end of GEMM


  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::SYMM(ESide side, EUplo uplo, const OrdinalType m, const OrdinalType n, const ScalarType alpha, const ScalarType* A, const OrdinalType lda, const ScalarType* B, const OrdinalType ldb, const ScalarType beta, ScalarType* C, const OrdinalType ldc) const
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    ScalarType zero = ScalarTraits<ScalarType>::zero();
    ScalarType one = ScalarTraits<ScalarType>::one();
    OrdinalType i, j, k, NRowA = m;
    ScalarType temp1, temp2;
    bool BadArgument = false;
    bool Upper = (EUploChar[uplo] == 'U');
    if (ESideChar[side] == 'R') { NRowA = n; }
    
    // Quick return.
    if ( (m==izero) || (n==izero) || ( (alpha==zero)&&(beta==one) ) ) { return; }
    if( m < 0 ) { 
      cout << "BLAS::SYMM Error: M == "<< m << endl;
      BadArgument = true; }
    if( n < 0 ) {
      cout << "BLAS::SYMM Error: N == "<< n << endl;
      BadArgument = true; }
    if( lda < NRowA ) {
      cout << "BLAS::SYMM Error: LDA == "<<lda<<endl;
      BadArgument = true; }
    if( ldb < m ) {
      cout << "BLAS::SYMM Error: LDB == "<<ldb<<endl;
      BadArgument = true; }
    if( ldc < m ) {
      cout << "BLAS::SYMM Error: LDC == "<<ldc<<endl;
      BadArgument = true; }

    if(!BadArgument) {

      // Only need to scale C and return.
      if (alpha == zero) {
	if (beta == zero ) {
	  for (j=izero; j<n; j++) {
	    for (i=izero; i<m; i++) {
	      C[j*ldc + i] = zero;
	    }
	  }
	} else {
	  for (j=izero; j<n; j++) {
	    for (i=izero; i<m; i++) {
	      C[j*ldc + i] *= beta;
	    }
	  }
	}
	return;
      }

      if ( ESideChar[side] == 'L') {
	// Compute C = alpha*A*B + beta*C

	if (Upper) {
	  // The symmetric part of A is stored in the upper triangular part of the matrix.
	  for (j=izero; j<n; j++) {
	    for (i=izero; i<m; i++) {
	      temp1 = alpha*B[j*ldb + i];
	      temp2 = zero;
	      for (k=izero; k<i; k++) {
		C[j*ldc + k] += temp1*A[i*lda + k];
		temp2 += B[j*ldb + k]*A[i*lda + k];
	      }
	      if (beta == zero) {
		C[j*ldc + i] = temp1*A[i*lda + i] + alpha*temp2;
	      } else {
		C[j*ldc + i] = beta*C[j*ldc + i] + temp1*A[i*lda + i] + alpha*temp2;
	      }
	    }
	  }
	} else {
	  // The symmetric part of A is stored in the lower triangular part of the matrix.
	  for (j=izero; j<n; j++) {
	    for (i=m-ione; i>-ione; i--) {
	      temp1 = alpha*B[j*ldb + i];
	      temp2 = zero;
	      for (k=i+ione; k<m; k++) {
		C[j*ldc + k] += temp1*A[i*lda + k];
		temp2 += B[j*ldb + k]*A[i*lda + k];
	      }
	      if (beta == zero) {
		C[j*ldc + i] = temp1*A[i*lda + i] + alpha*temp2;
	      } else {
		C[j*ldc + i] = beta*C[j*ldc + i] + temp1*A[i*lda + i] + alpha*temp2;
	      }
	    }
	  }
	}
      } else {
	// Compute C = alpha*B*A + beta*C.
	for (j=izero; j<n; j++) {
	  temp1 = alpha*A[j*lda + j];
	  if (beta == zero) {
	    for (i=izero; i<m; i++) {
	      C[j*ldc + i] = temp1*B[j*ldb + i];
	    }
	  } else {
	    for (i=izero; i<m; i++) {
	      C[j*ldc + i] = beta*C[j*ldc + i] + temp1*B[j*ldb + i];
	    }
	  }
	  for (k=izero; k<j; k++) {
	    if (Upper) {
	      temp1 = alpha*A[j*lda + k];
	    } else {
	      temp1 = alpha*A[k*lda + j];
	    }
	    for (i=izero; i<m; i++) {
	      C[j*ldc + i] += temp1*B[k*ldb + i];
	    }
	  }
	  for (k=j+ione; k<n; k++) {
	    if (Upper) {
	      temp1 = alpha*A[k*lda + j];
	    } else {
	      temp1 = alpha*A[j*lda + k];
	    }
	    for (i=izero; i<m; i++) {
	      C[j*ldc + i] += temp1*B[k*ldb + i];
	    }
	  }
	}
      } // end if (ESideChar[side]=='L') ...
    } // end if(!BadArgument) ...
} // end SYMM
  
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::TRMM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const OrdinalType m, const OrdinalType n, const ScalarType alpha, const ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb) const
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    ScalarType zero = ScalarTraits<ScalarType>::zero();
    ScalarType one = ScalarTraits<ScalarType>::one();
    OrdinalType i, j, k, NRowA = m;
    ScalarType temp;
    bool BadArgument = false;
    bool LSide = (ESideChar[side] == 'L');
    bool NoUnit = (EDiagChar[diag] == 'N');
    bool Upper = (EUploChar[uplo] == 'U');

    if(!LSide) { NRowA = n; }

    // Quick return.
    if (n==izero || m==izero) { return; }
    if( m < 0 ) {
      cout << "BLAS::TRMM Error: M == "<< m <<endl;
      BadArgument = true; }
    if( n < 0 ) {
      cout << "BLAS::TRMM Error: N == "<< n <<endl;
      BadArgument = true; }
    if( lda < NRowA ) {
      cout << "BLAS::TRMM Error: LDA == "<< lda << endl;
      BadArgument = true; }
    if( ldb < m ) {
      cout << "BLAS::TRMM Error: M == "<< ldb << endl;
      BadArgument = true; }

    if(!BadArgument) {

      // B only needs to be zeroed out.
      if( alpha == zero ) {
	for( j=izero; j<n; j++ ) {
	  for( i=izero; i<m; i++ ) {
	    B[j*ldb + i] = zero;
	  }
	}
	return;
      }
      
      //  Start the computations. 
      if ( LSide ) {
	// A is on the left side of B.
	
	if ( ETranspChar[transa]=='N' ) {
	  // Compute B = alpha*A*B

	  if ( Upper ) {
	    // A is upper triangular
	    for( j=izero; j<n; j++ ) {
	      for( k=izero; k<m; k++) {
		if ( B[j*ldb + k] != zero ) {
		  temp = alpha*B[j*ldb + k];
		  for( i=izero; i<k; i++ ) {
		    B[j*ldb + i] += temp*A[k*lda + i];
		  }
		  if ( NoUnit )
		    temp *=A[k*lda + k];
		  B[j*ldb + k] = temp;
		}
	      }
	    }
	  } else {
	    // A is lower triangular
	    for( j=izero; j<n; j++ ) {
	      for( k=m-ione; k>-ione; k-- ) {
		if( B[j*ldb + k] != zero ) {
		  temp = alpha*B[j*ldb + k];
		  B[j*ldb + k] = temp;
		  if ( NoUnit )
		    B[j*ldb + k] *= A[k*lda + k];
		  for( i=k+ione; i<m; i++ ) {
		    B[j*ldb + i] += temp*A[k*lda + i];
		  }
		}
	      }
	    }
	  }
	} else {
	  // Compute B = alpha*A'*B
	  if( Upper ) {
	    for( j=izero; j<n; j++ ) {
	      for( i=m-ione; i>-ione; i-- ) {
		temp = B[j*ldb + i];
		if( NoUnit )
		  temp *= A[i*lda + i];
		for( k=izero; k<i; k++ ) {
		  temp += A[i*lda + k]*B[j*ldb + k];
		}
		B[j*ldb + i] = alpha*temp;
	      }
	    }
	  } else {
	    for( j=izero; j<n; j++ ) {
	      for( i=izero; i<m; i++ ) {
		temp = B[j*ldb + i];
		if( NoUnit ) 
		  temp *= A[i*lda + i];
		for( k=i+ione; k<m; k++ ) {
		  temp += A[i*lda + k]*B[j*ldb + k];
		}
		B[j*ldb + i] = alpha*temp;
	      }
	    }
	  }
	}
      } else {
	// A is on the right hand side of B.
	
	if( ETranspChar[transa] == 'N' ) {
	  // Compute B = alpha*B*A

	  if( Upper ) {
	    // A is upper triangular.
	    for( j=n-ione; j>-ione; j-- ) {
	      temp = alpha;
	      if( NoUnit )
		temp *= A[j*lda + j];
	      for( i=izero; i<m; i++ ) {
		B[j*ldb + i] *= temp;
	      }
	      for( k=izero; k<j; k++ ) {
		if( A[j*lda + k] != zero ) {
		  temp = alpha*A[j*lda + k];
		  for( i=izero; i<m; i++ ) {
		    B[j*ldb + i] += temp*B[k*ldb + i];
		  }
		}
	      }
	    }
	  } else {
	    // A is lower triangular.
	    for( j=izero; j<n; j++ ) {
	      temp = alpha;
	      if( NoUnit )
		temp *= A[j*lda + j];
	      for( i=izero; i<m; i++ ) {
		B[j*ldb + i] *= temp;
	      }
	      for( k=j+ione; k<n; k++ ) {
		if( A[j*lda + k] != zero ) {
		  temp = alpha*A[j*lda + k];
		  for( i=izero; i<m; i++ ) {
		    B[j*ldb + i] += temp*B[k*ldb + i];
		  }
		}
	      }
	    }
	  }
	} else {
	  // Compute B = alpha*B*A'

	  if( Upper ) {
	    for( k=izero; k<n; k++ ) {
	      for( j=izero; j<k; j++ ) {
		if( A[k*lda + j] != zero ) {
		  temp = alpha*A[k*lda + j];
		  for( i=izero; i<m; i++ ) {
		    B[j*ldb + i] += temp*B[k*ldb + i];
		  }
		}
	      }
	      temp = alpha;
	      if( NoUnit ) 
		temp *= A[k*lda + k];
	      if( temp != one ) {
		for( i=izero; i<m; i++ ) {
		  B[k*ldb + i] *= temp;
		}
	      }
	    }
	  } else {
	    for( k=n-ione; k>-ione; k-- ) {
	      for( j=k+ione; j<n; j++ ) {
		if( A[k*lda + j] != zero ) {
		  temp = alpha*A[k*lda + j];
		  for( i=izero; i<m; i++ ) {
		    B[j*ldb + i] += temp*B[k*ldb + i];
		  }
		}
	      }
	      temp = alpha;
	      if( NoUnit )
		temp *= A[k*lda + k];
	      if( temp != one ) {
		for( i=izero; i<m; i++ ) {
		  B[k*ldb + i] *= temp;
		}
	      }
	    }
	  }
	} // end if( ETranspChar[transa] == 'N' ) ...
      } // end if ( LSide ) ...
    } // end if (!BadArgument)
  } // end TRMM
  
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::TRSM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const OrdinalType m, const OrdinalType n, const ScalarType alpha, const ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb) const
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    ScalarType zero = ScalarTraits<ScalarType>::zero();
    ScalarType one = ScalarTraits<ScalarType>::one();
    ScalarType temp;
    OrdinalType NRowA = m;
    bool BadArgument = false;
    bool NoUnit = (EDiagChar[diag]=='N');
    
    if (!(ESideChar[side] == 'L')) { NRowA = n; }

    // Quick return.
    if (n == izero || m == izero) { return; }
    if( m < izero ) {
      cout << "BLAS::TRSM Error: M == "<<m<<endl;
      BadArgument = true; }
    if( n < izero ) {
      cout << "BLAS::TRSM Error: N == "<<n<<endl;
      BadArgument = true; }
    if( lda < NRowA ) {
      cout << "BLAS::TRSM Error: LDA == "<<lda<<endl;
      BadArgument = true; }
    if( ldb < m ) {
      cout << "BLAS::TRSM Error: LDB == "<<ldb<<endl;
      BadArgument = true; }

    if(!BadArgument)
      {
	int i, j, k;
	// Set the solution to the zero vector.
	if(alpha == zero) {
	    for(j = izero; j < n; j++) {
	    	for( i = izero; i < m; i++) {
		    B[j*ldb+i] = zero;
	      	}
	    }
	}
	else 
	{ // Start the operations.
	    if(ESideChar[side] == 'L') {
		//
	    	// Perform computations for OP(A)*X = alpha*B	    
		//
		if(ETranspChar[transa] == 'N') {
		    //
		    //  Compute B = alpha*inv( A )*B
		    //
		    if(EUploChar[uplo] == 'U') { 
			// A is upper triangular.
			for(j = izero; j < n; j++) {
	    		    // Perform alpha*B if alpha is not 1.
	    		    if(alpha != one) {
	    	    		for( i = izero; i < m; i++) {
		    		    B[j*ldb+i] *= alpha;
		    		}
			    }
			    // Perform a backsolve for column j of B.
			    for(k = (m - ione); k > -ione; k--) {
				// If this entry is zero, we don't have to do anything.
				if (B[j*ldb + k] != zero) {
				    if (NoUnit) {
					B[j*ldb + k] /= A[k*lda + k];
				    }
				    for(i = izero; i < k; i++) {
					B[j*ldb + i] -= B[j*ldb + k] * A[k*lda + i];
				    }
				}
			    }
			}
		    }
		    else 
		    { // A is lower triangular.
                        for(j = izero; j < n; j++) {
                            // Perform alpha*B if alpha is not 1.
                            if(alpha != one) {
                                for( i = izero; i < m; i++) {
                                    B[j*ldb+i] *= alpha;
                                }
                            }
                            // Perform a forward solve for column j of B.
                            for(k = izero; k < m; k++) {
                                // If this entry is zero, we don't have to do anything.
                                if (B[j*ldb + k] != zero) {   
                                    if (NoUnit) {
                                        B[j*ldb + k] /= A[k*lda + k];
                                    }
                                    for(i = k+ione; i < m; i++) {
                                        B[j*ldb + i] -= B[j*ldb + k] * A[k*lda + i];
                                    }
                                }
                            }
                        }
		    } // end if (uplo == 'U')
		}  // if (transa =='N')	
	    	else { 
		    //
		    //  Compute B = alpha*inv( A' )*B
		    //
		    if(EUploChar[uplo] == 'U') { 
			// A is upper triangular.
			for(j = izero; j < n; j++) {
	    	    	    for( i = izero; i < m; i++) {
		    		temp = alpha*B[j*ldb+i];
			    	for(k = izero; k < i; k++) {
				    temp -= A[i*lda + k] * B[j*ldb + k];
				}
				if (NoUnit) {
				    temp /= A[i*lda + i];
				}
				B[j*ldb + i] = temp;
			    }
			}
		    }
		    else
		    { // A is lower triangular.
                        for(j = izero; j < n; j++) {
                            for(i = (m - ione) ; i > -ione; i--) {
                                temp = alpha*B[j*ldb+i];
                            	for(k = i+ione; k < m; k++) {
				    temp -= A[i*lda + k] * B[j*ldb + k];
				}
				if (NoUnit) {
				    temp /= A[i*lda + i];
				}
				B[j*ldb + i] = temp;
                            }
                        }
		    }
		}
	    }  // if (side == 'L')
	    else { 
	       // side == 'R'
	       //
	       // Perform computations for X*OP(A) = alpha*B	    
	       //
	      if (ETranspChar[transa] == 'N') {
		    //
		    //  Compute B = alpha*B*inv( A )
		    //
		    if(EUploChar[uplo] == 'U') { 
			// A is upper triangular.
	    		// Perform a backsolve for column j of B.
			for(j = izero; j < n; j++) {
	    		    // Perform alpha*B if alpha is not 1.
	    		    if(alpha != one) {
	    	    		for( i = izero; i < m; i++) {
		    		    B[j*ldb+i] *= alpha;
		    		}
			    }
			    for(k = izero; k < j; k++) {
				// If this entry is zero, we don't have to do anything.
				if (A[j*lda + k] != zero) {
				    for(i = izero; i < m; i++) {
					B[j*ldb + i] -= A[j*lda + k] * B[k*ldb + i];
				    }
				}
			    }
			    if (NoUnit) {
				temp = one/A[j*lda + j];
				for(i = izero; i < m; i++) {
				    B[j*ldb + i] *= temp;
				}
			    }
			}
		    }
		    else 
		    { // A is lower triangular.
                        for(j = (n - ione); j > -ione; j--) {
                            // Perform alpha*B if alpha is not 1.
                            if(alpha != one) {
                                for( i = izero; i < m; i++) {
                                    B[j*ldb+i] *= alpha;
                                }
                            }
                            // Perform a forward solve for column j of B.
                            for(k = j+ione; k < n; k++) {
                                // If this entry is zero, we don't have to do anything.
				if (A[j*lda + k] != zero) {
				    for(i = izero; i < m; i++) {
                                        B[j*ldb + i] -= A[j*lda + k] * B[k*ldb + i]; 
                                    }
                                } 
                            }
			    if (NoUnit) {
				temp = one/A[j*lda + j];
				for(i = izero; i < m; i++) {
				    B[j*ldb + i] *= temp;
				}
			    }			
                        }
		    } // end if (uplo == 'U')
		}  // if (transa =='N')	
	    	else { 
		    //
		    //  Compute B = alpha*B*inv( A' )
		    //
		    if(EUploChar[uplo] == 'U') { 
			// A is upper triangular.
			for(k = (n - ione); k > -ione; k--) {
			    if (NoUnit) {
				temp = one/A[k*lda + k];
	    	    	    	for(i = izero; i < m; i++) {
		    		    B[k*ldb + i] *= temp;
				}
			    }
			    for(j = izero; j < k; j++) {
				if (A[k*lda + j] != zero) {
				    temp = A[k*lda + j];
				    for(i = izero; i < m; i++) {
					B[j*ldb + i] -= temp*B[k*ldb + i];
				    }
				}
			    }
			    if (alpha != one) {
				for (i = izero; i < m; i++) {
				    B[k*ldb + i] *= alpha;
				}
			    }
			}
		    }
		    else
		    { // A is lower triangular.
			for(k = izero; k < n; k++) {
			    if (NoUnit) {
				temp = one/A[k*lda + k];
				for (i = izero; i < m; i++) {
				    B[k*ldb + i] *= temp;
				}
			    }
			    for(j = k+ione; j < n; j++) {
				if(A[k*lda + j] != zero) {
				    temp = A[k*lda + j];
				    for(i = izero; i < m; i++) {
					B[j*ldb + i] -= temp*B[k*ldb + i];
				    }
				}
			    }
			    if (alpha != one) {
				for (i = izero; i < m; i++) {
				    B[k*ldb + i] *= alpha;
				}
			    }
                        }
		    }
		}		
	    }
	}
    }
  }
  
#ifndef DOXYGEN_SHOULD_SKIP_THIS

  template<typename OrdinalType>
  class BLAS<OrdinalType, float>
  {    
  public:
    inline BLAS(void) {};
    inline BLAS(const BLAS& BLAS_source) {};
    inline virtual ~BLAS(void) {};
    void ROTG(float* da, float* db, float* c, float* s) const;
    float ASUM(const OrdinalType n, const float* x, const OrdinalType incx) const;
    void AXPY(const OrdinalType n, const float alpha, const float* x, const OrdinalType incx, float* y, const OrdinalType incy) const;
    void COPY(const OrdinalType n, const float* x, const OrdinalType incx, float* y, const OrdinalType incy) const;
    float DOT(const OrdinalType n, const float* x, const OrdinalType incx, const float* y, const OrdinalType incy) const;
    float NRM2(const OrdinalType n, const float* x, const OrdinalType incx) const;
    void SCAL(const OrdinalType n, const float alpha, float* x, const OrdinalType incx) const;
    OrdinalType IAMAX(const OrdinalType n, const float* x, const OrdinalType incx) const;
    void GEMV(ETransp trans, const OrdinalType m, const OrdinalType n, const float alpha, const float* A, const OrdinalType lda, const float* x, const OrdinalType incx, const float beta, float* y, const OrdinalType incy) const;
    void TRMV(EUplo uplo, ETransp trans, EDiag diag, const OrdinalType n, const float* A, const OrdinalType lda, float* x, const OrdinalType incx) const;
    void GER(const OrdinalType m, const OrdinalType n, const float alpha, const float* x, const OrdinalType incx, const float* y, const OrdinalType incy, float* A, const OrdinalType lda) const;
    void GEMM(ETransp transa, ETransp transb, const OrdinalType m, const OrdinalType n, const OrdinalType k, const float alpha, const float* A, const OrdinalType lda, const float* B, const OrdinalType ldb, const float beta, float* C, const OrdinalType ldc) const;
    void SYMM(ESide side, EUplo uplo, const OrdinalType m, const OrdinalType n, const float alpha, const float* A, const OrdinalType lda, const float *B, const OrdinalType ldb, const float beta, float *C, const OrdinalType ldc) const;
    void TRMM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const OrdinalType m, const OrdinalType n, const float alpha, const float* A, const OrdinalType lda, float* B, const OrdinalType ldb) const;
    void TRSM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const OrdinalType m, const OrdinalType n, const float alpha, const float* A, const OrdinalType lda, float* B, const OrdinalType ldb) const;
  };

  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::ROTG(float* da, float* db, float* c, float* s) const
  { SROTG_F77(da, db, c, s ); }

  template<typename OrdinalType>
  float BLAS<OrdinalType, float>::ASUM(const OrdinalType n, const float* x, const OrdinalType incx) const
  { return SASUM_F77(&n, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::AXPY(const OrdinalType n, const float alpha, const float* x, const OrdinalType incx, float* y, const OrdinalType incy) const
  { SAXPY_F77(&n, &alpha, x, &incx, y, &incy); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::COPY(const OrdinalType n, const float* x, const OrdinalType incx, float* y, const OrdinalType incy) const 
  { SCOPY_F77(&n, x, &incx, y, &incy); }
  
  template<typename OrdinalType>
  float BLAS<OrdinalType, float>::DOT(const OrdinalType n, const float* x, const OrdinalType incx, const float* y, const OrdinalType incy) const
  { return SDOT_F77(&n, x, &incx, y, &incy); }
  
  template<typename OrdinalType>
  OrdinalType BLAS<OrdinalType, float>::IAMAX(const OrdinalType n, const float* x, const OrdinalType incx) const
  { return ISAMAX_F77(&n, x, &incx); }

  template<typename OrdinalType>
  float BLAS<OrdinalType, float>::NRM2(const OrdinalType n, const float* x, const OrdinalType incx) const
  { return SNRM2_F77(&n, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::SCAL(const OrdinalType n, const float alpha, float* x, const OrdinalType incx) const
  { SSCAL_F77(&n, &alpha, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::GEMV(ETransp trans, const OrdinalType m, const OrdinalType n, const float alpha, const float* A, const OrdinalType lda, const float* x, const OrdinalType incx, const float beta, float* y, const OrdinalType incy) const
  { SGEMV_F77(&ETranspChar[trans], &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::GER(const OrdinalType m, const OrdinalType n, const float alpha, const float* x, const OrdinalType incx, const float* y, const OrdinalType incy, float* A, const OrdinalType lda) const
  { SGER_F77(&m, &n, &alpha, x, &incx, y, &incy, A, &lda); }

  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::TRMV(EUplo uplo, ETransp trans, EDiag diag, const OrdinalType n, const float* A, const OrdinalType lda, float* x, const OrdinalType incx) const
  { STRMV_F77(&EUploChar[uplo], &ETranspChar[trans], &EDiagChar[diag], &n, A, &lda, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::GEMM(ETransp transa, ETransp transb, const OrdinalType m, const OrdinalType n, const OrdinalType k, const float alpha, const float* A, const OrdinalType lda, const float* B, const OrdinalType ldb, const float beta, float* C, const OrdinalType ldc) const
  { SGEMM_F77(&ETranspChar[transa], &ETranspChar[transb], &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::SYMM(ESide side, EUplo uplo, const OrdinalType m, const OrdinalType n, const float alpha, const float* A, const OrdinalType lda, const float* B, const OrdinalType ldb, const float beta, float* C, const OrdinalType ldc) const
  { SSYMM_F77(&ESideChar[side], &EUploChar[uplo], &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::TRMM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const OrdinalType m, const OrdinalType n, const float alpha, const float* A, const OrdinalType lda, float* B, const OrdinalType ldb) const
  { STRMM_F77(&ESideChar[side], &EUploChar[uplo], &ETranspChar[transa], &EDiagChar[diag], &m, &n, &alpha, A, &lda, B, &ldb); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::TRSM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const OrdinalType m, const OrdinalType n, const float alpha, const float* A, const OrdinalType lda, float* B, const OrdinalType ldb) const
  { STRSM_F77(&ESideChar[side], &EUploChar[uplo], &ETranspChar[transa], &EDiagChar[diag], &m, &n, &alpha, A, &lda, B, &ldb); }


  template<typename OrdinalType>
  class BLAS<OrdinalType, double>
  {    
  public:
    inline BLAS(void) {};
    inline BLAS(const BLAS& BLAS_source) {};
    inline virtual ~BLAS(void) {};
    void ROTG(double* da, double* db, double* c, double* s) const;
    double ASUM(const OrdinalType n, const double* x, const OrdinalType incx) const;
    void AXPY(const OrdinalType n, const double alpha, const double* x, const OrdinalType incx, double* y, const OrdinalType incy) const;
    void COPY(const OrdinalType n, const double* x, const OrdinalType incx, double* y, const OrdinalType incy) const;
    double DOT(const OrdinalType n, const double* x, const OrdinalType incx, const double* y, const OrdinalType incy) const;
    double NRM2(const OrdinalType n, const double* x, const OrdinalType incx) const;
    void SCAL(const OrdinalType n, const double alpha, double* x, const OrdinalType incx) const;
    OrdinalType IAMAX(const OrdinalType n, const double* x, const OrdinalType incx) const;
    void GEMV(ETransp trans, const OrdinalType m, const OrdinalType n, const double alpha, const double* A, const OrdinalType lda, const double* x, const OrdinalType incx, const double beta, double* y, const OrdinalType incy) const;
    void TRMV(EUplo uplo, ETransp trans, EDiag diag, const OrdinalType n, const double* A, const OrdinalType lda, double* x, const OrdinalType incx) const;
    void GER(const OrdinalType m, const OrdinalType n, const double alpha, const double* x, const OrdinalType incx, const double* y, const OrdinalType incy, double* A, const OrdinalType lda) const;
    void GEMM(ETransp transa, ETransp transb, const OrdinalType m, const OrdinalType n, const OrdinalType k, const double alpha, const double* A, const OrdinalType lda, const double* B, const OrdinalType ldb, const double beta, double* C, const OrdinalType ldc) const;
    void SYMM(ESide side, EUplo uplo, const OrdinalType m, const OrdinalType n, const double alpha, const double* A, const OrdinalType lda, const double *B, const OrdinalType ldb, const double beta, double *C, const OrdinalType ldc) const;
    void TRMM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const OrdinalType m, const OrdinalType n, const double alpha, const double* A, const OrdinalType lda, double* B, const OrdinalType ldb) const;
    void TRSM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const OrdinalType m, const OrdinalType n, const double alpha, const double* A, const OrdinalType lda, double* B, const OrdinalType ldb) const;
  };
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::ROTG(double* da, double* db, double* c, double* s) const
  { DROTG_F77(da, db, c, s); }

  template<typename OrdinalType>
  double BLAS<OrdinalType, double>::ASUM(const OrdinalType n, const double* x, const OrdinalType incx) const
  { return DASUM_F77(&n, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::AXPY(const OrdinalType n, const double alpha, const double* x, const OrdinalType incx, double* y, const OrdinalType incy) const
  { DAXPY_F77(&n, &alpha, x, &incx, y, &incy); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::COPY(const OrdinalType n, const double* x, const OrdinalType incx, double* y, const OrdinalType incy) const
  { DCOPY_F77(&n, x, &incx, y, &incy); }
  
  template<typename OrdinalType>
  double BLAS<OrdinalType, double>::DOT(const OrdinalType n, const double* x, const OrdinalType incx, const double* y, const OrdinalType incy) const
  { return DDOT_F77(&n, x, &incx, y, &incy); }
  
  template<typename OrdinalType>
  OrdinalType BLAS<OrdinalType, double>::IAMAX(const OrdinalType n, const double* x, const OrdinalType incx) const
  { return IDAMAX_F77(&n, x, &incx); }

  template<typename OrdinalType>
  double BLAS<OrdinalType, double>::NRM2(const OrdinalType n, const double* x, const OrdinalType incx) const
  { return DNRM2_F77(&n, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::SCAL(const OrdinalType n, const double alpha, double* x, const OrdinalType incx) const
  { DSCAL_F77(&n, &alpha, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::GEMV(ETransp trans, const OrdinalType m, const OrdinalType n, const double alpha, const double* A, const OrdinalType lda, const double* x, const OrdinalType incx, const double beta, double* y, const OrdinalType incy) const
  { DGEMV_F77(&ETranspChar[trans], &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::GER(const OrdinalType m, const OrdinalType n, const double alpha, const double* x, const OrdinalType incx, const double* y, const OrdinalType incy, double* A, const OrdinalType lda) const
  { DGER_F77(&m, &n, &alpha, x, &incx, y, &incy, A, &lda); }

  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::TRMV(EUplo uplo, ETransp trans, EDiag diag, const OrdinalType n, const double* A, const OrdinalType lda, double* x, const OrdinalType incx) const
  { DTRMV_F77(&EUploChar[uplo], &ETranspChar[trans], &EDiagChar[diag], &n, A, &lda, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::GEMM(ETransp transa, ETransp transb, const OrdinalType m, const OrdinalType n, const OrdinalType k, const double alpha, const double* A, const OrdinalType lda, const double* B, const OrdinalType ldb, const double beta, double* C, const OrdinalType ldc) const
  { DGEMM_F77(&ETranspChar[transa], &ETranspChar[transb], &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::SYMM(ESide side, EUplo uplo, const OrdinalType m, const OrdinalType n, const double alpha, const double* A, const OrdinalType lda, const double *B, const OrdinalType ldb, const double beta, double *C, const OrdinalType ldc) const
  { DSYMM_F77(&ESideChar[side], &EUploChar[uplo], &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::TRMM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const OrdinalType m, const OrdinalType n, const double alpha, const double* A, const OrdinalType lda, double* B, const OrdinalType ldb) const
  { DTRMM_F77(&ESideChar[side], &EUploChar[uplo], &ETranspChar[transa], &EDiagChar[diag], &m, &n, &alpha, A, &lda, B, &ldb); }

  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::TRSM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const OrdinalType m, const OrdinalType n, const double alpha, const double* A, const OrdinalType lda, double* B, const OrdinalType ldb) const
  { DTRSM_F77(&ESideChar[side], &EUploChar[uplo], &ETranspChar[transa], &EDiagChar[diag], &m, &n, &alpha, A, &lda, B, &ldb); }
  
#ifdef HAVE_TEUCHOS_COMPLEX

  template<typename OrdinalType>
  class BLAS<OrdinalType, complex<float> >
  {    
  public:
    inline BLAS(void) {};
    inline BLAS(const BLAS& BLAS_source) {};
    inline virtual ~BLAS(void) {};
    void ROTG(complex<float>* da, complex<float>* db, complex<float>* c, complex<float>* s) const;
    complex<float> ASUM(const OrdinalType n, const complex<float>* x, const OrdinalType incx) const;
    void AXPY(const OrdinalType n, const complex<float> alpha, const complex<float>* x, const OrdinalType incx, complex<float>* y, const OrdinalType incy) const;
    void COPY(const OrdinalType n, const complex<float>* x, const OrdinalType incx, complex<float>* y, const OrdinalType incy) const;
    complex<float> DOT(const OrdinalType n, const complex<float>* x, const OrdinalType incx, const complex<float>* y, const OrdinalType incy) const;
    complex<float> NRM2(const OrdinalType n, const complex<float>* x, const OrdinalType incx) const;
    void SCAL(const OrdinalType n, const complex<float> alpha, complex<float>* x, const OrdinalType incx) const;
    OrdinalType IAMAX(const OrdinalType n, const complex<float>* x, const OrdinalType incx) const;
    void GEMV(ETransp trans, const OrdinalType m, const OrdinalType n, const complex<float> alpha, const complex<float>* A, const OrdinalType lda, const complex<float>* x, const OrdinalType incx, const complex<float> beta, complex<float>* y, const OrdinalType incy) const;
    void TRMV(EUplo uplo, ETransp trans, EDiag diag, const OrdinalType n, const complex<float>* A, const OrdinalType lda, complex<float>* x, const OrdinalType incx) const;
    void GER(const OrdinalType m, const OrdinalType n, const complex<float> alpha, const complex<float>* x, const OrdinalType incx, const complex<float>* y, const OrdinalType incy, complex<float>* A, const OrdinalType lda) const;
    void GEMM(ETransp transa, ETransp transb, const OrdinalType m, const OrdinalType n, const OrdinalType k, const complex<float> alpha, const complex<float>* A, const OrdinalType lda, const complex<float>* B, const OrdinalType ldb, const complex<float> beta, complex<float>* C, const OrdinalType ldc) const;
    void SYMM(ESide side, EUplo uplo, const OrdinalType m, const OrdinalType n, const complex<float> alpha, const complex<float>* A, const OrdinalType lda, const complex<float> *B, const OrdinalType ldb, const complex<float> beta, complex<float> *C, const OrdinalType ldc) const;
    void TRMM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const OrdinalType m, const OrdinalType n, const complex<float> alpha, const complex<float>* A, const OrdinalType lda, complex<float>* B, const OrdinalType ldb) const;
    void TRSM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const OrdinalType m, const OrdinalType n, const complex<float> alpha, const complex<float>* A, const OrdinalType lda, complex<float>* B, const OrdinalType ldb) const;
  };

  template<typename OrdinalType>
  void BLAS<OrdinalType, complex<float> >::ROTG(complex<float>* da, complex<float>* db, complex<float>* c, complex<float>* s) const
  { SROTG_F77(da, db, c, s ); }

  template<typename OrdinalType>
  complex<float> BLAS<OrdinalType, complex<float> >::ASUM(const OrdinalType n, const complex<float>* x, const OrdinalType incx) const
  { return SASUM_F77(&n, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, complex<float> >::AXPY(const OrdinalType n, const complex<float> alpha, const complex<float>* x, const OrdinalType incx, complex<float>* y, const OrdinalType incy) const
  { SAXPY_F77(&n, &alpha, x, &incx, y, &incy); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, complex<float> >::COPY(const OrdinalType n, const complex<float>* x, const OrdinalType incx, complex<float>* y, const OrdinalType incy) const
  { SCOPY_F77(&n, x, &incx, y, &incy); }
  
  template<typename OrdinalType>
  complex<float> BLAS<OrdinalType, complex<float> >::DOT(const OrdinalType n, const complex<float>* x, const OrdinalType incx, const complex<float>* y, const OrdinalType incy) const
  { return SDOT_F77(&n, x, &incx, y, &incy); }
  
  template<typename OrdinalType>
  OrdinalType BLAS<OrdinalType, complex<float> >::IAMAX(const OrdinalType n, const complex<float>* x, const OrdinalType incx) const
  { return ISAMAX_F77(&n, x, &incx); }

  template<typename OrdinalType>
  complex<float> BLAS<OrdinalType, complex<float> >::NRM2(const OrdinalType n, const complex<float>* x, const OrdinalType incx) const
  { return SNRM2_F77(&n, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, complex<float> >::SCAL(const OrdinalType n, const complex<float> alpha, complex<float>* x, const OrdinalType incx) const
  { SSCAL_F77(&n, &alpha, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, complex<float> >::GEMV(ETransp trans, const OrdinalType m, const OrdinalType n, const complex<float> alpha, const complex<float>* A, const OrdinalType lda, const complex<float>* x, const OrdinalType incx, const complex<float> beta, complex<float>* y, const OrdinalType incy) const
  { SGEMV_F77(&ETranspChar[trans], &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, complex<float> >::GER(const OrdinalType m, const OrdinalType n, const complex<float> alpha, const complex<float>* x, const OrdinalType incx, const complex<float>* y, const OrdinalType incy, complex<float>* A, const OrdinalType lda) const
  { SGER_F77(&m, &n, &alpha, x, &incx, y, &incy, A, &lda); }

  template<typename OrdinalType>
  void BLAS<OrdinalType, complex<float> >::TRMV(EUplo uplo, ETransp trans, EDiag diag, const OrdinalType n, const complex<float>* A, const OrdinalType lda, complex<float>* x, const OrdinalType incx) const
  { STRMV_F77(&EUploChar[uplo], &ETranspChar[trans], &EDiagChar[diag], &n, A, &lda, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, complex<float> >::GEMM(ETransp transa, ETransp transb, const OrdinalType m, const OrdinalType n, const OrdinalType k, const complex<float> alpha, const complex<float>* A, const OrdinalType lda, const complex<float>* B, const OrdinalType ldb, const complex<float> beta, complex<float>* C, const OrdinalType ldc) const
  { SGEMM_F77(&ETranspChar[transa], &ETranspChar[transb], &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); } 
 
  template<typename OrdinalType>
  void BLAS<OrdinalType, complex<float> >::SYMM(ESide side, EUplo uplo, const OrdinalType m, const OrdinalType n, const complex<float> alpha, const complex<float>* A, const OrdinalType lda, const complex<float>* B, const OrdinalType ldb, const complex<float> beta, complex<float>* C, const OrdinalType ldc) const
  { SSYMM_F77(&ESideChar[side], &EUploChar[uplo], &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, complex<float> >::TRMM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const OrdinalType m, const OrdinalType n, const complex<float> alpha, const complex<float>* A, const OrdinalType lda, complex<float>* B, const OrdinalType ldb) const
  { STRMM_F77(&ESideChar[side], &EUploChar[uplo], &ETranspChar[transa], &EDiagChar[diag], &m, &n, &alpha, A, &lda, B, &ldb); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, complex<float> >::TRSM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const OrdinalType m, const OrdinalType n, const complex<float> alpha, const complex<float>* A, const OrdinalType lda, complex<float>* B, const OrdinalType ldb) const
  { STRSM_F77(&ESideChar[side], &EUploChar[uplo], &ETranspChar[transa], &EDiagChar[diag], &m, &n, &alpha, A, &lda, B, &ldb); }


  template<typename OrdinalType>
  class BLAS<OrdinalType, complex<double> >
  {    
  public:
    inline BLAS(void) {};
    inline BLAS(const BLAS& BLAS_source) {};
    inline virtual ~BLAS(void) {};
    void ROTG(complex<double>* da, complex<double>* db, complex<double>* c, complex<double>* s) const;
    complex<double> ASUM(const OrdinalType n, const complex<double>* x, const OrdinalType incx) const;
    void AXPY(const OrdinalType n, const complex<double> alpha, const complex<double>* x, const OrdinalType incx, complex<double>* y, const OrdinalType incy) const;
    void COPY(const OrdinalType n, const complex<double>* x, const OrdinalType incx, complex<double>* y, const OrdinalType incy) const;
    complex<double> DOT(const OrdinalType n, const complex<double>* x, const OrdinalType incx, const complex<double>* y, const OrdinalType incy) const;
    complex<double> NRM2(const OrdinalType n, const complex<double>* x, const OrdinalType incx) const;
    void SCAL(const OrdinalType n, const complex<double> alpha, complex<double>* x, const OrdinalType incx) const;
    OrdinalType IAMAX(const OrdinalType n, const complex<double>* x, const OrdinalType incx) const;
    void GEMV(ETransp trans, const OrdinalType m, const OrdinalType n, const complex<double> alpha, const complex<double>* A, const OrdinalType lda, const complex<double>* x, const OrdinalType incx, const complex<double> beta, complex<double>* y, const OrdinalType incy) const;
    void TRMV(EUplo uplo, ETransp trans, EDiag diag, const OrdinalType n, const complex<double>* A, const OrdinalType lda, complex<double>* x, const OrdinalType incx) const;
    void GER(const OrdinalType m, const OrdinalType n, const complex<double> alpha, const complex<double>* x, const OrdinalType incx, const complex<double>* y, const OrdinalType incy, complex<double>* A, const OrdinalType lda) const;
    void GEMM(ETransp transa, ETransp transb, const OrdinalType m, const OrdinalType n, const OrdinalType k, const complex<double> alpha, const complex<double>* A, const OrdinalType lda, const complex<double>* B, const OrdinalType ldb, const complex<double> beta, complex<double>* C, const OrdinalType ldc) const;
    void SYMM(ESide side, EUplo uplo, const OrdinalType m, const OrdinalType n, const complex<double> alpha, const complex<double>* A, const OrdinalType lda, const complex<double> *B, const OrdinalType ldb, const complex<double> beta, complex<double> *C, const OrdinalType ldc) const;
    void TRMM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const OrdinalType m, const OrdinalType n, const complex<double> alpha, const complex<double>* A, const OrdinalType lda, complex<double>* B, const OrdinalType ldb) const;
    void TRSM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const OrdinalType m, const OrdinalType n, const complex<double> alpha, const complex<double>* A, const OrdinalType lda, complex<double>* B, const OrdinalType ldb) const;
  };
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, complex<double> >::ROTG(complex<double>* da, complex<double>* db, complex<double>* c, complex<double>* s) const
  { DROTG_F77(da, db, c, s); }

  template<typename OrdinalType>
  complex<double> BLAS<OrdinalType, complex<double> >::ASUM(const OrdinalType n, const complex<double>* x, const OrdinalType incx) const
  { return DASUM_F77(&n, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, complex<double> >::AXPY(const OrdinalType n, const complex<double> alpha, const complex<double>* x, const OrdinalType incx, complex<double>* y, const OrdinalType incy) const
  { DAXPY_F77(&n, &alpha, x, &incx, y, &incy); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, complex<double> >::COPY(const OrdinalType n, const complex<double>* x, const OrdinalType incx, complex<double>* y, const OrdinalType incy) const
  { DCOPY_F77(&n, x, &incx, y, &incy); }
  
  template<typename OrdinalType>
  complex<double> BLAS<OrdinalType, complex<double> >::DOT(const OrdinalType n, const complex<double>* x, const OrdinalType incx, const complex<double>* y, const OrdinalType incy) const
  { return DDOT_F77(&n, x, &incx, y, &incy); }
  
  template<typename OrdinalType>
  OrdinalType BLAS<OrdinalType, complex<double> >::IAMAX(const OrdinalType n, const complex<double>* x, const OrdinalType incx) const
  { return IDAMAX_F77(&n, x, &incx); }

  template<typename OrdinalType>
  complex<double> BLAS<OrdinalType, complex<double> >::NRM2(const OrdinalType n, const complex<double>* x, const OrdinalType incx) const
  { return DNRM2_F77(&n, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, complex<double> >::SCAL(const OrdinalType n, const complex<double> alpha, complex<double>* x, const OrdinalType incx) const
  { DSCAL_F77(&n, &alpha, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, complex<double> >::GEMV(ETransp trans, const OrdinalType m, const OrdinalType n, const complex<double> alpha, const complex<double>* A, const OrdinalType lda, const complex<double>* x, const OrdinalType incx, const complex<double> beta, complex<double>* y, const OrdinalType incy) const
  { DGEMV_F77(&ETranspChar[trans], &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, complex<double> >::GER(const OrdinalType m, const OrdinalType n, const complex<double> alpha, const complex<double>* x, const OrdinalType incx, const complex<double>* y, const OrdinalType incy, complex<double>* A, const OrdinalType lda) const
  { DGER_F77(&m, &n, &alpha, x, &incx, y, &incy, A, &lda); }

  template<typename OrdinalType>
  void BLAS<OrdinalType, complex<double> >::TRMV(EUplo uplo, ETransp trans, EDiag diag, const OrdinalType n, const complex<double>* A, const OrdinalType lda, complex<double>* x, const OrdinalType incx) const
  { DTRMV_F77(&EUploChar[uplo], &ETranspChar[trans], &EDiagChar[diag], &n, A, &lda, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, complex<double> >::GEMM(ETransp transa, ETransp transb, const OrdinalType m, const OrdinalType n, const OrdinalType k, const complex<double> alpha, const complex<double>* A, const OrdinalType lda, const complex<double>* B, const OrdinalType ldb, const complex<double> beta, complex<double>* C, const OrdinalType ldc) const
  { DGEMM_F77(&ETranspChar[transa], &ETranspChar[transb], &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, complex<double> >::SYMM(ESide side, EUplo uplo, const OrdinalType m, const OrdinalType n, const complex<double> alpha, const complex<double>* A, const OrdinalType lda, const complex<double> *B, const OrdinalType ldb, const complex<double> beta, complex<double> *C, const OrdinalType ldc) const
  { DSYMM_F77(&ESideChar[side], &EUploChar[uplo], &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, complex<double> >::TRMM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const OrdinalType m, const OrdinalType n, const complex<double> alpha, const complex<double>* A, const OrdinalType lda, complex<double>* B, const OrdinalType ldb) const
  { DTRMM_F77(&ESideChar[side], &EUploChar[uplo], &ETranspChar[transa], &EDiagChar[diag], &m, &n, &alpha, A, &lda, B, &ldb); }

  template<typename OrdinalType>
  void BLAS<OrdinalType, complex<double> >::TRSM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const OrdinalType m, const OrdinalType n, const complex<double> alpha, const complex<double>* A, const OrdinalType lda, complex<double>* B, const OrdinalType ldb) const
  { DTRSM_F77(&ESideChar[side], &EUploChar[uplo], &ETranspChar[transa], &EDiagChar[diag], &m, &n, &alpha, A, &lda, B, &ldb); }
  
#endif // HAVE_TEUCHOS_COMPLEX

#endif // DOXYGEN_SHOULD_SKIP_THIS

} // namespace Teuchos

#endif // _TEUCHOS_BLAS_HPP_
