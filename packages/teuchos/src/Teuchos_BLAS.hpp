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

#ifndef _TEUCHOS_BLAS_HPP_
#define _TEUCHOS_BLAS_HPP_

#include "Teuchos_BLAS_wrappers.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_OrdinalTraits.hpp"

namespace Teuchos
{
//! Teuchos::BLAS: The Templated BLAS Class.
/*! The Teuchos::BLAS class provides functionality similar to the BLAS
    (Basic Linear Algebra Subprograms).  The BLAS provide portable, high-
    performance implementations of kernels such as dense vectoer multiplication,
    dot products, dense matrix-vector multiplication and dense matrix-matrix
    multiplication.

    The standard BLAS interface is Fortran-specific.  Unfortunately, the 
    interface between C++ and Fortran is not standard across all computer
    platforms.  The Teuchos_BLAS class provides C++ bindings for the BLAS
    kernels in order to insulate the rest of Petra from the details of 
    C++ to Fortran translation.

    In addition to giving access the standard BLAS functionality.
    Teuchos::BLAS also provide functionality for any <scalarType> class that
    defines the +, - * and / operators.

    Teuchos::BLAS is a single memory image interface only.  This is appropriate 
    since the standard BLAS are only specified for serial execution 
    (or shared memory parallel).
*/

  template<typename OrdinalType, typename ScalarType>
  class BLAS
  {    
  public:
    //@{ \name Constructor/Destructor.

    //! Teuchos::BLAS empty constructor.
    inline BLAS(void) {};

    //! Teuchos::BLAS copy constructor.
    inline BLAS(const BLAS& BLAS_source) {};

    //! Teuchos::BLAS destructor.
    inline virtual ~BLAS(void) {};
    //@}

    //@{ \name Level 1 BLAS Routines.
    void ROTG(ScalarType da, ScalarType db, ScalarType c, ScalarType s);
    void SCAL(OrdinalType n, ScalarType alpha, ScalarType* x, OrdinalType incx);
    void COPY(OrdinalType n, ScalarType* x, OrdinalType incx, ScalarType* y, OrdinalType incy);
    void AXPY(OrdinalType n, ScalarType alpha, ScalarType* x, OrdinalType incx, ScalarType* y, OrdinalType incy);
    ScalarType ASUM(OrdinalType n, ScalarType* x, OrdinalType incx);
    ScalarType DOT(OrdinalType n, ScalarType* x, OrdinalType incx, ScalarType* y, OrdinalType incy);
    ScalarType NRM2(OrdinalType n, ScalarType* x, OrdinalType incx);
    OrdinalType IAMAX(OrdinalType n, ScalarType* x, OrdinalType incx);
    //@}

    //@{ \name Level 2 BLAS Routines.
    void GEMV(char Trans, OrdinalType m, OrdinalType n, ScalarType alpha, ScalarType* A, 
	      OrdinalType lda, ScalarType* x, OrdinalType incx, ScalarType beta, ScalarType* y, OrdinalType incy);
    void TRMV(char Uplo, char Trans, char Diag, OrdinalType n, ScalarType* A, 
	      OrdinalType lda, ScalarType* x, OrdinalType incx);
    void GER(OrdinalType m, OrdinalType n, ScalarType alpha, ScalarType* x, OrdinalType incx, 
	     ScalarType* y, OrdinalType incy, ScalarType* A, OrdinalType lda);
    //@}
    
    //@{ \name Level 3 BLAS Routines. 
    void GEMM(char TransA, char TransB, OrdinalType m, OrdinalType n, OrdinalType k, ScalarType alpha, ScalarType* A, OrdinalType lda, ScalarType* B, OrdinalType ldb, ScalarType beta, ScalarType* C, OrdinalType ldc);
    void SYMM(char Side, char Uplo, OrdinalType m, OrdinalType n, ScalarType alpha, ScalarType* A, OrdinalType lda, ScalarType* B, OrdinalType ldb, ScalarType beta, ScalarType* C, OrdinalType ldc);
    void TRMM(char Side, char Uplo, char TransA, char Diag, OrdinalType m, OrdinalType n,
                ScalarType alpha, ScalarType* A, OrdinalType lda, ScalarType* B, OrdinalType ldb);
    void TRSM(char Side, char Uplo, char TransA, char Diag, OrdinalType m, OrdinalType n,
                ScalarType alpha, ScalarType* A, OrdinalType lda, ScalarType* B, OrdinalType ldb);
    //@}
  };

//------------------------------------------------------------------------------------------
//      LEVEL 1 BLAS ROUTINES  
//------------------------------------------------------------------------------------------
    
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::ROTG(ScalarType da, ScalarType db, ScalarType c, ScalarType s)
  {
    cout<<"Not implemented right now!"<<endl;
  }

  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::SCAL(OrdinalType n, ScalarType alpha, ScalarType* x, OrdinalType incx)
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
  void BLAS<OrdinalType, ScalarType>::COPY(OrdinalType n, ScalarType* x, OrdinalType incx, ScalarType* y, OrdinalType incy)
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
  void BLAS<OrdinalType, ScalarType>::AXPY(OrdinalType n, ScalarType alpha, ScalarType* x, OrdinalType incx, ScalarType* y, OrdinalType incy)
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
  ScalarType BLAS<OrdinalType, ScalarType>::ASUM(OrdinalType n, ScalarType* x, OrdinalType incx)
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
  ScalarType BLAS<OrdinalType, ScalarType>::DOT(OrdinalType n, ScalarType* x, OrdinalType incx, ScalarType* y, OrdinalType incy)
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
  ScalarType BLAS<OrdinalType, ScalarType>::NRM2(OrdinalType n, ScalarType* x, OrdinalType incx)
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
  OrdinalType BLAS<OrdinalType, ScalarType>::IAMAX(OrdinalType n, ScalarType* x, OrdinalType incx)
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
  void BLAS<OrdinalType, ScalarType>::GEMV(char Trans, OrdinalType m, OrdinalType n, ScalarType alpha, ScalarType* A, OrdinalType lda, ScalarType* x, OrdinalType incx, ScalarType beta, ScalarType* y, OrdinalType incy)
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    ScalarType zero = ScalarTraits<ScalarType>::zero();
    ScalarType one = ScalarTraits<ScalarType>::one();
    bool BadArgument = false;

    // Quick return if there is nothing to do!
    if( m == izero || n == izero || ( alpha == zero && beta == one ) ){ return; }
    
    // Otherwise, we need to check the argument list.
    if(!((toupper(Trans) == 'N') || (toupper(Trans) == 'T') || (toupper(Trans) == 'C'))) {
        cout << "BLAS::GEMV Error: TRANS == " << Trans << endl;
      	BadArgument = true; 
    }
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
      if(toupper(Trans) == 'N') {
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

      if( Trans == 'N' ) {
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
 void BLAS<OrdinalType, ScalarType>::TRMV(char Uplo, char Trans, char Diag, OrdinalType n, ScalarType* A, OrdinalType lda, ScalarType* x, OrdinalType incx)
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    ScalarType zero = ScalarTraits<ScalarType>::zero();
    ScalarType one = ScalarTraits<ScalarType>::one();
    bool BadArgument = false;

    // Quick return if there is nothing to do!
    if( n == izero ){ return; }
    
    // Otherwise, we need to check the argument list.
    if(!((toupper(Uplo) == 'U') || (toupper(Uplo) == 'L'))) {
      cout << "BLAS::TRMV Error: UPLO == "<< Uplo << endl;
      BadArgument = true;
    }
    if(!((toupper(Trans) == 'N') || (toupper(Trans) == 'T') || (toupper(Trans) == 'C'))) {
      cout << "BLAS::TRMV Error: TRANS == " << Trans << endl;
      BadArgument = true; 
    }
    if(!((toupper(Diag) == 'U') || (toupper(Diag) == 'N'))) {
      cout << "BLAS::TRMV Error: DIAG == "<< Diag << endl;
      BadArgument = true;
    }
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
      bool NoUnit = (toupper(Diag) == 'N');

      // Set the starting pointer for the vector x if incx < 0.
      if (incx < izero) { kx = (-n+ione)*incx; }

      // Start the operations for a nontransposed triangular matrix 
      if (toupper(Trans) == 'N') {
	/* Compute x = A*x */
	if (toupper(Uplo) == 'U') {
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
	} /* if (toupper(Uplo)=='U') */
      } else { /* A is transposed/conjugated */
	/* Compute x = A'*x */
	if (toupper(Uplo)=='U') {
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
	} /* if (toupper(Uplo)=='U') */
      } /* if (toupper(Trans)=='N') */
    } /* if (!BadArgument) */
  } /* end TRMV */
        
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::GER(OrdinalType m, OrdinalType n, ScalarType alpha, ScalarType* x, OrdinalType incx, ScalarType* y, OrdinalType incy, ScalarType* A, OrdinalType lda)
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
  void BLAS<OrdinalType, ScalarType>::GEMM(char TransA, char TransB, OrdinalType m, OrdinalType n, OrdinalType k, ScalarType alpha, ScalarType* A, OrdinalType lda, ScalarType* B, OrdinalType ldb, ScalarType beta, ScalarType* C, OrdinalType ldc)
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    ScalarType zero = ScalarTraits<ScalarType>::zero();
    ScalarType one = ScalarTraits<ScalarType>::one();
    OrdinalType i, j, p;
    OrdinalType NRowA = m, NColA = k, NRowB = k;
    ScalarType temp;
    bool BadArgument = false;

    // Change dimensions of matrix if either matrix is transposed
    if( !(toupper(TransA)=='N') ) {
      NRowA = k;
      NColA = m;
    }
    if( !(toupper(TransB)=='N') ) {
      NRowB = n;
    }

    // Quick return if there is nothing to do!
    if( (m==izero) || (n==izero) || (((alpha==zero)||(k==izero)) && (beta==one)) ){ return; }
    
    // Otherwise, we need to check the argument list.
    if(!((toupper(TransA) == 'N') || (toupper(TransA) == 'T') || (toupper(TransA) == 'C'))) {
      cout << "BLAS::GEMM Error: TRANSA == " << TransA << endl;
      BadArgument = true; 
    }
    if(!((toupper(TransB) == 'N') || (toupper(TransB) == 'T') || (toupper(TransB) == 'C'))) {
      cout << "BLAS::GEMM Error: TRANSB == " << TransB << endl;
      BadArgument = true; 
    }
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
      if ( toupper(TransB)=='N' ) {
	if ( toupper(TransA)=='N' ) {
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
	if ( toupper(TransA)=='N' ) {
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
	} // end if (toupper(TransA)=='N') ...
      } // end if (toupper(TransB)=='N') ...
    } // end if (!BadArgument) ...
  } // end of GEMM


  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::SYMM(char Side, char Uplo, OrdinalType m, OrdinalType n, ScalarType alpha, ScalarType* A, OrdinalType lda, ScalarType* B, OrdinalType ldb, ScalarType beta, ScalarType* C, OrdinalType ldc)
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    ScalarType zero = ScalarTraits<ScalarType>::zero();
    ScalarType one = ScalarTraits<ScalarType>::one();
    OrdinalType i, j, k, NRowA = m;
    ScalarType temp1, temp2;
    bool BadArgument = false;
    bool Upper = (toupper(Uplo) == 'U');
    if (toupper(Side) == 'R') { NRowA = n; }
    
    // Quick return.
    if ( (m==izero) || (n==izero) || ( (alpha==zero)&&(beta==one) ) ) { return; }

    if(!((toupper(Side) == 'L') || (toupper(Side) == 'R'))) {
      cout << "BLAS::SYMM Error: SIDE == " << Side << endl;
      BadArgument = true; }
    if(!((toupper(Uplo) == 'L') || Upper)) {
      cout << "BLAS::SYMM Error: UPLO == " << Uplo << endl;
      BadArgument = true; }
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

      if ( toupper(Side) == 'L') {
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
      } // end if (toupper(Side)=='L') ...
    } // end if(!BadArgument) ...
} // end SYMM
  
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::TRMM(char Side, char Uplo, char TransA, char Diag, OrdinalType m, OrdinalType n, ScalarType alpha, ScalarType* A, OrdinalType lda, ScalarType* B, OrdinalType ldb)
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    ScalarType zero = ScalarTraits<ScalarType>::zero();
    ScalarType one = ScalarTraits<ScalarType>::one();
    OrdinalType i, j, k, NRowA = m;
    ScalarType temp;
    bool BadArgument = false;
    bool LSide = (toupper(Side) == 'L');
    bool NoUnit = (toupper(Diag) == 'N');
    bool Upper = (toupper(Uplo) == 'U');

    if(!LSide) { NRowA = n; }

    // Quick return.
    if (n==izero || m==izero) { return; }

    if(!( LSide || (toupper(Side) == 'R'))) {
      cout << "BLAS::TRMM Error: SIDE == " << Side << endl;
      BadArgument = true; }
    if(!((toupper(Uplo) == 'L') || Upper )) {
      cout << "BLAS::TRMM Error: UPLO == " << Uplo << endl;
      BadArgument = true; }
    if(!((toupper(TransA) == 'N') || (toupper(TransA) == 'T') || (toupper(TransA) == 'C'))) {
      cout << "BLAS::TRMM Error: TRANSA == " << TransA << endl;
      BadArgument = true; }
    if(!( NoUnit || (toupper(Diag) == 'U'))) {
      cout << "BLAS::TRMM Error: DIAG == " << Diag << endl;
      BadArgument = true; }
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
	
	if ( toupper(TransA)=='N' ) {
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
	
	if( toupper(TransA) == 'N' ) {
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
	} // end if( toupper(TransA) == 'N' ) ...
      } // end if ( LSide ) ...
    } // end if (!BadArgument)
  } // end TRMM
  
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::TRSM(char Side, char Uplo, char TransA, char Diag, OrdinalType m, OrdinalType n, ScalarType alpha, ScalarType* A, OrdinalType lda, ScalarType* B, OrdinalType ldb)
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    ScalarType zero = ScalarTraits<ScalarType>::zero();
    ScalarType one = ScalarTraits<ScalarType>::one();
    ScalarType temp;
    OrdinalType NRowA = m;
    bool BadArgument = false;
    bool NoUnit = (toupper(Diag)=='N');
    
    if (!(toupper(Side) == 'L')) { NRowA = n; }

    // Quick return.
    if (n == izero || m == izero) { return; }

    if(!((toupper(Side) == 'L') || (toupper(Side) == 'R'))) {
      cout << "BLAS::TRSM Error: SIDE == " << Side << endl;
      BadArgument = true; }
    if(!((toupper(Uplo) == 'L') || (toupper(Uplo) == 'U'))) {
      cout << "BLAS::TRSM Error: UPLO == " << Uplo << endl;
      BadArgument = true; }
    if(!((toupper(TransA) == 'N') || (toupper(TransA) == 'T') || (TransA == 'C'))) {
      cout << "BLAS::TRSM Error: TRANSA == " << TransA << endl;
      BadArgument = true; }
    if(!((toupper(Diag) == 'N') || (toupper(Diag) == 'U'))) {
      cout << "BLAS::TRSM Error: DIAG == " << Diag << endl;
      BadArgument = true; }
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
	    if(toupper(Side) == 'L') {
		//
	    	// Perform computations for OP(A)*X = alpha*B	    
		//
		if(toupper(TransA) == 'N') {
		    //
		    //  Compute B = alpha*inv( A )*B
		    //
		    if(toupper(Uplo) == 'U') { 
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
		    } // end if (Uplo == 'U')
		}  // if (TransA =='N')	
	    	else { 
		    //
		    //  Compute B = alpha*inv( A' )*B
		    //
		    if(Uplo == 'U') { 
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
	    }  // if (Side == 'L')
	    else { 
	       // Side == 'R'
	       //
	       // Perform computations for X*OP(A) = alpha*B	    
	       //
	      if (toupper(TransA) == 'N') {
		    //
		    //  Compute B = alpha*B*inv( A )
		    //
		    if(toupper(Uplo) == 'U') { 
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
		    } // end if (Uplo == 'U')
		}  // if (TransA =='N')	
	    	else { 
		    //
		    //  Compute B = alpha*B*inv( A' )
		    //
		    if(Uplo == 'U') { 
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
  
#if 0

  template<typename OrdinalType>
  class BLAS<OrdinalType, float>
  {    
  public:
    inline BLAS(void) {};
    inline BLAS(const BLAS& BLAS_source) {};
    inline virtual ~BLAS(void) {};
    void ROTG(float da, float db, float c, float s);
    float ASUM(OrdinalType, float*, OrdinalType);
    void AXPY(OrdinalType, float, float*, OrdinalType, float*, OrdinalType);
    void COPY(OrdinalType, float*, OrdinalType, float*, OrdinalType);
    float DOT(OrdinalType, float*, OrdinalType, float*, OrdinalType);
    OrdinalType IAMAX(OrdinalType, float*, OrdinalType);   
    float NRM2(OrdinalType, float*, OrdinalType);
    void SCAL(OrdinalType, float, float*, OrdinalType);
    void GEMV(char, OrdinalType, OrdinalType, float, float*, OrdinalType, float*, OrdinalType, float, float*, OrdinalType);
    void GER(OrdinalType, OrdinalType, float, float*, OrdinalType, float*, OrdinalType, float*, OrdinalType);
    void TRMV(char, char, char, OrdinalType, float*, OrdinalType, float*, OrdinalType);
    void GEMM(char, char, OrdinalType, OrdinalType, OrdinalType, float, float*, OrdinalType, float*, OrdinalType, float, float*, OrdinalType);
    void SYMM(char, char, OrdinalType, OrdinalType, float, float*, OrdinalType, float*, OrdinalType, float, float*, OrdinalType);
    void TRMM(char, char, char, char, OrdinalType, OrdinalType, float, float*, OrdinalType, float*, OrdinalType);
    void TRSM(char, char, char, char, OrdinalType, OrdinalType, float, float*, OrdinalType, float*, OrdinalType);
  };

  template<typename OrdinalType>
  float BLAS<OrdinalType, float>::ROTG(float da, float db, float c, float s);
  { SROTG_F77(&da, &db, &c, &s ); }

  template<typename OrdinalType>
  float BLAS<OrdinalType, float>::ASUM(OrdinalType n, float* x, OrdinalType incx)
  { return SASUM_F77(&n, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::AXPY(OrdinalType n, float alpha, float* x, OrdinalType incx, float* y, OrdinalType incy)
  { SAXPY_F77(&n, &alpha, x, &incx, y, &incy); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::COPY(OrdinalType n, float* x, OrdinalType incx, float* y, OrdinalType incy)
  { SCOPY_F77(&n, x, &incx, y, &incy); }
  
  template<typename OrdinalType>
  float BLAS<OrdinalType, float>::DOT(OrdinalType n, float* x, OrdinalType incx, float* y, OrdinalType incy)
  { return SDOT_F77(&n, x, &incx, y, &incy); }
  
  template<typename OrdinalType>
  OrdinalType BLAS<OrdinalType, float>::IAMAX(OrdinalType n, float* x, OrdinalType incx)
  { return ISAMAX_F77(&n, x, &incx); }

  template<typename OrdinalType>
  float BLAS<OrdinalType, float>::NRM2(OrdinalType n, float* x, OrdinalType incx)
  { return SNRM2_F77(&n, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::SCAL(OrdinalType n, float alpha, float* x, OrdinalType incx)
  { SSCAL_F77(&n, &alpha, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::GEMV(char trans, OrdinalType m, OrdinalType n, float alpha, float* A, OrdinalType lda, float* x, OrdinalType incx, float beta, float* y, OrdinalType incy)
  { SGEMV_F77(&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::GER(OrdinalType m, OrdinalType n, float alpha, float* x, OrdinalType incx, float* y, OrdinalType incy, float* A, OrdinalType lda)
  { SGER_F77(&m, &n, &alpha, x, &incx, y, &incy, A, &lda); }

  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::TRMV(char uplo, char trans, char diag, OrdinalType n, float* A, OrdinalType lda, float* x, OrdinalType incx)
  { STRMV_F77(&uplo, &trans, &diag, &n, A, &lda, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::GEMM(char transa, char transb, OrdinalType m, OrdinalType n, OrdinalType k, float alpha, float* A, OrdinalType lda, float* B, OrdinalType ldb, float beta, float* C, OrdinalType ldc)
  { SGEMM_F77(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::SYMM(char side, char uplo, OrdinalType m, OrdinalType n, float alpha, float* A, OrdinalType lda, float* B, OrdinalType ldb, float beta, float* C, OrdinalType ldc)
  { SSYMM_F77(&side, &uplo, &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::TRMM(char side, char uplo, char transa, char diag, OrdinalType m, OrdinalType n, float alpha, float* A, OrdinalType lda, float* B, OrdinalType ldb)
  { STRMM_F77(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::TRSM(char side, char uplo, char transa, char diag, OrdinalType m, OrdinalType n, float alpha, float* A, OrdinalType lda, float* B, OrdinalType ldb)
  { STRSM_F77(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb); }

#endif

  template<typename OrdinalType>
  class BLAS<OrdinalType, double>
  {    
  public:
    inline BLAS(void) {};
    inline BLAS(const BLAS& BLAS_source) {};
    inline virtual ~BLAS(void) {};
    void ROTG(double da, double db, double c, double s);
    double ASUM(OrdinalType n, double* x, OrdinalType incx);
    void AXPY(OrdinalType n, double alpha, double* x, OrdinalType incx, double* y, OrdinalType incy);
    void COPY(OrdinalType n, double* x, OrdinalType incx, double* y, OrdinalType incy);
    double DOT(OrdinalType n, double* x, OrdinalType incx, double* y, OrdinalType incy);
    double NRM2(OrdinalType n, double* x, OrdinalType incx);
    void SCAL(OrdinalType n, double alpha, double* x, OrdinalType incx);
    OrdinalType IAMAX(OrdinalType n, double* x, OrdinalType incx);
    void GEMV(char trans, OrdinalType m, OrdinalType n, double alpha, double* A, OrdinalType lda, double* x, OrdinalType incx, double beta, double* y, OrdinalType incy);
    void TRMV(char uplo, char trans, char diag, OrdinalType n, double* A, OrdinalType lda, double* x, OrdinalType incx);
    void GER(OrdinalType m, OrdinalType n, double alpha, double* x, OrdinalType incx, double* y, OrdinalType incy, double* A, OrdinalType lda);
    void GEMM(char transa, char transb, OrdinalType m, OrdinalType n, OrdinalType k, double alpha, double* A, OrdinalType lda, double* B, OrdinalType ldb, double beta, double* C, OrdinalType ldc);
    void SYMM(char side, char uplo, OrdinalType m, OrdinalType n, double alpha, double* A, OrdinalType lda, double *B, OrdinalType ldb, double beta, double *C, OrdinalType ldc);
    void TRMM(char side, char uplo, char transa, char diag, OrdinalType m, OrdinalType n, double alpha, double* A, OrdinalType lda, double* B, OrdinalType ldb);
    void TRSM(char side, char uplo, char transa, char diag, OrdinalType m, OrdinalType n, double alpha, double* A, OrdinalType lda, double* B, OrdinalType ldb);
  };
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::ROTG(double da, double db, double c, double s)
  { DROTG_F77(&da, &db, &c, &s); }

  template<typename OrdinalType>
  double BLAS<OrdinalType, double>::ASUM(OrdinalType n, double* x, OrdinalType incx)
  { return DASUM_F77(&n, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::AXPY(OrdinalType n, double alpha, double* x, OrdinalType incx, double* y, OrdinalType incy)
  { DAXPY_F77(&n, &alpha, x, &incx, y, &incy); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::COPY(OrdinalType n, double* x, OrdinalType incx, double* y, OrdinalType incy)
  { DCOPY_F77(&n, x, &incx, y, &incy); }
  
  template<typename OrdinalType>
  double BLAS<OrdinalType, double>::DOT(OrdinalType n, double* x, OrdinalType incx, double* y, OrdinalType incy)
  { return DDOT_F77(&n, x, &incx, y, &incy); }
  
  template<typename OrdinalType>
  OrdinalType BLAS<OrdinalType, double>::IAMAX(OrdinalType n, double* x, OrdinalType incx)
  { return IDAMAX_F77(&n, x, &incx); }

  template<typename OrdinalType>
  double BLAS<OrdinalType, double>::NRM2(OrdinalType n, double* x, OrdinalType incx)
  { return DNRM2_F77(&n, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::SCAL(OrdinalType n, double alpha, double* x, OrdinalType incx)
  { DSCAL_F77(&n, &alpha, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::GEMV(char trans, OrdinalType m, OrdinalType n, double alpha, double* A, OrdinalType lda, double* x, OrdinalType incx, double beta, double* y, OrdinalType incy)
  { DGEMV_F77(&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::GER(OrdinalType m, OrdinalType n, double alpha, double* x, OrdinalType incx, double* y, OrdinalType incy, double* A, OrdinalType lda)
  { DGER_F77(&m, &n, &alpha, x, &incx, y, &incy, A, &lda); }

  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::TRMV(char uplo, char trans, char diag, OrdinalType n, double* A, OrdinalType lda, double* x, OrdinalType incx)
  { DTRMV_F77(&uplo, &trans, &diag, &n, A, &lda, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::GEMM(char transa, char transb, OrdinalType m, OrdinalType n, OrdinalType k, double alpha, double* A, OrdinalType lda, double* B, OrdinalType ldb, double beta, double* C, OrdinalType ldc)
  { DGEMM_F77(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::SYMM(char side, char uplo, OrdinalType m, OrdinalType n, double alpha, double* A, OrdinalType lda, double *B, OrdinalType ldb, double beta, double *C, OrdinalType ldc)
  { DSYMM_F77(&side, &uplo, &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::TRMM(char side, char uplo, char transa, char diag, OrdinalType m, OrdinalType n, double alpha, double* A, OrdinalType lda, double* B, OrdinalType ldb)
  { DTRMM_F77(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb); }

  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::TRSM(char side, char uplo, char transa, char diag, OrdinalType m, OrdinalType n, double alpha, double* A, OrdinalType lda, double* B, OrdinalType ldb)
  { DTRSM_F77(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb); }
  
} // namespace Teuchos

#endif // _TEUCHOS_BLAS_HPP_
