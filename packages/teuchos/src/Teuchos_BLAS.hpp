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
    void SCAL(OrdinalType n, ScalarType alpha, ScalarType* x, OrdinalType incx);
    void COPY(OrdinalType n, ScalarType* x, OrdinalType incx, ScalarType* y, OrdinalType incy);
    void AXPY(OrdinalType n, ScalarType alpha, ScalarType* x, OrdinalType incx, ScalarType* y, OrdinalType incy);
    ScalarType ASUM(OrdinalType n, ScalarType* x, OrdinalType incx);
    ScalarType DOT(OrdinalType n, ScalarType* x, OrdinalType incx, ScalarType* y, OrdinalType incy);
    ScalarType NRM2(OrdinalType n, ScalarType* x, OrdinalType incx);
    OrdinalType IAMAX(OrdinalType n, ScalarType* x, OrdinalType incx);
    //@}

    //@{ \name Level 2 BLAS Routines.
    void GEMV(char Trans, OrdinalType m, OrdinalType n, ScalarType alpha, ScalarType* A, OrdinalType lda,
                ScalarType* x, OrdinalType incx, ScalarType beta, ScalarType* y, OrdinalType incy);
    void GER(int m, int n, ScalarType alpha, ScalarType* x, int incx, ScalarType* y, int incy,
                ScalarType* A, int lda);
    void TRMV(char Uplo, char Trans, char Diag, int n, ScalarType* A, int lda, ScalarType* x, int incx);
    //@}
    
    //@{ \name Level 3 BLAS Routines. 
    void GEMM(char TransA, char TransB, int m, int n, int k, ScalarType alpha, ScalarType* A, int lda,
                ScalarType* B, int ldb, ScalarType beta, ScalarType* C, int ldc);
    void SYMM(char Side, char Uplo, int m, int n, ScalarType alpha, ScalarType* A, int lda,
                ScalarType* B, int ldb, ScalarType beta, ScalarType* C, int ldc);
    void TRMM(char Side, char Uplo, char TransA, char Diag, int m, int n,
                ScalarType alpha, ScalarType* A, int lda, ScalarType* B, int ldb);
    void TRSM(char Side, char Uplo, char TransA, char Diag, int m, int n,
                ScalarType alpha, ScalarType* A, int lda, ScalarType* B, int ldb);
    //@}
  };

//------------------------------------------------------------------------------------------
//      LEVEL 1 BLAS ROUTINES  
//------------------------------------------------------------------------------------------
    
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
  }

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
  }

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
  }

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
  }
  
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
  }
  
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
  }
  
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
  }

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
    if( incx == 0 ) {
	cout << "BLAS::GEMV Error: INCX == 0"<< endl;
	BadArgument = true;
    }
    if( incy == 0 ) {
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
  }

    
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::GER(int m, int n, ScalarType alpha, ScalarType* x, int incx, ScalarType* y, int incy, ScalarType* A, int lda)
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    int i, j;
    for(i = 0; i < m; i++) {
      for(j = 0; j < n; j++) {
	A[(j * m) + i] += (alpha * x[i] * y[j]); }}
  }
  
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::TRMV(char Uplo, char Trans, char Diag, int n, ScalarType* A, int lda, ScalarType* x, int incx)
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    bool BadArgument = 0;
    if(!((Uplo == 'U') || (Uplo == 'L'))) {
      cout << "BLAS::TRMV Error: TRIANGULAR == " << Uplo << endl;
      BadArgument = 1; }
    if(!((Trans == 'N') || (Trans == 'T') || (Trans == 'C'))) {
      cout << "BLAS::TRMV Error: TRANSPOSE == " << Trans << endl;
      BadArgument = 1; }
    if(!((Diag == 'N') || (Diag == 'U'))) {
      cout << "BLAS::TRMV Error: UNITDIAGONAL == " << Diag << endl;
      BadArgument = 1; }
    if(!BadArgument) {
      ScalarType* temp = new ScalarType[n];
      int i, j;
      for(i = 0; i < n; i++) {
	temp[i] = ScalarTraits<ScalarType>::zero(); }
      if((Uplo == 'U') || ((Uplo == 'L') && ((Trans == 'T') || (Trans == 'C')))) {
	if(Diag == 'N') {
	  for(i = 0; i < n; i++) {
	    for(j = i; j < n; j++) {
	      temp[i] += A[(j * n) + i] * x[j]; }}}
	else if(Diag == 'U') {
	  for(i = 0; i < n; i++) {
	    for(j = i; j < n; j++) {
	      if(i == j) {
		temp[i] += x[j]; }
	      else {
		temp[i] += A[(j * n) + i] * x[j]; }}}}}
      else if((Uplo == 'L') || ((Uplo == 'U') && ((Trans == 'T') || (Trans == 'C')))) {
	if(Diag == 'N') {
	  for(i = 0; i < n; i++) {
	    for(j = i; j < n; j++) {
	      temp[i] += A[(j * n) + i] * x[j]; }}}
	else if(Diag == 'U') {
	  for(i = 0; i < n; i++) {
	    for(j = i; j < n; j++) {
	      if(i == j) {
		temp[i] += x[j]; }
	      else {
		temp[i] += A[(j * n) + i] * x[j]; }}}}}
      for(i = 0; i < n; i++) {
	x[i] = temp[i]; } 
      delete [] temp; }
  }
    
//------------------------------------------------------------------------------------------
//      LEVEL 3 BLAS ROUTINES
//------------------------------------------------------------------------------------------
        
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::GEMM(char TransA, char TransB, int m, int n, int p, ScalarType alpha, ScalarType* A, int lda, ScalarType* B, int ldb, ScalarType beta, ScalarType* C, int ldc)
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    bool BadArgument = 0;
    if(!((TransA == 'N') || (TransA == 'T') || (TransA == 'C'))) {
      cout << "BLAS::GEMM Error: TRANSA == " << TransA << endl;
      BadArgument = 1; }
    if(!((TransB == 'N') || (TransB == 'T') || (TransB == 'C'))) {
      cout << "BLAS::GEMM Error: TRANSB == " << TransB << endl;
      BadArgument = 1; }
    if(!BadArgument) {
      ScalarType zero = ScalarTraits<ScalarType>::zero();
      ScalarType one = ScalarTraits<ScalarType>::one();
      int i, j, k;
      if(beta == zero) {
	for(i = 0; i < (m * n); i++) {
	  C[i] = zero; }}
      else {
	for(i = 0; i < (m * n); i++) {
	  C[i] *= beta; }}
      if(alpha != zero) {
	if(alpha == one) {
	  if(TransA == 'N') {
	    if(TransB == 'N') {
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < p; k++) {
		    C[i + (j * m)] += (A[i + (k * m)] * B[k + (j * p)]); }}}}
	    else { // TransB == T || C
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < p; k++) {
		    C[i + (j * m)] += (A[i + (k * m)] * B[j + (k * n)]); }}}}}
	  else { // TransA == T || C
	    if(TransB == 'N') {
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < p; k++) {
		    C[i + (j * m)] += (A[k + (i * p)] * B[k + (j * p)]); }}}}
	    else { // TransB == T || C
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < p; k++) {
		    C[i + (j * m)] += (A[k + (i * p)] * B[j + (k * n)]); }}}}}}
	else { // alpha is neither 1 nor 0 
	  if(TransA == 'N') {
	    if(TransB == 'N') {
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < p; k++) {
		    C[i + (j * m)] += (alpha * A[i + (k * m)] * B[k + (j * p)]); }}}}
	    else { // TransB == T || C
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < p; k++) {
		    C[i + (j * m)] += (alpha * A[i + (k * m)] * B[j + (k * n)]); }}}}}
	  else { // TransA == T || C
	    if(TransB == 'N') {
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < p; k++) {
		    C[i + (j * m)] += (alpha * A[k + (i * p)] * B[k + (j * p)]); }}}}
	    else { // TransB == T || C
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < p; k++) {
		    C[i + (j * m)] += (alpha * A[k + (i * p)] * B[j + (k * n)]); }}}}}}}}
  }

  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::SYMM(char Side, char Uplo, int m, int n, ScalarType alpha, ScalarType* A, int lda, ScalarType* B, int ldb, ScalarType beta, ScalarType* C, int ldc)
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    ScalarType zero = ScalarTraits<ScalarType>::zero();
    ScalarType one = ScalarTraits<ScalarType>::one();
    bool BadArgument = 0;
    if(!((Side == 'L') || (Side == 'R'))) {
      cout << "BLAS::GEMM Error: SIDE == " << Side << endl;
      BadArgument = 1; }
    if(!((Uplo == 'L') || (Uplo == 'U'))) {
      cout << "BLAS::GEMM Error: UPLO == " << Uplo << endl;
      BadArgument = 1; }
    if(!BadArgument) {
      int i, j, k;
      if(beta == zero) {
	for(i = 0; i < (m * n); i++) {
	  C[i] = zero; }} 
      else {
	for(i = 0; i < (m * n); i++) {
	  C[i] *= beta; }}
      if(alpha != zero) {
	if(alpha == one) {
	  if(Side == 'L') {
	    if(Uplo == 'U') {
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < m; k++) {
		    if(k < i) {
		      C[i + (j * m)] += (A[k + (i * m)] * B[k + (j * m)]); }
		    else {
		      C[i + (j * m)] += (A[i + (k * m)] * B[k + (j * m)]); }}}}}
	    else { // Uplo == L
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < m; k++) {
		    if(k < i) {
		      C[i + (j * m)] += (A[i + (k * m)] * B[k + (j * m)]); }
		    else {
		      C[i + (j * m)] += (A[k + (i * m)] * B[k + (j * m)]); }}}}}}
	  else { // Side == R
	    if(Uplo == 'U') {
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < n; k++) {
		    if(k < j) {
		      C[i + (j * m)] += (B[i + (k * m)] * A[k + (j * n)]); }
		    else {
		      C[i + (j * m)] += (B[i + (k * m)] * A[j + (k * n)]); }}}}}
	    else { // Uplo == L
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < n; k++) {
		    if(k < j) {
		      C[i + (j * m)] += (B[i + (k * m)] * A[j + (k * n)]); } 
		    else {
		      C[i + (j * m)] += (B[i + (k * m)] * A[k + (j * n)]); }}}}}}}
	else { // alpha is neither 1 nor 0
	  if(Side == 'L') {
	    if(Uplo == 'U') {
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < m; k++) {
		    if(k < i) {
		      C[i + (j * m)] += (alpha * A[k + (i * m)] * B[k + (j * m)]); }
		    else {
		      C[i + (j * m)] += (alpha * A[i + (k * m)] * B[k + (j * m)]); }}}}}
	    else { // Uplo == L
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < m; k++) {
		    if(k < i) {
		      C[i + (j * m)] += (alpha * A[i + (k * m)] * B[k + (j * m)]); }
		    else {
		      C[i + (j * m)] += (alpha * A[k + (i * m)] * B[k + (j * m)]); }}}}}}
	  else { // Side == R
	    if(Uplo == 'U') {
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < n; k++) {
		    if(k < j) {
		      C[i + (j * m)] += (alpha * B[i + (k * m)] * A[k + (j * n)]); }
		    else {
		      C[i + (j * m)] += (alpha * B[i + (k * m)] * A[j + (k * n)]); }}}}}
	    else { // Uplo == L
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < n; k++) { 
		    if(k < j) { 
		      C[i + (j * m)] += (alpha * B[i + (k * m)] * A[j + (k * n)]); }
		    else {
		      C[i + (j * m)] += (alpha * B[i + (k * m)] * A[k + (j * n)]); }}}}}}}}}
  }
  
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::TRMM(char Side, char Uplo, char TransA, char Diag, int m, int n, ScalarType alpha, ScalarType* A, int lda, ScalarType* B, int ldb)
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    ScalarType zero = ScalarTraits<ScalarType>::zero();
    ScalarType one = ScalarTraits<ScalarType>::one();
    bool BadArgument = 0;
    if(!((Side == 'L') || (Side == 'R'))) {
      cout << "BLAS::GEMM Error: SIDE == " << Side << endl;
      BadArgument = 1; }
    if(!((Uplo == 'L') || (Uplo == 'U'))) {
      cout << "BLAS::GEMM Error: UPLO == " << Uplo << endl;
      BadArgument = 1; }
    if(!((TransA == 'N') || (TransA == 'T') || (TransA == 'C'))) {
      cout << "BLAS::GEMM Error: TRANSA == " << TransA << endl;
      BadArgument = 1; }
    if(!((Diag == 'N') || (Diag == 'U'))) {
      cout << "BLAS::GEMM Error: DIAG == " << Diag << endl;
      BadArgument = 1; }
    if(!BadArgument) {
      int i, j, k;
      if(alpha == zero) {
	for(i = 0; i < (m * n); i++) {
	  B[i] = zero; }}
      else {
	ScalarType* temp = new ScalarType[m * n];
	for(i = 0; i < (m * n); i++) {
	  temp[i] = B[i];
	  B[i] = zero; }
	if(alpha == one) {
	  if(Side == 'L') {
	    if(Uplo == 'U') {
	      if(TransA == 'N') {
		if(Diag == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = i; k < m; k++) {
			B[i + (j * m)] += (A[i + (k * m)] * temp[k + (j * m)]); }}}}
		else { // Diag == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = i; k < m; k++) {
			if(i == k) {
			  B[i + (j * m)] += temp[k + (j * m)]; }
			else {
			  B[i + (j * m)] += (A[i + (k * m)] * temp[k + (j * m)]); }}}}}}
	      else { // TransA == T or C
		if(Diag == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (i + 1); k++) {
			B[i + (j * m)] += (A[k + (i * m)] * temp[k + (j * m)]); }}}}
		else { // Diag == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (i + 1); k++) {
			if(i == k) {
			  B[i + (j * m)] += temp[k + (j * m)]; }
			else {
			  B[i + (j * m)] += (A[k + (i * m)] * temp[k + (j * m)]); }}}}}}}
	    else { // Uplo == L
	      if(TransA == 'N') {
		if(Diag == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (i + 1); k++) { 
			B[i + (j * m)] += (A[i + (k * m)] * temp[k + (j * m)]); }}}}
		else { // Diag == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (i + 1); k++) {
			if(i == k) {
			  B[i + (j * m)] += temp[k + (j * m)]; }
			else {
			  B[i + (j * m)] += (A[i + (k * m)] * temp[k + (j * m)]); }}}}}}
	      else { // TransA == T or C
		if(Diag == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = i; k < m; k++) {		    
			B[i + (j * m)] += (A[k + (i * m)] * temp[k + (j * m)]); }}}}
		else { // Diag == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = i; k < m; k++) {
			if(i == k) {
			  B[i + (j * m)] += temp[k + (j * m)]; }
			else {
			  B[i + (j * m)] += (A[k + (i * m)] * temp[k + (j * m)]); }}}}}}}}
	  else { // Side == R
	    if(Uplo == 'U') {
	      if(TransA == 'N') {
		if(Diag == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (j + 1); k++) {
			B[i + (j * m)] += (temp[i + (k * m)] * A[k + (j * n)]); }}}}
		else { // Diag == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (j + 1); k++) {
			if(j == k) {
			  B[i + (j * m)] += temp[i + (k * m)]; }
			else {
			  B[i + (j * m)] += (temp[i + (k * m)] * A[k + (j * n)]); }}}}}}
	      else { // TransA == T or C
		if(Diag == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = j; k < n; k++) {
			B[i + (j * m)] += (temp[i + (k * m)] * A[j + (k * n)]); }}}}
		else { // Diag == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = j; k < n; k++) {
			if(j == k) {
			  B[i + (j * m)] += temp[i + (k * m)]; }
			else {
			  B[i + (j * m)] += (temp[i + (k * m)] * A[j + (k * n)]); }}}}}}}
	    else { // Uplo == L
	      if(TransA == 'N') {
		if(Diag == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = j; k < n; k++) {
			B[i + (j * m)] += (temp[i + (k * m)] * A[k + (j * n)]); }}}}
		else { // Diag == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = j; k < n; k++) {
			if(j == k) {
			  B[i + (j * m)] += temp[i + (k * m)]; }
			else {
			  B[i + (j * m)] += (temp[i + (k * m)] * A[k + (j * n)]); }}}}}}
	      else { // TransA == T or C
		if(Diag == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (j + 1); k++) {
			B[i + (j * m)] += (temp[i + (k * m)] * A[j + (k * n)]); }}}}
		else { // Diag == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (j + 1); k++) {
			if(j == k) {
			  B[i + (j * m)] += temp[i + (k * m)]; } 
			else {
			  B[i + (j * m)] += (temp[i + (k * m)] * A[j + (k * n)]); }}}}}}}}}
	else // alpha is neither 0 nor 1
	  {
	  if(Side == 'L') {
	    if(Uplo == 'U') {
	      if(TransA == 'N') {
		if(Diag == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = i; k < m; k++) {
			B[i + (j * m)] += (alpha * A[i + (k * m)] * temp[k + (j * m)]); }}}}
		else { // Diag == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = i; k < m; k++) {
			if(i == k) {
			  B[i + (j * m)] += (alpha * temp[k + (j * m)]); }
			else {
			  B[i + (j * m)] += (alpha * A[i + (k * m)] * temp[k + (j * m)]); }}}}}}
	      else { // TransA == T or C
		if(Diag == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (i + 1); k++) {
			B[i + (j * m)] += (alpha * A[k + (i * m)] * temp[k + (j * m)]); }}}}
		else { // Diag == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (i + 1); k++) {
			if(i == k) {
			  B[i + (j * m)] += (alpha * temp[k + (j * m)]); }
			else {
			  B[i + (j * m)] += (alpha * A[k + (i * m)] * temp[k + (j * m)]); }}}}}}}
	    else { // Uplo == L
	      if(TransA == 'N') {
		if(Diag == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (i + 1); k++) { 
			B[i + (j * m)] += (alpha * A[i + (k * m)] * temp[k + (j * m)]); }}}}
		else { // Diag == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (i + 1); k++) {
			if(i == k) {
			  B[i + (j * m)] += (alpha * temp[k + (j * m)]); }
			else {
			  B[i + (j * m)] += (alpha * A[i + (k * m)] * temp[k + (j * m)]); }}}}}}
	      else { // TransA == T or C
		if(Diag == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = i; k < m; k++) {		    
			B[i + (j * m)] += (alpha * A[k + (i * m)] * temp[k + (j * m)]); }}}}
		else { // Diag == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = i; k < m; k++) {
			if(i == k) {
			  B[i + (j * m)] += (alpha * temp[k + (j * m)]); }
			else {
			  B[i + (j * m)] += (alpha * A[k + (i * m)] * temp[k + (j * m)]); }}}}}}}}
	  else { // Side == R
	    if(Uplo == 'U') {
	      if(TransA == 'N') {
		if(Diag == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (j + 1); k++) {
			B[i + (j * m)] += (alpha * temp[i + (k * m)] * A[k + (j * n)]); }}}}
		else { // Diag == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (j + 1); k++) {
			if(j == k) {
			  B[i + (j * m)] += (alpha * temp[i + (k * m)]); }
			else {
			  B[i + (j * m)] += (alpha * temp[i + (k * m)] * A[k + (j * n)]); }}}}}}
	      else { // TransA == T or C
		if(Diag == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = j; k < n; k++) {
			B[i + (j * m)] += (alpha * temp[i + (k * m)] * A[j + (k * n)]); }}}}
		else { // Diag == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = j; k < n; k++) {
			if(j == k) {
			  B[i + (j * m)] += (alpha * temp[i + (k * m)]); }
			else {
			  B[i + (j * m)] += (alpha * temp[i + (k * m)] * A[j + (k * n)]); }}}}}}}
	    else { // Uplo == L
	      if(TransA == 'N') {
		if(Diag == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = j; k < n; k++) {
			B[i + (j * m)] += (alpha * temp[i + (k * m)] * A[k + (j * n)]); }}}}
		else { // Diag == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = j; k < n; k++) {
			if(j == k) {
			  B[i + (j * m)] += (alpha * temp[i + (k * m)]); }
			else {
			  B[i + (j * m)] += (alpha * temp[i + (k * m)] * A[k + (j * n)]); }}}}}}
	      else { // TransA == T or C
		if(Diag == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (j + 1); k++) {
			B[i + (j * m)] += (alpha * temp[i + (k * m)] * A[j + (k * n)]); }}}}
		else { // Diag == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (j + 1); k++) {
			if(j == k) {
			  B[i + (j * m)] += (alpha * temp[i + (k * m)]); } 
			else {
			  B[i + (j * m)] += (alpha * temp[i + (k * m)] * A[j + (k * n)]); }}}}}}}}}
	delete [] temp; }}
  }
  
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::TRSM(char Side, char Uplo, char TransA, char Diag, int m, int n, ScalarType alpha, ScalarType* A, int lda, ScalarType* B, int ldb)
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    ScalarType zero = ScalarTraits<ScalarType>::zero();
    ScalarType one = ScalarTraits<ScalarType>::one();
    ScalarType temp;
    bool BadArgument = 0;
    if(!((Side == 'L') || (Side == 'R'))) {
      cout << "BLAS::GEMM Error: SIDE == " << Side << endl;
      BadArgument = 1; }
    if(!((Uplo == 'L') || (Uplo == 'U'))) {
      cout << "BLAS::GEMM Error: UPLO == " << Uplo << endl;
      BadArgument = 1; }
    if(!((TransA == 'N') || (TransA == 'T') || (TransA == 'C'))) {
      cout << "BLAS::GEMM Error: TRANSA == " << TransA << endl;
      BadArgument = 1; }
    if(!((Diag == 'N') || (Diag == 'U'))) {
      cout << "BLAS::GEMM Error: DIAG == " << Diag << endl;
      BadArgument = 1; }
    if(!BadArgument && n!=0)
      {
	int i, j, k;
	// Set the solution to the zero vector.
	if(alpha == zero) {
	    for(j = 0; j < n; j++) {
	    	for( i = 0; i < m; i++) {
		    B[j*ldb+i] = zero;
	      	}
	    }
	}
	else 
	{ // Start the operations.
	    if(Side == 'L') {
		//
	    	// Perform computations for OP(A)*X = alpha*B	    
		//
		if(TransA == 'N') {
		    //
		    //  Compute B = alpha*inv( A )*B
		    //
		    if(Uplo == 'U') { 
			// A is upper triangular.
			for(j = 0; j < n; j++) {
	    		    // Perform alpha*B if alpha is not 1.
	    		    if(alpha != one) {
	    	    		for( i = 0; i < m; i++) {
		    		    B[j*ldb+i] *= alpha;
		    		}
			    }
			    // Perform a backsolve for column j of B.
			    for(k = (m - 1); k > -1; k--) {
				// If this entry is zero, we don't have to do anything.
				if (B[j*ldb + k] != zero) {
				    if (Diag == 'N') {
					B[j*ldb + k] /= A[k*lda + k];
				    }
				    for(i = 0; i < k; i++) {
					B[j*ldb + i] -= B[j*ldb + k] * A[k*lda + i];
				    }
				}
			    }
			}
		    }
		    else 
		    { // A is lower triangular.
                        for(j = 0; j < n; j++) {
                            // Perform alpha*B if alpha is not 1.
                            if(alpha != one) {
                                for( i = 0; i < m; i++) {
                                    B[j*ldb+i] *= alpha;
                                }
                            }
                            // Perform a forward solve for column j of B.
                            for(k = 0; k < m; k++) {
                                // If this entry is zero, we don't have to do anything.
                                if (B[j*ldb + k] != zero) {   
                                    if (Diag == 'N') {
                                        B[j*ldb + k] /= A[k*lda + k];
                                    }
                                    for(i = k+1; i < m; i++) {
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
			for(j = 0; j < n; j++) {
	    	    	    for( i = 0; i < m; i++) {
		    		temp = alpha*B[j*ldb+i];
			    	for(k = 0; k < i; k++) {
				    temp -= A[i*lda + k] * B[j*ldb + k];
				}
				if (Diag == 'N') {
				    temp /= A[i*lda + i];
				}
				B[j*ldb + i] = temp;
			    }
			}
		    }
		    else
		    { // A is lower triangular.
                        for(j = 0; j < n; j++) {
                            for(i = (m - 1) ; i > -1; i--) {
                                temp = alpha*B[j*ldb+i];
                            	for(k = i+1; k < m; k++) {
				    temp -= A[i*lda + k] * B[j*ldb + k];
				}
				if (Diag == 'N') {
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
	      if (TransA == 'N') {
		    //
		    //  Compute B = alpha*B*inv( A )
		    //
		    if(Uplo == 'U') { 
			// A is upper triangular.
	    		// Perform a backsolve for column j of B.
			for(j = 0; j < n; j++) {
	    		    // Perform alpha*B if alpha is not 1.
	    		    if(alpha != one) {
	    	    		for( i = 0; i < m; i++) {
		    		    B[j*ldb+i] *= alpha;
		    		}
			    }
			    for(k = 0; k < j; k++) {
				// If this entry is zero, we don't have to do anything.
				if (A[j*lda + k] != zero) {
				    for(i = 0; i < m; i++) {
					B[j*ldb + i] -= A[j*lda + k] * B[k*ldb + i];
				    }
				}
			    }
			    if (Diag == 'N') {
				temp = one/A[j*lda + j];
				for(i = 0; i < m; i++) {
				    B[j*ldb + i] *= temp;
				}
			    }
			}
		    }
		    else 
		    { // A is lower triangular.
                        for(j = (n - 1); j > -1; j--) {
                            // Perform alpha*B if alpha is not 1.
                            if(alpha != one) {
                                for( i = 0; i < m; i++) {
                                    B[j*ldb+i] *= alpha;
                                }
                            }
                            // Perform a forward solve for column j of B.
                            for(k = j+1; k < n; k++) {
                                // If this entry is zero, we don't have to do anything.
				if (A[j*lda + k] != zero) {
				    for(i = 0; i < m; i++) {
                                        B[j*ldb + i] -= A[j*lda + k] * B[k*ldb + i]; 
                                    }
                                } 
                            }
			    if (Diag == 'N') {
				temp = one/A[j*lda + j];
				for(i = 0; i < m; i++) {
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
			for(k = (n - 1); k > -1; k--) {
			    if (Diag == 'N') {
				temp = one/A[k*lda + k];
	    	    	    	for(i = 0; i < m; i++) {
		    		    B[k*ldb + i] *= temp;
				}
			    }
			    for(j = 0; j < k; j++) {
				if (A[k*lda + j] != zero) {
				    temp = A[k*lda + j];
				    for(i = 0; i < m; i++) {
					B[j*ldb + i] -= temp*B[k*ldb + i];
				    }
				}
			    }
			    if (alpha != one) {
				for (i = 0; i < m; i++) {
				    B[k*ldb + i] *= alpha;
				}
			    }
			}
		    }
		    else
		    { // A is lower triangular.
			for(k = 0; k < n; k++) {
			    if (Diag == 'N') {
				temp = one/A[k*lda + k];
				for (i = 0; i < m; i++) {
				    B[k*ldb + i] *= temp;
				}
			    }
			    for(j = k+1; j < n; j++) {
				if(A[k*lda + j] != zero) {
				    temp = A[k*lda + j];
				    for(i = 0; i < m; i++) {
					B[j*ldb + i] -= temp*B[k*ldb + i];
				    }
				}
			    }
			    if (alpha != one) {
				for (i = 0; i < m; i++) {
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
