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
//             * All of the L2/L3 routines assume that the entire matrix is being used (that is, if A is mxn, LDA = m); they don't work on submatrices yet. This *should* be a reasonably trivial thing to fix, as well.
//          -- Removed warning messages for generic calls
// 08.08.03 -- TRSM now works for all cases where SIDE == L and DIAG == N. DIAG == U is implemented but does not work correctly; SIDE == R is not yet implemented.
// 08.14.03 -- TRSM now works for all cases and accepts (and uses) leading-dimension information.

#ifndef _TEUCHOS_BLAS_HPP_
#define _TEUCHOS_BLAS_HPP_

#include "Teuchos_BLAS_wrappers.hpp"
#include "Teuchos_ScalarTraits.hpp"

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
    inline BLAS(void) {};
    inline BLAS(const BLAS& BLAS_source) {};
    inline virtual ~BLAS(void) {};
    ScalarType ASUM(int, ScalarType*, int);
    void AXPY(int, ScalarType, ScalarType*, int, ScalarType*, int);
    void COPY(int, ScalarType*, int, ScalarType*, int);
    ScalarType DOT(int, ScalarType*, int, ScalarType*, int);
    int IAMAX(int, ScalarType*, int);    
    ScalarType NRM2(int, ScalarType*, int);
    void SCAL(int, ScalarType, ScalarType*, int);
    void GEMV(char, int, int, ScalarType, ScalarType*, int, ScalarType*, int, ScalarType, ScalarType*, int);
    void GER(int, int, ScalarType, ScalarType*, int, ScalarType*, int, ScalarType*, int);
    void TRMV(char, char, char, int, ScalarType*, int, ScalarType*, int);
    void GEMM(char, char, int, int, int, ScalarType, ScalarType*, int, ScalarType*, int, ScalarType, ScalarType*, int);
    void SYMM(char, char, int, int, ScalarType, ScalarType*, int, ScalarType*, int, ScalarType beta, ScalarType*, int);
    void TRMM(char, char, char, char, int, int, ScalarType, ScalarType*, int, ScalarType*, int);
    void TRSM(char, char, char, char, int, int, ScalarType, ScalarType*, int, ScalarType*, int);
    void XERBLA(char, int);
  };

  template<typename OrdinalType, typename ScalarType>
  ScalarType BLAS<OrdinalType, ScalarType>::ASUM(int n, ScalarType* x, int StepX)
  {
    ScalarType result = ScalarTraits<ScalarType>::zero();
    int i;
    for(i = 0; i < n; i += StepX)
      {
	result += ScalarTraits<ScalarType>::magnitude(x[i]);
      }
    return result;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::AXPY(int n, ScalarType alpha, ScalarType* x, int StepX, ScalarType* y, int StepY)
  {
    if((StepX == StepY) && (StepX > 0))
      {
	int i;
	for(i = 0; i < n; i += StepX)
	  {
	    y[i] += alpha * x[i];
	  }
      }
  }
  
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::COPY(int n, ScalarType* x, int StepX, ScalarType* y, int StepY)
  {
    if((StepX == StepY) && (StepX > 0))
      {
	int i;
	for(i = 0; i < n; i += StepX)
	  {
	    y[i] = x[i];
	  }
      }
  }
  
  template<typename OrdinalType, typename ScalarType>
  ScalarType BLAS<OrdinalType, ScalarType>::DOT(int n, ScalarType* x, int StepX, ScalarType* y, int StepY)
  {
    ScalarType result = ScalarTraits<ScalarType>::zero();
    if((StepX == StepY) && (StepX > 0))
      {
	int i;
	for(i = 0; i < n; i+= StepX)
	  {
	    result += x[i] * y[i];
	  }
      }
    return result;
  }
  
  template<typename OrdinalType, typename ScalarType>
  int BLAS<OrdinalType, ScalarType>::IAMAX(int n, ScalarType* x, int StepX)
  {
    int result = 0, i;
    for(i = 0 + StepX; i < n; i += StepX)
      {
	if(ScalarTraits<ScalarType>::magnitude(x[i]) > ScalarTraits<ScalarType>::magnitude(x[result]))
	  {
	    result = i;
	  }
      }
    return result + 1; // the BLAS I?AMAX functions return 1-indexed (Fortran-style) values
  }

  template<typename OrdinalType, typename ScalarType>
  ScalarType BLAS<OrdinalType, ScalarType>::NRM2(int n, ScalarType* x, int StepX)
  {
    ScalarType result = ScalarTraits<ScalarType>::zero();
    int i;
    for(i = 0; i < n; i += StepX)
      {
	result += x[i] * x[i];
      }
    result = ScalarTraits<ScalarType>::squareroot(result);
    return result;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::SCAL(int n, ScalarType alpha, ScalarType* x, int StepX)
  {
    int i;
    for(i = 0; i < n; i += StepX)
      {
	x[i] *= alpha;
      }
  }

  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::GEMV(char Transpose, int m, int n, ScalarType alpha, ScalarType* A, int LDA, ScalarType* x, int StepX, ScalarType beta, ScalarType* y, int StepY)
  {
    bool BadArgument = 0;
    if(!((Transpose == 'N') || (Transpose == 'T') || (Transpose == 'C'))) {
      cout << "BLAS::GEMM Error: TRANSA == " << Transpose << endl;
      BadArgument = 1; }
    if(!BadArgument) {
      ScalarType zero = ScalarTraits<ScalarType>::zero();
      ScalarType one = ScalarTraits<ScalarType>::one();
      int i, j;
      if(beta == zero) {
	for(i = 0; i < m; i++) {
	  y[i] = zero; }}
      else if(beta == one) {
	// do nothing
      }
      else {
	for(i = 0; i < m; i++) {
	  y[i] *= beta; }}
      if(Transpose == 'N') {
	if(alpha == zero) {
	  // do nothing
	}
	else if(alpha == one) {
	  for(i = 0; i < m; i++) {
	    for(j = 0; j < n; j++) {
	      y[i] += A[(j * m) + i] * x[j]; }}}
	else {
	  for(i = 0; i < m; i++) {
	    for(j = 0; j < n; j++) {
	      y[i] += alpha * (A[(j * m) + i] * x[j]); }}}}
      else if((Transpose == 'T') || (Transpose == 'C')) {
	if(alpha == zero) {
	  // do nothing 
	}
	else if(alpha == one) {
	  for(i = 0; i < m; i++) {
	    for(j = 0; j < n; j++) {
	      y[i] += (A[(i * n) + j] * x[j]); }}}
	else {
	  for(i = 0; i < m; i++) {
	    for(j = 0; j < n; j++) {
	      y[i] += alpha * (A[(i * n) + j] * x[j]); }}}}}
  }
    
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::GER(int m, int n, ScalarType alpha, ScalarType* x, int StepX, ScalarType* y, int StepY, ScalarType* A, int LDA)
  {
    int i, j;
    for(i = 0; i < m; i++) {
      for(j = 0; j < n; j++) {
	A[(j * m) + i] += (alpha * x[i] * y[j]); }}
  }
  
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::TRMV(char Triangular, char Transpose, char UnitDiagonal, int n, ScalarType* A, int LDA, ScalarType* x, int StepX)
  {
    bool BadArgument = 0;
    if(!((Triangular == 'U') || (Triangular == 'L'))) {
      cout << "BLAS::TRMV Error: TRIANGULAR == " << Triangular << endl;
      BadArgument = 1; }
    if(!((Transpose == 'N') || (Transpose == 'T') || (Transpose == 'C'))) {
      cout << "BLAS::TRMV Error: TRANSPOSE == " << Transpose << endl;
      BadArgument = 1; }
    if(!((UnitDiagonal == 'N') || (UnitDiagonal == 'U'))) {
      cout << "BLAS::TRMV Error: UNITDIAGONAL == " << UnitDiagonal << endl;
      BadArgument = 1; }
    if(!BadArgument) {
      ScalarType* temp = new ScalarType[n];
      int i, j;
      for(i = 0; i < n; i++) {
	temp[i] = ScalarTraits<ScalarType>::zero(); }
      if((Triangular == 'U') || ((Triangular == 'L') && ((Transpose == 'T') || (Transpose == 'C')))) {
	if(UnitDiagonal == 'N') {
	  for(i = 0; i < n; i++) {
	    for(j = i; j < n; j++) {
	      temp[i] += A[(j * n) + i] * x[j]; }}}
	else if(UnitDiagonal == 'U') {
	  for(i = 0; i < n; i++) {
	    for(j = i; j < n; j++) {
	      if(i == j) {
		temp[i] += x[j]; }
	      else {
		temp[i] += A[(j * n) + i] * x[j]; }}}}}
      else if((Triangular == 'L') || ((Triangular == 'U') && ((Transpose == 'T') || (Transpose == 'C')))) {
	if(UnitDiagonal == 'N') {
	  for(i = 0; i < n; i++) {
	    for(j = i; j < n; j++) {
	      temp[i] += A[(j * n) + i] * x[j]; }}}
	else if(UnitDiagonal == 'U') {
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
    
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::GEMM(char TransposeA, char TransposeB, int m, int n, int p, ScalarType alpha, ScalarType* A, int LDA, ScalarType* B, int LDB, ScalarType beta, ScalarType* C, int LDC)
  {
    bool BadArgument = 0;
    if(!((TransposeA == 'N') || (TransposeA == 'T') || (TransposeA == 'C'))) {
      cout << "BLAS::GEMM Error: TRANSA == " << TransposeA << endl;
      BadArgument = 1; }
    if(!((TransposeB == 'N') || (TransposeB == 'T') || (TransposeB == 'C'))) {
      cout << "BLAS::GEMM Error: TRANSB == " << TransposeB << endl;
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
	  if(TransposeA == 'N') {
	    if(TransposeB == 'N') {
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < p; k++) {
		    C[i + (j * m)] += (A[i + (k * m)] * B[k + (j * p)]); }}}}
	    else { // TransposeB == T || C
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < p; k++) {
		    C[i + (j * m)] += (A[i + (k * m)] * B[j + (k * n)]); }}}}}
	  else { // TransposeA == T || C
	    if(TransposeB == 'N') {
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < p; k++) {
		    C[i + (j * m)] += (A[k + (i * p)] * B[k + (j * p)]); }}}}
	    else { // TransposeB == T || C
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < p; k++) {
		    C[i + (j * m)] += (A[k + (i * p)] * B[j + (k * n)]); }}}}}}
	else { // alpha is neither 1 nor 0 
	  if(TransposeA == 'N') {
	    if(TransposeB == 'N') {
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < p; k++) {
		    C[i + (j * m)] += (alpha * A[i + (k * m)] * B[k + (j * p)]); }}}}
	    else { // TransposeB == T || C
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < p; k++) {
		    C[i + (j * m)] += (alpha * A[i + (k * m)] * B[j + (k * n)]); }}}}}
	  else { // TransposeA == T || C
	    if(TransposeB == 'N') {
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < p; k++) {
		    C[i + (j * m)] += (alpha * A[k + (i * p)] * B[k + (j * p)]); }}}}
	    else { // TransposeB == T || C
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < p; k++) {
		    C[i + (j * m)] += (alpha * A[k + (i * p)] * B[j + (k * n)]); }}}}}}}}
  }

  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::SYMM(char Side, char Triangular, int m, int n, ScalarType alpha, ScalarType* A, int LDA, ScalarType* B, int LDB, ScalarType beta, ScalarType* C, int LDC)
  {
    ScalarType zero = ScalarTraits<ScalarType>::zero();
    ScalarType one = ScalarTraits<ScalarType>::one();
    bool BadArgument = 0;
    if(!((Side == 'L') || (Side == 'R'))) {
      cout << "BLAS::GEMM Error: SIDE == " << Side << endl;
      BadArgument = 1; }
    if(!((Triangular == 'L') || (Triangular == 'U'))) {
      cout << "BLAS::GEMM Error: UPLO == " << Triangular << endl;
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
	    if(Triangular == 'U') {
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < m; k++) {
		    if(k < i) {
		      C[i + (j * m)] += (A[k + (i * m)] * B[k + (j * m)]); }
		    else {
		      C[i + (j * m)] += (A[i + (k * m)] * B[k + (j * m)]); }}}}}
	    else { // Triangular == L
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < m; k++) {
		    if(k < i) {
		      C[i + (j * m)] += (A[i + (k * m)] * B[k + (j * m)]); }
		    else {
		      C[i + (j * m)] += (A[k + (i * m)] * B[k + (j * m)]); }}}}}}
	  else { // Side == R
	    if(Triangular == 'U') {
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < n; k++) {
		    if(k < j) {
		      C[i + (j * m)] += (B[i + (k * m)] * A[k + (j * n)]); }
		    else {
		      C[i + (j * m)] += (B[i + (k * m)] * A[j + (k * n)]); }}}}}
	    else { // Triangular == L
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < n; k++) {
		    if(k < j) {
		      C[i + (j * m)] += (B[i + (k * m)] * A[j + (k * n)]); } 
		    else {
		      C[i + (j * m)] += (B[i + (k * m)] * A[k + (j * n)]); }}}}}}}
	else { // alpha is neither 1 nor 0
	  if(Side == 'L') {
	    if(Triangular == 'U') {
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < m; k++) {
		    if(k < i) {
		      C[i + (j * m)] += (alpha * A[k + (i * m)] * B[k + (j * m)]); }
		    else {
		      C[i + (j * m)] += (alpha * A[i + (k * m)] * B[k + (j * m)]); }}}}}
	    else { // Triangular == L
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < m; k++) {
		    if(k < i) {
		      C[i + (j * m)] += (alpha * A[i + (k * m)] * B[k + (j * m)]); }
		    else {
		      C[i + (j * m)] += (alpha * A[k + (i * m)] * B[k + (j * m)]); }}}}}}
	  else { // Side == R
	    if(Triangular == 'U') {
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < n; k++) {
		    if(k < j) {
		      C[i + (j * m)] += (alpha * B[i + (k * m)] * A[k + (j * n)]); }
		    else {
		      C[i + (j * m)] += (alpha * B[i + (k * m)] * A[j + (k * n)]); }}}}}
	    else { // Triangular == L
	      for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
		  for(k = 0; k < n; k++) { 
		    if(k < j) { 
		      C[i + (j * m)] += (alpha * B[i + (k * m)] * A[j + (k * n)]); }
		    else {
		      C[i + (j * m)] += (alpha * B[i + (k * m)] * A[k + (j * n)]); }}}}}}}}}
  }
  
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::TRMM(char Side, char Triangular, char TransposeA, char UnitDiagonal, int m, int n, ScalarType alpha, ScalarType* A, int LDA, ScalarType* B, int LDB)
  {
    ScalarType zero = ScalarTraits<ScalarType>::zero();
    ScalarType one = ScalarTraits<ScalarType>::one();
    bool BadArgument = 0;
    if(!((Side == 'L') || (Side == 'R'))) {
      cout << "BLAS::GEMM Error: SIDE == " << Side << endl;
      BadArgument = 1; }
    if(!((Triangular == 'L') || (Triangular == 'U'))) {
      cout << "BLAS::GEMM Error: UPLO == " << Triangular << endl;
      BadArgument = 1; }
    if(!((TransposeA == 'N') || (TransposeA == 'T') || (TransposeA == 'C'))) {
      cout << "BLAS::GEMM Error: TRANSA == " << TransposeA << endl;
      BadArgument = 1; }
    if(!((UnitDiagonal == 'N') || (UnitDiagonal == 'U'))) {
      cout << "BLAS::GEMM Error: DIAG == " << UnitDiagonal << endl;
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
	    if(Triangular == 'U') {
	      if(TransposeA == 'N') {
		if(UnitDiagonal == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = i; k < m; k++) {
			B[i + (j * m)] += (A[i + (k * m)] * temp[k + (j * m)]); }}}}
		else { // UnitDiagonal == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = i; k < m; k++) {
			if(i == k) {
			  B[i + (j * m)] += temp[k + (j * m)]; }
			else {
			  B[i + (j * m)] += (A[i + (k * m)] * temp[k + (j * m)]); }}}}}}
	      else { // TransposeA == T or C
		if(UnitDiagonal == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (i + 1); k++) {
			B[i + (j * m)] += (A[k + (i * m)] * temp[k + (j * m)]); }}}}
		else { // UnitDiagonal == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (i + 1); k++) {
			if(i == k) {
			  B[i + (j * m)] += temp[k + (j * m)]; }
			else {
			  B[i + (j * m)] += (A[k + (i * m)] * temp[k + (j * m)]); }}}}}}}
	    else { // Triangular == L
	      if(TransposeA == 'N') {
		if(UnitDiagonal == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (i + 1); k++) { 
			B[i + (j * m)] += (A[i + (k * m)] * temp[k + (j * m)]); }}}}
		else { // UnitDiagonal == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (i + 1); k++) {
			if(i == k) {
			  B[i + (j * m)] += temp[k + (j * m)]; }
			else {
			  B[i + (j * m)] += (A[i + (k * m)] * temp[k + (j * m)]); }}}}}}
	      else { // TransposeA == T or C
		if(UnitDiagonal == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = i; k < m; k++) {		    
			B[i + (j * m)] += (A[k + (i * m)] * temp[k + (j * m)]); }}}}
		else { // UnitDiagonal == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = i; k < m; k++) {
			if(i == k) {
			  B[i + (j * m)] += temp[k + (j * m)]; }
			else {
			  B[i + (j * m)] += (A[k + (i * m)] * temp[k + (j * m)]); }}}}}}}}
	  else { // Side == R
	    if(Triangular == 'U') {
	      if(TransposeA == 'N') {
		if(UnitDiagonal == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (j + 1); k++) {
			B[i + (j * m)] += (temp[i + (k * m)] * A[k + (j * n)]); }}}}
		else { // UnitDiagonal == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (j + 1); k++) {
			if(j == k) {
			  B[i + (j * m)] += temp[i + (k * m)]; }
			else {
			  B[i + (j * m)] += (temp[i + (k * m)] * A[k + (j * n)]); }}}}}}
	      else { // TransposeA == T or C
		if(UnitDiagonal == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = j; k < n; k++) {
			B[i + (j * m)] += (temp[i + (k * m)] * A[j + (k * n)]); }}}}
		else { // UnitDiagonal == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = j; k < n; k++) {
			if(j == k) {
			  B[i + (j * m)] += temp[i + (k * m)]; }
			else {
			  B[i + (j * m)] += (temp[i + (k * m)] * A[j + (k * n)]); }}}}}}}
	    else { // Triangular == L
	      if(TransposeA == 'N') {
		if(UnitDiagonal == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = j; k < n; k++) {
			B[i + (j * m)] += (temp[i + (k * m)] * A[k + (j * n)]); }}}}
		else { // UnitDiagonal == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = j; k < n; k++) {
			if(j == k) {
			  B[i + (j * m)] += temp[i + (k * m)]; }
			else {
			  B[i + (j * m)] += (temp[i + (k * m)] * A[k + (j * n)]); }}}}}}
	      else { // TransposeA == T or C
		if(UnitDiagonal == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (j + 1); k++) {
			B[i + (j * m)] += (temp[i + (k * m)] * A[j + (k * n)]); }}}}
		else { // UnitDiagonal == U
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
	    if(Triangular == 'U') {
	      if(TransposeA == 'N') {
		if(UnitDiagonal == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = i; k < m; k++) {
			B[i + (j * m)] += (alpha * A[i + (k * m)] * temp[k + (j * m)]); }}}}
		else { // UnitDiagonal == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = i; k < m; k++) {
			if(i == k) {
			  B[i + (j * m)] += (alpha * temp[k + (j * m)]); }
			else {
			  B[i + (j * m)] += (alpha * A[i + (k * m)] * temp[k + (j * m)]); }}}}}}
	      else { // TransposeA == T or C
		if(UnitDiagonal == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (i + 1); k++) {
			B[i + (j * m)] += (alpha * A[k + (i * m)] * temp[k + (j * m)]); }}}}
		else { // UnitDiagonal == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (i + 1); k++) {
			if(i == k) {
			  B[i + (j * m)] += (alpha * temp[k + (j * m)]); }
			else {
			  B[i + (j * m)] += (alpha * A[k + (i * m)] * temp[k + (j * m)]); }}}}}}}
	    else { // Triangular == L
	      if(TransposeA == 'N') {
		if(UnitDiagonal == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (i + 1); k++) { 
			B[i + (j * m)] += (alpha * A[i + (k * m)] * temp[k + (j * m)]); }}}}
		else { // UnitDiagonal == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (i + 1); k++) {
			if(i == k) {
			  B[i + (j * m)] += (alpha * temp[k + (j * m)]); }
			else {
			  B[i + (j * m)] += (alpha * A[i + (k * m)] * temp[k + (j * m)]); }}}}}}
	      else { // TransposeA == T or C
		if(UnitDiagonal == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = i; k < m; k++) {		    
			B[i + (j * m)] += (alpha * A[k + (i * m)] * temp[k + (j * m)]); }}}}
		else { // UnitDiagonal == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = i; k < m; k++) {
			if(i == k) {
			  B[i + (j * m)] += (alpha * temp[k + (j * m)]); }
			else {
			  B[i + (j * m)] += (alpha * A[k + (i * m)] * temp[k + (j * m)]); }}}}}}}}
	  else { // Side == R
	    if(Triangular == 'U') {
	      if(TransposeA == 'N') {
		if(UnitDiagonal == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (j + 1); k++) {
			B[i + (j * m)] += (alpha * temp[i + (k * m)] * A[k + (j * n)]); }}}}
		else { // UnitDiagonal == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (j + 1); k++) {
			if(j == k) {
			  B[i + (j * m)] += (alpha * temp[i + (k * m)]); }
			else {
			  B[i + (j * m)] += (alpha * temp[i + (k * m)] * A[k + (j * n)]); }}}}}}
	      else { // TransposeA == T or C
		if(UnitDiagonal == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = j; k < n; k++) {
			B[i + (j * m)] += (alpha * temp[i + (k * m)] * A[j + (k * n)]); }}}}
		else { // UnitDiagonal == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = j; k < n; k++) {
			if(j == k) {
			  B[i + (j * m)] += (alpha * temp[i + (k * m)]); }
			else {
			  B[i + (j * m)] += (alpha * temp[i + (k * m)] * A[j + (k * n)]); }}}}}}}
	    else { // Triangular == L
	      if(TransposeA == 'N') {
		if(UnitDiagonal == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = j; k < n; k++) {
			B[i + (j * m)] += (alpha * temp[i + (k * m)] * A[k + (j * n)]); }}}}
		else { // UnitDiagonal == U
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = j; k < n; k++) {
			if(j == k) {
			  B[i + (j * m)] += (alpha * temp[i + (k * m)]); }
			else {
			  B[i + (j * m)] += (alpha * temp[i + (k * m)] * A[k + (j * n)]); }}}}}}
	      else { // TransposeA == T or C
		if(UnitDiagonal == 'N') {
		  for(i = 0; i < m; i++) {
		    for(j = 0; j < n; j++) {
		      for(k = 0; k < (j + 1); k++) {
			B[i + (j * m)] += (alpha * temp[i + (k * m)] * A[j + (k * n)]); }}}}
		else { // UnitDiagonal == U
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
  void BLAS<OrdinalType, ScalarType>::TRSM(char Side, char Triangular, char TransposeA, char UnitDiagonal, int m, int n, ScalarType alpha, ScalarType* A, int lda, ScalarType* B, int ldb)
  {
    ScalarType zero = ScalarTraits<ScalarType>::zero();
    ScalarType one = ScalarTraits<ScalarType>::one();
    ScalarType temp;
    bool BadArgument = 0;
    if(!((Side == 'L') || (Side == 'R'))) {
      cout << "BLAS::GEMM Error: SIDE == " << Side << endl;
      BadArgument = 1; }
    if(!((Triangular == 'L') || (Triangular == 'U'))) {
      cout << "BLAS::GEMM Error: UPLO == " << Triangular << endl;
      BadArgument = 1; }
    if(!((TransposeA == 'N') || (TransposeA == 'T') || (TransposeA == 'C'))) {
      cout << "BLAS::GEMM Error: TRANSA == " << TransposeA << endl;
      BadArgument = 1; }
    if(!((UnitDiagonal == 'N') || (UnitDiagonal == 'U'))) {
      cout << "BLAS::GEMM Error: DIAG == " << UnitDiagonal << endl;
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
		if(TransposeA == 'N') {
		    //
		    //  Compute B = alpha*inv( A )*B
		    //
		    if(Triangular == 'U') { 
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
				    if (UnitDiagonal == 'N') {
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
                                    if (UnitDiagonal == 'N') {
                                        B[j*ldb + k] /= A[k*lda + k];
                                    }
                                    for(i = k+1; i < m; i++) {
                                        B[j*ldb + i] -= B[j*ldb + k] * A[k*lda + i];
                                    }
                                }
                            }
                        }
		    } // end if (Triangular == 'U')
		}  // if (TransposeA =='N')	
	    	else { 
		    //
		    //  Compute B = alpha*inv( A' )*B
		    //
		    if(Triangular == 'U') { 
			// A is upper triangular.
			for(j = 0; j < n; j++) {
	    	    	    for( i = 0; i < m; i++) {
		    		temp = alpha*B[j*ldb+i];
			    	for(k = 0; k < i; k++) {
				    temp -= A[i*lda + k] * B[j*ldb + k];
				}
				if (UnitDiagonal == 'N') {
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
				if (UnitDiagonal == 'N') {
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
	      if (TransposeA == 'N') {
		    //
		    //  Compute B = alpha*B*inv( A )
		    //
		    if(Triangular == 'U') { 
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
			    if (UnitDiagonal == 'N') {
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
			    if (UnitDiagonal == 'N') {
				temp = one/A[j*lda + j];
				for(i = 0; i < m; i++) {
				    B[j*ldb + i] *= temp;
				}
			    }			
                        }
		    } // end if (Triangular == 'U')
		}  // if (TransposeA =='N')	
	    	else { 
		    //
		    //  Compute B = alpha*B*inv( A' )
		    //
		    if(Triangular == 'U') { 
			// A is upper triangular.
			for(k = (n - 1); k > -1; k--) {
			    if (UnitDiagonal == 'N') {
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
			    if (UnitDiagonal == 'N') {
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
  
  template<typename OrdinalType, typename ScalarType>
  void BLAS<OrdinalType, ScalarType>::XERBLA(char xerbla_arg, int info)
  {
    std::cout << "Warning: default BLAS::XERBLA() doesn't do anything!" << std::endl;
  }

#if 0

  template<typename OrdinalType>
  class BLAS<OrdinalType, float>
  {    
  public:
    inline BLAS(void) {};
    inline BLAS(const BLAS& BLAS_source) {};
    inline virtual ~BLAS(void) {};
    float ASUM(int, float*, int);
    void AXPY(int, float, float*, int, float*, int);
    void COPY(int, float*, int, float*, int);
    float DOT(int, float*, int, float*, int);
    int IAMAX(int, float*, int);   
    float NRM2(int, float*, int);
    void SCAL(int, float, float*, int);
    void GEMV(char, int, int, float, float*, int, float*, int, float, float*, int);
    void GER(int, int, float, float*, int, float*, int, float*, int);
    void TRMV(char, char, char, int, float*, int, float*, int);
    void GEMM(char, char, int, int, int, float, float*, int, float*, int, float, float*, int);
    void SYMM(char, char, int, int, float, float*, int, float*, int, float, float*, int);
    void TRMM(char, char, char, char, int, int, float, float*, int, float*, int);
    void TRSM(char, char, char, char, int, int, float, float*, int, float*, int);
    void XERBLA(char, int);
  };

  template<typename OrdinalType>
  float BLAS<OrdinalType, float>::ASUM(int n, float* x, int incx)
  { return SASUM_F77(&n, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::AXPY(int n, float alpha, float* x, int incx, float* y, int incy)
  { SAXPY_F77(&n, &alpha, x, &incx, y, &incy); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::COPY(int n, float* x, int incx, float* y, int incy)
  { SCOPY_F77(&n, x, &incx, y, &incy); }
  
  template<typename OrdinalType>
  float BLAS<OrdinalType, float>::DOT(int n, float* x, int incx, float* y, int incy)
  { return SDOT_F77(&n, x, &incx, y, &incy); }
  
  template<typename OrdinalType>
  int BLAS<OrdinalType, float>::IAMAX(int n, float* x, int incx)
  { return ISAMAX_F77(&n, x, &incx); }

  template<typename OrdinalType>
  float BLAS<OrdinalType, float>::NRM2(int n, float* x, int incx)
  { return SNRM2_F77(&n, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::SCAL(int n, float alpha, float* x, int incx)
  { SSCAL_F77(&n, &alpha, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::GEMV(char trans, int m, int n, float alpha, float* A, int lda, float* x, int incx, float beta, float* y, int incy)
  { SGEMV_F77(&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::GER(int m, int n, float alpha, float* x, int incx, float* y, int incy, float* A, int lda)
  { SGER_F77(&m, &n, &alpha, x, &incx, y, &incy, A, &lda); }

  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::TRMV(char uplo, char trans, char diag, int n, float* A, int lda, float* x, int incx)
  { STRMV_F77(&uplo, &trans, &diag, &n, A, &lda, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::GEMM(char transa, char transb, int m, int n, int k, float alpha, float* A, int lda, float* B, int ldb, float beta, float* C, int ldc)
  { SGEMM_F77(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::SYMM(char side, char uplo, int m, int n, float alpha, float* A, int lda, float* B, int ldb, float beta, float* C, int ldc)
  { SSYMM_F77(&side, &uplo, &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::TRMM(char side, char uplo, char transa, char diag, int m, int n, float alpha, float* A, int lda, float* B, int ldb)
  { STRMM_F77(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::TRSM(char side, char uplo, char transa, char diag, int m, int n, float alpha, float* A, int lda, float* B, int ldb)
  { STRSM_F77(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb); }

  template<typename OrdinalType>
  void BLAS<OrdinalType, float>::XERBLA(char xerbla_arg, int info)
  { XERBLA_F77(&xerbla_arg, &info); }

#endif

  template<typename OrdinalType>
  class BLAS<OrdinalType, double>
  {    
  public:
    inline BLAS(void) {};
    inline BLAS(const BLAS& BLAS_source) {};
    inline virtual ~BLAS(void) {};
    double ASUM(int, double* x, int);
    void AXPY(int, double, double*, int, double*, int);
    void COPY(int, double*, int, double*, int);
    double DOT(int, double* x, int, double* y, int);
    double NRM2(int, double* x, int);
    void SCAL(int, double, double*, int);
    int IAMAX(int, double*, int);
    void GEMV(char, int, int, double, double*, int, double*, int, double, double*, int);
    void TRMV(char, char, char, int, double*, int, double*, int);
    void GER(int, int, double, double*, int, double*, int, double*, int);
    void GEMM(char, char, int, int, int, double, double*, int, double*, int, double, double*, int);
    void SYMM(char, char, int, int, double, double*, int, double*, int, double, double*, int);
    void TRMM(char, char, char, char, int, int, double, double*, int, double*, int);
    void TRSM(char, char, char, char, int, int, double, double*, int, double*, int);
    void XERBLA(char, int);
  };
  
  template<typename OrdinalType>
  double BLAS<OrdinalType, double>::ASUM(int n, double* x, int incx)
  { return DASUM_F77(&n, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::AXPY(int n, double alpha, double* x, int incx, double* y, int incy)
  { DAXPY_F77(&n, &alpha, x, &incx, y, &incy); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::COPY(int n, double* x, int incx, double* y, int incy)
  { DCOPY_F77(&n, x, &incx, y, &incy); }
  
  template<typename OrdinalType>
  double BLAS<OrdinalType, double>::DOT(int n, double* x, int incx, double* y, int incy)
  { return DDOT_F77(&n, x, &incx, y, &incy); }
  
  template<typename OrdinalType>
  int BLAS<OrdinalType, double>::IAMAX(int n, double* x, int incx)
  { return IDAMAX_F77(&n, x, &incx); }

  template<typename OrdinalType>
  double BLAS<OrdinalType, double>::NRM2(int n, double* x, int incx)
  { return DNRM2_F77(&n, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::SCAL(int n, double alpha, double* x, int incx)
  { DSCAL_F77(&n, &alpha, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::GEMV(char trans, int m, int n, double alpha, double* A, int lda, double* x, int incx, double beta, double* y, int incy)
  { DGEMV_F77(&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::GER(int m, int n, double alpha, double* x, int incx, double* y, int incy, double* A, int lda)
  { DGER_F77(&m, &n, &alpha, x, &incx, y, &incy, A, &lda); }

  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::TRMV(char uplo, char trans, char diag, int n, double* A, int lda, double* x, int incx)
  { DTRMV_F77(&uplo, &trans, &diag, &n, A, &lda, x, &incx); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::GEMM(char transa, char transb, int m, int n, int k, double alpha, double* A, int lda, double* B, int ldb, double beta, double* C, int ldc)
  { DGEMM_F77(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::SYMM(char side, char uplo, int m, int n, double alpha, double* A, int lda, double *B, int ldb, double beta, double *C, int ldc)
  { DSYMM_F77(&side, &uplo, &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::TRMM(char side, char uplo, char transa, char diag, int m, int n, double alpha, double* A, int lda, double* B, int ldb)
  { DTRMM_F77(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb); }

  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::TRSM(char side, char uplo, char transa, char diag, int m, int n, double alpha, double* A, int lda, double* B, int ldb)
  { DTRSM_F77(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb); }
  
  template<typename OrdinalType>
  void BLAS<OrdinalType, double>::XERBLA(char xerbla_arg, int info)
  { XERBLA_F77(&xerbla_arg, &info); }
  
} // namespace Teuchos

#endif // _TEUCHOS_BLAS_HPP_
