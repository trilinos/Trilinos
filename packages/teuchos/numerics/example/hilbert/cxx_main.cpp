// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

// Teuchos Example: Hilbert
  
// This example showcases the usage of BLAS generics with an arbitrary precision
// library -- ARPREC.
  
// Hilbert matrices are classical examples of ill-conditioned matrices. Cholesky
// factorization fails on higher-order Hilbert matrices because they lose their
// positive definiteness when represented with floating-point numbers. We have
// attempted to alleviate this problem with arbitrary precision.
  
// The example program compares two datatypes, scalar type 1 and scalar type 2,
// which can be customized below using #defines. Default types are mp_real (from
// ARPREC) and double. The mp_real datatype must be initialized with a maximum 
// precision value, also customizable below. (Default is 32.)
  
// For a given size n, an n-by-n Hilbert matrix H and a n-by-1 std::vector b are 
// constructed such that, if Hx* = b, the true solution x* is a one-std::vector.
// Cholesky factorization is attempted on H; if it fails, no further tests are
// attempted for that datatype. If it is successful, the approximate solution x~
// is computed with a pair of BLAS TRSM (triangular solve) calls. Then, the
// two-norm of (x* - x~) is computed with BLAS AXPY (std::vector update) and BLAS
// NRM2. The program output is of the form:

//     [size of Hilbert matrix]: [two-norm of (x* - x~)]
  
// Tests for scalar type 2 are performed before scalar type 1 because scalar
// type 2 fails at Cholesky factorization for much lower values of n if the
// mp_real precision is sufficiently large.
  
// Timing analysis still remains to be done for this example, which should be
// easily accomplished with the timing mechanisms native to Teuchos.

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_Version.hpp"
#include <typeinfo>

#ifdef HAVE_TEUCHOS_ARPREC
#include <arprec/mp_real.h>
#endif

#ifdef HAVE_TEUCHOS_GNU_MP
#include "gmp.h"
#include "gmpxx.h"
#endif


using namespace Teuchos;

#ifdef HAVE_TEUCHOS_ARPREC
#define SType1	   mp_real
#elif defined(HAVE_TEUCHOS_GNU_MP)
#define SType1     mpf_class
#else
#define SType1     double
#endif
#define SType2     double
#define OType	   int

template<typename TYPE>
void ConstructHilbertMatrix(TYPE*, int);

template<typename TYPE>
void ConstructHilbertSumVector(TYPE*, int);

template<typename TYPE1, typename TYPE2>
void ConvertHilbertMatrix(TYPE1*, TYPE2*, int);

template<typename TYPE1, typename TYPE2>
void ConvertHilbertSumVector(TYPE1*, TYPE2*, int);

#ifdef HAVE_TEUCHOS_ARPREC
template<>
void ConvertHilbertMatrix(mp_real*, double*, int);

template<>
void ConvertHilbertSumVector(mp_real*, double*, int);
#endif

#ifdef HAVE_TEUCHOS_GNU_MP
template<>
void ConvertHilbertMatrix(mpf_class*, double*, int);

template<>
void ConvertHilbertSumVector(mpf_class*, double*, int);
#endif

template<typename TYPE>
int Cholesky(TYPE*, int);

template<typename TYPE>
int Solve(int, TYPE*, TYPE*, TYPE*);

template<typename TYPE>
void PrintArrayAsVector(TYPE*, int);

template<typename TYPE>
void PrintArrayAsMatrix(TYPE*, int, int);

#ifdef HAVE_TEUCHOS_ARPREC
template<>
void PrintArrayAsVector(mp_real*, int);

template<>
void PrintArrayAsMatrix(mp_real*, int, int);
#endif

int main(int argc, char *argv[]) {

  std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;
  //
  // Create command line processor. 
  //
  Teuchos::CommandLineProcessor hilbertCLP(true, false);
  //
  // Set option for precision and verbosity  
  int precision = 32;
  hilbertCLP.setOption("precision", &precision, "Arbitrary precision");
  bool verbose = false;
  hilbertCLP.setOption("verbose", "quiet", &verbose, "Verbosity of example");
  //
  // Parse command line.
  hilbertCLP.parse( argc, argv );

#ifdef HAVE_TEUCHOS_ARPREC
  mp::mp_init( precision );
#endif

#ifdef HAVE_TEUCHOS_GNU_MP
  mpf_set_default_prec( precision );
  std::cout<< "The precision of the GNU MP variable is (in bits) : "<< mpf_get_default_prec() << std::endl;
#endif
  //
  // Keep track of valid datatypes
  //
  int compSType1 = 1;  // Perform cholesky factorization of matrices of SType1
  int compSType2 = 1;  // Perform cholesky factorization of matrices of SType2
  int convSType1 = 1;  // Perform cholesky factorization of matrices of SType1 (that were converted from SType2)
  int convSType2 = 1;  // Perform cholesky factorization of matrices of SType2 (that were converted from SType1)

  int n = 2;  // Initial dimension of hilbert matrix.
  //
  // Error in solution.
  //
  SType1 result1, result2_1;
  SType2 result2, result1_2;
  //
  // Create pointers to necessary matrices/vectors.
  //
  SType1 *H1=0, *b1=0;
  SType2 *H2=0, *b2=0;
  //
  while ( compSType1>0 || compSType2>0 || convSType1>0 || convSType2>0 ) {
    
    if (compSType1 > 0) {
      H1 = new SType1[ n*n ];
      b1 = new SType1[ n ];
      //
      // Construct problem.
      //
      ConstructHilbertMatrix< SType1 >(H1, n);
      ConstructHilbertSumVector< SType1 >(b1, n);
      //
      // Try to solve it.
      //
      compSType1 = Solve(n, H1, b1, &result1);
      if (compSType1 < 0 && verbose) 
	std::cout << typeid( result1 ).name() << " -- Cholesky factorization failed (negative diagonal) at row "<<-compSType1<< std::endl;
      //
      // Clean up always;
      delete [] H1; H1 = 0;
      delete [] b1; b1 = 0;
    }
    if (compSType2 > 0) {
      H2 = new SType2[ n*n ];
      b2 = new SType2[ n ];
      //
      // Construct problem.
      //
      ConstructHilbertMatrix< SType2 >(H2, n);
      ConstructHilbertSumVector< SType2 >(b2, n);
      //
      // Try to solve it.
      //
      compSType2 = Solve(n, H2, b2, &result2);
      if (compSType2 < 0 && verbose) 
	std::cout << typeid( result2 ).name() << " -- Cholesky factorization failed (negative diagonal) at row "<<-compSType2<< std::endl;
      //
      // Clean up always.
      delete [] H2; H2 = 0;
      delete [] b2; b2 = 0;      
    }
    if (convSType2 > 0) {
      //
      // Create and construct the problem in higher precision
      //
      if (!H1) H1 = new SType1[ n*n ];
      if (!b1) b1 = new SType1[ n ];
      ConstructHilbertMatrix( H1, n );
      ConstructHilbertSumVector( b1, n );
      //
      if (!H2) H2 = new SType2[ n*n ];
      if (!b2) b2 = new SType2[ n ];
      //
      // Convert the problem from SType1 to SType2 ( which should be of lesser precision )
      //
      ConvertHilbertMatrix(H1, H2, n);
      ConvertHilbertSumVector(b1, b2, n);
      //
      // Try to solve it.
      //
      convSType2 = Solve(n, H2, b2, &result1_2);
      if (convSType2 < 0 && verbose) 
	std::cout << typeid( result1_2 ).name() << " (converted) -- Cholesky factorization failed (negative diagonal) at row "<<-convSType2<< std::endl;
      //
      // Clean up
      //
      delete [] H2; H2 = 0;
      delete [] b2; b2 = 0;
      delete [] H1; H1 = 0;
      delete [] b1; b1 = 0;
    }
    if (convSType1 > 0) {
      //
      // Create and construct the problem in lower precision
      //
      if (!H2) H2 = new SType2[ n*n ];
      if (!b2) b2 = new SType2[ n ];
      ConstructHilbertMatrix(H2, n);
      ConstructHilbertSumVector(b2, n);
      //
      if (!H1) H1 = new SType1[ n*n ];
      if (!b1) b1 = new SType1[ n ];
      //
      // Convert the problem from SType2 to SType1 ( which should be of higher precision )
      //
      ConvertHilbertMatrix(H2, H1, n);
      ConvertHilbertSumVector(b2, b1, n);
      //
      // Try to solve it.
      //
      convSType1 = Solve(n, H1, b1, &result2_1);
      if (convSType1 < 0 && verbose) 
	std::cout << typeid( result2_1 ).name() << " (converted) -- Cholesky factorization failed (negative diagonal) at row "<<-convSType1<< std::endl;
      //
      // Clean up
      //
      delete [] H1; H1 = 0;
      delete [] b1; b1 = 0;
      delete [] H2; H2 = 0;
      delete [] b2; b2 = 0;
    }
    if (verbose && (compSType1>0 || compSType2>0 || convSType1>0 || convSType2>0) ) {
      std::cout << "***************************************************" << std::endl;
      std::cout << "Dimension of Hilbert Matrix : "<< n << std::endl;
      std::cout << "***************************************************" << std::endl;
      std::cout << "Datatype : Absolute error || x_hat - x ||"<< std::endl;
      std::cout << "---------------------------------------------------" << std::endl;
    }    
    if (compSType1>0 && verbose)
      std::cout << typeid( result1 ).name() << "\t : "<< result1 << std::endl;
    
    if (convSType1>0 && verbose)
      std::cout << typeid( result2_1 ).name() <<"(converted) : "<< result2_1 << std::endl;

    if (convSType2>0 && verbose)
      std::cout << typeid( result1_2 ).name() <<"(converted) : "<< result2_1 << std::endl;

    if (compSType2>0 && verbose) 
      std::cout << typeid( result2 ).name() << "\t : "<< result2 << std::endl;
    //
    // Increment counter.
    //
    n++;
  }

#ifdef HAVE_TEUCHOS_ARPREC
  mp::mp_finalize();
#endif

  return 0;
}

template<typename TYPE>
void ConstructHilbertMatrix(TYPE* A, int n) {
  TYPE scal1 = ScalarTraits<TYPE>::one();
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      A[i + (j * n)] = (scal1 / (i + j + scal1));
    }
  }
}

template<typename TYPE>
void ConstructHilbertSumVector(TYPE* x, int n) {
  TYPE scal0 = ScalarTraits<TYPE>::zero();
  TYPE scal1 = ScalarTraits<TYPE>::one();
  for(int i = 0; i < n; i++) {
    x[i] = scal0;
    for(int j = 0; j < n; j++) {
      x[i] += (scal1 / (i + j + scal1));
    }
  }  
}

template<typename TYPE1, typename TYPE2>
void ConvertHilbertMatrix(TYPE1* A, TYPE2* B, int n) {
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      B[i + (j * n)] = A[i + (j * n)];
    }
  }
}

template<typename TYPE1, typename TYPE2>
void ConvertHilbertSumVector(TYPE1* x, TYPE2* y, int n) {
  for(int i = 0; i < n; i++) {
    y[i] = x[i];
  }  
}

#ifdef HAVE_TEUCHOS_ARPREC
template<>
void ConvertHilbertMatrix(mp_real* A, double* B, int n) {
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      B[i + (j * n)] = dble( A[i + (j * n)] );
    }
  }
}

template<>
void ConvertHilbertSumVector(mp_real* x, double* y, int n) {
  for(int i = 0; i < n; i++) {
    y[i] = dble( x[i] );
  }  
}
#endif

#ifdef HAVE_TEUCHOS_GNU_MP
template<>
void ConvertHilbertMatrix(mpf_class* A, double* B, int n) {
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      B[i + (j * n)] = A[i + (j * n)].get_d();
    }
  }
}

template<>
void ConvertHilbertSumVector(mpf_class* x, double* y, int n) {
  for(int i = 0; i < n; i++) {
    y[i] = x[i].get_d();
  }  
}
#endif

template<typename TYPE>
int Cholesky(TYPE* A, int n) {
  TYPE scal0 = ScalarTraits<TYPE>::zero();
  for(int k = 0; k < n; k++) {
    for(int j = k + 1; j < n; j++) {
      TYPE alpha = A[k + (j * n)] / A[k + (k * n)];
      for(int i = j; i < n; i++) {
	A[j + (i * n)] -= (alpha * A[k + (i * n)]);
      }
    }
    if(A[k + (k * n)] <= scal0) {
      return -k;
    }
    TYPE beta = ScalarTraits<TYPE>::squareroot(A[k + (k * n)]);
    for(int i = k; i < n; i++) {
      A[k + (i * n)] /= beta;
    }
  }
  return 1;
}

template<typename TYPE>
int Solve(int n, TYPE* H, TYPE* b, TYPE* err) {
  TYPE scal0 = ScalarTraits<TYPE>::zero();
  TYPE scal1 = ScalarTraits<TYPE>::one();
  TYPE scalNeg1 = scal0 - scal1;
  BLAS<int, TYPE> blasObj;
  TYPE* x = new TYPE[n];
  for(int i = 0; i < n; i++) {
    x[i] = scal1;
  }
  int choleskySuccessful = Cholesky(H, n);
  if(choleskySuccessful > 0) {
    blasObj.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, Teuchos::NON_UNIT_DIAG, n, 1, scal1, H, n, b, n);
    blasObj.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, Teuchos::NON_UNIT_DIAG, n, 1, scal1, H, n, b, n);
    blasObj.AXPY(n, scalNeg1, x, 1, b, 1);
    *err = blasObj.NRM2(n, b, 1);
  }  
  delete[] x;
  return choleskySuccessful;
}

template<typename TYPE>
void PrintArrayAsVector(TYPE* x, int n) {
  std::cout << "[";
  for(int i = 0; i < n; i++) {
    std::cout << " " << x[i];
  }
  std::cout << " ]" << std::endl;
}

template<typename TYPE>
void PrintArrayAsMatrix(TYPE* a, int m, int n) {
  std::cout << "[";
  for(int i = 0; i < m; i++) {
    if(i != 0) {
      std::cout << " ";
    }
    std::cout << "[";
    for(int j = 0; j < n; j++) {
      std::cout << " " << a[i + (j * m)];
    }
    std::cout << " ]";
    if(i != (m - 1)) {
      std::cout << std::endl;
    }
  }
  std::cout << "]" << std::endl;
}

#ifdef HAVE_TEUCHOS_ARPREC
template<>
void PrintArrayAsVector(mp_real* x, int n) {
  std::cout << "[ ";
  for(int i = 0; i < n; i++) {
    if(i != 0) {
      std::cout << "  ";
    }
    std::cout << x[i];
  }
  std::cout << "]" << std::endl;
}

template<>
void PrintArrayAsMatrix(mp_real* a, int m, int n) {
  std::cout << "[";
  for(int i = 0; i < m; i++) {
    if(i != 0) {
      std::cout << " ";
    }
    std::cout << "[";
    for(int j = 0; j < n; j++) {
      if(j != 0) {
	std::cout << "  ";
      }
      std::cout << " " << a[i + (j * m)];
    }
    std::cout << " ]";
    if(i != (m - 1)) {
      std::cout << std::endl;
    }
  }
  std::cout << "]" << std::endl; 
}
#endif
