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
  
// For a given size n, an n-by-n Hilbert matrix H and a n-by-1 vector b are 
// constructed such that, if Hx* = b, the true solution x* is a one-vector.
// Cholesky factorization is attempted on H; if it fails, no further tests are
// attempted for that datatype. If it is successful, the approximate solution x~
// is computed with a pair of BLAS TRSM (triangular solve) calls. Then, the
// two-norm of (x* - x~) is computed with BLAS AXPY (vector update) and BLAS
// NRM2. The program output is of the form:

//     [size of Hilbert matrix]: [two-norm of (x* - x~)]
  
// Tests for scalar type 2 are performed before scalar type 1 because scalar
// type 2 fails at Cholesky factorization for much lower values of n if the
// mp_real precision is sufficiently large.
  
// Timing analysis still remains to be done for this example, which should be
// easily accomplished with the timing mechanisms native to Teuchos.

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_BLAS.hpp"

#ifdef HAVE_TEUCHOS_ARPREC
#include "mp/mpreal.h"
#endif

using namespace std;
using namespace Teuchos;

#ifdef HAVE_TEUCHOS_ARPREC
#define SType1	   mp_real
#define PRECISION  32
#else
#define SType1     double
#endif
#define SType2     double
#define OType	   int

template<typename TYPE>
void ConstructHilbertMatrix(TYPE*, int);

template<typename TYPE>
void ConstructHilbertSumVector(TYPE*, int);

template<typename TYPE>
bool Cholesky(TYPE*, int);

template<typename TYPE>
TYPE Solve(int, TYPE);

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

#ifdef HAVE_TEUCHOS_ARPREC
  mp::mp_init(PRECISION);
#endif

  cout << "SType2: " << endl;
  int n = 2;
  SType2 dummy2 = ScalarTraits<SType2>::zero();
  SType2 result2 = Solve(2, dummy2);
  while(result2 != 10000.0) {
    cout << n << ": " << result2 << endl;
    n++;
    result2 = Solve(n, dummy2);
  }
  cout << endl << endl << "SType1:" << endl;
  n = 2;
  SType1 dummy1 = ScalarTraits<SType1>::zero();
  SType1 result1 = Solve(2, dummy1);
  while(result1 != 10000.0) {
    cout << n << ": " << result1 << endl;
    n++;
    result1 = Solve(n, dummy1);
  }

#ifdef HAVE_TEUCHOS_ARPREC
  mp::mp_finalize();
#endif

  return 0;
}

template<typename TYPE>
void ConstructHilbertMatrix(TYPE* A, int n) {
  TYPE scal0 = ScalarTraits<TYPE>::zero();
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

template<typename TYPE>
bool Cholesky(TYPE* A, int n) {
  TYPE scal0 = ScalarTraits<TYPE>::zero();
  for(int k = 0; k < n; k++) {
    for(int j = k + 1; j < n; j++) {
      TYPE alpha = A[k + (j * n)] / A[k + (k * n)];
      for(int i = j; i < n; i++) {
	A[j + (i * n)] -= (alpha * A[k + (i * n)]);
      }
    }
    if(A[k + (k * n)] <= scal0) {
      return 0;
    }
    TYPE beta = ScalarTraits<TYPE>::squareroot(A[k + (k * n)]);
    for(int i = k; i < n; i++) {
      A[k + (i * n)] /= beta;
    }
  }
  return 1;
}

template<typename TYPE>
TYPE Solve(int n, TYPE dummy) {
  TYPE scal0 = ScalarTraits<TYPE>::zero();
  TYPE scal1 = ScalarTraits<TYPE>::one();
  TYPE scalNeg1 = scal0 - scal1;
  TYPE result = scal0;
  BLAS<int, TYPE> blasObj;
  TYPE* H = new TYPE[n*n];
  TYPE* b = new TYPE[n];
  TYPE* x = new TYPE[n];
  for(int i = 0; i < n; i++) {
    x[i] = scal1;
  }
  ConstructHilbertMatrix(H, n);
  bool choleskySuccessful = Cholesky(H, n);
  if(!choleskySuccessful) {
    return 10000.0;
  }
  else {
    ConstructHilbertSumVector(b, n);
    blasObj.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, Teuchos::NON_UNIT_DIAG, n, 1, scal1, H, n, b, n);
    blasObj.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, Teuchos::NON_UNIT_DIAG, n, 1, scal1, H, n, b, n);
    blasObj.AXPY(n, scalNeg1, x, 1, b, 1);
    result = blasObj.NRM2(n, b, 1);
    delete[] H;
    delete[] b;
    delete[] x;
    return result;
  }
}

template<typename TYPE>
void PrintArrayAsVector(TYPE* x, int n) {
  cout << "[";
  for(int i = 0; i < n; i++) {
    cout << " " << x[i];
  }
  cout << " ]" << endl;
}

template<typename TYPE>
void PrintArrayAsMatrix(TYPE* a, int m, int n) {
  cout << "[";
  for(int i = 0; i < m; i++) {
    if(i != 0) {
      cout << " ";
    }
    cout << "[";
    for(int j = 0; j < n; j++) {
      cout << " " << a[i + (j * m)];
    }
    cout << " ]";
    if(i != (m - 1)) {
      cout << endl;
    }
  }
  cout << "]" << endl;
}

#ifdef HAVE_TEUCHOS_ARPREC
template<>
void PrintArrayAsVector(mp_real* x, int n) {
  cout << "[ ";
  for(int i = 0; i < n; i++) {
    if(i != 0) {
      cout << "  ";
    }
    cout << x[i];
  }
  cout << "]" << endl;
}

template<>
void PrintArrayAsMatrix(mp_real* a, int m, int n) {
  cout << "[";
  for(int i = 0; i < m; i++) {
    if(i != 0) {
      cout << " ";
    }
    cout << "[";
    for(int j = 0; j < n; j++) {
      if(j != 0) {
	cout << "  ";
      }
      cout << " " << a[i + (j * m)];
    }
    cout << " ]";
    if(i != (m - 1)) {
      cout << endl;
    }
  }
  cout << "]" << endl; 
}
#endif
