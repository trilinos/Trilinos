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

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_BLAS.hpp"

#ifdef HAVE_TEUCHOS_ARPREC
#include "mp/mpreal.h"
#endif

using namespace std;
using namespace Teuchos;

#ifdef HAVE_TEUCHOS_ARPREC
#define SType1	   mp_real
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
  mp::mp_init(20);
#endif

  // SType1 dummy = ScalarTraits<SType1>::zero();
  // SType1 result = Solve(3, dummy);
  // cout << "result = " << result << endl;

  SType1* A = new SType1[3*3];
  SType1* b = new SType1[3];
  SType1* B = new SType1[3*2];

  A[0] = 1;
  A[1] = 0;
  A[2] = 0;
  A[3] = 4;
  A[4] = 2;
  A[5] = 0;
  A[6] = 5;
  A[7] = 6;
  A[8] = 3;

  b[0] = 13;
  b[1] = 20;
  b[2] = 12;

  B[0] = 7;
  B[1] = 7;
  B[2] = 3;
  B[3] = 19;
  B[4] = 16;
  B[5] = 6;

  BLAS<int, SType1> blasObj;

  PrintArrayAsMatrix(A, 3, 3);
  PrintArrayAsVector(b, 3);

  // void  TRSM (ESide side, EUplo uplo, ETransp transa, EDiag diag, const OrdinalType m, const OrdinalType n, const ScalarType alpha, const ScalarType *A, const OrdinalType lda, ScalarType *B, const OrdinalType ldb) const
  
  blasObj.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, Teuchos::NON_UNIT_DIAG, 3, 1, 1, A, 3, b, 3);

  cout << endl;

  PrintArrayAsVector(b, 3);

  cout << " *** " << endl << endl;

  PrintArrayAsMatrix(A, 3, 3);
  PrintArrayAsMatrix(B, 3, 2);
  
  blasObj.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, Teuchos::NON_UNIT_DIAG, 3, 2, 2, A, 3, B, 3);
 
  cout << endl;

  PrintArrayAsMatrix(B, 3, 2);

  delete[] A;
  delete[] b;
  delete[] B;

#ifdef HAVE_TEUCHOS_ARPREC
  mp::mp_finalize();
#endif

  return 0;
}

template<typename TYPE>
void ConstructHilbertMatrix(TYPE* A, int n) {
  TYPE sOne = ScalarTraits<TYPE>::one();
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      A[i + (j * n)] = (sOne / (i + j + sOne));
    }
  }
}

template<typename TYPE>
void ConstructHilbertSumVector(TYPE* x, int n) {
  TYPE sZero = ScalarTraits<TYPE>::zero();
  TYPE sOne = ScalarTraits<TYPE>::one();
  for(int i = 0; i < n; i++) {
    x[i] = sZero;
    for(int j = 0; j < n; j++) {
      x[i] += (sOne / (i + j + sOne));
    }
  }  
}

template<typename TYPE>
bool Cholesky(TYPE* A, int n) {
  TYPE sZero = ScalarTraits<TYPE>::zero();
  for(int k = 0; k < n; k++) {
    for(int j = k + 1; j < n; j++) {
      TYPE alpha = A[k + (j * n)] / A[k + (k * n)];
      for(int i = j; i < n; i++) {
	A[j + (i * n)] -= (alpha * A[k + (i * n)]);
      }
    }
    if(A[k + (k * n)] <= sZero) {
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
  BLAS<int, TYPE> blasObj;
  TYPE* H = new TYPE[n*n];
  TYPE* b = new TYPE[n];

  ConstructHilbertMatrix(H, n);
  ConstructHilbertSumVector(b, n);

  // void COPY (const OrdinalType n, const ScalarType *x, const OrdinalType incx, ScalarType *y, const OrdinalType incy) const -- Copy the vector x to the vector y. 

  PrintArrayAsMatrix(H, n, n);

  TYPE sOne = ScalarTraits<TYPE>::one();

  // void TRSM (ESide side, EUplo uplo, ETransp transa, EDiag diag, const OrdinalType m, const OrdinalType n, const ScalarType alpha, const ScalarType *A, const OrdinalType lda, ScalarType *B, const OrdinalType ldb) const

  PrintArrayAsVector(b, n);

  blasObj.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, Teuchos::NON_UNIT_DIAG, n, n, sOne, H, n, b, n);

  PrintArrayAsVector(b, n);

  delete[] H;
  delete[] b;
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

