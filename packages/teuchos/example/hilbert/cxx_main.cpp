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
#define SType1	   double
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
void PrintArrayAsVector(TYPE*, int);

template<>
void PrintArrayAsVector(mp_real*, int);

template<typename TYPE>
void PrintArrayAsMatrix(TYPE*, int, int);

template<>
void PrintArrayAsMatrix(mp_real*, int, int);

int main(int argc, char *argv[]) {

#ifdef HAVE_TEUCHOS_ARPREC
  mp::mp_init(40);
#endif

  int n = 3;
  SType1* A = new SType1[n*n];
  A[0] = 1;
  A[1] = 2;
  A[2] = 3;
  A[3] = 2;
  A[4] = 5;
  A[5] = 8;
  A[6] = 3;
  A[7] = 8;
  A[8] = 14;

  PrintArrayAsMatrix(A, n, n);

  Cholesky(A, n);

  PrintArrayAsMatrix(A, n, n);

  SType1* H5 = new SType1[5*5];
  ConstructHilbertMatrix(H5, 5);
  PrintArrayAsMatrix(H5, 5, 5);

  SType1* b = new SType1[5];
  ConstructHilbertSumVector(b, 5);
  PrintArrayAsVector(b, 5);

  delete[] A;
  delete[] H5;
  delete[] b;

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
void PrintArrayAsVector(TYPE* x, int n) {
  cout << "[";
  for(int i = 0; i < n; i++) {
    cout << " " << x[i];
  }
  cout << " ]" << endl;
}

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
