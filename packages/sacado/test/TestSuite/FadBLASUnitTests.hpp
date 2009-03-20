// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef FADBLASUNITTESTS_HPP
#define FADBLASUNITTESTS_HPP

// Sacado includes
#include "Sacado.hpp"
#include "Sacado_Fad_BLAS.hpp"
#include "Sacado_Random.hpp"

// Cppunit includes
#include <cppunit/extensions/HelperMacros.h>

#define COMPARE_VALUES(a, b) \
  CPPUNIT_ASSERT( std::abs(a-b) < this->tol_a + this->tol_r*std::abs(a) );

#define COMPARE_FADS(a, b)                              \
CPPUNIT_ASSERT(a.size() == b.size());			\
CPPUNIT_ASSERT(a.hasFastAccess() == b.hasFastAccess()); \
COMPARE_VALUES(a.val(), b.val());			\
for (int k=0; k<a.size(); k++) {			\
  COMPARE_VALUES(a.dx(k), b.dx(k));			\
  COMPARE_VALUES(a.fastAccessDx(k), b.fastAccessDx(k)); \
 }							\
 ;

#define COMPARE_FAD_VECTORS(X1, X2, n)		\
  CPPUNIT_ASSERT(X1.size() == std::size_t(n));	\
  CPPUNIT_ASSERT(X2.size() == std::size_t(n));	\
  for (unsigned int i=0; i<n; i++) {		\
    COMPARE_FADS(X1[i], X2[i]);			\
  }						\
  ;

// A class for testing differentiated BLAS operations for general Fad types
template <class FadType, class ScalarType>
class FadBLASUnitTests : public CppUnit::TestFixture {

  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;

  CPPUNIT_TEST_SUITE( FadBLASUnitTests );

  CPPUNIT_TEST(testSCAL1);
  CPPUNIT_TEST(testSCAL2);
  CPPUNIT_TEST(testSCAL3);
  CPPUNIT_TEST(testSCAL4);

  CPPUNIT_TEST(testCOPY1);
  CPPUNIT_TEST(testCOPY2);
  CPPUNIT_TEST(testCOPY3);
  CPPUNIT_TEST(testCOPY4);

  CPPUNIT_TEST(testAXPY1);
  CPPUNIT_TEST(testAXPY2);
  CPPUNIT_TEST(testAXPY3);
  CPPUNIT_TEST(testAXPY4);

  CPPUNIT_TEST(testDOT1);
  CPPUNIT_TEST(testDOT2);
  CPPUNIT_TEST(testDOT3);
  CPPUNIT_TEST(testDOT4);

  CPPUNIT_TEST(testNRM21);
  CPPUNIT_TEST(testNRM22);
  
  CPPUNIT_TEST(testGEMV1);
  CPPUNIT_TEST(testGEMV2);
  CPPUNIT_TEST(testGEMV3);
  CPPUNIT_TEST(testGEMV4);
  CPPUNIT_TEST(testGEMV5);
  CPPUNIT_TEST(testGEMV6);
  CPPUNIT_TEST(testGEMV7);
  CPPUNIT_TEST(testGEMV8);
  CPPUNIT_TEST(testGEMV9);

  CPPUNIT_TEST(testTRMV1);
  CPPUNIT_TEST(testTRMV2);
  CPPUNIT_TEST(testTRMV3);
  CPPUNIT_TEST(testTRMV4);

  CPPUNIT_TEST(testGER1);
  CPPUNIT_TEST(testGER2);
  CPPUNIT_TEST(testGER3);
  CPPUNIT_TEST(testGER4);
  CPPUNIT_TEST(testGER5);
  CPPUNIT_TEST(testGER6);
  CPPUNIT_TEST(testGER7);

  CPPUNIT_TEST(testGEMM1);
  CPPUNIT_TEST(testGEMM2);
  CPPUNIT_TEST(testGEMM3);
  CPPUNIT_TEST(testGEMM4);
  CPPUNIT_TEST(testGEMM5);
  CPPUNIT_TEST(testGEMM6);
  CPPUNIT_TEST(testGEMM7);
  CPPUNIT_TEST(testGEMM8);
  CPPUNIT_TEST(testGEMM9);
  CPPUNIT_TEST(testGEMM10);

  CPPUNIT_TEST(testSYMM1);
  CPPUNIT_TEST(testSYMM2);
  CPPUNIT_TEST(testSYMM3);
  CPPUNIT_TEST(testSYMM4);
  CPPUNIT_TEST(testSYMM5);
  CPPUNIT_TEST(testSYMM6);
  CPPUNIT_TEST(testSYMM7);
  CPPUNIT_TEST(testSYMM8);
  CPPUNIT_TEST(testSYMM9);

  CPPUNIT_TEST(testTRMM1);
  CPPUNIT_TEST(testTRMM2);
  CPPUNIT_TEST(testTRMM3);
  CPPUNIT_TEST(testTRMM4);
  CPPUNIT_TEST(testTRMM5);
  CPPUNIT_TEST(testTRMM6);
  CPPUNIT_TEST(testTRMM7);

  CPPUNIT_TEST(testTRSM1);
  CPPUNIT_TEST(testTRSM2);
  CPPUNIT_TEST(testTRSM3);
  CPPUNIT_TEST(testTRSM4);
  CPPUNIT_TEST(testTRSM5);
  CPPUNIT_TEST(testTRSM6);
  CPPUNIT_TEST(testTRSM7);

  CPPUNIT_TEST_SUITE_END();

public:

  FadBLASUnitTests();

  FadBLASUnitTests(int m, int n, int l, int ndot, 
		   double absolute_tolerance, double relative_tolerance);

  void setUp();
  void tearDown();

  void testSCAL1();
  void testSCAL2();
  void testSCAL3();
  void testSCAL4();

  void testCOPY1();
  void testCOPY2();
  void testCOPY3();
  void testCOPY4();

  void testAXPY1();
  void testAXPY2();
  void testAXPY3();
  void testAXPY4();

  void testDOT1();
  void testDOT2();
  void testDOT3();
  void testDOT4();

  void testNRM21();
  void testNRM22();

  void testGEMV1();
  void testGEMV2();
  void testGEMV3();
  void testGEMV4();
  void testGEMV5();
  void testGEMV6();
  void testGEMV7();
  void testGEMV8();
  void testGEMV9();

  void testTRMV1();
  void testTRMV2();
  void testTRMV3();
  void testTRMV4();

  void testGER1();
  void testGER2();
  void testGER3();
  void testGER4();
  void testGER5();
  void testGER6();
  void testGER7();

  void testGEMM1();
  void testGEMM2();
  void testGEMM3();
  void testGEMM4();
  void testGEMM5();
  void testGEMM6();
  void testGEMM7();
  void testGEMM8();
  void testGEMM9();
  void testGEMM10();

  void testSYMM1();
  void testSYMM2();
  void testSYMM3();
  void testSYMM4();
  void testSYMM5();
  void testSYMM6();
  void testSYMM7();
  void testSYMM8();
  void testSYMM9();

  void testTRMM1();
  void testTRMM2();
  void testTRMM3();
  void testTRMM4();
  void testTRMM5();
  void testTRMM6();
  void testTRMM7();

  void testTRSM1();
  void testTRSM2();
  void testTRSM3();
  void testTRSM4();
  void testTRSM5();
  void testTRSM6();
  void testTRSM7();

protected:

  // Random number generator
  Sacado::Random urand;

  // Number of matrix rows
  unsigned int m;

  // Number of matrix columns
  unsigned int n;

  // Number of matrix columns for level 3 blas
  unsigned int l;

  // Number of derivative components
  unsigned int ndot;

  // Tolerances to which fad objects should be the same
  double tol_a, tol_r;

}; // class FadBLASUnitTests

template <class FadType, class ScalarType>
FadBLASUnitTests<FadType,ScalarType>::
FadBLASUnitTests() :
  urand(0.0, 1.0), m(5), n(6), l(4), ndot(7), tol_a(1.0e-12), tol_r(1.0e-12) {}

template <class FadType, class ScalarType>
FadBLASUnitTests<FadType,ScalarType>::
FadBLASUnitTests(int m_, int n_, int l_, int ndot_, double absolute_tolerance, 
		 double relative_tolerance) :
  urand(0.0, 1.0), 
  m(m_),
  n(n_),
  l(l_),
  ndot(ndot_), 
  tol_a(absolute_tolerance), 
  tol_r(relative_tolerance) {}

template <class FadType, class ScalarType>
void FadBLASUnitTests<FadType,ScalarType>::
setUp() {}

template <class FadType, class ScalarType>
void FadBLASUnitTests<FadType,ScalarType>::
tearDown() {}

// Tests all arguments
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testSCAL1() {
  VectorType x1(m,ndot), x2(m,ndot), x3(m,ndot);
  for (unsigned int i=0; i<m; i++) {
    ScalarType val = urand.number();
    x1[i] = FadType(ndot, val);
    x2[i] = FadType(ndot, val);
    x3[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = urand.number();
      x1[i].fastAccessDx(k) = val;
      x2[i].fastAccessDx(k) = val;
      x3[i].fastAccessDx(k) = val;
    }
  }
  FadType alpha(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.SCAL(m, alpha, &x1[0], 1);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.SCAL(m, alpha, &x2[0], 1);

  COMPARE_FAD_VECTORS(x1, x2, m);

  unsigned int sz = m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.SCAL(m, alpha, &x3[0], 1);

  COMPARE_FAD_VECTORS(x1, x3, m);
}

// Tests non-unit inc
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testSCAL2() {
  unsigned int incx = 2;
  VectorType x1(m*incx,ndot), x2(m*incx,ndot), x3(m*incx,ndot);
  for (unsigned int i=0; i<m*incx; i++) {
    ScalarType val = urand.number();
    x1[i] = FadType(ndot, val);
    x2[i] = FadType(ndot, val);
    x3[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = urand.number();
      x1[i].fastAccessDx(k) = val;
      x2[i].fastAccessDx(k) = val;
      x3[i].fastAccessDx(k) = val;
    }
  }
  FadType alpha(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.SCAL(m, alpha, &x1[0], incx);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.SCAL(m, alpha, &x2[0], incx);

  COMPARE_FAD_VECTORS(x1, x2, m*incx);

  unsigned int sz = m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.SCAL(m, alpha, &x3[0], incx);

  COMPARE_FAD_VECTORS(x1, x3, m*incx);
}

// Tests constant alpha
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testSCAL3() {
  VectorType x1(m,ndot), x2(m,ndot), x3(m,ndot);
  for (unsigned int i=0; i<m; i++) {
    ScalarType val = urand.number();
    x1[i] = FadType(ndot, val);
    x2[i] = FadType(ndot, val);
    x3[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = urand.number();
      x1[i].fastAccessDx(k) = val;
      x2[i].fastAccessDx(k) = val;
      x3[i].fastAccessDx(k) = val;
    }
  }
  ScalarType alpha = urand.number();

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.SCAL(m, alpha, &x1[0], 1);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.SCAL(m, alpha, &x2[0], 1);

  COMPARE_FAD_VECTORS(x1, x2, m);

  unsigned int sz = m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.SCAL(m, alpha, &x3[0], 1);

  COMPARE_FAD_VECTORS(x1, x3, m);
}

// Tests constant x
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testSCAL4() {
  VectorType x1(m,ndot), x2(m,ndot), x3(m,ndot);
  for (unsigned int i=0; i<m; i++) {
    ScalarType val = urand.number();
    x1[i] = val;
    x2[i] = val;
    x3[i] = val;
  }
  FadType alpha = FadType(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++)
    alpha.fastAccessDx(k) = urand.number();

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.SCAL(m, alpha, &x1[0], 1);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.SCAL(m, alpha, &x2[0], 1);

  COMPARE_FAD_VECTORS(x1, x2, m);

  unsigned int sz = m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.SCAL(m, alpha, &x3[0], 1);

  COMPARE_FAD_VECTORS(x1, x3, m);
}

// Tests all arguments
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testCOPY1() {
  VectorType x(m,ndot), y1(m,ndot), y2(m,ndot), y3(m,ndot);
  for (unsigned int i=0; i<m; i++) {
    x[i] = FadType(ndot, urand.number());
    ScalarType val = urand.number();
    y1[i] = FadType(ndot, val);
    y2[i] = FadType(ndot, val);
    y3[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      x[i].fastAccessDx(k) = urand.number();
      val = urand.number();
      y1[i].fastAccessDx(k) = val;
      y2[i].fastAccessDx(k) = val;
      y3[i].fastAccessDx(k) = val;
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.COPY(m, &x[0], 1, &y1[0], 1);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.COPY(m, &x[0], 1, &y2[0], 1);

  COMPARE_FAD_VECTORS(y1, y2, m);

  unsigned int sz = 2*m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.COPY(m, &x[0], 1, &y3[0], 1);

  COMPARE_FAD_VECTORS(y1, y3, m);
}

// Tests non unit inc
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testCOPY2() {
  unsigned int incx = 2;
  unsigned int incy = 3;
  VectorType x(m*incx,ndot), y1(m*incy,ndot), y2(m*incy,ndot), y3(m*incy,ndot);
  for (unsigned int i=0; i<m*incx; i++) {
    x[i] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++) {
      x[i].fastAccessDx(k) = urand.number();
    }
  }
  for (unsigned int i=0; i<m*incy; i++) {
    ScalarType val = urand.number();
    y1[i] = FadType(ndot, val);
    y2[i] = FadType(ndot, val);
    y3[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = urand.number();
      y1[i].fastAccessDx(k) = val;
      y2[i].fastAccessDx(k) = val;
      y3[i].fastAccessDx(k) = val;
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.COPY(m, &x[0], incx, &y1[0], incy);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.COPY(m, &x[0], incx, &y2[0], incy);

  COMPARE_FAD_VECTORS(y1, y2, m*incy);

  unsigned int sz = 2*m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.COPY(m, &x[0], incx, &y3[0], incy);

  COMPARE_FAD_VECTORS(y1, y3, m*incy);
}

// Tests x constant
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testCOPY3() {
  VectorType x(m,ndot), y1(m,ndot), y2(m,ndot), y3(m,ndot);
  for (unsigned int i=0; i<m; i++) {
    x[i] = urand.number();
  }
  for (unsigned int i=0; i<m; i++) {
    ScalarType val = urand.number();
    y1[i] = FadType(ndot, val);
    y2[i] = FadType(ndot, val);
    y3[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = urand.number();
      y1[i].fastAccessDx(k) = val;
      y2[i].fastAccessDx(k) = val;
      y3[i].fastAccessDx(k) = val;
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.COPY(m, &x[0], 1, &y1[0], 1);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.COPY(m, &x[0], 1, &y2[0], 1);

  COMPARE_FAD_VECTORS(y1, y2, m);

  unsigned int sz = 2*m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.COPY(m, &x[0], 1, &y3[0], 1);

  COMPARE_FAD_VECTORS(y1, y3, m);
}

// Tests y constant
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testCOPY4() {
  VectorType x(m,ndot), y1(m,ndot), y2(m,ndot), y3(m,ndot);
  for (unsigned int i=0; i<m; i++) {
    x[i] = FadType(ndot, urand.number());
    ScalarType val = urand.number();
    y1[i] = val;
    y2[i] = val;
    y3[i] = val;
    for (unsigned int k=0; k<ndot; k++) {
      x[i].fastAccessDx(k) = urand.number();
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.COPY(m, &x[0], 1, &y1[0], 1);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.COPY(m, &x[0], 1, &y2[0], 1);

  COMPARE_FAD_VECTORS(y1, y2, m);

  unsigned int sz = 2*m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.COPY(m, &x[0], 1, &y3[0], 1);

  COMPARE_FAD_VECTORS(y1, y3, m);
}

// Tests all arguments
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testAXPY1() {
  VectorType x(m,ndot), y1(m,ndot), y2(m,ndot), y3(m,ndot);
  for (unsigned int i=0; i<m; i++) {
    x[i] = FadType(ndot, urand.number());
    ScalarType val = urand.number();
    y1[i] = FadType(ndot, val);
    y2[i] = FadType(ndot, val);
    y3[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      x[i].fastAccessDx(k) = urand.number();
      val = urand.number();
      y1[i].fastAccessDx(k) = val;
      y2[i].fastAccessDx(k) = val;
      y3[i].fastAccessDx(k) = val;
    }
  }
  FadType alpha(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) 
    alpha.fastAccessDx(k) = urand.number();
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.AXPY(m, alpha, &x[0], 1, &y1[0], 1);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.AXPY(m, alpha, &x[0], 1, &y2[0], 1);

  COMPARE_FAD_VECTORS(y1, y2, m);

  unsigned int sz = 2*m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.AXPY(m, alpha, &x[0], 1, &y3[0], 1);

  COMPARE_FAD_VECTORS(y1, y3, m);
}

// Tests non unit inc
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testAXPY2() {
  unsigned int incx = 2;
  unsigned int incy = 3;
  VectorType x(m*incx,ndot), y1(m*incy,ndot), y2(m*incy,ndot), y3(m*incy,ndot);
  for (unsigned int i=0; i<m*incx; i++) {
    x[i] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++) {
      x[i].fastAccessDx(k) = urand.number();
    }
  }
  for (unsigned int i=0; i<m*incy; i++) {
    ScalarType val = urand.number();
    y1[i] = FadType(ndot, val);
    y2[i] = FadType(ndot, val);
    y3[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = urand.number();
      y1[i].fastAccessDx(k) = val;
      y2[i].fastAccessDx(k) = val;
      y3[i].fastAccessDx(k) = val;
    }
  }
  FadType alpha(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) 
    alpha.fastAccessDx(k) = urand.number();
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.AXPY(m, alpha, &x[0], incx, &y1[0], incy);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.AXPY(m, alpha, &x[0], incx, &y2[0], incy);

  COMPARE_FAD_VECTORS(y1, y2, m*incy);

  unsigned int sz = 2*m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.AXPY(m, alpha, &x[0], incx, &y3[0], incy);

  COMPARE_FAD_VECTORS(y1, y3, m*incy);
}

// Tests x constant
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testAXPY3() {
  VectorType x(m,ndot), y1(m,ndot), y2(m,ndot), y3(m,ndot), y4(m,ndot);
  std::vector<ScalarType> xx(m);
  for (unsigned int i=0; i<m; i++) {
    xx[i] = urand.number();
    x[i] = xx[i];
    ScalarType val = urand.number();
    y1[i] = FadType(ndot, val);
    y2[i] = FadType(ndot, val);
    y3[i] = FadType(ndot, val);
    y4[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = urand.number();
      y1[i].fastAccessDx(k) = val;
      y2[i].fastAccessDx(k) = val;
      y3[i].fastAccessDx(k) = val;
      y4[i].fastAccessDx(k) = val;
    }
  }
  FadType alpha(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) 
    alpha.fastAccessDx(k) = urand.number();
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.AXPY(m, alpha, &x[0], 1, &y1[0], 1);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.AXPY(m, alpha, &x[0], 1, &y2[0], 1);

  COMPARE_FAD_VECTORS(y1, y2, m);

  unsigned int sz = m*(1+ndot)+m;
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.AXPY(m, alpha, &x[0], 1, &y3[0], 1);

  COMPARE_FAD_VECTORS(y1, y3, m);

  sacado_blas.AXPY(m, alpha, &xx[0], 1, &y4[0], 1);

  COMPARE_FAD_VECTORS(y1, y4, m);
}

// Tests y constant
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testAXPY4() {
  VectorType x(m,ndot), y1(m,ndot), y2(m,ndot), y3(m,ndot);
  for (unsigned int i=0; i<m; i++) {
    x[i] = FadType(ndot, urand.number());
    ScalarType val = urand.number();
    y1[i] = val;
    y2[i] = val;
    y3[i] = val;
    for (unsigned int k=0; k<ndot; k++) {
      x[i].fastAccessDx(k) = urand.number();
    }
  }
  FadType alpha(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) 
    alpha.fastAccessDx(k) = urand.number();
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.AXPY(m, alpha, &x[0], 1, &y1[0], 1);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.AXPY(m, alpha, &x[0], 1, &y2[0], 1);

  COMPARE_FAD_VECTORS(y1, y2, m);

  unsigned int sz = 2*m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.AXPY(m, alpha, &x[0], 1, &y3[0], 1);

  COMPARE_FAD_VECTORS(y1, y3, m);
}

// Tests all arguments
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testDOT1() {
  VectorType X(m,ndot), Y(m,ndot);
  for (unsigned int i=0; i<m; i++) {
    X[i] = FadType(ndot, urand.number());
    Y[i] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++) {
      X[i].fastAccessDx(k) = urand.number();
      Y[i].fastAccessDx(k) = urand.number();
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  FadType z1 = teuchos_blas.DOT(m, &X[0], 1, &Y[0], 1);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  FadType z2 = sacado_blas.DOT(m, &X[0], 1, &Y[0], 1);

  COMPARE_FADS(z1, z2);

  unsigned int sz = 2*m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  FadType z3 = sacado_blas2.DOT(m, &X[0], 1, &Y[0], 1);

  COMPARE_FADS(z1, z3);
}

// Tests non-unit inc
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testDOT2() {
  unsigned int incx = 2;
  unsigned int incy = 3;
  VectorType X(m*incx,ndot), Y(m*incy,ndot);
  for (unsigned int i=0; i<m*incx; i++) {
    X[i] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++) {
      X[i].fastAccessDx(k) = urand.number();
    }
  }
  for (unsigned int i=0; i<m*incy; i++) {
    Y[i] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++) {
      Y[i].fastAccessDx(k) = urand.number();
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  FadType z1 = teuchos_blas.DOT(m, &X[0], incx, &Y[0], incy);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  FadType z2 = sacado_blas.DOT(m, &X[0], incx, &Y[0], incy);

  COMPARE_FADS(z1, z2);

  unsigned int sz = 2*m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  FadType z3 = sacado_blas2.DOT(m, &X[0], incx, &Y[0], incy);

  COMPARE_FADS(z1, z3);
}

// Tests X constant
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testDOT3() {
  VectorType X(m,0), Y(m,ndot);
  std::vector<ScalarType> x(m);
  for (unsigned int i=0; i<m; i++) {
    x[i] = urand.number();
    X[i] = x[i];
    Y[i] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++) {
      Y[i].fastAccessDx(k) = urand.number();
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  FadType z1 = teuchos_blas.DOT(m, &X[0], 1, &Y[0], 1);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  FadType z2 = sacado_blas.DOT(m, &X[0], 1, &Y[0], 1);

  COMPARE_FADS(z1, z2);

  unsigned int sz = 2*m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  FadType z3 = sacado_blas2.DOT(m, &X[0], 1, &Y[0], 1);

  COMPARE_FADS(z1, z3);

  FadType z4 = sacado_blas.DOT(m, &x[0], 1, &Y[0], 1);

  COMPARE_FADS(z1, z4);
}

// Tests Y constant
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testDOT4() {
  VectorType X(m,ndot), Y(m,0);
  std::vector<ScalarType> y(m);
  for (unsigned int i=0; i<m; i++) {
    X[i] = FadType(ndot, urand.number());
    y[i] = urand.number();
    Y[i] = y[i];
    for (unsigned int k=0; k<ndot; k++) {
      X[i].fastAccessDx(k) = urand.number();
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  FadType z1 = teuchos_blas.DOT(m, &X[0], 1, &Y[0], 1);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  FadType z2 = sacado_blas.DOT(m, &X[0], 1, &Y[0], 1);

  COMPARE_FADS(z1, z2);

  unsigned int sz = 2*m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  FadType z3 = sacado_blas2.DOT(m, &X[0], 1, &Y[0], 1);

  COMPARE_FADS(z1, z3);

  FadType z4 = sacado_blas.DOT(m, &X[0], 1, &y[0], 1);

  COMPARE_FADS(z1, z4);
}

// Tests all arguments
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testNRM21() {
  VectorType X(m,ndot);
  for (unsigned int i=0; i<m; i++) {
    X[i] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++) {
      X[i].fastAccessDx(k) = urand.number();
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  typename Teuchos::ScalarTraits<FadType>::magnitudeType z1 = 
    teuchos_blas.NRM2(m, &X[0], 1);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  typename Teuchos::ScalarTraits<FadType>::magnitudeType z2 = 
    sacado_blas.NRM2(m, &X[0], 1);

  COMPARE_FADS(z1, z2);

  unsigned int sz = m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  typename Teuchos::ScalarTraits<FadType>::magnitudeType z3 = 
    sacado_blas2.NRM2(m, &X[0], 1);

  COMPARE_FADS(z1, z3);
}

// Tests non-unit inc
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testNRM22() {
  unsigned int incx = 2;
  VectorType X(m*incx,ndot);
  for (unsigned int i=0; i<m*incx; i++) {
    X[i] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++) {
      X[i].fastAccessDx(k) = urand.number();
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  typename Teuchos::ScalarTraits<FadType>::magnitudeType z1 = 
    teuchos_blas.NRM2(m, &X[0], incx);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  typename Teuchos::ScalarTraits<FadType>::magnitudeType z2 = 
    sacado_blas.NRM2(m, &X[0], incx);

  COMPARE_FADS(z1, z2);

  unsigned int sz = m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  typename Teuchos::ScalarTraits<FadType>::magnitudeType z3 = 
    sacado_blas2.NRM2(m, &X[0], incx);

  COMPARE_FADS(z1, z3);
}

// Tests all arguments
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testGEMV1() {
  VectorType A(m*n,ndot), B(n,ndot), C1(m,ndot), C2(m,ndot), C3(m,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*m].fastAccessDx(k) = urand.number();
    }
    B[j] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++)
      B[j].fastAccessDx(k) = urand.number();
  }
  FadType alpha(ndot, urand.number());
  FadType beta(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
    beta.fastAccessDx(k) = urand.number();
  }

  for (unsigned int i=0; i<m; i++) {
    ScalarType val = urand.number();
    C1[i] = FadType(ndot, val);
    C2[i] = FadType(ndot, val);
    C3[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = urand.number();
      C1[i].fastAccessDx(k) = val;
      C2[i].fastAccessDx(k) = val;
      C3[i].fastAccessDx(k) = val;
    } 
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1, 
		    beta, &C1[0], 1);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1, 
		   beta, &C2[0], 1);

  COMPARE_FAD_VECTORS(C1, C2, m);

  unsigned int sz = m*n*(1+ndot) + n*(1+ndot) + m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1, 
		    beta, &C3[0], 1);

  COMPARE_FAD_VECTORS(C1, C3, m);
}

// Tests non-unit inc and different lda
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testGEMV2() {
  unsigned int lda = m+3;
  unsigned int incb = 2;
  unsigned int incc = 3;
  VectorType A(lda*n,ndot), B(n*incb,ndot), C1(m*incc,ndot), C2(m*incc,ndot), 
    C3(m*incc,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<lda; i++) {
      A[i+j*lda] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*lda].fastAccessDx(k) = urand.number();
    }
  }
  for (unsigned int j=0; j<n*incb; j++) {
    B[j] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++)
      B[j].fastAccessDx(k) = urand.number();
  }
  FadType alpha(ndot, urand.number());
  FadType beta(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
    beta.fastAccessDx(k) = urand.number();
  }

  for (unsigned int i=0; i<m*incc; i++) {
    ScalarType val = urand.number();
    C1[i] = FadType(ndot, val);
    C2[i] = FadType(ndot, val);
    C3[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = urand.number();
      C1[i].fastAccessDx(k) = val;
      C2[i].fastAccessDx(k) = val;
      C3[i].fastAccessDx(k) = val;
    } 
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], lda, &B[0], incb, 
		    beta, &C1[0], incc);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], lda, &B[0], incb, 
		   beta, &C2[0], incc);

  COMPARE_FAD_VECTORS(C1, C2, m*incc);

  unsigned int sz = m*n*(1+ndot) + n*(1+ndot) + m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], lda, &B[0], incb, 
		    beta, &C3[0], incc);

  COMPARE_FAD_VECTORS(C1, C3, m*incc);
}

// Tests transpose with all arguments
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testGEMV3() {
  VectorType A(m*n,ndot), B(m,ndot), C1(n,ndot), C2(n,ndot), C3(n,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*m].fastAccessDx(k) = urand.number();
    }
  }
  for (unsigned int j=0; j<m; j++) {
    B[j] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++)
      B[j].fastAccessDx(k) = urand.number();
  }
  FadType alpha(ndot, urand.number());
  FadType beta(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
    beta.fastAccessDx(k) = urand.number();
  }

  for (unsigned int i=0; i<n; i++) {
    ScalarType val = urand.number();
    C1[i] = FadType(ndot, val);
    C2[i] = FadType(ndot, val);
    C3[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = urand.number();
      C1[i].fastAccessDx(k) = val;
      C2[i].fastAccessDx(k) = val;
      C3[i].fastAccessDx(k) = val;
    } 
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMV(Teuchos::TRANS, m, n, alpha, &A[0], m, &B[0], 1, 
		    beta, &C1[0], 1);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.GEMV(Teuchos::TRANS, m, n, alpha, &A[0], m, &B[0], 1, 
		   beta, &C2[0], 1);

  COMPARE_FAD_VECTORS(C1, C2, n);

  unsigned int sz = m*n*(1+ndot) + n*(1+ndot) + m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMV(Teuchos::TRANS, m, n, alpha, &A[0], m, &B[0], 1, 
		    beta, &C3[0], 1);

  COMPARE_FAD_VECTORS(C1, C3, n);
}

// Tests transpose with non-unit inc and different lda
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testGEMV4() {
  unsigned int lda = m+3;
  unsigned int incb = 2;
  unsigned int incc = 3;
  VectorType A(lda*n,ndot), B(m*incb,ndot), C1(n*incc,ndot), C2(n*incc,ndot), 
    C3(n*incc,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<lda; i++) {
      A[i+j*lda] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*lda].fastAccessDx(k) = urand.number();
    }
  }
  for (unsigned int j=0; j<m*incb; j++) {
    B[j] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++)
      B[j].fastAccessDx(k) = urand.number();
  }
  FadType alpha(ndot, urand.number());
  FadType beta(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
    beta.fastAccessDx(k) = urand.number();
  }

  for (unsigned int i=0; i<n*incc; i++) {
    ScalarType val = urand.number();
    C1[i] = FadType(ndot, val);
    C2[i] = FadType(ndot, val);
    C3[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = urand.number();
      C1[i].fastAccessDx(k) = val;
      C2[i].fastAccessDx(k) = val;
      C3[i].fastAccessDx(k) = val;
    } 
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMV(Teuchos::TRANS, m, n, alpha, &A[0], lda, &B[0], incb, 
		    beta, &C1[0], incc);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.GEMV(Teuchos::TRANS, m, n, alpha, &A[0], lda, &B[0], incb, 
		   beta, &C2[0], incc);

  COMPARE_FAD_VECTORS(C1, C2, n*incc);

  unsigned int sz = m*n*(1+ndot) + n*(1+ndot) + m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMV(Teuchos::TRANS, m, n, alpha, &A[0], lda, &B[0], incb, 
  		    beta, &C3[0], incc);

  COMPARE_FAD_VECTORS(C1, C3, n*incc);
}

// Tests with constant C
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testGEMV5() {
  VectorType A(m*n,ndot), B(n,ndot), C1(m,ndot), C2(m,ndot), C3(m,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*m].fastAccessDx(k) = urand.number();
    }
    B[j] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++)
      B[j].fastAccessDx(k) = urand.number();
  }
  FadType alpha(ndot, urand.number());
  FadType beta(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
    beta.fastAccessDx(k) = urand.number();
  }

  for (unsigned int i=0; i<m; i++) {
    ScalarType val = urand.number();
    C1[i] = val;
    C2[i] = val;
    C3[i] = val;
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1, 
		    beta, &C1[0], 1);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1, 
		   beta, &C2[0], 1);

  COMPARE_FAD_VECTORS(C1, C2, m);

  unsigned int sz = m*n*(1+ndot) + n*(1+ndot) + m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1, 
		    beta, &C3[0], 1);

  COMPARE_FAD_VECTORS(C1, C3, m);
}

// Tests with constant alpha, beta
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testGEMV6() {
  VectorType A(m*n,ndot), B(n,ndot), C1(m,ndot), C2(m,ndot), C3(m,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*m].fastAccessDx(k) = urand.number();
    }
    B[j] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++)
      B[j].fastAccessDx(k) = urand.number();
  }
  ScalarType alpha = urand.number();
  ScalarType beta = urand.number();

  for (unsigned int i=0; i<m; i++) {
    ScalarType val = urand.number();
    C1[i] = FadType(ndot, val);
    C2[i] = FadType(ndot, val);
    C3[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = urand.number();
      C1[i].fastAccessDx(k) = val;
      C2[i].fastAccessDx(k) = val;
      C3[i].fastAccessDx(k) = val;
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1, 
		    beta, &C1[0], 1);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1, 
		   beta, &C2[0], 1);

  COMPARE_FAD_VECTORS(C1, C2, m);

  unsigned int sz = m*n*(1+ndot) + n*(1+ndot) + m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1, 
		    beta, &C3[0], 1);

  COMPARE_FAD_VECTORS(C1, C3, m);
}

// Tests wth constant B
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testGEMV7() {
  VectorType A(m*n,ndot), B(n,0), C1(m,ndot), C2(m,ndot), C3(m,ndot), 
    C4(m,ndot);
  std::vector<ScalarType> b(n);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*m].fastAccessDx(k) = urand.number();
    }
    b[j] = urand.number();
    B[j] = b[j];
  }
  FadType alpha(ndot, urand.number());
  FadType beta(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
    beta.fastAccessDx(k) = urand.number();
  }

  for (unsigned int i=0; i<m; i++) {
    ScalarType val = urand.number();
    C1[i] = FadType(ndot, val);
    C2[i] = FadType(ndot, val);
    C3[i] = FadType(ndot, val);
    C4[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = urand.number();
      C1[i].fastAccessDx(k) = val;
      C2[i].fastAccessDx(k) = val;
      C3[i].fastAccessDx(k) = val;
      C4[i].fastAccessDx(k) = val;
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1, 
		    beta, &C1[0], 1);
  
  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1, 
		   beta, &C2[0], 1);
  
  COMPARE_FAD_VECTORS(C1, C2, m);
  
  unsigned int sz = m*n*(1+ndot) + n + m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1, 
		    beta, &C3[0], 1);
  
  COMPARE_FAD_VECTORS(C1, C3, m);

  sacado_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &b[0], 1, 
		   beta, &C4[0], 1);
  
  COMPARE_FAD_VECTORS(C1, C4, m);
}

// Tests with constant A
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testGEMV8() {
  VectorType A(m*n,0), B(n,ndot), C1(m,ndot), C2(m,ndot), C3(m,ndot),
    C4(m,ndot);
  std::vector<ScalarType> a(m*n);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      a[i+j*m] = urand.number();
      A[i+j*m] = a[i+j*m];
    }
    B[j] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++)
      B[j].fastAccessDx(k) = urand.number();
  }
  FadType alpha(ndot, urand.number());
  FadType beta(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
    beta.fastAccessDx(k) = urand.number();
  }
  
  for (unsigned int i=0; i<m; i++) {
    ScalarType val = urand.number();
    C1[i] = FadType(ndot, val);
    C2[i] = FadType(ndot, val);
    C3[i] = FadType(ndot, val);
    C4[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = urand.number();
      C1[i].fastAccessDx(k) = val;
      C2[i].fastAccessDx(k) = val;
      C3[i].fastAccessDx(k) = val;
      C4[i].fastAccessDx(k) = val;
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1, 
		    beta, &C1[0], 1);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1, 
		   beta, &C2[0], 1);

  COMPARE_FAD_VECTORS(C1, C2, m);

  unsigned int sz = m*n* + n*(1+ndot) + m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1, 
		    beta, &C3[0], 1);
  
  COMPARE_FAD_VECTORS(C1, C3, m);

  sacado_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &a[0], m, &B[0], 1, 
		   beta, &C4[0], 1);

  COMPARE_FAD_VECTORS(C1, C4, m);
}

// Tests with constant A, B
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testGEMV9() {
  VectorType A(m*n,0), B(n,ndot), C1(m,ndot), C2(m,ndot), C3(m,ndot),
    C4(m,ndot);
  std::vector<ScalarType> a(m*n), b(n);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      a[i+j*m] = urand.number();
      A[i+j*m] = a[i+j*m];
    }
    b[j] = urand.number();
    B[j] = b[j];
  }
  FadType alpha(ndot, urand.number());
  FadType beta(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
    beta.fastAccessDx(k) = urand.number();
  }

  for (unsigned int i=0; i<m; i++) {
    ScalarType val = urand.number();
    C1[i] = FadType(ndot, val);
    C2[i] = FadType(ndot, val);
    C3[i] = FadType(ndot, val);
    C4[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = urand.number();
      C1[i].fastAccessDx(k) = val;
      C2[i].fastAccessDx(k) = val;
      C3[i].fastAccessDx(k) = val;
      C4[i].fastAccessDx(k) = val;
    } 
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1, 
		    beta, &C1[0], 1);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1, 
		   beta, &C2[0], 1);

  COMPARE_FAD_VECTORS(C1, C2, m);

  unsigned int sz = m*n* + n*(1+ndot) + m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1, 
		    beta, &C3[0], 1);
  
  COMPARE_FAD_VECTORS(C1, C3, m);

  sacado_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &a[0], m, &b[0], 1, 
		   beta, &C4[0], 1);

  COMPARE_FAD_VECTORS(C1, C4, m);
}

// Tests all arguments
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testTRMV1() {
  VectorType A(n*n,ndot), x1(n,ndot), x2(n,ndot), x3(n,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<n; i++) {
      A[i+j*n] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*n].fastAccessDx(k) = urand.number();
    }
    ScalarType val = urand.number();
    x1[j] = FadType(ndot, val);
    x2[j] = FadType(ndot, val);
    x3[j] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = urand.number();
      x1[j].fastAccessDx(k) = val;
      x2[j].fastAccessDx(k) = val;
      x3[j].fastAccessDx(k) = val;
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x1[0], 1);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x2[0], 1);

  COMPARE_FAD_VECTORS(x1, x2, n);

  unsigned int sz = n*n*(1+ndot) + n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x3[0], 1);

  COMPARE_FAD_VECTORS(x1, x3, n);

  teuchos_blas.TRMV(Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x1[0], 1);
  sacado_blas.TRMV(Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x2[0], 1);
  sacado_blas2.TRMV(Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x3[0], 1);
  COMPARE_FAD_VECTORS(x1, x2, n);
  COMPARE_FAD_VECTORS(x1, x3, n);

  teuchos_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x1[0], 1);
  sacado_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x2[0], 1);
  sacado_blas2.TRMV(Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x3[0], 1);
  COMPARE_FAD_VECTORS(x1, x2, n);
  COMPARE_FAD_VECTORS(x1, x3, n);

  for (unsigned int i=0; i<n; i++) {
    A[i*n+i].val() = 1.0;
    for (unsigned int k=0; k<ndot; k++)
      A[i*n+i].fastAccessDx(k) = 0.0;
  }
  teuchos_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, n, &A[0], n, &x1[0], 1);
  sacado_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, n, &A[0], n, &x2[0], 1);
  sacado_blas2.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, n, &A[0], n, &x3[0], 1);
  COMPARE_FAD_VECTORS(x1, x2, n);
  COMPARE_FAD_VECTORS(x1, x3, n);
}

// Tests non unit inc, different lda
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testTRMV2() {
  unsigned int lda = n+3;
  unsigned int incx = 2;
  VectorType A(lda*n,ndot), x1(n*incx,ndot), x2(n*incx,ndot), x3(n*incx,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<lda; i++) {
      A[i+j*lda] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*lda].fastAccessDx(k) = urand.number();
    }
  }
  for (unsigned int j=0; j<n*incx; j++) {
    ScalarType val = urand.number();
    x1[j] = FadType(ndot, val);
    x2[j] = FadType(ndot, val);
    x3[j] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = urand.number();
      x1[j].fastAccessDx(k) = val;
      x2[j].fastAccessDx(k) = val;
      x3[j].fastAccessDx(k) = val;
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], lda, &x1[0], incx);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], lda, &x2[0], incx);

  COMPARE_FAD_VECTORS(x1, x2, n*incx);

  unsigned int sz = n*n*(1+ndot) + n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], lda, &x3[0], incx);

  COMPARE_FAD_VECTORS(x1, x3, n*incx);

  teuchos_blas.TRMV(Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], lda, &x1[0], incx);
  sacado_blas.TRMV(Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], lda, &x2[0], incx);
  sacado_blas2.TRMV(Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], lda, &x3[0], incx);
  COMPARE_FAD_VECTORS(x1, x2, n*incx);
  COMPARE_FAD_VECTORS(x1, x3, n*incx);

  teuchos_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], lda, &x1[0], incx);
  sacado_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], lda, &x2[0], incx);
  sacado_blas2.TRMV(Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], lda, &x3[0], incx);
  COMPARE_FAD_VECTORS(x1, x2, n*incx);
  COMPARE_FAD_VECTORS(x1, x3, n*incx);

  for (unsigned int i=0; i<n; i++) {
    A[i*lda+i].val() = 1.0;
    for (unsigned int k=0; k<ndot; k++)
      A[i*lda+i].fastAccessDx(k) = 0.0;
  }
  teuchos_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, n, &A[0], lda, &x1[0], incx);
  sacado_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, n, &A[0], lda, &x2[0], incx);
  sacado_blas2.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, n, &A[0], lda, &x3[0], incx);
  COMPARE_FAD_VECTORS(x1, x2, n*incx);
  COMPARE_FAD_VECTORS(x1, x3, n*incx);
}

// Tests A constant
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testTRMV3() {
  VectorType A(n*n,ndot), x1(n,ndot), x2(n,ndot), x3(n,ndot), x4(n,ndot), 
    x5(n,ndot);
  std::vector<ScalarType> a(n*n);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<n; i++) {
      a[i+j*n] = urand.number();
      A[i+j*n] = a[i+j*n];
    }
    ScalarType val = urand.number();
    x1[j] = FadType(ndot, val);
    x2[j] = FadType(ndot, val);
    x3[j] = FadType(ndot, val);
    x4[j] = FadType(ndot, val);
    x5[j] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = urand.number();
      x1[j].fastAccessDx(k) = val;
      x2[j].fastAccessDx(k) = val;
      x3[j].fastAccessDx(k) = val;
      x4[j].fastAccessDx(k) = val;
      x5[j].fastAccessDx(k) = val;
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x1[0], 1);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x2[0], 1);

  COMPARE_FAD_VECTORS(x1, x2, n);

  unsigned int sz = n*n+n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x3[0], 1);

  COMPARE_FAD_VECTORS(x1, x3, n);

  sacado_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		   Teuchos::NON_UNIT_DIAG, n, &a[0], n, &x4[0], 1);

  COMPARE_FAD_VECTORS(x1, x4, n);

  sacado_blas2.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x5[0], 1);

  COMPARE_FAD_VECTORS(x1, x5, n);

  teuchos_blas.TRMV(Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x1[0], 1);
  sacado_blas.TRMV(Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x2[0], 1);
  sacado_blas2.TRMV(Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x3[0], 1);
  sacado_blas.TRMV(Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &a[0], n, &x4[0], 1);
  sacado_blas2.TRMV(Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &a[0], n, &x5[0], 1);
  COMPARE_FAD_VECTORS(x1, x2, n);
  COMPARE_FAD_VECTORS(x1, x3, n);
  COMPARE_FAD_VECTORS(x1, x4, n);
  COMPARE_FAD_VECTORS(x1, x5, n);

  teuchos_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x1[0], 1);
  sacado_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x2[0], 1);
  sacado_blas2.TRMV(Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x3[0], 1);
  sacado_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &a[0], n, &x4[0], 1);
  sacado_blas2.TRMV(Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &a[0], n, &x5[0], 1);
  COMPARE_FAD_VECTORS(x1, x2, n);
  COMPARE_FAD_VECTORS(x1, x3, n);
  COMPARE_FAD_VECTORS(x1, x4, n);
  COMPARE_FAD_VECTORS(x1, x5, n);

  for (unsigned int i=0; i<n; i++) {
    A[i*n+i].val() = 1.0;
    for (unsigned int k=0; k<ndot; k++)
      A[i*n+i].fastAccessDx(k) = 0.0;
  }
  teuchos_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, n, &A[0], n, &x1[0], 1);
  sacado_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, n, &A[0], n, &x2[0], 1);
  sacado_blas2.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, n, &A[0], n, &x3[0], 1);
  sacado_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, n, &a[0], n, &x4[0], 1);
  sacado_blas2.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, n, &a[0], n, &x5[0], 1);
  COMPARE_FAD_VECTORS(x1, x2, n);
  COMPARE_FAD_VECTORS(x1, x3, n);
  COMPARE_FAD_VECTORS(x1, x4, n);
  COMPARE_FAD_VECTORS(x1, x5, n);
}

// Tests x constant
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testTRMV4() {
  VectorType A(n*n,ndot), x1(n,ndot), x2(n,ndot), x3(n,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<n; i++) {
      A[i+j*n] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*n].fastAccessDx(k) = urand.number();
    }
    ScalarType val = urand.number();
    x1[j] = val;
    x2[j] = val;
    x3[j] = val;
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x1[0], 1);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x2[0], 1);

  COMPARE_FAD_VECTORS(x1, x2, n);

  unsigned int sz = n*n*(1+ndot) + n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x3[0], 1);

  COMPARE_FAD_VECTORS(x1, x3, n);

  teuchos_blas.TRMV(Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x1[0], 1);
  sacado_blas.TRMV(Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x2[0], 1);
  sacado_blas2.TRMV(Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x3[0], 1);
  COMPARE_FAD_VECTORS(x1, x2, n);
  COMPARE_FAD_VECTORS(x1, x3, n);

  teuchos_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x1[0], 1);
  sacado_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x2[0], 1);
  sacado_blas2.TRMV(Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x3[0], 1);
  COMPARE_FAD_VECTORS(x1, x2, n);
  COMPARE_FAD_VECTORS(x1, x3, n);

  for (unsigned int i=0; i<n; i++) {
    A[i*n+i].val() = 1.0;
    for (unsigned int k=0; k<ndot; k++)
      A[i*n+i].fastAccessDx(k) = 0.0;
  }
  teuchos_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, n, &A[0], n, &x1[0], 1);
  sacado_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, n, &A[0], n, &x2[0], 1);
  sacado_blas2.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, n, &A[0], n, &x3[0], 1);
  COMPARE_FAD_VECTORS(x1, x2, n);
  COMPARE_FAD_VECTORS(x1, x3, n);
}

// Tests all arguments
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testGER1() {

  // GER is apparently not defined in the BLAS for complex types
  if (Teuchos::ScalarTraits<ScalarType>::isComplex)
    return; 

  VectorType A1(m*n,ndot), A2(m*n,ndot), A3(m*n,ndot), x(m,ndot), y(n,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = urand.number();
      A1[i+j*m] = FadType(ndot, val);
      A2[i+j*m] = FadType(ndot, val);
      A3[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
      	val = urand.number();
      	A1[i+j*m].fastAccessDx(k) = val;
      	A2[i+j*m].fastAccessDx(k) = val;
      	A3[i+j*m].fastAccessDx(k) = val;
      }
    }
  }
  for (unsigned int i=0; i<m; i++) {
    x[i] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++)
      x[i].fastAccessDx(k) = urand.number();
  }
  for (unsigned int i=0; i<n; i++) {
    y[i] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++)
      y[i].fastAccessDx(k) = urand.number();
  }
  FadType alpha(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A1[0], m);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A2[0], m);

  COMPARE_FAD_VECTORS(A1, A2, m*n);

  unsigned int sz = m*n*(1+ndot) + n*(1+ndot) + m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A3[0], m);

  COMPARE_FAD_VECTORS(A1, A3, m*n);
}

// Tests non unit inc, different lda
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testGER2() {

  // GER is apparently not defined in the BLAS for complex types
  if (Teuchos::ScalarTraits<ScalarType>::isComplex)
    return; 

  unsigned int lda = m+3;
  unsigned int incx = 2;
  unsigned int incy = 3;
  VectorType A1(lda*n,ndot), A2(lda*n,ndot), A3(lda*n,ndot), x(m*incx,ndot), 
    y(n*incy,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<lda; i++) {
      ScalarType val = urand.number();
      A1[i+j*lda] = FadType(ndot, val);
      A2[i+j*lda] = FadType(ndot, val);
      A3[i+j*lda] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
      	val = urand.number();
      	A1[i+j*lda].fastAccessDx(k) = val;
      	A2[i+j*lda].fastAccessDx(k) = val;
      	A3[i+j*lda].fastAccessDx(k) = val;
      }
    }
  }
  for (unsigned int i=0; i<m*incx; i++) {
    x[i] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++)
      x[i].fastAccessDx(k) = urand.number();
  }
  for (unsigned int i=0; i<n*incy; i++) {
    y[i] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++)
      y[i].fastAccessDx(k) = urand.number();
  }
  FadType alpha(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GER(m, n, alpha, &x[0], incx, &y[0], incy, &A1[0], lda);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.GER(m, n, alpha, &x[0], incx, &y[0], incy, &A2[0], lda);

  COMPARE_FAD_VECTORS(A1, A2, lda*n);

  unsigned int sz = m*n*(1+ndot) + n*(1+ndot) + m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GER(m, n, alpha, &x[0], incx, &y[0], incy, &A3[0], lda);

  COMPARE_FAD_VECTORS(A1, A3, lda*n);
}

// Tests constant alpha
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testGER3() {

  // GER is apparently not defined in the BLAS for complex types
  if (Teuchos::ScalarTraits<ScalarType>::isComplex)
    return; 

  VectorType A1(m*n,ndot), A2(m*n,ndot), A3(m*n,ndot), x(m,ndot), y(n,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = urand.number();
      A1[i+j*m] = FadType(ndot, val);
      A2[i+j*m] = FadType(ndot, val);
      A3[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
      	val = urand.number();
      	A1[i+j*m].fastAccessDx(k) = val;
      	A2[i+j*m].fastAccessDx(k) = val;
      	A3[i+j*m].fastAccessDx(k) = val;
      }
    }
  }
  for (unsigned int i=0; i<m; i++) {
    x[i] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++)
      x[i].fastAccessDx(k) = urand.number();
  }
  for (unsigned int i=0; i<n; i++) {
    y[i] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++)
      y[i].fastAccessDx(k) = urand.number();
  }
  ScalarType alpha = urand.number();
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A1[0], m);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A2[0], m);

  COMPARE_FAD_VECTORS(A1, A2, m*n);

  unsigned int sz = m*n*(1+ndot) + n*(1+ndot) + m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A3[0], m);

  COMPARE_FAD_VECTORS(A1, A3, m*n);
}

// Tests constant x
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testGER4() {

  // GER is apparently not defined in the BLAS for complex types
  if (Teuchos::ScalarTraits<ScalarType>::isComplex)
    return; 

  VectorType A1(m*n,ndot), A2(m*n,ndot), A3(m*n,ndot), A4(m*n,ndot), 
    A5(m*n,ndot), x(m,ndot), y(n,ndot);
  std::vector<ScalarType> xx(m);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = urand.number();
      A1[i+j*m] = FadType(ndot, val);
      A2[i+j*m] = FadType(ndot, val);
      A3[i+j*m] = FadType(ndot, val);
      A4[i+j*m] = FadType(ndot, val);
      A5[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
      	val = urand.number();
      	A1[i+j*m].fastAccessDx(k) = val;
      	A2[i+j*m].fastAccessDx(k) = val;
      	A3[i+j*m].fastAccessDx(k) = val;
	A4[i+j*m].fastAccessDx(k) = val;
	A5[i+j*m].fastAccessDx(k) = val;
      }
    }
  }
  for (unsigned int i=0; i<m; i++) {
    xx[i] = urand.number();
    x[i] = xx[i];
  }
  for (unsigned int i=0; i<n; i++) {
    y[i] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++)
      y[i].fastAccessDx(k) = urand.number();
  }
  FadType alpha(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A1[0], m);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A2[0], m);

  COMPARE_FAD_VECTORS(A1, A2, m*n);

  unsigned int sz = m*n*(1+ndot) + n*(1+ndot) + m;
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A3[0], m);

  COMPARE_FAD_VECTORS(A1, A3, m*n);

  sacado_blas.GER(m, n, alpha, &xx[0], 1, &y[0], 1, &A4[0], m);

  COMPARE_FAD_VECTORS(A1, A4, m*n);

  sacado_blas2.GER(m, n, alpha, &xx[0], 1, &y[0], 1, &A5[0], m);

  COMPARE_FAD_VECTORS(A1, A5, m*n);
}

// Tests constant y
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testGER5() {

  // GER is apparently not defined in the BLAS for complex types
  if (Teuchos::ScalarTraits<ScalarType>::isComplex)
    return; 

  VectorType A1(m*n,ndot), A2(m*n,ndot), A3(m*n,ndot), A4(m*n,ndot), 
    A5(m*n,ndot), x(m,ndot), y(n,ndot);
  std::vector<ScalarType> yy(n);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = urand.number();
      A1[i+j*m] = FadType(ndot, val);
      A2[i+j*m] = FadType(ndot, val);
      A3[i+j*m] = FadType(ndot, val);
      A4[i+j*m] = FadType(ndot, val);
      A5[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
      	val = urand.number();
      	A1[i+j*m].fastAccessDx(k) = val;
      	A2[i+j*m].fastAccessDx(k) = val;
      	A3[i+j*m].fastAccessDx(k) = val;
	A4[i+j*m].fastAccessDx(k) = val;
	A5[i+j*m].fastAccessDx(k) = val;
      }
    }
  }
  for (unsigned int i=0; i<m; i++) {
    x[i] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++)
      x[i].fastAccessDx(k) = urand.number();
  }
  for (unsigned int i=0; i<n; i++) {
    yy[i] = urand.number();
    y[i] = yy[i];
  }
  FadType alpha(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A1[0], m);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A2[0], m);

  COMPARE_FAD_VECTORS(A1, A2, m*n);

  unsigned int sz = m*n*(1+ndot) + m*(1+ndot) + n;
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A3[0], m);

  COMPARE_FAD_VECTORS(A1, A3, m*n);

  sacado_blas.GER(m, n, alpha, &x[0], 1, &yy[0], 1, &A4[0], m);

  COMPARE_FAD_VECTORS(A1, A4, m*n);

  sacado_blas2.GER(m, n, alpha, &x[0], 1, &yy[0], 1, &A5[0], m);

  COMPARE_FAD_VECTORS(A1, A5, m*n);
}

// Tests constant x and y
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testGER6() {

  // GER is apparently not defined in the BLAS for complex types
  if (Teuchos::ScalarTraits<ScalarType>::isComplex)
    return; 

  VectorType A1(m*n,ndot), A2(m*n,ndot), A3(m*n,ndot), A4(m*n,ndot), 
    A5(m*n,ndot), x(m,ndot), y(n,ndot);
  std::vector<ScalarType> xx(n), yy(n);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = urand.number();
      A1[i+j*m] = FadType(ndot, val);
      A2[i+j*m] = FadType(ndot, val);
      A3[i+j*m] = FadType(ndot, val);
      A4[i+j*m] = FadType(ndot, val);
      A5[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
      	val = urand.number();
      	A1[i+j*m].fastAccessDx(k) = val;
      	A2[i+j*m].fastAccessDx(k) = val;
      	A3[i+j*m].fastAccessDx(k) = val;
	A4[i+j*m].fastAccessDx(k) = val;
	A5[i+j*m].fastAccessDx(k) = val;
      }
    }
  }
  for (unsigned int i=0; i<m; i++) {
    xx[i] = urand.number();
    x[i] = xx[i];
  }
  for (unsigned int i=0; i<n; i++) {
    yy[i] = urand.number();
    y[i] = yy[i];
  }
  FadType alpha(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A1[0], m);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A2[0], m);

  COMPARE_FAD_VECTORS(A1, A2, m*n);

  unsigned int sz = m*n*(1+ndot) + m + n;
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A3[0], m);

  COMPARE_FAD_VECTORS(A1, A3, m*n);

  sacado_blas.GER(m, n, alpha, &xx[0], 1, &yy[0], 1, &A4[0], m);

  COMPARE_FAD_VECTORS(A1, A4, m*n);

  sacado_blas2.GER(m, n, alpha, &xx[0], 1, &yy[0], 1, &A5[0], m);

  COMPARE_FAD_VECTORS(A1, A5, m*n);
}

// Tests constant A
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testGER7() {

  // GER is apparently not defined in the BLAS for complex types
  if (Teuchos::ScalarTraits<ScalarType>::isComplex)
    return; 

  VectorType A1(m*n,ndot), A2(m*n,ndot), A3(m*n,ndot), x(m,ndot), y(n,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = urand.number();
      A1[i+j*m] = val;
      A2[i+j*m] = val;
      A3[i+j*m] = val;
    }
  }
  for (unsigned int i=0; i<m; i++) {
    x[i] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++)
      x[i].fastAccessDx(k) = urand.number();
  }
  for (unsigned int i=0; i<n; i++) {
    y[i] = FadType(ndot, urand.number());
    for (unsigned int k=0; k<ndot; k++)
      y[i].fastAccessDx(k) = urand.number();
  }
  FadType alpha(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A1[0], m);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A2[0], m);

  COMPARE_FAD_VECTORS(A1, A2, m*n);

  unsigned int sz = m*n*(1+ndot) + n*(1+ndot) + m*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A3[0], m);

  COMPARE_FAD_VECTORS(A1, A3, m*n);
}

// Tests all arguments
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testGEMM1() {
  VectorType A(m*l,ndot), B(l*n,ndot), C1(m*n,ndot), C2(m*n,ndot), C3(m*n,ndot);
  for (unsigned int j=0; j<l; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*m].fastAccessDx(k) = urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<l; i++) {
      B[i+j*l] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
	B[i+j*l].fastAccessDx(k) = urand.number();
    }
  }
  FadType alpha(ndot, urand.number());
  FadType beta(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
    beta.fastAccessDx(k) = urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = urand.number();
      C1[i+j*m] = FadType(ndot, val);
      C2[i+j*m] = FadType(ndot, val);
      C3[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
	val = urand.number();
	C1[i+j*m].fastAccessDx(k) = val;
	C2[i+j*m].fastAccessDx(k) = val;
	C3[i+j*m].fastAccessDx(k) = val;
      } 
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], l, beta, &C1[0], m);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], l, beta, &C2[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);

  unsigned int sz = m*l*(1+ndot) + l*n*(1+ndot) + m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], l, beta, &C3[0], m);

  COMPARE_FAD_VECTORS(C1, C3, m*n);

  // transa
  teuchos_blas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], l, beta, &C1[0], m);
  sacado_blas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], l, beta, &C2[0], m);
  sacado_blas2.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], l, beta, &C3[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);
  COMPARE_FAD_VECTORS(C1, C3, m*n);

  // transb
  teuchos_blas.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], n, beta, &C1[0], m);
  sacado_blas.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], n, beta, &C2[0], m);
  sacado_blas2.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], n, beta, &C3[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);
  COMPARE_FAD_VECTORS(C1, C3, m*n);

  // transa and transb
  teuchos_blas.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], n, beta, &C1[0], m);
  sacado_blas.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], n, beta, &C2[0], m);
  sacado_blas2.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], n, beta, &C3[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);
  COMPARE_FAD_VECTORS(C1, C3, m*n);
}

// Tests different lda, ldb, ldc
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testGEMM2() {
  unsigned int lda = m+4;
  unsigned int ldb = l+4;
  unsigned int ldc = m+5;
  VectorType A(lda*l,ndot), B(ldb*n,ndot), C1(ldc*n,ndot), C2(ldc*n,ndot), 
    C3(ldc*n,ndot);
  for (unsigned int j=0; j<l; j++) {
    for (unsigned int i=0; i<lda; i++) {
      A[i+j*lda] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*lda].fastAccessDx(k) = urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<ldb; i++) {
      B[i+j*ldb] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
	B[i+j*ldb].fastAccessDx(k) = urand.number();
    }
  }
  FadType alpha(ndot, urand.number());
  FadType beta(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
    beta.fastAccessDx(k) = urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<ldc; i++) {
      ScalarType val = urand.number();
      C1[i+j*ldc] = FadType(ndot, val);
      C2[i+j*ldc] = FadType(ndot, val);
      C3[i+j*ldc] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
	val = urand.number();
	C1[i+j*ldc].fastAccessDx(k) = val;
	C2[i+j*ldc].fastAccessDx(k) = val;
	C3[i+j*ldc].fastAccessDx(k) = val;
      } 
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], lda, &B[0], ldb, beta, &C1[0], ldc);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], lda, &B[0], ldb, beta, &C2[0], ldc);

  COMPARE_FAD_VECTORS(C1, C2, ldc*n);

  unsigned int sz = m*l*(1+ndot) + l*n*(1+ndot) + m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], lda, &B[0], ldb, beta, &C3[0], ldc);

  COMPARE_FAD_VECTORS(C1, C3, ldc*n);
}

// Tests different lda, ldb, ldc with transa
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testGEMM3() {
  unsigned int lda = l+3;
  unsigned int ldb = l+4;
  unsigned int ldc = m+5;
  VectorType A(lda*m,ndot), B(ldb*n,ndot), C1(ldc*n,ndot), C2(ldc*n,ndot), 
    C3(ldc*n,ndot);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<lda; i++) {
      A[i+j*lda] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*lda].fastAccessDx(k) = urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<ldb; i++) {
      B[i+j*ldb] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
	B[i+j*ldb].fastAccessDx(k) = urand.number();
    }
  }
  FadType alpha(ndot, urand.number());
  FadType beta(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
    beta.fastAccessDx(k) = urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<ldc; i++) {
      ScalarType val = urand.number();
      C1[i+j*ldc] = FadType(ndot, val);
      C2[i+j*ldc] = FadType(ndot, val);
      C3[i+j*ldc] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
	val = urand.number();
	C1[i+j*ldc].fastAccessDx(k) = val;
	C2[i+j*ldc].fastAccessDx(k) = val;
	C3[i+j*ldc].fastAccessDx(k) = val;
      } 
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], lda, &B[0], ldb, beta, &C1[0], ldc);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], lda, &B[0], ldb, beta, &C2[0], ldc);

  COMPARE_FAD_VECTORS(C1, C2, ldc*n);

  unsigned int sz = m*l*(1+ndot) + l*n*(1+ndot) + m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], lda, &B[0], ldb, beta, &C3[0], ldc);

  COMPARE_FAD_VECTORS(C1, C3, ldc*n);
}

// Tests different lda, ldb, ldc with transb
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testGEMM4() {
  unsigned int lda = m+4;
  unsigned int ldb = n+4;
  unsigned int ldc = m+5;
  VectorType A(lda*l,ndot), B(ldb*l,ndot), C1(ldc*n,ndot), C2(ldc*n,ndot), 
    C3(ldc*n,ndot);
  for (unsigned int j=0; j<l; j++) {
    for (unsigned int i=0; i<lda; i++) {
      A[i+j*lda] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*lda].fastAccessDx(k) = urand.number();
    }
  }
  for (unsigned int j=0; j<l; j++) {
    for (unsigned int i=0; i<ldb; i++) {
      B[i+j*ldb] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
	B[i+j*ldb].fastAccessDx(k) = urand.number();
    }
  }
  FadType alpha(ndot, urand.number());
  FadType beta(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
    beta.fastAccessDx(k) = urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<ldc; i++) {
      ScalarType val = urand.number();
      C1[i+j*ldc] = FadType(ndot, val);
      C2[i+j*ldc] = FadType(ndot, val);
      C3[i+j*ldc] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
	val = urand.number();
	C1[i+j*ldc].fastAccessDx(k) = val;
	C2[i+j*ldc].fastAccessDx(k) = val;
	C3[i+j*ldc].fastAccessDx(k) = val;
      } 
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], lda, &B[0], ldb, beta, &C1[0], ldc);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], lda, &B[0], ldb, beta, &C2[0], ldc);

  COMPARE_FAD_VECTORS(C1, C2, ldc*n);

  unsigned int sz = m*l*(1+ndot) + l*n*(1+ndot) + m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], lda, &B[0], ldb, beta, &C3[0], ldc);

  COMPARE_FAD_VECTORS(C1, C3, ldc*n);
}

// Tests different lda, ldb, ldc with transa and transb
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testGEMM5() {
  unsigned int lda = l+3;
  unsigned int ldb = n+4;
  unsigned int ldc = m+5;
  VectorType A(lda*m,ndot), B(ldb*l,ndot), C1(ldc*n,ndot), C2(ldc*n,ndot), 
    C3(ldc*n,ndot);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<lda; i++) {
      A[i+j*lda] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*lda].fastAccessDx(k) = urand.number();
    }
  }
  for (unsigned int j=0; j<l; j++) {
    for (unsigned int i=0; i<ldb; i++) {
      B[i+j*ldb] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
	B[i+j*ldb].fastAccessDx(k) = urand.number();
    }
  }
  FadType alpha(ndot, urand.number());
  FadType beta(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
    beta.fastAccessDx(k) = urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<ldc; i++) {
      ScalarType val = urand.number();
      C1[i+j*ldc] = FadType(ndot, val);
      C2[i+j*ldc] = FadType(ndot, val);
      C3[i+j*ldc] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
	val = urand.number();
	C1[i+j*ldc].fastAccessDx(k) = val;
	C2[i+j*ldc].fastAccessDx(k) = val;
	C3[i+j*ldc].fastAccessDx(k) = val;
      } 
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], lda, &B[0], ldb, beta, &C1[0], ldc);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], lda, &B[0], ldb, beta, &C2[0], ldc);

  COMPARE_FAD_VECTORS(C1, C2, ldc*n);

  unsigned int sz = m*l*(1+ndot) + l*n*(1+ndot) + m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], lda, &B[0], ldb, beta, &C3[0], ldc);

  COMPARE_FAD_VECTORS(C1, C3, ldc*n);
}

// Tests with constant C
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testGEMM6() {
  VectorType A(m*l,ndot), B(l*n,ndot), C1(m*n,ndot), C2(m*n,ndot), C3(m*n,ndot);
  for (unsigned int j=0; j<l; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*m].fastAccessDx(k) = urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<l; i++) {
      B[i+j*l] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
	B[i+j*l].fastAccessDx(k) = urand.number();
    }
  }
  FadType alpha(ndot, urand.number());
  FadType beta(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
    beta.fastAccessDx(k) = urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = urand.number();
      C1[i+j*m] = val;
      C2[i+j*m] = val;
      C3[i+j*m] = val;
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], l, beta, &C1[0], m);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], l, beta, &C2[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);

  unsigned int sz = m*l*(1+ndot) + l*n*(1+ndot) + m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], l, beta, &C3[0], m);

  COMPARE_FAD_VECTORS(C1, C3, m*n);

  // transa
  teuchos_blas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], l, beta, &C1[0], m);
  sacado_blas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], l, beta, &C2[0], m);
  sacado_blas2.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], l, beta, &C3[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);
  COMPARE_FAD_VECTORS(C1, C3, m*n);

  // transb
  teuchos_blas.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], n, beta, &C1[0], m);
  sacado_blas.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], n, beta, &C2[0], m);
  sacado_blas2.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], n, beta, &C3[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);
  COMPARE_FAD_VECTORS(C1, C3, m*n);

  // transa and transb
  teuchos_blas.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], n, beta, &C1[0], m);
  sacado_blas.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], n, beta, &C2[0], m);
  sacado_blas2.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], n, beta, &C3[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);
  COMPARE_FAD_VECTORS(C1, C3, m*n);
}

// Tests with constant alpha, beta
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testGEMM7() {
  VectorType A(m*l,ndot), B(l*n,ndot), C1(m*n,ndot), C2(m*n,ndot), C3(m*n,ndot);
  for (unsigned int j=0; j<l; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*m].fastAccessDx(k) = urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<l; i++) {
      B[i+j*l] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
	B[i+j*l].fastAccessDx(k) = urand.number();
    }
  }
  ScalarType alpha = urand.number();
  ScalarType beta = urand.number();

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = urand.number();
      C1[i+j*m] = FadType(ndot, val);
      C2[i+j*m] = FadType(ndot, val);
      C3[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
	val = urand.number();
	C1[i+j*m].fastAccessDx(k) = val;
	C2[i+j*m].fastAccessDx(k) = val;
	C3[i+j*m].fastAccessDx(k) = val;
      }
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], l, beta, &C1[0], m);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], l, beta, &C2[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);

  unsigned int sz = m*l*(1+ndot) + l*n*(1+ndot) + m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], l, beta, &C3[0], m);

  COMPARE_FAD_VECTORS(C1, C3, m*n);

  // transa
  teuchos_blas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], l, beta, &C1[0], m);
  sacado_blas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], l, beta, &C2[0], m);
  sacado_blas2.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], l, beta, &C3[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);
  COMPARE_FAD_VECTORS(C1, C3, m*n);

  // transb
  teuchos_blas.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], n, beta, &C1[0], m);
  sacado_blas.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], n, beta, &C2[0], m);
  sacado_blas2.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], n, beta, &C3[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);
  COMPARE_FAD_VECTORS(C1, C3, m*n);

  // transa and transb
  teuchos_blas.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], n, beta, &C1[0], m);
  sacado_blas.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], n, beta, &C2[0], m);
  sacado_blas2.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], n, beta, &C3[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);
  COMPARE_FAD_VECTORS(C1, C3, m*n);
}

// Tests with constant A
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testGEMM8() {
  VectorType A(m*l,ndot), B(l*n,ndot), C1(m*n,ndot), C2(m*n,ndot), C3(m*n,ndot),
    C4(m*n,ndot), C5(m*n,ndot);
  std::vector<ScalarType> a(m*l);
  for (unsigned int j=0; j<l; j++) {
    for (unsigned int i=0; i<m; i++) {
      a[i+j*m] = urand.number();
      A[i+j*m] = a[i+j*m];
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<l; i++) {
      B[i+j*l] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
	B[i+j*l].fastAccessDx(k) = urand.number();
    }
  }
  FadType alpha(ndot, urand.number());
  FadType beta(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
    beta.fastAccessDx(k) = urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = urand.number();
      C1[i+j*m] = FadType(ndot, val);
      C2[i+j*m] = FadType(ndot, val);
      C3[i+j*m] = FadType(ndot, val);
      C4[i+j*m] = FadType(ndot, val);
      C5[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
	val = urand.number();
	C1[i+j*m].fastAccessDx(k) = val;
	C2[i+j*m].fastAccessDx(k) = val;
	C3[i+j*m].fastAccessDx(k) = val;
	C4[i+j*m].fastAccessDx(k) = val;
	C5[i+j*m].fastAccessDx(k) = val;
      }
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], l, beta, &C1[0], m);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], l, beta, &C2[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);

  unsigned int sz = m*l + l*n*(1+ndot) + m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], l, beta, &C3[0], m);

  COMPARE_FAD_VECTORS(C1, C3, m*n);

  sacado_blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &a[0], m, &B[0], l, beta, &C4[0], m);

  COMPARE_FAD_VECTORS(C1, C4, m*n);

  sacado_blas2.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &a[0], m, &B[0], l, beta, &C5[0], m);

  COMPARE_FAD_VECTORS(C1, C5, m*n);

  // transa
  teuchos_blas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], l, beta, &C1[0], m);
  sacado_blas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], l, beta, &C2[0], m);
  sacado_blas2.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], l, beta, &C3[0], m);
  sacado_blas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &a[0], l, &B[0], l, beta, &C4[0], m);
  sacado_blas2.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &a[0], l, &B[0], l, beta, &C5[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);
  COMPARE_FAD_VECTORS(C1, C3, m*n);
  COMPARE_FAD_VECTORS(C1, C4, m*n);
  COMPARE_FAD_VECTORS(C1, C5, m*n);

  // transb
  teuchos_blas.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], n, beta, &C1[0], m);
  sacado_blas.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], n, beta, &C2[0], m);
  sacado_blas2.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], n, beta, &C3[0], m);
  sacado_blas.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &a[0], m, &B[0], n, beta, &C4[0], m);
  sacado_blas2.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &a[0], m, &B[0], n, beta, &C5[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);
  COMPARE_FAD_VECTORS(C1, C3, m*n);
  COMPARE_FAD_VECTORS(C1, C4, m*n);
  COMPARE_FAD_VECTORS(C1, C5, m*n);

  // transa and transb
  teuchos_blas.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], n, beta, &C1[0], m);
  sacado_blas.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], n, beta, &C2[0], m);
  sacado_blas2.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], n, beta, &C3[0], m);
  sacado_blas.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &a[0], l, &B[0], n, beta, &C4[0], m);
  sacado_blas2.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &a[0], l, &B[0], n, beta, &C5[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);
  COMPARE_FAD_VECTORS(C1, C3, m*n);
  COMPARE_FAD_VECTORS(C1, C4, m*n);
  COMPARE_FAD_VECTORS(C1, C5, m*n);
}

// Tests with constant B
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testGEMM9() {
  VectorType A(m*l,ndot), B(l*n,ndot), C1(m*n,ndot), C2(m*n,ndot), C3(m*n,ndot),
    C4(m*n,ndot), C5(m*n,ndot);
  std::vector<ScalarType> b(l*n);
  for (unsigned int j=0; j<l; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*m].fastAccessDx(k) = urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<l; i++) {
      b[i+j*l] = urand.number();
      B[i+j*l] = b[i+j*l];
    }
  }
  FadType alpha(ndot, urand.number());
  FadType beta(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
    beta.fastAccessDx(k) = urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = urand.number();
      C1[i+j*m] = FadType(ndot, val);
      C2[i+j*m] = FadType(ndot, val);
      C3[i+j*m] = FadType(ndot, val);
      C4[i+j*m] = FadType(ndot, val);
      C5[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
    	val = urand.number();
    	C1[i+j*m].fastAccessDx(k) = val;
    	C2[i+j*m].fastAccessDx(k) = val;
    	C3[i+j*m].fastAccessDx(k) = val;
    	C4[i+j*m].fastAccessDx(k) = val;
    	C5[i+j*m].fastAccessDx(k) = val;
      } 
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], l, beta, &C1[0], m);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], l, beta, &C2[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);

  unsigned int sz = m*l*(1+ndot) + l*n + m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], l, beta, &C3[0], m);

  COMPARE_FAD_VECTORS(C1, C3, m*n);

  sacado_blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], m, &b[0], l, beta, &C4[0], m);

  COMPARE_FAD_VECTORS(C1, C4, m*n);

  sacado_blas2.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], m, &b[0], l, beta, &C5[0], m);

  COMPARE_FAD_VECTORS(C1, C5, m*n);

  // transa
  teuchos_blas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], l, beta, &C1[0], m);
  sacado_blas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], l, beta, &C2[0], m);
  sacado_blas2.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], l, beta, &C3[0], m);
  sacado_blas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], l, &b[0], l, beta, &C4[0], m);
  sacado_blas2.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], l, &b[0], l, beta, &C5[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);
  COMPARE_FAD_VECTORS(C1, C3, m*n);
  COMPARE_FAD_VECTORS(C1, C4, m*n);
  COMPARE_FAD_VECTORS(C1, C5, m*n);

  // transb
  teuchos_blas.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], n, beta, &C1[0], m);
  sacado_blas.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], n, beta, &C2[0], m);
  sacado_blas2.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], n, beta, &C3[0], m);
  sacado_blas.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], m, &b[0], n, beta, &C4[0], m);
  sacado_blas2.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], m, &b[0], n, beta, &C5[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);
  COMPARE_FAD_VECTORS(C1, C3, m*n);
  COMPARE_FAD_VECTORS(C1, C4, m*n);
  COMPARE_FAD_VECTORS(C1, C5, m*n);

  // transa and transb
  teuchos_blas.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], n, beta, &C1[0], m);
  sacado_blas.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], n, beta, &C2[0], m);
  sacado_blas2.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], n, beta, &C3[0], m);
  sacado_blas.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], l, &b[0], n, beta, &C4[0], m);
  sacado_blas2.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], l, &b[0], n, beta, &C5[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);
  COMPARE_FAD_VECTORS(C1, C3, m*n);
  COMPARE_FAD_VECTORS(C1, C4, m*n);
  COMPARE_FAD_VECTORS(C1, C5, m*n);
}

// Tests with constant A and B
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testGEMM10() {
  VectorType A(m*l,ndot), B(l*n,ndot), C1(m*n,ndot), C2(m*n,ndot), C3(m*n,ndot),
    C4(m*n,ndot), C5(m*n,ndot);
  std::vector<ScalarType> a(m*l), b(l*n);
  for (unsigned int j=0; j<l; j++) {
    for (unsigned int i=0; i<m; i++) {
      a[i+j*m] = urand.number();
      A[i+j*m] = a[i+j*m];
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<l; i++) {
      b[i+j*l] = urand.number();
      B[i+j*l] = b[i+j*l];
    }
  }
  FadType alpha(ndot, urand.number());
  FadType beta(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
    beta.fastAccessDx(k) = urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = urand.number();
      C1[i+j*m] = FadType(ndot, val);
      C2[i+j*m] = FadType(ndot, val);
      C3[i+j*m] = FadType(ndot, val);
      C4[i+j*m] = FadType(ndot, val);
      C5[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
    	val = urand.number();
    	C1[i+j*m].fastAccessDx(k) = val;
    	C2[i+j*m].fastAccessDx(k) = val;
    	C3[i+j*m].fastAccessDx(k) = val;
    	C4[i+j*m].fastAccessDx(k) = val;
    	C5[i+j*m].fastAccessDx(k) = val;
      } 
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], l, beta, &C1[0], m);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], l, beta, &C2[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);

  unsigned int sz = m*l + l*n + m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], l, beta, &C3[0], m);

  COMPARE_FAD_VECTORS(C1, C3, m*n);

  sacado_blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &a[0], m, &b[0], l, beta, &C4[0], m);

  COMPARE_FAD_VECTORS(C1, C4, m*n);

  sacado_blas2.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &a[0], m, &b[0], l, beta, &C5[0], m);

  COMPARE_FAD_VECTORS(C1, C5, m*n);

  // transa
  teuchos_blas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], l, beta, &C1[0], m);
  sacado_blas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], l, beta, &C2[0], m);
  sacado_blas2.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], l, beta, &C3[0], m);
  sacado_blas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &a[0], l, &b[0], l, beta, &C4[0], m);
  sacado_blas2.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha, 
		    &a[0], l, &b[0], l, beta, &C5[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);
  COMPARE_FAD_VECTORS(C1, C3, m*n);
  COMPARE_FAD_VECTORS(C1, C4, m*n);
  COMPARE_FAD_VECTORS(C1, C5, m*n);

  // transb
  teuchos_blas.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], n, beta, &C1[0], m);
  sacado_blas.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], n, beta, &C2[0], m);
  sacado_blas2.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], m, &B[0], n, beta, &C3[0], m);
  sacado_blas.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &a[0], m, &b[0], n, beta, &C4[0], m);
  sacado_blas2.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &a[0], m, &b[0], n, beta, &C5[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);
  COMPARE_FAD_VECTORS(C1, C3, m*n);
  COMPARE_FAD_VECTORS(C1, C4, m*n);
  COMPARE_FAD_VECTORS(C1, C5, m*n);

  // transa and transb
  teuchos_blas.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], n, beta, &C1[0], m);
  sacado_blas.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], n, beta, &C2[0], m);
  sacado_blas2.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &A[0], l, &B[0], n, beta, &C3[0], m);
  sacado_blas.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &a[0], l, &b[0], n, beta, &C4[0], m);
  sacado_blas2.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha, 
		    &a[0], l, &b[0], n, beta, &C5[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);
  COMPARE_FAD_VECTORS(C1, C3, m*n);
  COMPARE_FAD_VECTORS(C1, C4, m*n);
  COMPARE_FAD_VECTORS(C1, C5, m*n);
}

// Tests all arguments, left side
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testSYMM1() {

  // SYMM is apparently not defined in the BLAS for complex types
  if (Teuchos::ScalarTraits<ScalarType>::isComplex)
    return; 

  VectorType A(m*m,ndot), B(m*n,ndot), C1(m*n,ndot), C2(m*n,ndot), C3(m*n,ndot);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*m].fastAccessDx(k) = urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      B[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
	B[i+j*m].fastAccessDx(k) = urand.number();
    }
  }
  FadType alpha(ndot, urand.number());
  FadType beta(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
    beta.fastAccessDx(k) = urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = urand.number();
      C1[i+j*m] = FadType(ndot, val);
      C2[i+j*m] = FadType(ndot, val);
      C3[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
	val = urand.number();
	C1[i+j*m].fastAccessDx(k) = val;
	C2[i+j*m].fastAccessDx(k) = val;
	C3[i+j*m].fastAccessDx(k) = val;
      } 
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C1[0], m);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C2[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);

  unsigned int sz = m*m*(1+ndot) + 2*m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C3[0], m);

  COMPARE_FAD_VECTORS(C1, C3, m*n);

  // lower tri
  teuchos_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C1[0], m);
  sacado_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C2[0], m);
  sacado_blas2.SYMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C3[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);
  COMPARE_FAD_VECTORS(C1, C3, m*n);
}

// Tests all arguments, right side
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testSYMM2() {

  // SYMM is apparently not defined in the BLAS for complex types
  if (Teuchos::ScalarTraits<ScalarType>::isComplex)
    return; 

  VectorType A(n*n,ndot), B(m*n,ndot), C1(m*n,ndot), C2(m*n,ndot), C3(m*n,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<n; i++) {
      A[i+j*n] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*n].fastAccessDx(k) = urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      B[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
	B[i+j*m].fastAccessDx(k) = urand.number();
    }
  }
  FadType alpha(ndot, urand.number());
  FadType beta(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
    beta.fastAccessDx(k) = urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = urand.number();
      C1[i+j*m] = FadType(ndot, val);
      C2[i+j*m] = FadType(ndot, val);
      C3[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
	val = urand.number();
	C1[i+j*m].fastAccessDx(k) = val;
	C2[i+j*m].fastAccessDx(k) = val;
	C3[i+j*m].fastAccessDx(k) = val;
      } 
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.SYMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &A[0], n, &B[0], m, beta, &C1[0], m);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.SYMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &A[0], n, &B[0], m, beta, &C2[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);

  unsigned int sz = n*n*(1+ndot) + 2*m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.SYMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &A[0], n, &B[0], m, beta, &C3[0], m);

  COMPARE_FAD_VECTORS(C1, C3, m*n);

  // lower tri
  teuchos_blas.SYMM(Teuchos::RIGHT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &A[0], n, &B[0], m, beta, &C1[0], m);
  sacado_blas.SYMM(Teuchos::RIGHT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &A[0], n, &B[0], m, beta, &C2[0], m);
  sacado_blas2.SYMM(Teuchos::RIGHT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &A[0], n, &B[0], m, beta, &C3[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);
  COMPARE_FAD_VECTORS(C1, C3, m*n);
}

// Tests different lda, ldb, ldc, left side
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testSYMM3() {

  // SYMM is apparently not defined in the BLAS for complex types
  if (Teuchos::ScalarTraits<ScalarType>::isComplex)
    return; 

  unsigned int lda = m+4;
  unsigned int ldb = m+5;
  unsigned int ldc = m+6;
  VectorType A(lda*m,ndot), B(ldb*n,ndot), C1(ldc*n,ndot), C2(ldc*n,ndot), 
    C3(ldc*n,ndot);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<lda; i++) {
      A[i+j*lda] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*lda].fastAccessDx(k) = urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<ldb; i++) {
      B[i+j*ldb] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
	B[i+j*ldb].fastAccessDx(k) = urand.number();
    }
  }
  FadType alpha(ndot, urand.number());
  FadType beta(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
    beta.fastAccessDx(k) = urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<ldc; i++) {
      ScalarType val = urand.number();
      C1[i+j*ldc] = FadType(ndot, val);
      C2[i+j*ldc] = FadType(ndot, val);
      C3[i+j*ldc] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
	val = urand.number();
	C1[i+j*ldc].fastAccessDx(k) = val;
	C2[i+j*ldc].fastAccessDx(k) = val;
	C3[i+j*ldc].fastAccessDx(k) = val;
      } 
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &A[0], lda, &B[0], ldb, beta, &C1[0], ldc);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &A[0], lda, &B[0], ldb, beta, &C2[0], ldc);

  COMPARE_FAD_VECTORS(C1, C2, ldc*n);

  unsigned int sz = m*m*(1+ndot) + 2*m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &A[0], lda, &B[0], ldb, beta, &C3[0], ldc);

  COMPARE_FAD_VECTORS(C1, C3, ldc*n);

  // lower tri
  teuchos_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &A[0], lda, &B[0], ldb, beta, &C1[0], ldc);
  sacado_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &A[0], lda, &B[0], ldb, beta, &C2[0], ldc);
  sacado_blas2.SYMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &A[0], lda, &B[0], ldb, beta, &C3[0], ldc);

  COMPARE_FAD_VECTORS(C1, C2, ldc*n);
  COMPARE_FAD_VECTORS(C1, C3, ldc*n);
}

// Tests different lda, ldb, ldc, right side
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testSYMM4() {

  // SYMM is apparently not defined in the BLAS for complex types
  if (Teuchos::ScalarTraits<ScalarType>::isComplex)
    return; 

  unsigned int lda = n+4;
  unsigned int ldb = m+5;
  unsigned int ldc = m+6;
  VectorType A(lda*n,ndot), B(ldb*n,ndot), C1(ldc*n,ndot), C2(ldc*n,ndot), 
    C3(ldc*n,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<lda; i++) {
      A[i+j*lda] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*lda].fastAccessDx(k) = urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<ldb; i++) {
      B[i+j*ldb] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
	B[i+j*ldb].fastAccessDx(k) = urand.number();
    }
  }
  FadType alpha(ndot, urand.number());
  FadType beta(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
    beta.fastAccessDx(k) = urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<ldc; i++) {
      ScalarType val = urand.number();
      C1[i+j*ldc] = FadType(ndot, val);
      C2[i+j*ldc] = FadType(ndot, val);
      C3[i+j*ldc] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
	val = urand.number();
	C1[i+j*ldc].fastAccessDx(k) = val;
	C2[i+j*ldc].fastAccessDx(k) = val;
	C3[i+j*ldc].fastAccessDx(k) = val;
      } 
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.SYMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &A[0], lda, &B[0], ldb, beta, &C1[0], ldc);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.SYMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &A[0], lda, &B[0], ldb, beta, &C2[0], ldc);

  COMPARE_FAD_VECTORS(C1, C2, ldc*n);

  unsigned int sz = n*n*(1+ndot) + 2*m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.SYMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &A[0], lda, &B[0], ldb, beta, &C3[0], ldc);

  COMPARE_FAD_VECTORS(C1, C3, ldc*n);

  // lower tri
  teuchos_blas.SYMM(Teuchos::RIGHT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &A[0], lda, &B[0], ldb, beta, &C1[0], ldc);
  sacado_blas.SYMM(Teuchos::RIGHT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &A[0], lda, &B[0], ldb, beta, &C2[0], ldc);
  sacado_blas2.SYMM(Teuchos::RIGHT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &A[0], lda, &B[0], ldb, beta, &C3[0], ldc);

  COMPARE_FAD_VECTORS(C1, C2, ldc*n);
  COMPARE_FAD_VECTORS(C1, C3, ldc*n);
}

// Tests with constant C
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testSYMM5() {

  // SYMM is apparently not defined in the BLAS for complex types
  if (Teuchos::ScalarTraits<ScalarType>::isComplex)
    return; 

  VectorType A(m*m,ndot), B(m*n,ndot), C1(m*n,ndot), C2(m*n,ndot), C3(m*n,ndot);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*m].fastAccessDx(k) = urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      B[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
	B[i+j*m].fastAccessDx(k) = urand.number();
    }
  }
  FadType alpha(ndot, urand.number());
  FadType beta(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
    beta.fastAccessDx(k) = urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = urand.number();
      C1[i+j*m] = val;
      C2[i+j*m] = val;
      C3[i+j*m] = val;
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C1[0], m);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C2[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);

  unsigned int sz = m*m*(1+ndot) + 2*m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C3[0], m);

  COMPARE_FAD_VECTORS(C1, C3, m*n);

  // lower tri
  teuchos_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C1[0], m);
  sacado_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C2[0], m);
  sacado_blas2.SYMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C3[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);
  COMPARE_FAD_VECTORS(C1, C3, m*n);
}

// Tests with constant alpha, beta
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testSYMM6() {

  // SYMM is apparently not defined in the BLAS for complex types
  if (Teuchos::ScalarTraits<ScalarType>::isComplex)
    return; 

  VectorType A(m*m,ndot), B(m*n,ndot), C1(m*n,ndot), C2(m*n,ndot), C3(m*n,ndot);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*m].fastAccessDx(k) = urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      B[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
	B[i+j*m].fastAccessDx(k) = urand.number();
    }
  }
  ScalarType alpha = urand.number();
  ScalarType beta = urand.number();

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = urand.number();
      C1[i+j*m] = FadType(ndot, val);
      C2[i+j*m] = FadType(ndot, val);
      C3[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
	val = urand.number();
	C1[i+j*m].fastAccessDx(k) = val;
	C2[i+j*m].fastAccessDx(k) = val;
	C3[i+j*m].fastAccessDx(k) = val;
      } 
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C1[0], m);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C2[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);

  unsigned int sz = m*m*(1+ndot) + 2*m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C3[0], m);

  COMPARE_FAD_VECTORS(C1, C3, m*n);

  // lower tri
  teuchos_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C1[0], m);
  sacado_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C2[0], m);
  sacado_blas2.SYMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C3[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);
  COMPARE_FAD_VECTORS(C1, C3, m*n);
}

// Tests with constant A
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testSYMM7() {

  // SYMM is apparently not defined in the BLAS for complex types
  if (Teuchos::ScalarTraits<ScalarType>::isComplex)
    return; 

  VectorType A(m*m,ndot), B(m*n,ndot), C1(m*n,ndot), C2(m*n,ndot), C3(m*n,ndot),
    C4(m*n,ndot), C5(m*n,ndot);
  std::vector<ScalarType> a(m*m);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<m; i++) {
      a[i+j*m] = urand.number();
      A[i+j*m] = a[i+j*m];
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      B[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
	B[i+j*m].fastAccessDx(k) = urand.number();
    }
  }
  FadType alpha(ndot, urand.number());
  FadType beta(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
    beta.fastAccessDx(k) = urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = urand.number();
      C1[i+j*m] = FadType(ndot, val);
      C2[i+j*m] = FadType(ndot, val);
      C3[i+j*m] = FadType(ndot, val);
      C4[i+j*m] = FadType(ndot, val);
      C5[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
	val = urand.number();
	C1[i+j*m].fastAccessDx(k) = val;
	C2[i+j*m].fastAccessDx(k) = val;
	C3[i+j*m].fastAccessDx(k) = val;
	C4[i+j*m].fastAccessDx(k) = val;
	C5[i+j*m].fastAccessDx(k) = val;
      } 
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C1[0], m);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C2[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);

  unsigned int sz = m*m + 2*m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C3[0], m);

  COMPARE_FAD_VECTORS(C1, C3, m*n);

  sacado_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &a[0], m, &B[0], m, beta, &C4[0], m);

  COMPARE_FAD_VECTORS(C1, C4, m*n);

  sacado_blas2.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &a[0], m, &B[0], m, beta, &C5[0], m);

  COMPARE_FAD_VECTORS(C1, C5, m*n);

  // lower tri
  teuchos_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C1[0], m);
  sacado_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C2[0], m);
  sacado_blas2.SYMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C3[0], m);
  sacado_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &a[0], m, &B[0], m, beta, &C4[0], m);
  sacado_blas2.SYMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &a[0], m, &B[0], m, beta, &C5[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);
  COMPARE_FAD_VECTORS(C1, C3, m*n);
  COMPARE_FAD_VECTORS(C1, C4, m*n);
  COMPARE_FAD_VECTORS(C1, C5, m*n);
}

// Tests with constant B
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testSYMM8() {

  // SYMM is apparently not defined in the BLAS for complex types
  if (Teuchos::ScalarTraits<ScalarType>::isComplex)
    return; 

  VectorType A(m*m,ndot), B(m*n,ndot), C1(m*n,ndot), C2(m*n,ndot), C3(m*n,ndot),
    C4(m*n,ndot), C5(m*n,ndot);
  std::vector<ScalarType> b(m*n);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*m].fastAccessDx(k) = urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      b[i+j*m] = urand.number();
      B[i+j*m] = b[i+j*m];
    }
  }
  FadType alpha(ndot, urand.number());
  FadType beta(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
    beta.fastAccessDx(k) = urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = urand.number();
      C1[i+j*m] = FadType(ndot, val);
      C2[i+j*m] = FadType(ndot, val);
      C3[i+j*m] = FadType(ndot, val);
      C4[i+j*m] = FadType(ndot, val);
      C5[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
	val = urand.number();
	C1[i+j*m].fastAccessDx(k) = val;
	C2[i+j*m].fastAccessDx(k) = val;
	C3[i+j*m].fastAccessDx(k) = val;
	C4[i+j*m].fastAccessDx(k) = val;
	C5[i+j*m].fastAccessDx(k) = val;
      } 
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C1[0], m);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C2[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);

  unsigned int sz = m*m*(1+ndot) + m*n*(2+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C3[0], m);

  COMPARE_FAD_VECTORS(C1, C3, m*n);

  sacado_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &A[0], m, &b[0], m, beta, &C4[0], m);

  COMPARE_FAD_VECTORS(C1, C4, m*n);

  sacado_blas2.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &A[0], m, &b[0], m, beta, &C5[0], m);

  COMPARE_FAD_VECTORS(C1, C5, m*n);

  // lower tri
  teuchos_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C1[0], m);
  sacado_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C2[0], m);
  sacado_blas2.SYMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C3[0], m);
  sacado_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &A[0], m, &b[0], m, beta, &C4[0], m);
  sacado_blas2.SYMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &A[0], m, &b[0], m, beta, &C5[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);
  COMPARE_FAD_VECTORS(C1, C3, m*n);
  COMPARE_FAD_VECTORS(C1, C4, m*n);
  COMPARE_FAD_VECTORS(C1, C5, m*n);
}

// Tests with constant A and B
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testSYMM9() {

  // SYMM is apparently not defined in the BLAS for complex types
  if (Teuchos::ScalarTraits<ScalarType>::isComplex)
    return; 

  VectorType A(m*m,ndot), B(m*n,ndot), C1(m*n,ndot), C2(m*n,ndot), C3(m*n,ndot),
    C4(m*n,ndot), C5(m*n,ndot);
  std::vector<ScalarType> a(m*m), b(m*n);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<m; i++) {
      a[i+j*m] = urand.number();
      A[i+j*m] = a[i+j*m];
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      b[i+j*m] = urand.number();
      B[i+j*m] = b[i+j*m];
    }
  }
  FadType alpha(ndot, urand.number());
  FadType beta(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
    beta.fastAccessDx(k) = urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = urand.number();
      C1[i+j*m] = FadType(ndot, val);
      C2[i+j*m] = FadType(ndot, val);
      C3[i+j*m] = FadType(ndot, val);
      C4[i+j*m] = FadType(ndot, val);
      C5[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
	val = urand.number();
	C1[i+j*m].fastAccessDx(k) = val;
	C2[i+j*m].fastAccessDx(k) = val;
	C3[i+j*m].fastAccessDx(k) = val;
	C4[i+j*m].fastAccessDx(k) = val;
	C5[i+j*m].fastAccessDx(k) = val;
      } 
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C1[0], m);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C2[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);

  unsigned int sz = m*m + m*n*(2+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C3[0], m);

  COMPARE_FAD_VECTORS(C1, C3, m*n);

  sacado_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &a[0], m, &b[0], m, beta, &C4[0], m);

  COMPARE_FAD_VECTORS(C1, C4, m*n);

  sacado_blas2.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha, 
		    &a[0], m, &b[0], m, beta, &C5[0], m);

  COMPARE_FAD_VECTORS(C1, C5, m*n);

  // lower tri
  teuchos_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C1[0], m);
  sacado_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C2[0], m);
  sacado_blas2.SYMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &A[0], m, &B[0], m, beta, &C3[0], m);
  sacado_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &a[0], m, &b[0], m, beta, &C4[0], m);
  sacado_blas2.SYMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, m, n, alpha, 
		    &a[0], m, &b[0], m, beta, &C5[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);
  COMPARE_FAD_VECTORS(C1, C3, m*n);
  COMPARE_FAD_VECTORS(C1, C4, m*n);
  COMPARE_FAD_VECTORS(C1, C5, m*n);
}

// Tests all arguments, left side
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testTRMM1() {
  VectorType A(m*m,ndot), B1(m*n,ndot), B2(m*n,ndot), B3(m*n,ndot);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*m].fastAccessDx(k) = urand.number();
    }
  }
  FadType alpha(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = urand.number();
      B1[i+j*m] = FadType(ndot, val);
      B2[i+j*m] = FadType(ndot, val);
      B3[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
      	val = urand.number();
      	B1[i+j*m].fastAccessDx(k) = val;
      	B2[i+j*m].fastAccessDx(k) = val;
      	B3[i+j*m].fastAccessDx(k) = val;
      } 
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);

  COMPARE_FAD_VECTORS(B1, B2, m*n);

  unsigned int sz = m*m*(1+ndot) + m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);

  COMPARE_FAD_VECTORS(B1, B3, m*n);

  teuchos_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);
  sacado_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);
  sacado_blas2.TRMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);

  teuchos_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);
  sacado_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);
  sacado_blas2.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);

  for (unsigned int i=0; i<m; i++) {
    A[i*m+i].val() = 1.0;
    for (unsigned int k=0; k<ndot; k++)
      A[i*m+i].fastAccessDx(k) = 0.0;
  }
  teuchos_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);
  sacado_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);
  sacado_blas2.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);
}

// Tests all arguments, right side
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testTRMM2() {
  VectorType A(n*n,ndot), B1(m*n,ndot), B2(m*n,ndot), B3(m*n,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<n; i++) {
      A[i+j*n] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*n].fastAccessDx(k) = urand.number();
    }
  }
  FadType alpha(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = urand.number();
      B1[i+j*m] = FadType(ndot, val);
      B2[i+j*m] = FadType(ndot, val);
      B3[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
      	val = urand.number();
      	B1[i+j*m].fastAccessDx(k) = val;
      	B2[i+j*m].fastAccessDx(k) = val;
      	B3[i+j*m].fastAccessDx(k) = val;
      } 
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], n, &B1[0], m);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.TRMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], n, &B2[0], m);

  COMPARE_FAD_VECTORS(B1, B2, m*n);

  unsigned int sz = n*n*(1+ndot) + m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.TRMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], n, &B3[0], m);

  COMPARE_FAD_VECTORS(B1, B3, m*n);

  teuchos_blas.TRMM(Teuchos::RIGHT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], n, &B1[0], m);
  sacado_blas.TRMM(Teuchos::RIGHT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], n, &B2[0], m);
  sacado_blas2.TRMM(Teuchos::RIGHT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], n, &B3[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);

  teuchos_blas.TRMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], n, &B1[0], m);
  sacado_blas.TRMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], n, &B2[0], m);
  sacado_blas2.TRMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], n, &B3[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);

  for (unsigned int i=0; i<n; i++) {
    A[i*n+i].val() = 1.0;
    for (unsigned int k=0; k<ndot; k++)
      A[i*n+i].fastAccessDx(k) = 0.0;
  }
  teuchos_blas.TRMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], n, &B1[0], m);
  sacado_blas.TRMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], n, &B2[0], m);
  sacado_blas2.TRMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], n, &B3[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);
}

// Tests all arguments, left side, different lda, ldb
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testTRMM3() {
  unsigned int lda = m+4;
  unsigned int ldb = m+5;
  VectorType A(lda*m,ndot), B1(ldb*n,ndot), B2(ldb*n,ndot), B3(ldb*n,ndot);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<lda; i++) {
      A[i+j*lda] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*lda].fastAccessDx(k) = urand.number();
    }
  }
  FadType alpha(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<ldb; i++) {
      ScalarType val = urand.number();
      B1[i+j*ldb] = FadType(ndot, val);
      B2[i+j*ldb] = FadType(ndot, val);
      B3[i+j*ldb] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
      	val = urand.number();
      	B1[i+j*ldb].fastAccessDx(k) = val;
      	B2[i+j*ldb].fastAccessDx(k) = val;
      	B3[i+j*ldb].fastAccessDx(k) = val;
      } 
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B1[0], ldb);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B2[0], ldb);

  COMPARE_FAD_VECTORS(B1, B2, ldb*n);

  unsigned int sz = m*m*(1+ndot) + m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B3[0], ldb);

  COMPARE_FAD_VECTORS(B1, B3, ldb*n);

  teuchos_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B1[0], ldb);
  sacado_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B2[0], ldb);
  sacado_blas2.TRMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B3[0], ldb);
  COMPARE_FAD_VECTORS(B1, B2, ldb*n);
  COMPARE_FAD_VECTORS(B1, B3, ldb*n);

  teuchos_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B1[0], ldb);
  sacado_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B2[0], ldb);
  sacado_blas2.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B3[0], ldb);
  COMPARE_FAD_VECTORS(B1, B2, ldb*n);
  COMPARE_FAD_VECTORS(B1, B3, ldb*n);

  for (unsigned int i=0; i<m; i++) {
    A[i*lda+i].val() = 1.0;
    for (unsigned int k=0; k<ndot; k++)
      A[i*lda+i].fastAccessDx(k) = 0.0;
  }
  teuchos_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], lda, &B1[0], ldb);
  sacado_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], lda, &B2[0], ldb);
  sacado_blas2.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], lda, &B3[0], ldb);
  COMPARE_FAD_VECTORS(B1, B2, ldb*n);
  COMPARE_FAD_VECTORS(B1, B3, ldb*n);
}

// Tests all arguments, right side, different lda, ldb
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testTRMM4() {
  unsigned int lda = n+4;
  unsigned int ldb = m+5;
  VectorType A(lda*n,ndot), B1(ldb*n,ndot), B2(ldb*n,ndot), B3(ldb*n,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<lda; i++) {
      A[i+j*lda] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*lda].fastAccessDx(k) = urand.number();
    }
  }
  FadType alpha(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<ldb; i++) {
      ScalarType val = urand.number();
      B1[i+j*ldb] = FadType(ndot, val);
      B2[i+j*ldb] = FadType(ndot, val);
      B3[i+j*ldb] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
      	val = urand.number();
      	B1[i+j*ldb].fastAccessDx(k) = val;
      	B2[i+j*ldb].fastAccessDx(k) = val;
      	B3[i+j*ldb].fastAccessDx(k) = val;
      } 
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B1[0], ldb);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.TRMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		   Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B2[0], ldb);

  COMPARE_FAD_VECTORS(B1, B2, ldb*n);

  unsigned int sz = n*n*(1+ndot) + m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.TRMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B3[0], ldb);

  COMPARE_FAD_VECTORS(B1, B3, ldb*n);

  teuchos_blas.TRMM(Teuchos::RIGHT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B1[0], ldb);
  sacado_blas.TRMM(Teuchos::RIGHT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		   Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B2[0], ldb);
  sacado_blas2.TRMM(Teuchos::RIGHT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B3[0], ldb);
  COMPARE_FAD_VECTORS(B1, B2, ldb*n);
  COMPARE_FAD_VECTORS(B1, B3, ldb*n);

  teuchos_blas.TRMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B1[0], ldb);
  sacado_blas.TRMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		   Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B2[0], ldb);
  sacado_blas2.TRMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B3[0], ldb);
  COMPARE_FAD_VECTORS(B1, B2, ldb*n);
  COMPARE_FAD_VECTORS(B1, B3, ldb*n);

  for (unsigned int i=0; i<n; i++) {
    A[i*lda+i].val() = 1.0;
    for (unsigned int k=0; k<ndot; k++)
      A[i*lda+i].fastAccessDx(k) = 0.0;
  }
  teuchos_blas.TRMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], lda, &B1[0], ldb);
  sacado_blas.TRMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], lda, &B2[0], ldb);
  sacado_blas2.TRMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], lda, &B3[0], ldb);
  COMPARE_FAD_VECTORS(B1, B2, ldb*n);
  COMPARE_FAD_VECTORS(B1, B3, ldb*n);
}

// Tests constant alpha
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testTRMM5() {
  VectorType A(m*m,ndot), B1(m*n,ndot), B2(m*n,ndot), B3(m*n,ndot);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*m].fastAccessDx(k) = urand.number();
    }
  }
  ScalarType alpha = urand.number();

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = urand.number();
      B1[i+j*m] = FadType(ndot, val);
      B2[i+j*m] = FadType(ndot, val);
      B3[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
      	val = urand.number();
      	B1[i+j*m].fastAccessDx(k) = val;
      	B2[i+j*m].fastAccessDx(k) = val;
      	B3[i+j*m].fastAccessDx(k) = val;
      } 
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);

  COMPARE_FAD_VECTORS(B1, B2, m*n);

  unsigned int sz = m*m*(1+ndot) + m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);

  COMPARE_FAD_VECTORS(B1, B3, m*n);

  teuchos_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);
  sacado_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);
  sacado_blas2.TRMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);

  teuchos_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);
  sacado_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);
  sacado_blas2.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);

  for (unsigned int i=0; i<m; i++) {
    A[i*m+i].val() = 1.0;
    for (unsigned int k=0; k<ndot; k++)
      A[i*m+i].fastAccessDx(k) = 0.0;
  }
  teuchos_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);
  sacado_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);
  sacado_blas2.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);
}

// Tests constant B
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testTRMM6() {
  VectorType A(m*m,ndot), B1(m*n,ndot), B2(m*n,ndot), B3(m*n,ndot);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*m].fastAccessDx(k) = urand.number();
    }
  }
  FadType alpha(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = urand.number();
      B1[i+j*m] = val;
      B2[i+j*m] = val;
      B3[i+j*m] = val;
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);

  COMPARE_FAD_VECTORS(B1, B2, m*n);

  unsigned int sz = m*m*(1+ndot) + m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);

  COMPARE_FAD_VECTORS(B1, B3, m*n);

  teuchos_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);
  sacado_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);
  sacado_blas2.TRMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);

  teuchos_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);
  sacado_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);
  sacado_blas2.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);

  for (unsigned int i=0; i<m; i++) {
    A[i*m+i].val() = 1.0;
    for (unsigned int k=0; k<ndot; k++)
      A[i*m+i].fastAccessDx(k) = 0.0;
  }
  teuchos_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);
  sacado_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);
  sacado_blas2.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);
}

// Tests constant A
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testTRMM7() {
  VectorType A(m*m,ndot), B1(m*n,ndot), B2(m*n,ndot), B3(m*n,ndot), 
    B4(m*n,ndot), B5(m*n,ndot);
  std::vector<ScalarType> a(m*m);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<m; i++) {
      a[i+j*m] = urand.number();
      A[i+j*m] = a[i+j*m];
    }
  }
  FadType alpha(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = urand.number();
      B1[i+j*m] = FadType(ndot, val);
      B2[i+j*m] = FadType(ndot, val);
      B3[i+j*m] = FadType(ndot, val);
      B4[i+j*m] = FadType(ndot, val);
      B5[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
      	val = urand.number();
      	B1[i+j*m].fastAccessDx(k) = val;
      	B2[i+j*m].fastAccessDx(k) = val;
      	B3[i+j*m].fastAccessDx(k) = val;
	B4[i+j*m].fastAccessDx(k) = val;
	B5[i+j*m].fastAccessDx(k) = val;
      } 
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);

  COMPARE_FAD_VECTORS(B1, B2, m*n);

  unsigned int sz = m*m + m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);

  COMPARE_FAD_VECTORS(B1, B3, m*n);

  sacado_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &a[0], m, &B4[0], m);

  COMPARE_FAD_VECTORS(B1, B4, m*n);

  sacado_blas2.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &a[0], m, &B5[0], m);

  COMPARE_FAD_VECTORS(B1, B5, m*n);

  teuchos_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);
  sacado_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);
  sacado_blas2.TRMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);
  sacado_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &a[0], m, &B4[0], m);
  sacado_blas2.TRMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &a[0], m, &B5[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);
  COMPARE_FAD_VECTORS(B1, B4, m*n);
  COMPARE_FAD_VECTORS(B1, B5, m*n);

  teuchos_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);
  sacado_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);
  sacado_blas2.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);
  sacado_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &a[0], m, &B4[0], m);
  sacado_blas2.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &a[0], m, &B5[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);
  COMPARE_FAD_VECTORS(B1, B4, m*n);
  COMPARE_FAD_VECTORS(B1, B5, m*n);

  for (unsigned int i=0; i<m; i++) {
    A[i*m+i].val() = 1.0;
    for (unsigned int k=0; k<ndot; k++)
      A[i*m+i].fastAccessDx(k) = 0.0;
  }
  teuchos_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);
  sacado_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);
  sacado_blas2.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);
  sacado_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &a[0], m, &B4[0], m);
  sacado_blas2.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &a[0], m, &B5[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);
  COMPARE_FAD_VECTORS(B1, B4, m*n);
  COMPARE_FAD_VECTORS(B1, B5, m*n);
}

// Tests all arguments, left side
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testTRSM1() {
  VectorType A(m*m,ndot), B1(m*n,ndot), B2(m*n,ndot), B3(m*n,ndot);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<m; i++) {
      //A[i+j*m] = urand.number();
      A[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*m].fastAccessDx(k) = urand.number();
    }
  }
  FadType alpha(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
  }
  //ScalarType alpha = urand.number();

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = urand.number();
      // B1[i+j*m] = val;
      // B2[i+j*m] = val;
      // B3[i+j*m] = val;
      B1[i+j*m] = FadType(ndot, val);
      B2[i+j*m] = FadType(ndot, val);
      B3[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
      	val = urand.number();
      	B1[i+j*m].fastAccessDx(k) = val;
      	B2[i+j*m].fastAccessDx(k) = val;
      	B3[i+j*m].fastAccessDx(k) = val;
      } 
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);

  COMPARE_FAD_VECTORS(B1, B2, m*n);

  unsigned int sz = m*m*(1+ndot) + m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);

  COMPARE_FAD_VECTORS(B1, B3, m*n);

  teuchos_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
  		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);
  sacado_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
  		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);
  sacado_blas2.TRSM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
  		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);

  teuchos_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
  		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);
  sacado_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
  		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);
  sacado_blas2.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
  		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);

  for (unsigned int i=0; i<m; i++) {
    A[i*m+i].val() = 1.0;
    for (unsigned int k=0; k<ndot; k++)
      A[i*m+i].fastAccessDx(k) = 0.0;
  }
  teuchos_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
  		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);
  sacado_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
  		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);
  sacado_blas2.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
  		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);
}

// Tests all arguments, right side
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testTRSM2() {
  VectorType A(n*n,ndot), B1(m*n,ndot), B2(m*n,ndot), B3(m*n,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<n; i++) {
      A[i+j*n] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*n].fastAccessDx(k) = urand.number();
    }
  }
  FadType alpha(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = urand.number();
      B1[i+j*m] = FadType(ndot, val);
      B2[i+j*m] = FadType(ndot, val);
      B3[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
      	val = urand.number();
      	B1[i+j*m].fastAccessDx(k) = val;
      	B2[i+j*m].fastAccessDx(k) = val;
      	B3[i+j*m].fastAccessDx(k) = val;
      } 
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRSM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], n, &B1[0], m);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.TRSM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], n, &B2[0], m);

  COMPARE_FAD_VECTORS(B1, B2, m*n);

  unsigned int sz = n*n*(1+ndot) + m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.TRSM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], n, &B3[0], m);

  COMPARE_FAD_VECTORS(B1, B3, m*n);

  teuchos_blas.TRSM(Teuchos::RIGHT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], n, &B1[0], m);
  sacado_blas.TRSM(Teuchos::RIGHT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], n, &B2[0], m);
  sacado_blas2.TRSM(Teuchos::RIGHT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], n, &B3[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);

  teuchos_blas.TRSM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], n, &B1[0], m);
  sacado_blas.TRSM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], n, &B2[0], m);
  sacado_blas2.TRSM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], n, &B3[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);

  for (unsigned int i=0; i<n; i++) {
    A[i*n+i].val() = 1.0;
    for (unsigned int k=0; k<ndot; k++)
      A[i*n+i].fastAccessDx(k) = 0.0;
  }
  teuchos_blas.TRSM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], n, &B1[0], m);
  sacado_blas.TRSM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], n, &B2[0], m);
  sacado_blas2.TRSM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], n, &B3[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);
}

// Tests all arguments, left side, different lda, ldb
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testTRSM3() {
  unsigned int lda = m+4;
  unsigned int ldb = m+5;
  VectorType A(lda*m,ndot), B1(ldb*n,ndot), B2(ldb*n,ndot), B3(ldb*n,ndot);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<lda; i++) {
      A[i+j*lda] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*lda].fastAccessDx(k) = urand.number();
    }
  }
  FadType alpha(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<ldb; i++) {
      ScalarType val = urand.number();
      B1[i+j*ldb] = FadType(ndot, val);
      B2[i+j*ldb] = FadType(ndot, val);
      B3[i+j*ldb] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
      	val = urand.number();
      	B1[i+j*ldb].fastAccessDx(k) = val;
      	B2[i+j*ldb].fastAccessDx(k) = val;
      	B3[i+j*ldb].fastAccessDx(k) = val;
      } 
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B1[0], ldb);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B2[0], ldb);

  COMPARE_FAD_VECTORS(B1, B2, ldb*n);

  unsigned int sz = m*m*(1+ndot) + m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B3[0], ldb);

  COMPARE_FAD_VECTORS(B1, B3, ldb*n);

  teuchos_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B1[0], ldb);
  sacado_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B2[0], ldb);
  sacado_blas2.TRSM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B3[0], ldb);
  COMPARE_FAD_VECTORS(B1, B2, ldb*n);
  COMPARE_FAD_VECTORS(B1, B3, ldb*n);

  teuchos_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B1[0], ldb);
  sacado_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B2[0], ldb);
  sacado_blas2.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B3[0], ldb);
  COMPARE_FAD_VECTORS(B1, B2, ldb*n);
  COMPARE_FAD_VECTORS(B1, B3, ldb*n);

  for (unsigned int i=0; i<m; i++) {
    A[i*lda+i].val() = 1.0;
    for (unsigned int k=0; k<ndot; k++)
      A[i*lda+i].fastAccessDx(k) = 0.0;
  }
  teuchos_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], lda, &B1[0], ldb);
  sacado_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], lda, &B2[0], ldb);
  sacado_blas2.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], lda, &B3[0], ldb);
  COMPARE_FAD_VECTORS(B1, B2, ldb*n);
  COMPARE_FAD_VECTORS(B1, B3, ldb*n);
}

// Tests all arguments, right side, different lda, ldb
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testTRSM4() {
  unsigned int lda = n+4;
  unsigned int ldb = m+5;
  VectorType A(lda*n,ndot), B1(ldb*n,ndot), B2(ldb*n,ndot), B3(ldb*n,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<lda; i++) {
      A[i+j*lda] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*lda].fastAccessDx(k) = urand.number();
    }
  }
  FadType alpha(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<ldb; i++) {
      ScalarType val = urand.number();
      B1[i+j*ldb] = FadType(ndot, val);
      B2[i+j*ldb] = FadType(ndot, val);
      B3[i+j*ldb] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
      	val = urand.number();
      	B1[i+j*ldb].fastAccessDx(k) = val;
      	B2[i+j*ldb].fastAccessDx(k) = val;
      	B3[i+j*ldb].fastAccessDx(k) = val;
      } 
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRSM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B1[0], ldb);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.TRSM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		   Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B2[0], ldb);

  COMPARE_FAD_VECTORS(B1, B2, ldb*n);

  unsigned int sz = n*n*(1+ndot) + m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.TRSM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B3[0], ldb);

  COMPARE_FAD_VECTORS(B1, B3, ldb*n);

  teuchos_blas.TRSM(Teuchos::RIGHT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B1[0], ldb);
  sacado_blas.TRSM(Teuchos::RIGHT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		   Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B2[0], ldb);
  sacado_blas2.TRSM(Teuchos::RIGHT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B3[0], ldb);
  COMPARE_FAD_VECTORS(B1, B2, ldb*n);
  COMPARE_FAD_VECTORS(B1, B3, ldb*n);

  teuchos_blas.TRSM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B1[0], ldb);
  sacado_blas.TRSM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		   Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B2[0], ldb);
  sacado_blas2.TRSM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B3[0], ldb);
  COMPARE_FAD_VECTORS(B1, B2, ldb*n);
  COMPARE_FAD_VECTORS(B1, B3, ldb*n);

  for (unsigned int i=0; i<n; i++) {
    A[i*lda+i].val() = 1.0;
    for (unsigned int k=0; k<ndot; k++)
      A[i*lda+i].fastAccessDx(k) = 0.0;
  }
  teuchos_blas.TRSM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], lda, &B1[0], ldb);
  sacado_blas.TRSM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], lda, &B2[0], ldb);
  sacado_blas2.TRSM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], lda, &B3[0], ldb);
  COMPARE_FAD_VECTORS(B1, B2, ldb*n);
  COMPARE_FAD_VECTORS(B1, B3, ldb*n);
}

// Tests constant alpha
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testTRSM5() {
  VectorType A(m*m,ndot), B1(m*n,ndot), B2(m*n,ndot), B3(m*n,ndot);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*m].fastAccessDx(k) = urand.number();
    }
  }
  ScalarType alpha = urand.number();

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = urand.number();
      B1[i+j*m] = FadType(ndot, val);
      B2[i+j*m] = FadType(ndot, val);
      B3[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
      	val = urand.number();
      	B1[i+j*m].fastAccessDx(k) = val;
      	B2[i+j*m].fastAccessDx(k) = val;
      	B3[i+j*m].fastAccessDx(k) = val;
      } 
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);

  COMPARE_FAD_VECTORS(B1, B2, m*n);

  unsigned int sz = m*m*(1+ndot) + m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);

  COMPARE_FAD_VECTORS(B1, B3, m*n);

  teuchos_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);
  sacado_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);
  sacado_blas2.TRSM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);

  teuchos_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);
  sacado_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);
  sacado_blas2.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);

  for (unsigned int i=0; i<m; i++) {
    A[i*m+i].val() = 1.0;
    for (unsigned int k=0; k<ndot; k++)
      A[i*m+i].fastAccessDx(k) = 0.0;
  }
  teuchos_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);
  sacado_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);
  sacado_blas2.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);
}

// Tests constant B
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testTRSM6() {
  VectorType A(m*m,ndot), B1(m*n,ndot), B2(m*n,ndot), B3(m*n,ndot);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, urand.number());
      for (unsigned int k=0; k<ndot; k++)
      	A[i+j*m].fastAccessDx(k) = urand.number();
    }
  }
  FadType alpha(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = urand.number();
      B1[i+j*m] = val;
      B2[i+j*m] = val;
      B3[i+j*m] = val;
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);

  COMPARE_FAD_VECTORS(B1, B2, m*n);

  unsigned int sz = m*m*(1+ndot) + m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);

  COMPARE_FAD_VECTORS(B1, B3, m*n);

  teuchos_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);
  sacado_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);
  sacado_blas2.TRSM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);

  teuchos_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);
  sacado_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);
  sacado_blas2.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);

  for (unsigned int i=0; i<m; i++) {
    A[i*m+i].val() = 1.0;
    for (unsigned int k=0; k<ndot; k++)
      A[i*m+i].fastAccessDx(k) = 0.0;
  }
  teuchos_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);
  sacado_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);
  sacado_blas2.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);
}

// Tests constant A
template <class FadType, class ScalarType>
void
FadBLASUnitTests<FadType,ScalarType>::
testTRSM7() {
  VectorType A(m*m,ndot), B1(m*n,ndot), B2(m*n,ndot), B3(m*n,ndot), 
    B4(m*n,ndot), B5(m*n,ndot);
  std::vector<ScalarType> a(m*m);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<m; i++) {
      a[i+j*m] = urand.number();
      A[i+j*m] = a[i+j*m];
    }
  }
  FadType alpha(ndot, urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = urand.number();
      B1[i+j*m] = FadType(ndot, val);
      B2[i+j*m] = FadType(ndot, val);
      B3[i+j*m] = FadType(ndot, val);
      B4[i+j*m] = FadType(ndot, val);
      B5[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
      	val = urand.number();
      	B1[i+j*m].fastAccessDx(k) = val;
      	B2[i+j*m].fastAccessDx(k) = val;
      	B3[i+j*m].fastAccessDx(k) = val;
	B4[i+j*m].fastAccessDx(k) = val;
	B5[i+j*m].fastAccessDx(k) = val;
      } 
    }
  }
  
  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);

  Sacado::Fad::BLAS<int,FadType> sacado_blas;
  sacado_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);

  COMPARE_FAD_VECTORS(B1, B2, m*n);

  unsigned int sz = m*m + m*n*(1+ndot);
  Sacado::Fad::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);

  COMPARE_FAD_VECTORS(B1, B3, m*n);

  sacado_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &a[0], m, &B4[0], m);

  COMPARE_FAD_VECTORS(B1, B4, m*n);

  sacado_blas2.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &a[0], m, &B5[0], m);

  COMPARE_FAD_VECTORS(B1, B5, m*n);

  teuchos_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);
  sacado_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);
  sacado_blas2.TRSM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);
  sacado_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &a[0], m, &B4[0], m);
  sacado_blas2.TRSM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &a[0], m, &B5[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);
  COMPARE_FAD_VECTORS(B1, B4, m*n);
  COMPARE_FAD_VECTORS(B1, B5, m*n);

  teuchos_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);
  sacado_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);
  sacado_blas2.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);
  sacado_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &a[0], m, &B4[0], m);
  sacado_blas2.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, 
		    Teuchos::NON_UNIT_DIAG, m, n, alpha, &a[0], m, &B5[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);
  COMPARE_FAD_VECTORS(B1, B4, m*n);
  COMPARE_FAD_VECTORS(B1, B5, m*n);

  for (unsigned int i=0; i<m; i++) {
    A[i*m+i].val() = 1.0;
    for (unsigned int k=0; k<ndot; k++)
      A[i*m+i].fastAccessDx(k) = 0.0;
  }
  teuchos_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);
  sacado_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);
  sacado_blas2.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &A[0], m, &B3[0], m);
  sacado_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &a[0], m, &B4[0], m);
  sacado_blas2.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
		    Teuchos::UNIT_DIAG, m, n, alpha, &a[0], m, &B5[0], m);
  COMPARE_FAD_VECTORS(B1, B2, m*n);
  COMPARE_FAD_VECTORS(B1, B3, m*n);
  COMPARE_FAD_VECTORS(B1, B4, m*n);
  COMPARE_FAD_VECTORS(B1, B5, m*n);
}

#undef COMPARE_VALUES
#undef COMPARE_FADS
#undef COMPARE_FAD_VECTORS

#endif // FADBLASUNITTESTS_HPP
