// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef FADBLASUNITTESTS_HPP
#define FADBLASUNITTESTS_HPP

// Sacado includes
#include "Sacado_No_Kokkos.hpp"
#include "Sacado_Fad_BLAS.hpp"
#include "Sacado_Random.hpp"

// gtest includes
#include <gtest/gtest.h>

#include "GTestUtils.hpp"

#define COMPARE_FAD_VECTORS(X1, X2, n)          \
  ASSERT_TRUE(X1.size() == std::size_t(n));  \
  ASSERT_TRUE(X2.size() == std::size_t(n));  \
  for (unsigned int i=0; i<n; i++) {            \
    COMPARE_FADS(X1[i], X2[i]);                 \
  }                                             \
  ;

// A class for testing differentiated BLAS operations for general Fad types
template <class FadType>
class FadBLASUnitTests : public ::testing::Test {
protected:
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;

  // Random number generator
  Sacado::Random<ScalarType> urand;

  // Real random number generator for derivative components
  Sacado::Random<double> real_urand;

  // Number of matrix rows
  unsigned int m_;

  // Number of matrix columns
  unsigned int n_;

  // Number of matrix columns for level 3 blas
  unsigned int l_;

  // Number of derivative components
  unsigned int ndot_;

  // Tolerances to which fad objects should be the same
  double tol_a, tol_r;

  FadType fad;

  FadBLASUnitTests() :
    urand(), real_urand(), m_(5), n_(6), l_(4), ndot_(7),
    tol_a(1.0e-11), tol_r(1.0e-11) {}

}; // class FadBLASUnitTests

TYPED_TEST_SUITE_P(FadBLASUnitTests);

// Tests all arguments
TYPED_TEST_P(FadBLASUnitTests, testSCAL1) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto ndot = this->ndot_;

  VectorType x1(m,ndot), x2(m,ndot), x3(m,ndot);
  for (unsigned int i=0; i<m; i++) {
    ScalarType val = this->urand.number();
    x1[i] = FadType(ndot, val);
    x2[i] = FadType(ndot, val);
    x3[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = this->urand.number();
      x1[i].fastAccessDx(k) = val;
      x2[i].fastAccessDx(k) = val;
      x3[i].fastAccessDx(k) = val;
    }
  }
  FadType alpha(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.SCAL(m, alpha, &x1[0], 1);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.SCAL(m, alpha, &x2[0], 1);

  COMPARE_FAD_VECTORS(x1, x2, m);

  unsigned int sz = m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.SCAL(m, alpha, &x3[0], 1);

  COMPARE_FAD_VECTORS(x1, x3, m);
}

// Tests non-unit inc
TYPED_TEST_P(FadBLASUnitTests, testSCAL2) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto ndot = this->ndot_;

  unsigned int incx = 2;
  VectorType x1(m*incx,ndot), x2(m*incx,ndot), x3(m*incx,ndot);
  for (unsigned int i=0; i<m*incx; i++) {
    ScalarType val = this->urand.number();
    x1[i] = FadType(ndot, val);
    x2[i] = FadType(ndot, val);
    x3[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = this->urand.number();
      x1[i].fastAccessDx(k) = val;
      x2[i].fastAccessDx(k) = val;
      x3[i].fastAccessDx(k) = val;
    }
  }
  FadType alpha(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.SCAL(m, alpha, &x1[0], incx);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.SCAL(m, alpha, &x2[0], incx);

  COMPARE_FAD_VECTORS(x1, x2, m*incx);

  unsigned int sz = m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.SCAL(m, alpha, &x3[0], incx);

  COMPARE_FAD_VECTORS(x1, x3, m*incx);
}

// Tests constant alpha
TYPED_TEST_P(FadBLASUnitTests, testSCAL3) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto ndot = this->ndot_;

  VectorType x1(m,ndot), x2(m,ndot), x3(m,ndot);
  for (unsigned int i=0; i<m; i++) {
    ScalarType val = this->urand.number();
    x1[i] = FadType(ndot, val);
    x2[i] = FadType(ndot, val);
    x3[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = this->urand.number();
      x1[i].fastAccessDx(k) = val;
      x2[i].fastAccessDx(k) = val;
      x3[i].fastAccessDx(k) = val;
    }
  }
  ScalarType alpha = this->urand.number();

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.SCAL(m, alpha, &x1[0], 1);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.SCAL(m, alpha, &x2[0], 1);

  COMPARE_FAD_VECTORS(x1, x2, m);

  unsigned int sz = m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.SCAL(m, alpha, &x3[0], 1);

  COMPARE_FAD_VECTORS(x1, x3, m);
}

// Tests constant x
TYPED_TEST_P(FadBLASUnitTests, testSCAL4) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto ndot = this->ndot_;

  VectorType x1(m,ndot), x2(m,ndot), x3(m,ndot);
  for (unsigned int i=0; i<m; i++) {
    ScalarType val = this->urand.number();
    x1[i] = val;
    x2[i] = val;
    x3[i] = val;
  }
  FadType alpha = FadType(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++)
    alpha.fastAccessDx(k) = this->urand.number();

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.SCAL(m, alpha, &x1[0], 1);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.SCAL(m, alpha, &x2[0], 1);

  COMPARE_FAD_VECTORS(x1, x2, m);

  unsigned int sz = m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.SCAL(m, alpha, &x3[0], 1);

  COMPARE_FAD_VECTORS(x1, x3, m);
}

// Tests all arguments
TYPED_TEST_P(FadBLASUnitTests, testCOPY1) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto ndot = this->ndot_;

  VectorType x(m,ndot), y1(m,ndot), y2(m,ndot), y3(m,ndot);
  for (unsigned int i=0; i<m; i++) {
    x[i] = FadType(ndot, this->urand.number());
    ScalarType val = this->urand.number();
    y1[i] = FadType(ndot, val);
    y2[i] = FadType(ndot, val);
    y3[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      x[i].fastAccessDx(k) = this->urand.number();
      val = this->urand.number();
      y1[i].fastAccessDx(k) = val;
      y2[i].fastAccessDx(k) = val;
      y3[i].fastAccessDx(k) = val;
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.COPY(m, &x[0], 1, &y1[0], 1);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.COPY(m, &x[0], 1, &y2[0], 1);

  COMPARE_FAD_VECTORS(y1, y2, m);

  unsigned int sz = 2*m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.COPY(m, &x[0], 1, &y3[0], 1);

  COMPARE_FAD_VECTORS(y1, y3, m);
}

// Tests non unit inc
TYPED_TEST_P(FadBLASUnitTests, testCOPY2) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto ndot = this->ndot_;

  unsigned int incx = 2;
  unsigned int incy = 3;
  VectorType x(m*incx,ndot), y1(m*incy,ndot), y2(m*incy,ndot), y3(m*incy,ndot);
  for (unsigned int i=0; i<m*incx; i++) {
    x[i] = FadType(ndot, this->urand.number());
    for (unsigned int k=0; k<ndot; k++) {
      x[i].fastAccessDx(k) = this->urand.number();
    }
  }
  for (unsigned int i=0; i<m*incy; i++) {
    ScalarType val = this->urand.number();
    y1[i] = FadType(ndot, val);
    y2[i] = FadType(ndot, val);
    y3[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = this->urand.number();
      y1[i].fastAccessDx(k) = val;
      y2[i].fastAccessDx(k) = val;
      y3[i].fastAccessDx(k) = val;
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.COPY(m, &x[0], incx, &y1[0], incy);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.COPY(m, &x[0], incx, &y2[0], incy);

  COMPARE_FAD_VECTORS(y1, y2, m*incy);

  unsigned int sz = 2*m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.COPY(m, &x[0], incx, &y3[0], incy);

  COMPARE_FAD_VECTORS(y1, y3, m*incy);
}

// Tests x constant
TYPED_TEST_P(FadBLASUnitTests, testCOPY3) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto ndot = this->ndot_;

  VectorType x(m,ndot), y1(m,ndot), y2(m,ndot), y3(m,ndot);
  for (unsigned int i=0; i<m; i++) {
    x[i] = this->urand.number();
  }
  for (unsigned int i=0; i<m; i++) {
    ScalarType val = this->urand.number();
    y1[i] = FadType(ndot, val);
    y2[i] = FadType(ndot, val);
    y3[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = this->urand.number();
      y1[i].fastAccessDx(k) = val;
      y2[i].fastAccessDx(k) = val;
      y3[i].fastAccessDx(k) = val;
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.COPY(m, &x[0], 1, &y1[0], 1);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.COPY(m, &x[0], 1, &y2[0], 1);

  COMPARE_FAD_VECTORS(y1, y2, m);

  unsigned int sz = 2*m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.COPY(m, &x[0], 1, &y3[0], 1);

  COMPARE_FAD_VECTORS(y1, y3, m);
}

// Tests y constant
TYPED_TEST_P(FadBLASUnitTests, testCOPY4) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto ndot = this->ndot_;

  VectorType x(m,ndot), y1(m,ndot), y2(m,ndot), y3(m,ndot);
  for (unsigned int i=0; i<m; i++) {
    x[i] = FadType(ndot, this->urand.number());
    ScalarType val = this->urand.number();
    y1[i] = val;
    y2[i] = val;
    y3[i] = val;
    for (unsigned int k=0; k<ndot; k++) {
      x[i].fastAccessDx(k) = this->urand.number();
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.COPY(m, &x[0], 1, &y1[0], 1);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.COPY(m, &x[0], 1, &y2[0], 1);

  COMPARE_FAD_VECTORS(y1, y2, m);

  unsigned int sz = 2*m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.COPY(m, &x[0], 1, &y3[0], 1);

  COMPARE_FAD_VECTORS(y1, y3, m);
}

// Tests all arguments
TYPED_TEST_P(FadBLASUnitTests, testAXPY1) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto ndot = this->ndot_;

  VectorType x(m,ndot), y1(m,ndot), y2(m,ndot), y3(m,ndot);
  for (unsigned int i=0; i<m; i++) {
    x[i] = FadType(ndot, this->urand.number());
    ScalarType val = this->urand.number();
    y1[i] = FadType(ndot, val);
    y2[i] = FadType(ndot, val);
    y3[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      x[i].fastAccessDx(k) = this->urand.number();
      val = this->urand.number();
      y1[i].fastAccessDx(k) = val;
      y2[i].fastAccessDx(k) = val;
      y3[i].fastAccessDx(k) = val;
    }
  }
  FadType alpha(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++)
    alpha.fastAccessDx(k) = this->urand.number();

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.AXPY(m, alpha, &x[0], 1, &y1[0], 1);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.AXPY(m, alpha, &x[0], 1, &y2[0], 1);

  COMPARE_FAD_VECTORS(y1, y2, m);

  unsigned int sz = 2*m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.AXPY(m, alpha, &x[0], 1, &y3[0], 1);

  COMPARE_FAD_VECTORS(y1, y3, m);
}

// Tests non unit inc
TYPED_TEST_P(FadBLASUnitTests, testAXPY2) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto ndot = this->ndot_;

  unsigned int incx = 2;
  unsigned int incy = 3;
  VectorType x(m*incx,ndot), y1(m*incy,ndot), y2(m*incy,ndot), y3(m*incy,ndot);
  for (unsigned int i=0; i<m*incx; i++) {
    x[i] = FadType(ndot, this->urand.number());
    for (unsigned int k=0; k<ndot; k++) {
      x[i].fastAccessDx(k) = this->urand.number();
    }
  }
  for (unsigned int i=0; i<m*incy; i++) {
    ScalarType val = this->urand.number();
    y1[i] = FadType(ndot, val);
    y2[i] = FadType(ndot, val);
    y3[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = this->urand.number();
      y1[i].fastAccessDx(k) = val;
      y2[i].fastAccessDx(k) = val;
      y3[i].fastAccessDx(k) = val;
    }
  }
  FadType alpha(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++)
    alpha.fastAccessDx(k) = this->urand.number();

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.AXPY(m, alpha, &x[0], incx, &y1[0], incy);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.AXPY(m, alpha, &x[0], incx, &y2[0], incy);

  COMPARE_FAD_VECTORS(y1, y2, m*incy);

  unsigned int sz = 2*m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.AXPY(m, alpha, &x[0], incx, &y3[0], incy);

  COMPARE_FAD_VECTORS(y1, y3, m*incy);
}

// Tests x constant
TYPED_TEST_P(FadBLASUnitTests, testAXPY3) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto ndot = this->ndot_;

  VectorType x(m,ndot), y1(m,ndot), y2(m,ndot), y3(m,ndot), y4(m,ndot);
  std::vector<ScalarType> xx(m);
  for (unsigned int i=0; i<m; i++) {
    xx[i] = this->urand.number();
    x[i] = xx[i];
    ScalarType val = this->urand.number();
    y1[i] = FadType(ndot, val);
    y2[i] = FadType(ndot, val);
    y3[i] = FadType(ndot, val);
    y4[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = this->urand.number();
      y1[i].fastAccessDx(k) = val;
      y2[i].fastAccessDx(k) = val;
      y3[i].fastAccessDx(k) = val;
      y4[i].fastAccessDx(k) = val;
    }
  }
  FadType alpha(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++)
    alpha.fastAccessDx(k) = this->urand.number();

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.AXPY(m, alpha, &x[0], 1, &y1[0], 1);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.AXPY(m, alpha, &x[0], 1, &y2[0], 1);

  COMPARE_FAD_VECTORS(y1, y2, m);

  unsigned int sz = m*(1+ndot)+m;
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.AXPY(m, alpha, &x[0], 1, &y3[0], 1);

  COMPARE_FAD_VECTORS(y1, y3, m);

  sacado_blas.AXPY(m, alpha, &xx[0], 1, &y4[0], 1);

  COMPARE_FAD_VECTORS(y1, y4, m);
}

// Tests y constant
TYPED_TEST_P(FadBLASUnitTests, testAXPY4) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto ndot = this->ndot_;

  VectorType x(m,ndot), y1(m,ndot), y2(m,ndot), y3(m,ndot);
  for (unsigned int i=0; i<m; i++) {
    x[i] = FadType(ndot, this->urand.number());
    ScalarType val = this->urand.number();
    y1[i] = val;
    y2[i] = val;
    y3[i] = val;
    for (unsigned int k=0; k<ndot; k++) {
      x[i].fastAccessDx(k) = this->urand.number();
    }
  }
  FadType alpha(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++)
    alpha.fastAccessDx(k) = this->urand.number();

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.AXPY(m, alpha, &x[0], 1, &y1[0], 1);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.AXPY(m, alpha, &x[0], 1, &y2[0], 1);

  COMPARE_FAD_VECTORS(y1, y2, m);

  unsigned int sz = 2*m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.AXPY(m, alpha, &x[0], 1, &y3[0], 1);

  COMPARE_FAD_VECTORS(y1, y3, m);
}

// Tests all arguments
TYPED_TEST_P(FadBLASUnitTests, testDOT1) {
  typedef decltype(this->fad) FadType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto ndot = this->ndot_;

  VectorType X(m,ndot), Y(m,ndot);
  for (unsigned int i=0; i<m; i++) {
    X[i] = FadType(ndot, this->real_urand.number());
    Y[i] = FadType(ndot, this->real_urand.number());
    for (unsigned int k=0; k<ndot; k++) {
      X[i].fastAccessDx(k) = this->real_urand.number();
      Y[i].fastAccessDx(k) = this->real_urand.number();
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  FadType z1 = teuchos_blas.DOT(m, &X[0], 1, &Y[0], 1);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  FadType z2 = sacado_blas.DOT(m, &X[0], 1, &Y[0], 1);

  COMPARE_FADS(z1, z2);

  unsigned int sz = 2*m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  FadType z3 = sacado_blas2.DOT(m, &X[0], 1, &Y[0], 1);

  COMPARE_FADS(z1, z3);
}

// Tests non-unit inc
TYPED_TEST_P(FadBLASUnitTests, testDOT2) {
  typedef decltype(this->fad) FadType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto ndot = this->ndot_;

  unsigned int incx = 2;
  unsigned int incy = 3;
  VectorType X(m*incx,ndot), Y(m*incy,ndot);
  for (unsigned int i=0; i<m*incx; i++) {
    X[i] = FadType(ndot, this->real_urand.number());
    for (unsigned int k=0; k<ndot; k++) {
      X[i].fastAccessDx(k) = this->real_urand.number();
    }
  }
  for (unsigned int i=0; i<m*incy; i++) {
    Y[i] = FadType(ndot, this->real_urand.number());
    for (unsigned int k=0; k<ndot; k++) {
      Y[i].fastAccessDx(k) = this->real_urand.number();
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  FadType z1 = teuchos_blas.DOT(m, &X[0], incx, &Y[0], incy);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  FadType z2 = sacado_blas.DOT(m, &X[0], incx, &Y[0], incy);

  COMPARE_FADS(z1, z2);

  unsigned int sz = 2*m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  FadType z3 = sacado_blas2.DOT(m, &X[0], incx, &Y[0], incy);

  COMPARE_FADS(z1, z3);
}

// Tests X constant
TYPED_TEST_P(FadBLASUnitTests, testDOT3) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto ndot = this->ndot_;

  VectorType X(m,0), Y(m,ndot);
  std::vector<ScalarType> x(m);
  for (unsigned int i=0; i<m; i++) {
    x[i] = this->urand.number();
    X[i] = x[i];
    Y[i] = FadType(ndot, this->real_urand.number());
    for (unsigned int k=0; k<ndot; k++) {
      Y[i].fastAccessDx(k) = this->real_urand.number();
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  FadType z1 = teuchos_blas.DOT(m, &X[0], 1, &Y[0], 1);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  FadType z2 = sacado_blas.DOT(m, &X[0], 1, &Y[0], 1);

  COMPARE_FADS(z1, z2);

  unsigned int sz = 2*m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  FadType z3 = sacado_blas2.DOT(m, &X[0], 1, &Y[0], 1);

  COMPARE_FADS(z1, z3);

  FadType z4 = sacado_blas.DOT(m, &x[0], 1, &Y[0], 1);

  COMPARE_FADS(z1, z4);
}

// Tests Y constant
TYPED_TEST_P(FadBLASUnitTests, testDOT4) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto ndot = this->ndot_;

  VectorType X(m,ndot), Y(m,0);
  std::vector<ScalarType> y(m);
  for (unsigned int i=0; i<m; i++) {
    X[i] = FadType(ndot, this->real_urand.number());
    y[i] = this->urand.number();
    Y[i] = y[i];
    for (unsigned int k=0; k<ndot; k++) {
      X[i].fastAccessDx(k) = this->real_urand.number();
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  FadType z1 = teuchos_blas.DOT(m, &X[0], 1, &Y[0], 1);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  FadType z2 = sacado_blas.DOT(m, &X[0], 1, &Y[0], 1);

  COMPARE_FADS(z1, z2);

  unsigned int sz = 2*m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  FadType z3 = sacado_blas2.DOT(m, &X[0], 1, &Y[0], 1);

  COMPARE_FADS(z1, z3);

  FadType z4 = sacado_blas.DOT(m, &X[0], 1, &y[0], 1);

  COMPARE_FADS(z1, z4);
}

// Tests all arguments
TYPED_TEST_P(FadBLASUnitTests, testNRM21) {
  typedef decltype(this->fad) FadType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto ndot = this->ndot_;

  VectorType X(m,ndot);
  for (unsigned int i=0; i<m; i++) {
    X[i] = FadType(ndot, this->real_urand.number());
    for (unsigned int k=0; k<ndot; k++) {
      X[i].fastAccessDx(k) = this->real_urand.number();
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  typename Teuchos::ScalarTraits<FadType>::magnitudeType z1 =
    teuchos_blas.NRM2(m, &X[0], 1);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  typename Teuchos::ScalarTraits<FadType>::magnitudeType z2 =
    sacado_blas.NRM2(m, &X[0], 1);

  COMPARE_FADS(z1, z2);

  unsigned int sz = m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  typename Teuchos::ScalarTraits<FadType>::magnitudeType z3 =
    sacado_blas2.NRM2(m, &X[0], 1);

  COMPARE_FADS(z1, z3);
}

// Tests non-unit inc
TYPED_TEST_P(FadBLASUnitTests, testNRM22) {
  typedef decltype(this->fad) FadType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto ndot = this->ndot_;

  unsigned int incx = 2;
  VectorType X(m*incx,ndot);
  for (unsigned int i=0; i<m*incx; i++) {
    X[i] = FadType(ndot, this->real_urand.number());
    for (unsigned int k=0; k<ndot; k++) {
      X[i].fastAccessDx(k) = this->real_urand.number();
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  typename Teuchos::ScalarTraits<FadType>::magnitudeType z1 =
    teuchos_blas.NRM2(m, &X[0], incx);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  typename Teuchos::ScalarTraits<FadType>::magnitudeType z2 =
    sacado_blas.NRM2(m, &X[0], incx);

  COMPARE_FADS(z1, z2);

  unsigned int sz = m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  typename Teuchos::ScalarTraits<FadType>::magnitudeType z3 =
    sacado_blas2.NRM2(m, &X[0], incx);

  COMPARE_FADS(z1, z3);
}

// Tests all arguments
TYPED_TEST_P(FadBLASUnitTests, testGEMV1) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;

  VectorType A(m*n,ndot), B(n,ndot), C1(m,ndot), C2(m,ndot), C3(m,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*m].fastAccessDx(k) = this->urand.number();
    }
    B[j] = FadType(ndot, this->urand.number());
    for (unsigned int k=0; k<ndot; k++)
      B[j].fastAccessDx(k) = this->urand.number();
  }
  FadType alpha(ndot, this->urand.number());
  FadType beta(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
    beta.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int i=0; i<m; i++) {
    ScalarType val = this->urand.number();
    C1[i] = FadType(ndot, val);
    C2[i] = FadType(ndot, val);
    C3[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = this->urand.number();
      C1[i].fastAccessDx(k) = val;
      C2[i].fastAccessDx(k) = val;
      C3[i].fastAccessDx(k) = val;
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1,
                    beta, &C1[0], 1);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1,
                   beta, &C2[0], 1);

  COMPARE_FAD_VECTORS(C1, C2, m);

  unsigned int sz = m*n*(1+ndot) + n*(1+ndot) + m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1,
                    beta, &C3[0], 1);

  COMPARE_FAD_VECTORS(C1, C3, m);
}

// Tests non-unit inc and different lda
TYPED_TEST_P(FadBLASUnitTests, testGEMV2) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;

  unsigned int lda = m+3;
  unsigned int incb = 2;
  unsigned int incc = 3;
  VectorType A(lda*n,ndot), B(n*incb,ndot), C1(m*incc,ndot), C2(m*incc,ndot),
    C3(m*incc,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<lda; i++) {
      A[i+j*lda] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*lda].fastAccessDx(k) = this->urand.number();
    }
  }
  for (unsigned int j=0; j<n*incb; j++) {
    B[j] = FadType(ndot, this->urand.number());
    for (unsigned int k=0; k<ndot; k++)
      B[j].fastAccessDx(k) = this->urand.number();
  }
  FadType alpha(ndot, this->urand.number());
  FadType beta(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
    beta.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int i=0; i<m*incc; i++) {
    ScalarType val = this->urand.number();
    C1[i] = FadType(ndot, val);
    C2[i] = FadType(ndot, val);
    C3[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = this->urand.number();
      C1[i].fastAccessDx(k) = val;
      C2[i].fastAccessDx(k) = val;
      C3[i].fastAccessDx(k) = val;
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], lda, &B[0], incb,
                    beta, &C1[0], incc);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], lda, &B[0], incb,
                   beta, &C2[0], incc);

  COMPARE_FAD_VECTORS(C1, C2, m*incc);

  unsigned int sz = m*n*(1+ndot) + n*(1+ndot) + m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], lda, &B[0], incb,
                    beta, &C3[0], incc);

  COMPARE_FAD_VECTORS(C1, C3, m*incc);
}

// Tests transpose with all arguments
TYPED_TEST_P(FadBLASUnitTests, testGEMV3) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;

  VectorType A(m*n,ndot), B(m,ndot), C1(n,ndot), C2(n,ndot), C3(n,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*m].fastAccessDx(k) = this->urand.number();
    }
  }
  for (unsigned int j=0; j<m; j++) {
    B[j] = FadType(ndot, this->urand.number());
    for (unsigned int k=0; k<ndot; k++)
      B[j].fastAccessDx(k) = this->urand.number();
  }
  FadType alpha(ndot, this->urand.number());
  FadType beta(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
    beta.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int i=0; i<n; i++) {
    ScalarType val = this->urand.number();
    C1[i] = FadType(ndot, val);
    C2[i] = FadType(ndot, val);
    C3[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = this->urand.number();
      C1[i].fastAccessDx(k) = val;
      C2[i].fastAccessDx(k) = val;
      C3[i].fastAccessDx(k) = val;
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMV(Teuchos::TRANS, m, n, alpha, &A[0], m, &B[0], 1,
                    beta, &C1[0], 1);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.GEMV(Teuchos::TRANS, m, n, alpha, &A[0], m, &B[0], 1,
                   beta, &C2[0], 1);

  COMPARE_FAD_VECTORS(C1, C2, n);

  unsigned int sz = m*n*(1+ndot) + n*(1+ndot) + m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMV(Teuchos::TRANS, m, n, alpha, &A[0], m, &B[0], 1,
                    beta, &C3[0], 1);

  COMPARE_FAD_VECTORS(C1, C3, n);
}

// Tests transpose with non-unit inc and different lda
TYPED_TEST_P(FadBLASUnitTests, testGEMV4) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;

  unsigned int lda = m+3;
  unsigned int incb = 2;
  unsigned int incc = 3;
  VectorType A(lda*n,ndot), B(m*incb,ndot), C1(n*incc,ndot), C2(n*incc,ndot),
    C3(n*incc,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<lda; i++) {
      A[i+j*lda] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*lda].fastAccessDx(k) = this->urand.number();
    }
  }
  for (unsigned int j=0; j<m*incb; j++) {
    B[j] = FadType(ndot, this->urand.number());
    for (unsigned int k=0; k<ndot; k++)
      B[j].fastAccessDx(k) = this->urand.number();
  }
  FadType alpha(ndot, this->urand.number());
  FadType beta(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
    beta.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int i=0; i<n*incc; i++) {
    ScalarType val = this->urand.number();
    C1[i] = FadType(ndot, val);
    C2[i] = FadType(ndot, val);
    C3[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = this->urand.number();
      C1[i].fastAccessDx(k) = val;
      C2[i].fastAccessDx(k) = val;
      C3[i].fastAccessDx(k) = val;
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMV(Teuchos::TRANS, m, n, alpha, &A[0], lda, &B[0], incb,
                    beta, &C1[0], incc);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.GEMV(Teuchos::TRANS, m, n, alpha, &A[0], lda, &B[0], incb,
                   beta, &C2[0], incc);

  COMPARE_FAD_VECTORS(C1, C2, n*incc);

  unsigned int sz = m*n*(1+ndot) + n*(1+ndot) + m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMV(Teuchos::TRANS, m, n, alpha, &A[0], lda, &B[0], incb,
                    beta, &C3[0], incc);

  COMPARE_FAD_VECTORS(C1, C3, n*incc);
}

// Tests with constant C
TYPED_TEST_P(FadBLASUnitTests, testGEMV5) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;

  VectorType A(m*n,ndot), B(n,ndot), C1(m,ndot), C2(m,ndot), C3(m,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*m].fastAccessDx(k) = this->urand.number();
    }
    B[j] = FadType(ndot, this->urand.number());
    for (unsigned int k=0; k<ndot; k++)
      B[j].fastAccessDx(k) = this->urand.number();
  }
  FadType alpha(ndot, this->urand.number());
  FadType beta(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
    beta.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int i=0; i<m; i++) {
    ScalarType val = this->urand.number();
    C1[i] = val;
    C2[i] = val;
    C3[i] = val;
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1,
                    beta, &C1[0], 1);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1,
                   beta, &C2[0], 1);

  COMPARE_FAD_VECTORS(C1, C2, m);

  unsigned int sz = m*n*(1+ndot) + n*(1+ndot) + m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1,
                    beta, &C3[0], 1);

  COMPARE_FAD_VECTORS(C1, C3, m);
}

// Tests with constant alpha, beta
TYPED_TEST_P(FadBLASUnitTests, testGEMV6) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;

  VectorType A(m*n,ndot), B(n,ndot), C1(m,ndot), C2(m,ndot), C3(m,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*m].fastAccessDx(k) = this->urand.number();
    }
    B[j] = FadType(ndot, this->urand.number());
    for (unsigned int k=0; k<ndot; k++)
      B[j].fastAccessDx(k) = this->urand.number();
  }
  ScalarType alpha = this->urand.number();
  ScalarType beta = this->urand.number();

  for (unsigned int i=0; i<m; i++) {
    ScalarType val = this->urand.number();
    C1[i] = FadType(ndot, val);
    C2[i] = FadType(ndot, val);
    C3[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = this->urand.number();
      C1[i].fastAccessDx(k) = val;
      C2[i].fastAccessDx(k) = val;
      C3[i].fastAccessDx(k) = val;
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1,
                    beta, &C1[0], 1);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1,
                   beta, &C2[0], 1);

  COMPARE_FAD_VECTORS(C1, C2, m);

  unsigned int sz = m*n*(1+ndot) + n*(1+ndot) + m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1,
                    beta, &C3[0], 1);

  COMPARE_FAD_VECTORS(C1, C3, m);
}

// Tests wth constant B
TYPED_TEST_P(FadBLASUnitTests, testGEMV7) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;

  VectorType A(m*n,ndot), B(n,0), C1(m,ndot), C2(m,ndot), C3(m,ndot),
    C4(m,ndot);
  std::vector<ScalarType> b(n);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*m].fastAccessDx(k) = this->urand.number();
    }
    b[j] = this->urand.number();
    B[j] = b[j];
  }
  FadType alpha(ndot, this->urand.number());
  FadType beta(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
    beta.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int i=0; i<m; i++) {
    ScalarType val = this->urand.number();
    C1[i] = FadType(ndot, val);
    C2[i] = FadType(ndot, val);
    C3[i] = FadType(ndot, val);
    C4[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = this->urand.number();
      C1[i].fastAccessDx(k) = val;
      C2[i].fastAccessDx(k) = val;
      C3[i].fastAccessDx(k) = val;
      C4[i].fastAccessDx(k) = val;
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1,
                    beta, &C1[0], 1);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1,
                   beta, &C2[0], 1);

  COMPARE_FAD_VECTORS(C1, C2, m);

  unsigned int sz = m*n*(1+ndot) + n + m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1,
                    beta, &C3[0], 1);

  COMPARE_FAD_VECTORS(C1, C3, m);

  sacado_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &b[0], 1,
                   beta, &C4[0], 1);

  COMPARE_FAD_VECTORS(C1, C4, m);
}

// Tests with constant A
TYPED_TEST_P(FadBLASUnitTests, testGEMV8) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;

  VectorType A(m*n,0), B(n,ndot), C1(m,ndot), C2(m,ndot), C3(m,ndot),
    C4(m,ndot);
  std::vector<ScalarType> a(m*n);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      a[i+j*m] = this->urand.number();
      A[i+j*m] = a[i+j*m];
    }
    B[j] = FadType(ndot, this->urand.number());
    for (unsigned int k=0; k<ndot; k++)
      B[j].fastAccessDx(k) = this->urand.number();
  }
  FadType alpha(ndot, this->urand.number());
  FadType beta(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
    beta.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int i=0; i<m; i++) {
    ScalarType val = this->urand.number();
    C1[i] = FadType(ndot, val);
    C2[i] = FadType(ndot, val);
    C3[i] = FadType(ndot, val);
    C4[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = this->urand.number();
      C1[i].fastAccessDx(k) = val;
      C2[i].fastAccessDx(k) = val;
      C3[i].fastAccessDx(k) = val;
      C4[i].fastAccessDx(k) = val;
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1,
                    beta, &C1[0], 1);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1,
                   beta, &C2[0], 1);

  COMPARE_FAD_VECTORS(C1, C2, m);

  unsigned int sz = m*n* + n*(1+ndot) + m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1,
                    beta, &C3[0], 1);

  COMPARE_FAD_VECTORS(C1, C3, m);

  sacado_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &a[0], m, &B[0], 1,
                   beta, &C4[0], 1);

  COMPARE_FAD_VECTORS(C1, C4, m);
}

// Tests with constant A, B
TYPED_TEST_P(FadBLASUnitTests, testGEMV9) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;

  VectorType A(m*n,0), B(n,ndot), C1(m,ndot), C2(m,ndot), C3(m,ndot),
    C4(m,ndot);
  std::vector<ScalarType> a(m*n), b(n);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      a[i+j*m] = this->urand.number();
      A[i+j*m] = a[i+j*m];
    }
    b[j] = this->urand.number();
    B[j] = b[j];
  }
  FadType alpha(ndot, this->urand.number());
  FadType beta(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
    beta.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int i=0; i<m; i++) {
    ScalarType val = this->urand.number();
    C1[i] = FadType(ndot, val);
    C2[i] = FadType(ndot, val);
    C3[i] = FadType(ndot, val);
    C4[i] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = this->urand.number();
      C1[i].fastAccessDx(k) = val;
      C2[i].fastAccessDx(k) = val;
      C3[i].fastAccessDx(k) = val;
      C4[i].fastAccessDx(k) = val;
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1,
                    beta, &C1[0], 1);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1,
                   beta, &C2[0], 1);

  COMPARE_FAD_VECTORS(C1, C2, m);

  unsigned int sz = m*n* + n*(1+ndot) + m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMV(Teuchos::NO_TRANS, m, n, alpha, &A[0], m, &B[0], 1,
                    beta, &C3[0], 1);

  COMPARE_FAD_VECTORS(C1, C3, m);

  sacado_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, &a[0], m, &b[0], 1,
                   beta, &C4[0], 1);

  COMPARE_FAD_VECTORS(C1, C4, m);
}

// Tests all arguments
TYPED_TEST_P(FadBLASUnitTests, testTRMV1) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto n = this->n_;
  auto ndot = this->ndot_;

  VectorType A(n*n,ndot), x1(n,ndot), x2(n,ndot), x3(n,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<n; i++) {
      A[i+j*n] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*n].fastAccessDx(k) = this->urand.number();
    }
    ScalarType val = this->urand.number();
    x1[j] = FadType(ndot, val);
    x2[j] = FadType(ndot, val);
    x3[j] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = this->urand.number();
      x1[j].fastAccessDx(k) = val;
      x2[j].fastAccessDx(k) = val;
      x3[j].fastAccessDx(k) = val;
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x1[0], 1);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x2[0], 1);

  COMPARE_FAD_VECTORS(x1, x2, n);

  unsigned int sz = n*n*(1+ndot) + n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testTRMV2) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto n = this->n_;
  auto ndot = this->ndot_;

  unsigned int lda = n+3;
  unsigned int incx = 2;
  VectorType A(lda*n,ndot), x1(n*incx,ndot), x2(n*incx,ndot), x3(n*incx,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<lda; i++) {
      A[i+j*lda] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*lda].fastAccessDx(k) = this->urand.number();
    }
  }
  for (unsigned int j=0; j<n*incx; j++) {
    ScalarType val = this->urand.number();
    x1[j] = FadType(ndot, val);
    x2[j] = FadType(ndot, val);
    x3[j] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = this->urand.number();
      x1[j].fastAccessDx(k) = val;
      x2[j].fastAccessDx(k) = val;
      x3[j].fastAccessDx(k) = val;
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, n, &A[0], lda, &x1[0], incx);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, n, &A[0], lda, &x2[0], incx);

  COMPARE_FAD_VECTORS(x1, x2, n*incx);

  unsigned int sz = n*n*(1+ndot) + n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testTRMV3) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto n = this->n_;
  auto ndot = this->ndot_;

  VectorType A(n*n,ndot), x1(n,ndot), x2(n,ndot), x3(n,ndot), x4(n,ndot),
    x5(n,ndot);
  std::vector<ScalarType> a(n*n);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<n; i++) {
      a[i+j*n] = this->urand.number();
      A[i+j*n] = a[i+j*n];
    }
    ScalarType val = this->urand.number();
    x1[j] = FadType(ndot, val);
    x2[j] = FadType(ndot, val);
    x3[j] = FadType(ndot, val);
    x4[j] = FadType(ndot, val);
    x5[j] = FadType(ndot, val);
    for (unsigned int k=0; k<ndot; k++) {
      val = this->urand.number();
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

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x2[0], 1);

  COMPARE_FAD_VECTORS(x1, x2, n);

  unsigned int sz = n*n+n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testTRMV4) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto n = this->n_;
  auto ndot = this->ndot_;

  VectorType A(n*n,ndot), x1(n,ndot), x2(n,ndot), x3(n,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<n; i++) {
      A[i+j*n] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*n].fastAccessDx(k) = this->urand.number();
    }
    ScalarType val = this->urand.number();
    x1[j] = val;
    x2[j] = val;
    x3[j] = val;
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x1[0], 1);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, n, &A[0], n, &x2[0], 1);

  COMPARE_FAD_VECTORS(x1, x2, n);

  unsigned int sz = n*n*(1+ndot) + n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testGER1) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;


  // GER is apparently not defined in the BLAS for complex types
  if (Teuchos::ScalarTraits<ScalarType>::isComplex)
    return;

  VectorType A1(m*n,ndot), A2(m*n,ndot), A3(m*n,ndot), x(m,ndot), y(n,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = this->urand.number();
      A1[i+j*m] = FadType(ndot, val);
      A2[i+j*m] = FadType(ndot, val);
      A3[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
        A1[i+j*m].fastAccessDx(k) = val;
        A2[i+j*m].fastAccessDx(k) = val;
        A3[i+j*m].fastAccessDx(k) = val;
      }
    }
  }
  for (unsigned int i=0; i<m; i++) {
    x[i] = FadType(ndot, this->urand.number());
    for (unsigned int k=0; k<ndot; k++)
      x[i].fastAccessDx(k) = this->urand.number();
  }
  for (unsigned int i=0; i<n; i++) {
    y[i] = FadType(ndot, this->urand.number());
    for (unsigned int k=0; k<ndot; k++)
      y[i].fastAccessDx(k) = this->urand.number();
  }
  FadType alpha(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A1[0], m);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A2[0], m);

  COMPARE_FAD_VECTORS(A1, A2, m*n);

  unsigned int sz = m*n*(1+ndot) + n*(1+ndot) + m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A3[0], m);

  COMPARE_FAD_VECTORS(A1, A3, m*n);
}

// Tests non unit inc, different lda
TYPED_TEST_P(FadBLASUnitTests, testGER2) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;


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
      ScalarType val = this->urand.number();
      A1[i+j*lda] = FadType(ndot, val);
      A2[i+j*lda] = FadType(ndot, val);
      A3[i+j*lda] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
        A1[i+j*lda].fastAccessDx(k) = val;
        A2[i+j*lda].fastAccessDx(k) = val;
        A3[i+j*lda].fastAccessDx(k) = val;
      }
    }
  }
  for (unsigned int i=0; i<m*incx; i++) {
    x[i] = FadType(ndot, this->urand.number());
    for (unsigned int k=0; k<ndot; k++)
      x[i].fastAccessDx(k) = this->urand.number();
  }
  for (unsigned int i=0; i<n*incy; i++) {
    y[i] = FadType(ndot, this->urand.number());
    for (unsigned int k=0; k<ndot; k++)
      y[i].fastAccessDx(k) = this->urand.number();
  }
  FadType alpha(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GER(m, n, alpha, &x[0], incx, &y[0], incy, &A1[0], lda);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.GER(m, n, alpha, &x[0], incx, &y[0], incy, &A2[0], lda);

  COMPARE_FAD_VECTORS(A1, A2, lda*n);

  unsigned int sz = m*n*(1+ndot) + n*(1+ndot) + m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GER(m, n, alpha, &x[0], incx, &y[0], incy, &A3[0], lda);

  COMPARE_FAD_VECTORS(A1, A3, lda*n);
}

// Tests constant alpha
TYPED_TEST_P(FadBLASUnitTests, testGER3) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;


  // GER is apparently not defined in the BLAS for complex types
  if (Teuchos::ScalarTraits<ScalarType>::isComplex)
    return;

  VectorType A1(m*n,ndot), A2(m*n,ndot), A3(m*n,ndot), x(m,ndot), y(n,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = this->urand.number();
      A1[i+j*m] = FadType(ndot, val);
      A2[i+j*m] = FadType(ndot, val);
      A3[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
        A1[i+j*m].fastAccessDx(k) = val;
        A2[i+j*m].fastAccessDx(k) = val;
        A3[i+j*m].fastAccessDx(k) = val;
      }
    }
  }
  for (unsigned int i=0; i<m; i++) {
    x[i] = FadType(ndot, this->urand.number());
    for (unsigned int k=0; k<ndot; k++)
      x[i].fastAccessDx(k) = this->urand.number();
  }
  for (unsigned int i=0; i<n; i++) {
    y[i] = FadType(ndot, this->urand.number());
    for (unsigned int k=0; k<ndot; k++)
      y[i].fastAccessDx(k) = this->urand.number();
  }
  ScalarType alpha = this->urand.number();

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A1[0], m);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A2[0], m);

  COMPARE_FAD_VECTORS(A1, A2, m*n);

  unsigned int sz = m*n*(1+ndot) + n*(1+ndot) + m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A3[0], m);

  COMPARE_FAD_VECTORS(A1, A3, m*n);
}

// Tests constant x
TYPED_TEST_P(FadBLASUnitTests, testGER4) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;


  // GER is apparently not defined in the BLAS for complex types
  if (Teuchos::ScalarTraits<ScalarType>::isComplex)
    return;

  VectorType A1(m*n,ndot), A2(m*n,ndot), A3(m*n,ndot), A4(m*n,ndot),
    A5(m*n,ndot), x(m,ndot), y(n,ndot);
  std::vector<ScalarType> xx(m);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = this->urand.number();
      A1[i+j*m] = FadType(ndot, val);
      A2[i+j*m] = FadType(ndot, val);
      A3[i+j*m] = FadType(ndot, val);
      A4[i+j*m] = FadType(ndot, val);
      A5[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
        A1[i+j*m].fastAccessDx(k) = val;
        A2[i+j*m].fastAccessDx(k) = val;
        A3[i+j*m].fastAccessDx(k) = val;
        A4[i+j*m].fastAccessDx(k) = val;
        A5[i+j*m].fastAccessDx(k) = val;
      }
    }
  }
  for (unsigned int i=0; i<m; i++) {
    xx[i] = this->urand.number();
    x[i] = xx[i];
  }
  for (unsigned int i=0; i<n; i++) {
    y[i] = FadType(ndot, this->urand.number());
    for (unsigned int k=0; k<ndot; k++)
      y[i].fastAccessDx(k) = this->urand.number();
  }
  FadType alpha(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A1[0], m);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A2[0], m);

  COMPARE_FAD_VECTORS(A1, A2, m*n);

  unsigned int sz = m*n*(1+ndot) + n*(1+ndot) + m;
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A3[0], m);

  COMPARE_FAD_VECTORS(A1, A3, m*n);

  sacado_blas.GER(m, n, alpha, &xx[0], 1, &y[0], 1, &A4[0], m);

  COMPARE_FAD_VECTORS(A1, A4, m*n);

  sacado_blas2.GER(m, n, alpha, &xx[0], 1, &y[0], 1, &A5[0], m);

  COMPARE_FAD_VECTORS(A1, A5, m*n);
}

// Tests constant y
TYPED_TEST_P(FadBLASUnitTests, testGER5) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;


  // GER is apparently not defined in the BLAS for complex types
  if (Teuchos::ScalarTraits<ScalarType>::isComplex)
    return;

  VectorType A1(m*n,ndot), A2(m*n,ndot), A3(m*n,ndot), A4(m*n,ndot),
    A5(m*n,ndot), x(m,ndot), y(n,ndot);
  std::vector<ScalarType> yy(n);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = this->urand.number();
      A1[i+j*m] = FadType(ndot, val);
      A2[i+j*m] = FadType(ndot, val);
      A3[i+j*m] = FadType(ndot, val);
      A4[i+j*m] = FadType(ndot, val);
      A5[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
        A1[i+j*m].fastAccessDx(k) = val;
        A2[i+j*m].fastAccessDx(k) = val;
        A3[i+j*m].fastAccessDx(k) = val;
        A4[i+j*m].fastAccessDx(k) = val;
        A5[i+j*m].fastAccessDx(k) = val;
      }
    }
  }
  for (unsigned int i=0; i<m; i++) {
    x[i] = FadType(ndot, this->urand.number());
    for (unsigned int k=0; k<ndot; k++)
      x[i].fastAccessDx(k) = this->urand.number();
  }
  for (unsigned int i=0; i<n; i++) {
    yy[i] = this->urand.number();
    y[i] = yy[i];
  }
  FadType alpha(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A1[0], m);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A2[0], m);

  COMPARE_FAD_VECTORS(A1, A2, m*n);

  unsigned int sz = m*n*(1+ndot) + m*(1+ndot) + n;
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A3[0], m);

  COMPARE_FAD_VECTORS(A1, A3, m*n);

  sacado_blas.GER(m, n, alpha, &x[0], 1, &yy[0], 1, &A4[0], m);

  COMPARE_FAD_VECTORS(A1, A4, m*n);

  sacado_blas2.GER(m, n, alpha, &x[0], 1, &yy[0], 1, &A5[0], m);

  COMPARE_FAD_VECTORS(A1, A5, m*n);
}

// Tests constant x and y
TYPED_TEST_P(FadBLASUnitTests, testGER6) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;


  // GER is apparently not defined in the BLAS for complex types
  if (Teuchos::ScalarTraits<ScalarType>::isComplex)
    return;

  VectorType A1(m*n,ndot), A2(m*n,ndot), A3(m*n,ndot), A4(m*n,ndot),
    A5(m*n,ndot), x(m,ndot), y(n,ndot);
  std::vector<ScalarType> xx(n), yy(n);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = this->urand.number();
      A1[i+j*m] = FadType(ndot, val);
      A2[i+j*m] = FadType(ndot, val);
      A3[i+j*m] = FadType(ndot, val);
      A4[i+j*m] = FadType(ndot, val);
      A5[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
        A1[i+j*m].fastAccessDx(k) = val;
        A2[i+j*m].fastAccessDx(k) = val;
        A3[i+j*m].fastAccessDx(k) = val;
        A4[i+j*m].fastAccessDx(k) = val;
        A5[i+j*m].fastAccessDx(k) = val;
      }
    }
  }
  for (unsigned int i=0; i<m; i++) {
    xx[i] = this->urand.number();
    x[i] = xx[i];
  }
  for (unsigned int i=0; i<n; i++) {
    yy[i] = this->urand.number();
    y[i] = yy[i];
  }
  FadType alpha(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A1[0], m);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A2[0], m);

  COMPARE_FAD_VECTORS(A1, A2, m*n);

  unsigned int sz = m*n*(1+ndot) + m + n;
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A3[0], m);

  COMPARE_FAD_VECTORS(A1, A3, m*n);

  sacado_blas.GER(m, n, alpha, &xx[0], 1, &yy[0], 1, &A4[0], m);

  COMPARE_FAD_VECTORS(A1, A4, m*n);

  sacado_blas2.GER(m, n, alpha, &xx[0], 1, &yy[0], 1, &A5[0], m);

  COMPARE_FAD_VECTORS(A1, A5, m*n);
}

// Tests constant A
TYPED_TEST_P(FadBLASUnitTests, testGER7) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;


  // GER is apparently not defined in the BLAS for complex types
  if (Teuchos::ScalarTraits<ScalarType>::isComplex)
    return;

  VectorType A1(m*n,ndot), A2(m*n,ndot), A3(m*n,ndot), x(m,ndot), y(n,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = this->urand.number();
      A1[i+j*m] = val;
      A2[i+j*m] = val;
      A3[i+j*m] = val;
    }
  }
  for (unsigned int i=0; i<m; i++) {
    x[i] = FadType(ndot, this->urand.number());
    for (unsigned int k=0; k<ndot; k++)
      x[i].fastAccessDx(k) = this->urand.number();
  }
  for (unsigned int i=0; i<n; i++) {
    y[i] = FadType(ndot, this->urand.number());
    for (unsigned int k=0; k<ndot; k++)
      y[i].fastAccessDx(k) = this->urand.number();
  }
  FadType alpha(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A1[0], m);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A2[0], m);

  COMPARE_FAD_VECTORS(A1, A2, m*n);

  unsigned int sz = m*n*(1+ndot) + n*(1+ndot) + m*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GER(m, n, alpha, &x[0], 1, &y[0], 1, &A3[0], m);

  COMPARE_FAD_VECTORS(A1, A3, m*n);
}

// Tests all arguments
TYPED_TEST_P(FadBLASUnitTests, testGEMM1) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto l = this->l_;
  auto ndot = this->ndot_;

  VectorType A(m*l,ndot), B(l*n,ndot), C1(m*n,ndot), C2(m*n,ndot), C3(m*n,ndot);
  for (unsigned int j=0; j<l; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*m].fastAccessDx(k) = this->urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<l; i++) {
      B[i+j*l] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        B[i+j*l].fastAccessDx(k) = this->urand.number();
    }
  }
  FadType alpha(ndot, this->urand.number());
  FadType beta(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
    beta.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = this->urand.number();
      C1[i+j*m] = FadType(ndot, val);
      C2[i+j*m] = FadType(ndot, val);
      C3[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
        C1[i+j*m].fastAccessDx(k) = val;
        C2[i+j*m].fastAccessDx(k) = val;
        C3[i+j*m].fastAccessDx(k) = val;
      }
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha,
                    &A[0], m, &B[0], l, beta, &C1[0], m);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha,
                    &A[0], m, &B[0], l, beta, &C2[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);

  unsigned int sz = m*l*(1+ndot) + l*n*(1+ndot) + m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testGEMM2) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto l = this->l_;
  auto ndot = this->ndot_;

  unsigned int lda = m+4;
  unsigned int ldb = l+4;
  unsigned int ldc = m+5;
  VectorType A(lda*l,ndot), B(ldb*n,ndot), C1(ldc*n,ndot), C2(ldc*n,ndot),
    C3(ldc*n,ndot);
  for (unsigned int j=0; j<l; j++) {
    for (unsigned int i=0; i<lda; i++) {
      A[i+j*lda] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*lda].fastAccessDx(k) = this->urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<ldb; i++) {
      B[i+j*ldb] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        B[i+j*ldb].fastAccessDx(k) = this->urand.number();
    }
  }
  FadType alpha(ndot, this->urand.number());
  FadType beta(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
    beta.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<ldc; i++) {
      ScalarType val = this->urand.number();
      C1[i+j*ldc] = FadType(ndot, val);
      C2[i+j*ldc] = FadType(ndot, val);
      C3[i+j*ldc] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
        C1[i+j*ldc].fastAccessDx(k) = val;
        C2[i+j*ldc].fastAccessDx(k) = val;
        C3[i+j*ldc].fastAccessDx(k) = val;
      }
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha,
                    &A[0], lda, &B[0], ldb, beta, &C1[0], ldc);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha,
                    &A[0], lda, &B[0], ldb, beta, &C2[0], ldc);

  COMPARE_FAD_VECTORS(C1, C2, ldc*n);

  unsigned int sz = m*l*(1+ndot) + l*n*(1+ndot) + m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha,
                    &A[0], lda, &B[0], ldb, beta, &C3[0], ldc);

  COMPARE_FAD_VECTORS(C1, C3, ldc*n);
}

// Tests different lda, ldb, ldc with transa
TYPED_TEST_P(FadBLASUnitTests, testGEMM3) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto l = this->l_;
  auto ndot = this->ndot_;

  unsigned int lda = l+3;
  unsigned int ldb = l+4;
  unsigned int ldc = m+5;
  VectorType A(lda*m,ndot), B(ldb*n,ndot), C1(ldc*n,ndot), C2(ldc*n,ndot),
    C3(ldc*n,ndot);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<lda; i++) {
      A[i+j*lda] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*lda].fastAccessDx(k) = this->urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<ldb; i++) {
      B[i+j*ldb] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        B[i+j*ldb].fastAccessDx(k) = this->urand.number();
    }
  }
  FadType alpha(ndot, this->urand.number());
  FadType beta(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
    beta.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<ldc; i++) {
      ScalarType val = this->urand.number();
      C1[i+j*ldc] = FadType(ndot, val);
      C2[i+j*ldc] = FadType(ndot, val);
      C3[i+j*ldc] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
        C1[i+j*ldc].fastAccessDx(k) = val;
        C2[i+j*ldc].fastAccessDx(k) = val;
        C3[i+j*ldc].fastAccessDx(k) = val;
      }
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha,
                    &A[0], lda, &B[0], ldb, beta, &C1[0], ldc);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha,
                    &A[0], lda, &B[0], ldb, beta, &C2[0], ldc);

  COMPARE_FAD_VECTORS(C1, C2, ldc*n);

  unsigned int sz = m*l*(1+ndot) + l*n*(1+ndot) + m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, l, alpha,
                    &A[0], lda, &B[0], ldb, beta, &C3[0], ldc);

  COMPARE_FAD_VECTORS(C1, C3, ldc*n);
}

// Tests different lda, ldb, ldc with transb
TYPED_TEST_P(FadBLASUnitTests, testGEMM4) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto l = this->l_;
  auto ndot = this->ndot_;

  unsigned int lda = m+4;
  unsigned int ldb = n+4;
  unsigned int ldc = m+5;
  VectorType A(lda*l,ndot), B(ldb*l,ndot), C1(ldc*n,ndot), C2(ldc*n,ndot),
    C3(ldc*n,ndot);
  for (unsigned int j=0; j<l; j++) {
    for (unsigned int i=0; i<lda; i++) {
      A[i+j*lda] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*lda].fastAccessDx(k) = this->urand.number();
    }
  }
  for (unsigned int j=0; j<l; j++) {
    for (unsigned int i=0; i<ldb; i++) {
      B[i+j*ldb] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        B[i+j*ldb].fastAccessDx(k) = this->urand.number();
    }
  }
  FadType alpha(ndot, this->urand.number());
  FadType beta(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
    beta.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<ldc; i++) {
      ScalarType val = this->urand.number();
      C1[i+j*ldc] = FadType(ndot, val);
      C2[i+j*ldc] = FadType(ndot, val);
      C3[i+j*ldc] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
        C1[i+j*ldc].fastAccessDx(k) = val;
        C2[i+j*ldc].fastAccessDx(k) = val;
        C3[i+j*ldc].fastAccessDx(k) = val;
      }
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha,
                    &A[0], lda, &B[0], ldb, beta, &C1[0], ldc);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha,
                    &A[0], lda, &B[0], ldb, beta, &C2[0], ldc);

  COMPARE_FAD_VECTORS(C1, C2, ldc*n);

  unsigned int sz = m*l*(1+ndot) + l*n*(1+ndot) + m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, m, n, l, alpha,
                    &A[0], lda, &B[0], ldb, beta, &C3[0], ldc);

  COMPARE_FAD_VECTORS(C1, C3, ldc*n);
}

// Tests different lda, ldb, ldc with transa and transb
TYPED_TEST_P(FadBLASUnitTests, testGEMM5) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto l = this->l_;
  auto ndot = this->ndot_;

  unsigned int lda = l+3;
  unsigned int ldb = n+4;
  unsigned int ldc = m+5;
  VectorType A(lda*m,ndot), B(ldb*l,ndot), C1(ldc*n,ndot), C2(ldc*n,ndot),
    C3(ldc*n,ndot);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<lda; i++) {
      A[i+j*lda] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*lda].fastAccessDx(k) = this->urand.number();
    }
  }
  for (unsigned int j=0; j<l; j++) {
    for (unsigned int i=0; i<ldb; i++) {
      B[i+j*ldb] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        B[i+j*ldb].fastAccessDx(k) = this->urand.number();
    }
  }
  FadType alpha(ndot, this->urand.number());
  FadType beta(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
    beta.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<ldc; i++) {
      ScalarType val = this->urand.number();
      C1[i+j*ldc] = FadType(ndot, val);
      C2[i+j*ldc] = FadType(ndot, val);
      C3[i+j*ldc] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
        C1[i+j*ldc].fastAccessDx(k) = val;
        C2[i+j*ldc].fastAccessDx(k) = val;
        C3[i+j*ldc].fastAccessDx(k) = val;
      }
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha,
                    &A[0], lda, &B[0], ldb, beta, &C1[0], ldc);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha,
                    &A[0], lda, &B[0], ldb, beta, &C2[0], ldc);

  COMPARE_FAD_VECTORS(C1, C2, ldc*n);

  unsigned int sz = m*l*(1+ndot) + l*n*(1+ndot) + m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
  sacado_blas2.GEMM(Teuchos::TRANS, Teuchos::TRANS, m, n, l, alpha,
                    &A[0], lda, &B[0], ldb, beta, &C3[0], ldc);

  COMPARE_FAD_VECTORS(C1, C3, ldc*n);
}

// Tests with constant C
TYPED_TEST_P(FadBLASUnitTests, testGEMM6) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto l = this->l_;
  auto ndot = this->ndot_;

  VectorType A(m*l,ndot), B(l*n,ndot), C1(m*n,ndot), C2(m*n,ndot), C3(m*n,ndot);
  for (unsigned int j=0; j<l; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*m].fastAccessDx(k) = this->urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<l; i++) {
      B[i+j*l] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        B[i+j*l].fastAccessDx(k) = this->urand.number();
    }
  }
  FadType alpha(ndot, this->urand.number());
  FadType beta(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
    beta.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = this->urand.number();
      C1[i+j*m] = val;
      C2[i+j*m] = val;
      C3[i+j*m] = val;
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha,
                    &A[0], m, &B[0], l, beta, &C1[0], m);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha,
                    &A[0], m, &B[0], l, beta, &C2[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);

  unsigned int sz = m*l*(1+ndot) + l*n*(1+ndot) + m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testGEMM7) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto l = this->l_;
  auto ndot = this->ndot_;

  VectorType A(m*l,ndot), B(l*n,ndot), C1(m*n,ndot), C2(m*n,ndot), C3(m*n,ndot);
  for (unsigned int j=0; j<l; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*m].fastAccessDx(k) = this->urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<l; i++) {
      B[i+j*l] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        B[i+j*l].fastAccessDx(k) = this->urand.number();
    }
  }
  ScalarType alpha = this->urand.number();
  ScalarType beta = this->urand.number();

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = this->urand.number();
      C1[i+j*m] = FadType(ndot, val);
      C2[i+j*m] = FadType(ndot, val);
      C3[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
        C1[i+j*m].fastAccessDx(k) = val;
        C2[i+j*m].fastAccessDx(k) = val;
        C3[i+j*m].fastAccessDx(k) = val;
      }
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha,
                    &A[0], m, &B[0], l, beta, &C1[0], m);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha,
                    &A[0], m, &B[0], l, beta, &C2[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);

  unsigned int sz = m*l*(1+ndot) + l*n*(1+ndot) + m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testGEMM8) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto l = this->l_;
  auto ndot = this->ndot_;

  VectorType A(m*l,ndot), B(l*n,ndot), C1(m*n,ndot), C2(m*n,ndot), C3(m*n,ndot),
    C4(m*n,ndot), C5(m*n,ndot);
  std::vector<ScalarType> a(m*l);
  for (unsigned int j=0; j<l; j++) {
    for (unsigned int i=0; i<m; i++) {
      a[i+j*m] = this->urand.number();
      A[i+j*m] = a[i+j*m];
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<l; i++) {
      B[i+j*l] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        B[i+j*l].fastAccessDx(k) = this->urand.number();
    }
  }
  FadType alpha(ndot, this->urand.number());
  FadType beta(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
    beta.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = this->urand.number();
      C1[i+j*m] = FadType(ndot, val);
      C2[i+j*m] = FadType(ndot, val);
      C3[i+j*m] = FadType(ndot, val);
      C4[i+j*m] = FadType(ndot, val);
      C5[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
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

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha,
                    &A[0], m, &B[0], l, beta, &C2[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);

  unsigned int sz = m*l + l*n*(1+ndot) + m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testGEMM9) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto l = this->l_;
  auto ndot = this->ndot_;

  VectorType A(m*l,ndot), B(l*n,ndot), C1(m*n,ndot), C2(m*n,ndot), C3(m*n,ndot),
    C4(m*n,ndot), C5(m*n,ndot);
  std::vector<ScalarType> b(l*n);
  for (unsigned int j=0; j<l; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*m].fastAccessDx(k) = this->urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<l; i++) {
      b[i+j*l] = this->urand.number();
      B[i+j*l] = b[i+j*l];
    }
  }
  FadType alpha(ndot, this->urand.number());
  FadType beta(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
    beta.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = this->urand.number();
      C1[i+j*m] = FadType(ndot, val);
      C2[i+j*m] = FadType(ndot, val);
      C3[i+j*m] = FadType(ndot, val);
      C4[i+j*m] = FadType(ndot, val);
      C5[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
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

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha,
                    &A[0], m, &B[0], l, beta, &C2[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);

  unsigned int sz = m*l*(1+ndot) + l*n + m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testGEMM10) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto l = this->l_;
  auto ndot = this->ndot_;

  VectorType A(m*l,ndot), B(l*n,ndot), C1(m*n,ndot), C2(m*n,ndot), C3(m*n,ndot),
    C4(m*n,ndot), C5(m*n,ndot);
  std::vector<ScalarType> a(m*l), b(l*n);
  for (unsigned int j=0; j<l; j++) {
    for (unsigned int i=0; i<m; i++) {
      a[i+j*m] = this->urand.number();
      A[i+j*m] = a[i+j*m];
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<l; i++) {
      b[i+j*l] = this->urand.number();
      B[i+j*l] = b[i+j*l];
    }
  }
  FadType alpha(ndot, this->urand.number());
  FadType beta(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
    beta.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = this->urand.number();
      C1[i+j*m] = FadType(ndot, val);
      C2[i+j*m] = FadType(ndot, val);
      C3[i+j*m] = FadType(ndot, val);
      C4[i+j*m] = FadType(ndot, val);
      C5[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
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

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, l, alpha,
                    &A[0], m, &B[0], l, beta, &C2[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);

  unsigned int sz = m*l + l*n + m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testSYMM1) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;


  // SYMM is apparently not defined in the BLAS for complex types
  if (Teuchos::ScalarTraits<ScalarType>::isComplex)
    return;

  VectorType A(m*m,ndot), B(m*n,ndot), C1(m*n,ndot), C2(m*n,ndot), C3(m*n,ndot);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*m].fastAccessDx(k) = this->urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      B[i+j*m] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        B[i+j*m].fastAccessDx(k) = this->urand.number();
    }
  }
  FadType alpha(ndot, this->urand.number());
  FadType beta(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
    beta.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = this->urand.number();
      C1[i+j*m] = FadType(ndot, val);
      C2[i+j*m] = FadType(ndot, val);
      C3[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
        C1[i+j*m].fastAccessDx(k) = val;
        C2[i+j*m].fastAccessDx(k) = val;
        C3[i+j*m].fastAccessDx(k) = val;
      }
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha,
                    &A[0], m, &B[0], m, beta, &C1[0], m);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha,
                    &A[0], m, &B[0], m, beta, &C2[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);

  unsigned int sz = m*m*(1+ndot) + 2*m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testSYMM2) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;


  // SYMM is apparently not defined in the BLAS for complex types
  if (Teuchos::ScalarTraits<ScalarType>::isComplex)
    return;

  VectorType A(n*n,ndot), B(m*n,ndot), C1(m*n,ndot), C2(m*n,ndot), C3(m*n,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<n; i++) {
      A[i+j*n] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*n].fastAccessDx(k) = this->urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      B[i+j*m] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        B[i+j*m].fastAccessDx(k) = this->urand.number();
    }
  }
  FadType alpha(ndot, this->urand.number());
  FadType beta(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
    beta.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = this->urand.number();
      C1[i+j*m] = FadType(ndot, val);
      C2[i+j*m] = FadType(ndot, val);
      C3[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
        C1[i+j*m].fastAccessDx(k) = val;
        C2[i+j*m].fastAccessDx(k) = val;
        C3[i+j*m].fastAccessDx(k) = val;
      }
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.SYMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, m, n, alpha,
                    &A[0], n, &B[0], m, beta, &C1[0], m);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.SYMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, m, n, alpha,
                    &A[0], n, &B[0], m, beta, &C2[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);

  unsigned int sz = n*n*(1+ndot) + 2*m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testSYMM3) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;


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
      A[i+j*lda] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*lda].fastAccessDx(k) = this->urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<ldb; i++) {
      B[i+j*ldb] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        B[i+j*ldb].fastAccessDx(k) = this->urand.number();
    }
  }
  FadType alpha(ndot, this->urand.number());
  FadType beta(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
    beta.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<ldc; i++) {
      ScalarType val = this->urand.number();
      C1[i+j*ldc] = FadType(ndot, val);
      C2[i+j*ldc] = FadType(ndot, val);
      C3[i+j*ldc] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
        C1[i+j*ldc].fastAccessDx(k) = val;
        C2[i+j*ldc].fastAccessDx(k) = val;
        C3[i+j*ldc].fastAccessDx(k) = val;
      }
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha,
                    &A[0], lda, &B[0], ldb, beta, &C1[0], ldc);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha,
                    &A[0], lda, &B[0], ldb, beta, &C2[0], ldc);

  COMPARE_FAD_VECTORS(C1, C2, ldc*n);

  unsigned int sz = m*m*(1+ndot) + 2*m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testSYMM4) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;


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
      A[i+j*lda] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*lda].fastAccessDx(k) = this->urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<ldb; i++) {
      B[i+j*ldb] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        B[i+j*ldb].fastAccessDx(k) = this->urand.number();
    }
  }
  FadType alpha(ndot, this->urand.number());
  FadType beta(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
    beta.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<ldc; i++) {
      ScalarType val = this->urand.number();
      C1[i+j*ldc] = FadType(ndot, val);
      C2[i+j*ldc] = FadType(ndot, val);
      C3[i+j*ldc] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
        C1[i+j*ldc].fastAccessDx(k) = val;
        C2[i+j*ldc].fastAccessDx(k) = val;
        C3[i+j*ldc].fastAccessDx(k) = val;
      }
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.SYMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, m, n, alpha,
                    &A[0], lda, &B[0], ldb, beta, &C1[0], ldc);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.SYMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, m, n, alpha,
                    &A[0], lda, &B[0], ldb, beta, &C2[0], ldc);

  COMPARE_FAD_VECTORS(C1, C2, ldc*n);

  unsigned int sz = n*n*(1+ndot) + 2*m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testSYMM5) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;


  // SYMM is apparently not defined in the BLAS for complex types
  if (Teuchos::ScalarTraits<ScalarType>::isComplex)
    return;

  VectorType A(m*m,ndot), B(m*n,ndot), C1(m*n,ndot), C2(m*n,ndot), C3(m*n,ndot);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*m].fastAccessDx(k) = this->urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      B[i+j*m] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        B[i+j*m].fastAccessDx(k) = this->urand.number();
    }
  }
  FadType alpha(ndot, this->urand.number());
  FadType beta(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
    beta.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = this->urand.number();
      C1[i+j*m] = val;
      C2[i+j*m] = val;
      C3[i+j*m] = val;
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha,
                    &A[0], m, &B[0], m, beta, &C1[0], m);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha,
                    &A[0], m, &B[0], m, beta, &C2[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);

  unsigned int sz = m*m*(1+ndot) + 2*m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testSYMM6) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;


  // SYMM is apparently not defined in the BLAS for complex types
  if (Teuchos::ScalarTraits<ScalarType>::isComplex)
    return;

  VectorType A(m*m,ndot), B(m*n,ndot), C1(m*n,ndot), C2(m*n,ndot), C3(m*n,ndot);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*m].fastAccessDx(k) = this->urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      B[i+j*m] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        B[i+j*m].fastAccessDx(k) = this->urand.number();
    }
  }
  ScalarType alpha = this->urand.number();
  ScalarType beta = this->urand.number();

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = this->urand.number();
      C1[i+j*m] = FadType(ndot, val);
      C2[i+j*m] = FadType(ndot, val);
      C3[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
        C1[i+j*m].fastAccessDx(k) = val;
        C2[i+j*m].fastAccessDx(k) = val;
        C3[i+j*m].fastAccessDx(k) = val;
      }
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha,
                    &A[0], m, &B[0], m, beta, &C1[0], m);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha,
                    &A[0], m, &B[0], m, beta, &C2[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);

  unsigned int sz = m*m*(1+ndot) + 2*m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testSYMM7) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;


  // SYMM is apparently not defined in the BLAS for complex types
  if (Teuchos::ScalarTraits<ScalarType>::isComplex)
    return;

  VectorType A(m*m,ndot), B(m*n,ndot), C1(m*n,ndot), C2(m*n,ndot), C3(m*n,ndot),
    C4(m*n,ndot), C5(m*n,ndot);
  std::vector<ScalarType> a(m*m);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<m; i++) {
      a[i+j*m] = this->urand.number();
      A[i+j*m] = a[i+j*m];
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      B[i+j*m] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        B[i+j*m].fastAccessDx(k) = this->urand.number();
    }
  }
  FadType alpha(ndot, this->urand.number());
  FadType beta(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
    beta.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = this->urand.number();
      C1[i+j*m] = FadType(ndot, val);
      C2[i+j*m] = FadType(ndot, val);
      C3[i+j*m] = FadType(ndot, val);
      C4[i+j*m] = FadType(ndot, val);
      C5[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
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

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha,
                    &A[0], m, &B[0], m, beta, &C2[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);

  unsigned int sz = m*m + 2*m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testSYMM8) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;


  // SYMM is apparently not defined in the BLAS for complex types
  if (Teuchos::ScalarTraits<ScalarType>::isComplex)
    return;

  VectorType A(m*m,ndot), B(m*n,ndot), C1(m*n,ndot), C2(m*n,ndot), C3(m*n,ndot),
    C4(m*n,ndot), C5(m*n,ndot);
  std::vector<ScalarType> b(m*n);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*m].fastAccessDx(k) = this->urand.number();
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      b[i+j*m] = this->urand.number();
      B[i+j*m] = b[i+j*m];
    }
  }
  FadType alpha(ndot, this->urand.number());
  FadType beta(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
    beta.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = this->urand.number();
      C1[i+j*m] = FadType(ndot, val);
      C2[i+j*m] = FadType(ndot, val);
      C3[i+j*m] = FadType(ndot, val);
      C4[i+j*m] = FadType(ndot, val);
      C5[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
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

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha,
                    &A[0], m, &B[0], m, beta, &C2[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);

  unsigned int sz = m*m*(1+ndot) + m*n*(2+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testSYMM9) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;


  // SYMM is apparently not defined in the BLAS for complex types
  if (Teuchos::ScalarTraits<ScalarType>::isComplex)
    return;

  VectorType A(m*m,ndot), B(m*n,ndot), C1(m*n,ndot), C2(m*n,ndot), C3(m*n,ndot),
    C4(m*n,ndot), C5(m*n,ndot);
  std::vector<ScalarType> a(m*m), b(m*n);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<m; i++) {
      a[i+j*m] = this->urand.number();
      A[i+j*m] = a[i+j*m];
    }
  }
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      b[i+j*m] = this->urand.number();
      B[i+j*m] = b[i+j*m];
    }
  }
  FadType alpha(ndot, this->urand.number());
  FadType beta(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
    beta.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = this->urand.number();
      C1[i+j*m] = FadType(ndot, val);
      C2[i+j*m] = FadType(ndot, val);
      C3[i+j*m] = FadType(ndot, val);
      C4[i+j*m] = FadType(ndot, val);
      C5[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
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

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.SYMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, m, n, alpha,
                    &A[0], m, &B[0], m, beta, &C2[0], m);

  COMPARE_FAD_VECTORS(C1, C2, m*n);

  unsigned int sz = m*m + m*n*(2+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testTRMM1) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;

  VectorType A(m*m,ndot), B1(m*n,ndot), B2(m*n,ndot), B3(m*n,ndot);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*m].fastAccessDx(k) = this->urand.number();
    }
  }
  FadType alpha(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = this->urand.number();
      B1[i+j*m] = FadType(ndot, val);
      B2[i+j*m] = FadType(ndot, val);
      B3[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
        B1[i+j*m].fastAccessDx(k) = val;
        B2[i+j*m].fastAccessDx(k) = val;
        B3[i+j*m].fastAccessDx(k) = val;
      }
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);

  COMPARE_FAD_VECTORS(B1, B2, m*n);

  unsigned int sz = m*m*(1+ndot) + m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testTRMM2) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;

  VectorType A(n*n,ndot), B1(m*n,ndot), B2(m*n,ndot), B3(m*n,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<n; i++) {
      A[i+j*n] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*n].fastAccessDx(k) = this->urand.number();
    }
  }
  FadType alpha(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = this->urand.number();
      B1[i+j*m] = FadType(ndot, val);
      B2[i+j*m] = FadType(ndot, val);
      B3[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
        B1[i+j*m].fastAccessDx(k) = val;
        B2[i+j*m].fastAccessDx(k) = val;
        B3[i+j*m].fastAccessDx(k) = val;
      }
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], n, &B1[0], m);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.TRMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], n, &B2[0], m);

  COMPARE_FAD_VECTORS(B1, B2, m*n);

  unsigned int sz = n*n*(1+ndot) + m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testTRMM3) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;

  unsigned int lda = m+4;
  unsigned int ldb = m+5;
  VectorType A(lda*m,ndot), B1(ldb*n,ndot), B2(ldb*n,ndot), B3(ldb*n,ndot);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<lda; i++) {
      A[i+j*lda] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*lda].fastAccessDx(k) = this->urand.number();
    }
  }
  FadType alpha(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<ldb; i++) {
      ScalarType val = this->urand.number();
      B1[i+j*ldb] = FadType(ndot, val);
      B2[i+j*ldb] = FadType(ndot, val);
      B3[i+j*ldb] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
        B1[i+j*ldb].fastAccessDx(k) = val;
        B2[i+j*ldb].fastAccessDx(k) = val;
        B3[i+j*ldb].fastAccessDx(k) = val;
      }
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B1[0], ldb);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B2[0], ldb);

  COMPARE_FAD_VECTORS(B1, B2, ldb*n);

  unsigned int sz = m*m*(1+ndot) + m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testTRMM4) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;

  unsigned int lda = n+4;
  unsigned int ldb = m+5;
  VectorType A(lda*n,ndot), B1(ldb*n,ndot), B2(ldb*n,ndot), B3(ldb*n,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<lda; i++) {
      A[i+j*lda] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*lda].fastAccessDx(k) = this->urand.number();
    }
  }
  FadType alpha(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<ldb; i++) {
      ScalarType val = this->urand.number();
      B1[i+j*ldb] = FadType(ndot, val);
      B2[i+j*ldb] = FadType(ndot, val);
      B3[i+j*ldb] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
        B1[i+j*ldb].fastAccessDx(k) = val;
        B2[i+j*ldb].fastAccessDx(k) = val;
        B3[i+j*ldb].fastAccessDx(k) = val;
      }
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B1[0], ldb);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.TRMM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                   Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B2[0], ldb);

  COMPARE_FAD_VECTORS(B1, B2, ldb*n);

  unsigned int sz = n*n*(1+ndot) + m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testTRMM5) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;

  VectorType A(m*m,ndot), B1(m*n,ndot), B2(m*n,ndot), B3(m*n,ndot);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*m].fastAccessDx(k) = this->urand.number();
    }
  }
  ScalarType alpha = this->urand.number();

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = this->urand.number();
      B1[i+j*m] = FadType(ndot, val);
      B2[i+j*m] = FadType(ndot, val);
      B3[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
        B1[i+j*m].fastAccessDx(k) = val;
        B2[i+j*m].fastAccessDx(k) = val;
        B3[i+j*m].fastAccessDx(k) = val;
      }
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);

  COMPARE_FAD_VECTORS(B1, B2, m*n);

  unsigned int sz = m*m*(1+ndot) + m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testTRMM6) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;

  VectorType A(m*m,ndot), B1(m*n,ndot), B2(m*n,ndot), B3(m*n,ndot);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*m].fastAccessDx(k) = this->urand.number();
    }
  }
  FadType alpha(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = this->urand.number();
      B1[i+j*m] = val;
      B2[i+j*m] = val;
      B3[i+j*m] = val;
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);

  COMPARE_FAD_VECTORS(B1, B2, m*n);

  unsigned int sz = m*m*(1+ndot) + m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testTRMM7) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;

  VectorType A(m*m,ndot), B1(m*n,ndot), B2(m*n,ndot), B3(m*n,ndot),
    B4(m*n,ndot), B5(m*n,ndot);
  std::vector<ScalarType> a(m*m);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<m; i++) {
      a[i+j*m] = this->urand.number();
      A[i+j*m] = a[i+j*m];
    }
  }
  FadType alpha(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = this->urand.number();
      B1[i+j*m] = FadType(ndot, val);
      B2[i+j*m] = FadType(ndot, val);
      B3[i+j*m] = FadType(ndot, val);
      B4[i+j*m] = FadType(ndot, val);
      B5[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
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

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);

  COMPARE_FAD_VECTORS(B1, B2, m*n);

  unsigned int sz = m*m + m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testTRSM1) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;

  VectorType A(m*m,ndot), B1(m*n,ndot), B2(m*n,ndot), B3(m*n,ndot);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<m; i++) {
      //A[i+j*m] = this->urand.number();
      A[i+j*m] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*m].fastAccessDx(k) = this->urand.number();
    }
  }
  FadType alpha(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
  }
  //ScalarType alpha = this->urand.number();

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = this->urand.number();
      // B1[i+j*m] = val;
      // B2[i+j*m] = val;
      // B3[i+j*m] = val;
      B1[i+j*m] = FadType(ndot, val);
      B2[i+j*m] = FadType(ndot, val);
      B3[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
        B1[i+j*m].fastAccessDx(k) = val;
        B2[i+j*m].fastAccessDx(k) = val;
        B3[i+j*m].fastAccessDx(k) = val;
      }
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);

  COMPARE_FAD_VECTORS(B1, B2, m*n);

  unsigned int sz = m*m*(1+ndot) + m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testTRSM2) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;

  VectorType A(n*n,ndot), B1(m*n,ndot), B2(m*n,ndot), B3(m*n,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<n; i++) {
      A[i+j*n] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*n].fastAccessDx(k) = this->urand.number();
    }
  }
  FadType alpha(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = this->urand.number();
      B1[i+j*m] = FadType(ndot, val);
      B2[i+j*m] = FadType(ndot, val);
      B3[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
        B1[i+j*m].fastAccessDx(k) = val;
        B2[i+j*m].fastAccessDx(k) = val;
        B3[i+j*m].fastAccessDx(k) = val;
      }
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRSM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], n, &B1[0], m);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.TRSM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], n, &B2[0], m);

  COMPARE_FAD_VECTORS(B1, B2, m*n);

  unsigned int sz = n*n*(1+ndot) + m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testTRSM3) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;

  unsigned int lda = m+4;
  unsigned int ldb = m+5;
  VectorType A(lda*m,ndot), B1(ldb*n,ndot), B2(ldb*n,ndot), B3(ldb*n,ndot);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<lda; i++) {
      A[i+j*lda] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*lda].fastAccessDx(k) = this->urand.number();
    }
  }
  FadType alpha(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<ldb; i++) {
      ScalarType val = this->urand.number();
      B1[i+j*ldb] = FadType(ndot, val);
      B2[i+j*ldb] = FadType(ndot, val);
      B3[i+j*ldb] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
        B1[i+j*ldb].fastAccessDx(k) = val;
        B2[i+j*ldb].fastAccessDx(k) = val;
        B3[i+j*ldb].fastAccessDx(k) = val;
      }
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B1[0], ldb);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B2[0], ldb);

  COMPARE_FAD_VECTORS(B1, B2, ldb*n);

  unsigned int sz = m*m*(1+ndot) + m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testTRSM4) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;

  unsigned int lda = n+4;
  unsigned int ldb = m+5;
  VectorType A(lda*n,ndot), B1(ldb*n,ndot), B2(ldb*n,ndot), B3(ldb*n,ndot);
  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<lda; i++) {
      A[i+j*lda] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*lda].fastAccessDx(k) = this->urand.number();
    }
  }
  FadType alpha(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<ldb; i++) {
      ScalarType val = this->urand.number();
      B1[i+j*ldb] = FadType(ndot, val);
      B2[i+j*ldb] = FadType(ndot, val);
      B3[i+j*ldb] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
        B1[i+j*ldb].fastAccessDx(k) = val;
        B2[i+j*ldb].fastAccessDx(k) = val;
        B3[i+j*ldb].fastAccessDx(k) = val;
      }
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRSM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B1[0], ldb);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.TRSM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                   Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], lda, &B2[0], ldb);

  COMPARE_FAD_VECTORS(B1, B2, ldb*n);

  unsigned int sz = n*n*(1+ndot) + m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testTRSM5) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;

  VectorType A(m*m,ndot), B1(m*n,ndot), B2(m*n,ndot), B3(m*n,ndot);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*m].fastAccessDx(k) = this->urand.number();
    }
  }
  ScalarType alpha = this->urand.number();

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = this->urand.number();
      B1[i+j*m] = FadType(ndot, val);
      B2[i+j*m] = FadType(ndot, val);
      B3[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
        B1[i+j*m].fastAccessDx(k) = val;
        B2[i+j*m].fastAccessDx(k) = val;
        B3[i+j*m].fastAccessDx(k) = val;
      }
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);

  COMPARE_FAD_VECTORS(B1, B2, m*n);

  unsigned int sz = m*m*(1+ndot) + m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testTRSM6) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;

  VectorType A(m*m,ndot), B1(m*n,ndot), B2(m*n,ndot), B3(m*n,ndot);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<m; i++) {
      A[i+j*m] = FadType(ndot, this->urand.number());
      for (unsigned int k=0; k<ndot; k++)
        A[i+j*m].fastAccessDx(k) = this->urand.number();
    }
  }
  FadType alpha(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = this->urand.number();
      B1[i+j*m] = val;
      B2[i+j*m] = val;
      B3[i+j*m] = val;
    }
  }

  Teuchos::BLAS<int,FadType> teuchos_blas;
  teuchos_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B1[0], m);

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);

  COMPARE_FAD_VECTORS(B1, B2, m*n);

  unsigned int sz = m*m*(1+ndot) + m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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
TYPED_TEST_P(FadBLASUnitTests, testTRSM7) {
  typedef decltype(this->fad) FadType;
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Sacado::Fad::Vector<unsigned int,FadType> VectorType;
  auto m = this->m_;
  auto n = this->n_;
  auto ndot = this->ndot_;

  VectorType A(m*m,ndot), B1(m*n,ndot), B2(m*n,ndot), B3(m*n,ndot),
    B4(m*n,ndot), B5(m*n,ndot);
  std::vector<ScalarType> a(m*m);
  for (unsigned int j=0; j<m; j++) {
    for (unsigned int i=0; i<m; i++) {
      a[i+j*m] = this->urand.number();
      A[i+j*m] = a[i+j*m];
    }
  }
  FadType alpha(ndot, this->urand.number());
  for (unsigned int k=0; k<ndot; k++) {
    alpha.fastAccessDx(k) = this->urand.number();
  }

  for (unsigned int j=0; j<n; j++) {
    for (unsigned int i=0; i<m; i++) {
      ScalarType val = this->urand.number();
      B1[i+j*m] = FadType(ndot, val);
      B2[i+j*m] = FadType(ndot, val);
      B3[i+j*m] = FadType(ndot, val);
      B4[i+j*m] = FadType(ndot, val);
      B5[i+j*m] = FadType(ndot, val);
      for (unsigned int k=0; k<ndot; k++) {
        val = this->urand.number();
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

  Teuchos::BLAS<int,FadType> sacado_blas(false);
  sacado_blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, m, n, alpha, &A[0], m, &B2[0], m);

  COMPARE_FAD_VECTORS(B1, B2, m*n);

  unsigned int sz = m*m + m*n*(1+ndot);
  Teuchos::BLAS<int,FadType> sacado_blas2(false,false,sz);
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

REGISTER_TYPED_TEST_SUITE_P(
  FadBLASUnitTests,
  testSCAL1,
  testSCAL2,
  testSCAL3,
  testSCAL4,

  testCOPY1,
  testCOPY2,
  testCOPY3,
  testCOPY4,

  testAXPY1,
  testAXPY2,
  testAXPY3,
  testAXPY4,

  testDOT1,
  testDOT2,
  testDOT3,
  testDOT4,

  testNRM21,
  testNRM22,

  testGEMV1,
  testGEMV2,
  testGEMV3,
  testGEMV4,
  testGEMV5,
  testGEMV6,
  testGEMV7,
  testGEMV8,
  testGEMV9,

  testTRMV1,
  testTRMV2,
  testTRMV3,
  testTRMV4,

  testGER1,
  testGER2,
  testGER3,
  testGER4,
  testGER5,
  testGER6,
  testGER7,

  testGEMM1,
  testGEMM2,
  testGEMM3,
  testGEMM4,
  testGEMM5,
  testGEMM6,
  testGEMM7,
  testGEMM8,
  testGEMM9,
  testGEMM10,

  testSYMM1,
  testSYMM2,
  testSYMM3,
  testSYMM4,
  testSYMM5,
  testSYMM6,
  testSYMM7,
  testSYMM8,
  testSYMM9,

  testTRMM1,
  testTRMM2,
  testTRMM3,
  testTRMM4,
  testTRMM5,
  testTRMM6,
  testTRMM7,

  testTRSM1,
  testTRSM2,
  testTRSM3,
  testTRSM4,
  testTRSM5,
  testTRSM6,
  testTRSM7);


#endif // FADBLASUNITTESTS_HPP
