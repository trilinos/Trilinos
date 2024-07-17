// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Stokhos.hpp"
#include "Stokhos_UnitTestHelpers.hpp"
#ifdef HAVE_STOKHOS_FORUQTK
#include "Stokhos_gaussq.h"
#endif

namespace LegendreBasisUnitTest {

  // Common setup for unit tests
  template <typename OrdinalType, typename ValueType>
  struct UnitTestSetup {
    ValueType rtol, atol;
    OrdinalType p;
    Stokhos::LegendreBasis<OrdinalType,ValueType> basis;
    
    UnitTestSetup() : rtol(1e-12), atol(1e-10), p(10), basis(p,true) {}
    
  };

  UnitTestSetup<int,double> setup;

  TEUCHOS_UNIT_TEST( Stokhos_NormalizedLegendreBasis, Order ) {
    int order = setup.basis.order();
    TEUCHOS_TEST_EQUALITY(order, setup.p, out, success);
  }

  TEUCHOS_UNIT_TEST( Stokhos_NormalizedLegendreBasis, Size ) {
    int size = setup.basis.size();
    TEUCHOS_TEST_EQUALITY(size, setup.p+1, out, success);
  }

  TEUCHOS_UNIT_TEST( Stokhos_NormalizedLegendreBasis, Norm ) {
    const Teuchos::Array<double>& n1 = setup.basis.norm_squared();
    Teuchos::Array<double> n2(setup.p+1);
    for (int i=0; i<=setup.p; i++)
      n2[i] = 1.0;
    success = Stokhos::compareArrays(n1, "n1", n2, "n2", setup.rtol, setup.atol,
				     out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_NormalizedLegendreBasis, Norm2 ) {
    Teuchos::Array<double> n1(setup.p+1);
    Teuchos::Array<double> n2(setup.p+1);
    for (int i=0; i<=setup.p; i++) {
      n1[i] = setup.basis.norm_squared(i);
      n2[i] = 1.0;
    }
    success = Stokhos::compareArrays(n1, "n1", n2, "n2", setup.rtol, setup.atol,
				     out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_NormalizedLegendreBasis, QuadNorm ) {
    const Teuchos::Array<double>& n1 = setup.basis.norm_squared();
    Teuchos::Array<double> n2(setup.p+1);
    Teuchos::Array<double> x, w;
    Teuchos::Array< Teuchos::Array<double> > v;
    setup.basis.getQuadPoints(2*setup.p, x, w, v);
    int nqp = w.size();
    for (int i=0; i<=setup.p; i++) {
      n2[i] = 0;
      for (int j=0; j<nqp; j++)
	n2[i] += w[j]*v[j][i]*v[j][i];
    }
    success = Stokhos::compareArrays(n1, "n1", n2, "n2", setup.rtol, setup.atol,
				     out);
  }

  // Tests basis is orthogonal
  TEUCHOS_UNIT_TEST( Stokhos_NormalizedLegendreBasis, QuadOrthog ) {
    const Teuchos::Array<double>& norms = setup.basis.norm_squared();
    Teuchos::Array<double> x, w;
    Teuchos::Array< Teuchos::Array<double> > v;
    setup.basis.getQuadPoints(2*setup.p, x, w, v);
    int nqp = w.size();
    Teuchos::SerialDenseMatrix<int,double> mat(setup.p+1, setup.p+1);
    for (int i=0; i<=setup.p; i++) {
      for (int j=0; j<=setup.p; j++) {
	for (int k=0; k<nqp; k++)
	  mat(i,j) += w[k]*v[k][i]*v[k][j];
	mat(i,j) /= std::sqrt(norms[i]*norms[j]);
      }
      mat(i,i) -= 1.0;
    }
    success = mat.normInf() < setup.atol;
    if (!success) {
      out << "\n Error, mat.normInf() < atol = " << mat.normInf() 
	  << " < " << setup.atol << ": failed!\n";
      out << "mat = " << printMat(mat) << std::endl;
    }
  }

  TEUCHOS_UNIT_TEST( Stokhos_NormalizedLegendreBasis, TripleProduct ) {
    Teuchos::RCP< Stokhos::Dense3Tensor<int, double> > Cijk = 
      setup.basis.computeTripleProductTensor();
    
    Teuchos::Array<double> x, w;
    Teuchos::Array< Teuchos::Array<double> > v;
    setup.basis.getQuadPoints(3*setup.p, x, w, v);

    success = true;
    for (int i=0; i<=setup.p; i++) {
      for (int j=0; j<=setup.p; j++) {
  	for (int k=0; k<=setup.p; k++) {
  	  double c = 0.0;
  	  int nqp = w.size();
  	  for (int qp=0; qp<nqp; qp++)
  	    c += w[qp]*v[qp][i]*v[qp][j]*v[qp][k];
  	  std::stringstream ss;
  	  ss << "Cijk(" << i << "," << j << "," << k << ")";
  	  success = success && 
  	    Stokhos::compareValues((*Cijk)(i,j,k), ss.str(), c, "c",
  				   setup.rtol, setup.atol, out);
  	}
      }
    }
  }

  TEUCHOS_UNIT_TEST( Stokhos_NormalizedLegendreBasis, DerivDoubleProduct ) {
    Teuchos::RCP< Teuchos::SerialDenseMatrix<int, double> > Bij = 
      setup.basis.computeDerivDoubleProductTensor();
    
    Teuchos::Array<double> x, w;
    Teuchos::Array< Teuchos::Array<double> > v, val, deriv;
    setup.basis.getQuadPoints(3*setup.p, x, w, v);
    int nqp = w.size();
    val.resize(nqp);
    deriv.resize(nqp);
    for (int i=0; i<nqp; i++) {
      val[i].resize(setup.p+1);
      deriv[i].resize(setup.p+1);
      setup.basis.evaluateBasesAndDerivatives(x[i], val[i], deriv[i]);
    }

    success = true;
    for (int i=0; i<=setup.p; i++) {
      for (int j=0; j<=setup.p; j++) {
  	double b = 0.0;	
  	for (int qp=0; qp<nqp; qp++)
  	  b += w[qp]*deriv[qp][i]*val[qp][j];
  	std::stringstream ss;
  	ss << "Bij(" << i << "," << j << ")";
  	success = success && 
  	  Stokhos::compareValues((*Bij)(i,j), ss.str(), b, "b",
  				 setup.rtol, setup.atol, out);
      }
    }
  }

  // TEUCHOS_UNIT_TEST( Stokhos_NormalizedLegendreBasis, DerivDoubleProduct2 ) {
  //   Teuchos::RCP< Teuchos::SerialDenseMatrix<int, double> > Bij = 
  //     setup.basis.computeDerivDoubleProductTensor();
  //   const Teuchos::Array<double>& n = setup.basis.norm_squared();
    
  //   Teuchos::Array< Teuchos::Array<double> > deriv_coeffs(setup.p+1);
  //   deriv_coeffs[0].resize(setup.p+1,0.0);
  //   if (setup.p >= 1) {
  //     deriv_coeffs[1].resize(setup.p+1,0.0);
  //     deriv_coeffs[1][0] = 1.0;
  //   }
  //   if (setup.p >= 2) {
  //     deriv_coeffs[2].resize(setup.p+1,0.0);
  //     deriv_coeffs[2][1] = 3.0;
  //   }
  //   for (int k=3; k<=setup.p; k++) {
  //     deriv_coeffs[k].resize(setup.p+1,0.0);
  //     deriv_coeffs[k][0] = 1.0/3.0*deriv_coeffs[k-1][1];
  //     for (int i=1; i<=k-3; i++)
  // 	deriv_coeffs[k][i] = i/(2.0*i - 1.0)*deriv_coeffs[k-1][i-1] + 
  // 	  (i+1.0)/(2.0*i + 3.0)*deriv_coeffs[k-1][i+1];
  //     deriv_coeffs[k][k-2] = (k-2.0)/(2.0*k-5.0)*deriv_coeffs[k-1][k-3];
  //     deriv_coeffs[k][k-1] = (k-1.0)/(2.0*k-3.0)*deriv_coeffs[k-1][k-2] + k;
  //   }

  //   success = true;
  //   for (int i=0; i<=setup.p; i++) {
  //     for (int j=0; j<=setup.p; j++) {
  // 	double b = deriv_coeffs[i][j]*n[j];
  // 	std::stringstream ss;
  // 	ss << "Bij(" << i << "," << j << ")";
  // 	success = success && 
  // 	  Stokhos::compareValues((*Bij)(i,j), ss.str(), b, "b",
  // 				 setup.rtol, setup.atol, out);
  //     }
  //   }
  // }

  TEUCHOS_UNIT_TEST( Stokhos_NormalizedLegendreBasis, EvaluateBases ) {
    double x = 0.1234;
    Teuchos::Array<double> v1(setup.p+1), v2(setup.p+1);
    setup.basis.evaluateBases(x, v1);

    double b1, b2;
    b2 = std::sqrt(1.0/3.0);
    b1 = b2;
    v2[0] = 1.0;
    if (setup.p >= 1)
      v2[1] = x/b1;
    for (int i=2; i<=setup.p; i++) {
      b1 = std::sqrt((i*i)/((2.0*i+1.0)*(2.0*i-1.0)));
      v2[i] = (x*v2[i-1] - b2*v2[i-2]) / b1;
      b2 = b1;
    }
    success = Stokhos::compareArrays(v1, "v1", v2, "v2", setup.rtol, setup.atol,
  				     out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_NormalizedLegendreBasis, EvaluateBasesAndDerivatives ) {
    double x = 0.1234;
    Teuchos::Array<double> val1(setup.p+1), deriv1(setup.p+1), 
      val2(setup.p+1), deriv2(setup.p+1);
    setup.basis.evaluateBasesAndDerivatives(x, val1, deriv1);

    // evaluate bases and derivatives using formula:
    // P_0(x) = 1
    // P_1(x) = sqrt(3)*x
    // P_i(x) = (x*P_{i-1}(x) - b[i-1]*P_{i-2}(x))/b[i], i=2,3,...
    // where b[i] = std::sqrt((i*i)/((2.0*i+1.0)*(2.0*i-1.0)))
    // P_0'(x) = 0
    // P_1'(x) = sqrt(3)
    // P_i'(x) = (P_{i-1}(x) + x*P_{i-1}'(x) - b[i-1]*P_{i-2}')/b[i], i=2,3,...
    double b1, b2;
    b2 = std::sqrt(1.0/3.0);
    b1 = b2;
    val2[0] = 1.0;
    if (setup.p >= 1)
      val2[1] = x/b1;
    for (int i=2; i<=setup.p; i++) {
      b1 = std::sqrt((i*i)/((2.0*i+1.0)*(2.0*i-1.0)));
      val2[i] = (x*val2[i-1] - b2*val2[i-2])/b1;
      b2 = b1;
    }

    b2 = std::sqrt(1.0/3.0);
    b1 = b2;
    deriv2[0] = 0.0;
    if (setup.p >= 1)
      deriv2[1] = 1.0/b1;
    for (int i=2; i<=setup.p; i++) {
      b1 = std::sqrt((i*i)/((2.0*i+1.0)*(2.0*i-1.0)));
      deriv2[i] = (val2[i-1] + x*deriv2[i-1] - b2*deriv2[i-2])/b1;
      b2 = b1;
    }
    success = Stokhos::compareArrays(val1, "val1", val2, "val2", 
  				     setup.rtol, setup.atol, out);
    success = success && 
      Stokhos::compareArrays(deriv1, "deriv1", deriv2, "deriv2", 
  			     setup.rtol, setup.atol, out);

    
  }

  TEUCHOS_UNIT_TEST( Stokhos_NormalizedLegendreBasis, Evaluate ) {
    double x = 0.1234;
    Teuchos::Array<double> v1(setup.p+1), v2(setup.p+1);
    for (int i=0; i<=setup.p; i++)
      v1[i] = setup.basis.evaluate(x, i);

    double b1, b2;
    b2 = std::sqrt(1.0/3.0);
    b1 = b2;
    v2[0] = 1.0;
    if (setup.p >= 1)
      v2[1] = x/b1;
    for (int i=2; i<=setup.p; i++) {
      b1 = std::sqrt((i*i)/((2.0*i+1.0)*(2.0*i-1.0)));
      v2[i] = (x*v2[i-1] - b2*v2[i-2])/b1;
      b2 = b1;
    }
    success = Stokhos::compareArrays(v1, "v1", v2, "v2", setup.rtol, setup.atol,
  				     out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_NormalizedLegendreBasis, Recurrence ) {
    Teuchos::Array<double> a1(setup.p+1), b1(setup.p+1), c1(setup.p+1), 
      d1(setup.p+1);
    Teuchos::Array<double> a2(setup.p+1), b2(setup.p+1), c2(setup.p+1), 
      d2(setup.p+1);
    setup.basis.getRecurrenceCoefficients(a1, b1, c1, d1);

    // compute coefficients using formula
    a2[0] = 0.0; b2[0] = 1.0; c2[0] = 1.0; d2[0] = 1.0;
    for (int i=1; i<=setup.p; i++) {
      a2[i] = 0.0;
      b2[i] = std::sqrt((i*i)/((2.0*i+1.0)*(2.0*i-1.0))); 
      c2[i] = 1.0;
      d2[i] = b2[i];
    }
    success = true;
    success = success && 
      Stokhos::compareArrays(a1, "a1", a2, "a2", setup.rtol, setup.atol, out);
    success = success && 
      Stokhos::compareArrays(b1, "b1", b2, "b2", setup.rtol, setup.atol, out);
    success = success && 
      Stokhos::compareArrays(c1, "c1", c2, "c2", setup.rtol, setup.atol, out);
    success = success && 
      Stokhos::compareArrays(d1, "d1", d2, "d2", setup.rtol, setup.atol, out);
  }

  // Tests alpha coefficients satisfy 
  // alpha_k = delta_k * (t*psi_k,psi_k)/(psi_k,psi_k)
  TEUCHOS_UNIT_TEST( Stokhos_NormalizedLegendreBasis, RecurrenceAlpha ) {
    Teuchos::Array<double> a1(setup.p+1), b1(setup.p+1), c1(setup.p+1), 
      d1(setup.p+1);
    setup.basis.getRecurrenceCoefficients(a1, b1, c1, d1);

    Teuchos::Array<double> a2(setup.p+1);
    const Teuchos::Array<double>& n1 = setup.basis.norm_squared();
    Teuchos::Array<double> x, w;
    Teuchos::Array< Teuchos::Array<double> > v;
    setup.basis.getQuadPoints(2*setup.p, x, w, v);
    int nqp = w.size();
    for (int i=0; i<=setup.p; i++) {
      a2[i] = 0;
      for (int j=0; j<nqp; j++)
	a2[i] += w[j]*x[j]*v[j][i]*v[j][i];
      a2[i] = a2[i]*c1[i]/n1[i];
    }
    success = Stokhos::compareArrays(a1, "a1", a2, "a2", setup.rtol, setup.atol,
				     out);
  }

  // Tests beta coefficients satisfy 
  // beta_k = 
  //    gamma_k * delta_k/delta_{k-1} * (psi_k,psi_k)/(psi_{k-1},psi_{k-1})
  TEUCHOS_UNIT_TEST( Stokhos_NormalizedLegendreBasis, RecurrenceBeta ) {
    Teuchos::Array<double> a1(setup.p+1), b1(setup.p+1), c1(setup.p+1), 
      d1(setup.p+1);
    setup.basis.getRecurrenceCoefficients(a1, b1, c1, d1);

    Teuchos::Array<double> b2(setup.p+1);
    const Teuchos::Array<double>& n1 = setup.basis.norm_squared();
    b2[0] = b1[0];
    for (int i=1; i<=setup.p; i++) {
      b2[i] = d1[i]*(c1[i]/c1[i-1])*(n1[i]/n1[i-1]);
    }
    success = Stokhos::compareArrays(b1, "b1", b2, "b2", setup.rtol, setup.atol,
				     out);
  }

#ifdef HAVE_STOKHOS_DAKOTA
  TEUCHOS_UNIT_TEST( Stokhos_NormalizedLegendreBasis, QuadPointsDakota ) {
    int n = static_cast<int>(std::ceil((2*setup.p+1)/2.0));
    Teuchos::Array<double> x1, w1;
    Teuchos::Array< Teuchos::Array<double> > v1;
    setup.basis.getQuadPoints(2*setup.p, x1, w1, v1);

    Teuchos::Array<double> x2(n), w2(n);
    Teuchos::Array< Teuchos::Array<double> > v2(n);
    webbur::legendre_compute(n, &x2[0], &w2[0]);

    for (int i=0; i<n; i++) {
      w2[i] *= 0.5;  // measure = 1/2
      v2[i].resize(setup.p+1);
      setup.basis.evaluateBases(x2[i], v2[i]);
    }
    success = true;
    success = success && 
      Stokhos::compareArrays(x1, "x1", x2, "x2", setup.rtol, setup.atol, out);
    success = success && 
      Stokhos::compareArrays(w1, "w1", w2, "w2", setup.rtol, setup.atol, out);
    for (int i=0; i<n; i++) {
      std::stringstream ss1, ss2;
      ss1 << "v1[" << i << "]";
      ss2 << "v2[" << i << "]";
      success = success && 
	Stokhos::compareArrays(v1[i], ss1.str(), v2[i], ss2.str(), 
			       setup.rtol, setup.atol, out);
    }
  }
#endif

#ifdef HAVE_STOKHOS_FORUQTK
  TEUCHOS_UNIT_TEST( Stokhos_NormalizedLegendreBasis, QuadPointsForUQTK ) {
    int n = static_cast<int>(std::ceil((setup.p+1)/2.0));
    Teuchos::Array<double> x1, w1;
    Teuchos::Array< Teuchos::Array<double> > v1;
    setup.basis.getQuadPoints(setup.p, x1, w1, v1);

    Teuchos::Array<double> x2(n), w2(n);
    Teuchos::Array< Teuchos::Array<double> > v2(n);
    int kind = 1;
    int kpts = 0;
    double endpts[2] = {0.0, 0.0};
    Teuchos::Array<double> b(n);
    double alpha = 0.0;
    double beta = 0.0;
    GAUSSQ_F77(&kind, &n, &alpha, &beta, &kpts, endpts, &b[0], &x2[0], &w2[0]);

    for (int i=0; i<n; i++) {
      w2[i] *= 0.5;  // measure = 1/2
      v2[i].resize(setup.p+1);
      setup.basis.evaluateBases(x2[i], v2[i]);
    }
    success = true;
    success = success && 
      Stokhos::compareArrays(x1, "x1", x2, "x2", setup.rtol, setup.atol, out);
    success = success && 
      Stokhos::compareArrays(w1, "w1", w2, "w2", setup.rtol, setup.atol, out);
    for (int i=0; i<n; i++) {
      std::stringstream ss1, ss2;
      ss1 << "v1[" << i << "]";
      ss2 << "v2[" << i << "]";
      success = success && 
	Stokhos::compareArrays(v1[i], ss1.str(), v2[i], ss2.str(), 
			       setup.rtol, setup.atol, out);
    }
  }
#endif

}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
