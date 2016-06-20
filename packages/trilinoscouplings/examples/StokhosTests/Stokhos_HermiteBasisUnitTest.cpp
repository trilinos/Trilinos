// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
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

namespace HermiteBasisUnitTest {

  // Common setup for unit tests
  template <typename OrdinalType, typename ValueType>
  struct UnitTestSetup {
    ValueType rtol, atol;
    OrdinalType p;
    Stokhos::HermiteBasis<OrdinalType,ValueType> basis;
    
    UnitTestSetup() : rtol(1e-12), atol(1e-7), p(10), basis(p) {}
    
  };

  UnitTestSetup<int,double> setup;

  TEUCHOS_UNIT_TEST( Stokhos_HermiteBasis, Order ) {
    int order = setup.basis.order();
    TEUCHOS_TEST_EQUALITY(order, setup.p, out, success);
  }

  TEUCHOS_UNIT_TEST( Stokhos_HermiteBasis, Size ) {
    int size = setup.basis.size();
    TEUCHOS_TEST_EQUALITY(size, setup.p+1, out, success);
  }

  TEUCHOS_UNIT_TEST( Stokhos_HermiteBasis, Norm ) {
    const Teuchos::Array<double>& n1 = setup.basis.norm_squared();
    Teuchos::Array<double> n2(setup.p+1);
    n2[0] = 1.0;
    for (int i=1; i<=setup.p; i++)
      n2[i] = i*n2[i-1];
    success = Stokhos::compareArrays(n1, "n1", n2, "n2", setup.rtol, setup.atol,
				     out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_HermiteBasis, Norm2 ) {
    Teuchos::Array<double> n1(setup.p+1);
    Teuchos::Array<double> n2(setup.p+1);
    n1[0] = setup.basis.norm_squared(0);
    n2[0] = 1.0;
    for (int i=1; i<=setup.p; i++) {
      n1[i] = setup.basis.norm_squared(i);
      n2[i] = i*n2[i-1];
    }
    success = Stokhos::compareArrays(n1, "n1", n2, "n2", setup.rtol, setup.atol,
				     out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_HermiteBasis, QuadNorm ) {
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
  TEUCHOS_UNIT_TEST( Stokhos_HermiteBasis, QuadOrthog ) {
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
      out << "mat = " << mat << std::endl;
    }
  }

  TEUCHOS_UNIT_TEST( Stokhos_HermiteBasis, TripleProduct ) {
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

  TEUCHOS_UNIT_TEST( Stokhos_HermiteBasis, DerivDoubleProduct ) {
    Teuchos::RCP< Teuchos::SerialDenseMatrix<int, double> > Bij = 
      setup.basis.computeDerivDoubleProductTensor();
    
    Teuchos::Array<double> x, w;
    Teuchos::Array< Teuchos::Array<double> > v, val, deriv;
    setup.basis.getQuadPoints(2*setup.p, x, w, v);
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

  TEUCHOS_UNIT_TEST( Stokhos_HermiteBasis, DerivDoubleProduct2 ) {
    Teuchos::RCP< Teuchos::SerialDenseMatrix<int, double> > Bij = 
      setup.basis.computeDerivDoubleProductTensor();
    const Teuchos::Array<double>& n = setup.basis.norm_squared();

    success = true;
    for (int i=0; i<=setup.p; i++) {
      for (int j=0; j<=setup.p; j++) {
	double b = 0.0;
	if (j == i-1)
	  b = i*n[j];
	std::stringstream ss;
	ss << "Bij(" << i << "," << j << ")";
	success = success && 
	  Stokhos::compareValues((*Bij)(i,j), ss.str(), b, "b",
				 setup.rtol, setup.atol, out);
      }
    }
  }

  TEUCHOS_UNIT_TEST( Stokhos_HermiteBasis, EvaluateBases ) {
    double x = 0.1234;
    Teuchos::Array<double> v1(setup.p+1), v2(setup.p+1);
    setup.basis.evaluateBases(x, v1);

    // evaluate bases using formula
    // He_0(x) = 1
    // He_1(x) = x
    // He_i(x) = x*He_{i-1}(x) - (i-1)*He_{i-2}(x), i=2,3,...
    v2[0] = 1.0;
    if (setup.p >= 1)
      v2[1] = x;
    for (int i=2; i<=setup.p; i++)
      v2[i] = x*v2[i-1] - (i-1.0)*v2[i-2];
    success = Stokhos::compareArrays(v1, "v1", v2, "v2", setup.rtol, setup.atol,
				     out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_HermiteBasis, EvaluateBasesAndDerivatives ) {
    double x = 0.1234;
    Teuchos::Array<double> val1(setup.p+1), deriv1(setup.p+1), 
      val2(setup.p+1), deriv2(setup.p+1);
    setup.basis.evaluateBasesAndDerivatives(x, val1, deriv1);

    // evaluate bases and derivatives using formula:
    // He_0(x) = 1
    // He_1(x) = x
    // He_i(x) = x*He_{i-1}(x) - (i-1)*He_{i-2}(x), i=2,3,...
    val2[0] = 1.0;
    if (setup.p >= 1)
      val2[1] = x;
    for (int i=2; i<=setup.p; i++)
      val2[i] = x*val2[i-1] - (i-1.0)*val2[i-2];

    deriv2[0] = 0.0;
    if (setup.p >= 1)
      deriv2[1] = 1.0;
    for (int i=2; i<=setup.p; i++)
      deriv2[i] = i*val2[i-1];
    success = Stokhos::compareArrays(val1, "val1", val2, "val2", 
				     setup.rtol, setup.atol, out);
    success = success && 
      Stokhos::compareArrays(deriv1, "deriv1", deriv2, "deriv2", 
			     setup.rtol, setup.atol, out);

    
  }

  TEUCHOS_UNIT_TEST( Stokhos_HermiteBasis, Evaluate ) {
    double x = 0.1234;
    Teuchos::Array<double> v1(setup.p+1), v2(setup.p+1);
    for (int i=0; i<=setup.p; i++)
      v1[i] = setup.basis.evaluate(x, i);

    // evaluate bases using formula
    // He_0(x) = 1
    // He_1(x) = x
    // He_i(x) = x*He_{i-1}(x) - (i-1)*He_{i-2}(x), i=2,3,...
    v2[0] = 1.0;
    if (setup.p >= 1)
      v2[1] = x;
    for (int i=2; i<=setup.p; i++)
      v2[i] = x*v2[i-1] - (i-1.0)*v2[i-2];
    success = Stokhos::compareArrays(v1, "v1", v2, "v2", setup.rtol, setup.atol,
				     out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_HermiteBasis, Recurrence ) {
    Teuchos::Array<double> a1(setup.p+1), b1(setup.p+1), c1(setup.p+1), 
      d1(setup.p+1);
    Teuchos::Array<double> a2(setup.p+1), b2(setup.p+1), c2(setup.p+1), 
      d2(setup.p+1);
    setup.basis.getRecurrenceCoefficients(a1, b1, c1, d1);

    // compute coefficients using formula
    a2[0] = 0.0; b2[0] = 1.0; c2[0] = 1.0; d2[0] = 1.0;
    for (int i=1; i<=setup.p; i++) {
      a2[i] = 0.0;
      b2[i] = i;
      c2[i] = 1.0;
      d2[i] = 1.0;
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
  TEUCHOS_UNIT_TEST( Stokhos_HermiteBasis, RecurrenceAlpha ) {
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
  TEUCHOS_UNIT_TEST( Stokhos_HermiteBasis, RecurrenceBeta ) {
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
  TEUCHOS_UNIT_TEST( Stokhos_HermiteBasis, QuadPointsDakota ) {
    int n = static_cast<int>(std::ceil((setup.p+1)/2.0));
    Teuchos::Array<double> x1, w1;
    Teuchos::Array< Teuchos::Array<double> > v1;
    setup.basis.getQuadPoints(setup.p, x1, w1, v1);

    Teuchos::Array<double> x2(n), w2(n);
    Teuchos::Array< Teuchos::Array<double> > v2(n);
    webbur::hermite_compute(n, &x2[0], &w2[0]);

    for (int i=0; i<n; i++) {
      w2[i] *= 0.5/std::sqrt(std::atan(1.0)); // 1/sqrt(pi)
      x2[i] *= std::sqrt(2.0);
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
  TEUCHOS_UNIT_TEST( Stokhos_HermiteBasis, QuadPointsForUQTK ) {
    int n = static_cast<int>(std::ceil((setup.p+1)/2.0));
    Teuchos::Array<double> x1, w1;
    Teuchos::Array< Teuchos::Array<double> > v1;
    setup.basis.getQuadPoints(setup.p, x1, w1, v1);

    Teuchos::Array<double> x2(n), w2(n);
    Teuchos::Array< Teuchos::Array<double> > v2(n);
    int kind = 4;
    int kpts = 0;
    double endpts[2] = {0.0, 0.0};
    Teuchos::Array<double> b(n);
    double alpha = 0.0;
    double beta = 0.0;
    GAUSSQ_F77(&kind, &n, &alpha, &beta, &kpts, endpts, &b[0], &x2[0], &w2[0]);

    for (int i=0; i<n; i++) {
      w2[i] *= 0.5/std::sqrt(std::atan(1.0)); // 1/sqrt(pi)
      x2[i] *= std::sqrt(2.0);
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
