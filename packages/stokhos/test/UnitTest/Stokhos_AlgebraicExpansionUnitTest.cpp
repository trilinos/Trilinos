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

namespace AlgebraicExpansionUnitTest {

  // Common setup for unit tests
  template <typename OrdinalType, typename ValueType>
  struct UnitTestSetup {
    ValueType rtol, atol;
    ValueType crtol, catol;
    OrdinalType sz;
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<OrdinalType,ValueType> > basis;
    Teuchos::RCP<const Stokhos::Quadrature<OrdinalType,ValueType> > quad;
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk, Cijk_linear;
    Teuchos::RCP< Stokhos::AlgebraicOrthogPolyExpansion<OrdinalType,ValueType> > exp, exp_linear;
    Teuchos::RCP< Stokhos::QuadOrthogPolyExpansion<OrdinalType,ValueType> > qexp;
    Stokhos::OrthogPolyApprox<OrdinalType,ValueType> x, y, u, u2, cx, cu, cu2, sx, su, su2;
    ValueType a;
    
    UnitTestSetup() {
      rtol = 1e-4;
      atol = 1e-5;
      crtol = 1e-12;
      catol = 1e-12;
      a = 3.1;
      const OrdinalType d = 2;
      const OrdinalType p = 7;
      
      // Create product basis
      Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<OrdinalType,ValueType> > > bases(d);
      for (OrdinalType i=0; i<d; i++)
	bases[i] = 
	  Teuchos::rcp(new Stokhos::LegendreBasis<OrdinalType,ValueType>(p));
      basis =
	Teuchos::rcp(new Stokhos::CompletePolynomialBasis<OrdinalType,ValueType>(bases));

      // Tensor product quadrature
      quad = 
	Teuchos::rcp(new Stokhos::TensorProductQuadrature<OrdinalType,ValueType>(basis));

      // Triple product tensor
      Cijk = basis->computeTripleProductTensor();
      Cijk_linear = basis->computeLinearTripleProductTensor();
      
      // Algebraic expansion
      exp = 
	Teuchos::rcp(new Stokhos::AlgebraicOrthogPolyExpansion<OrdinalType,ValueType>(basis, Cijk));
      exp_linear = 
	Teuchos::rcp(new Stokhos::AlgebraicOrthogPolyExpansion<OrdinalType,ValueType>(basis, Cijk_linear));

      // Quadrature expansion
      qexp = 
	Teuchos::rcp(new Stokhos::QuadOrthogPolyExpansion<OrdinalType,ValueType>(basis, Cijk, quad));
      
      // Create approximation
      sz = basis->size();
      x.reset(basis);
      y.reset(basis);
      u.reset(basis); 
      u2.reset(basis);
      cx.reset(basis, 1);
      x.term(0, 0) = 1.0;
      cx.term(0, 0) = a;
      cu.reset(basis);
      cu2.reset(basis, 1);
      sx.reset(basis, d+1);
      su.reset(basis, d+1);
      su2.reset(basis, d+1);
      for (OrdinalType i=0; i<d; i++) {
	x.term(i, 1) = 0.1;
	sx.term(i, 1) = 0.0;
      }
      y.term(0, 0) = 2.0;
      for (OrdinalType i=0; i<d; i++)
	y.term(i, 1) = 0.25;
    }

    template <class Func>
    void computePCE1(Stokhos::OrthogPolyApprox<OrdinalType,ValueType>& c,
		     const Stokhos::OrthogPolyApprox<OrdinalType,ValueType>& a) 
    {
      // Quadrature data
      const Teuchos::Array<ValueType>& weights = quad->getQuadWeights();
      const Teuchos::Array< Teuchos::Array<ValueType> >& points = 
	quad->getQuadPoints();
      const Teuchos::Array< Teuchos::Array<ValueType> >& values = 
	quad->getBasisAtQuadPoints();
      OrdinalType nqp = weights.size();

      // Initialize
      for (OrdinalType i=0; i<c.size(); i++)
	c[i] = 0.0;
      
      // Compute PCE via quadrature
      Func func;
      for (OrdinalType k=0; k<nqp; k++) {
	ValueType val = a.evaluate(points[k], values[k]);
	val = func(val);
	for (int i=0; i<c.size(); i++)
	  c[i] += weights[k]*val*values[k][i] / basis->norm_squared(i);
      }
    }

    template <class Func>
    void computePCE2(Stokhos::OrthogPolyApprox<OrdinalType,ValueType>& c,
		     const Stokhos::OrthogPolyApprox<OrdinalType,ValueType>& a,
		     const Stokhos::OrthogPolyApprox<OrdinalType,ValueType>& b) 
    {
      // Quadrature data
      const Teuchos::Array<ValueType>& weights = quad->getQuadWeights();
      const Teuchos::Array< Teuchos::Array<ValueType> >& points = 
	quad->getQuadPoints();
      const Teuchos::Array< Teuchos::Array<ValueType> >& values = 
	quad->getBasisAtQuadPoints();
      OrdinalType nqp = weights.size();

      // Initialize
      for (OrdinalType i=0; i<c.size(); i++)
	c[i] = 0.0;
      
      // Compute PCE via quadrature
      Func func;
      for (OrdinalType k=0; k<nqp; k++) {
	ValueType val1 = a.evaluate(points[k], values[k]);
	ValueType val2 = b.evaluate(points[k], values[k]);
	ValueType val = func(val1, val2);
	for (int i=0; i<c.size(); i++)
	  c[i] += weights[k]*val*values[k][i] / basis->norm_squared(i);
      }
    }

    template <class Func>
    void computePCE2LC(
		  Stokhos::OrthogPolyApprox<OrdinalType,ValueType>& c,
		  ValueType a,
		  const Stokhos::OrthogPolyApprox<OrdinalType,ValueType>& b) 
    {
      // Quadrature data
      const Teuchos::Array<ValueType>& weights = quad->getQuadWeights();
      const Teuchos::Array< Teuchos::Array<ValueType> >& points = 
	quad->getQuadPoints();
      const Teuchos::Array< Teuchos::Array<ValueType> >& values = 
	quad->getBasisAtQuadPoints();
      OrdinalType nqp = weights.size();

      // Initialize
      for (OrdinalType i=0; i<c.size(); i++)
	c[i] = 0.0;
      
      // Compute PCE via quadrature
      Func func;
      for (OrdinalType k=0; k<nqp; k++) {
	ValueType val2 = b.evaluate(points[k], values[k]);
	ValueType val = func(a, val2);
	for (int i=0; i<c.size(); i++)
	  c[i] += weights[k]*val*values[k][i] / basis->norm_squared(i);
      }
    }

    template <class Func>
    void computePCE2RC(
		    Stokhos::OrthogPolyApprox<OrdinalType,ValueType>& c,
		    const Stokhos::OrthogPolyApprox<OrdinalType,ValueType>& a,
		    ValueType b) 
    {
      // Quadrature data
      const Teuchos::Array<ValueType>& weights = quad->getQuadWeights();
      const Teuchos::Array< Teuchos::Array<ValueType> >& points = 
	quad->getQuadPoints();
      const Teuchos::Array< Teuchos::Array<ValueType> >& values = 
	quad->getBasisAtQuadPoints();
      OrdinalType nqp = weights.size();

      // Initialize
      for (OrdinalType i=0; i<c.size(); i++)
	c[i] = 0.0;
      
      // Compute PCE via quadrature
      Func func;
      for (OrdinalType k=0; k<nqp; k++) {
	ValueType val1 = a.evaluate(points[k], values[k]);
	ValueType val = func(val1, b);
	for (int i=0; i<c.size(); i++)
	  c[i] += weights[k]*val*values[k][i] / basis->norm_squared(i);
      }
    }
    
  };

  UnitTestSetup<int,double> setup;

  struct UMinusFunc { 
    double operator() (double a) const { return -a; } 
  };
  
  struct ASinhFunc { 
    double operator() (double a) const { 
      return std::log(a+std::sqrt(a*a+1.0)); 
    } 
  };
  struct ACoshFunc { 
    double operator() (double a) const { 
      return std::log(a+std::sqrt(a*a-1.0)); 
    } 
  };
  struct ATanhFunc { 
    double operator() (double a) const { 
      return 0.5*std::log((1.0+a)/(1.0-a)); 
    } 
  };

  struct PlusFunc { 
    double operator() (double a, double b) const { return a + b; } 
  };
  struct MinusFunc { 
    double operator() (double a, double b) const { return a - b; } 
  };
  struct TimesFunc { 
    double operator() (double a, double b) const { return a * b; } 
  };
  struct DivideFunc { 
    double operator() (double a, double b) const { return a / b; } 
  };

  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, UMinus ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.exp->unaryMinus(setup.u, v);
    setup.computePCE1<UMinusFunc>(setup.u2, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, ExpConst ) {
    setup.exp->exp(setup.cu, setup.cx);
    setup.cu2[0] = std::exp(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, LogConst ) {
    setup.exp->log(setup.cu, setup.cx);
    setup.cu2[0] = std::log(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, Log10Const ) {
    setup.exp->log10(setup.cu, setup.cx);
    setup.cu2[0] = std::log10(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, SqrtConst ) {
    setup.exp->sqrt(setup.cu, setup.cx);
    setup.cu2[0] = std::sqrt(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, CbrtConst ) {
    setup.exp->cbrt(setup.cu, setup.cx);
    setup.cu2[0] = std::cbrt(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, SinConst ) {
    setup.exp->sin(setup.cu, setup.cx);
    setup.cu2[0] = std::sin(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, CosConst ) {
    setup.exp->cos(setup.cu, setup.cx);
    setup.cu2[0] = std::cos(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, TanConst ) {
    setup.exp->tan(setup.cu, setup.cx);
    setup.cu2[0] = std::tan(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, SinhConst ) {
    setup.exp->sinh(setup.cu, setup.cx);
    setup.cu2[0] = std::sinh(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, CoshConst ) {
    setup.exp->cosh(setup.cu, setup.cx);
    setup.cu2[0] = std::cosh(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, TanhConst ) {
    setup.exp->tanh(setup.cu, setup.cx);
    setup.cu2[0] = std::tanh(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, ASinConst ) {
    setup.exp->asin(setup.cu, setup.cx);
    setup.cu2[0] = std::asin(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, ACosConst ) {
    setup.exp->acos(setup.cu, setup.cx);
    setup.cu2[0] = std::acos(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, ATanConst ) {
    setup.exp->atan(setup.cu, setup.cx);
    setup.cu2[0] = std::atan(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, ASinhConst ) {
    ASinhFunc f;
    setup.exp->asinh(setup.cu, setup.cx);
    setup.cu2[0] = f(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, ACoshConst ) {
    ACoshFunc f;
    setup.exp->acosh(setup.cu, setup.cx);
    setup.cu2[0] = f(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, ATanhConst ) {
    ATanhFunc f;
    setup.exp->atanh(setup.cu, setup.cx);
    setup.cu2[0] = f(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, Plus ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.qexp->cos(w, setup.y);
    setup.exp->plus(setup.u, v, w);
    setup.computePCE2<PlusFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, PlusLC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.exp->plus(setup.u, setup.a, v);
    setup.computePCE2LC<PlusFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, PlusRC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.exp->plus(setup.u, v, setup.a);
    setup.computePCE2RC<PlusFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, PlusCC ) {
    setup.exp->plus(setup.cu, setup.cx, setup.cx);
    setup.computePCE2<PlusFunc>(setup.cu2, setup.cx, setup.cx);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, PlusLC2 ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.exp->plus(setup.u, setup.cx, v);
    setup.computePCE2LC<PlusFunc>(setup.u2, setup.cx[0], v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, PlusRC2 ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.exp->plus(setup.u, v, setup.cx);
    setup.computePCE2RC<PlusFunc>(setup.u2, v, setup.cx[0]);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, PlusResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.qexp->sin(v, setup.x);
    setup.qexp->cos(w, setup.y);
    setup.exp->plus(ru, v, w);
    setup.computePCE2<PlusFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, PlusLCResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.qexp->sin(v, setup.x);
    setup.exp->plus(ru, setup.a, v);
    setup.computePCE2LC<PlusFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, PlusRCResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.qexp->sin(v, setup.x);
    setup.exp->plus(ru, v, setup.a);
    setup.computePCE2RC<PlusFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, PlusLS ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.exp->plus(setup.u, setup.sx, v);
    setup.computePCE2<PlusFunc>(setup.u2, setup.sx, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, PlusRS ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.exp->plus(setup.u, v, setup.sx);
    setup.computePCE2<PlusFunc>(setup.u2, v, setup.sx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, PlusLSRC ) {
    setup.exp->plus(setup.su, setup.sx, setup.a);
    setup.computePCE2RC<PlusFunc>(setup.su2, setup.sx, setup.a);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, PlusRSLC ) {
    setup.exp->plus(setup.su, setup.a, setup.sx);
    setup.computePCE2LC<PlusFunc>(setup.su2, setup.a, setup.sx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, PlusLSRC2 ) {
    setup.exp->plus(setup.su, setup.sx, setup.cx);
    setup.computePCE2<PlusFunc>(setup.su2, setup.sx, setup.cx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, PlusRSLC2 ) {
    setup.exp->plus(setup.su, setup.cx, setup.sx);
    setup.computePCE2<PlusFunc>(setup.su2, setup.cx, setup.sx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, Minus ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.qexp->cos(w, setup.y);
    setup.exp->minus(setup.u, v, w);
    setup.computePCE2<MinusFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, MinusLC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.exp->minus(setup.u, setup.a, v);
    setup.computePCE2LC<MinusFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, MinusRC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.exp->minus(setup.u, v, setup.a);
    setup.computePCE2RC<MinusFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, MinusCC ) {
    setup.exp->minus(setup.cu, setup.cx, setup.cx);
    setup.computePCE2<MinusFunc>(setup.cu2, setup.cx, setup.cx);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, MinusLC2 ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.exp->minus(setup.u, setup.cx, v);
    setup.computePCE2LC<MinusFunc>(setup.u2, setup.cx[0], v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, MinusRC2 ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.exp->minus(setup.u, v, setup.cx);
    setup.computePCE2RC<MinusFunc>(setup.u2, v, setup.cx[0]);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, MinusResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.qexp->sin(v, setup.x);
    setup.qexp->cos(w, setup.y);
    setup.exp->minus(ru, v, w);
    setup.computePCE2<MinusFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, MinusLCResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.qexp->sin(v, setup.x);
    setup.exp->minus(ru, setup.a, v);
    setup.computePCE2LC<MinusFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, MinusRCResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.qexp->sin(v, setup.x);
    setup.exp->minus(ru, v, setup.a);
    setup.computePCE2RC<MinusFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, MinusLS ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.exp->minus(setup.u, setup.sx, v);
    setup.computePCE2<MinusFunc>(setup.u2, setup.sx, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, MinusRS ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.exp->minus(setup.u, v, setup.sx);
    setup.computePCE2<MinusFunc>(setup.u2, v, setup.sx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, MinusLSRC ) {
    setup.exp->minus(setup.su, setup.sx, setup.a);
    setup.computePCE2RC<MinusFunc>(setup.su2, setup.sx, setup.a);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, MinusRSLC ) {
    setup.exp->minus(setup.su, setup.a, setup.sx);
    setup.computePCE2LC<MinusFunc>(setup.su2, setup.a, setup.sx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, MinusLSRC2 ) {
    setup.exp->minus(setup.su, setup.sx, setup.cx);
    setup.computePCE2<MinusFunc>(setup.su2, setup.sx, setup.cx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, MinusRSLC2 ) {
    setup.exp->minus(setup.su, setup.cx, setup.sx);
    setup.computePCE2<MinusFunc>(setup.su2, setup.cx, setup.sx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, Times ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.qexp->cos(w, setup.y);
    setup.exp->times(setup.u, v, w);
    setup.computePCE2<TimesFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, TimesLC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.exp->times(setup.u, setup.a, v);
    setup.computePCE2LC<TimesFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, TimesRC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.exp->times(setup.u, v, setup.a);
    setup.computePCE2RC<TimesFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, TimesCC ) {
    setup.exp->times(setup.cu, setup.cx, setup.cx);
    setup.computePCE2<TimesFunc>(setup.cu2, setup.cx, setup.cx);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, TimesLC2 ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.exp->times(setup.u, setup.cx, v);
    setup.computePCE2LC<TimesFunc>(setup.u2, setup.cx[0], v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, TimesRC2 ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.exp->times(setup.u, v, setup.cx);
    setup.computePCE2RC<TimesFunc>(setup.u2, v, setup.cx[0]);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, TimesResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.qexp->sin(v, setup.x);
    setup.qexp->cos(w, setup.y);
    setup.exp->times(ru, v, w);
    setup.computePCE2<TimesFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, TimesLCResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.qexp->sin(v, setup.x);
    setup.exp->times(ru, setup.a, v);
    setup.computePCE2LC<TimesFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, TimesRCResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.qexp->sin(v, setup.x);
    setup.exp->times(ru, v, setup.a);
    setup.computePCE2RC<TimesFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, TimesLS ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.exp->times(setup.u, setup.sx, v);
    setup.computePCE2<TimesFunc>(setup.u2, setup.sx, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, TimesRS ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.exp->times(setup.u, v, setup.sx);
    setup.computePCE2<TimesFunc>(setup.u2, v, setup.sx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, TimesLSLinear ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.exp_linear->times(setup.u, setup.sx, v);
    setup.computePCE2<TimesFunc>(setup.u2, setup.sx, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, TimesRSLinear ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.exp_linear->times(setup.u, v, setup.sx);
    setup.computePCE2<TimesFunc>(setup.u2, v, setup.sx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, TimesLSRC ) {
    setup.exp->times(setup.su, setup.sx, setup.a);
    setup.computePCE2RC<TimesFunc>(setup.su2, setup.sx, setup.a);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, TimesRSLC ) {
    setup.exp->times(setup.su, setup.a, setup.sx);
    setup.computePCE2LC<TimesFunc>(setup.su2, setup.a, setup.sx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, TimesLSRC2 ) {
    setup.exp->times(setup.su, setup.sx, setup.cx);
    setup.computePCE2<TimesFunc>(setup.su2, setup.sx, setup.cx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, TimesRSLC2 ) {
    setup.exp->times(setup.su, setup.cx, setup.sx);
    setup.computePCE2<TimesFunc>(setup.su2, setup.cx, setup.sx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }

  
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, DivideRC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.exp->divide(setup.u, v, setup.a);
    setup.computePCE2RC<DivideFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, DivideCC ) {
    setup.exp->divide(setup.cu, setup.cx, setup.cx);
    setup.computePCE2<DivideFunc>(setup.cu2, setup.cx, setup.cx);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, DivideRC2 ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.exp->divide(setup.u, v, setup.cx);
    setup.computePCE2RC<DivideFunc>(setup.u2, v, setup.cx[0]);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, DivideRCResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.qexp->sin(v, setup.x);
    setup.exp->divide(ru, v, setup.a);
    setup.computePCE2RC<DivideFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(ru, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, DivideLSRC ) {
    setup.exp->divide(setup.su, setup.sx, setup.a);
    setup.computePCE2RC<DivideFunc>(setup.su2, setup.sx, setup.a);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, DivideLSRC2 ) {
    setup.exp->divide(setup.su, setup.sx, setup.cx);
    setup.computePCE2<DivideFunc>(setup.su2, setup.sx, setup.cx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, PowConst ) {
    setup.exp->pow(setup.cu, setup.cx, setup.cx);
    setup.cu2[0] = std::pow(setup.cx[0], setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, PlusEqual ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.qexp->cos(setup.u, setup.x);
    setup.computePCE2<PlusFunc>(setup.u2, setup.u, v);
    setup.exp->plusEqual(setup.u, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, PlusEqualC ) {
    setup.qexp->cos(setup.u, setup.x);
    setup.computePCE2RC<PlusFunc>(setup.u2, setup.u, setup.a);
    setup.exp->plusEqual(setup.u, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, PlusEqualC2 ) {
    setup.qexp->cos(setup.u, setup.x);
    setup.computePCE2RC<PlusFunc>(setup.u2, setup.u, setup.cx[0]);
    setup.exp->plusEqual(setup.u, setup.cx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, PlusEqualResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.qexp->sin(v, setup.x);
    setup.exp->plusEqual(ru, v);
    success = Stokhos::comparePCEs(ru, "ru", v, "v", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, PlusEqualS ) {
    setup.qexp->cos(setup.u, setup.x);
    setup.computePCE2<PlusFunc>(setup.u2, setup.u, setup.sx);
    setup.exp->plusEqual(setup.u, setup.sx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, PlusEqualSC ) {
    setup.su = setup.sx;
    setup.computePCE2RC<PlusFunc>(setup.su2, setup.su, setup.a);
    setup.exp->plusEqual(setup.su, setup.a);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, PlusEqualSC2 ) {
    setup.su = setup.sx;
    setup.computePCE2<PlusFunc>(setup.su2, setup.su, setup.cx);
    setup.exp->plusEqual(setup.su, setup.cx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
  				   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, MinusEqual ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.qexp->cos(setup.u, setup.x);
    setup.computePCE2<MinusFunc>(setup.u2, setup.u, v);
    setup.exp->minusEqual(setup.u, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, MinusEqualC ) {
    setup.qexp->cos(setup.u, setup.x);
    setup.computePCE2RC<MinusFunc>(setup.u2, setup.u, setup.a);
    setup.exp->minusEqual(setup.u, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, MinusEqualC2 ) {
    setup.qexp->cos(setup.u, setup.x);
    setup.computePCE2RC<MinusFunc>(setup.u2, setup.u, setup.cx[0]);
    setup.exp->minusEqual(setup.u, setup.cx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, MinusEqualResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.qexp->sin(v, setup.x);
    setup.exp->minusEqual(ru, v);
    setup.exp->unaryMinus(v, v);
    success = Stokhos::comparePCEs(ru, "ru", v, "v", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, MinusEqualS ) {
    setup.qexp->cos(setup.u, setup.x);
    setup.computePCE2<MinusFunc>(setup.u2, setup.u, setup.sx);
    setup.exp->minusEqual(setup.u, setup.sx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, MinusEqualSC ) {
    setup.su = setup.sx;
    setup.computePCE2RC<MinusFunc>(setup.su2, setup.su, setup.a);
    setup.exp->minusEqual(setup.su, setup.a);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, MinusEqualSC2 ) {
    setup.su = setup.sx;
    setup.computePCE2<MinusFunc>(setup.su2, setup.su, setup.cx);
    setup.exp->minusEqual(setup.su, setup.cx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
  				   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, TimesEqual ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.qexp->cos(setup.u, setup.x);
    setup.computePCE2<TimesFunc>(setup.u2, setup.u, v);
    setup.exp->timesEqual(setup.u, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, TimesEqualC ) {
    setup.qexp->cos(setup.u, setup.x);
    setup.computePCE2RC<TimesFunc>(setup.u2, setup.u, setup.a);
    setup.exp->timesEqual(setup.u, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, TimesEqualC2 ) {
    setup.qexp->cos(setup.u, setup.x);
    setup.computePCE2RC<TimesFunc>(setup.u2, setup.u, setup.cx[0]);
    setup.exp->timesEqual(setup.u, setup.cx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, TimesEqualResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.qexp->sin(v, setup.x);
    setup.su = setup.sx;
    setup.computePCE2<TimesFunc>(setup.u2, setup.su, v);
    setup.exp->timesEqual(setup.su, v);
    success = Stokhos::comparePCEs(setup.su, "su", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, TimesEqualS ) {
    setup.qexp->cos(setup.u, setup.x);
    setup.computePCE2<TimesFunc>(setup.u2, setup.u, setup.sx);
    setup.exp->timesEqual(setup.u, setup.sx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, TimesEqualSLinear ) {
    setup.qexp->cos(setup.u, setup.x);
    setup.computePCE2<TimesFunc>(setup.u2, setup.u, setup.sx);
    setup.exp_linear->timesEqual(setup.u, setup.sx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, TimesEqualSC ) {
    setup.su = setup.sx;
    setup.computePCE2RC<TimesFunc>(setup.su2, setup.su, setup.a);
    setup.exp->timesEqual(setup.su, setup.a);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, TimesEqualSC2 ) {
    setup.su = setup.sx;
    setup.computePCE2<TimesFunc>(setup.su2, setup.su, setup.cx);
    setup.exp->timesEqual(setup.su, setup.cx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
  				   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, DivideEqualC ) {
    setup.qexp->cos(setup.u, setup.x);
    setup.computePCE2RC<DivideFunc>(setup.u2, setup.u, setup.a);
    setup.exp->divideEqual(setup.u, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, DivideEqualC2 ) {
    setup.qexp->cos(setup.u, setup.x);
    setup.computePCE2RC<DivideFunc>(setup.u2, setup.u, setup.cx[0]);
    setup.exp->divideEqual(setup.u, setup.cx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, DivideEqualSC ) {
    setup.su = setup.sx;
    setup.computePCE2RC<DivideFunc>(setup.su2, setup.su, setup.a);
    setup.exp->divideEqual(setup.su, setup.a);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_AlgebraicExpansion, DivideEqualSC2 ) {
    setup.su = setup.sx;
    setup.computePCE2<DivideFunc>(setup.su2, setup.su, setup.cx);
    setup.exp->divideEqual(setup.su, setup.cx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
  				   setup.rtol, setup.atol, out);
  }

  // Not testing atan2(), max(), min(), abs(), fabs() since these are
  // not smooth functions

}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
