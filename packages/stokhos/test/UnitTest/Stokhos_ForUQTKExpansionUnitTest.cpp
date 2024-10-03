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

namespace ForUQTKExpansionUnitTest {

  // Common setup for unit tests
  template <typename OrdinalType, typename ValueType>
  struct UnitTestSetup {
    ValueType rtol, atol;
    ValueType crtol, catol;
    OrdinalType sz;
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<OrdinalType,ValueType> > basis;
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk;
    Teuchos::RCP<const Stokhos::Quadrature<OrdinalType,ValueType> > quad;
    Teuchos::RCP< Stokhos::ForUQTKOrthogPolyExpansion<OrdinalType,ValueType> > exp;
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

      // Triple product tensor
      Cijk = basis->computeTripleProductTensor();
      
      // Tensor product quadrature
      quad = 
	Teuchos::rcp(new Stokhos::TensorProductQuadrature<OrdinalType,ValueType>(basis));
      
      // Quadrature expansion
      exp = 
	Teuchos::rcp(new Stokhos::ForUQTKOrthogPolyExpansion<OrdinalType,ValueType>(basis, Cijk, Stokhos::ForUQTKOrthogPolyExpansion<OrdinalType,ValueType>::INTEGRATION));
      
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
  struct ExpFunc { 
    double operator() (double a) const { return std::exp(a); } 
  };
  struct LogFunc { 
    double operator() (double a) const { return std::log(a); } 
  };
  struct Log10Func { 
    double operator() (double a) const { return std::log10(a); } 
  };
  struct SqrtFunc { 
    double operator() (double a) const { return std::sqrt(a); } 
  };
  struct CbrtFunc { 
    double operator() (double a) const { return std::sqrt(a); } 
  };
  struct SinFunc { 
    double operator() (double a) const { return std::sin(a); } 
  };
  struct CosFunc { 
    double operator() (double a) const { return std::cos(a); } 
  };
  struct TanFunc { 
    double operator() (double a) const { return std::tan(a); } 
  };
  struct SinhFunc { 
    double operator() (double a) const { return std::sinh(a); } 
  };
  struct CoshFunc { 
    double operator() (double a) const { return std::cosh(a); } 
  };
  struct TanhFunc { 
    double operator() (double a) const { return std::tanh(a); } 
  };
  struct ASinFunc { 
    double operator() (double a) const { return std::asin(a); } 
  };
  struct ACosFunc { 
    double operator() (double a) const { return std::acos(a); } 
  };
  struct ATanFunc { 
    double operator() (double a) const { return std::atan(a); } 
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
  struct PowFunc { 
    double operator() (double a, double b) const { return std::pow(a,b); } 
  };

  /*
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, UMinus ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->unaryMinus(setup.u, v);
    setup.computePCE1<UMinusFunc>(setup.u2, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, Exp ) {
    setup.exp->exp(setup.u, setup.x);
    setup.computePCE1<ExpFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, Log ) {
    setup.exp->log(setup.u, setup.x);
    setup.computePCE1<LogFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, Log10 ) {
    setup.exp->log10(setup.u, setup.x);
    setup.computePCE1<Log10Func>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, Sqrt ) {
    setup.exp->sqrt(setup.u, setup.x);
    setup.computePCE1<SqrtFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, Cbrt ) {
    setup.exp->cbrt(setup.u, setup.x);
    setup.computePCE1<CbrtFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, Sinh ) {
    setup.exp->sinh(setup.u, setup.x);
    setup.computePCE1<SinhFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, Cosh ) {
    setup.exp->cosh(setup.u, setup.x);
    setup.computePCE1<CoshFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, Tanh ) {
    setup.exp->tanh(setup.u, setup.x);
    setup.computePCE1<TanhFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, Plus ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->cosh(w, setup.y);
    setup.exp->plus(setup.u, v, w);
    setup.computePCE2<PlusFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, PlusLC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->plus(setup.u, setup.a, v);
    setup.computePCE2LC<PlusFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, PlusRC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->plus(setup.u, v, setup.a);
    setup.computePCE2RC<PlusFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, Minus ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->cosh(w, setup.y);
    setup.exp->minus(setup.u, v, w);
    setup.computePCE2<MinusFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, MinusLC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->minus(setup.u, setup.a, v);
    setup.computePCE2LC<MinusFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, MinusRC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->minus(setup.u, v, setup.a);
    setup.computePCE2RC<MinusFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, Times ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->cosh(w, setup.y);
    setup.exp->times(setup.u, v, w);
    setup.computePCE2<TimesFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, TimesLC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->times(setup.u, setup.a, v);
    setup.computePCE2LC<TimesFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, TimesRC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->times(setup.u, v, setup.a);
    setup.computePCE2RC<TimesFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, Divide ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->cosh(w, setup.y);
    setup.exp->divide(setup.u, v, w);
    setup.computePCE2<DivideFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, DivideLC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->divide(setup.u, setup.a, v);
    setup.computePCE2LC<DivideFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, DivideRC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->divide(setup.u, v, setup.a);
    setup.computePCE2RC<DivideFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, Pow ) {
    setup.exp->pow(setup.u, setup.x, setup.y);
    setup.computePCE2<PowFunc>(setup.u2, setup.x, setup.y);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, PowLC ) {
    setup.exp->pow(setup.u, setup.a, setup.y);
    setup.computePCE2LC<PowFunc>(setup.u2, setup.a, setup.y);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, PowRC ) {
    setup.exp->pow(setup.u, setup.x, setup.a);
    setup.computePCE2RC<PowFunc>(setup.u2, setup.x, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, PlusEqual ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->cosh(setup.u, setup.x);
    setup.computePCE2<PlusFunc>(setup.u2, setup.u, v);
    setup.exp->plusEqual(setup.u, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, PlusEqualC ) {
    setup.exp->cosh(setup.u, setup.x);
    setup.computePCE2RC<PlusFunc>(setup.u2, setup.u, setup.a);
    setup.exp->plusEqual(setup.u, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, MinusEqual ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->cosh(setup.u, setup.x);
    setup.computePCE2<MinusFunc>(setup.u2, setup.u, v);
    setup.exp->minusEqual(setup.u, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, MinusEqualC ) {
    setup.exp->cosh(setup.u, setup.x);
    setup.computePCE2RC<MinusFunc>(setup.u2, setup.u, setup.a);
    setup.exp->minusEqual(setup.u, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, TimesEqual ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->cosh(setup.u, setup.x);
    setup.computePCE2<TimesFunc>(setup.u2, setup.u, v);
    setup.exp->timesEqual(setup.u, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, TimesEqualC ) {
    setup.exp->cosh(setup.u, setup.x);
    setup.computePCE2RC<TimesFunc>(setup.u2, setup.u, setup.a);
    setup.exp->timesEqual(setup.u, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, DivideEqual ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->cosh(setup.u, setup.x);
    setup.computePCE2<DivideFunc>(setup.u2, setup.u, v);
    setup.exp->divideEqual(setup.u, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, DivideEqualC ) {
    setup.exp->cosh(setup.u, setup.x);
    setup.computePCE2RC<DivideFunc>(setup.u2, setup.u, setup.a);
    setup.exp->divideEqual(setup.u, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  */

  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, UMinus ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->unaryMinus(setup.u, v);
    setup.computePCE1<UMinusFunc>(setup.u2, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, Exp ) {
    setup.exp->exp(setup.u, setup.x);
    setup.computePCE1<ExpFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, ExpConst ) {
    setup.exp->exp(setup.cu, setup.cx);
    setup.cu2[0] = std::exp(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, ExpResize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->exp(ru, setup.x);
    setup.computePCE1<ExpFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, Log ) {
    setup.exp->log(setup.u, setup.x);
    setup.computePCE1<LogFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, LogConst ) {
    setup.exp->log(setup.cu, setup.cx);
    setup.cu2[0] = std::log(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, LogResize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->log(ru, setup.x);
    setup.computePCE1<LogFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, Log10 ) {
    setup.exp->log10(setup.u, setup.x);
    setup.computePCE1<Log10Func>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, Log10Const ) {
    setup.exp->log10(setup.cu, setup.cx);
    setup.cu2[0] = std::log10(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, Log10Resize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->log10(ru, setup.x);
    setup.computePCE1<Log10Func>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, Sqrt ) {
    setup.exp->sqrt(setup.u, setup.x);
    setup.computePCE1<SqrtFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, SqrtConst ) {
    setup.exp->sqrt(setup.cu, setup.cx);
    setup.cu2[0] = std::sqrt(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, SqrtResize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sqrt(ru, setup.x);
    setup.computePCE1<SqrtFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
   TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, Cbrt ) {
    setup.exp->cbrt(setup.u, setup.x);
    setup.computePCE1<CbrtFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, CbrtConst ) {
    setup.exp->cbrt(setup.cu, setup.cx);
    setup.cu2[0] = std::cbrt(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, CbrtResize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->cbrt(ru, setup.x);
    setup.computePCE1<CbrtFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, SinConst ) {
    setup.exp->sin(setup.cu, setup.cx);
    setup.cu2[0] = std::sin(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, CosConst ) {
    setup.exp->cos(setup.cu, setup.cx);
    setup.cu2[0] = std::cos(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, TanConst ) {
    setup.exp->tan(setup.cu, setup.cx);
    setup.cu2[0] = std::tan(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, Sinh ) {
    setup.exp->sinh(setup.u, setup.x);
    setup.computePCE1<SinhFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, SinhConst ) {
    setup.exp->sinh(setup.cu, setup.cx);
    setup.cu2[0] = std::sinh(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, SinhResize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sinh(ru, setup.x);
    setup.computePCE1<SinhFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, Cosh ) {
    setup.exp->cosh(setup.u, setup.x);
    setup.computePCE1<CoshFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, CoshConst ) {
    setup.exp->cosh(setup.cu, setup.cx);
    setup.cu2[0] = std::cosh(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, CoshResize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->cosh(ru, setup.x);
    setup.computePCE1<CoshFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, Tanh ) {
    setup.exp->tanh(setup.u, setup.x);
    setup.computePCE1<TanhFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, TanhConst ) {
    setup.exp->tanh(setup.cu, setup.cx);
    setup.cu2[0] = std::tanh(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, TanhResize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->tanh(ru, setup.x);
    setup.computePCE1<TanhFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, ASinConst ) {
    setup.exp->asin(setup.cu, setup.cx);
    setup.cu2[0] = std::asin(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, ACosConst ) {
    setup.exp->acos(setup.cu, setup.cx);
    setup.cu2[0] = std::acos(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, ATanConst ) {
    setup.exp->atan(setup.cu, setup.cx);
    setup.cu2[0] = std::atan(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, ASinhConst ) {
    ASinhFunc f;
    setup.exp->asinh(setup.cu, setup.cx);
    setup.cu2[0] = f(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, ACoshConst ) {
    ACoshFunc f;
    setup.exp->acosh(setup.cu, setup.cx);
    setup.cu2[0] = f(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, ATanhConst ) {
    ATanhFunc f;
    setup.exp->atanh(setup.cu, setup.cx);
    setup.cu2[0] = f(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, Plus ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->cosh(w, setup.y);
    setup.exp->plus(setup.u, v, w);
    setup.computePCE2<PlusFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, PlusLC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->plus(setup.u, setup.a, v);
    setup.computePCE2LC<PlusFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, PlusRC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->plus(setup.u, v, setup.a);
    setup.computePCE2RC<PlusFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, PlusCC ) {
    setup.exp->plus(setup.cu, setup.cx, setup.cx);
    setup.computePCE2<PlusFunc>(setup.cu2, setup.cx, setup.cx);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, PlusLC2 ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->plus(setup.u, setup.cx, v);
    setup.computePCE2LC<PlusFunc>(setup.u2, setup.cx[0], v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, PlusRC2 ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->plus(setup.u, v, setup.cx);
    setup.computePCE2RC<PlusFunc>(setup.u2, v, setup.cx[0]);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, PlusResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sinh(v, setup.x);
    setup.exp->cosh(w, setup.y);
    setup.exp->plus(ru, v, w);
    setup.computePCE2<PlusFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, PlusLCResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sinh(v, setup.x);
    setup.exp->plus(ru, setup.a, v);
    setup.computePCE2LC<PlusFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, PlusRCResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sinh(v, setup.x);
    setup.exp->plus(ru, v, setup.a);
    setup.computePCE2RC<PlusFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, PlusLS ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->plus(setup.u, setup.sx, v);
    setup.computePCE2<PlusFunc>(setup.u2, setup.sx, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, PlusRS ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->plus(setup.u, v, setup.sx);
    setup.computePCE2<PlusFunc>(setup.u2, v, setup.sx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, PlusLSRC ) {
    setup.exp->plus(setup.su, setup.sx, setup.a);
    setup.computePCE2RC<PlusFunc>(setup.su2, setup.sx, setup.a);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, PlusRSLC ) {
    setup.exp->plus(setup.su, setup.a, setup.sx);
    setup.computePCE2LC<PlusFunc>(setup.su2, setup.a, setup.sx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, PlusLSRC2 ) {
    setup.exp->plus(setup.su, setup.sx, setup.cx);
    setup.computePCE2<PlusFunc>(setup.su2, setup.sx, setup.cx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, PlusRSLC2 ) {
    setup.exp->plus(setup.su, setup.cx, setup.sx);
    setup.computePCE2<PlusFunc>(setup.su2, setup.cx, setup.sx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, Minus ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->cosh(w, setup.y);
    setup.exp->minus(setup.u, v, w);
    setup.computePCE2<MinusFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, MinusLC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->minus(setup.u, setup.a, v);
    setup.computePCE2LC<MinusFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, MinusRC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->minus(setup.u, v, setup.a);
    setup.computePCE2RC<MinusFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, MinusCC ) {
    setup.exp->minus(setup.cu, setup.cx, setup.cx);
    setup.computePCE2<MinusFunc>(setup.cu2, setup.cx, setup.cx);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, MinusLC2 ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->minus(setup.u, setup.cx, v);
    setup.computePCE2LC<MinusFunc>(setup.u2, setup.cx[0], v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, MinusRC2 ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->minus(setup.u, v, setup.cx);
    setup.computePCE2RC<MinusFunc>(setup.u2, v, setup.cx[0]);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, MinusResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sinh(v, setup.x);
    setup.exp->cosh(w, setup.y);
    setup.exp->minus(ru, v, w);
    setup.computePCE2<MinusFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, MinusLCResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sinh(v, setup.x);
    setup.exp->minus(ru, setup.a, v);
    setup.computePCE2LC<MinusFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, MinusRCResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sinh(v, setup.x);
    setup.exp->minus(ru, v, setup.a);
    setup.computePCE2RC<MinusFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, MinusLS ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->minus(setup.u, setup.sx, v);
    setup.computePCE2<MinusFunc>(setup.u2, setup.sx, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, MinusRS ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->minus(setup.u, v, setup.sx);
    setup.computePCE2<MinusFunc>(setup.u2, v, setup.sx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, MinusLSRC ) {
    setup.exp->minus(setup.su, setup.sx, setup.a);
    setup.computePCE2RC<MinusFunc>(setup.su2, setup.sx, setup.a);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, MinusRSLC ) {
    setup.exp->minus(setup.su, setup.a, setup.sx);
    setup.computePCE2LC<MinusFunc>(setup.su2, setup.a, setup.sx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, MinusLSRC2 ) {
    setup.exp->minus(setup.su, setup.sx, setup.cx);
    setup.computePCE2<MinusFunc>(setup.su2, setup.sx, setup.cx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, MinusRSLC2 ) {
    setup.exp->minus(setup.su, setup.cx, setup.sx);
    setup.computePCE2<MinusFunc>(setup.su2, setup.cx, setup.sx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, Times ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->cosh(w, setup.y);
    setup.exp->times(setup.u, v, w);
    setup.computePCE2<TimesFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, TimesLC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->times(setup.u, setup.a, v);
    setup.computePCE2LC<TimesFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, TimesRC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->times(setup.u, v, setup.a);
    setup.computePCE2RC<TimesFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, TimesCC ) {
    setup.exp->times(setup.cu, setup.cx, setup.cx);
    setup.computePCE2<TimesFunc>(setup.cu2, setup.cx, setup.cx);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, TimesLC2 ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->times(setup.u, setup.cx, v);
    setup.computePCE2LC<TimesFunc>(setup.u2, setup.cx[0], v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, TimesRC2 ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->times(setup.u, v, setup.cx);
    setup.computePCE2RC<TimesFunc>(setup.u2, v, setup.cx[0]);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, TimesResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sinh(v, setup.x);
    setup.exp->cosh(w, setup.y);
    setup.exp->times(ru, v, w);
    setup.computePCE2<TimesFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, TimesLCResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sinh(v, setup.x);
    setup.exp->times(ru, setup.a, v);
    setup.computePCE2LC<TimesFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, TimesRCResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sinh(v, setup.x);
    setup.exp->times(ru, v, setup.a);
    setup.computePCE2RC<TimesFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, Divide ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->exp(w, setup.y);
    setup.exp->divide(setup.u, v, w);
    setup.computePCE2<DivideFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, DivideLC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->divide(setup.u, setup.a, v);
    setup.computePCE2LC<DivideFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, DivideRC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->divide(setup.u, v, setup.a);
    setup.computePCE2RC<DivideFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, DivideCC ) {
    setup.exp->divide(setup.cu, setup.cx, setup.cx);
    setup.computePCE2<DivideFunc>(setup.cu2, setup.cx, setup.cx);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, DivideRC2 ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->divide(setup.u, v, setup.cx);
    setup.computePCE2RC<DivideFunc>(setup.u2, v, setup.cx[0]);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
   TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, DivideResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sinh(v, setup.x);
    setup.exp->exp(w, setup.y);
    setup.exp->divide(ru, v, w);
    setup.computePCE2<DivideFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, DivideLCResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sinh(v, setup.x);
    setup.exp->divide(ru, setup.a, v);
    setup.computePCE2LC<DivideFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, DivideRCResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sinh(v, setup.x);
    setup.exp->divide(ru, v, setup.a);
    setup.computePCE2RC<DivideFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(ru, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, Pow ) {
    setup.exp->pow(setup.u, setup.x, setup.y);
    setup.computePCE2<PowFunc>(setup.u2, setup.x, setup.y);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, PowLC ) {
    setup.exp->pow(setup.u, setup.a, setup.y);
    setup.computePCE2LC<PowFunc>(setup.u2, setup.a, setup.y);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, PowRC ) {
    setup.exp->pow(setup.u, setup.x, setup.a);
    setup.computePCE2RC<PowFunc>(setup.u2, setup.x, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, PowCC ) {
    setup.exp->pow(setup.cu, setup.cx, setup.cx);
    setup.computePCE2<PowFunc>(setup.cu2, setup.cx, setup.cx);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, PowLC2 ) {
    setup.exp->pow(setup.u, setup.cx, setup.y);
    setup.computePCE2LC<PowFunc>(setup.u2, setup.cx[0], setup.y);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, PowRC2 ) {
    setup.exp->pow(setup.u, setup.x, setup.cx);
    setup.computePCE2RC<PowFunc>(setup.u2, setup.x, setup.cx[0]);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, PowResize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->pow(ru, setup.x, setup.y);
    setup.computePCE2<PowFunc>(setup.u2, setup.x, setup.y);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, PowLCResize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->pow(ru, setup.a, setup.y);
    setup.computePCE2LC<PowFunc>(setup.u2, setup.a, setup.y);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, PowRCResize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->pow(ru, setup.x, setup.a);
    setup.computePCE2RC<PowFunc>(setup.u2, setup.x, setup.a);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, PlusEqual ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->cosh(setup.u, setup.x);
    setup.computePCE2<PlusFunc>(setup.u2, setup.u, v);
    setup.exp->plusEqual(setup.u, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, PlusEqualC ) {
    setup.exp->cosh(setup.u, setup.x);
    setup.computePCE2RC<PlusFunc>(setup.u2, setup.u, setup.a);
    setup.exp->plusEqual(setup.u, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, PlusEqualC2 ) {
    setup.exp->cosh(setup.u, setup.x);
    setup.computePCE2RC<PlusFunc>(setup.u2, setup.u, setup.cx[0]);
    setup.exp->plusEqual(setup.u, setup.cx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, PlusEqualResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sinh(v, setup.x);
    setup.exp->plusEqual(ru, v);
    success = Stokhos::comparePCEs(ru, "ru", v, "v", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, PlusEqualS ) {
    setup.exp->cosh(setup.u, setup.x);
    setup.computePCE2<PlusFunc>(setup.u2, setup.u, setup.sx);
    setup.exp->plusEqual(setup.u, setup.sx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, PlusEqualSC ) {
    setup.su = setup.sx;
    setup.computePCE2RC<PlusFunc>(setup.su2, setup.su, setup.a);
    setup.exp->plusEqual(setup.su, setup.a);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, PlusEqualSC2 ) {
    setup.su = setup.sx;
    setup.computePCE2<PlusFunc>(setup.su2, setup.su, setup.cx);
    setup.exp->plusEqual(setup.su, setup.cx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
  				   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, MinusEqual ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->cosh(setup.u, setup.x);
    setup.computePCE2<MinusFunc>(setup.u2, setup.u, v);
    setup.exp->minusEqual(setup.u, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, MinusEqualC ) {
    setup.exp->cosh(setup.u, setup.x);
    setup.computePCE2RC<MinusFunc>(setup.u2, setup.u, setup.a);
    setup.exp->minusEqual(setup.u, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, MinusEqualC2 ) {
    setup.exp->cosh(setup.u, setup.x);
    setup.computePCE2RC<MinusFunc>(setup.u2, setup.u, setup.cx[0]);
    setup.exp->minusEqual(setup.u, setup.cx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, MinusEqualResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sinh(v, setup.x);
    setup.exp->minusEqual(ru, v);
    setup.exp->unaryMinus(v, v);
    success = Stokhos::comparePCEs(ru, "ru", v, "v", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, MinusEqualS ) {
    setup.exp->cosh(setup.u, setup.x);
    setup.computePCE2<MinusFunc>(setup.u2, setup.u, setup.sx);
    setup.exp->minusEqual(setup.u, setup.sx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, MinusEqualSC ) {
    setup.su = setup.sx;
    setup.computePCE2RC<MinusFunc>(setup.su2, setup.su, setup.a);
    setup.exp->minusEqual(setup.su, setup.a);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, MinusEqualSC2 ) {
    setup.su = setup.sx;
    setup.computePCE2<MinusFunc>(setup.su2, setup.su, setup.cx);
    setup.exp->minusEqual(setup.su, setup.cx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
  				   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, TimesEqual ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->cosh(setup.u, setup.x);
    setup.computePCE2<TimesFunc>(setup.u2, setup.u, v);
    setup.exp->timesEqual(setup.u, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, TimesEqualC ) {
    setup.exp->cosh(setup.u, setup.x);
    setup.computePCE2RC<TimesFunc>(setup.u2, setup.u, setup.a);
    setup.exp->timesEqual(setup.u, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, TimesEqualC2 ) {
    setup.exp->cosh(setup.u, setup.x);
    setup.computePCE2RC<TimesFunc>(setup.u2, setup.u, setup.cx[0]);
    setup.exp->timesEqual(setup.u, setup.cx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, TimesEqualResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.su = setup.sx;
    setup.computePCE2<TimesFunc>(setup.u2, setup.su, v);
    setup.exp->timesEqual(setup.su, v);
    success = Stokhos::comparePCEs(setup.su, "su", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, DivideEqual ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.exp->cosh(setup.u, setup.x);
    setup.computePCE2<DivideFunc>(setup.u2, setup.u, v);
    setup.exp->divideEqual(setup.u, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, DivideEqualC ) {
    setup.exp->cosh(setup.u, setup.x);
    setup.computePCE2RC<DivideFunc>(setup.u2, setup.u, setup.a);
    setup.exp->divideEqual(setup.u, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, DivideEqualC2 ) {
    setup.exp->cosh(setup.u, setup.x);
    setup.computePCE2RC<DivideFunc>(setup.u2, setup.u, setup.cx[0]);
    setup.exp->divideEqual(setup.u, setup.cx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ForUQTKExpansion, DivideEqualResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sinh(v, setup.x);
    setup.su = setup.sx;
    setup.computePCE2<DivideFunc>(setup.u2, setup.su, v);
    setup.exp->divideEqual(setup.su, v);
    success = Stokhos::comparePCEs(setup.su, "su", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }

  // Not testing atan2(), max(), min(), abs(), fabs() since these are
  // not smooth functions

}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
