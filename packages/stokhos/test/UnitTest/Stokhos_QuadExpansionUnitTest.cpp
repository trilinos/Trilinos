// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2008) Sandia Corporation
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

namespace QuadExpansionUnitTest {

  // Common setup for unit tests
  template <typename OrdinalType, typename ValueType>
  struct UnitTestSetup {
    ValueType rtol, atol;
    OrdinalType sz;
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<OrdinalType,ValueType> > basis;
    Teuchos::RCP<const Stokhos::Quadrature<OrdinalType,ValueType> > quad;
    Teuchos::RCP< Stokhos::QuadOrthogPolyExpansion<OrdinalType,ValueType> > exp;
    Stokhos::OrthogPolyApprox<OrdinalType,ValueType> x, y, u, u2;
    ValueType a;
    
    UnitTestSetup() {
      rtol = 1e-4;
      atol = 1e-5;
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
      
    // Quadrature expansion
      exp = 
	Teuchos::rcp(new Stokhos::QuadOrthogPolyExpansion<OrdinalType,ValueType>(basis, quad));
      
      // Create approximation
      sz = basis->size();
      x.reset(basis);
      y.reset(basis);
      u.reset(basis); 
      u2.reset(basis);
      x.term(0, 0) = 1.0;
      for (OrdinalType i=0; i<d; i++)
	x.term(i, 1) = 0.1;
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
      for (OrdinalType i=0; i<sz; i++)
	c[i] = 0.0;
      
      // Compute PCE via quadrature
      Func func;
      for (OrdinalType k=0; k<nqp; k++) {
	ValueType val = a.evaluate(points[k], values[k]);
	val = func(val);
	for (int i=0; i<sz; i++)
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
      for (OrdinalType i=0; i<sz; i++)
	c[i] = 0.0;
      
      // Compute PCE via quadrature
      Func func;
      for (OrdinalType k=0; k<nqp; k++) {
	ValueType val1 = a.evaluate(points[k], values[k]);
	ValueType val2 = b.evaluate(points[k], values[k]);
	ValueType val = func(val1, val2);
	for (int i=0; i<sz; i++)
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
      for (OrdinalType i=0; i<sz; i++)
	c[i] = 0.0;
      
      // Compute PCE via quadrature
      Func func;
      for (OrdinalType k=0; k<nqp; k++) {
	ValueType val2 = b.evaluate(points[k], values[k]);
	ValueType val = func(a, val2);
	for (int i=0; i<sz; i++)
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
      for (OrdinalType i=0; i<sz; i++)
	c[i] = 0.0;
      
      // Compute PCE via quadrature
      Func func;
      for (OrdinalType k=0; k<nqp; k++) {
	ValueType val1 = a.evaluate(points[k], values[k]);
	ValueType val = func(val1, b);
	for (int i=0; i<sz; i++)
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

  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, UMinus ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
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
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, Sin ) {
    setup.exp->sin(setup.u, setup.x);
    setup.computePCE1<SinFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, Cos ) {
    setup.exp->cos(setup.u, setup.x);
    setup.computePCE1<CosFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, Tan ) {
    setup.exp->tan(setup.u, setup.x);
    setup.computePCE1<TanFunc>(setup.u2, setup.x);
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
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, ASin ) {
    setup.exp->asin(setup.u, setup.x);
    setup.computePCE1<ASinFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, ACos ) {
    setup.exp->acos(setup.u, setup.x);
    setup.computePCE1<ACosFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, ATan ) {
    setup.exp->atan(setup.u, setup.x);
    setup.computePCE1<ATanFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, ASinh ) {
    setup.exp->asinh(setup.u, setup.x);
    setup.computePCE1<ASinhFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, ACosh ) {
    setup.exp->acosh(setup.u, setup.x);
    setup.computePCE1<ACoshFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, ATanh ) {
    setup.exp->atanh(setup.u, setup.x);
    setup.computePCE1<ATanhFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, Plus ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->cos(w, setup.y);
    setup.exp->plus(setup.u, v, w);
    setup.computePCE2<PlusFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, PlusLC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->plus(setup.u, setup.a, v);
    setup.computePCE2LC<PlusFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, PlusRC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->plus(setup.u, v, setup.a);
    setup.computePCE2RC<PlusFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, Minus ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->cos(w, setup.y);
    setup.exp->minus(setup.u, v, w);
    setup.computePCE2<MinusFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, MinusLC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->minus(setup.u, setup.a, v);
    setup.computePCE2LC<MinusFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, MinusRC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->minus(setup.u, v, setup.a);
    setup.computePCE2RC<MinusFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, Times ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->cos(w, setup.y);
    setup.exp->times(setup.u, v, w);
    setup.computePCE2<TimesFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, TimesLC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->times(setup.u, setup.a, v);
    setup.computePCE2LC<TimesFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, TimesRC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->times(setup.u, v, setup.a);
    setup.computePCE2RC<TimesFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, Divide ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->exp(w, setup.y);
    setup.exp->divide(setup.u, v, w);
    setup.computePCE2<DivideFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, DivideLC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->divide(setup.u, setup.a, v);
    setup.computePCE2LC<DivideFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, DivideRC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
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
    setup.exp->sin(v, setup.x);
    setup.exp->cos(setup.u, setup.x);
    setup.computePCE2<PlusFunc>(setup.u2, setup.u, v);
    setup.exp->plusEqual(setup.u, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, PlusEqualC ) {
    setup.exp->cos(setup.u, setup.x);
    setup.computePCE2RC<PlusFunc>(setup.u2, setup.u, setup.a);
    setup.exp->plusEqual(setup.u, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, MinusEqual ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->cos(setup.u, setup.x);
    setup.computePCE2<MinusFunc>(setup.u2, setup.u, v);
    setup.exp->minusEqual(setup.u, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, MinusEqualC ) {
    setup.exp->cos(setup.u, setup.x);
    setup.computePCE2RC<MinusFunc>(setup.u2, setup.u, setup.a);
    setup.exp->minusEqual(setup.u, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, TimesEqual ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->cos(setup.u, setup.x);
    setup.computePCE2<TimesFunc>(setup.u2, setup.u, v);
    setup.exp->timesEqual(setup.u, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, TimesEqualC ) {
    setup.exp->cos(setup.u, setup.x);
    setup.computePCE2RC<TimesFunc>(setup.u2, setup.u, setup.a);
    setup.exp->timesEqual(setup.u, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, DivideEqual ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->cos(setup.u, setup.x);
    setup.computePCE2<DivideFunc>(setup.u2, setup.u, v);
    setup.exp->divideEqual(setup.u, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_QuadExpansion, DivideEqualC ) {
    setup.exp->cos(setup.u, setup.x);
    setup.computePCE2RC<DivideFunc>(setup.u2, setup.u, setup.a);
    setup.exp->divideEqual(setup.u, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }

  // Not testing atan2(), max(), min(), abs(), fabs() since these are
  // not smooth functions

}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
