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
#include "Stokhos_PseudoSpectralOrthogPolyExpansion.hpp"

namespace PseudoSpectralExpansionUnitTest {

  // Common setup for unit tests
  template <typename OrdinalType, typename ValueType>
  struct UnitTestSetup {
    //typedef Stokhos::TensorProductBasis<OrdinalType,ValueType> product_basis_type;
    typedef Stokhos::SmolyakBasis<OrdinalType,ValueType> product_basis_type;
    ValueType rtol, atol;
    ValueType crtol, catol;
    OrdinalType sz;
    Teuchos::RCP<const product_basis_type> basis;
    Teuchos::RCP<const Stokhos::Quadrature<OrdinalType,ValueType> > quad;
    Teuchos::RCP<const Stokhos::PseudoSpectralOperator<OrdinalType,ValueType> > ps_op;
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk, Cijk_linear;
    Teuchos::RCP< Stokhos::PseudoSpectralOrthogPolyExpansion<OrdinalType,ValueType> > exp, exp_linear;
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

      Stokhos::TotalOrderIndexSet<OrdinalType> coeff_index_set(d, p);
      Teuchos::Array< Stokhos::IdentityGrowthRule<OrdinalType> > coeff_growth(d);
      basis =
	Teuchos::rcp(new product_basis_type(bases, coeff_index_set, coeff_growth));
      
      // Tensor product quadrature
      quad = 
	Teuchos::rcp(new Stokhos::TensorProductQuadrature<OrdinalType,ValueType>(basis));

      // Tensor product pseudospectral operator
      Teuchos::Array< Stokhos::EvenGrowthRule<OrdinalType> > point_growth(d);
      ps_op = 
	Teuchos::rcp(new Stokhos::SmolyakPseudoSpectralOperator<OrdinalType,ValueType>(*basis, point_growth, true, true));

      // Triple product tensor
      Cijk = basis->computeTripleProductTensor(basis->size());
      Cijk_linear = basis->computeTripleProductTensor(basis->dimension()+1);
      
      // Quadrature expansion
      exp = 
	Teuchos::rcp(new Stokhos::PseudoSpectralOrthogPolyExpansion<OrdinalType,ValueType>(basis, Cijk, ps_op));
      exp_linear = 
	Teuchos::rcp(new Stokhos::PseudoSpectralOrthogPolyExpansion<OrdinalType,ValueType>(basis, Cijk_linear, ps_op));
      
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

  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, UMinus ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->unaryMinus(setup.u, v);
    setup.computePCE1<UMinusFunc>(setup.u2, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, Exp ) {
    setup.exp->exp(setup.u, setup.x);
    setup.computePCE1<ExpFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, ExpConst ) {
    setup.exp->exp(setup.cu, setup.cx);
    setup.cu2[0] = std::exp(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, ExpResize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->exp(ru, setup.x);
    setup.computePCE1<ExpFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, Log ) {
    setup.exp->log(setup.u, setup.x);
    setup.computePCE1<LogFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, LogConst ) {
    setup.exp->log(setup.cu, setup.cx);
    setup.cu2[0] = std::log(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, LogResize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->log(ru, setup.x);
    setup.computePCE1<LogFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, Log10 ) {
    setup.exp->log10(setup.u, setup.x);
    setup.computePCE1<Log10Func>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, Log10Const ) {
    setup.exp->log10(setup.cu, setup.cx);
    setup.cu2[0] = std::log10(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, Log10Resize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->log10(ru, setup.x);
    setup.computePCE1<Log10Func>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, Sqrt ) {
    setup.exp->sqrt(setup.u, setup.x);
    setup.computePCE1<SqrtFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, SqrtConst ) {
    setup.exp->sqrt(setup.cu, setup.cx);
    setup.cu2[0] = std::sqrt(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, SqrtResize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sqrt(ru, setup.x);
    setup.computePCE1<SqrtFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, Sin ) {
    setup.exp->sin(setup.u, setup.x);
    setup.computePCE1<SinFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, SinConst ) {
    setup.exp->sin(setup.cu, setup.cx);
    setup.cu2[0] = std::sin(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, SinResize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sin(ru, setup.x);
    setup.computePCE1<SinFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, Cos ) {
    setup.exp->cos(setup.u, setup.x);
    setup.computePCE1<CosFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, CosConst ) {
    setup.exp->cos(setup.cu, setup.cx);
    setup.cu2[0] = std::cos(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, CosResize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->cos(ru, setup.x);
    setup.computePCE1<CosFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, Tan ) {
    setup.exp->tan(setup.u, setup.x);
    setup.computePCE1<TanFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, TanConst ) {
    setup.exp->tan(setup.cu, setup.cx);
    setup.cu2[0] = std::tan(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, TanResize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->tan(ru, setup.x);
    setup.computePCE1<TanFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, Sinh ) {
    setup.exp->sinh(setup.u, setup.x);
    setup.computePCE1<SinhFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, SinhConst ) {
    setup.exp->sinh(setup.cu, setup.cx);
    setup.cu2[0] = std::sinh(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, SinhResize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sinh(ru, setup.x);
    setup.computePCE1<SinhFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, Cosh ) {
    setup.exp->cosh(setup.u, setup.x);
    setup.computePCE1<CoshFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, CoshConst ) {
    setup.exp->cosh(setup.cu, setup.cx);
    setup.cu2[0] = std::cosh(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, CoshResize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->cosh(ru, setup.x);
    setup.computePCE1<CoshFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, Tanh ) {
    setup.exp->tanh(setup.u, setup.x);
    setup.computePCE1<TanhFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, TanhConst ) {
    setup.exp->tanh(setup.cu, setup.cx);
    setup.cu2[0] = std::tanh(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, TanhResize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->tanh(ru, setup.x);
    setup.computePCE1<TanhFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, ASin ) {
    setup.exp->asin(setup.u, setup.x);
    setup.computePCE1<ASinFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, ASinConst ) {
    setup.exp->asin(setup.cu, setup.cx);
    setup.cu2[0] = std::asin(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, ASinResize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->asin(ru, setup.x);
    setup.computePCE1<ASinFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, ACos ) {
    setup.exp->acos(setup.u, setup.x);
    setup.computePCE1<ACosFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, ACosConst ) {
    setup.exp->acos(setup.cu, setup.cx);
    setup.cu2[0] = std::acos(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, ACosResize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->acos(ru, setup.x);
    setup.computePCE1<ACosFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, ATan ) {
    setup.exp->atan(setup.u, setup.x);
    setup.computePCE1<ATanFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, ATanConst ) {
    setup.exp->atan(setup.cu, setup.cx);
    setup.cu2[0] = std::atan(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, ATanResize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->atan(ru, setup.x);
    setup.computePCE1<ATanFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, ASinh ) {
    setup.exp->asinh(setup.u, setup.x);
    setup.computePCE1<ASinhFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, ASinhConst ) {
    ASinhFunc f;
    setup.exp->asinh(setup.cu, setup.cx);
    setup.cu2[0] = f(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, ASinhResize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->asinh(ru, setup.x);
    setup.computePCE1<ASinhFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, ACosh ) {
    setup.exp->acosh(setup.u, setup.x);
    setup.computePCE1<ACoshFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, ACoshConst ) {
    ACoshFunc f;
    setup.exp->acosh(setup.cu, setup.cx);
    setup.cu2[0] = f(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, ACoshResize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->acosh(ru, setup.x);
    setup.computePCE1<ACoshFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, ATanh ) {
    setup.exp->atanh(setup.u, setup.x);
    setup.computePCE1<ATanhFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, ATanhConst ) {
    ATanhFunc f;
    setup.exp->atanh(setup.cu, setup.cx);
    setup.cu2[0] = f(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, ATanhResize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->atanh(ru, setup.x);
    setup.computePCE1<ATanhFunc>(setup.u2, setup.x);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, Plus ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->cos(w, setup.y);
    setup.exp->plus(setup.u, v, w);
    setup.computePCE2<PlusFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PlusLC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->plus(setup.u, setup.a, v);
    setup.computePCE2LC<PlusFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PlusRC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->plus(setup.u, v, setup.a);
    setup.computePCE2RC<PlusFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PlusCC ) {
    setup.exp->plus(setup.cu, setup.cx, setup.cx);
    setup.computePCE2<PlusFunc>(setup.cu2, setup.cx, setup.cx);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PlusLC2 ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->plus(setup.u, setup.cx, v);
    setup.computePCE2LC<PlusFunc>(setup.u2, setup.cx[0], v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PlusRC2 ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->plus(setup.u, v, setup.cx);
    setup.computePCE2RC<PlusFunc>(setup.u2, v, setup.cx[0]);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PlusResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sin(v, setup.x);
    setup.exp->cos(w, setup.y);
    setup.exp->plus(ru, v, w);
    setup.computePCE2<PlusFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PlusLCResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sin(v, setup.x);
    setup.exp->plus(ru, setup.a, v);
    setup.computePCE2LC<PlusFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PlusRCResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sin(v, setup.x);
    setup.exp->plus(ru, v, setup.a);
    setup.computePCE2RC<PlusFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PlusLS ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->plus(setup.u, setup.sx, v);
    setup.computePCE2<PlusFunc>(setup.u2, setup.sx, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PlusRS ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->plus(setup.u, v, setup.sx);
    setup.computePCE2<PlusFunc>(setup.u2, v, setup.sx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PlusLSRC ) {
    setup.exp->plus(setup.su, setup.sx, setup.a);
    setup.computePCE2RC<PlusFunc>(setup.su2, setup.sx, setup.a);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PlusRSLC ) {
    setup.exp->plus(setup.su, setup.a, setup.sx);
    setup.computePCE2LC<PlusFunc>(setup.su2, setup.a, setup.sx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PlusLSRC2 ) {
    setup.exp->plus(setup.su, setup.sx, setup.cx);
    setup.computePCE2<PlusFunc>(setup.su2, setup.sx, setup.cx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PlusRSLC2 ) {
    setup.exp->plus(setup.su, setup.cx, setup.sx);
    setup.computePCE2<PlusFunc>(setup.su2, setup.cx, setup.sx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, Minus ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->cos(w, setup.y);
    setup.exp->minus(setup.u, v, w);
    setup.computePCE2<MinusFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, MinusLC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->minus(setup.u, setup.a, v);
    setup.computePCE2LC<MinusFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, MinusRC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->minus(setup.u, v, setup.a);
    setup.computePCE2RC<MinusFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, MinusCC ) {
    setup.exp->minus(setup.cu, setup.cx, setup.cx);
    setup.computePCE2<MinusFunc>(setup.cu2, setup.cx, setup.cx);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, MinusLC2 ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->minus(setup.u, setup.cx, v);
    setup.computePCE2LC<MinusFunc>(setup.u2, setup.cx[0], v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, MinusRC2 ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->minus(setup.u, v, setup.cx);
    setup.computePCE2RC<MinusFunc>(setup.u2, v, setup.cx[0]);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, MinusResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sin(v, setup.x);
    setup.exp->cos(w, setup.y);
    setup.exp->minus(ru, v, w);
    setup.computePCE2<MinusFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, MinusLCResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sin(v, setup.x);
    setup.exp->minus(ru, setup.a, v);
    setup.computePCE2LC<MinusFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, MinusRCResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sin(v, setup.x);
    setup.exp->minus(ru, v, setup.a);
    setup.computePCE2RC<MinusFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, MinusLS ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->minus(setup.u, setup.sx, v);
    setup.computePCE2<MinusFunc>(setup.u2, setup.sx, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, MinusRS ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->minus(setup.u, v, setup.sx);
    setup.computePCE2<MinusFunc>(setup.u2, v, setup.sx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, MinusLSRC ) {
    setup.exp->minus(setup.su, setup.sx, setup.a);
    setup.computePCE2RC<MinusFunc>(setup.su2, setup.sx, setup.a);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, MinusRSLC ) {
    setup.exp->minus(setup.su, setup.a, setup.sx);
    setup.computePCE2LC<MinusFunc>(setup.su2, setup.a, setup.sx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, MinusLSRC2 ) {
    setup.exp->minus(setup.su, setup.sx, setup.cx);
    setup.computePCE2<MinusFunc>(setup.su2, setup.sx, setup.cx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, MinusRSLC2 ) {
    setup.exp->minus(setup.su, setup.cx, setup.sx);
    setup.computePCE2<MinusFunc>(setup.su2, setup.cx, setup.sx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, Times ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->cos(w, setup.y);
    setup.exp->times(setup.u, v, w);
    setup.computePCE2<TimesFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, TimesLC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->times(setup.u, setup.a, v);
    setup.computePCE2LC<TimesFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, TimesRC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->times(setup.u, v, setup.a);
    setup.computePCE2RC<TimesFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, TimesCC ) {
    setup.exp->times(setup.cu, setup.cx, setup.cx);
    setup.computePCE2<TimesFunc>(setup.cu2, setup.cx, setup.cx);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, TimesLC2 ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->times(setup.u, setup.cx, v);
    setup.computePCE2LC<TimesFunc>(setup.u2, setup.cx[0], v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, TimesRC2 ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->times(setup.u, v, setup.cx);
    setup.computePCE2RC<TimesFunc>(setup.u2, v, setup.cx[0]);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, TimesResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sin(v, setup.x);
    setup.exp->cos(w, setup.y);
    setup.exp->times(ru, v, w);
    setup.computePCE2<TimesFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, TimesLCResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sin(v, setup.x);
    setup.exp->times(ru, setup.a, v);
    setup.computePCE2LC<TimesFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, TimesRCResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sin(v, setup.x);
    setup.exp->times(ru, v, setup.a);
    setup.computePCE2RC<TimesFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, TimesLS ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->times(setup.u, setup.sx, v);
    setup.computePCE2<TimesFunc>(setup.u2, setup.sx, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, TimesRS ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->times(setup.u, v, setup.sx);
    setup.computePCE2<TimesFunc>(setup.u2, v, setup.sx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, TimesLSLinear ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp_linear->times(setup.u, setup.sx, v);
    setup.computePCE2<TimesFunc>(setup.u2, setup.sx, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, TimesRSLinear ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp_linear->times(setup.u, v, setup.sx);
    setup.computePCE2<TimesFunc>(setup.u2, v, setup.sx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, TimesLSRC ) {
    setup.exp->times(setup.su, setup.sx, setup.a);
    setup.computePCE2RC<TimesFunc>(setup.su2, setup.sx, setup.a);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, TimesRSLC ) {
    setup.exp->times(setup.su, setup.a, setup.sx);
    setup.computePCE2LC<TimesFunc>(setup.su2, setup.a, setup.sx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, TimesLSRC2 ) {
    setup.exp->times(setup.su, setup.sx, setup.cx);
    setup.computePCE2<TimesFunc>(setup.su2, setup.sx, setup.cx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, TimesRSLC2 ) {
    setup.exp->times(setup.su, setup.cx, setup.sx);
    setup.computePCE2<TimesFunc>(setup.su2, setup.cx, setup.sx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, Divide ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->exp(w, setup.y);
    setup.exp->divide(setup.u, v, w);
    setup.computePCE2<DivideFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, DivideLC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->divide(setup.u, setup.a, v);
    setup.computePCE2LC<DivideFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, DivideRC ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->divide(setup.u, v, setup.a);
    setup.computePCE2RC<DivideFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, DivideCC ) {
    setup.exp->divide(setup.cu, setup.cx, setup.cx);
    setup.computePCE2<DivideFunc>(setup.cu2, setup.cx, setup.cx);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, DivideLC2 ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->divide(setup.u, setup.cx, v);
    setup.computePCE2LC<DivideFunc>(setup.u2, setup.cx[0], v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, DivideRC2 ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->divide(setup.u, v, setup.cx);
    setup.computePCE2RC<DivideFunc>(setup.u2, v, setup.cx[0]);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
   TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, DivideResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis), w(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sin(v, setup.x);
    setup.exp->exp(w, setup.y);
    setup.exp->divide(ru, v, w);
    setup.computePCE2<DivideFunc>(setup.u2, v, w);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, DivideLCResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sin(v, setup.x);
    setup.exp->divide(ru, setup.a, v);
    setup.computePCE2LC<DivideFunc>(setup.u2, setup.a, v);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, DivideRCResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sin(v, setup.x);
    setup.exp->divide(ru, v, setup.a);
    setup.computePCE2RC<DivideFunc>(setup.u2, v, setup.a);
    success = Stokhos::comparePCEs(ru, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, DivideLS ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->divide(setup.u, setup.sx, v);
    setup.computePCE2<DivideFunc>(setup.u2, setup.sx, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, DivideRS ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->divide(setup.u, v, setup.sx);
    setup.computePCE2<DivideFunc>(setup.u2, v, setup.sx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, DivideLSRC ) {
    setup.exp->divide(setup.su, setup.sx, setup.a);
    setup.computePCE2RC<DivideFunc>(setup.su2, setup.sx, setup.a);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, DivideRSLC ) {
    setup.exp->divide(setup.u, setup.a, setup.sx);
    setup.computePCE2LC<DivideFunc>(setup.u2, setup.a, setup.sx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, DivideLSRC2 ) {
    setup.exp->divide(setup.su, setup.sx, setup.cx);
    setup.computePCE2<DivideFunc>(setup.su2, setup.sx, setup.cx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, DivideRSLC2 ) {
    setup.exp->divide(setup.u, setup.cx, setup.sx);
    setup.computePCE2<DivideFunc>(setup.u2, setup.cx, setup.sx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, Pow ) {
    setup.exp->pow(setup.u, setup.x, setup.y);
    setup.computePCE2<PowFunc>(setup.u2, setup.x, setup.y);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PowLC ) {
    setup.exp->pow(setup.u, setup.a, setup.y);
    setup.computePCE2LC<PowFunc>(setup.u2, setup.a, setup.y);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PowRC ) {
    setup.exp->pow(setup.u, setup.x, setup.a);
    setup.computePCE2RC<PowFunc>(setup.u2, setup.x, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PowCC ) {
    setup.exp->pow(setup.cu, setup.cx, setup.cx);
    setup.computePCE2<PowFunc>(setup.cu2, setup.cx, setup.cx);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PowLC2 ) {
    setup.exp->pow(setup.u, setup.cx, setup.y);
    setup.computePCE2LC<PowFunc>(setup.u2, setup.cx[0], setup.y);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PowRC2 ) {
    setup.exp->pow(setup.u, setup.x, setup.cx);
    setup.computePCE2RC<PowFunc>(setup.u2, setup.x, setup.cx[0]);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PowResize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->pow(ru, setup.x, setup.y);
    setup.computePCE2<PowFunc>(setup.u2, setup.x, setup.y);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PowLCResize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->pow(ru, setup.a, setup.y);
    setup.computePCE2LC<PowFunc>(setup.u2, setup.a, setup.y);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PowRCResize ) {
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->pow(ru, setup.x, setup.a);
    setup.computePCE2RC<PowFunc>(setup.u2, setup.x, setup.a);
    success = Stokhos::comparePCEs(ru, "ru", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PowLS ) {
    setup.exp->pow(setup.u, setup.sx, setup.y);
    setup.computePCE2<PowFunc>(setup.u2, setup.sx, setup.y);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PowRS ) {
    setup.exp->pow(setup.u, setup.x, setup.sx);
    setup.computePCE2<PowFunc>(setup.u2, setup.x, setup.sx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PowLSRC ) {
    setup.exp->pow(setup.u, setup.sx, setup.a);
    setup.computePCE2RC<PowFunc>(setup.u2, setup.sx, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PowRSLC ) {
    setup.exp->pow(setup.u, setup.a, setup.sx);
    setup.computePCE2LC<PowFunc>(setup.u2, setup.a, setup.sx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PowLSRC2 ) {
    setup.exp->pow(setup.u, setup.sx, setup.cx);
    setup.computePCE2<PowFunc>(setup.u2, setup.sx, setup.cx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PowRSLC2 ) {
    setup.exp->pow(setup.u, setup.cx, setup.sx);
    setup.computePCE2<PowFunc>(setup.u2, setup.cx, setup.sx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PlusEqual ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->cos(setup.u, setup.x);
    setup.computePCE2<PlusFunc>(setup.u2, setup.u, v);
    setup.exp->plusEqual(setup.u, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PlusEqualC ) {
    setup.exp->cos(setup.u, setup.x);
    setup.computePCE2RC<PlusFunc>(setup.u2, setup.u, setup.a);
    setup.exp->plusEqual(setup.u, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PlusEqualC2 ) {
    setup.exp->cos(setup.u, setup.x);
    setup.computePCE2RC<PlusFunc>(setup.u2, setup.u, setup.cx[0]);
    setup.exp->plusEqual(setup.u, setup.cx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PlusEqualResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sin(v, setup.x);
    setup.exp->plusEqual(ru, v);
    success = Stokhos::comparePCEs(ru, "ru", v, "v", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PlusEqualS ) {
    setup.exp->cos(setup.u, setup.x);
    setup.computePCE2<PlusFunc>(setup.u2, setup.u, setup.sx);
    setup.exp->plusEqual(setup.u, setup.sx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PlusEqualSC ) {
    setup.su = setup.sx;
    setup.computePCE2RC<PlusFunc>(setup.su2, setup.su, setup.a);
    setup.exp->plusEqual(setup.su, setup.a);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, PlusEqualSC2 ) {
    setup.su = setup.sx;
    setup.computePCE2<PlusFunc>(setup.su2, setup.su, setup.cx);
    setup.exp->plusEqual(setup.su, setup.cx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
  				   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, MinusEqual ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->cos(setup.u, setup.x);
    setup.computePCE2<MinusFunc>(setup.u2, setup.u, v);
    setup.exp->minusEqual(setup.u, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, MinusEqualC ) {
    setup.exp->cos(setup.u, setup.x);
    setup.computePCE2RC<MinusFunc>(setup.u2, setup.u, setup.a);
    setup.exp->minusEqual(setup.u, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, MinusEqualC2 ) {
    setup.exp->cos(setup.u, setup.x);
    setup.computePCE2RC<MinusFunc>(setup.u2, setup.u, setup.cx[0]);
    setup.exp->minusEqual(setup.u, setup.cx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, MinusEqualResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    Stokhos::OrthogPolyApprox<int, double> ru(setup.basis, 0);
    setup.exp->sin(v, setup.x);
    setup.exp->minusEqual(ru, v);
    setup.exp->unaryMinus(v, v);
    success = Stokhos::comparePCEs(ru, "ru", v, "v", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, MinusEqualS ) {
    setup.exp->cos(setup.u, setup.x);
    setup.computePCE2<MinusFunc>(setup.u2, setup.u, setup.sx);
    setup.exp->minusEqual(setup.u, setup.sx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, MinusEqualSC ) {
    setup.su = setup.sx;
    setup.computePCE2RC<MinusFunc>(setup.su2, setup.su, setup.a);
    setup.exp->minusEqual(setup.su, setup.a);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, MinusEqualSC2 ) {
    setup.su = setup.sx;
    setup.computePCE2<MinusFunc>(setup.su2, setup.su, setup.cx);
    setup.exp->minusEqual(setup.su, setup.cx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
  				   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, TimesEqual ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->cos(setup.u, setup.x);
    setup.computePCE2<TimesFunc>(setup.u2, setup.u, v);
    setup.exp->timesEqual(setup.u, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, TimesEqualC ) {
    setup.exp->cos(setup.u, setup.x);
    setup.computePCE2RC<TimesFunc>(setup.u2, setup.u, setup.a);
    setup.exp->timesEqual(setup.u, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, TimesEqualC2 ) {
    setup.exp->cos(setup.u, setup.x);
    setup.computePCE2RC<TimesFunc>(setup.u2, setup.u, setup.cx[0]);
    setup.exp->timesEqual(setup.u, setup.cx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, TimesEqualResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.su = setup.sx;
    setup.computePCE2<TimesFunc>(setup.u2, setup.su, v);
    setup.exp->timesEqual(setup.su, v);
    success = Stokhos::comparePCEs(setup.su, "su", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, TimesEqualS ) {
    setup.exp->cos(setup.u, setup.x);
    setup.computePCE2<TimesFunc>(setup.u2, setup.u, setup.sx);
    setup.exp->timesEqual(setup.u, setup.sx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, TimesEqualSLinear ) {
    setup.exp->cos(setup.u, setup.x);
    setup.computePCE2<TimesFunc>(setup.u2, setup.u, setup.sx);
    setup.exp_linear->timesEqual(setup.u, setup.sx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, TimesEqualSC ) {
    setup.su = setup.sx;
    setup.computePCE2RC<TimesFunc>(setup.su2, setup.su, setup.a);
    setup.exp->timesEqual(setup.su, setup.a);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, TimesEqualSC2 ) {
    setup.su = setup.sx;
    setup.computePCE2<TimesFunc>(setup.su2, setup.su, setup.cx);
    setup.exp->timesEqual(setup.su, setup.cx);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
  				   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, DivideEqual ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.exp->cos(setup.u, setup.x);
    setup.computePCE2<DivideFunc>(setup.u2, setup.u, v);
    setup.exp->divideEqual(setup.u, v);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, DivideEqualC ) {
    setup.exp->cos(setup.u, setup.x);
    setup.computePCE2RC<DivideFunc>(setup.u2, setup.u, setup.a);
    setup.exp->divideEqual(setup.u, setup.a);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, DivideEqualC2 ) {
    setup.exp->cos(setup.u, setup.x);
    setup.computePCE2RC<DivideFunc>(setup.u2, setup.u, setup.cx[0]);
    setup.exp->divideEqual(setup.u, setup.cx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, DivideEqualResize ) {
    Stokhos::OrthogPolyApprox<int,double> v(setup.basis);
    setup.exp->sin(v, setup.x);
    setup.su = setup.sx;
    setup.computePCE2<DivideFunc>(setup.u2, setup.su, v);
    setup.exp->divideEqual(setup.su, v);
    success = Stokhos::comparePCEs(setup.su, "su", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, DivideEqualS ) {
    setup.exp->cos(setup.u, setup.x);
    setup.computePCE2<DivideFunc>(setup.u2, setup.u, setup.sx);
    setup.exp->divideEqual(setup.u, setup.sx);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, DivideEqualSC ) {
    setup.su = setup.sx;
    setup.computePCE2RC<DivideFunc>(setup.su2, setup.su, setup.a);
    setup.exp->divideEqual(setup.su, setup.a);
    success = Stokhos::comparePCEs(setup.su, "su", setup.su2, "su2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PSExpansion, DivideEqualSC2 ) {
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
