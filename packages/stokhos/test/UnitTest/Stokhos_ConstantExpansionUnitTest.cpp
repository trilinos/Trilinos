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

namespace ConstantExpansionUnitTest {

  // Common setup for unit tests
  template <typename OrdinalType, typename ValueType>
  struct UnitTestSetup {
    ValueType rtol, atol;
    ValueType crtol, catol;
    OrdinalType sz;
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<OrdinalType,ValueType> > basis;
    Teuchos::RCP< Stokhos::ConstantOrthogPolyExpansion<OrdinalType,ValueType> > exp;
    Stokhos::OrthogPolyApprox<OrdinalType,ValueType> cx, cu, cu2;
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
      
      // Constant expansion
      exp = 
	Teuchos::rcp(new Stokhos::ConstantOrthogPolyExpansion<OrdinalType,ValueType>());
      
      // Create approximation
      cx.reset(basis, 1);
      cx.term(0, 0) = a;
      cu.reset(basis, 1);
      cu2.reset(basis, 1);
    }
    
  };

  UnitTestSetup<int,double> setup;
  
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

  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, UMinus ) {
    setup.exp->unaryMinus(setup.cu, setup.cx);
    setup.cu2[0] = -setup.cx[0];
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, ExpConst ) {
    setup.exp->exp(setup.cu, setup.cx);
    setup.cu2[0] = std::exp(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, LogConst ) {
    setup.exp->log(setup.cu, setup.cx);
    setup.cu2[0] = std::log(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, Log10Const ) {
    setup.exp->log10(setup.cu, setup.cx);
    setup.cu2[0] = std::log10(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, SqrtConst ) {
    setup.exp->sqrt(setup.cu, setup.cx);
    setup.cu2[0] = std::sqrt(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, CbrtConst ) {
    setup.exp->cbrt(setup.cu, setup.cx);
    setup.cu2[0] = std::cbrt(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, SinConst ) {
    setup.exp->sin(setup.cu, setup.cx);
    setup.cu2[0] = std::sin(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, CosConst ) {
    setup.exp->cos(setup.cu, setup.cx);
    setup.cu2[0] = std::cos(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, TanConst ) {
    setup.exp->tan(setup.cu, setup.cx);
    setup.cu2[0] = std::tan(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, SinhConst ) {
    setup.exp->sinh(setup.cu, setup.cx);
    setup.cu2[0] = std::sinh(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, CoshConst ) {
    setup.exp->cosh(setup.cu, setup.cx);
    setup.cu2[0] = std::cosh(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, TanhConst ) {
    setup.exp->tanh(setup.cu, setup.cx);
    setup.cu2[0] = std::tanh(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, ASinConst ) {
    setup.exp->asin(setup.cu, setup.cx);
    setup.cu2[0] = std::asin(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, ACosConst ) {
    setup.exp->acos(setup.cu, setup.cx);
    setup.cu2[0] = std::acos(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, ATanConst ) {
    setup.exp->atan(setup.cu, setup.cx);
    setup.cu2[0] = std::atan(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, ASinhConst ) {
    ASinhFunc f;
    setup.exp->asinh(setup.cu, setup.cx);
    setup.cu2[0] = f(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, ACoshConst ) {
    ACoshFunc f;
    setup.exp->acosh(setup.cu, setup.cx);
    setup.cu2[0] = f(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, ATanhConst ) {
    ATanhFunc f;
    setup.exp->atanh(setup.cu, setup.cx);
    setup.cu2[0] = f(setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, Plus ) {
    setup.exp->plus(setup.cu, setup.cx, setup.cx);
    setup.cu2[0] = setup.cx[0] + setup.cx[0];
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, PlusLC ) {
    setup.exp->plus(setup.cu, setup.a, setup.cx);
    setup.cu2[0] = setup.a + setup.cx[0];
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, PlusRC ) {
    setup.exp->plus(setup.cu, setup.cx, setup.a);
    setup.cu2[0] = setup.cx[0] + setup.a;
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, Minus ) {
    setup.exp->minus(setup.cu, setup.cx, setup.cx);
    setup.cu2[0] = setup.cx[0] - setup.cx[0];
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, MinusLC ) {
    setup.exp->minus(setup.cu, setup.a, setup.cx);
    setup.cu2[0] = setup.a - setup.cx[0];
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, MinusRC ) {
    setup.exp->minus(setup.cu, setup.cx, setup.a);
    setup.cu2[0] = setup.cx[0] - setup.a;
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.rtol, setup.atol, out);
  }
  
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, Times ) {
    setup.exp->times(setup.cu, setup.cx, setup.cx);
    setup.cu2[0] = setup.cx[0] * setup.cx[0];
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, TimesLC ) {
    setup.exp->times(setup.cu, setup.a, setup.cx);
    setup.cu2[0] = setup.a * setup.cx[0];
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, TimesRC ) {
    setup.exp->times(setup.cu, setup.cx, setup.a);
    setup.cu2[0] = setup.cx[0] * setup.a;
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, Divide ) {
    setup.exp->divide(setup.cu, setup.cx, setup.cx);
    setup.cu2[0] = setup.cx[0] / setup.cx[0];
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, DivideLC ) {
    setup.exp->divide(setup.cu, setup.a, setup.cx);
    setup.cu2[0] = setup.a / setup.cx[0];
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, DivideRC ) {
    setup.exp->divide(setup.cu, setup.cx, setup.a);
    setup.cu2[0] = setup.cx[0] / setup.a;
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, PowConst ) {
    setup.exp->pow(setup.cu, setup.cx, setup.cx);
    setup.cu2[0] = std::pow(setup.cx[0], setup.cx[0]);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
				   setup.crtol, setup.catol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, PlusEqualC ) {
    setup.cu = setup.cx;
    setup.cu2[0] = setup.cu[0] + setup.a;
    setup.exp->plusEqual(setup.cu, setup.a);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, PlusEqualC2 ) {
    setup.cu = setup.cx;
    setup.cu2[0] = setup.cu[0] + setup.cx[0];
    setup.exp->plusEqual(setup.cu, setup.cx);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
  				   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, MinusEqualC ) {
    setup.cu = setup.cx;
    setup.cu2[0] = setup.cu[0] - setup.a;
    setup.exp->minusEqual(setup.cu, setup.a);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, MinusEqualC2 ) {
    setup.cu = setup.cx;
    setup.cu2[0] = setup.cu[0] - setup.cx[0];
    setup.exp->minusEqual(setup.cu, setup.cx);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
  				   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, TimesEqualC ) {
    setup.cu = setup.cx;
    setup.cu2[0] = setup.cu[0] * setup.a;
    setup.exp->timesEqual(setup.cu, setup.a);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, TimesEqualC2 ) {
    setup.cu = setup.cx;
    setup.cu2[0] = setup.cu[0] * setup.cx[0];
    setup.exp->timesEqual(setup.cu, setup.cx);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
  				   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, DivideEqualC ) {
    setup.cu = setup.cx;
    setup.cu2[0] = setup.cu[0] / setup.a;
    setup.exp->divideEqual(setup.cu, setup.a);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
  				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_ConstantExpansion, DivideEqualC2 ) {
    setup.cu = setup.cx;
    setup.cu2[0] = setup.cu[0] / setup.cx[0];
    setup.exp->divideEqual(setup.cu, setup.cx);
    success = Stokhos::comparePCEs(setup.cu, "cu", setup.cu2, "cu2", 
  				   setup.rtol, setup.atol, out);
  }

  // Not testing atan2(), max(), min(), abs(), fabs() since these are
  // not smooth functions

}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
