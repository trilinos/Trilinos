// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TestingHelpers.hpp"

#include "Sacado.hpp"

const int N = 10;

// Check whether the safe_sqrt() function works as expected
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SafeSqrt, SafeSqrt, AD )
{
  typedef AD ad_type;

  success = true;

  // Check non-zero value
  ad_type x(N, 1.5);
  for (int i=0; i<N; ++i)
    x.fastAccessDx(i) = 2.0;
  ad_type y = safe_sqrt(x);
  ad_type z = sqrt(x);
  TEST_EQUALITY( y.val(), z.val() );
  for (int i=0; i<N; ++i)
    TEST_EQUALITY( y.dx(i), z.dx(i) );

  // Check zero value
  x.val() = 0.0;
  y = safe_sqrt(x);
  TEST_EQUALITY_CONST( y.val(), 0.0 );
  for (int i=0; i<N; ++i)
    TEST_EQUALITY_CONST( y.dx(i), 0.0 );

  // Check double
  double a = 1.5;
  double b = Sacado::safe_sqrt(a);
  double c = std::sqrt(a);
  TEST_EQUALITY( b, c );
}

typedef Sacado::Fad::DFad<double> Fad_DFadType;
typedef Sacado::Fad::SLFad<double,N> Fad_SLFadType;
typedef Sacado::Fad::SFad<double,N> Fad_SFadType;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( SafeSqrt, SafeSqrt, Fad_DFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( SafeSqrt, SafeSqrt, Fad_SLFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( SafeSqrt, SafeSqrt, Fad_SFadType )

typedef Sacado::ELRFad::DFad<double> ELRFad_DFadType;
typedef Sacado::ELRFad::SLFad<double,N> ELRFad_SLFadType;
typedef Sacado::ELRFad::SFad<double,N> ELRFad_SFadType;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( SafeSqrt, SafeSqrt, ELRFad_DFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( SafeSqrt, SafeSqrt, ELRFad_SLFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( SafeSqrt, SafeSqrt, ELRFad_SFadType )

typedef Sacado::CacheFad::DFad<double> CacheFad_DFadType;
typedef Sacado::CacheFad::SLFad<double,N> CacheFad_SLFadType;
typedef Sacado::CacheFad::SFad<double,N> CacheFad_SFadType;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( SafeSqrt, SafeSqrt, CacheFad_DFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( SafeSqrt, SafeSqrt, CacheFad_SLFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( SafeSqrt, SafeSqrt, CacheFad_SFadType )

typedef Sacado::ELRCacheFad::DFad<double> ELRCacheFad_DFadType;
typedef Sacado::ELRCacheFad::SLFad<double,N> ELRCacheFad_SLFadType;
typedef Sacado::ELRCacheFad::SFad<double,N> ELRCacheFad_SFadType;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( SafeSqrt, SafeSqrt, ELRCacheFad_DFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( SafeSqrt, SafeSqrt, ELRCacheFad_SLFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( SafeSqrt, SafeSqrt, ELRCacheFad_SFadType )

#if defined(SACADO_ENABLE_NEW_DESIGN) && !defined(SACADO_NEW_FAD_DESIGN_IS_DEFAULT)
typedef Sacado::Fad::Exp::DFad<double> ExpFad_DFadType;
typedef Sacado::Fad::Exp::SLFad<double,N> ExpFad_SLFadType;
typedef Sacado::Fad::Exp::SFad<double,N> ExpFad_SFadType;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( SafeSqrt, SafeSqrt, ExpFad_DFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( SafeSqrt, SafeSqrt, ExpFad_SLFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( SafeSqrt, SafeSqrt, ExpFad_SFadType )
#endif

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
