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

// Size used for all Fad types
const int global_fad_size = 10;

// Check whether the ternary operator works as expected
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Ternary, Ternary, AD )
{
  typedef AD ad_type;
  typedef typename ad_type::value_type value_type;

  success = true;

  ad_type x(global_fad_size, value_type(1.5));
  for (int i=0; i<global_fad_size; ++i)
    x.fastAccessDx(i) = 2.0;

  ad_type y = x > 0 ? -x : x;

  TEST_EQUALITY_CONST( y.val(), -1.5 );
  for (int i=0; i<global_fad_size; ++i)
    TEST_EQUALITY_CONST( y.dx(i), -2.0 );
}

typedef Sacado::Fad::DFad<double> Fad_DFadType;
typedef Sacado::Fad::SLFad<double,global_fad_size> Fad_SLFadType;
typedef Sacado::Fad::SFad<double,global_fad_size> Fad_SFadType;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Ternary, Ternary, Fad_DFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Ternary, Ternary, Fad_SLFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Ternary, Ternary, Fad_SFadType )

typedef Sacado::ELRFad::DFad<double> ELRFad_DFadType;
typedef Sacado::ELRFad::SLFad<double,global_fad_size> ELRFad_SLFadType;
typedef Sacado::ELRFad::SFad<double,global_fad_size> ELRFad_SFadType;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Ternary, Ternary, ELRFad_DFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Ternary, Ternary, ELRFad_SLFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Ternary, Ternary, ELRFad_SFadType )

typedef Sacado::CacheFad::DFad<double> CacheFad_DFadType;
typedef Sacado::CacheFad::SLFad<double,global_fad_size> CacheFad_SLFadType;
typedef Sacado::CacheFad::SFad<double,global_fad_size> CacheFad_SFadType;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Ternary, Ternary, CacheFad_DFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Ternary, Ternary, CacheFad_SLFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Ternary, Ternary, CacheFad_SFadType )

typedef Sacado::ELRCacheFad::DFad<double> ELRCacheFad_DFadType;
typedef Sacado::ELRCacheFad::SLFad<double,global_fad_size> ELRCacheFad_SLFadType;
typedef Sacado::ELRCacheFad::SFad<double,global_fad_size> ELRCacheFad_SFadType;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Ternary, Ternary, ELRCacheFad_DFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Ternary, Ternary, ELRCacheFad_SLFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Ternary, Ternary, ELRCacheFad_SFadType )

#if defined(SACADO_ENABLE_NEW_DESIGN) && !defined(SACADO_NEW_FAD_DESIGN_IS_DEFAULT)
typedef Sacado::Fad::Exp::DFad<double> ExpFad_DFadType;
typedef Sacado::Fad::Exp::SLFad<double,global_fad_size> ExpFad_SLFadType;
typedef Sacado::Fad::Exp::SFad<double,global_fad_size> ExpFad_SFadType;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Ternary, Ternary, ExpFad_DFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Ternary, Ternary, ExpFad_SLFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Ternary, Ternary, ExpFad_SFadType )
#endif

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
