// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include <vector>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TestingHelpers.hpp"

#include "Sacado.hpp"

template <typename T>
struct container {
  T x;
};

// Test that checks ConditionalReturnType used in the return type deduction for
// conditional operations uses SFINAE correctly to remove the Fad overloads
// when the compiler searches for unrelated overloads of !=
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ConditionalReturnType, Fad, FAD )
{
  std::vector< container<FAD> > x(1);
  std::vector< container<FAD> > y(x);

  // Test is really to check whether the code compiles, so if it compiles,
  // it passes
  success = true;
}

const int global_fad_size = 10;
typedef Sacado::Fad::Exp::DFad<double> Fad_DFadType;
typedef Sacado::Fad::Exp::SFad<double,global_fad_size> Fad_SFadType;
typedef Sacado::Fad::Exp::SLFad<double,global_fad_size> Fad_SLFadType;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(ConditionalReturnType, Fad, Fad_DFadType)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(ConditionalReturnType, Fad, Fad_SFadType)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(ConditionalReturnType, Fad, Fad_SLFadType)

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
