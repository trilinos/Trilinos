// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
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
