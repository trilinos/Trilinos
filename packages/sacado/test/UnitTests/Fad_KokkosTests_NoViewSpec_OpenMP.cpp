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
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

// Disable view specializations
#define SACADO_DISABLE_FAD_VIEW_SPEC

#include "Fad_KokkosTests.hpp"

#include "Kokkos_Core.hpp"

// Instantiate tests for OpenMP device
#define SACADO_TEST_DFAD 1
using Kokkos::OpenMP;
VIEW_FAD_TESTS_D( OpenMP )

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  Kokkos::initialize(argc,argv);
  Kokkos::print_configuration(std::cout);

  int res = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  // Finalize OpenMP
  Kokkos::finalize();

  return res;
}
