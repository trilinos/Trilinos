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

#include "Fad_KokkosTests.hpp"

#include "Kokkos_Core.hpp"

// Instantiate tests for Serial device
using Kokkos::Serial;
VIEW_FAD_TESTS_D( Serial )

// Add a unit test verifying something from Albany compiles
TEUCHOS_UNIT_TEST(Kokkos_View_Fad, DynRankMauroDeepCopy )
{
  Kokkos::DynRankView<Sacado::Fad::DFad<double>,Kokkos::Serial> v1(
    "v1", 3, 5);
  Kokkos::DynRankView<Sacado::Fad::DFad<double>,Kokkos::LayoutRight,Kokkos::HostSpace> v2("v2", 3 , 5);

  Kokkos::deep_copy(v1, v2);

  // We're just verifying this compiles
  success = true;
}

// Tests assignment between ViewFad and Fad of different type
TEUCHOS_UNIT_TEST(Kokkos_View_Fad, MixedViewFadTypes )
{
  const int StaticDim = 5;
  typedef Sacado::Fad::SFad<double,StaticDim> FadType1;
  typedef Sacado::Fad::SLFad<double,StaticDim> FadType2;
  typedef Kokkos::View<FadType1*,Kokkos::Serial> ViewType;

  const size_t num_rows = 11;
  const size_t fad_size = StaticDim;
  ViewType v("v", num_rows, fad_size+1);

  FadType2 x = 0.0;
  x = v(1);

  const size_t sz = x.size();
  TEUCHOS_TEST_EQUALITY(sz, fad_size, out, success);
}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Initialize serial
  Kokkos::initialize(argc,argv);

  int res = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  Kokkos::finalize();

  return res;
}
