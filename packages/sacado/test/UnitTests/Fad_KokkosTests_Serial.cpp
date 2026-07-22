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

#define SACADO_TEST_DFAD 1
#include "Fad_KokkosTests.hpp"

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
