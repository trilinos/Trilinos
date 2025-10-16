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

#include "Kokkos_Core.hpp"
#include "Sacado.hpp"

// Tests special size alignment for SFad on Cuda is correct
TEUCHOS_UNIT_TEST(Kokkos_View_CornerCases, ViewOfUnmanagedView)
{
  using fad_t = Sacado::Fad::SFad<double, 3>;
  using inner_t = Kokkos::View<fad_t*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using outer_t = Kokkos::View<inner_t*, Kokkos::HostSpace>;

  Kokkos::View<fad_t*, Kokkos::HostSpace> data("data", 9, 4);
  TEUCHOS_TEST_EQUALITY(data.use_count(), 1, out, success);
  {
    outer_t A("A", 10);
 
    for(int i=0 ; i<10; i++) A(i) = data;
    TEUCHOS_TEST_EQUALITY(data.use_count(), 1, out, success);
  }
  TEUCHOS_TEST_EQUALITY(data.use_count(), 1, out, success);
  data = Kokkos::View<fad_t*, Kokkos::HostSpace>();
  TEUCHOS_TEST_EQUALITY(data.use_count(), 0, out, success);
}

TEUCHOS_UNIT_TEST(Kokkos_View_CornerCases, SubviewOfUnmanagedView)
{
  using fad_t = Sacado::Fad::SFad<double, 3>;
  using view_um_t = Kokkos::View<fad_t*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  Kokkos::View<fad_t*, Kokkos::HostSpace> data("data", 9, 4);
  TEUCHOS_TEST_EQUALITY(data.use_count(), 1, out, success);

  view_um_t um_data = data;
  TEUCHOS_TEST_EQUALITY(data.use_count(), 1, out, success);
  TEUCHOS_TEST_EQUALITY(data.data(), um_data.data(), out, success);
  TEUCHOS_TEST_EQUALITY(um_data.extent_int(0), 9, out, success);

  view_um_t sub = Kokkos::subview(um_data, Kokkos::pair(1,8));
  TEUCHOS_TEST_EQUALITY(sub.extent_int(0), 7, out, success);
}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Initialize serial
  Kokkos::initialize(argc,argv);

  int res = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  Kokkos::finalize();

  return res;
}
