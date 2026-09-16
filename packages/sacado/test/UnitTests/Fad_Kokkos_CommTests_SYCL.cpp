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

#include "Sacado.hpp"
#include "Fad_CommTests.hpp"

Sacado::Random<double> rnd;

typedef int Ordinal;
typedef Sacado::Fad::SLFad<double,10> Fad_SLFadType;
typedef Sacado::Fad::SFad<double,5> Fad_SFadType;

// DFad is omitted:  it needs device-side allocation, as for HIP.
FAD_KOKKOS_COMM_TESTS_SYCL(Fad_SLFadType, Fad_SLFad)
FAD_KOKKOS_COMM_TESTS_SYCL(Fad_SFadType, Fad_SFad)

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Initialize SYCL
  Kokkos::InitializationSettings init_args;
  init_args.set_device_id(0);
  Kokkos::initialize( init_args );
  Kokkos::print_configuration(std::cout);

  int ret = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  // Finalize SYCL
  Kokkos::finalize();

  return ret;
}
