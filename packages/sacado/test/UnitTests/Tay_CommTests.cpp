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

#include "Sacado_No_Kokkos.hpp"
#include "Sacado_Tay_CacheTaylor.hpp"
#include "Tay_CommTests.hpp"

typedef int Ordinal;
Sacado::Random<double> rnd;
TAY_COMM_TESTS(Sacado::Tay::Taylor<double>, Taylor)
TAY_COMM_TESTS(Sacado::Tay::CacheTaylor<double>, CacheTaylor)

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
