// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_GlobalMPISession.hpp"
#include "Kokkos_Core.hpp"
#include "Teuchos_UnitTestRepository.hpp"

/**
 * Initialize both Teuchos MPI Comm and Kokkos.
 * Run the tests and return test suite result code.
 */
int main( int argc, char* argv[] )
{
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);

    // Initialize Kokkos
    Kokkos::initialize(argc,argv);

    // For unit tests to reduce the result when MPI is used
    Teuchos::UnitTestRepository::setGloballyReduceTestResult(true);

    // Run all registered unit tests
    int res = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

    // Finalize Kokkos
    Kokkos::finalize();

    return res;
}
