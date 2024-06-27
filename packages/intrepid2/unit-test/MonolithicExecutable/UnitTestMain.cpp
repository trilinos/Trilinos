// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   UnitTestMain.cpp
    \brief  Main for Teuchos unit tests.
 */

#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Teuchos_StackedTimer.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_DefaultComm.hpp"

#include "Kokkos_Core.hpp"

#include <fstream>

int main( int argc, char* argv[] )
{
  // Note that the dtor for GlobalMPISession will call Kokkos::finalize_all() but does not call Kokkos::initialize()...
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Kokkos::initialize(argc,argv);
  Teuchos::UnitTestRepository::setGloballyReduceTestResult(true);
  
  int result = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
  
  return result;
}
