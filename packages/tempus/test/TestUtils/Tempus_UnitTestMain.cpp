//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Tempus_UnitTestMainUtils.hpp"

int main(int argc, char* argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Tempus_Test::enableFPE(true);
  Teuchos::UnitTestRepository::setGloballyReduceTestResult(true);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
