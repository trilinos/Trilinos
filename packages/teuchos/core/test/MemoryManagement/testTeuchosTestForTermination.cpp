// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TestForException.hpp"

int main(int argc, char* argv[])
{

  using Teuchos::GlobalMPISession;
  GlobalMPISession mpiSession(&argc,&argv);

  Teuchos::CommandLineProcessor clp;
  int terminate_on_procid = 0;
  clp.setOption("terminate-on-procid", &terminate_on_procid);
  (void)clp.parse(argc, argv);

  TEUCHOS_TEST_FOR_TERMINATION(
    GlobalMPISession::getRank() == terminate_on_procid,
    "Bingo, we are terminating on procid == "
    "terminate_on_procid = "<<GlobalMPISession::getRank()<<"!"
    );

  return 1; // Will never be called!

}
