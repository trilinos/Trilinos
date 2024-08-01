// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_Assert.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "SimpleThrowFunctions.hpp"


int main(int argc, char* argv[])
{

  using Teuchos::CommandLineProcessor;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  CommandLineProcessor clp(false); // Don't throw exceptions

  bool enableStackTrace = true;

  clp.setOption("enable-stacktrace", "no-enable-stacktrace", &enableStackTrace,
    "Determine if stacktracing is shown or not on exception throw." );

  CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);

  if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) {
    std::cerr << "\nEnd Result: TEST FAILED" << std::endl;
    return parse_return;
  }

  Teuchos::TestForException_setEnableStacktrace(enableStackTrace);

  bool success = true;

  try {
    func_b();
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success)
  return (success ? 0 : 1);
}
