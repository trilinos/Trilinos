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
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_RCPNode.hpp"

#include "TestClasses.hpp"

int main(int argc, char* argv[])
{

  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  const Teuchos::RCP<Teuchos::FancyOStream> out
    = Teuchos::VerboseObjectBase::getDefaultOStream();

  bool printActiveRcpNodesOnExit = true;
  bool success = true;

  try {

    // Read in commandline arguments
    Teuchos::CommandLineProcessor clp;
    clp.setOption("print-active-rcp-nodes-on-exit", "no-print-active-rcp-nodes-on-exit",
      &printActiveRcpNodesOnExit);
    (void)clp.parse(argc, argv);

    // Create a circular reference that will result in dangling active
    // RCPNodes at shutdown.
    RCP<A> a = rcp(new A);
    RCP<C> c = rcp(new C);
    a->set_C(c);
    c->set_A(a);

  } // try
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

  if(success) {
    *out << "\nEnd Result: TEST PASSED" << std::endl;
  }

  if (!printActiveRcpNodesOnExit) {
    Teuchos::RCPNodeTracer::setPrintActiveRcpNodesOnExit(false);
  }
  // ABOVE: If not set to 'fasle', we leave it at the default which should be
  // true to make sure that the default is true!

  return ( success ? 0 : 1 );
}
