/*
// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
*/

#include "Teuchos_TestForException.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_GlobalMPISession.hpp"


void func_a()
{
  TEST_FOR_EXCEPTION(true, std::logic_error, "This is an exception I throw!");
}


void func_b()
{
  func_a();
}


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
  
  TestForException_setEnableStacktrace(enableStackTrace);

  bool success = true;
  
  try {
    func_b();
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success)
  return (success ? 0 : 1);
}
