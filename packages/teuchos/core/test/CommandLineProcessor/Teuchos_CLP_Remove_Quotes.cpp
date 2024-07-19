// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"


int main( int argc, char* argv[] )
{
  using Teuchos::CommandLineProcessor;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  bool verbose = true;
  bool parse_successful = true;

  std::string test_arg[4];
  std::string test_title[4];
  test_arg[0] = "--method=Name_of_Method";       test_title[0] = "without Quotes";
  test_arg[1] = "--method=\"Name_of_Method\"";   test_title[1] = "with Quotes";
  test_arg[2] = "--method=\"Name_of_Method";     test_title[2] = "with Leading Quote";
  test_arg[3] = "--method=Name_of_Method\"";     test_title[3] = "with Trailing Quote";

  for (int i=0; i<4; i++) {
    try {
      if (verbose)
        std::cout << "Test "<<i<<" :  CLP - remove_quotes() " << test_title[i] << std::endl;

      argv[1] = &(test_arg[i][0]);
      CommandLineProcessor CLP(true, true);  // Recognize all options AND throw exceptions

      std::string method_name = "";

      CLP.setOption("method", &method_name, "Name of Method");
      CLP.parse(argc, argv);

      if (verbose)
        std::cout << "Test "<<i<<" :  CLP - remove_quotes() " << test_title[i] << ": ";
      if (method_name != "Name_of_Method")
      {
        parse_successful = false;
        if (verbose) std::cout << "FAILED" << std::endl;
      }
      else
        if (verbose) std::cout << "PASSED" << std::endl;
    }
    catch( CommandLineProcessor::UnrecognizedOption &excpt ) {
      if(verbose)
        std::cout << "*** Caught EXPECTED standard exception : " << excpt.what() << std::endl
                  << "Test "<<i<<" :  CLP - remove_quotes() " << test_title[i] << ": PASSED" << std::endl;
    }
    catch( ... ) {
      if(verbose)
        std::cout << "*** Caught UNEXPECTED unknown exception" << std::endl
                  << "Test "<<i<<" :  CLP - remove_quotes() " << test_title[i] << ": FAILED" << std::endl;
      parse_successful = false;  // No exceptions should be thrown for this command line processor.
    }
  }

  // Return whether the command line processor tests passed.
  if (parse_successful) {
    std::cout << "End Result: TEST PASSED" << std::endl;
    return 0;
  }
  else {
    std::cout << "End Result: TEST FAILED" << std::endl;
    return 1;
  }
}
