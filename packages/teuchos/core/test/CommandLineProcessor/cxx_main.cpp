// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Version.hpp"

int main( int argc, char* argv[] )
{

  using Teuchos::CommandLineProcessor;

  bool verbose = true;  
  bool parse_successful = true;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // First create tests for a command line processor that doesn't throw exceptions.
  try {
    // Read options from the commandline
    CommandLineProcessor  clp(false, false); // Don't throw exceptions

    double rel_proc_speed = 1e-5; // Should 
    clp.setOption( "rel-proc-speed", &rel_proc_speed, "Relative processor speed (try around 1.0 for timing)." );

    int size = 1;
    clp.setOption( "size", &size, "Size of memory blocks created." );

#ifdef HAVE_TEUCHOS_LONG_LONG_INT
    long long longLongOption = 42;
    clp.setOption( "longLongOption", &longLongOption, "An option of type long long." );
#endif // HAVE_TEUCHOS_LONG_LONG_INT

    // Parse the current input, which should return succesful.
    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if (verbose)
      std::cout << "Test 1:  CommandLineProcessor - No exceptions - All extra options ignored: ";
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL )
    {
      parse_successful = false;
      if (verbose) std::cout << "FAILED" << std::endl;
    }
    else
      if (verbose) std::cout << "PASSED" << std::endl;

    // Add a new option that is required
    int num = 1;
    clp.setOption( "num", &num, "Number of memory blocks created (required option).", true );
   
    // Now parse with this new option (which should not be passed in on the command line) 
    parse_return = clp.parse(argc,argv);
    if (verbose)
      std::cout << "Test 2:  CommandLineProcessor - No exceptions - All extra options ignored - 1 required: ";
    if( parse_return != CommandLineProcessor::PARSE_ERROR )
    {
      parse_successful = false;
      if (verbose) std::cout << "FAILED" << std::endl;
    }
    else
      if (verbose) std::cout << "PASSED" << std::endl;

  }
  catch( ... ) {
    if(verbose)
      std::cerr << "*** Caught UNEXPECTED unknown exception\n";
    parse_successful = false;  // No exceptions should be thrown for this command line processor.
  }

  // Next create tests for a command line processor that does throw exceptions.
  // Read options from the commandline
  try {
    CommandLineProcessor  clp2(true, false); // Throw exceptions

    clp2.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not." );

    double rel_proc_speed = 1e-5; // Should 
    clp2.setOption( "rel-proc-speed", &rel_proc_speed, "Relative processor speed (try around 1.0 for timing)." );

    int size = 1;
    clp2.setOption( "size", &size, "Size of memory blocks created." );

    // Add a new option that is required
    int num = 1;
    clp2.setOption( "num", &num, "Number of memory blocks created (required option).", true );
  
    // Parse the argument line and see if we get an exception thrown 
    clp2.parse(argc,argv);
  }
  catch( CommandLineProcessor::ParseError &excpt ) {
    if(verbose)
      std::cout << "*** Caught EXPECTED standard exception : " << excpt.what() << std::endl
                << "Test 3:  CommandLineProcessor - Throw exceptions - All extra options ignored - 1 required: PASSED" << std::endl;
  }
  catch( ... ) {
    if(verbose)
      std::cout << "*** Caught UNEXPECTED unknown exception" << std::endl
                << "Test 3:  CommandLineProcessor - Throw exceptions - All extra options ignored - 1 required: FAILED" << std::endl;
    parse_successful = false;  // No exceptions should be thrown for this command line processor.
  }

  // Next create tests for a command line processor that doesn't throw exceptions, and doesn't recognize all options.
  try {
    CommandLineProcessor  clp3(false, true); // Don't recognize all options

    // Parse the current input, which should not be successful because the test is run with "--verbose" argument.
    CommandLineProcessor::EParseCommandLineReturn parse_return = clp3.parse(argc,argv);
    if (verbose)
      std::cout << "Test 4 :  CommandLineProcessor - No exceptions - Extra options not recognized: ";
    if( parse_return != CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION )
    {
      parse_successful = false;
      if (verbose) std::cout << "FAILED" << std::endl;
    }
    else
      if (verbose) std::cout << "PASSED" << std::endl;

    // Now add the verbose option back in and add a required option.
    clp3.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not." );

    int num = 1;
    clp3.setOption( "num", &num, "Number of memory blocks created (required option).", true );

    parse_return = clp3.parse(argc,argv);
    if (verbose)
      std::cout << "Test 5 :  CommandLineProcessor - No exceptions - Extra options not recognized - 1 required: ";
    if( parse_return != CommandLineProcessor::PARSE_ERROR )
    {
      parse_successful = false;
      if (verbose) std::cout << "FAILED" << std::endl;
    }
    else
      if (verbose) std::cout << "PASSED" << std::endl;
  } 
  catch( ... ) {
    if(verbose)
      std::cerr << "*** Caught UNEXPECTED unknown exception" << std::endl;
    parse_successful = false;  // No exceptions should be thrown for this command line processor.
  }
 
  // Next create tests for a command line processor that doesn't throw exceptions, and doesn't recognize all options.
  try {
    if (verbose)
      std::cout << "Test 6 :  CommandLineProcessor - Throw exceptions - Extra options not recognized: ";

    CommandLineProcessor  clp4(true, true); // Don't recognize all options AND throw exceptions (default mode)

    // Parse the current input, which should not be successful because the test is run with "--verbose" argument.
    clp4.parse(argc,argv);
  } 
  catch( CommandLineProcessor::UnrecognizedOption &excpt ) {
    if(verbose)
      std::cout << "*** Caught EXPECTED standard exception : " << excpt.what() << std::endl
                << "Test 6:  CommandLineProcessor - Throw exceptions - Extra options not recognized: PASSED" << std::endl;
  }
  catch( ... ) {
    if(verbose)
      std::cout << "*** Caught UNEXPECTED unknown exception" << std::endl
                << "Test 5:  CommandLineProcessor - Throw exceptions - Extra options not recognized: FAILED" << std::endl;
    parse_successful = false;  // No exceptions should be thrown for this command line processor.
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
