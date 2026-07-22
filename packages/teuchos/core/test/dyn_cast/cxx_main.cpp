// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Version.hpp"

class A { public: virtual ~A(){} };
class B : public A { public: void f(bool verbose) { if(verbose) std::cout << "\nB::f() called!\n"; } };
class C : public A {};

int main( int argc, char* argv[] )
{

  using Teuchos::CommandLineProcessor;

  bool verbose = true;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  try {

    // Read options from the commandline
    CommandLineProcessor  clp(false); // Don't throw exceptions
    clp.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not." );
    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    if(verbose) std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;

    if(verbose) std::cout
      << "\n*******************************************"
      << "\n*** Basic test of Teuchos::dyn_cast<>() ***"
      << "\n*******************************************\n";
    B b;
    A &a = b;
    try {
      if(verbose) std::cout << "\nTrying: dynamic_cast<C&>(a); [Should throw a std::bad_cast std::exception with very bad error message]\n";
      C& c_from_an_a = dynamic_cast<C&>(a);
      (void) c_from_an_a; // forestall "unused variable" compiler warnings
    }
    catch( const std::bad_cast &e ) {
      if(verbose) std::cout << "\nCaught std::bad_cast std::exception e where e.what() = \"" << e.what() << "\"\n";
    }
    try {
      if(verbose) std::cout << "\nTrying: Teuchos::dyn_cast<C>(a); [Should throw a std::bad_cast std::exception with a very good error message]\n";
      Teuchos::dyn_cast<C>(a);
    }
    catch( const std::bad_cast &e ) {
      if(verbose) std::cout << "\nCaught std::bad_cast std::exception e where e.what() = \"" << e.what() << "\"\n";
    }
    if(verbose) std::cout << "\nTrying:  Teuchos::dyn_cast<B>(a).f(); [Should succeed and print \"B::f() called\"]\n";
    Teuchos::dyn_cast<B>(a).f(verbose);
    if(verbose) std::cout << "\nAll tests check out!\n";
  }
  catch( const std::exception &excpt ) {
    if(verbose)
      std::cerr << "*** Caught standard std::exception : " << excpt.what() << std::endl;
    return 1;
  }
  catch( ... ) {
    if(verbose)
      std::cerr << "*** Caught an unknown std::exception\n";
    return 1;
  }
        return 0;
}
