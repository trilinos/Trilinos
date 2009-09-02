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
      dynamic_cast<C&>(a);
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
