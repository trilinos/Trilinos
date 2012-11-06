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
