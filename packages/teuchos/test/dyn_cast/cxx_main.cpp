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
#include "Teuchos_Version.hpp"

class A { public: virtual ~A(){} };
class B : public A { public: void f() { std::cout << "\nB::f() called!\n"; } };
class C : public A {};

int main( int argc, char* argv[] ) {

  std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;

  std::cout
    << "\n*******************************************"
    << "\n*** Basic test of Teuchos::dyn_cast<>() ***"
    << "\n*******************************************\n";
  B b;
  A &a = b;
  try {
    std::cout << "\nTrying: dynamic_cast<C&>(a); [Should throw a std::bad_cast exception with very bad error message]\n";
    dynamic_cast<C&>(a);
  }
  catch( const std::bad_cast &e ) {
    std::cout << "\nCaught std::bad_cast exception e where e.what() = \"" << e.what() << "\"\n";
  }
  try {
    std::cout << "\nTrying: Teuchos::dyn_cast<C>(a); [Should throw a std::bad_cast exception with a very good error message]\n";
    Teuchos::dyn_cast<C>(a);
  }
  catch( const std::bad_cast &e ) {
    std::cout << "\nCaught std::bad_cast exception e where e.what() = \"" << e.what() << "\"\n";
  }
	std::cout << "\nTrying:  Teuchos::dyn_cast<B>(a).f(); [Should succeed and print \"B::f() called\"]\n";
	Teuchos::dyn_cast<B>(a).f();
	std::cout << "\nAll tests check out!\n";
	return 0;
}
