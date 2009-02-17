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

#include "Teuchos_exit.h"
#include "Teuchos_StandardCatchMacros.hpp"
#include "some_c_func.h"

int main( int argc, char* argv[] ) {

  bool result, success = true;

  result = true;
	try {
    std::cerr << "\nCall a C function that call TEUCHOS_EXIT(...) ...\n";
    some_c_func();
	}
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, result);
  if (result) success = false;
  
  result = true;
  try {
    std::cerr << "\nRaise an std::exception with TEUCHOS_MSG_EXIT(...) right here in C++ code with a message ...\n";
    TEUCHOS_MSG_EXIT("This std::exception is raised from C++ code!",1);
	}
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, result);
  if (result) success = false;
  
  if(success)
    std::cerr << "\nEnd Result: TEST PASSED" << std::endl;	
  
  return ( success ? 0 : 1 );
  
}
