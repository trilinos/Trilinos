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

#include "Teuchos_TestForException.hpp"


namespace { int throwNumber = 0; }


void TestForException_incrThrowNumber()
{
  ++throwNumber;
}


int TestForException_getThrowNumber()
{
  return throwNumber;
}


void TestForException_break( const std::string &errorMsg )
{
  int break_on_me;
  break_on_me = errorMsg.length(); // Use errMsg to avoid compiler warning.
  // Above is just some statement for the debugger to break on.  Note: now is
  // a good time to examine the stack trace and look at the error message in
  // 'errorMsg' to see what happened.  In GDB just type 'where' or you can go
  // up by typing 'up' and moving up in the stack trace to see where you are
  // and how you got to this point in the code where you are throwning this
  // exception!  Typing in a 'p errorMsg' will show you what the error message
  // is.  Also, you should consider adding a conditional breakpoint in this
  // function based on a specific value of 'throwNumber' if the exception you
  // want to examine is not the first exception thrown.
}
