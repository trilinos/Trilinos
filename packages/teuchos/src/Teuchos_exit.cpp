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
#include "Teuchos_TestForException.hpp"
#include "Teuchos_VerboseObject.hpp"

void Teuchos_exit_helper(
  const char file[],
  int line,
  const char msg[],
  int error_code
  )
{
  std::ostringstream omsg;
  omsg << file << ":" << line << ": error code = " << error_code;
  if(msg)
    omsg << ": " << msg;
  const std::string &omsgstr = omsg.str();
  TestForException_break(omsgstr); // Allows us to set a breakpoint!
#ifdef HAVE_TEUCHOS_C_EXCEPTIONS
  throw std::logic_error(omsgstr);
#else
  *Teuchos::VerboseObjectBase::getDefaultOStream()
    << omsg.str() << "\n" << std::flush;
  std::exit(error_code);
#endif
}
