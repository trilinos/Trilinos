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

#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_toString.hpp"


namespace {


// Make sure this is initialized whenever needed before main starts!
bool& showTestFailureLocationImpl()
{
  static bool showTFL = false;
  return showTFL;
}


} // namespace


const std::string
Teuchos::passfail_with_location(const bool result,
  const std::string &file, const int lineNumber)
{
  std::string rtn = passfail(result);
  if (!result && showTestFailureLocation()) {
    rtn += " ==> "+file+":"+toString(lineNumber);
  }
  return rtn;
}


void Teuchos::showTestFailureLocation(bool showTFL)
{
  showTestFailureLocationImpl() = showTFL;
}


bool Teuchos::showTestFailureLocation()
{
  return showTestFailureLocationImpl();
}
