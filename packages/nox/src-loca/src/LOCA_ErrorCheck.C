// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "LOCA_ErrorCheck.H"     // class definition

#include "LOCA_Utils.H"          // for printing utilities

LOCA::ErrorCheck::ErrorCheck() 
{
 
}

LOCA::ErrorCheck::~ErrorCheck()
{

}


void LOCA::ErrorCheck::throwError(const string& callingFunction, 
				  const string& message, 
				  const string& throwLabel)
{
  if (LOCA::Utils::doPrint(LOCA::Utils::Error)) {
    cout << "************************" << "\n";
    cout << "ERROR: " << callingFunction << "\n";
    if (message != "") 
      cout << message << "\n";       
    cout << "************************" << endl;
  }
  throw (throwLabel.c_str());
  return;    
}

void LOCA::ErrorCheck::printWarning(const string& callingFunction, 
				    const string& message)
{
  if (LOCA::Utils::doPrint(LOCA::Utils::Warning)) {
    cout << "WARNING: " << callingFunction << " - ";
    if (message != "") 
      cout << message << endl;     
  }  
  return;    
}

void
LOCA::ErrorCheck::checkReturnType(
			     const NOX::Abstract::Group::ReturnType& status,
			     const string& callingFunction)
{
  if (status == NOX::Abstract::Group::Ok)
    return;
  else if (status == NOX::Abstract::Group::Failed || 
	   status == NOX::Abstract::Group::NotDefined ||
	   status == NOX::Abstract::Group::BadDependency)
    checkReturnType(status, LOCA::ErrorCheck::ThrowError, callingFunction);
  else if (status == NOX::Abstract::Group::NotConverged)
    checkReturnType(status, LOCA::ErrorCheck::PrintWarning, 
		       callingFunction);
  else 
    throwError("LOCA::ErrorCheck::checkReturnType", "Unknown status");
}

void LOCA::ErrorCheck::checkReturnType(
			       const NOX::Abstract::Group::ReturnType& status,
			       const ActionType& action,
			       const string& callingFunction,
			       const string& message)
{
  if (status != NOX::Abstract::Group::Ok) {
    
    if (action == ThrowError) {
      const string messageWithReturnType = message + "\n" + "Return Type = " + 
	                                   getReturnTypeString(status);

      throwError(callingFunction, messageWithReturnType);
    }
    else if (action == PrintWarning) {
      const string messageWithReturnType = message + "\n" + "Return Type = " + 
	                                   getReturnTypeString(status);

      printWarning(callingFunction, messageWithReturnType);
    }
    else {
      printWarning("LOCA::ErrorCheck::checkReturnType", 
		   "Unknown ActionType!");
    }

  }

  return;
}

NOX::Abstract::Group::ReturnType 
LOCA::ErrorCheck::combineReturnTypes(
			      const NOX::Abstract::Group::ReturnType& status1,
			      const NOX::Abstract::Group::ReturnType& status2)
{
  if (status1 == NOX::Abstract::Group::NotDefined ||
	   status2 == NOX::Abstract::Group::NotDefined)
    return NOX::Abstract::Group::NotDefined;
  else if (status1 == NOX::Abstract::Group::BadDependency ||
	   status2 == NOX::Abstract::Group::BadDependency)
    return NOX::Abstract::Group::BadDependency;
  else if (status1 == NOX::Abstract::Group::Failed ||
	   status2 == NOX::Abstract::Group::Failed)
    return NOX::Abstract::Group::Failed;
  else if (status1 == NOX::Abstract::Group::NotConverged ||
	   status2 == NOX::Abstract::Group::NotConverged)
    return NOX::Abstract::Group::NotConverged;
  else
    return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::ErrorCheck::combineAndCheckReturnTypes(
			      const NOX::Abstract::Group::ReturnType& status1,
			      const NOX::Abstract::Group::ReturnType& status2,
			      const string& callingFunction)
{
  NOX::Abstract::Group::ReturnType status3 = 
    combineReturnTypes(status1, status2);
  checkReturnType(status3, callingFunction);
  return status3;
}

string LOCA::ErrorCheck::getReturnTypeString(
			       NOX::Abstract::Group::ReturnType status)
{
  if (status == NOX::Abstract::Group::Ok)
    return "Ok";
  else if (status == NOX::Abstract::Group::NotDefined)
    return "NotDefined";
  else if (status == NOX::Abstract::Group::BadDependency)
    return "BadDependency";
  else if (status == NOX::Abstract::Group::NotConverged)
    return "NotConverged";
  else if (status == NOX::Abstract::Group::Failed)
    return "Failed";

  // If we got ot here the return type is bad
  return "<Unknown Return Type>";  
}
