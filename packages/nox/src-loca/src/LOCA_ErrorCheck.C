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
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
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


void LOCA::ErrorCheck::throwError(string callingFunction, 
				  string message, 
				  string throwLabel) const
{
  if (LOCA::Utils::doPrint(LOCA::Utils::Error)) {
    cout << "************************" << "\n";
    cout << "ERROR in " << callingFunction << "\n";
    if (message != "") 
      cout << message << "\n";       
    cout << "************************" << endl;
  }  
  throw throwLabel;
  return;    
}

void LOCA::ErrorCheck::printWarning(string callingFunction, 
				    string message) const
{
  if (LOCA::Utils::doPrint(LOCA::Utils::Warning)) {
    cout << "WARNING: " << callingFunction << " - ";
    if (message != "") 
      cout << message << endl;     
  }  
  return;    
}

void LOCA::ErrorCheck::checkReturnType(NOX::Abstract::Group::ReturnType status,
				       ActionType action,
				       string callingFunction,
				       string message) const
{
  if (status != NOX::Abstract::Group::Ok) {
    
    if (action == ThrowError) {
      string messageWithReturnType = message + "\n" + "Return Type = " + 
	                             getReturnTypeString(status);

      throwError(callingFunction, messageWithReturnType);
    }
    else if (action == PrintWarning) {
      string messageWithReturnType = message + "\n" + "Return Type = " + 
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

string LOCA::ErrorCheck::getReturnTypeString(NOX::Abstract::Group::ReturnType status) const
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
