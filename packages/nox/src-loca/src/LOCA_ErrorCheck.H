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

#ifndef LOCA_ERRORCHECK_H
#define LOCA_ERRORCHECK_H

#include "NOX_Common.H"          // for string
#include "NOX_Abstract_Group.H"  // for ReturnType

namespace LOCA { 
  
  //! \brief An Error checking algorithm for NOX/LOCA routines.  
  /*! This object will check the return types on objects and print a warning or throw an error if appropriate
  */

  class ErrorCheck  {

  public:

    //! Constructor
    ErrorCheck();

    //! Destructor
    ~ErrorCheck();

    //! Type for the action to be taken if the ReturnType is not Ok.
    enum ActionType {ThrowError, PrintWarning};

    //! Generic call to throw that prints info to the screen.
    void throwError(string callingFunction = "<Unknown Method>",
		    string message = "",
		    string throwLabel = "LOCA Error") const;

    //! Generic printing algorithm for sending warnings to the screen.
    void printWarning(string callingFunction = "<Unknown Method>",
		      string message = "") const;

    //! Checks the return type for the NOX::AbstractGroup and may throw an error or print a warning to the screen based on the ActionType requested.
    void checkReturnType(NOX::Abstract::Group::ReturnType status,
			 ActionType action,
			 string callingFunction = "<Unknown Method>",
			 string message = "") const;

  private:

    string getReturnTypeString(NOX::Abstract::Group::ReturnType status) const;

  };
} // namespace LOCA

#endif
