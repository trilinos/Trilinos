// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef LOCA_ERRORCHECK_H
#define LOCA_ERRORCHECK_H

#include "Teuchos_RCP.hpp"

#include "NOX_Common.H"          // for string
#include "NOX_Abstract_Group.H"  // for ReturnType

// forward declarations
namespace LOCA {
  class GlobalData;
}

namespace LOCA { 
  
  //! \brief An Error checking algorithm for NOX/LOCA routines.  
  /*! 
   * This object will check the return types on objects and print a warning 
   * or throw an error if appropriate
   */

  class ErrorCheck  {

  public:

    //! Constructor
    ErrorCheck(const Teuchos::RCP<LOCA::GlobalData>& global_data);

    //! Destructor
    ~ErrorCheck();

    //! Type for the action to be taken if the ReturnType is not Ok.
    enum ActionType {
      ThrowError,  //!< Throw an error
      PrintWarning //!< Print a warning
    };

    //! Generic call to throw that prints info to the screen.
    void throwError(const string& callingFunction = "<Unknown Method>",
		    const string& message = "",
		    const string& throwLabel = "LOCA Error");

    //! Generic printing algorithm for sending warnings to the screen.
    void printWarning(const string& callingFunction = "<Unknown Method>",
		      const string& message = "");

    //! Checks the supplied return type and performs an appropriate action
    /*!
     * This routine performs the following actions depending on the value
     * of \em status
     * <ul> 
     * <li> NOX::Abstract::Group::Ok -- nothing
     * <li> NOX::Abstract::Group::Failed -- print message and throw a string
     * <li> NOX::Abstract::Group::NotDefined -- print message and throw a string
     * <li> NOX::Abstract::Group::BadDependency -- print message and throw a string
     * <li> NOX::Abstract::Group::NotConverged -- print a warning message
     * </ul>
     */
    void checkReturnType(
		  const NOX::Abstract::Group::ReturnType& status,
		  const string& callingFunction = string("<Unknown Method>"));

    /*! 
     * \brief Checks the return type for the NOX::AbstractGroup and may 
     * throw an error or print a warning to the screen based on the 
     * ActionType requested.
     */
    void checkReturnType(
		  const NOX::Abstract::Group::ReturnType& status,
		  const ActionType& action,
		  const string& callingFunction = string("<Unknown Method>"),
		  const string& message = string(""));

    //! Combines two return types
    /*!
     * If either return type is NOX::Abstract::Group::NotDefined, returns 
     * NotDefined.  Otherwise if either is BadDependcy, returns BadDependency, 
     * if either is Failed, returns Failed, if either is NotConverged,
     * returns NotConverged, and otherwise returns Ok.
     */
    NOX::Abstract::Group::ReturnType 
    combineReturnTypes(const NOX::Abstract::Group::ReturnType& status1,
		       const NOX::Abstract::Group::ReturnType& status2);

    //! Combines two return types and checks the first
    /*!
     * First combines \em status1 and \em status2 using combineReturnTypes()
     * and checks the first using checkReturnType().  
     */
    NOX::Abstract::Group::ReturnType 
    combineAndCheckReturnTypes(
		  const NOX::Abstract::Group::ReturnType& status1,
		  const NOX::Abstract::Group::ReturnType& status2,
		  const string& callingFunction = string("<Unknown Method>"));

  protected:

    //! Converts the NOX::Abstract::Group::ReturnType to a string
    string getReturnTypeString(NOX::Abstract::Group::ReturnType status);

  protected:

    //! Global data
    Teuchos::RCP<LOCA::GlobalData> globalData;

  };
} // namespace LOCA

#endif
