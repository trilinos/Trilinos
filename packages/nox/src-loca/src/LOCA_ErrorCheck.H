// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_ERRORCHECK_H
#define LOCA_ERRORCHECK_H

#include "Teuchos_RCP.hpp"

#include "NOX_Common.H"          // for std::string
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
    void throwError(const std::string& callingFunction = "<Unknown Method>",
            const std::string& message = "",
            const std::string& throwLabel = "LOCA Error");

    //! Generic printing algorithm for sending warnings to the screen.
    void printWarning(const std::string& callingFunction = "<Unknown Method>",
              const std::string& message = "");

    //! Checks the supplied return type and performs an appropriate action
    /*!
     * This routine performs the following actions depending on the value
     * of \em status
     * <ul>
     * <li> NOX::Abstract::Group::Ok -- nothing
     * <li> NOX::Abstract::Group::Failed -- print message and throw a std::string
     * <li> NOX::Abstract::Group::NotDefined -- print message and throw a std::string
     * <li> NOX::Abstract::Group::BadDependency -- print message and throw a std::string
     * <li> NOX::Abstract::Group::NotConverged -- print a warning message
     * </ul>
     */
    void checkReturnType(
          const NOX::Abstract::Group::ReturnType& status,
          const std::string& callingFunction = std::string("<Unknown Method>"));

    /*!
     * \brief Checks the return type for the NOX::AbstractGroup and may
     * throw an error or print a warning to the screen based on the
     * ActionType requested.
     */
    void checkReturnType(
          const NOX::Abstract::Group::ReturnType& status,
          const ActionType& action,
          const std::string& callingFunction = std::string("<Unknown Method>"),
          const std::string& message = std::string(""));

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
          const std::string& callingFunction = std::string("<Unknown Method>"));

  protected:

    //! Converts the NOX::Abstract::Group::ReturnType to a std::string
    std::string getReturnTypeString(NOX::Abstract::Group::ReturnType status);

  protected:

    //! Global data
    Teuchos::RCP<LOCA::GlobalData> globalData;

  };
} // namespace LOCA

#endif
