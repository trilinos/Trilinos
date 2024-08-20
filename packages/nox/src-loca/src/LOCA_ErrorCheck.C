// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_ErrorCheck.H"     // class definition
#include "LOCA_GlobalData.H"
#include "NOX_Utils.H"          // for printing utilities
#include "Teuchos_Assert.hpp"

LOCA::ErrorCheck::ErrorCheck(
          const Teuchos::RCP<LOCA::GlobalData>& global_data) :
  globalData(global_data)
{
}

LOCA::ErrorCheck::~ErrorCheck()
{
}


void LOCA::ErrorCheck::throwError(const std::string& callingFunction,
                  const std::string& message,
                  const std::string& throwLabel)
{
  if (globalData->locaUtils->isPrintType(NOX::Utils::Error)) {
    globalData->locaUtils->err() << "************************" << "\n";
    globalData->locaUtils->err() << "ERROR: " << callingFunction << "\n";
    if (message != "")
      globalData->locaUtils->err() << message << "\n";
    globalData->locaUtils->err() << "************************" << std::endl;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
    "ERROR: " << callingFunction << "\n"
    << "ThrowLabel: " << throwLabel << "\n"
    << message << "\n");
}

void LOCA::ErrorCheck::printWarning(const std::string& callingFunction,
                    const std::string& message)
{
  if (globalData->locaUtils->isPrintType(NOX::Utils::Warning)) {
    globalData->locaUtils->out() << "WARNING: " << callingFunction << " - ";
    if (message != "")
      globalData->locaUtils->out() << message << std::endl;
  }
  return;
}

void
LOCA::ErrorCheck::checkReturnType(
                 const NOX::Abstract::Group::ReturnType& status,
                 const std::string& callingFunction)
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
                   const std::string& callingFunction,
                   const std::string& message)
{
  if (status != NOX::Abstract::Group::Ok) {

    if (action == ThrowError) {
      const std::string messageWithReturnType = message + "\n" + "Return Type = " +
                                       getReturnTypeString(status);

      throwError(callingFunction, messageWithReturnType);
    }
    else if (action == PrintWarning) {
      const std::string messageWithReturnType = message + "\n" + "Return Type = " +
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
                  const std::string& callingFunction)
{
  NOX::Abstract::Group::ReturnType status3 =
    combineReturnTypes(status1, status2);
  checkReturnType(status3, callingFunction);
  return status3;
}

std::string LOCA::ErrorCheck::getReturnTypeString(NOX::Abstract::Group::ReturnType status)
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
