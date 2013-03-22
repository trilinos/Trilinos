// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
  return;    
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
