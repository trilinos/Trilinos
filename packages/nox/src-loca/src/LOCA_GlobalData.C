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

#include "LOCA_GlobalData.H"
#include "NOX_Utils.H"
#include "Teuchos_ParameterList.hpp"
#include "LOCA_ErrorCheck.H"
#include "LOCA_Factory.H"
#include "LOCA_Parameter_SublistParser.H"

LOCA::GlobalData::GlobalData(
           const Teuchos::RCP<NOX::Utils>& loca_utils,
           const Teuchos::RCP<LOCA::ErrorCheck>& loca_error_check,
           const Teuchos::RCP<LOCA::Factory>& loca_factory) :
  locaUtils(loca_utils),
  locaErrorCheck(loca_error_check),
  locaFactory(loca_factory),
  parsedParams()
{
}

LOCA::GlobalData::~GlobalData()
{
}

Teuchos::RCP<LOCA::GlobalData>
LOCA::createGlobalData(
          const Teuchos::RCP<Teuchos::ParameterList>& paramList,
          const Teuchos::RCP<LOCA::Abstract::Factory>& userFactory)
{
  // Create a global data object with null data fields
  Teuchos::RCP<LOCA::GlobalData> globalData =
    Teuchos::rcp(new LOCA::GlobalData(Teuchos::null,
                      Teuchos::null,
                      Teuchos::null));

  // Create utils
  globalData->locaUtils =
    Teuchos::rcp(new NOX::Utils(paramList->sublist("NOX").sublist("Printing")));

  // Create error check
  globalData->locaErrorCheck =
    Teuchos::rcp(new LOCA::ErrorCheck(globalData));

  // Create factory
  if (userFactory != Teuchos::null)
    globalData->locaFactory = Teuchos::rcp(new LOCA::Factory(globalData,
                                 userFactory));
  else
    globalData->locaFactory = Teuchos::rcp(new LOCA::Factory(globalData));

  // Parse parameter list
  globalData->parsedParams =
    Teuchos::rcp(new Parameter::SublistParser(globalData));
  globalData->parsedParams->parseSublists(paramList);

  return globalData;
}

void
LOCA::destroyGlobalData(
            const Teuchos::RCP<LOCA::GlobalData>& globalData)
{
  globalData->locaUtils = Teuchos::null;
  globalData->locaErrorCheck = Teuchos::null;
  globalData->locaFactory = Teuchos::null;
  globalData->parsedParams = Teuchos::null;
}
