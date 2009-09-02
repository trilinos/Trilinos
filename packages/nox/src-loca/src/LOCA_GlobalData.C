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
