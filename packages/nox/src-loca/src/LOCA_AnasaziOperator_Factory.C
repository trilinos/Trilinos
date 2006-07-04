// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//                  LOCA Continuation Algorithm Package
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov) or Eric Phipps
// (etphipp@sandia.gov), Sandia National Laboratories.
//
// ************************************************************************
//@HEADER

#include "Teuchos_ParameterList.hpp"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

#include "LOCA_AnasaziOperator_Factory.H"
#include "LOCA_AnasaziOperator_AbstractStrategy.H"
#include "LOCA_AnasaziOperator_JacobianInverse.H"
#include "LOCA_AnasaziOperator_ShiftInvert.H"
#include "LOCA_AnasaziOperator_Cayley.H"

LOCA::AnasaziOperator::Factory::Factory(
	        const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data) : 
  globalData(global_data)
{
}

LOCA::AnasaziOperator::Factory::~Factory()
{
}

Teuchos::RefCountPtr<LOCA::AnasaziOperator::AbstractStrategy>
LOCA::AnasaziOperator::Factory::create(
       const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RefCountPtr<Teuchos::ParameterList>& eigenParams,
       const Teuchos::RefCountPtr<Teuchos::ParameterList>& solverParams,
       const Teuchos::RefCountPtr<NOX::Abstract::Group>& grp)
{
  string methodName = "LOCA::AnasaziOperator::Factory::create()";
  Teuchos::RefCountPtr<LOCA::AnasaziOperator::AbstractStrategy> strategy;

  // Get name of strategy
  const string& name = strategyName(*eigenParams);

  if (name == "Jacobian Inverse")
    strategy = 
      Teuchos::rcp(new LOCA::AnasaziOperator::JacobianInverse(globalData,
							      topParams,
							      eigenParams,
							      solverParams,
							      grp));
  else if (name == "Shift-Invert") {
    Teuchos::RefCountPtr<LOCA::TimeDependent::AbstractGroup> tdGrp = 
      Teuchos::rcp_dynamic_cast<LOCA::TimeDependent::AbstractGroup>(grp);
    if (tdGrp == Teuchos::null)
      globalData->locaErrorCheck->throwError(
	methodName, 
	std::string("Group argument for Shift-Invert Anasazi operator ") + 
	std::string("strategy must be a LOCA::TimeDependent::AbstractGroup."));
    strategy = 
      Teuchos::rcp(new LOCA::AnasaziOperator::ShiftInvert(globalData,
							  topParams,
							  eigenParams,
							  solverParams,
							  tdGrp));
  }
  else if (name == "Cayley") {
    Teuchos::RefCountPtr<LOCA::TimeDependent::AbstractGroup> tdGrp = 
      Teuchos::rcp_dynamic_cast<LOCA::TimeDependent::AbstractGroup>(grp);
    if (tdGrp == Teuchos::null)
      globalData->locaErrorCheck->throwError(
	methodName, 
	std::string("Group argument for Shift-Invert Anasazi operator ") + 
	std::string("strategy must be a LOCA::TimeDependent::AbstractGroup."));
    strategy = 
      Teuchos::rcp(new LOCA::AnasaziOperator::Cayley(globalData,
						     topParams,
						     eigenParams,
						     solverParams,
						     tdGrp));
  }
  else if (name == "User-Defined") {

    // Get name of user-defined strategy
    string userDefinedName = 
      eigenParams->get("Operator User-Defined Name", "???");
    if ((*eigenParams).INVALID_TEMPLATE_QUALIFIER
	isType< Teuchos::RefCountPtr<LOCA::AnasaziOperator::AbstractStrategy> >(userDefinedName))
      strategy = (*eigenParams).INVALID_TEMPLATE_QUALIFIER
	get< Teuchos::RefCountPtr<LOCA::AnasaziOperator::AbstractStrategy> >(userDefinedName);
    else
       globalData->locaErrorCheck->throwError(
				       methodName,
				       "Cannot find user-defined strategy: " + 
				       userDefinedName);
  }
  else
    globalData->locaErrorCheck->throwError(
				      methodName,
				      "Invalid Anasazi operator strategy: " + 
				      name);

  return strategy;
}

const string&
LOCA::AnasaziOperator::Factory::strategyName(
				  Teuchos::ParameterList& eigenParams) const
{
  return eigenParams.get("Operator", "Jacobian Inverse");
}
