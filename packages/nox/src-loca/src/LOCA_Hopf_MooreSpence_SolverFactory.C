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

#include "LOCA_Hopf_MooreSpence_SolverFactory.H"
#include "LOCA_Hopf_MooreSpence_SolverStrategy.H"
#include "LOCA_Hopf_MooreSpence_SalingerBordering.H"

LOCA::Hopf::MooreSpence::SolverFactory::SolverFactory(
	        const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data) : 
  globalData(global_data)
{
}

LOCA::Hopf::MooreSpence::SolverFactory::~SolverFactory()
{
}

Teuchos::RefCountPtr<LOCA::Hopf::MooreSpence::SolverStrategy>
LOCA::Hopf::MooreSpence::SolverFactory::create(
       const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RefCountPtr<Teuchos::ParameterList>& solverParams)
{
  string methodName = 
    "LOCA::Hopf::MooreSpence::SolverFactory::create()";
  Teuchos::RefCountPtr<LOCA::Hopf::MooreSpence::SolverStrategy> strategy;

  // Get name of strategy
  const string& name = strategyName(*solverParams);

  if (name == "Salinger Bordering")
    strategy = 
      Teuchos::rcp(new LOCA::Hopf::MooreSpence::SalingerBordering(
								globalData,
								topParams,
								solverParams));

  else if (name == "User-Defined") {

    // Get name of user-defined strategy
    string userDefinedName = solverParams->get("User-Defined Name", "???");
    if ((*solverParams).INVALID_TEMPLATE_QUALIFIER
	isType< Teuchos::RefCountPtr<LOCA::Hopf::MooreSpence::SolverStrategy> >(userDefinedName))
      strategy = (*solverParams).INVALID_TEMPLATE_QUALIFIER
	get< Teuchos::RefCountPtr<LOCA::Hopf::MooreSpence::SolverStrategy> >(userDefinedName);
    else
       globalData->locaErrorCheck->throwError(
				       methodName,
				       "Cannot find user-defined strategy: " + 
				       userDefinedName);
  }
  else
    globalData->locaErrorCheck->throwError(
				      methodName,
				      "Invalid bordered solver strategy: " + 
				      name);

  return strategy;
}

const string&
LOCA::Hopf::MooreSpence::SolverFactory::strategyName(
				  Teuchos::ParameterList& solverParams) const
{
  return solverParams.get("Solver Method", "Salinger Bordering");
}
