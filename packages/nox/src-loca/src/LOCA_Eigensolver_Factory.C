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

#include "NOX_Parameter_List.H"
#include "LOCA_Eigensolver_Factory.H"
#include "LOCA_Eigensolver_DefaultStrategy.H"
#include "LOCA_Eigensolver_AnasaziStrategy.H"
#include "LOCA_ErrorCheck.H"

Teuchos::RefCountPtr<LOCA::Eigensolver::AbstractStrategy>
LOCA::Eigensolver::Factory::createStrategy(
	     const Teuchos::RefCountPtr<NOX::Parameter::List>& eigenParams,
	     const Teuchos::RefCountPtr<NOX::Parameter::List>& solverParams)
{
  string methodName = "LOCA::Eigensolver::Factory::createStrategy()";
  Teuchos::RefCountPtr<LOCA::Eigensolver::AbstractStrategy> strategy;

  // Get name of strategy
  string name = eigenParams->getParameter("Method", "Default");

  if (name == "Default")
    strategy = 
      Teuchos::rcp(new LOCA::Eigensolver::DefaultStrategy(eigenParams,
							  solverParams));
  else if (name == "Anasazi")
    strategy = 
      Teuchos::rcp(new LOCA::Eigensolver::AnasaziStrategy(eigenParams,
							  solverParams));
  else if (name == "User-Defined") {

    // Get name of user-defined strategy
    string userDefinedName = eigenParams->getParameter("User-Defined Name",
						       "???");
    if ((*eigenParams).template 
	  isParameterRcp<LOCA::Eigensolver::AbstractStrategy>(userDefinedName))
      strategy = (*eigenParams).template 
	getRcpParameter<LOCA::Eigensolver::AbstractStrategy>(userDefinedName);
    else
       LOCA::ErrorCheck::throwError(methodName,
				    "Cannot find user-defined strategy: " + 
				    userDefinedName);
  }
  else
    LOCA::ErrorCheck::throwError(methodName,
				 "Invalid eigensolver strategy: " + 
				 name);

  return strategy;
}
