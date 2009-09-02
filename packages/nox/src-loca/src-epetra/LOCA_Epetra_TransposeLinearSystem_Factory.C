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

#include "Teuchos_ParameterList.hpp"

#include "NOX_Common.H" // for NOX_Config.h

#include "LOCA_Epetra_TransposeLinearSystem_Factory.H"
#include "LOCA_Epetra_TransposeLinearSystem_AbstractStrategy.H"
#include "LOCA_Epetra_TransposeLinearSystem_TransposePreconditioner.H"
#ifdef HAVE_NOX_EPETRAEXT
#include "LOCA_Epetra_TransposeLinearSystem_ExplicitTranspose.H"
#endif
#include "LOCA_Epetra_TransposeLinearSystem_LeftPreconditioning.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

LOCA::Epetra::TransposeLinearSystem::Factory::Factory(
		  const Teuchos::RCP<LOCA::GlobalData>& global_data) : 
  globalData(global_data)
{
}

LOCA::Epetra::TransposeLinearSystem::Factory::~Factory()
{
}

Teuchos::RCP<LOCA::Epetra::TransposeLinearSystem::AbstractStrategy>
LOCA::Epetra::TransposeLinearSystem::Factory::create(
		const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
		const Teuchos::RCP<NOX::Epetra::LinearSystem>& linsys)
{
  string methodName = "LOCA::Epetra::TransposeLinearSystem::Factory::create()";
  Teuchos::RCP<LOCA::Epetra::TransposeLinearSystem::AbstractStrategy> strategy;

  // Get name of strategy
  const string& name = strategyName(*solverParams);

  if (name == "Transpose Preconditioner")
    strategy = 
      Teuchos::rcp(new LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner(globalData, solverParams, linsys));

#ifdef HAVE_NOX_EPETRAEXT
  else if (name == "Explicit Transpose") {
    strategy = 
      Teuchos::rcp(new LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose(globalData, solverParams, linsys));
  }
#endif

  else if (name == "Left Preconditioning") {
    strategy = 
      Teuchos::rcp(new LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning(globalData, solverParams, linsys));
  }

  else if (name == "User-Defined") {

    // Get name of user-defined strategy
    string userDefinedName = solverParams->get("User-Defined Name",
							"???");
    if ((*solverParams).INVALID_TEMPLATE_QUALIFIER
	isType< Teuchos::RCP<LOCA::Epetra::TransposeLinearSystem::AbstractStrategy> >(userDefinedName))
      strategy = (*solverParams).INVALID_TEMPLATE_QUALIFIER
	get< Teuchos::RCP<LOCA::Epetra::TransposeLinearSystem::AbstractStrategy> >(userDefinedName);
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
LOCA::Epetra::TransposeLinearSystem::Factory::strategyName(
				  Teuchos::ParameterList& solverParams) const
{
  return solverParams.get("Transpose Solver Method", 
				   "Transpose Preconditioner");
}
