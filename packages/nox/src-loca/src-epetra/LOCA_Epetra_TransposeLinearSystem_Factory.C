// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
  std::string methodName = "LOCA::Epetra::TransposeLinearSystem::Factory::create()";
  Teuchos::RCP<LOCA::Epetra::TransposeLinearSystem::AbstractStrategy> strategy;

  // Get name of strategy
  const std::string& name = strategyName(*solverParams);

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
    std::string userDefinedName = solverParams->get("User-Defined Name",
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

const std::string&
LOCA::Epetra::TransposeLinearSystem::Factory::strategyName(
                  Teuchos::ParameterList& solverParams) const
{
  return solverParams.get("Transpose Solver Method",
                   "Transpose Preconditioner");
}
