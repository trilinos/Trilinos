// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ParameterList.hpp"

#include "LOCA_Tpetra_Factory.hpp"
#include "LOCA_Parameter_SublistParser.H"
#include "LOCA_BorderedSolver_TpetraHouseholder.hpp"

LOCA::Tpetra::Factory::Factory() :
  globalData()
{
}

LOCA::Tpetra::Factory::~Factory()
{
}

void
LOCA::Tpetra::Factory::init(
           const Teuchos::RCP<LOCA::GlobalData>& global_data)
{
  globalData = global_data;
}

bool
LOCA::Tpetra::Factory::createBorderedSolverStrategy(
       const std::string& strategyName,
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
       Teuchos::RCP<LOCA::BorderedSolver::AbstractStrategy>& strategy)
{
  // Instantiate Householder strategy if requested
  if (strategyName == "Householder") {
    strategy =
      Teuchos::rcp(new LOCA::BorderedSolver::TpetraHouseholder(globalData,
                                                               topParams,
                                                               solverParams));
    return true;
  }
  else
    return false;
}
