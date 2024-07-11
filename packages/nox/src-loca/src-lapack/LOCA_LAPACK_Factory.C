// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ParameterList.hpp"

#include "LOCA_LAPACK_Factory.H"
#include "LOCA_Parameter_SublistParser.H"
#include "LOCA_BorderedSolver_LAPACKDirectSolve.H"
#include "LOCA_Eigensolver_DGGEVStrategy.H"

LOCA::LAPACK::Factory::Factory() :
  globalData()
{
}

LOCA::LAPACK::Factory::~Factory()
{
}

void
LOCA::LAPACK::Factory::init(
           const Teuchos::RCP<LOCA::GlobalData>& global_data)
{
  globalData = global_data;
}

bool
LOCA::LAPACK::Factory::createBorderedSolverStrategy(
       const std::string& strategyName,
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
       Teuchos::RCP<LOCA::BorderedSolver::AbstractStrategy>& strategy)
{
  // Instantiate DGGEV strategy if requested
  if (strategyName == "LAPACK Direct Solve") {
    strategy =
      Teuchos::rcp(new LOCA::BorderedSolver::LAPACKDirectSolve(globalData,
                                   topParams,
                                   solverParams));
    return true;
  }
  else
    return false;
}

bool
LOCA::LAPACK::Factory::createEigensolverStrategy(
         const std::string& strategyName,
     const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     const Teuchos::RCP<Teuchos::ParameterList>& eigenParams,
     Teuchos::RCP<LOCA::Eigensolver::AbstractStrategy>& strategy)
{
  // Instantiate DGGEV strategy if requested
  if (strategyName == "DGGEV") {
    strategy =
      Teuchos::rcp(new LOCA::Eigensolver::DGGEVStrategy(globalData,
                            topParams,
                            eigenParams));
    return true;
  }
  else
    return false;
}

