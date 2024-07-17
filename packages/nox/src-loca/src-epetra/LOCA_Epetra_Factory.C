// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ParameterList.hpp"

#include "LOCA_Epetra_Factory.H"
#include "LOCA_Parameter_SublistParser.H"
#include "LOCA_BorderedSolver_EpetraHouseholder.H"
#include "LOCA_BorderedSolver_EpetraAugmented.H"
#ifdef HAVE_NOX_EPETRAEXT
#ifdef HAVE_MPI
#include "LOCA_Epetra_AnasaziOperator_Floquet.H"
#endif
#endif

LOCA::Epetra::Factory::Factory() :
  globalData()
{
}

LOCA::Epetra::Factory::~Factory()
{
}

void
LOCA::Epetra::Factory::init(
           const Teuchos::RCP<LOCA::GlobalData>& global_data)
{
  globalData = global_data;
}

bool
LOCA::Epetra::Factory::createBorderedSolverStrategy(
       const std::string& strategyName,
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
       Teuchos::RCP<LOCA::BorderedSolver::AbstractStrategy>& strategy)
{
  // Instantiate Householder strategy if requested
  if (strategyName == "Householder") {
    strategy =
      Teuchos::rcp(new LOCA::BorderedSolver::EpetraHouseholder(globalData,
                                   topParams,
                                   solverParams));
    return true;
  }
  // Instantiate augmented strategy if requested
  else if (strategyName == "Augmented") {
    strategy =
      Teuchos::rcp(new LOCA::BorderedSolver::EpetraAugmented(globalData,
                                 topParams,
                                 solverParams));
    return true;
  }
  else
    return false;
}

bool
LOCA::Epetra::Factory::createAnasaziOperatorStrategy(
       const std::string& strategyName,
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RCP<Teuchos::ParameterList>& eigenParams,
       const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
       const Teuchos::RCP<NOX::Abstract::Group>& grp,
       Teuchos::RCP<LOCA::AnasaziOperator::AbstractStrategy>& strategy)
{
#ifdef HAVE_NOX_EPETRAEXT
#ifdef HAVE_MPI
if (strategyName == "Floquet") {

    strategy =
      Teuchos::rcp(new LOCA::Epetra::AnasaziOperator::Floquet(globalData,
                                topParams, eigenParams, solverParams, grp));
    return true;
  }
 else
#endif
#endif
    return false;
}


