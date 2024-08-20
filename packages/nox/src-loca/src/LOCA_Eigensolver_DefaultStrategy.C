// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_Eigensolver_DefaultStrategy.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"


LOCA::Eigensolver::DefaultStrategy::DefaultStrategy(
    const Teuchos::RCP<LOCA::GlobalData>& global_data,
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& /* topParams */,
    const Teuchos::RCP<Teuchos::ParameterList>& /* eigenParams */) :
  globalData(global_data)
{
}

LOCA::Eigensolver::DefaultStrategy::~DefaultStrategy()
{
}

NOX::Abstract::Group::ReturnType
LOCA::Eigensolver::DefaultStrategy::computeEigenvalues(
         NOX::Abstract::Group& /* group */,
         Teuchos::RCP< std::vector<double> >& /* evals_r */,
         Teuchos::RCP< std::vector<double> >& /* evals_i */,
         Teuchos::RCP< NOX::Abstract::MultiVector >& /* evecs_r */,
             Teuchos::RCP< NOX::Abstract::MultiVector >& /* evecs_i */)
{
  // Print a warning that this eigensolver strategy doesn't do anything
  globalData->locaErrorCheck->printWarning(
   "LOCA::Eigensolver::DefaultStrategy::computeEigenvalues()",
   "\nThe default Eigensolver strategy does not compute eigenvalues.\nSet the \"Method\" parameter of the \"Eigensolver\" sublist to chose an \neigensolver method.");
  return NOX::Abstract::Group::Ok;
}
