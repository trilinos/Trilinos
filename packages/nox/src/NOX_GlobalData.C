// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_GlobalData.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_MeritFunction_SumOfSquares.H"
#include "NOX_SolverStats.hpp"

NOX::GlobalData::
GlobalData(const Teuchos::RCP<Teuchos::ParameterList>& noxParams)
{ this->initialize(noxParams); }

NOX::GlobalData::
GlobalData(const Teuchos::RCP<NOX::Utils>& utils,
           const Teuchos::RCP<NOX::MeritFunction::Generic>& mf) :
  utilsPtr(utils),
  meritFunctionPtr(mf)
{
  solverStatsPtr = Teuchos::rcp(new NOX::SolverStats);
}

NOX::GlobalData::~GlobalData() {}

void NOX::GlobalData::
initialize(const Teuchos::RCP<Teuchos::ParameterList>& noxParams)
{
  paramListPtr = noxParams;
  utilsPtr = Teuchos::rcp(new NOX::Utils(noxParams->sublist("Printing")));
  if (is_null(solverStatsPtr))
    solverStatsPtr = Teuchos::rcp(new NOX::SolverStats);

  Teuchos::ParameterList& so = noxParams->sublist("Solver Options");

  if (so.isType< Teuchos::RCP<NOX::MeritFunction::Generic> >("User Defined Merit Function")) {
    meritFunctionPtr = so.get< Teuchos::RCP<NOX::MeritFunction::Generic> >("User Defined Merit Function");
  }

  // PL validator sets a default null RCP. If it is null, allocate
  // a concrete default.
  if (is_null(meritFunctionPtr))
    meritFunctionPtr = Teuchos::rcp(new NOX::MeritFunction::SumOfSquares(utilsPtr)); 
}

Teuchos::RCP<NOX::Utils> NOX::GlobalData::getUtils() const
{ return utilsPtr; }

Teuchos::RCP<NOX::MeritFunction::Generic>
NOX::GlobalData::getMeritFunction() const
{ return meritFunctionPtr; }

Teuchos::RCP<Teuchos::ParameterList>
NOX::GlobalData::getNoxParameterList() const
{ return paramListPtr; }

Teuchos::RCP<const NOX::SolverStats>
NOX::GlobalData::getSolverStatistics() const
{ return solverStatsPtr; }

Teuchos::RCP<NOX::SolverStats>
NOX::GlobalData::getNonConstSolverStatistics()
{ return solverStatsPtr; }
