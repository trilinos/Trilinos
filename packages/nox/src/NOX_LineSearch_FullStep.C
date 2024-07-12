// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_LineSearch_FullStep.H" // class definition

#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_GlobalData.H"

using namespace NOX;
using namespace NOX::LineSearch;

FullStep::FullStep(const Teuchos::RCP<NOX::GlobalData>& /* gd */,
           Teuchos::ParameterList& params)
{
  Teuchos::ParameterList& p = params.sublist("Full Step");
  stepSize = p.get("Full Step", 1.0);
}

FullStep::~FullStep()
{

}

bool FullStep::reset(const Teuchos::RCP<NOX::GlobalData>& /* gd */,
             Teuchos::ParameterList& params)
{
  Teuchos::ParameterList& p = params.sublist("Full Step");
  stepSize = p.get("Full Step", 1.0);
  return true;
}

bool FullStep::compute(Abstract::Group& grp, double& step,
               const Abstract::Vector& dir,
               const Solver::Generic& s)
{
  step = stepSize;
  const Abstract::Group& oldGrp = s.getPreviousSolutionGroup();
  grp.computeX(oldGrp, dir, step);
  return true;
}

