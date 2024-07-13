// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Common.H"
#include "NOX_Direction_SteepestDescent.H" // class definition
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Utils.H"
#include "NOX_MeritFunction_Generic.H"
#include "NOX_GlobalData.H"

NOX::Direction::SteepestDescent::
SteepestDescent(const Teuchos::RCP<NOX::GlobalData>& gd,
        Teuchos::ParameterList& params)
{
  reset(gd, params);
}

NOX::Direction::SteepestDescent::~SteepestDescent()
{

}

bool NOX::Direction::SteepestDescent::
reset(const Teuchos::RCP<NOX::GlobalData>& gd,
      Teuchos::ParameterList& params)
{
  globalDataPtr = gd;
  utils = gd->getUtils();
  meritFuncPtr = gd->getMeritFunction();

  Teuchos::ParameterList& p = params.sublist("Steepest Descent");

  const std::string tmp = p.get("Scaling Type", "2-Norm");
  if (tmp == "2-Norm")
    scaleType = NOX::Direction::SteepestDescent::TwoNorm;
  else if (tmp == "F 2-Norm")
    scaleType = NOX::Direction::SteepestDescent::FunctionTwoNorm;
  else if (tmp == "Quadratic Model Min")
    scaleType = NOX::Direction::SteepestDescent::QuadMin;
  else if (tmp == "None")
    scaleType = NOX::Direction::SteepestDescent::None;
  else {
    utils->out() << "NOX::Direction::SteepestDescent::reset - Invalid choice "
         << "\"" << tmp << "\" for \"Scaling Type\"" << std::endl;
    throw std::runtime_error("NOX Error");
  }

 return true;
}

bool NOX::Direction::SteepestDescent::compute(Abstract::Vector& dir,
                 Abstract::Group& soln,
                 const Solver::Generic& /* solver */)
{
  NOX::Abstract::Group::ReturnType status;

  // Compute F at current solution
  status = soln.computeF();
  if (status != NOX::Abstract::Group::Ok)
    throwError("compute", "Unable to compute F");

  // Compute Jacobian at current solution
  status = soln.computeJacobian();
  if (status != NOX::Abstract::Group::Ok)
    throwError("compute", "Unable to compute Jacobian");

  // Scale
  switch (scaleType) {

  case NOX::Direction::SteepestDescent::TwoNorm:

    meritFuncPtr->computeGradient(soln, dir);
    dir.scale(-1.0/dir.norm());
    break;

  case NOX::Direction::SteepestDescent::FunctionTwoNorm:

    meritFuncPtr->computeGradient(soln, dir);
    dir.scale(-1.0/soln.getNormF());
    break;

  case NOX::Direction::SteepestDescent::QuadMin:

    meritFuncPtr->computeQuadraticMinimizer(soln, dir);

    break;

  case NOX::Direction::SteepestDescent::None:

    meritFuncPtr->computeGradient(soln, dir);
    dir.scale( -1.0 );
    break;

  default:

    throwError("compute", "Invalid scaleType");

  }

  return true;
}

bool NOX::Direction::SteepestDescent::compute(Abstract::Vector& dir,
                 Abstract::Group& soln,
                 const Solver::LineSearchBased& solver)
{
  return NOX::Direction::Generic::compute( dir, soln, solver );
}

void NOX::Direction::SteepestDescent::throwError(const std::string& functionName,
                         const std::string& errorMsg)
{
    if (utils->isPrintType(Utils::Error))
      utils->err() << "NOX::Direction::SteepestDescent::" << functionName
       << " - " << errorMsg << std::endl;
    throw std::runtime_error("NOX Error");
}
