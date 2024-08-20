// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_LineSearch_NonlinearCG.H"

#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Utils.H"
#include "NOX_GlobalData.H"
#include "NOX_StatusTest_FiniteValue.H"

NOX::LineSearch::NonlinearCG::
NonlinearCG(const Teuchos::RCP<NOX::GlobalData>& gd,
        Teuchos::ParameterList& params) :
  finiteValueTester( Teuchos::rcp(new StatusTest::FiniteValue()) )
{
  reset(gd, params);
}

NOX::LineSearch::NonlinearCG::~NonlinearCG()
{

}

bool NOX::LineSearch::NonlinearCG::
reset(const Teuchos::RCP<NOX::GlobalData>& gd,
      Teuchos::ParameterList& /* params */)
{
  utils = gd->getUtils();
  //Teuchos::ParameterList& p = params.sublist("NonlinearCG");
  return true;
}

bool NOX::LineSearch::NonlinearCG::compute(Abstract::Group& newgrp,
                   double& step,
                   const Abstract::Vector& dir,
                   const Solver::Generic& s)
{
  if (utils->isPrintType(NOX::Utils::InnerIteration))
  {
    utils->out() << "\n" << NOX::Utils::fill(72) << "\n"
        << "-- NonlinearCG Line Search -- \n";
  }

  const Abstract::Group& oldgrp = s.getPreviousSolutionGroup();

  // Perform single-step linesearch

  // Note that the following could be wrapped with a while loop to allow
  // iterations to be attempted

  double numerator = oldgrp.getF().innerProduct(dir);
  double denominator = computeDirectionalDerivative(dir, oldgrp).innerProduct(dir);

  if( finiteValueTester->finiteNumberTest(step) )
  {
    utils->out() << "NOX::LineSearch::NonlinearCG::compute "
         << "- step value is NaN or Inf. " << std::endl;
    throw std::runtime_error("NOX Error");
  }

  step = - numerator / denominator;
  newgrp.computeX(oldgrp, dir, step);
  newgrp.computeF();

  double checkOrthogonality = fabs( newgrp.getF().innerProduct(dir) );

  if (utils->isPrintType(Utils::InnerIteration)) {
    utils->out() << std::setw(3) << "1" << ":";
    utils->out() << " step = " << utils->sciformat(step);
    utils->out() << " orth = " << utils->sciformat(checkOrthogonality);
    utils->out() << "\n" << NOX::Utils::fill(72) << "\n" << std::endl;
  }

  return true;
}


NOX::Abstract::Vector& NOX::LineSearch::NonlinearCG::
computeDirectionalDerivative(const Abstract::Vector& dir,
                 const Abstract::Group& grp)
{
  // Allocate space for vecPtr and grpPtr if necessary
  if (Teuchos::is_null(vecPtr))
    vecPtr = dir.clone(ShapeCopy);
  if (Teuchos::is_null(grpPtr))
    grpPtr = grp.clone(ShapeCopy);

  // Check that F exists
  if (!grp.isF())
  {
    utils->out() << "NOX::LineSearch::NonlinearCG::computeDirectionalDerivative "
         << "- Invalid F" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  // Compute the perturbation parameter
  double lambda = 1.0e-6;
  double denominator = dir.norm();

  // Don't divide by zero
  if (denominator == 0.0)
    denominator = 1.0;

  double eta = lambda * (lambda + grp.getX().norm() / denominator);

  // Don't divide by zero
  if (eta == 0.0)
    eta = 1.0e-6;

  // Perturb the solution vector
  vecPtr->update(eta, dir, 1.0, grp.getX(), 0.0);

  // Compute the new F --> F(x + eta * dir)
  grpPtr->setX(*vecPtr);
  grpPtr->computeF();

  // Compute Js = (F(x + eta * dir) - F(x))/eta
  vecPtr->update(-1.0/eta, grp.getF(), 1.0/eta, grpPtr->getF(), 0.0);

  return(*vecPtr);
}
