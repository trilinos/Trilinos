// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_LineSearch_Utils_Slope.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_GlobalData.H"

NOX::LineSearch::Utils::Slope::Slope() {}

NOX::LineSearch::Utils::Slope::
Slope(const Teuchos::RCP<NOX::GlobalData>& gd) :
  utils(*(gd->getUtils()))
{

}

NOX::LineSearch::Utils::Slope::~Slope()
{

}

void NOX::LineSearch::Utils::Slope::
reset(const Teuchos::RCP<NOX::GlobalData>& gd)
{
  utils = *(gd->getUtils());
}

double NOX::LineSearch::Utils::Slope::
computeSlope(const Abstract::Vector& dir, const Abstract::Group& grp)
{
   if (grp.isGradient())
     return(dir.innerProduct(grp.getGradient()));

  // Allocate space for vecPtr if necessary
   if (Teuchos::is_null(vecPtr))
     vecPtr = dir.clone(ShapeCopy);

  // v = J * dir
  NOX::Abstract::Group::ReturnType status = grp.applyJacobian(dir,*vecPtr);

  if (status != NOX::Abstract::Group::Ok)
  {
    utils.out() << "NOX::LineSearch::Utils::Slope::computeSlope -  Unable to apply Jacobian!" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  // Check that F exists
  if (!grp.isF())
  {
    utils.out() << "NOX::LineSearch::Utils::Slope::computeSlope - Invalid F" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  // Return <v, F> = F' * J * dir = <J'F, dir> = <g, dir>
  return(vecPtr->innerProduct(grp.getF()));
}

double NOX::LineSearch::Utils::Slope::
computeSlopeWithOutJac(const Abstract::Vector& dir,
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
    utils.out() << "NOX::LineSearch::Utils::Slope::computeSlope - Invalid F" << std::endl;
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

  return(vecPtr->innerProduct(grp.getF()));
}
