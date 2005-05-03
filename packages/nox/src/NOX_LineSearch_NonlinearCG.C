//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifdef WITH_PRERELEASE

#include "NOX_LineSearch_NonlinearCG.H"

#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Parameter_List.H"
#include "NOX_Utils.H"

using namespace NOX;
using namespace NOX::LineSearch;

NonlinearCG::NonlinearCG(const NOX::Utils& u, Parameter::List& params) :
  utils(u),
  vecPtr(0),
  grpPtr(0)
{
  reset(params);
}

NonlinearCG::~NonlinearCG()
{

}

bool NonlinearCG::reset(Parameter::List& params)
{ 
  //NOX::Parameter::List& p = params.sublist("NonlinearCG");
  return true;
}

bool NonlinearCG::compute(Abstract::Group& newgrp, 
		     double& step, 
		     const Abstract::Vector& dir,
		     const Solver::Generic& s) 
{
  if (utils.isPrintProcessAndType(NOX::Utils::InnerIteration))
  {
    cout << "\n" << NOX::Utils::fill(72) << "\n" << "-- NonlinearCG Line Search -- \n";
  }

  const Abstract::Group& oldgrp = s.getPreviousSolutionGroup();

  // Perform single-step linesearch

  // Note that the following could be wrapped with a while loop to allow
  // iterations to be attempted 

  double numerator = oldgrp.getF().dot(dir);
  double denominator = computeDirectionalDerivative(dir, oldgrp).dot(dir);

  step = - numerator / denominator;
  newgrp.computeX(oldgrp, dir, step);
  newgrp.computeF(); 

  double checkOrthogonality = fabs( newgrp.getF().dot(dir) ); 

  if (utils.isPrintProcessAndType(Utils::InnerIteration)) {
    cout << setw(3) << "1" << ":";
    cout << " step = " << utils.sciformat(step);
    cout << " orth = " << utils.sciformat(checkOrthogonality);
    cout << "\n" << NOX::Utils::fill(72) << "\n" << endl;
  }
  
  return true;
}


NOX::Abstract::Vector& NonlinearCG::computeDirectionalDerivative(
                                const Abstract::Vector& dir,
                                const Abstract::Group& grp)
{
  // Allocate space for vecPtr and grpPtr if necessary
  if (vecPtr == 0)
    vecPtr = dir.clone(ShapeCopy);
  if (grpPtr == 0)
    grpPtr = grp.clone(ShapeCopy);

  // Check that F exists
  if (!grp.isF())
  {
    cout << "NOX::LineSearch::NonlinearCG::computeDirectionalDerivative "
         << "- Invalid F" << endl;
    throw "NOX Error";
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

#endif
