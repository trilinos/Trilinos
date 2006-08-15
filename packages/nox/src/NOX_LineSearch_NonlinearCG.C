//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "NOX_LineSearch_NonlinearCG.H"

#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Utils.H"
#include "NOX_GlobalData.H"

NOX::LineSearch::NonlinearCG::
NonlinearCG(const Teuchos::RefCountPtr<NOX::GlobalData>& gd,
	    Teuchos::ParameterList& params)
{
  reset(gd, params);
}

NOX::LineSearch::NonlinearCG::~NonlinearCG()
{

}

bool NOX::LineSearch::NonlinearCG::
reset(const Teuchos::RefCountPtr<NOX::GlobalData>& gd,
      Teuchos::ParameterList& params)
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
  double denominator = 
    computeDirectionalDerivative(dir, oldgrp).innerProduct(dir);

  step = - numerator / denominator;
  newgrp.computeX(oldgrp, dir, step);
  newgrp.computeF(); 

  double checkOrthogonality = fabs( newgrp.getF().innerProduct(dir) ); 

  if (utils->isPrintType(Utils::InnerIteration)) {
    utils->out() << setw(3) << "1" << ":";
    utils->out() << " step = " << utils->sciformat(step);
    utils->out() << " orth = " << utils->sciformat(checkOrthogonality);
    utils->out() << "\n" << NOX::Utils::fill(72) << "\n" << endl;
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
