// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#include "NOX_LineSearch_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Utils.H"
#include "NOX_Parameter_List.H"

using namespace NOX;
using namespace NOX::LineSearch;

Common::Common(const NOX::Utils& u, NOX::Parameter::List& lineSearchParams) :
  utils(u),
  paramsPtr(&lineSearchParams),
  vecPtr(NULL)
{
  reset(lineSearchParams);
}

Common::~Common()
{
  delete vecPtr;
}

bool Common::reset(NOX::Parameter::List& lineSearchParams)
{
  paramsPtr = &lineSearchParams;
  totalNumLineSearchCalls = 0;
  totalNumNonTrivialLineSearches = 0;
  totalNumFailedLineSearches = 0;
  totalNumIterations = 0;
  return true;
}

void Common::printStep(int n, double step, double oldf, double newf, const string s) const
{
  if (utils.isPrintProcessAndType(Utils::InnerIteration)) 
  {
    cout << setw(3) << n << ":";
    cout << Utils::fill(1,' ') << "step = " << utils.sciformat(step);
    cout << Utils::fill(1,' ') << "oldf = " << utils.sciformat(sqrt(2. * oldf));
    cout << Utils::fill(1,' ') << "newf = " << utils.sciformat(sqrt(2. * newf));
    if (!s.empty()) 
    {
      cout << " " << s << "\n";
      cout << Utils::fill(72);
    }
    cout << endl;
  }
}

double Common::computeSlope(const Abstract::Vector& dir, const Abstract::Group& grp) 
{
   if (grp.isGradient()) 
     return(dir.dot(grp.getGradient()));

  // Allocate space for vecPtr if necessary
  if (vecPtr == NULL) 
    vecPtr = dir.clone(ShapeCopy);

  // v = J * dir
  NOX::Abstract::Group::ReturnType status = grp.applyJacobian(dir,*vecPtr);
  
  if (status != NOX::Abstract::Group::Ok) 
  {
    cout << "NOX::LineSearch::Common::computeSlope -  Unable to apply Jacobian!" << endl;
    throw "NOX Error";
  }

  // Check that F exists
  if (!grp.isF()) 
  {
    cout << "NOX::LineSearch::Common::computeSlope - Invalid F" << endl;
    throw "NOX Error";
  }

  // Return <v, F> = F' * J * dir = <J'F, dir> = <g, dir>
  return(vecPtr->dot(grp.getF()));
}

double Common::computeSlopeWithOutJac(const Abstract::Vector& dir, const Abstract::Group& grp) 
{
  // Allocate space for vecPtr and grpPtr if necessary
  if (vecPtr == 0) 
    vecPtr = dir.clone(ShapeCopy);
  if (grpPtr == 0)
    grpPtr = grp.clone(ShapeCopy);

  // Check that F exists
  if (!grp.isF()) 
  {
    cout << "NOX::LineSearch::Common::computeSlope - Invalid F" << endl;
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
  *vecPtr = grp.getX();
  vecPtr->update(eta, dir, 1.0);

  // Compute the new F --> F(x + eta * dir)
  grpPtr->setX(*vecPtr);  
  grpPtr->computeF();

  // Compute Js = (F(x + eta * dir) - F(x))/eta
  *vecPtr = grpPtr->getF();
  vecPtr->update(-1.0, grp.getF(), 1.0);
  vecPtr->scale(1.0/eta);
  
  return(vecPtr->dot(grp.getF()));
}

bool Common::setCommonDataValues() 
{
  NOX::Parameter::List& outputList = paramsPtr->sublist("Output");
  outputList.setParameter("Total Number of Line Search Calls", totalNumLineSearchCalls);
  outputList.setParameter("Total Number of Non-trivial Line Searches", totalNumNonTrivialLineSearches);
  outputList.setParameter("Total Number of Failed Line Searches", totalNumFailedLineSearches);
  outputList.setParameter("Total Number of Line Search Inner Iterations", totalNumIterations);
  return true;
}
