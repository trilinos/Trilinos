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

#include "NOX_LineSearch_Backtrack.H" // class definition

#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Parameter_List.H"
#include "NOX_Utils.H"

NOX::LineSearch::Backtrack::Backtrack(const NOX::Utils& u, NOX::Parameter::List& params) :
  utils(u)
{
  reset(params);
}

NOX::LineSearch::Backtrack::~Backtrack()
{

}

bool NOX::LineSearch::Backtrack::reset(NOX::Parameter::List& params)
{ 
  NOX::Parameter::List& p = params.sublist("Backtrack");

  minStep = p.getParameter("Minimum Step", 1.0e-12);
  defaultStep = p.getParameter("Default Step", 1.0);
  recoveryStep = p.getParameter("Recovery Step", defaultStep);
  maxIters = p.getParameter("Max Iters", 100);

  const string tmp = p.getParameter("Decrease Condition", "Max Norm");
  
  if (tmp == "Max Norm")
    normType = NOX::Abstract::Vector::MaxNorm;
  else if (tmp == "Two Norm")
    normType = NOX::Abstract::Vector::TwoNorm;
  else 
  {
    cout << "NOX::LineSearch::Backtrack::reset - Invalid choice \"" << tmp 
	 << "\" for \"Decrease Condition\"" << endl;
    throw "NOX Error";
  }

  return true;
}

double NOX::LineSearch::Backtrack::getNormF(const NOX::Abstract::Group& grp) const
{
  return (normType == NOX::Abstract::Vector::MaxNorm) ? 
    grp.getF().norm(normType) : grp.getNormF();
}

bool NOX::LineSearch::Backtrack::compute(NOX::Abstract::Group& grp, double& step, 
					 const NOX::Abstract::Vector& dir,
					 const NOX::Solver::Generic& s)
{
  const Abstract::Group& oldGrp = s.getPreviousSolutionGroup();
  double oldF = getNormF(oldGrp);
  double newF;
  bool isFailed = false;

  step = defaultStep;
  grp.computeX(oldGrp, dir, step);

  NOX::Abstract::Group::ReturnType rtype;

  rtype = grp.computeF();    
  if (rtype != NOX::Abstract::Group::Ok)
  {
    cerr << "NOX::LineSearch::BackTrack::compute - Unable to compute F" << endl;
    throw "NOX Error";
  }

  newF = getNormF(grp);
  int nIters = 1;

  if (utils.isPrintProcessAndType(Utils::InnerIteration)) 
  {
   cout << "\n" << Utils::fill(72) << "\n" << "-- Backtrack Line Search -- \n";
  }

  while ((newF >= oldF) && (!isFailed)) 
  {

    if (utils.isPrintProcessAndType(Utils::InnerIteration)) 
    {
      cout << setw(3) << nIters << ":";
      cout << " step = " << utils.sciformat(step);
      cout << " oldF = " << utils.sciformat(oldF);
      cout << " newF = " << utils.sciformat(newF);
      cout << endl;
    }

    nIters ++;
    step = step * 0.5;

    if ((step < minStep) || (nIters > maxIters)) 
    {
      isFailed = true;
      step = recoveryStep;
    }

    grp.computeX(oldGrp, dir, step);

    rtype = grp.computeF();    
    if (rtype != NOX::Abstract::Group::Ok)
    {
      cerr << "NOX::LineSearch::BackTrack::compute - Unable to compute F" << endl;
      throw "NOX Error";
    }

    newF = getNormF(grp);
  } 

  if (utils.isPrintProcessAndType(Utils::InnerIteration)) 
  {
    cout << setw(3) << nIters << ":";
    cout << " step = " << utils.sciformat(step);
    cout << " oldF = " << utils.sciformat(oldF);
    cout << " newF = " << utils.sciformat(newF);
    if (isFailed)
      cout << " (USING RECOVERY STEP!)" << endl;
    else
      cout << " (STEP ACCEPTED!)" << endl;
    cout << Utils::fill(72) << "\n" << endl;
  }

  return (!isFailed);
}

