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
#include "NOX_Parameter_List.H"
#include "NOX_Utils.H"

using namespace NOX;
using namespace NOX::LineSearch;

Backtrack::Backtrack(Parameter::List& params) 
{
  reset(params);
}

Backtrack::~Backtrack()
{

}

bool Backtrack::reset(Parameter::List& params)
{ 
  minstep = params.getParameter("Minimum Step", 1.0e-12);
  defaultstep = params.getParameter("Default Step", 1.0);
  recoverystep = params.getParameter("Recovery Step", defaultstep);
  maxiters = params.getParameter("Max Iters", 100);

  const string tmp = params.getParameter("Decrease Condition", "Max Norm");
  
  if (tmp == "Max Norm")
    normtype = NOX::Abstract::Vector::INF;
  else if (tmp == "Two Norm")
    normtype = NOX::Abstract::Vector::TWO;
  else {
    cout << "NOX::LineSearch::Backtrack::reset - Invalid choice \"" << tmp 
	 << "\" for \"Decrease Condition\"" << endl;
    throw "NOX Error";
  }

  return true;
}

double Backtrack::getNormRHS(const Abstract::Group& grp) const
{
  return (normtype == NOX::Abstract::Vector::INF) ? 
    grp.getRHS().norm(normtype) : grp.getNormRHS();
}

bool Backtrack::operator()(Abstract::Group& newgrp, double& step, 
			 const Abstract::Group& oldgrp, const Abstract::Vector& dir) 
{
  double oldf = getNormRHS(oldgrp);
  double newf;
  bool isfailed = false;

  step = defaultstep;
  newgrp.computeX(oldgrp, dir, step);
  newgrp.computeRHS();    
  newf = getNormRHS(newgrp);
  int niters = 1;

  if (Utils::doPrint(Utils::InnerIteration)) {
   cout << "\n" << Utils::fill(72) << "\n" << "-- Backtrack Line Search -- \n";
  }
  while ((newf >= oldf) && (!isfailed)) {

    if (Utils::doPrint(Utils::InnerIteration)) {
      cout << setw(3) << niters << ":";
      cout << " step = " << Utils::sci(step);
      cout << " oldf = " << Utils::sci(oldf);
      cout << " newf = " << Utils::sci(newf);
      cout << endl;
    }

    niters ++;
    step = step * 0.5;

    if ((step < minstep) || (niters > maxiters)) {
      isfailed = true;
      step = recoverystep;
    }

    newgrp.computeX(oldgrp, dir, step);
    newgrp.computeRHS();    
    newf = getNormRHS(newgrp);
  } 

  if (Utils::doPrint(Utils::InnerIteration)) {
    cout << setw(3) << niters << ":";
    cout << " step = " << Utils::sci(step);
    cout << " oldf = " << Utils::sci(oldf);
    cout << " newf = " << Utils::sci(newf);
    if (isfailed)
      cout << " (USING RECOVERY STEP!)" << endl;
    else
      cout << " (STEP ACCEPTED!)" << endl;
    cout << Utils::fill(72) << "\n" << endl;
  }

  return (!isfailed);
}

