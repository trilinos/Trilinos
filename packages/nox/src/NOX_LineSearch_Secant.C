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

#include "NOX_LineSearch_Secant.H"

#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Parameter_List.H"
#include "NOX_Utils.H"

using namespace NOX;
using namespace NOX::LineSearch;

Secant::Secant(Parameter::List& params) 
{
  reset(params);
}

Secant::~Secant()
{

}

bool Secant::reset(Parameter::List& params)
{ 
  minstep = params.getParameter("Minimum Step", 1.0e-12);
  defaultstep = params.getParameter("Default Step", 1.0);
  recoverystep = params.getParameter("Recovery Step", defaultstep);
  maxiters = params.getParameter("Max Iters", 20);
  return true;
}

bool Secant::compute(Abstract::Group& newgrp, 
		     double& step, 
		     const Abstract::Vector& dir,
		     const Solver::Generic& s) 
{

  const Abstract::Group& oldgrp = s.getPreviousSolutionGroup();

  double oldf = 0.5*oldgrp.getNormF()*oldgrp.getNormF();  
  double oldfprime = dir.dot(oldgrp.getF()); 

  step = 1.0; // Could use different, user specified initial step
  newgrp.computeX(oldgrp, dir, step);
  newgrp.computeF(); // Assumed gradient direction for this linesearch
  double newf = 0.5*newgrp.getNormF()*newgrp.getNormF();  

  int niters = 0;

  if (Utils::doPrint(Utils::InnerIteration)) {
    cout << "\n" << Utils::fill(72) << "\n" << "-- Secant Line Search -- \n";
    cout << setw(3) << niters << ":";
    cout << " step = " << Utils::sci(step);
    cout << " oldf = " << Utils::sci(sqrt(2.*oldf));
    cout << " newf = " << Utils::sci(sqrt(2.*newf));
    //cout << endl;
  }

  relStepChange = step; // tolerance is hard-coded for now, RH
  double oldstep = 0.;

  while ( niters<maxiters && relStepChange>1.e-6 ) { 

    niters++;

    oldstep = step;
    step = - step * oldfprime/(dir.dot(newgrp.getF())-oldfprime);
    newgrp.computeX(oldgrp, dir, step);
    newgrp.computeF();
    newf = 0.5*newgrp.getNormF()*newgrp.getNormF();  

    if (Utils::doPrint(Utils::InnerIteration)) {
      cout << endl;
      cout << setw(3) << niters << ":";
      cout << " step = " << Utils::sci(step);
      cout << " oldf = " << Utils::sci(sqrt(2.*oldf));
      cout << " newf = " << Utils::sci(sqrt(2.*newf));
      //cout << endl;
    }

    // Update relative change in step size used in convergence check
    if( fabs(oldstep)<minstep ) 
      relStepChange = fabs(step - oldstep);
    else 
      relStepChange = fabs(step - oldstep) / oldstep;

  } // end while loop


  // Now check acceptability of computed step

  if( step<minstep) { 
    step = recoverystep;
    newgrp.computeX(oldgrp, dir, step);
    newgrp.computeF(); 
    newf = 0.5*newgrp.getNormF()*newgrp.getNormF();  
    niters++;
    if (Utils::doPrint(Utils::InnerIteration)) {
      cout << endl;
      cout << setw(3) << niters << ":";
      cout << " step = " << Utils::sci(step);
      cout << " oldf = " << Utils::sci(sqrt(2.*oldf));
      cout << " newf = " << Utils::sci(sqrt(2.*newf));
      cout << " (USING RECOVERY STEP!)" << endl;
      cout << endl;
    }
  }
  else 
    if (Utils::doPrint(Utils::InnerIteration)) 
      cout << " (STEP ACCEPTED!)" << endl;

  return true;
}

