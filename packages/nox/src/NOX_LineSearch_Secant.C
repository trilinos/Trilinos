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
  maxiters = params.getParameter("Max Iters", 100);
  return true;
}

bool Secant::compute(Abstract::Group& newgrp, double& step, 
			 const Abstract::Group& oldgrp, const Abstract::Vector& dir) 
{

  double oldf = 0.5*oldgrp.getNormF()*oldgrp.getNormF();  
  double fmin = 10.0*oldf; // Allowable proximity to oldf for bestStep
  double newf = 0.;
  bool isfailed = false;

  alpha = 1.e-5; // Used in a backward difference approximation for
		 // initialization of the numerical hessian
  newgrp.computeX(oldgrp, dir, -alpha);
  newgrp.computeF(); // Assumed gradient direction for this linesearch
  newf = 0.5*newgrp.getNormF()*newgrp.getNormF();  

  double etaOld = dir.dot(newgrp.getF());
  double eta;

  int niters = 1;

  step = 0.;
  newgrp.computeX(oldgrp, dir, step);

  if (Utils::doPrint(Utils::InnerIteration)) {
    cout << "\n" << Utils::fill(72) << "\n" << "-- Secant Line Search -- \n";
  }

  while ((abs(alpha)>1.e-8) && (niters<=maxiters)) { 

    newgrp.computeF();
    newf = 0.5*newgrp.getNormF()*newgrp.getNormF();  

    if (Utils::doPrint(Utils::InnerIteration)) {
      cout << setw(3) << niters << ":";
      cout << " alpha = " << Utils::sci(alpha);
      cout << " step = " << Utils::sci(step);
      cout << " oldf = " << Utils::sci(sqrt(2.*oldf));
      cout << " newf = " << Utils::sci(sqrt(2.*newf));
      cout << endl;
    }

    eta = dir.dot(newgrp.getF());

    if(newf < fmin) {
      bestStep = step;
      fmin = newf;
    }

    alpha = alpha*(eta/(etaOld-eta));
    step += alpha;
    niters++;

    newgrp.computeX(oldgrp, dir, step);
    etaOld = eta;

    //   Bounds on step length could be enforced here

  } // end while loop


  newgrp.computeF();
  newf = 0.5*newgrp.getNormF()*newgrp.getNormF();  

  if ((newf < oldf) && (abs(step)>minstep)) {
    if (Utils::doPrint(Utils::InnerIteration)) {
        cout << setw(3) << niters << ":";
        cout << " alpha = " << Utils::sci(alpha);
        cout << " step = " << Utils::sci(step);
        cout << " oldf = " << Utils::sci(sqrt(2.*oldf));
        cout << " newf = " << Utils::sci(sqrt(2.*newf));
        cout << " (STEP ACCEPTED!)" << endl;
        cout << Utils::fill(72) << "\n" << endl;
    }
  }

  else {
    step = bestStep; // Could also use Recovery step here
    newgrp.computeX(oldgrp, dir, step);
    newgrp.computeF();
    newf = 0.5*newgrp.getNormF()*newgrp.getNormF();  
    if (Utils::doPrint(Utils::InnerIteration)) {
      cout << Utils::fill(5,' ') << "alpha = " << Utils::sci(alpha);
      cout << Utils::fill(5,' ') << "step = " << Utils::sci(step);
      cout << Utils::fill(1,' ') << "oldf = " << Utils::sci(sqrt(2.*oldf));
      cout << Utils::fill(1,' ') << "newf = " << Utils::sci(sqrt(2.*newf));
      cout << " (USING BEST CURRENT STEP!)" << endl;
      cout << Utils::fill(72) << "\n" << endl;
    }
    // isfailed = true; // Relaxed severity of failure for now
    return(!isfailed);
  }
  return true;
}

