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

#include "NOX_Linesearch_Halving.H"

#include "NOX_Utils.H"		// for static doPrint function
#include <iomanip>		// for setw

using namespace NOX;
using namespace NOX::Linesearch;

Halving::Halving(const Parameter::List& params) 
{
  reset(params);
}

Halving::~Halving()
{

}

void Halving::reset(const Parameter::List& params)
{ 
  minstep = params.getParameter("Minimum Step", 1.0e-12);
  defaultstep = params.getParameter("Default Step", 1.0);
  recoverystep = params.getParameter("Recovery Step", defaultstep);
}

bool Halving::operator()(Abstract::Group& newgrp, double& step, 
			 const Abstract::Group& oldgrp, const Abstract::Vector& dir) 
{
  double oldf = oldgrp.getNormRHS();
  double newf;
  bool isfailed = false;

  step = defaultstep;
  newgrp.computeX(oldgrp, dir, step);
  newgrp.computeRHS();    
  newf = newgrp.getNormRHS();

  if (Utils::doPrint(1)) {
   cout << "\n" << Utils::fill(72) << "\n" << " -- Interval Halving Line Search -- \n";
  }
  while ((newf >= oldf) && (!isfailed)) {

    if (Utils::doPrint(1)) {
      cout << Utils::fill(5,' ') << "step = " << Utils::sci(step);
      cout << Utils::fill(1,' ') << "oldf = " << Utils::sci(oldf);
      cout << Utils::fill(1,' ') << "newf = " << Utils::sci(newf);
      cout << endl;
    }

    step = step * 0.5;

    if (step < minstep) {
      isfailed = true;
      step = recoverystep;
    }

    newgrp.computeX(oldgrp, dir, step);
    newgrp.computeRHS();    
    newf = newgrp.getNormRHS();
  } 

  if (Utils::doPrint(1)) {
    cout << Utils::fill(5,' ') << "step = " << Utils::sci(step);
    cout << Utils::fill(1,' ') << "oldf = " << Utils::sci(oldf);
    cout << Utils::fill(1,' ') << "newf = " << Utils::sci(newf);
    cout << endl;
    if (isfailed)
      cout << "--Linesearch Failed!--" << endl;
    else
      cout << "--Step Accepted!--" << endl;
    cout << Utils::fill(72) << "\n" << endl;
  }

  return (!isfailed);
}

