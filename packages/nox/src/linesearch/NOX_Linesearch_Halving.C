// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

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

