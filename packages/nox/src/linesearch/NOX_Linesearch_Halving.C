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
}

bool Halving::operator()(Abstract::Group& newgrp, double& step, 
			 const Abstract::Group& oldgrp, const Abstract::Vector& dir) const
{
  double oldf = oldgrp.getNormRHS();
  double newf;
  int precision = Utils::precision;

  step = defaultstep;
  newgrp.computeX(oldgrp, dir, step);
  newgrp.computeRHS();    
  newf = newgrp.getNormRHS();

  if (Utils::doPrint(1)) {
    cout.setf(ios::scientific);
    cout.precision(precision);
    cout << "\n" << Utils::stars << "Interval Halving Line Search ->\n";
  }
  while (newf >= oldf) {

    if (Utils::doPrint(1)) {
      cout << " step = " << setw(precision + 6) << step
	   << " oldf = " << setw(precision + 6) << oldf
	   << " newf = " << setw(precision + 6) << newf
	   << endl;
    }

    step = step * 0.5;

    if (step < minstep)
      return false;

    newgrp.computeX(oldgrp, dir, step);
    newgrp.computeRHS();    
    newf = newgrp.getNormRHS();
  } 

  
  if (Utils::doPrint(1)) {
    cout << " step = " << setw(precision + 6) << step
	 << " oldf = " << setw(precision + 6) << oldf
	 << " newf = " << setw(precision + 6) << newf
	 << " (STEP ACCEPTED!)"
	 << "\n" 
	 << Utils::stars << endl;

    cout.unsetf(ios::scientific);
  }

  return true;
}

