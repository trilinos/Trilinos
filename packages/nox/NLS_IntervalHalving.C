// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NLS_IntervalHalving.H"
#include "NLS_Utilities.H"	// for static doPrint function
#include <iomanip>		// for setw

NLS_IntervalHalving::NLS_IntervalHalving(const NLS_ParameterList& params) :
  minstep(1.0e-12)
{
  reset(params);
}

NLS_IntervalHalving::~NLS_IntervalHalving()
{

}

void NLS_IntervalHalving::reset(const NLS_ParameterList& params)
{
  minstep = params.getParameter("Minimum Step", minstep);
}

bool NLS_IntervalHalving::operator()(NLS_Group& newgrp, double& step, 
				     const NLS_Group& oldgrp, const NLS_Vector& dir) const
{
  double oldf = oldgrp.getNormRHS();
  double newf;
  int precision = NLS_Utilities::precision;

  newgrp.computeX(oldgrp, dir, step);
  newgrp.computeRHS();    
  newf = newgrp.getNormRHS();

  if (NLS_Utilities::doPrint(1)) {
    cout.setf(ios::scientific);
    cout.precision(precision);
    cout << "\n" << NLS_Utilities::stars << "Interval Halving Line Search ->\n";
  }
  while (newf >= oldf) {

    if (NLS_Utilities::doPrint(1)) {
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

  
  if (NLS_Utilities::doPrint(1)) {
    cout << " step = " << setw(precision + 6) << step
	 << " oldf = " << setw(precision + 6) << oldf
	 << " newf = " << setw(precision + 6) << newf
	 << " (STEP ACCEPTED!)"
	 << "\n" 
	 << NLS_Utilities::stars << endl;

    cout.unsetf(ios::scientific);
  }

  return true;
}

