// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NLS_IntervalHalving.H"

NLS_IntervalHalving::NLS_IntervalHalving(const NLS_ParameterList& params) :
  maxstep(1.0),
  minstep(1.0e-12)
{
  reset(params);
}

NLS_IntervalHalving::~NLS_IntervalHalving()
{

}

void NLS_IntervalHalving::reset(const NLS_ParameterList& params)
{
  maxstep = params.getParameter("Maximum Step", maxstep);
  minstep = params.getParameter("Minimum Step", minstep);
  params.unused();
}

bool NLS_IntervalHalving::search(const NLS_Group& oldgrp, const NLS_Vector& dir, NLS_Group& newgrp) const
{
  double oldf = oldgrp.getRHS().norm();
  double newf;
  double step = maxstep;

  do {

    if (step < minstep)
      return false;

    newgrp.computeX(oldgrp, dir, step);
    newf = newgrp.computeRHS().norm();
    step = step * 0.5;

  } while (oldf < newf);

  return true;
}

