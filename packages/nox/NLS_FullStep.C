// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NLS_FullStep.H"

NLS_FullStep::NLS_FullStep(const NLS_ParameterList& params) :
  step(1.0)
{
  reset(params);
}

NLS_FullStep::~NLS_FullStep()
{

}

void NLS_FullStep::reset(const NLS_ParameterList& params)
{
  step = params.getParameter("Step", step);
}

bool NLS_FullStep::search(const NLS_Group& oldgrp, const NLS_Vector& dir, NLS_Group& newgrp) const
{
  newgrp.computeX(oldgrp, dir, step);
  newgrp.computeRHS();
  return true;
}

