// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NLS_FullStep.H"

NLS_FullStep::NLS_FullStep(const NLS_ParameterList& params) 
{
  reset(params);
}

NLS_FullStep::~NLS_FullStep()
{

}

void NLS_FullStep::reset(const NLS_ParameterList& params)
{
}

bool NLS_FullStep::operator()(NLS_Group& newgrp, double& step, 
			      const NLS_Group& oldgrp, const NLS_Vector& dir) const
{
  newgrp.computeX(oldgrp, dir, step);
  return true;
}

