// $Id$ 
// $Source$ 

// NOX: An Object-Oriented Nonlinear Solver Package
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NOX_Linesearch_Fullstep.H"

using namespace NOX;
using namespace NOX::Linesearch;

FullStep::FullStep(const Parameter::List& params) 
{
  reset(params);
}

FullStep::~FullStep()
{

}

void FullStep::reset(const Parameter::List& params)
{
  defaultstep = params.getParameter("Default Step", 1.0);
}

bool FullStep::operator()(Abstract::Group& newgrp, double& step, 
			  const Abstract::Group& oldgrp, const Abstract::Vector& dir) const
{
  step = defaultstep;
  newgrp.computeX(oldgrp, dir, step);
  return true;
}

