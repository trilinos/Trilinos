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

#include "NOX_LineSearch_FullStep.H" // class definition

#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Parameter_List.H"
#include "NOX_Utils.H"


using namespace NOX;
using namespace NOX::LineSearch;

FullStep::FullStep(Parameter::List& params) 
{
  reset(params);
}

FullStep::~FullStep()
{

}

bool FullStep::reset(Parameter::List& params)
{
  NOX::Parameter::List& p = params.sublist("Full Step");
  fullstep = p.getParameter("Full Step", 1.0);
  return true;
}

bool FullStep::compute(Abstract::Group& grp, double& step, 
		       const Abstract::Vector& dir,
		       const Solver::Generic& s)
{
  step = fullstep;
  const Abstract::Group& oldGrp = s.getPreviousSolutionGroup();
  grp.computeX(oldGrp, dir, step);
  return true;
}

