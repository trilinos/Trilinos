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

#include "LOCA_Predictor_Secant.H"
#include "LOCA_Continuation_Group.H"

LOCA::Predictor::Secant::Secant(NOX::Parameter::List& params) 
{
}

LOCA::Predictor::Secant::~Secant()
{
}

NOX::Abstract::Group::ReturnType 
LOCA::Predictor::Secant::reset(NOX::Parameter::List& params) 
{
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::Predictor::Secant::compute(LOCA::Continuation::Group& prevGroup,
				 LOCA::Continuation::Group& curGroup,
				 LOCA::Continuation::Vector& result) 
{
  NOX::Abstract::Group::ReturnType res;

  res = curGroup.computeSecant();
  if (res != LOCA::Abstract::Group::Ok)
    return res;
  
  result = curGroup.getPredictorDirection();

  return res;
}
