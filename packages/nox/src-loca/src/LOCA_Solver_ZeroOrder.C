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

#include "LOCA_Solver_ZeroOrder.H"	// class definition
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Common.H"
#include "NOX_Utils.H"

using namespace LOCA;
using namespace LOCA::Solver;
using namespace NOX;

/* Some compilers (in particular the SGI and ASCI Red - TFLOP) 
 * fail to find the max and min function.  Therfore we redefine them 
 * here. 
 */ 
#ifdef max
#undef max
#endif

#define max(a,b) ((a)>(b)) ? (a) : (b);

#ifdef min
#undef min
#endif

#define min(a,b) ((a)<(b)) ? (a) : (b);

ZeroOrder::ZeroOrder(Abstract::Group& xgrp, Status::Test& t, const Parameter::List& p)
{
  solverPtr = new NOX::Solver::Manager(xgrp, t, p);
}

bool ZeroOrder::reset(Abstract::Group& xgrp, Status::Test& t, const Parameter::List& p) 
{
  return solverPtr->reset(xgrp, t, p);
}

ZeroOrder::~ZeroOrder() 
{
  delete solverPtr;
}

Status::StatusType ZeroOrder::getStatus()
{
  return solverPtr->getStatus();;
}

Status::StatusType ZeroOrder::iterate()
{
  return solverPtr->iterate();
}

Status::StatusType ZeroOrder::solve()
{
  return solverPtr->solve();
}

const Abstract::Group& ZeroOrder::getSolutionGroup() const
{
  return solverPtr->getSolutionGroup();
}

const Abstract::Group& ZeroOrder::getPreviousSolutionGroup() const
{
  return solverPtr->getPreviousSolutionGroup();
}

int ZeroOrder::getNumIterations() const
{
  return solverPtr->getNumIterations();
}

const Parameter::List& ZeroOrder::getOutputParameters() const
{
  return solverPtr->getOutputParameters();
}

bool ZeroOrder::getStepLength(double& stepLength)
{
  // How do we compute this?
  // Should this even be here?
  return true;
}
