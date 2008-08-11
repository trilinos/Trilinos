// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "NOX_StatusTest_Stagnation.H" // class definition
#include "NOX_Common.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Abstract_Group.H"

NOX::StatusTest::Stagnation::Stagnation(int maxSteps_, double tolerance_) :
  maxSteps(maxSteps_),
  numSteps(0),
  lastIteration(-1),
  tolerance(tolerance_),
  convRate(1.0),
  status(NOX::StatusTest::Unevaluated)
{
    
}

NOX::StatusTest::Stagnation::~Stagnation()
{
}

NOX::StatusTest::StatusType 
NOX::StatusTest::Stagnation::
checkStatus(const Solver::Generic& problem,
	    NOX::StatusTest::CheckType checkType)
{
  status = Unconverged;

  // This test should ignore the checkType!  This test must be run
  // each iteration because it triggers after a set number of
  // iterations.

  // First time through we don't do anything
  int niters = problem.getNumIterations(); 
  if (niters == 0) {
    lastIteration = 0;
    numSteps = 0;
    return Unconverged;
  }

  // Make sure we have not already counted the last nonlinear iteration.
  // This protects against multiple calls to checkStatus() in between 
  // nonlinear iterations.
  bool isCounted = false;
  if (niters == lastIteration) {
    isCounted = true;
  }
  else
    lastIteration = niters;

  // Compute the convergence rate and set counter appropriately
  if (!isCounted) {

    convRate = problem.getSolutionGroup().getNormF() / 
               problem.getPreviousSolutionGroup().getNormF();
    
    if (convRate >= tolerance)
      numSteps ++;
    else
      numSteps = 0;
   
  }

  if (numSteps >= maxSteps)
    status = Failed;

  return status;
}

NOX::StatusTest::StatusType NOX::StatusTest::Stagnation::getStatus() const
{
  return status;
}

ostream& NOX::StatusTest::Stagnation::print(ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  stream << "Stagnation Count = " << numSteps << " < " << maxSteps << "\n";

  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << "             (convergence rate = " << convRate << ")";
  stream << endl;
 return stream;
}


int NOX::StatusTest::Stagnation::getMaxNumSteps() const
{
  return maxSteps;
}

int NOX::StatusTest::Stagnation::getCurrentNumSteps() const
{
  return numSteps;
}

double NOX::StatusTest::Stagnation::getTolerance() const
{
  return tolerance;
}

double NOX::StatusTest::Stagnation::getConvRate() const
{
  return convRate;
}
