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

#include "NOX_StatusTest_NormF.H"
#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"

using namespace NOX;
using namespace NOX::StatusTest;

NormF::NormF(double tolerance, Abstract::Vector::NormType ntype, ScaleType stype) :
  status(Unconverged),
  normType(ntype),
  scaleType(stype),
  toleranceType(Absolute),
  specifiedTolerance(tolerance),
  initialTolerance(1.0),
  trueTolerance(tolerance)
{
}

NormF::NormF(double tolerance, ScaleType stype) :
  status(Unconverged),
  normType(NOX::Abstract::Vector::TwoNorm),
  scaleType(stype),
  toleranceType(Absolute),
  specifiedTolerance(tolerance),
  initialTolerance(1.0),
  trueTolerance(tolerance)
{
}

NormF::NormF(Abstract::Group& initialGuess, double tolerance, Abstract::Vector::NormType ntype, ScaleType stype) :
  status(Unconverged),
  normType(ntype),
  scaleType(stype),
  toleranceType(Relative),
  specifiedTolerance(tolerance),
  initialTolerance(0.0),
  trueTolerance(0.0)
{
  initialGuess.computeF();
  initialTolerance = initialGuess.getNormF();
  trueTolerance = specifiedTolerance / initialTolerance;
}


NormF::NormF(Abstract::Group& initialGuess, double tolerance, ScaleType stype) :
  status(Unconverged),
  normType(Abstract::Vector::TwoNorm),
  scaleType(stype),
  toleranceType(Relative),
  specifiedTolerance(tolerance),
  initialTolerance(0.0),
  trueTolerance(0.0)
{
  initialGuess.computeF();
  initialTolerance = initialGuess.getNormF();
  trueTolerance = specifiedTolerance / initialTolerance;
}

NormF::~NormF()
{
}

StatusType NormF::checkStatus(const Solver::Generic& problem)
{
  status = Unconverged;
  const Abstract::Group& grp = problem.getSolutionGroup();
  int n = grp.getX().length();

  switch (normType) {
    
  case NOX::Abstract::Vector::TwoNorm:
    normF = grp.getNormF();
    if (scaleType == Scaled)
      normF /= sqrt(1.0 * n);
    break;

  default:
    normF = grp.getF().norm(normType);
    if (scaleType == Scaled)
      normF /= n;
    break;

  }

  if (normF < trueTolerance)
    status = Converged;

  return status;
}

StatusType NormF::getStatus() const
{
  return status;
}

ostream& NormF::print(ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  if (toleranceType == Absolute) {
    stream << "Absolute F-Norm = " << Utils::sci(normF) << " < " 
	   << Utils::sci(trueTolerance);
  }
  else {
    stream << "Relative F-Norm = " << Utils::sci(normF) << " < " 
	   << Utils::sci(trueTolerance);
  }
  stream << endl;
  return stream;
}
