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

NOX::StatusTest::NormF::NormF(double tolerance, NOX::Abstract::Vector::NormType ntype, ScaleType stype) :
  status(Unconverged),
  normType(ntype),
  scaleType(stype),
  toleranceType(Absolute),
  specifiedTolerance(tolerance),
  initialTolerance(1.0),
  trueTolerance(tolerance)
{
}

NOX::StatusTest::NormF::NormF(double tolerance, ScaleType stype) :
  status(Unconverged),
  normType(NOX::Abstract::Vector::TwoNorm),
  scaleType(stype),
  toleranceType(Absolute),
  specifiedTolerance(tolerance),
  initialTolerance(1.0),
  trueTolerance(tolerance)
{
}

NOX::StatusTest::NormF::NormF(NOX::Abstract::Group& initialGuess, double tolerance, 
			      NOX::Abstract::Vector::NormType ntype, 
			      NOX::StatusTest::NormF::ScaleType stype) :
  status(Unconverged),
  normType(ntype),
  scaleType(stype),
  toleranceType(Relative),
  specifiedTolerance(tolerance),
  initialTolerance(0.0),
  trueTolerance(0.0)
{
  relativeSetup(initialGuess);
}


NOX::StatusTest::NormF::NormF(NOX::Abstract::Group& initialGuess, double tolerance, ScaleType stype) :
  status(Unconverged),
  normType(NOX::Abstract::Vector::TwoNorm),
  scaleType(stype),
  toleranceType(Relative),
  specifiedTolerance(tolerance),
  initialTolerance(0.0),
  trueTolerance(0.0)
{
  relativeSetup(initialGuess);
}

NOX::StatusTest::NormF::~NormF()
{
}

double NOX::StatusTest::NormF::computeNorm(const NOX::Abstract::Group& grp)
{
  double norm;
  int n = grp.getX().length();

  switch (normType) 
  {
    
  case NOX::Abstract::Vector::TwoNorm:
    norm = grp.getNormF();
    if (scaleType == Scaled)
      norm /= sqrt(1.0 * n);
    break;

  default:
    norm = grp.getF().norm(normType);
    if (scaleType == Scaled)
      norm /= n;
    break;

  }

  return norm;
}


void NOX::StatusTest::NormF::relativeSetup(NOX::Abstract::Group& initialGuess)
{
  NOX::Abstract::Group::ReturnType rtype;
  rtype = initialGuess.computeF();
  if (rtype != NOX::Abstract::Group::Ok) 
  {
    cerr << "NOX::StatusTest::NormF::NormF - Unable to compute F" << endl;
    throw "NOX Error";
  }
    
  initialTolerance = computeNorm(initialGuess); 
  trueTolerance = specifiedTolerance / initialTolerance;
}


NOX::StatusTest::StatusType NOX::StatusTest::NormF::checkStatus(const NOX::Solver::Generic& problem)
{
  normF = computeNorm( problem.getSolutionGroup() );

  if (normF < trueTolerance)
    status = Converged;
  else
    status = Unconverged;

  return status;
}

NOX::StatusTest::StatusType NOX::StatusTest::NormF::getStatus() const
{
  return status;
}

ostream& NOX::StatusTest::NormF::print(ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  stream << "F-Norm = " << Utils::sciformat(normF,3);
  stream << " < " << Utils::sciformat(trueTolerance, 3);
  stream << "\n";

  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << setw(13) << " ";
  stream << "(";

  if (scaleType == Scaled)
    stream << "Length-Scaled";
  else
    stream << "Unscaled";

  stream << " ";

  if (normType == NOX::Abstract::Vector::TwoNorm)
    stream << "Two-Norm";
  else if (normType == NOX::Abstract::Vector::OneNorm)
    stream << "One-Norm";
  else if (normType == NOX::Abstract::Vector::MaxNorm)
    stream << "Max-Norm";
  
  stream << ", ";

  if (toleranceType == Absolute) 
    stream << "Absolute Tolerance";
  else 
    stream << "Relative Tolerance";

  stream << ")";

  stream << endl;

  return stream;
}


double NOX::StatusTest::NormF::getNormF() const
{
  return normF;
}

double NOX::StatusTest::NormF::getTrueTolerance() const
{
  return trueTolerance;
}

double NOX::StatusTest::NormF::getSpecifiedTolerance() const
{
  return specifiedTolerance;
}

double NOX::StatusTest::NormF::getInitialTolerance() const
{
  return initialTolerance;
}
