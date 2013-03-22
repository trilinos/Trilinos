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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

#include "NOX_StatusTest_NormF.H"
#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"

NOX::StatusTest::NormF::
NormF(double tolerance, 
      NOX::Abstract::Vector::NormType ntype, ScaleType stype, 
      const NOX::Utils* u) :
  status(Unevaluated),
  normType(ntype),
  scaleType(stype),
  toleranceType(Absolute),
  specifiedTolerance(tolerance),
  initialTolerance(1.0),
  trueTolerance(tolerance),
  normF(0.0)
{
  if (u != NULL)
    utils = *u;
}

NOX::StatusTest::NormF::
NormF(double tolerance, ScaleType stype, 
      const NOX::Utils* u) :
  status(Unevaluated),
  normType(NOX::Abstract::Vector::TwoNorm),
  scaleType(stype),
  toleranceType(Absolute),
  specifiedTolerance(tolerance),
  initialTolerance(1.0),
  trueTolerance(tolerance),
  normF(0.0)
{
  if (u != NULL)
    utils = *u;
}

NOX::StatusTest::NormF::
NormF(NOX::Abstract::Group& initialGuess, double tolerance, 
      NOX::Abstract::Vector::NormType ntype, 
      NOX::StatusTest::NormF::ScaleType stype, 
      const NOX::Utils* u) :
  status(Unevaluated),
  normType(ntype),
  scaleType(stype),
  toleranceType(Relative),
  specifiedTolerance(tolerance),
  initialTolerance(0.0),
  trueTolerance(0.0),
  normF(0.0)
{
  if (u != NULL)
    utils = *u;

  relativeSetup(initialGuess);
}


NOX::StatusTest::NormF::
NormF(NOX::Abstract::Group& initialGuess, double tolerance, ScaleType stype, 
      const NOX::Utils* u) :
  status(Unevaluated),
  normType(NOX::Abstract::Vector::TwoNorm),
  scaleType(stype),
  toleranceType(Relative),
  specifiedTolerance(tolerance),
  initialTolerance(0.0),
  trueTolerance(0.0),
  normF(0.0)
{
  if (u != NULL)
    utils = *u;

  relativeSetup(initialGuess);
}

NOX::StatusTest::NormF::~NormF()
{
}

void NOX::StatusTest::NormF::relativeSetup(NOX::Abstract::Group& initialGuess)
{
  NOX::Abstract::Group::ReturnType rtype;
  rtype = initialGuess.computeF();
  if (rtype != NOX::Abstract::Group::Ok) 
  {
    utils.err() << "NOX::StatusTest::NormF::NormF - Unable to compute F" 
		<< std::endl;
    throw "NOX Error";
  }
    
  initialTolerance = computeNorm(initialGuess); 
  trueTolerance = specifiedTolerance * initialTolerance;
}

void NOX::StatusTest::NormF::reset(double tolerance)
{
  specifiedTolerance = tolerance;
  
  if (toleranceType == Absolute)
    trueTolerance = tolerance;
  else
    trueTolerance = specifiedTolerance * initialTolerance;
}

void NOX::StatusTest::NormF::reset(NOX::Abstract::Group& initialGuess,
				   double tolerance)
{
  specifiedTolerance = tolerance;
  relativeSetup(initialGuess);
}

double NOX::StatusTest::NormF::computeNorm(const NOX::Abstract::Group& grp)
{
  if (!grp.isF())
    return -1.0;

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


NOX::StatusTest::StatusType NOX::StatusTest::NormF::
checkStatus(const NOX::Solver::Generic& problem,
	    NOX::StatusTest::CheckType checkType)
{
  if (checkType == NOX::StatusTest::None)
  {
    normF = 0.0;
    status = Unevaluated;
  }
  else
  {
    normF = computeNorm( problem.getSolutionGroup() );
    status = ((normF != -1) && (normF < trueTolerance)) ? Converged : Unconverged;
  }

  return status;
}

NOX::StatusTest::StatusType NOX::StatusTest::NormF::getStatus() const
{
  return status;
}

std::ostream& NOX::StatusTest::NormF::print(std::ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  stream << "F-Norm = " << Utils::sciformat(normF,3);
  stream << " < " << Utils::sciformat(trueTolerance, 3);
  stream << "\n";

  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << std::setw(13) << " ";
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

  stream << std::endl;

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
