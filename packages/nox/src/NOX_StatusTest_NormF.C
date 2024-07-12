// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
    throw std::runtime_error("NOX Error");
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
  NOX::size_type n = grp.getX().length();

  switch (normType)
  {

  case NOX::Abstract::Vector::TwoNorm:
    norm = grp.getNormF();
    if (scaleType == Scaled)
      norm /= sqrt(1.0 * static_cast<double>(n));
    break;

  default:
    norm = grp.getF().norm(normType);
    if (scaleType == Scaled)
      norm /= static_cast<double>(n);
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
