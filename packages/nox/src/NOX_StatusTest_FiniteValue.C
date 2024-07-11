// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_StatusTest_FiniteValue.H" // class definition
#include "NOX_Common.H"  // for std::string class
#include "NOX_Solver_Generic.H"
#include "NOX_Abstract_Group.H"

NOX::StatusTest::FiniteValue::
FiniteValue(VectorType v, NOX::Abstract::Vector::NormType n) :
  vectorType(v),
  vectorTypeLabel("?"),
  normType(n),
  normTypeLabel("?"),
  status(Unevaluated),
  result(-1),
  normValue(-1.0)
{
  // Set the vector type label for printing
  if (vectorType == FVector)
    vectorTypeLabel = "F";
  else
    vectorTypeLabel = "Solution";

  // Set the norm type label for printing
  if (normType == NOX::Abstract::Vector::TwoNorm)
    normTypeLabel = "Two-Norm";
  else if (normType == NOX::Abstract::Vector::OneNorm)
    normTypeLabel = "One-Norm";
  else
    normTypeLabel = "Max-Norm";

}

NOX::StatusTest::FiniteValue::~FiniteValue()
{
}

NOX::StatusTest::StatusType NOX::StatusTest::FiniteValue::
checkStatus(const Solver::Generic& problem,
        NOX::StatusTest::CheckType checkType)
{
  // Reset the check
  normValue = -1.0;
  const NOX::Abstract::Group& grp = problem.getSolutionGroup();

  switch (checkType)
  {
  case NOX::StatusTest::Complete:
  case NOX::StatusTest::Minimal:


    if (vectorType == FVector)
    {
      if (normType == NOX::Abstract::Vector::TwoNorm)
    normValue = grp.getNormF();  // More efficient than recomputing norm
      else
    normValue = grp.getF().norm(normType);
    }
    else
      normValue = grp.getX().norm(normType);

    result = finiteNumberTest(normValue);

    status = (result == 0) ? Unconverged : Failed;
    break;

  case NOX::StatusTest::None:
  default:
    result = 1;
    status = Unevaluated;
    break;
  }

  return status;
}

NOX::StatusTest::StatusType NOX::StatusTest::FiniteValue::getStatus() const
{
  return status;
}

std::ostream& NOX::StatusTest::FiniteValue::print(std::ostream& stream, int indent) const
{


  // Set the correct label for the check result
  std::string label = "Unknown";
  if (result == 0)
    label = "Finite";
  else if (result == -1)
    label = "NaN";
  else if (result == -2)
    label = "Infinite";

  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  stream << "Finite Number Check (" << normTypeLabel;
  stream << " " << vectorTypeLabel;
  stream << ") = ";
  stream << label;
  //stream << " (" << normValue << ")";
  stream << std::endl;

  return stream;
}

int NOX::StatusTest::FiniteValue::finiteNumberTest(double x) const
{
  if (NOX_isnan(x))
    return -1;

  if (NOX_isinf(x))
    return -2;

  return 0;
}

bool NOX::StatusTest::FiniteValue::NOX_isnan(double x) const
{
#if defined(HAVE_NAN_SUPPORT)
  if (isnan(x))
#elif defined(FINITE_VALUE_HAVE_STD_ISNAN)
  if (std::isnan(x))
#elif defined(FINITE_VALUE_HAVE_GLOBAL_ISNAN)
  if (isnan(x))
#else
  if (x != x)
#endif
    return true;

  return false;
}

bool NOX::StatusTest::FiniteValue::NOX_isinf(double x) const
{
#ifdef HAVE_INF_SUPPORT
  if (isinf(x))
#elif defined(FINITE_VALUE_HAVE_STD_ISINF)
  if (std::isinf(x))
#elif defined(FINITE_VALUE_HAVE_GLOBAL_ISINF)
  if (isinf(x))
#else
  // Use IEEE 754 definition: Inf * 0 = NaN
  double z = 0.0 * x;
  if (NOX_isnan(z))
#endif
    return true;

  return false;
}
