// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Common.H"
#include "NOX_StatusTest_RelativeNormF.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"
#include "Teuchos_Assert.hpp"
#include "Teuchos_as.hpp"

NOX::StatusTest::RelativeNormF::
RelativeNormF(double in_tolerance, bool in_scale_norms_by_vector_length,
              const NOX::Utils* u, NOX::Abstract::Vector::NormType ntype) :
  status(NOX::StatusTest::Unevaluated),
  tolerance(in_tolerance),
  normF_0(0.0),
  normF(0.0),
  scale_norms_by_vector_length(in_scale_norms_by_vector_length),
  norm_type(ntype)
{
  if (u != NULL)
    utils = *u;
}

NOX::StatusTest::StatusType NOX::StatusTest::RelativeNormF::
checkStatus(const NOX::Solver::Generic& problem,
            NOX::StatusTest::CheckType checkType)
{
  // On initial iteration, compute initial norm F
  if (problem.getNumIterations() == 0) {
    normF_0 = problem.getSolutionGroup().getF().norm(norm_type);

    if (scale_norms_by_vector_length)
      normF_0 /= std::sqrt(Teuchos::as<double>(problem.getSolutionGroup().getF().length()));
  }

  if (checkType == NOX::StatusTest::None)
  {
    normF = -1.0;
    status = Unevaluated;
  }
  else
    {
      normF = problem.getSolutionGroup().getF().norm(norm_type);

      if (scale_norms_by_vector_length)
        normF /= std::sqrt(Teuchos::as<double>(problem.getSolutionGroup().getF().length()));

      status = (normF < tolerance * normF_0) ? Converged : Unconverged;
  }

  return status;
}

NOX::StatusTest::StatusType NOX::StatusTest::RelativeNormF::getStatus() const
{
  return status;
}

std::ostream& NOX::StatusTest::RelativeNormF::print(std::ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  stream << "F-Norm = " << Utils::sciformat(normF,3);
  stream << " < " << Utils::sciformat(tolerance * normF_0, 3);
  stream << " (" << Utils::sciformat(tolerance, 3);
  stream << " * " << Utils::sciformat(normF_0, 3);
  stream << ")" << std::endl;

  return stream;
}
