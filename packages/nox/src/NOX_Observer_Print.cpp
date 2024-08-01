// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Observer_Print.hpp"
#include "NOX_Solver_Generic.H"
#include "NOX_Abstract_Group.H"
#include "NOX_SolverStats.hpp"
#include "NOX_Utils.H"
#include "NOX_Solver_LineSearchBased.H"

NOX::ObserverPrint::ObserverPrint(const Teuchos::RCP<NOX::Utils>& os) :
  os_(os)
{}

void NOX::ObserverPrint::runPreIterate(const NOX::Solver::Generic& solver)
{
  if (solver.getNumIterations() == 0) {
    auto& os = os_->out();
    auto original_flags = os.flags();

    os.setf(std::ios::left);
    os.width(5);
    os << "N";

    os.width(13);
    os << "Status";

    os.setf(std::ios::left);
    os.width(14);
    os << "||F||";

    const auto* is_linesearch = dynamic_cast<const NOX::Solver::LineSearchBased*>(&solver);
    if (is_linesearch) {
      os.setf(std::ios::left);
      os.width(14);
      os << "Step Size";
    }

    os.setf(std::ios::left);
    os.width(11);
    os << "Linear Its";

    os.setf(std::ios::left);
    os.width(14);
    os << "Achieved Tol";

    os << std::endl;

    os.flags(original_flags);
    this->printStep(solver);
  }
}

void NOX::ObserverPrint::runPostIterate(const NOX::Solver::Generic& solver)
{
  this->printStep(solver);
}

void NOX::ObserverPrint::printStep(const NOX::Solver::Generic& solver)
{
  const auto& stats = *solver.getSolverStatistics();
  auto& os = os_->out();
  auto original_flags = os.flags();
  const int precision = 6;

  os.width(5);
  os.setf(std::ios::left);
  os << stats.numNonlinearIterations;

  os.width(13);
  if (solver.getStatus() == NOX::StatusTest::Unconverged)
    os << "Unconverged";
  else if (solver.getStatus() == NOX::StatusTest::Converged)
    os << "Converged!";
  else if (solver.getStatus() == NOX::StatusTest::Failed)
    os << "Failed!";

  os.width(14);
  os.precision(precision);
  os.setf(std::ios::left|std::ios::scientific);
  auto& grp = solver.getSolutionGroup();
  if (!grp.isF())
    const_cast<NOX::Abstract::Group&>(grp).computeF();
  os << grp.getNormF();

  const auto* is_linesearch = dynamic_cast<const NOX::Solver::LineSearchBased*>(&solver);
  if (is_linesearch) {
    os.width(14);
    os.setf(std::ios::left|std::ios::scientific);
    os.precision(precision);
    os << is_linesearch->getStepSize();
  }

  os.width(11);
  os.setf(std::ios::left);
  os.precision(precision);
  os << stats.linearSolve.lastLinearSolve_NumIterations;

  os.width(14);
  os.setf(std::ios::left|std::ios::scientific);
  os.precision(precision);
  os << stats.linearSolve.lastLinearSolve_AchievedTolerance;

  os << std::endl;

  os.flags(original_flags);
}
