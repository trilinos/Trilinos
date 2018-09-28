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
//@HEADER

#include "NOX_Observer_Print.hpp"
#include "NOX_Solver_Generic.H"
#include "NOX_Abstract_Group.H"
#include "NOX_SolverStats.hpp"
#include "NOX_Utils.H"

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

    os.width(12);
    os << "Status";

    os.setf(std::ios::left);
    os.width(14);
    os << "||F||";

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

  os.width(12);
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
  if (not grp.isF())
    const_cast<NOX::Abstract::Group&>(grp).computeF();
  os << grp.getNormF();

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
