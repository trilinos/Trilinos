// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#include <MatrixBased_LOWS.hpp>
#include "Piro_config.hpp"
#include <Teuchos_TestForException.hpp>
#ifdef HAVE_PIRO_STRATIMIKOS
  #include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#endif

MatrixBased_LOWS::
    MatrixBased_LOWS(
        const Teuchos::RCP<Thyra::LinearOpBase<double>> &mat) : mat_(mat) {}

MatrixBased_LOWS::
    ~MatrixBased_LOWS() {}

Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
MatrixBased_LOWS::
    domain() const
{
  return mat_->domain();
}

Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
MatrixBased_LOWS::
    range() const
{
  return mat_->range();
}

Teuchos::RCP<Thyra::LinearOpBase<double>>
MatrixBased_LOWS::
    getMatrix()
{
  return mat_;
}

void
MatrixBased_LOWS::
    initializeSolver(Teuchos::RCP<Teuchos::ParameterList> solverParamList)
{
  #ifdef HAVE_PIRO_STRATIMIKOS
    std::string solverType = solverParamList->get<std::string>("Linear Solver Type");
    Stratimikos::DefaultLinearSolverBuilder strat;
    strat.setParameterList(solverParamList);
    auto lows_factory = strat.createLinearSolveStrategy(solverType);
    solver_ = lows_factory->createOp();
    Thyra::initializeOp<double>(*lows_factory, mat_, solver_.ptr(), Thyra::SUPPORT_SOLVE_FORWARD_ONLY);
  #else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error! MatrixBased_LOWS::initializeSolver, Implementation requires Stratimikos.\n");
  #endif
}

bool
MatrixBased_LOWS::
    opSupportedImpl(Thyra::EOpTransp M_trans) const
{
  return mat_->opSupported(M_trans);
}

void
MatrixBased_LOWS::
    applyImpl(const Thyra::EOpTransp M_trans,
              const Thyra::MultiVectorBase<double> &X,
              const Teuchos::Ptr<Thyra::MultiVectorBase<double>> &Y,
              const double alpha,
              const double beta) const
{
  mat_->apply(M_trans, X, Y, alpha, beta);
}

Thyra::SolveStatus<double>
MatrixBased_LOWS::
    solveImpl(
        const Thyra::EOpTransp transp,
        const Thyra::MultiVectorBase<double> &B,
        const Teuchos::Ptr<Thyra::MultiVectorBase<double>> &X,
        const Teuchos::Ptr<const Thyra::SolveCriteria<double>> solveCriteria) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(solver_), std::runtime_error, "Error! MatrixBased_LOWS::solveImpl, Solver not initialized, call initializeSolver first.\n");
  return solver_->solve(transp, B, X, solveCriteria);
}