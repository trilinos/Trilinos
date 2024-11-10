// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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