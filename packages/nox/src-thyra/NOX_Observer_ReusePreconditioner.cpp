// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#include "NOX_Observer_ReusePreconditioner.hpp"
#include "NOX_Solver_Generic.H"
#include "NOX_SolverStats.hpp"
#include "NOX_Thyra_Group.H"
#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"

NOX::ObserverReusePreconditioner::ObserverReusePreconditioner()
  : update_at_start_of_nonlinear_solve_(false),
    update_after_n_nonlinear_solves_(false),
    num_nonlinear_solves_for_update_(-1),
    reset_nonlinear_solve_count_on_failed_solve_(false),
    num_nonlinear_solves_count_(0),
    update_after_n_linear_iterations_(false),
    num_linear_iterations_for_update_(5),
    update_if_stalled_(false),
    max_count_for_stall_(-1),
    max_linear_iterations_for_stall_(-1),
    stall_count_(0),
    iterations_since_last_update_(0),
    update_preconditioner_count_(0)
{}

void
NOX::ObserverReusePreconditioner::
setOperatorsAndFactory(const Teuchos::RCP<::Thyra::PreconditionerBase<Scalar>>& precOperator,
                       const Teuchos::RCP<::Thyra::PreconditionerFactoryBase<Scalar>>& precFactory)
{
  precOperator_ = precOperator;
  precFactory_ = precFactory;
}

void
NOX::ObserverReusePreconditioner::updateAtStartOfSolve()
{
  update_at_start_of_nonlinear_solve_ = true;
}

void
NOX::ObserverReusePreconditioner::
updateAfterNNonlinearSolves(const int num_nonlinear_solves_for_update,
                            const bool reset_nonlinear_solve_count_on_failed_solve)
{
  update_after_n_nonlinear_solves_ = true;
  num_nonlinear_solves_for_update_ = num_nonlinear_solves_for_update;
  reset_nonlinear_solve_count_on_failed_solve_ = reset_nonlinear_solve_count_on_failed_solve;
}

void
NOX::ObserverReusePreconditioner::
updateAfterNIterations(const int num_iterations_for_update)
{
  update_after_n_linear_iterations_ = true;
  num_linear_iterations_for_update_ = num_iterations_for_update;
}

void
NOX::ObserverReusePreconditioner::
updateOnLinearSolverStall(const int max_linear_iterations,
                          const int max_count)
{
  update_if_stalled_ = true;
  max_linear_iterations_for_stall_ = max_linear_iterations;
  max_count_for_stall_ = max_count;
}

void NOX::ObserverReusePreconditioner::runPreSolve(const NOX::Solver::Generic& solver)
{
  // First time through, grab the operator and tell the group that the
  // user will control preconditioning.
  if (precOperator_.is_null()) {
    const auto& const_group = solver.getSolutionGroup();
    auto& group = const_cast<NOX::Abstract::Group&>(const_group);
    auto& thyra_group = dynamic_cast<NOX::Thyra::Group&>(group);

    auto lows_factory = thyra_group.getLinearOpWithSolveFactory();
    TEUCHOS_ASSERT(nonnull(lows_factory));

    // Creating the preconditioner objects follows the logic of the
    // NOX::Thyra::Group constructors for consistency.
    precFactory_ = thyra_group.getPreconditionerFactory();
    if (precFactory_.is_null()) {
      precFactory_ = lows_factory->getPreconditionerFactory();
    }

    precOperator_ = thyra_group.getNonconstPreconditioner();
    if (precOperator_.is_null()) {
      if (nonnull(precFactory_)) {
        precOperator_ = precFactory_->createPrec();
      } else if (thyra_group.getModel()->createOutArgs().supports( ::Thyra::ModelEvaluatorBase::OUT_ARG_W_prec)) {
        // precOperator_ = thyra_group.getModel()->create_W_prec();
      }
    }

    thyra_group.takeControlOfPreconditionerUpdates(precOperator_);

    TEUCHOS_TEST_FOR_EXCEPTION(precFactory_.is_null(),std::runtime_error,"ERROR: the ReusePreconditioner observer is registered but there is no preconditioner available!");
    TEUCHOS_TEST_FOR_EXCEPTION(precOperator_.is_null(),std::runtime_error,"ERROR: the ReusePreconditioner observer is registered but there is no preconditioner available!");

    // First time through, create the preconditioner
    // this->updatePreconditioner(thyra_group);
  }
}

void NOX::ObserverReusePreconditioner::runPreIterate(const NOX::Solver::Generic& solver)
{
  TEUCHOS_ASSERT(this->isInitialized());

  bool update_the_preconditioner = false;

  if (update_at_start_of_nonlinear_solve_ && (solver.getNumIterations() == 0) ) {
    update_the_preconditioner = true;
  }

  if (update_after_n_nonlinear_solves_ && (solver.getNumIterations() == 0) ) {
    if (num_nonlinear_solves_count_ % num_nonlinear_solves_for_update_ == 0) {
      update_the_preconditioner = true;
    }
  }

  if (update_after_n_linear_iterations_) {
    if (solver.getNumIterations() % num_linear_iterations_for_update_ == 0) {
      update_the_preconditioner = true;
    }
  }

  if (update_if_stalled_) {
    const int num_linear_iterations = solver.getSolverStatistics()->linearSolve.lastLinearSolve_NumIterations;
    const bool is_converged = solver.getSolverStatistics()->linearSolve.lastLinearSolve_Converged;

    if (num_linear_iterations >= max_linear_iterations_for_stall_) {
      ++stall_count_;
    } else {  // reset count on a successful linear solve
      stall_count_ = 0;
    }

    if ( (!is_converged) || (stall_count_ >= max_count_for_stall_) )
      update_the_preconditioner = true;
  }

  if (update_the_preconditioner) {
    const auto& const_group = solver.getSolutionGroup();
    auto& group = const_cast<NOX::Abstract::Group&>(const_group);
    auto& thyraGroup = dynamic_cast<NOX::Thyra::Group&>(group);
    this->updatePreconditioner(thyraGroup);
  }
}

void NOX::ObserverReusePreconditioner::runPostIterate(const NOX::Solver::Generic& solver)
{++iterations_since_last_update_;}

void NOX::ObserverReusePreconditioner::runPostSolve(const NOX::Solver::Generic& solver)
{
  ++num_nonlinear_solves_count_;

  if (reset_nonlinear_solve_count_on_failed_solve_) {
    if (solver.getStatus() == NOX::StatusTest::Failed) {
      num_nonlinear_solves_count_ = 0;
    }
  }
}

size_t NOX::ObserverReusePreconditioner::getNumPreconditionerUpdates() const
{return update_preconditioner_count_;}

size_t NOX::ObserverReusePreconditioner::getNumNonlinearSolvesCount() const
{return num_nonlinear_solves_count_;}

void NOX::ObserverReusePreconditioner::updatePreconditioner(NOX::Thyra::Group& group)
{
  TEUCHOS_ASSERT(this->isInitialized());

  if (!group.isJacobian())
    group.computeJacobian();

  auto jacOperator = group.getScaledJacobianOperator();

  ::Thyra::initializePrec<double>(*precFactory_,jacOperator,precOperator_.ptr());

  group.unscaleJacobianOperator();

  stall_count_ = 0;
  iterations_since_last_update_ = 0;

  ++update_preconditioner_count_;
}

bool NOX::ObserverReusePreconditioner::isInitialized() const
{
  return ( (update_at_start_of_nonlinear_solve_ ||
            update_after_n_nonlinear_solves_ ||
            update_after_n_linear_iterations_ ||
            update_if_stalled_) &&
           nonnull(precOperator_) && nonnull(precFactory_) );
}
