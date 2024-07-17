// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef NOX_OBSERVER_REUSE_PRECONDITIONER_HPP
#define NOX_OBSERVER_REUSE_PRECONDITIONER_HPP

#include "NOX_Config.h"
#include "NOX_Observer.hpp"
#include "Teuchos_RCP.hpp"

namespace Thyra {
  template<typename T> class PreconditionerBase;
  template<typename T> class PreconditionerFactoryBase;
}

namespace NOX {

  namespace Thyra {
    class Group;
  }

  /** \brief Observer that controls when to update the preconditioner for the Thyra interface.

      The preconditioner can be updated at the start of each new
      nonlinear solve and/or when convergence stalls out.
   */
  class ObserverReusePreconditioner : public NOX::Observer {
  public:
    using Scalar = double;

    ObserverReusePreconditioner();

    /** @name Initialization methods */
    //@{

    /// Setup the preconditioner objects used to update the preconditioner
    void setOperatorsAndFactory(const Teuchos::RCP<::Thyra::PreconditionerBase<Scalar>>& precOperator,
                                const Teuchos::RCP<::Thyra::PreconditionerFactoryBase<Scalar>>& precFactory);

    /** \brief Enables updating the preconditioner at the start of each nonlinear solve */
    void updateAtStartOfSolve();

    /** \brief Enables updating of the preconditioner after a set number of nonlinear solves.

        This is intended to reuse a preconditioner across all stages
        of an RK method in a single time step. The parameter should
        be set to the number of RK stages. If a nonlinear solve fails

        @param [in] num_nonlinear_solves_for_update (int) Updates the preconditioner after this number of nonlinear iterations.
        @param [in] reset_nonlinear_solve_count_on_faile_solve (bool) If set to true, when a nonlinear solve fails, the nonlinear solve count will be reset. When a nonlinear solve fails in the middle of an RK stage, we assume a new time step will start for the next nonlinear solve.
     */
    void updateAfterNNonlinearSolves(const int num_nonlinear_solves_for_update,
                                     const bool reset_nonlinear_solve_count_on_failed_solve = true);

    /** \brief Enables updating of the preconditioner after a set number of nonlinear iterations.

        @param [in] num_iterations_for_update (int) Updates the preconditioner after this number of nonlinear iterations.
     */
    void updateAfterNIterations(const int num_iterations_for_update);

    /** \brief Enables updating of preconditioner if the observer detects a stall or failure in the linear solver.

        This algorithm tries to assess a stalled computation due to
        reusing the preconditioner. It will always recompute for a
        failed linear solve. It will also recompute if the last
        max_count number of iterations each had linear solves that
        took more iterations than max_linear_iterations.

        @param [in] max_linear_iterations (int) Declare a stalled iteraiton if the number of linear solver iterations is above this value
        @param [in] max_count (int) Recompute the preconditioner after this many stalled iterations
     */
    void updateOnLinearSolverStall(const int max_linear_iterations,
                                   const int max_count);

    //@}

    /** @name Methods derived from NOX::Observer */
    //@{
    void runPreSolve(const NOX::Solver::Generic& solver) override;
    void runPreIterate(const NOX::Solver::Generic& solver) override;
    void runPostIterate(const NOX::Solver::Generic& solver) override;
    void runPostSolve(const NOX::Solver::Generic& solver) override;
    //@}

    /** @name Query methods used in unit testing */
    //@{
    /// Return the number of times the preconditioner has been updated.
    size_t getNumPreconditionerUpdates() const;
    /// Return the number of nonlinear solves that have been run.
    size_t getNumNonlinearSolvesCount() const;
    //@}

  private:
    void updatePreconditioner(NOX::Thyra::Group& group);
    bool isInitialized() const;

  private:
    bool update_at_start_of_nonlinear_solve_;

    bool update_after_n_nonlinear_solves_;
    int num_nonlinear_solves_for_update_;
    bool reset_nonlinear_solve_count_on_failed_solve_;
    size_t num_nonlinear_solves_count_;

    bool update_after_n_linear_iterations_;
    int num_linear_iterations_for_update_;

    bool update_if_stalled_;
    int max_count_for_stall_;
    int max_linear_iterations_for_stall_;
    int stall_count_;

    int iterations_since_last_update_;

    size_t update_preconditioner_count_;

    Teuchos::RCP<::Thyra::PreconditionerBase<Scalar>> precOperator_;
    Teuchos::RCP<::Thyra::PreconditionerFactoryBase<Scalar>> precFactory_;
  };

}

#endif
