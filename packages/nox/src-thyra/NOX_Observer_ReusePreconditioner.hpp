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

    /// Setup the preconditioner objects used to update the preconditioner
    void setOperatorsAndFactory(const Teuchos::RCP<::Thyra::PreconditionerBase<Scalar>>& precOperator,
                                const Teuchos::RCP<::Thyra::PreconditionerFactoryBase<Scalar>>& precFactory);

    /** \brief Enables updating the preconditioner at the start of each nonlinear solve */
    void updateAtStartOfSolve();

    /** \brief Enables updating of the preconditioner after a set number of nonlinear iteraitons.

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

    // Derived from NOX::Observer
    void runPreSolve(const NOX::Solver::Generic& solver) override;
    void runPreIterate(const NOX::Solver::Generic& solver) override;
    void runPostIterate(const NOX::Solver::Generic& solver) override;

  private:
    void updatePreconditioner(NOX::Thyra::Group& group);
    bool isInitialized() const;

  private:
    bool update_at_start_of_solve_;

    bool update_after_n_iterations_;
    int num_iterations_for_update_;

    bool update_if_stalled_;
    int max_count_for_stall_;
    int max_linear_iterations_for_stall_;
    int stall_count_;

    int iterations_since_last_update_;

    Teuchos::RCP<::Thyra::PreconditionerBase<Scalar>> precOperator_;
    Teuchos::RCP<::Thyra::PreconditionerFactoryBase<Scalar>> precFactory_;
  };

}

#endif
