// @HEADER
// @HEADER

#ifndef THYRA_NONLINEARSOLVER_NOX
#define THYRA_NONLINEARSOLVER_NOX

#include "Teuchos_RCP.hpp"
#include "Thyra_NonlinearSolverBase.hpp"

namespace NOX {
  class SolverManager;
  namespace Thyra {
    class Group;
  }
  namespace StatusTest {
    class Generic;
  }
  namespace Solver {
    class Generic;
  }
}

namespace Thyra {


/** \brief Concrete nonlinear solver for NOX.
 *
 * This class implemets a NOX nonlinear solver of type <double>.
 */
class NOXNonlinearSolver : public NonlinearSolverBase<double> {
public:

  NOXNonlinearSolver();
  ~NOXNonlinearSolver();

  /** @name Overridden from ParameterListAcceptor */
  //@{

  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);
  Teuchos::RCP<Teuchos::ParameterList> getParameterList();
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
  Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  //@}

  /** @name Overridden from NonlinearSolverBase */
  //@{

  void setModel(const Teuchos::RCP< const ModelEvaluator<double> > &model);
  Teuchos::RCP< const ModelEvaluator<double> > getModel() const;
  SolveStatus<double> solve(VectorBase<double> *x,
			    const SolveCriteria<double> *solveCriteria,
			    VectorBase<double> *delta);
  Teuchos::RCP< LinearOpWithSolveBase<double> > get_nonconst_W();
  Teuchos::RCP< const LinearOpWithSolveBase<double> > get_W() const;
  
  //@}

private:

  //! Called to rebuild the solver if a new parameter list is set.
  void resetSolver();

  //! Builds status tests - first looks for parameter list to use, otherwise builds a default set of status tests.
  Teuchos::RCP<NOX::StatusTest::Generic> 
  buildStatusTests(Teuchos::ParameterList& p);

private:

  Teuchos::RCP<Teuchos::ParameterList> param_list_;
  Teuchos::RCP<Teuchos::ParameterList> valid_param_list_;
  Teuchos::RCP<const ModelEvaluator<double> > model_;
  Teuchos::RCP<LinearOpWithSolveBase<double> > J_;

  Teuchos::RCP<NOX::Thyra::Group> nox_group_;
  Teuchos::RCP<NOX::StatusTest::Generic> status_test_;
  Teuchos::RCP<NOX::Solver::Generic> solver_;

};



} // namespace Thyra


#endif // THYRA_NONLINEARSOLVER_NOX
