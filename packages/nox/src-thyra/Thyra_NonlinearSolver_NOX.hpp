// @HEADER
// @HEADER

#ifndef THYRA_NONLINEARSOLVER_NOX
#define THYRA_NONLINEARSOLVER_NOX


#include "Thyra_NonlinearSolverBase.hpp"

namespace NOX {
  class SolverManager;
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

  //! Called to rebuild the solver if a new parameter list 
  void resetSolver();

private:

  Teuchos::RCP<Teuchos::ParameterList> param_list_;
  Teuchos::RCP<Teuchos::ParameterList> valid_param_list_;
  Teuchos::RCP<const ModelEvaluator<double> > model_;
  Teuchos::RCP<LinearOpWithSolveBase<double> > J_;

  Teuchos::RCP<NOX::SolverManager> solver_;

};



} // namespace Thyra


#endif // THYRA_NONLINEARSOLVER_NOX
