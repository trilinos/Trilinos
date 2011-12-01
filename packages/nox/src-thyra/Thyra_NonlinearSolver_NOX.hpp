// @HEADER
// @HEADER

#ifndef THYRA_NONLINEARSOLVER_NOX
#define THYRA_NONLINEARSOLVER_NOX

#include "Teuchos_RCP.hpp"
#include "Thyra_NonlinearSolverBase.hpp"
#include "NOX_PrePostOperator_RowSumScaling.H"

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

  /** \brief . */
  void setParameterList(RCP<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  RCP<Teuchos::ParameterList> getNonconstParameterList();
  /** \brief . */
  RCP<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  RCP<const Teuchos::ParameterList> getParameterList() const;
  /** \brief . */
  RCP<const Teuchos::ParameterList> getValidParameters() const;

  //@}

  /** @name Overridden from NonlinearSolverBase */
  //@{

  /** \brief . */
  void setModel(const RCP< const ModelEvaluator<double> > &model);
  /** \brief . */
  RCP< const ModelEvaluator<double> > getModel() const;
  /** \brief . */
  SolveStatus<double> solve(
    VectorBase<double> *x,
    const SolveCriteria<double> *solveCriteria,
    VectorBase<double> *delta
    );
  /** \brief . */
  RCP<const VectorBase<double> > get_current_x() const;
  /** \brief . */
  bool is_W_current() const;
  /** \brief . */
  RCP< LinearOpWithSolveBase<double> >
  get_nonconst_W(const bool forceUpToDate);
  /** \brief . */
  RCP< const LinearOpWithSolveBase<double> > get_W() const;
  
  //@}

  RCP<const NOX::Solver::Generic> getNOXSolver() const;

private:

  //! Called to rebuild the solver if a new parameter list is set.
  void resetSolver();

  //! Builds status tests - first looks for parameter list to use, otherwise builds a default set of status tests.
  RCP<NOX::StatusTest::Generic> 
  buildStatusTests(Teuchos::ParameterList& p);

  void
  validateAndParseThyraGroupOptions(Teuchos::ParameterList& thyra_group_sublist);

  /** If row sum scaling is enabled, allocates scaling vector and
      creates and sets the nox rePostOperator in the main parmaeter
      list.
  */
  void setupRowSumScalingObjects();

private:

  RCP<Teuchos::ParameterList> param_list_;
  RCP<Teuchos::ParameterList> valid_param_list_;
  RCP<const ModelEvaluator<double> > model_;

  RCP<NOX::Thyra::Group> nox_group_;
  RCP<NOX::StatusTest::Generic> status_test_;
  RCP<NOX::Solver::Generic> solver_;

  // Options
  std::string function_scaling_;
  bool do_row_sum_scaling_;
  NOX::RowSumScaling::ENOX_WhenToUpdateScaling when_to_update_;
  Teuchos::RCP< ::Thyra::VectorBase<double> > scaling_vector_;

};



} // namespace Thyra


#endif // THYRA_NONLINEARSOLVER_NOX
