// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_NONLINEARSOLVER_NOX
#define THYRA_NONLINEARSOLVER_NOX

#include "Teuchos_RCP.hpp"
#include "Thyra_NonlinearSolverBase.hpp"
#include "NOX_PrePostOperator_RowSumScaling.H"

namespace NOX {
  class SolverManager;
  namespace Abstract {
    class Group;
  }
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
    const SolveCriteria<double> *solveCriteria = nullptr,
    VectorBase<double> *delta = nullptr
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

  /** \brief . */
  RCP< const LinearOpBase<double> >
  get_W_op() const;
  /** \brief . */
  RCP< LinearOpBase<double> >
  get_nonconst_W_op(const bool forceUpToDate);

  /** \brief . */
  RCP< const PreconditionerBase<double> >
  get_prec_op() const;
  /** \brief . */
  RCP< PreconditionerBase<double> >
  get_nonconst_prec_op();

  //@}

  /** \brief . */
  void setBasePoint(const ModelEvaluatorBase::InArgs<double> &modelInArgs);

  /** \brief Users can optionally set the preconditioner

      \param[in] precOp The preconditioner operator
      \param[in] precFactory Optional preconditioner factory
      \param[in] updatePreconditioner Optional flag, if true the Group will automatically update the preconditioner using precFactory or model evalautor
   */
  void setPrecOp(const Teuchos::RCP< ::Thyra::PreconditionerBase<double>>& precOp,
                 const Teuchos::RCP< ::Thyra::PreconditionerFactoryBase<double>>& precFactory = Teuchos::null,
                 const bool updatePreconditioner = true);

  /** \brief Power user function to set the NOX solution Group. Normally, this will be created automatically but there are exceptional use cases where we need to set manually.

      \param[in] group The nox group with the initial guess
   */
  void setGroup(const Teuchos::RCP<NOX::Abstract::Group>& group);

  RCP<const NOX::Solver::Generic> getNOXSolver() const;

  //! Called to rebuild the solver if a new parameter list is set.
  void resetSolver();

private:

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

  // Returns the underlying NOX::Thyra::Group. First tries to dynamic
  // cast the group input param to a thyra group. If that fails, it
  // will try to pull out a nested group and recursively call this
  // function to eventually get a thyra group. This is needed when we
  // use groups that use a nesting approach as opposed to inheritance
  // (LOCA does this alot).
  Teuchos::RCP<NOX::Thyra::Group> getThyraGroupNonConst(const Teuchos::RCP<NOX::Abstract::Group>& group);

  // Returns the underlying NOX::Thyra::Group. First tries to dynamic
  // cast the group input param to a thyra group. If that fails, it
  // will try to pull out a nested group and recursively call this
  // function to eventually get a thyra group. This is needed when we
  // use groups that use a nesting approach as opposed to inheritance
  // (LOCA does this alot).
  Teuchos::RCP<const NOX::Thyra::Group> getThyraGroupConst(const Teuchos::RCP<const NOX::Abstract::Group>& group) const;

private:

  RCP<Teuchos::ParameterList> param_list_;
  RCP<Teuchos::ParameterList> valid_param_list_;
  RCP<const ModelEvaluator<double> > model_;
  ModelEvaluatorBase::InArgs<double> basePoint_;

  //! Note that this is not a thyra group. To support loca groups that
  //! use nesting, the nox_group is the base class. When we need the
  //! underlying Thyra group, the getThyraGroup calls can pull out the
  //! nested group.
  RCP<NOX::Abstract::Group> nox_group_;
  RCP<NOX::StatusTest::Generic> status_test_;
  RCP<NOX::Solver::Generic> solver_;

  // Options
  std::string function_scaling_;
  bool do_row_sum_scaling_;
  NOX::RowSumScaling::ENOX_WhenToUpdateScaling when_to_update_;
  Teuchos::RCP< ::Thyra::VectorBase<double> > scaling_vector_;
  Teuchos::RCP< ::Thyra::VectorBase<double> > right_scaling_vector_;
  Teuchos::RCP< ::Thyra::VectorBase<double> > inv_right_scaling_vector_;
  bool rightScalingFirst_;

  bool rebuild_solver_;

  Teuchos::RCP< ::Thyra::PreconditionerBase<double>> precOp_;
  Teuchos::RCP< ::Thyra::PreconditionerFactoryBase<double>> precFactory_;
  bool updatePreconditioner_;

  RCP<NOX::Abstract::Group> user_defined_nox_group_;

  bool use_base_point_;
};



} // namespace Thyra


#endif // THYRA_NONLINEARSOLVER_NOX
