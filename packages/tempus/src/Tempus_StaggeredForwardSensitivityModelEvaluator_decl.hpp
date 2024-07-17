//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StaggeredForwardSensitivityModelEvaluator_decl_hpp
#define Tempus_StaggeredForwardSensitivityModelEvaluator_decl_hpp

#include "Thyra_StateFuncModelEvaluatorBase.hpp"
#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"
#include "NOX_Thyra.H"

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_SensitivityModelEvaluatorBase.hpp"

namespace Tempus {

/** \brief Transform a ModelEvaluator's sensitivity equations to its residual
 *
 * This class wraps a given ModelEvalutor encapsulating f(x,p) and creates
 * a new "residual" for the forward sensitivity equations:
 *        F(X) = (df/dx)(x,p) * X + df/dp(x,p) = 0
 * where X = dx/dp (transient terms supressed for simplicity).  This model
 * evaluator can then be handed to a regular (non)linear solver to compute X.
 * Note that even though these equations are linear in X, it is not necessarily
 * the case that the underlying model evaluator accurately computes df/dx in its
 * evaluation of W.  Therefore this model evaluator can optionally reinterpret
 * the model's df/dp out-arg as (df/dx)(x,p) * dx/dp + df/dp(x,p) where dx/dp is
 * passed as another parameter (product) vector (encapsulated in the
 * Thyra::DefaultMultiVectorProductVector).  This is not standard model
 * evaluator behavior, but is useful for models where W is only an approximation
 * to df/dx and/or the model is capable of directly computing
 * (df/dx)(x,p) * dx/dp + df/dp(x,p).
 *
 * This model evaluator differes from CombinedForwardSensitivityModelEvaluator
 * in that it doesn't include the state equations in the residual, just the
 * sensitivity equations.  Therefore it provides methods to set the state
 * solution vector (x) and time derivatives (x_dot, x_dot_dot) for use in
 * evaluating the senstivity residual.  It also provides methods for setting
 * the linear operator (W a.k.a. alpha*df/dx + beta*df/dx_dot) and its
 * preconditioner in cases where they can be reused from the state model
 * evaluations.
 */
template <typename Scalar>
class StaggeredForwardSensitivityModelEvaluator
  : public Thyra::StateFuncModelEvaluatorBase<Scalar>,
    public SensitivityModelEvaluatorBase<Scalar> {
 public:
  typedef Thyra::VectorBase<Scalar> Vector;
  typedef Thyra::MultiVectorBase<Scalar> MultiVector;

  //! Constructor
  /*!
   * The optionally supplied parameter list supports the following options:
   * <ul>
   *   <li> "Use DfDp as Tangent" (default:  false) Reinterpret the df/dp
   *        out-arg as the tangent vector (df/dx)(x,p) * dx/dp + df/dp(x,p)
   *        as described above.  If it is false, this implementation will
   *        compute the tangent through several calls to the models
   *        evalModel().
   *   <li> "Sensitivity Parameter Index" (default: 0) Model evaluator
   *        parameter index for which sensitivities will be computed.
   *   <li> "Response Function Index", (default: 0) The model evaluator
   *        response index for which sensitivities will be computed if
   *        requested.
   *   <li> "Sensitivity X Tangent Index" (default: 1) If "Use DfDp as Tangent"
   *        is true, the model evaluator parameter index for passing dx/dp
   *        as a Thyra::DefaultMultiVectorProductVector.
   *   <li> "Sensitivity X-Dot Tangent Index" (default: 2) If
   *        "Use DfDp as Tangent" is true, the model evaluator parameter index
   *        for passing dx_dot/dp as a Thyra::DefaultMultiVectorProductVector.
   *   <li> "Sensitivity X-Dot-Dot Tangent Index" (default: 3) If
   *        "Use DfDp as Tangent" is true, the model evaluator parameter index
   *        for passing dx_dot_dot/dp as a
   *        Thyra::DefaultMultiVectorProductVector (if the model supports
   *        x_dot_dot).
   * </ul>
   */
  StaggeredForwardSensitivityModelEvaluator(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >&
          sens_residual_model,
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >&
          sens_solve_model,
      const bool is_pseudotransient,
      const Teuchos::RCP<const Teuchos::ParameterList>& pList = Teuchos::null,
      const Teuchos::RCP<MultiVector>& dxdp_init              = Teuchos::null,
      const Teuchos::RCP<MultiVector>& dx_dotdp_init          = Teuchos::null,
      const Teuchos::RCP<MultiVector>& dx_dotdot_dp_init      = Teuchos::null);

  /** \name Public functions overridden from SensitivityModelEvaulator. */
  //@{

  //! Get the underlying model 'f'
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > getForwardModel() const
  {
    return model_;
  }

  //! Set solution history from forward state evaluation (for interpolation)
  void setForwardSolutionHistory(
      const Teuchos::RCP<const Tempus::SolutionHistory<Scalar> >& sh);

  //! Set solution state from forward state evaluation (for frozen state)
  virtual void setForwardSolutionState(
      const Teuchos::RCP<const Tempus::SolutionState<Scalar> >& s);

  //! Set the solver of the underlying model if you want to reuse it
  virtual void setSolver(
      const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
      const bool force_W_update)
  {
    auto nox_solver =
        Teuchos::rcp_dynamic_cast<Thyra::NOXNonlinearSolver>(solver, true);
    lo_ = nox_solver->get_nonconst_W_op(force_W_update);
    po_ = nox_solver->get_nonconst_prec_op();
  }

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_space(int p) const;

  Teuchos::RCP<const Teuchos::Array<std::string> > get_p_names(int p) const;

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_space() const;

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_g_space(int j) const;

  Teuchos::ArrayView<const std::string> get_g_names(int j) const;

  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_W_op() const;

  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_DgDx_dot_op(int j) const;

  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_DgDx_op(int j) const;

  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_DgDp_op(int j, int l) const;

  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
  get_W_factory() const;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;

  //@}

  static Teuchos::RCP<const Teuchos::ParameterList> getValidParameters();

 private:
  typedef Thyra::DefaultMultiVectorProductVectorSpace<Scalar> DMVPVS;

  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;

  void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> prototypeInArgs_;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> prototypeOutArgs_;

  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > model_;
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > sens_residual_model_;
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > sens_solve_model_;
  Teuchos::RCP<MultiVector> dxdp_init_;
  Teuchos::RCP<MultiVector> dx_dotdp_init_;
  Teuchos::RCP<MultiVector> dx_dotdotdp_init_;
  int p_index_;
  int g_index_;
  int x_tangent_index_;
  int xdot_tangent_index_;
  int xdotdot_tangent_index_;
  bool use_dfdp_as_tangent_;
  bool use_dgdp_as_tangent_;

  int num_param_;
  int num_response_;
  int g_offset_;
  Teuchos::RCP<const DMVPVS> dxdp_space_;
  Teuchos::RCP<const DMVPVS> dfdp_space_;
  Teuchos::RCP<const DMVPVS> dgdp_space_;
  Teuchos::RCP<const Tempus::SolutionHistory<Scalar> > sh_;

  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > lo_;
  Teuchos::RCP<Thyra::PreconditionerBase<Scalar> > po_;

  bool is_pseudotransient_;
  mutable bool mass_matrix_is_computed_;
  mutable bool jacobian_matrix_is_computed_;
  mutable bool acceleration_matrix_is_computed_;
  mutable bool residual_sensitivity_is_computed_;
  mutable Teuchos::RCP<Thyra::LinearOpBase<Scalar> > my_dfdx_;
  mutable Teuchos::RCP<Thyra::LinearOpBase<Scalar> > my_dfdxdot_;
  mutable Teuchos::RCP<Thyra::LinearOpBase<Scalar> > my_dfdxdotdot_;
  mutable Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > my_dfdp_;
  mutable Teuchos::RCP<Thyra::LinearOpBase<Scalar> > my_dgdx_;
  mutable Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > my_dgdx_mv_;
  mutable Teuchos::RCP<const Tempus::SolutionState<Scalar> > forward_state_;
  mutable Teuchos::RCP<Tempus::SolutionState<Scalar> > nc_forward_state_;
  mutable Scalar t_interp_;
};

}  // namespace Tempus

#endif
