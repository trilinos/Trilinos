//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_AdjointSensitivityModelEvaluator_decl_hpp
#define Tempus_AdjointSensitivityModelEvaluator_decl_hpp

#include "Thyra_StateFuncModelEvaluatorBase.hpp"
#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"

namespace Tempus {

/** \brief ModelEvaluator for forming adjoint sensitivity equations. */
/**
 * This class wraps a given ModelEvalutor encapsulating f(x_dot,x,p) and creates
 * a new "residual" for the adjoint sensitivity equations:
 *        F(y) = d/dt( df/dx_dot^T*y ) - df/dx^T*y + dg/dx^T = 0
 * where y is the adjoint variable.  Nominally y is a multi-vector belonging
 * to the vector space of the residual f with number of columns equal to the
 * dimension of g.  To satisfy the model evaluator interface, y is converted
 * to a product vector formed by its columns.  This is the form of the adjoint
 * equations suitable for computing pseudotransientsensitivities of a response
 * function g(x(T),p) or adjoint sensitivities for a time-integrated response
 * function int_0^T g(x(t),p).
 *
 * To compute adjoint sensitivities, the equations f(x_dot,x,p) = 0 are
 * integrated forward in time to some final time T, with the adjoint equations
 * above integrated backward in time starting at T.  Since the tempus
 * integrator can only integrate forward in time, we define tau = T - t and
 * transform the adjoint equations to:
 *        F(y) = d/dtau( df/dx_dot^T*y ) + df/dx^T*y - dg/dx^T = 0.
 * The initial conditions for y are y(T) = 0. The sensitivity of g(T) is then
 *        int_0^T(dg/dp(T) - df/dp^T*y)dt + dx/dp(0)^T*df/dx_dot(0)^T*y(0)
 * for transient adjoint sensitivites and is
 *        dg/dp(T) - df/dp^T*y
 * for pseudotransient.
 *
 * This model evaluator supports both implicit and explict forms of the
 * governing equations.  However it assumes df/dx_dot is constant with respect
 * to x_dot, x, and t so that the adjoint equations become
 *        F(y) = df/dxdot^T*y_dot + df/dx^T*y - dg/dx^T = 0.
 * Moving beyond this assumption will require modifications to the steppers
 * in how they generate time-derivative terms.
 */
template <typename Scalar>
class AdjointSensitivityModelEvaluator
  : public Thyra::StateFuncModelEvaluatorBase<Scalar> {
 public:
  typedef Thyra::VectorBase<Scalar> Vector;
  typedef Thyra::MultiVectorBase<Scalar> MultiVector;

  //! Constructor
  /*!
   * t_final is the final integration time used to implement the time
   * transformation described above.  num_adjoint is the number of adjoint
   * variables defining the number of columns in y, which is determined by
   * the number of elements in the vector space for the response g.
   * The optionally supplied parameter list supports the following options:
   * <ul>
   *    <li> "Sensitivity Parameter Index", (default: 0) The model evaluator
   *          parameter index for which sensitivities will be computed.
   *    <li> "Response Function Index", (default: 0) The model evaluator
   *         response index for which sensitivities will be computed.
   *    <li> "Mass Matrix Is Constant" (default: true) Whether the mass matrix
   *         df/dx_dot is a constant matrix.  As describe above, this is
   *         currently required to be true.
   *    <li> "Mass Matrix Is Identity" (default: false) Whether the mass matrix
   *         is the identity matrix, in which some computations can be skipped.
   * </ul>
   */
  AdjointSensitivityModelEvaluator(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >&
          adjoint_residual_model,
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >&
          adjoint_solve_model,
      const Scalar& t_init, const Scalar& t_final,
      const bool is_pseudotransient,
      const Teuchos::RCP<const Teuchos::ParameterList>& pList = Teuchos::null);

  //! Get the underlying model 'f'
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > getModel() const
  {
    return model_;
  }

  //! Get the underlying adjoint residual model
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > getAdjointResidualModel()
      const
  {
    return adjoint_residual_model_;
  }

  //! Get the underlying adjoint solve model
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > getAdjointSolveModel()
      const
  {
    return adjoint_solve_model_;
  }

  //! Set the final time from the forward evaluation
  void setFinalTime(const Scalar t_final);

  //! Set solution history from forward evaluation
  void setForwardSolutionHistory(
      const Teuchos::RCP<const Tempus::SolutionHistory<Scalar> >& sh);

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_space(int p) const;

  Teuchos::RCP<const Teuchos::Array<std::string> > get_p_names(int p) const;

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_space() const;

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_g_space(int j) const;

  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_W_op() const;

  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
  get_W_factory() const;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;

  //@}

  static Teuchos::RCP<const Teuchos::ParameterList> getValidParameters();

 private:
  typedef Thyra::DefaultMultiVectorProductVectorSpace<Scalar> DMVPVS;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;

  void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> prototypeInArgs_;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> prototypeOutArgs_;

  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > model_;
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > adjoint_residual_model_;
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > adjoint_solve_model_;

  Teuchos::RCP<const DMVPVS> adjoint_space_;
  Teuchos::RCP<const DMVPVS> residual_space_;
  Teuchos::RCP<const DMVPVS> response_space_;
  Teuchos::RCP<const Tempus::SolutionHistory<Scalar> > sh_;
  Scalar t_init_;
  Scalar t_final_;
  bool is_pseudotransient_;
  bool mass_matrix_is_constant_;
  bool mass_matrix_is_identity_;
  int p_index_;
  int g_index_;
  int num_adjoint_;

  mutable bool mass_matrix_is_computed_;
  mutable bool jacobian_matrix_is_computed_;
  mutable bool response_gradient_is_computed_;
  mutable Teuchos::RCP<Thyra::VectorBase<Scalar> > my_x_dot_;
  mutable Teuchos::RCP<Thyra::LinearOpBase<Scalar> > my_dfdx_;
  mutable Teuchos::RCP<Thyra::LinearOpBase<Scalar> > my_dfdxdot_;
  mutable Teuchos::RCP<Thyra::LinearOpBase<Scalar> > my_dfdp_op_;
  mutable Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > my_dfdp_mv_;
  mutable Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > my_dgdx_mv_;
  mutable Teuchos::RCP<Tempus::SolutionState<Scalar> > forward_state_;
  mutable Scalar t_interp_;
};

}  // namespace Tempus

#endif
