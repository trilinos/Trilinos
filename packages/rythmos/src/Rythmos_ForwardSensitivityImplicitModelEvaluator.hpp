//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef RYTHMOS_FORWARD_SENSITIVITY_IMPLICIT_MODEL_EVALUATOR_HPP
#define RYTHMOS_FORWARD_SENSITIVITY_IMPLICIT_MODEL_EVALUATOR_HPP


#include "Rythmos_IntegratorBase.hpp"
#include "Rythmos_ForwardSensitivityModelEvaluatorBase.hpp"
#include "Rythmos_SolverAcceptingStepperBase.hpp"
#include "Rythmos_SingleResidualModelEvaluator.hpp"
#include "Thyra_ModelEvaluator.hpp" // Interface
#include "Thyra_StateFuncModelEvaluatorBase.hpp" // Implementation
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_PhysicallyBlockedLinearOpWithSolveBase.hpp" // Interface
#include "Thyra_DefaultBlockedTriangularLinearOpWithSolve.hpp" // Implementation
#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"
#include "Thyra_DefaultMultiVectorLinearOpWithSolve.hpp"
#include "Teuchos_implicit_cast.hpp"
#include "Teuchos_Assert.hpp"


namespace Rythmos {


/** \brief Forward sensitivity transient <tt>ModelEvaluator</tt> subclass.
 *
 * This class provides a very general implemenation of a linear forward
 * sensitivity model evaluator for a DAE.
 *
 * \section Rythmos_ForwardSensitivityModelEvaluatorBase_intro_sec Introduction
 *
 * The form of the parameterized state equation is:

 \verbatim

   f(x_dot(t),x(t),{p_l},t) = 0, over t = [t0,tf]

   x(t0) = x_init(p)

 \endverbatim

 * As shown above, the parameters are assumed to be steady state and can enter
 * through the intial condition and/or through the DAE equation itself.
 *
 * The forward sensitivity equations written in multi-vector form are:

 \verbatim

   F_sens(S_dot,S,t) = d(f)/d(x_dot)*S_dot + d(f)/d(x)*S + d(f)/d(p) = 0, over t = [t0,tf]

   S(t0) = d(x_init)/d(p)

 \endverbatim

 * where <tt>S</tt> is a multi-vector with <tt>np</tt> columns where each
 * column <tt>S(:,j) = d(x)/d(p_j)</tt> is the sensitivity of <tt>x(t)</tt>
 * with respect to the <tt>p_j</tt> parameter.  The sensitivity parameter
 * subvector <tt>p</tt> here is really just one of the parameter subvectors in
 * the state equation.  This index of the parameter subvector for which the
 * sensitivity equations are written for is given by <tt>p_index</tt>.  Note
 * that above <tt>d(f)/d(x_dot)</tt>, <tt>d(f)/d(x)</tt> and <tt>d(f)/d(p</tt>
 * are all evaluated at the solution to the state equation
 * <tt>(x_dot(t),x(t),t)</tt> and are not functions of <tt>S_dot</tt> or
 * <tt>S</tt>.
 *
 * Since the model evaluator interface must be expressed in vector form, the
 * multi-vector form of the forward sensitivity equations is flattened out
 * into:

 \verbatim

   f_sens(s_bar_dot(t),s_bar(t),{p_l},t) = 0, over t = [t0,tf]

   s_bar(t0) = s_bar_init

 \endverbatim

 * where

 \verbatim

   s_bar = [ S(:,0); S(:,0); ...; S(:,np-1) ]

            [ d(f)/d(x_dot)*S_dot(:,0) + d(f)/d(x)*S(:,0) + d(f)/d(p(0))          ]
            [ d(f)/d(x_dot)*S_dot(:,1) + d(f)/d(x)*S(:,1) + d(f)/d(p(1))          ]
   f_sens = [ ...                                                                 ]
            [ d(f)/d(x_dot)*S_dot(:,np-1) + d(f)/d(x)*S(:,np-1) + d(f)/d(p(np-1)) ]

   s_bar_init = [ d(x_init)/d(p(0)); d(x_init)/d(p(1)); ...; d(x_init)/d(p(np-1)) ]

 \endverbatim

 * The product vector <tt>s_bar</tt> is represented as a specialized
 * <tt>Thyra::ProductVectorBase</tt> subclass object with <tt>np</tt> "blocks"
 * in terms of a single <tt>Thyra::MultiVectorBase</tt> object (which has
 * <tt>np</tt> columns).
 *
 * \section Rythmos_ForwardSensitivityModelEvaluatorBase_details_sec Implementation Details
 *
 * In order to achieve high performance, the representation of the forward
 * sensitivity equations and the computation of the timestep for the
 * sensitivities must reuse as much from the state solve as possible.  Here we
 * are most concerned about implicit time stepping methods which compute and
 * use the composite matrix

 \verbatim
   W = alpha*d(f)/d(x_dot) + beta*d(f)/d(x)
 \endverbatim

 * which is the Jacobian for the nonlinear timestep equation in many methods.
 *
 * First, let us consider how to represent the forward sensitivity equations
 * using a precomputed LOWSB object

 \verbatim
   W_tilde = coeff_x_dot * d(f)/d(x_dot) + coeff_x * d(f)/d(x_dot)
 \endverbatim

 * computed at some point <tt>(x_dot_tilde,x_tilde,t_tilde,...)</tt>.
 *
 * Here is how to evaluate the forward sensitivity equations <tt>F_sens</tt>
 * using <tt>W_tilde</tt>:

 \verbatim

    F_sens = d(f)/d(x_dot)*S_dot
             + ( (1/coeff_x) * ( W_tilde - coeff_x_dot * d(f)/d(x_dot) ) ) * S
             + d(f)/d(p)

    ==>

    F_sens = (1/coeff_x) * W_tilde * S
             + d(f)/d(x_dot) * ( S_dot - (coeff_x_dot/coeff_x)*S )
             + d(f)/d(p)

 \endverbatim

 * Above, we see that in order to compute the multi-vector form of the
 * residual for the forward sensitivity equations <tt>F_sens</tt>, we must be
 * able to compute something like:

 \verbatim

    d(f)/d(x_dot) * V + d(f)/d(p)

 \endverbatim

 * This can be done by computing <tt>d(f)/d(x_dot)</tt> as <tt>W</tt> with
 * <tt>alpha=1.0</tt> and <tt>beta=0.0</tt> by calling the underlying state
 * model.  Or, a special sensitivity computation could be added to the model
 * evaluator that would compute a generalization of:

 \verbatim

    F_sens_var = x_dot_scalar * d(f)/d(x_dot) * V_x_dot
                   + x_scalar * d(f)/d(x) * V_x
                   + p_scalar * d(f)/d(p) * V_p

 \endverbatim

 * We could then compute what we need using <tt>x_dot_scalar=1.0</tt>,
 * <tt>x_scalar=0.0</tt>, <tt>p_scalar=1.0</tt>, <tt>V_x_dot=V</tt>, and
 * <tt>V_p=I</tt>.  For explicit time-stepping methods, this operation could
 * compute the needed sensitivity in one shot.  Such an addition to the
 * <tt>Thyra::ModelEvaluator</tt> interface would be handled through additions
 * to the InArgs and OutArgs classes and the details of which are yet to be
 * worked out.
 *
 * Up to this point, the only assumption that we have made about the time
 * stepping algorithm is that the sensitivities are only needed and only
 * represented at a single point <tt>(x_dot_tilde,x_tilde,t_tilde,...)</tt>.
 * It is up to the client to ensure that this is indeed the case and it will
 * be the case for many different types of timestepping methods.
 *
 * In order to get maximum reuse out of a precomputed <tt>W_tilde</tt> matrix
 * we also have to put in special logic for the evaluation of <tt>W_hat</tt>
 * with respect to the forward sensitivity residual.  When using a precomputed
 * <tt>W_tilde</tt> from the state timestep computation, the matrix
 * <tt>W_hat</tt> for the sensitivity residual becomes:

 \verbatim

    W_hat = alpha * d(f)/d(x_dot)
            + beta * (1/coeff_x) * ( W_tilde - coeff_x_dot * d(f)/d(x_dot) )

 \endverbatim

 * The special logic is that when <tt>alpha=coeff_x_dot</tt> and
 * <tt>beta=coeff_x</tt>, then <tt>W_hat</tt> above simplifies to
 * <tt>W_hat=W_tilde</tt> as:


 \verbatim

    W_hat = coeff_x_dot * d(f)/d(x_dot)
            + coeff_x * (1/coeff_x) * ( W_tilde - coeff_x_dot * d(f)/d(x_dot) )

    ==>

    W_hat = coeff_x_dot * d(f)/d(x_dot) - coeff_x_dot * d(f)/d(x_dot)
            + W_tilde

    ==>

    W_hat = W_tilde

 \endverbatim

 * For all of the timestepping methods currently implemented in Rythmos at the
 * time of this writing, <tt>alpha=coeff_x_dot</tt> and <tt>beta=coeff_x</tt>
 * will always be true and this is checked for in the code implementation.
 * Because of this, only a single <tt>W</tt> LOWSB object will be given to a
 * client and only it will be returned.  Any other values of <tt>alpha</tt>
 * and <tt>beta</tt> requested by the client will thown exceptions.  In the
 * future, other values of <tt>alpha</tt> and <tt>beta</tt> might be allowed
 * but for now this would only be an error.
 *
 * Note that there are also a few simplifications in cases where single
 * residual timestepping methods are used.  In the case of BDF methods (of
 * which backward Euler is one type), we have:

 \verbatim

    S_dot = coeff_x_dot * S_tilde + B_x_dot

    S_dot = coeff_x * S_tilde
    
 \endverbatim

 * In this type of method we see that :

 \verbatim

    S_dot - (coeff_x_dot/coeff_x) * S

    ==> 

    coeff_x_dot * S_tilde + B_x_dot - (coeff_x_dot/coeff_x) * ( coeff_x * S_tilde )

    ==> 

    coeff_x_dot * S_tilde - coeff_x_dot *  S_tilde + B_x_dot

    ==>

    B_x_dot
    
 \endverbatim

 * Therefore, in this type of method, the term involving
 * <tt>d(f)/d(x_dot)</tt> becomes:

 \verbatim

    d(f)/d(x_dot) * ( S_dot - (coeff_x_dot/coeff_x)*S )

    ==>

    d(f)/d(x_dot) * B_x_dot

 \endverbatim
 
 * and is independent of the unknown quantity <tt>S_tilde</tt>.  What this
 * means is that if the residual for the sensitivity equaitions is to be
 * computed multiple times for different values of <tt>S_tilde</tt>, the term
 * <tt>d(f)/d(x_dot) * B_x_dot</tt> need only be computed once and can then be
 * reused each time.
 *
 * ToDo: Finish documention!
 */
template<class Scalar>
class ForwardSensitivityImplicitModelEvaluator
  : virtual public ForwardSensitivityModelEvaluatorBase<Scalar>
{
public:

  /** \name Constructors/Intializers/Accessors */
  //@{

  /** \brief . */
  ForwardSensitivityImplicitModelEvaluator();

  /** \brief Intialize the with the model structure.
   *
   * \param  stateModel
   *           [in,persisting] The ModelEvaluator that defines the
   *           parameterized state model <tt>f(x_dot,x,p)</tt>.
   * \param  p_index
   *           [in] The index of the parameter subvector in <tt>stateModel</tt>
   *           for which sensitivities will be computed for.
   *
   * This function only intializes the spaces etc. needed to define structure
   * of the problem.  <tt>*this</tt> model object is not fully initialized at
   * this point in that <tt>evalModel()</tt> will not work yet and will thrown
   * exceptions if called.  The function <tt>initalizeState()</tt> must be
   * called later in order to fully initalize the model.
   */
  void initializeStructure(
    const RCP<const Thyra::ModelEvaluator<Scalar> > &stateModel,
    const int p_index
    );
  
  /** \brief . */
  RCP<const Thyra::ModelEvaluator<Scalar> >
  getStateModel() const;

  /** \brief . */
  RCP<Thyra::ModelEvaluator<Scalar> >
  getNonconstStateModel() const;
  
  /** \brief . */
  int get_p_index() const;

  /** \brief Set the state integrator that will be used to get x and x_dot at
   * various time points.
   *
   * \param stateIntegrator [in,persisting] The integrator that will deliver x
   * and xdot at any valid time t.  This integrator must be completely step up
   * and ready to return points from stateIntegrator->getFwdPoints(...).
   *
   * \param stateBasePoint [in,persisting] Contains all of the base point
   * information for an evaluation except for <tt>t</tt>, <tt>x</tt>, and
   * <tt>x_dot</tt>.  For exmaple, this will contain any parameters or extra
   * input that defines the point.
   *
   * With a state integrator is set here, the client need not call
   * <tt>initializePointState()</tt> for individual time points!
   */
  void setStateIntegrator(
    const RCP<IntegratorBase<Scalar> > &stateIntegrator,
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &stateBasePoint
    );

  /** \brief Initialize full state for a single point in time.
   *
   * \param stateBasePoint [in] The base point
   * <tt>(x_dot_tilde,x_tilde,t_tilde,...)</tt> for which the sensitivities
   * will be computed and represented at time <tt>t_tilde</tt>.  Any
   * sensitivities that are needed must be computed at this same time point.
   * The value of <tt>alpha</tt> and <tt>beta</tt> here are ignored.
   *
   * \param W_tilde [in,persisting] The composite state derivative matrix
   * computed at the base point <tt>stateBasePoint</tt> with
   * <tt>alpha=coeff_x_dot</tt> and <tt>beta=coeff_x</tt>.  This argument can
   * be null, in which case it can be computed internally if needed or not at
   * all.
   *
   * \param coeff_x_dot [in] The value of <tt>alpha</tt> for which
   * <tt>W_tilde</tt> was computed.
   *
   * \param coeff_x [in] The value of <tt>beta</tt> for which <tt>W_tilde</tt>
   * was computed.
   *
   * \param DfDx_dot [in] The matrix <tt>d(f)/d(x_dot)</tt> computed at
   * <tt>stateBasePoint</tt> if available.  If this argument is null, then it
   * will be computed internally if needed.  The default value is
   * <tt>Teuchos::null</tt>.
   *
   * \param DfDp [in] The matrix <tt>d(f)/d(p(p_index))</tt> computed at
   * <tt>stateBasePoint</tt> if available.  If this argument is null, then it
   * will be computed internally if needed.  The default value is
   * <tt>Teuchos::null</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>nonnull(W_tilde)</tt>
   * </ul>
   *
   * This function must be called after <tt>intializeStructure()</tt> and
   * before <tt>evalModel()</tt> is called.  After this function is called,
   * <tt>*this</tt> model is fully initialized and ready to compute the
   * requested outputs.
   */
  void initializePointState(
      Ptr<StepperBase<Scalar> > stateStepper,
      bool forceUpToDateW
      );

//  /** \brief Initialize full state for a single point in time.
//   *
//   * \param stateBasePoint [in] The base point
//   * <tt>(x_dot_tilde,x_tilde,t_tilde,...)</tt> for which the sensitivities
//   * will be computed and represented at time <tt>t_tilde</tt>.  Any
//   * sensitivities that are needed must be computed at this same time point.
//   * The value of <tt>alpha</tt> and <tt>beta</tt> here are ignored.
//   *
//   * \param W_tilde [in,persisting] The composite state derivative matrix
//   * computed at the base point <tt>stateBasePoint</tt> with
//   * <tt>alpha=coeff_x_dot</tt> and <tt>beta=coeff_x</tt>.  This argument can
//   * be null, in which case it can be computed internally if needed or not at
//   * all.
//   *
//   * \param coeff_x_dot [in] The value of <tt>alpha</tt> for which
//   * <tt>W_tilde</tt> was computed.
//   *
//   * \param coeff_x [in] The value of <tt>beta</tt> for which <tt>W_tilde</tt>
//   * was computed.
//   *
//   * \param DfDx_dot [in] The matrix <tt>d(f)/d(x_dot)</tt> computed at
//   * <tt>stateBasePoint</tt> if available.  If this argument is null, then it
//   * will be computed internally if needed.  The default value is
//   * <tt>Teuchos::null</tt>.
//   *
//   * \param DfDp [in] The matrix <tt>d(f)/d(p(p_index))</tt> computed at
//   * <tt>stateBasePoint</tt> if available.  If this argument is null, then it
//   * will be computed internally if needed.  The default value is
//   * <tt>Teuchos::null</tt>.
//   *
//   * <b>Preconditions:</b><ul>
//   * <li><tt>!is_null(W_tilde)</tt>
//   * </ul>
//   *
//   * This function must be called after <tt>intializeStructure()</tt> and
//   * before <tt>evalModel()</tt> is called.  After this function is called,
//   * <tt>*this</tt> model is fully initialized and ready to compute the
//   * requested outputs.
//   */
//  void initializePointState(
//    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &stateBasePoint,
//    const RCP<const Thyra::LinearOpWithSolveBase<Scalar> > &W_tilde,
//    const Scalar &coeff_x_dot,
//    const Scalar &coeff_x,
//    const RCP<const Thyra::LinearOpBase<Scalar> > &DfDx_dot = Teuchos::null,
//    const RCP<const Thyra::MultiVectorBase<Scalar> > &DfDp = Teuchos::null
//    );

  // 2007/05/22: rabartl: ToDo: Add an InterpolationBufferBase
  // stateInterpBuffer object to the initailizeState(...) function that can be
  // used to get x and x_dot at different points in time t.  Then, modify the
  // logic to recompute all of the needed matrices if t != t_base (as passed
  // in through stateBasePoint).  The values of x(t) and xdot(t) can then be
  // gotten from the stateInterpBuffer object!
  
  //@}

  /** \name Public functions overridden from ForwardSensitivityModelEvaluatorBase. */
  //@{

  /** \brief . */
  void initializeStructureInitCondOnly(
    const RCP<const Thyra::ModelEvaluator<Scalar> >& stateModel,
    const RCP<const Thyra::VectorSpaceBase<Scalar> >& p_space
    );
  
  /** \brief . */
  RCP<const Thyra::DefaultMultiVectorProductVectorSpace<Scalar> >
  get_s_bar_space() const;
  
  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_sens_space() const;

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;
  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_space() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief . */
  RCP<Thyra::LinearOpWithSolveBase<Scalar> > create_W() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

  //@}

  /** \name Deprecated. */
  //@{

  /** \brief Deprecated. */
  void initializeState(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &stateBasePoint,
    const RCP<const Thyra::LinearOpWithSolveBase<Scalar> > &W_tilde,
    const Scalar &coeff_x_dot,
    const Scalar &coeff_x,
    const RCP<const Thyra::LinearOpBase<Scalar> > &DfDx_dot = Teuchos::null,
    const RCP<const Thyra::MultiVectorBase<Scalar> > &DfDp = Teuchos::null
    );

  //@}

private:

  /** \name Private functions overridden from ModelEvaulatorDefaultBase. */
  //@{

  /** \brief . */
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  /** \brief . */
  void evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
    ) const;

  //@}

private:

  // /////////////////////////
  // Private data members

  RCP<const Thyra::ModelEvaluator<Scalar> > stateModel_;
  int p_index_;
  RCP<const Thyra::VectorSpaceBase<Scalar> > p_space_;
  int np_;

  RCP<IntegratorBase<Scalar> > stateIntegrator_;

  RCP<const Thyra::DefaultMultiVectorProductVectorSpace<Scalar> > s_bar_space_;
  RCP<const Thyra::DefaultMultiVectorProductVectorSpace<Scalar> > f_sens_space_;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues_;

  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar> stateBasePoint_;

  mutable RCP<const Thyra::LinearOpWithSolveBase<Scalar> > W_tilde_;
  mutable Scalar coeff_x_dot_;
  mutable Scalar coeff_x_;
  mutable RCP<const Thyra::LinearOpBase<Scalar> > DfDx_dot_;
  mutable RCP<const Thyra::MultiVectorBase<Scalar> > DfDp_;

  mutable RCP<Thyra::LinearOpWithSolveBase<Scalar> > W_tilde_compute_;
  mutable RCP<Thyra::LinearOpBase<Scalar> > DfDx_dot_compute_;
  mutable RCP<Thyra::MultiVectorBase<Scalar> > DfDp_compute_;

  // /////////////////////////
  // Private member functions

  bool hasStateFuncParams() const { return p_index_ >= 0; }
  
  void initializeStructureCommon(
    const RCP<const Thyra::ModelEvaluator<Scalar> > &stateModel,
    const int p_index,
    const RCP<const Thyra::VectorSpaceBase<Scalar> > &p_space
    );

  void wrapNominalValuesAndBounds();

  void computeDerivativeMatrices(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &point
    ) const;
  
};


// /////////////////////////////////
// Implementations


// Constructors/Intializers/Accessors

template<class Scalar>
RCP<ForwardSensitivityImplicitModelEvaluator<Scalar> > forwardSensitivityImplicitModelEvaluator()
{
  return rcp(new ForwardSensitivityImplicitModelEvaluator<Scalar>());
}

template<class Scalar>
ForwardSensitivityImplicitModelEvaluator<Scalar>::ForwardSensitivityImplicitModelEvaluator()
  : p_index_(0), np_(-1)
{}


template<class Scalar>
RCP<const Thyra::ModelEvaluator<Scalar> >
ForwardSensitivityImplicitModelEvaluator<Scalar>::getStateModel() const
{
  return stateModel_;
}


template<class Scalar>
RCP<Thyra::ModelEvaluator<Scalar> >
ForwardSensitivityImplicitModelEvaluator<Scalar>::getNonconstStateModel() const
{
  return Teuchos::null;
}


template<class Scalar>
int ForwardSensitivityImplicitModelEvaluator<Scalar>::get_p_index() const
{
  return p_index_;
}


template<class Scalar>
void ForwardSensitivityImplicitModelEvaluator<Scalar>::setStateIntegrator(
  const RCP<IntegratorBase<Scalar> > &stateIntegrator,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &stateBasePoint
  )
{
  stateIntegrator_ = stateIntegrator.assert_not_null();
  stateBasePoint_ = stateBasePoint;
}

template<class Scalar>
void ForwardSensitivityImplicitModelEvaluator<Scalar>::initializePointState(
      Ptr<StepperBase<Scalar> > stateStepper,
      bool forceUpToDateW
      )
{
  TEUCHOS_ASSERT( Teuchos::nonnull(stateStepper) );
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();

  const StepStatus<Scalar> stateStepStatus = stateStepper->getStepStatus();
  TEST_FOR_EXCEPTION(
      stateStepStatus.stepStatus != STEP_STATUS_CONVERGED, std::logic_error,
      "Error, the status should be converged since a positive step size was returned!"
      );
  Scalar curr_t = stateStepStatus.time;

  // Get both x and x_dot since these are needed compute other derivative
  // objects at these points.
  RCP<const Thyra::VectorBase<Scalar> > x, x_dot;
  get_x_and_x_dot(*stateStepper,curr_t,&x,&x_dot);
      
  stateBasePoint_ = stateStepper->getInitialCondition();
  stateBasePoint_.set_x_dot( x_dot );
  stateBasePoint_.set_x( x );
  stateBasePoint_.set_t( curr_t );

  // Grab the SingleResidualModel that was used to compute the state timestep.
  // From this, we can get the constants that where used to compute W!
  RCP<SolverAcceptingStepperBase<Scalar> > 
    sasStepper = Teuchos::rcp_dynamic_cast<SolverAcceptingStepperBase<Scalar> >(
        Teuchos::rcpFromRef(*stateStepper),true
        );
  RCP<Thyra::NonlinearSolverBase<Scalar> > 
    stateTimeStepSolver = sasStepper->getNonconstSolver();
  RCP<const SingleResidualModelEvaluatorBase<Scalar> >
    singleResidualModel
    = Teuchos::rcp_dynamic_cast<const SingleResidualModelEvaluatorBase<Scalar> >(
        stateTimeStepSolver->getModel(), true
        );
  const Scalar
    coeff_x_dot = singleResidualModel->get_coeff_x_dot(),
    coeff_x = singleResidualModel->get_coeff_x();

  // Get W (and force an update if not up to date already)

  using Teuchos::as;
  if ((as<int>(verbLevel) >= as<int>(Teuchos::VERB_MEDIUM)) && forceUpToDateW) {
    *out << "\nForcing an update of W at the converged state timestep ...\n";
  }

  RCP<Thyra::LinearOpWithSolveBase<Scalar> >
    W_tilde = stateTimeStepSolver->get_nonconst_W(forceUpToDateW);

  TEST_FOR_EXCEPTION(
      is_null(W_tilde), std::logic_error,
      "Error, the W from the state time step must be non-null!"
      );

#ifdef RYTHMOS_DEBUG
  TEST_FOR_EXCEPTION(
    is_null(stateModel_), std::logic_error,
    "Error, you must call intializeStructure(...) before you call initializeState(...)"
    );
  TEST_FOR_EXCEPT( is_null(stateBasePoint_.get_x()) );
  TEST_FOR_EXCEPT( is_null(stateBasePoint_.get_x_dot()) );
  TEST_FOR_EXCEPT( is_null(stateBasePoint_.get_p(p_index_)) );
  // What about the other parameter values?  We really can't say anything
  // about these and we can't check them.  They can be null just fine.
  if (!is_null(W_tilde)) {
    THYRA_ASSERT_VEC_SPACES("",*W_tilde->domain(),*stateModel_->get_x_space());
    THYRA_ASSERT_VEC_SPACES("",*W_tilde->range(),*stateModel_->get_f_space());
  }
#endif

  W_tilde_ = W_tilde;
  coeff_x_dot_ = coeff_x_dot;
  coeff_x_ = coeff_x;
  DfDx_dot_ = Teuchos::null;
  DfDp_ = Teuchos::null;

  wrapNominalValuesAndBounds();

}


template<class Scalar>
void ForwardSensitivityImplicitModelEvaluator<Scalar>::initializeState(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &stateBasePoint,
  const RCP<const Thyra::LinearOpWithSolveBase<Scalar> > &W_tilde,
  const Scalar &coeff_x_dot,
  const Scalar &coeff_x,
  const RCP<const Thyra::LinearOpBase<Scalar> > &DfDx_dot,
  const RCP<const Thyra::MultiVectorBase<Scalar> > &DfDp
  )
{

  typedef Thyra::ModelEvaluatorBase MEB;

#ifdef RYTHMOS_DEBUG
  TEST_FOR_EXCEPTION(
    is_null(stateModel_), std::logic_error,
    "Error, you must call intializeStructure(...) before you call initializeState(...)"
    );
  TEST_FOR_EXCEPT( is_null(stateBasePoint.get_x()) );
  TEST_FOR_EXCEPT( is_null(stateBasePoint.get_x_dot()) );
  if (hasStateFuncParams()) {
    TEST_FOR_EXCEPT( is_null(stateBasePoint.get_p(p_index_)) );
  }
  // What about the other parameter values?  We really can't say anything
  // about these and we can't check them.  They can be null just fine.
  if (nonnull(W_tilde)) {
    THYRA_ASSERT_VEC_SPACES("", *W_tilde->domain(), *stateModel_->get_x_space());
    THYRA_ASSERT_VEC_SPACES("", *W_tilde->range(), *stateModel_->get_f_space());
  }
  if (nonnull(DfDx_dot)) {
    THYRA_ASSERT_VEC_SPACES("", *DfDx_dot->domain(), *stateModel_->get_x_space());
    THYRA_ASSERT_VEC_SPACES("", *DfDx_dot->range(), *stateModel_->get_f_space());
  }
  if (nonnull(DfDp)) {
    TEUCHOS_ASSERT(hasStateFuncParams());
    THYRA_ASSERT_VEC_SPACES("", *DfDp->domain(), *p_space_);
    THYRA_ASSERT_VEC_SPACES("", *DfDp->range(), *stateModel_->get_f_space());
  }
#endif

  stateBasePoint_ = stateBasePoint;

  // Set whatever derivatives where passed in.  If an input in null, then the
  // member will be null and the null linear operators will be computed later
  // just in time.

  W_tilde_ = W_tilde;
  coeff_x_dot_ = coeff_x_dot;
  coeff_x_ = coeff_x;
  DfDx_dot_ = DfDx_dot;
  DfDp_ = DfDp;

  wrapNominalValuesAndBounds();

}


// Public functions overridden from ForwardSensitivityModelEvaluatorBase


template<class Scalar>
void ForwardSensitivityImplicitModelEvaluator<Scalar>::initializeStructure(
  const RCP<const Thyra::ModelEvaluator<Scalar> >& stateModel,
  const int p_index
  )
{
  initializeStructureCommon(stateModel, p_index, Teuchos::null);
}


template<class Scalar>
void ForwardSensitivityImplicitModelEvaluator<Scalar>::initializeStructureInitCondOnly(
  const RCP<const Thyra::ModelEvaluator<Scalar> >& stateModel,
  const RCP<const Thyra::VectorSpaceBase<Scalar> >& p_space
  )
{
  initializeStructureCommon(stateModel, -1, p_space);
}


template<class Scalar>
RCP<const Thyra::DefaultMultiVectorProductVectorSpace<Scalar> >
ForwardSensitivityImplicitModelEvaluator<Scalar>::get_s_bar_space() const
{
  return s_bar_space_;
}


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
ForwardSensitivityImplicitModelEvaluator<Scalar>::get_p_sens_space() const
{
  return p_space_;
}


// Public functions overridden from ModelEvaulator


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
ForwardSensitivityImplicitModelEvaluator<Scalar>::get_x_space() const
{
  return s_bar_space_;
}


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
ForwardSensitivityImplicitModelEvaluator<Scalar>::get_f_space() const
{
  return f_sens_space_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
ForwardSensitivityImplicitModelEvaluator<Scalar>::getNominalValues() const
{
  return nominalValues_;
}


template<class Scalar>
RCP<Thyra::LinearOpWithSolveBase<Scalar> >
ForwardSensitivityImplicitModelEvaluator<Scalar>::create_W() const
{
  return Thyra::multiVectorLinearOpWithSolve<Scalar>();
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
ForwardSensitivityImplicitModelEvaluator<Scalar>::createInArgs() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgs<Scalar> stateModelInArgs = stateModel_->createInArgs();
  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports( MEB::IN_ARG_x_dot,
    stateModelInArgs.supports(MEB::IN_ARG_x_dot) );
  inArgs.setSupports( MEB::IN_ARG_x );
  inArgs.setSupports( MEB::IN_ARG_t );
  inArgs.setSupports( MEB::IN_ARG_alpha,
    stateModelInArgs.supports(MEB::IN_ARG_alpha) );
  inArgs.setSupports( MEB::IN_ARG_beta,
    stateModelInArgs.supports(MEB::IN_ARG_beta) );
  return inArgs;
}


// Private functions overridden from ModelEvaulatorDefaultBase


template<class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
ForwardSensitivityImplicitModelEvaluator<Scalar>::createOutArgsImpl() const
{

  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::OutArgs<Scalar> stateModelOutArgs = stateModel_->createOutArgs();
  MEB::OutArgsSetup<Scalar> outArgs;

  outArgs.setModelEvalDescription(this->description());

  outArgs.setSupports(MEB::OUT_ARG_f);

  if (stateModelOutArgs.supports(MEB::OUT_ARG_W) ) {
    outArgs.setSupports(MEB::OUT_ARG_W);
    outArgs.set_W_properties(stateModelOutArgs.get_W_properties());
  }

  return outArgs;

}


template<class Scalar>
void ForwardSensitivityImplicitModelEvaluator<Scalar>::evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  ) const
{

  using Teuchos::rcp_dynamic_cast;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::VerboseObjectTempState<Thyra::ModelEvaluatorBase> VOTSME;

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_GEN_BEGIN(
    "ForwardSensitivityImplicitModelEvaluator", inArgs, outArgs, Teuchos::null );

  //
  // Update the derivative matrices if they are not already updated for the
  // given time!.
  //
  
  {
#ifdef ENABLE_RYTHMOS_TIMERS
    TEUCHOS_FUNC_TIME_MONITOR(
      "Rythmos:ForwardSensitivityImplicitModelEvaluator::evalModel: computeMatrices");
#endif
    computeDerivativeMatrices(inArgs);
  }

  //
  // InArgs
  //

  RCP<const Thyra::DefaultMultiVectorProductVector<Scalar> >
    s_bar = rcp_dynamic_cast<const Thyra::DefaultMultiVectorProductVector<Scalar> >(
      inArgs.get_x().assert_not_null(), true
      );
  RCP<const Thyra::DefaultMultiVectorProductVector<Scalar> >
    s_bar_dot = rcp_dynamic_cast<const Thyra::DefaultMultiVectorProductVector<Scalar> >(
      inArgs.get_x_dot().assert_not_null(), true
      );
  const Scalar
    alpha = inArgs.get_alpha();
  const Scalar
    beta = inArgs.get_beta();

  RCP<const Thyra::MultiVectorBase<Scalar> >
    S = s_bar->getMultiVector();
  RCP<const Thyra::MultiVectorBase<Scalar> >
    S_dot = s_bar_dot->getMultiVector();
  
  //
  // OutArgs
  //

  RCP<Thyra::DefaultMultiVectorProductVector<Scalar> >
    f_sens = rcp_dynamic_cast<Thyra::DefaultMultiVectorProductVector<Scalar> >(
      outArgs.get_f(), true
      );
  RCP<Thyra::DefaultMultiVectorLinearOpWithSolve<Scalar> >
    W_sens = rcp_dynamic_cast<Thyra::DefaultMultiVectorLinearOpWithSolve<Scalar> >(
      outArgs.get_W(), true
      );

  RCP<Thyra::MultiVectorBase<Scalar> >
    F_sens = f_sens->getNonconstMultiVector().assert_not_null();

  //
  // Compute the requested functions
  //

  if(nonnull(F_sens)) {

#ifdef ENABLE_RYTHMOS_TIMERS
    TEUCHOS_FUNC_TIME_MONITOR(
      "Rythmos:ForwardSensitivityImplicitModelEvaluator::evalModel: computeSens");
#endif

    // S_diff =  -(coeff_x_dot/coeff_x)*S + S_dot
    RCP<Thyra::MultiVectorBase<Scalar> >
      S_diff = createMembers( stateModel_->get_x_space(), np_ );
    V_StVpV( &*S_diff, Scalar(-coeff_x_dot_/coeff_x_), *S, *S_dot );
    // F_sens = (1/coeff_x) * W_tilde * S
    Thyra::apply(
      *W_tilde_, Thyra::NOTRANS,
      *S, &*F_sens,
      Scalar(1.0/coeff_x_), ST::zero()
      );
    // F_sens += d(f)/d(x_dot) * S_diff
    Thyra::apply(
      *DfDx_dot_, Thyra::NOTRANS,
      *S_diff, &*F_sens,
      ST::one(), ST::one()
      );
    // F_sens += d(f)/d(p)
    if (hasStateFuncParams())
      Vp_V( &*F_sens, *DfDp_ );
  }
  
  if(nonnull(W_sens)) {
    TEST_FOR_EXCEPTION(
      alpha != coeff_x_dot_, std::logic_error,
      "Error, alpha="<<alpha<<" != coeff_x_dot="<<coeff_x_dot_
      <<" with difference = "<<(alpha-coeff_x_dot_)<<"!" );
    TEST_FOR_EXCEPTION(
      beta != coeff_x_, std::logic_error,
      "Error, beta="<<beta<<" != coeff_x="<<coeff_x_
      <<" with difference = "<<(beta-coeff_x_)<<"!" );
    W_sens->initialize( W_tilde_, s_bar_space_, f_sens_space_ );
    
  }

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_END();

}


// private


template<class Scalar>
void ForwardSensitivityImplicitModelEvaluator<Scalar>::initializeStructureCommon(
  const RCP<const Thyra::ModelEvaluator<Scalar> >& stateModel,
  const int p_index,
  const RCP<const Thyra::VectorSpaceBase<Scalar> >& p_space
  )
{

  typedef Thyra::ModelEvaluatorBase MEB;

  //
  // Validate input
  //

  TEUCHOS_ASSERT(nonnull(stateModel));
  TEUCHOS_ASSERT(p_index >= 0 || nonnull(p_space));
  if (p_index >= 0) {
    TEST_FOR_EXCEPTION(
      !( 0 <= p_index && p_index < stateModel->Np() ), std::logic_error,
      "Error, p_index does not fall in the range [0,"<<(stateModel->Np()-1)<<"]!" );
    // ToDo: Validate support for DfDp!
  }
  else {
    TEUCHOS_ASSERT_EQUALITY(p_index, -1);
  }

  //
  // Set the input objects
  //

  stateModel_ = stateModel;

  if (p_index >= 0) {
    p_index_ = p_index;
    p_space_ = stateModel_->get_p_space(p_index);
  }
  else {
    p_index_ = -1;
    p_space_ = p_space;
  }

  np_ = p_space_->dim();

  //
  // Create the structure of the model
  //

  s_bar_space_ = Thyra::multiVectorProductVectorSpace(
    stateModel_->get_x_space(), np_
    );

  f_sens_space_ = Thyra::multiVectorProductVectorSpace(
    stateModel_->get_f_space(), np_
    );

  nominalValues_ = this->createInArgs();

  //
  // Wipe out matrix storage
  //

  stateBasePoint_ = MEB::InArgs<Scalar>();
  W_tilde_ = Teuchos::null;
  coeff_x_dot_ = 0.0;
  coeff_x_ = 0.0;
  DfDx_dot_ = Teuchos::null;
  DfDp_ = Teuchos::null;
  W_tilde_compute_ = Teuchos::null;
  DfDx_dot_compute_ = Teuchos::null;
  DfDp_compute_ = Teuchos::null;

}


template<class Scalar>
void ForwardSensitivityImplicitModelEvaluator<Scalar>::wrapNominalValuesAndBounds()
{

  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::ModelEvaluatorBase MEB;

  // nominalValues_.clear(); // ToDo: Implement this!

  nominalValues_.set_t(stateModel_->getNominalValues().get_t());
  
  // 2007/05/22: rabartl: Note: Currently there is not much of a reason to set
  // an initial condition here since the initial condition for the
  // sensitivities is really being set in the ForwardSensitivityStepper
  // object!  In the future, a more general use of this class might benefit
  // from setting the initial condition here.

}


template<class Scalar>
void ForwardSensitivityImplicitModelEvaluator<Scalar>::computeDerivativeMatrices(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &point
  ) const
{

  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::VerboseObjectTempState<MEB> VOTSME;

  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();

  const Scalar
    t_base = stateBasePoint_.get_t(),
    t = point.get_t();

  TEUCHOS_ASSERT_EQUALITY( t , t_base );

  if (is_null(W_tilde_)) {
    TEST_FOR_EXCEPT_MSG(true, "ToDo: compute W_tilde from scratch!");
  }
  
  if ( is_null(DfDx_dot_) || is_null(DfDp_) ) {

    MEB::InArgs<Scalar> inArgs = stateBasePoint_;
    MEB::OutArgs<Scalar> outArgs = stateModel_->createOutArgs();
    
    Teuchos::RCP<Thyra::LinearOpBase<Scalar> > DfDx_dot_compute;
    if (is_null(DfDx_dot_)) {
      DfDx_dot_compute = stateModel_->create_W_op();
      inArgs.set_alpha(1.0);
      inArgs.set_beta(0.0);
      outArgs.set_W_op(DfDx_dot_compute);
    }

    Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > DfDp_compute;
    if (is_null(DfDp_) && hasStateFuncParams()) {
      DfDp_compute = Thyra::create_DfDp_mv(
        *stateModel_, p_index_,
        MEB::DERIV_MV_BY_COL
        ).getMultiVector();
      outArgs.set_DfDp(
        p_index_,
        MEB::Derivative<Scalar>(DfDp_compute,MEB::DERIV_MV_BY_COL)
        );
    }
    
    VOTSME stateModel_outputTempState(stateModel_,out,verbLevel);
    stateModel_->evalModel(inArgs,outArgs);
    if (nonnull(DfDx_dot_compute))
      DfDx_dot_ = DfDx_dot_compute;
    if (nonnull(DfDp_compute))
      DfDp_ = DfDp_compute;
  
  }

  TEST_FOR_EXCEPT_MSG( nonnull(stateIntegrator_),
    "ToDo: Update for using the stateIntegrator!" );

/* 2007/12/11: rabartl: ToDo: Update the code below to work for the general
 * case.  It does not work in its current form!
 */

/*

  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::VerboseObjectTempState<MEB> VOTSME;

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  
  const Scalar t = point.get_t();

  //
  // A) Set up the base point at time t and determine what needs to be
  // computed here.
  //

  bool update_W_tilde = false;
  bool update_DfDx_dot = false;
  bool update_DfDp = false;

  if (nonnull(stateIntegrator_)) {
    // Get x and x_dot at t and flag to recompute all matrices!
    RCP<const Thyra::VectorBase<Scalar> > x, x_dot;
    get_fwd_x_and_x_dot( *stateIntegrator_, t, &x, &x_dot );
    stateBasePoint_.set_t(t);
    stateBasePoint_.set_x(x);
    stateBasePoint_.set_x_dot(x_dot);
    update_W_tilde = true;
    update_DfDx_dot = true;
    update_DfDp = true;
  }
  else {
    // The time point should be the same so only compute matrices that have
    // not already been computed.
    TEUCHOS_ASSERT_EQUALITY( t , stateBasePoint_.get_t() );
    if (is_null(W_tilde_))
      update_W_tilde = true;
    if (is_null(DfDx_dot_))
      update_DfDx_dot = true;
    if (is_null(DfDp_))
      update_DfDx_dot = true;
  }

  //
  // B) Compute DfDx_dot and DfDp at the same time if needed
  //

  if ( update_DfDx_dot || update_DfDp ) {

    // B.1) Allocate storage for matrices.  Note: All of these matrices are
    // needed so any of these that is null must be coputed!

    if ( is_null(DfDx_dot_) || is_null(DfDx_dot_compute_) )
      DfDx_dot_compute_ = stateModel_->create_W_op();

    if ( is_null(DfDp_) || is_null(DfDp_compute_) )
      DfDp_compute_ = Thyra::create_DfDp_mv(
        *stateModel_,p_index_,
        MEB::DERIV_MV_BY_COL
        ).getMultiVector();

    // B.2) Setup the InArgs and OutArgs

    MEB::InArgs<Scalar> inArgs = stateBasePoint_;
    MEB::OutArgs<Scalar> outArgs = stateModel_->createOutArgs();
    
    if (update_DfDx_dot) {
      inArgs.set_alpha(1.0);
      inArgs.set_beta(0.0);
      outArgs.set_W_op(DfDx_dot_compute_);
    }
    // 2007/12/07: rabartl: ToDo: I need to change the structure of the
    // derivative properties in OutArgs to keep track of whether DfDx_dot is
    // constant or not separate from W.  Right now I can't do this.  This is
    // one reason to add a new DfDx_dot output object to OutArgs.  A better
    // solution is a more general data structure giving the dependence of f()
    // and g[j]() on x, x_dot, p[l], t, etc ...

    if (update_DfDp) {
      outArgs.set_DfDp(
        p_index_,
        MEB::Derivative<Scalar>(DfDp_compute_, MEB::DERIV_MV_BY_COL)
        );
    }

    // B.3) Evaluate the outputs
    
    VOTSME stateModel_outputTempState(stateModel_,out,verbLevel);
    stateModel_->evalModel(inArgs,outArgs);

    // B.4) Set the outputs

    if (nonnull(DfDx_dot_compute_))
      DfDx_dot_ = DfDx_dot_compute_;

    if (nonnull(DfDp_compute_))
      DfDp_ = DfDp_compute_;
  
  }

  //
  // C) Compute W_tilde separately if needed
  //

  if ( update_W_tilde ) {

    // C.1) Allocate storage for matrices.  Note: All of these matrices are
    // needed so any of these that is null must be coputed!

    if ( is_null(W_tilde_) || is_null(W_tilde_compute_) )
      W_tilde_compute_ = stateModel_->create_W();

    // C.2) Setup the InArgs and OutArgs

    MEB::InArgs<Scalar> inArgs = stateBasePoint_;
    MEB::OutArgs<Scalar> outArgs = stateModel_->createOutArgs();
    
    if (is_null(W_tilde_)) {
      coeff_x_dot_ = point.get_alpha();
      coeff_x_ = point.get_beta();
      inArgs.set_alpha(coeff_x_dot_);
      inArgs.set_beta(coeff_x_);
      outArgs.set_W(W_tilde_compute_);
    }

    // C.3) Evaluate the outputs
    
    VOTSME stateModel_outputTempState(stateModel_,out,verbLevel);
    stateModel_->evalModel(inArgs,outArgs);

    // B.4) Set the outputs

    if (nonnull(W_tilde_compute_))
      W_tilde_ = W_tilde_compute_;

  }

  // 2007/12/07: rabartl: Note: Above, we see an example of where it would be
  // good to have a separate OutArg for DfDx_dot (as a LOB object) so that we
  // can compute W and DfDx_dot (and any other output quantitiy) at the same
  // time.

*/


}


} // namespace Rythmos


#endif // RYTHMOS_FORWARD_SENSITIVITY_IMPLICIT_MODEL_EVALUATOR_HPP
