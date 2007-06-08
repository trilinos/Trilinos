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

#ifndef RYTHMOS_FORWARD_SENSITIVITY_MODEL_EVALUATOR_HPP
#define RYTHMOS_FORWARD_SENSITIVITY_MODEL_EVALUATOR_HPP


#include "Thyra_ModelEvaluator.hpp" // Interface
#include "Thyra_StateFuncModelEvaluatorBase.hpp" // Implementation
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_PhysicallyBlockedLinearOpWithSolveBase.hpp" // Interface
#include "Thyra_DefaultBlockedTriangularLinearOpWithSolve.hpp" // Implementation
#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"
#include "Thyra_DefaultMultiVectorLinearOpWithSolve.hpp"
#include "Teuchos_implicit_cast.hpp"


namespace Rythmos {


/** \brief Forward sensitivity transient <tt>ModelEvaluator</tt> subclass.
 *
 * This class provides a very general implemenation of a forward sensitivity
 * model evaluator for a DAE.
 *
 * \section Rythmos_ForwardSensitivityModelEvaluator_intro_sec Introduction
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
 * \section Rythmos_ForwardSensitivityModelEvaluator_details_sec Implementation Details
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

    x_dot_scalar * d(f)/d(x_dot) * V_x_dot
    + x_scalar * d(f)/d(x) * V_x
    + p_scalar * d(f)/d(p) * V_p

 \endverbatim

 * We could then compute what we need using <tt>x_dot_scalar=1.0</tt>,
 * <tt>x_scalar=0.0</tt>, <tt>p_scalar=1.0</tt>, <tt>V_x_dot=V</tt>, and
 * <tt>V_p=I</tt>.  For explicit time-stepping methods, this is all that would
 * be needed.  Such an addition to the <tt>Thyra::ModelEvaluator</tt>
 * interface would be handled through additions to the InArgs and OutArgs
 * classes and the details of which are yet to be worked out.
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

 * Therefore, in this types of method, the term involving
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
class ForwardSensitivityModelEvaluator
  : virtual public Thyra::ModelEvaluator<Scalar>, // Public interface
    virtual protected Thyra::StateFuncModelEvaluatorBase<Scalar> // Protected implementation
{
public:

  /** \name Constructors/Intializers/Accessors */
  //@{

  /** \brief . */
  ForwardSensitivityModelEvaluator();

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
    const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &stateModel,
    const int p_index
    );
  
  /** \brief . */
  Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> >
  getStateModel() const;
  
  /** \brief . */
  int get_p_index() const;

  /** \brief Initialize full state.
   *
   * \param  stateBasePoint
   *           [in] The base point <tt>(x_dot_tilde,x_tilde,t_tilde,...)</tt>
   *           for which the sensitivities will be computed and represented
   *           at.  Any sensitivities that are needed are computed at this
   *           point.  The value of <tt>alpha</tt> and <tt>beta</tt> here are
   *           ignored.
   * \param  W_tilde
   *           [in,persisting] The composite state derivative matrix computed
   *           at the base point <tt>stateBasePoint</tt> with
   *           <tt>alpha=coeff_x_dot</tt> and <tt>beta=coeff_x</tt>.  This
   *           argument can be null, in which case it can be computed
   *           internally if needed or not at all.
   * \param  coeff_x_dot
   *           [in] The value of <tt>alpha</tt> for which <tt>W_tilde</tt> was
   *           computed.
   * \param  coeff_x
   *           [in] The value of <tt>beta</tt> for which <tt>W_tilde</tt> was
   *           computed.
   * \param  DfDx_dot
   *           [in] The matrix <tt>d(f)/d(x_dot)</tt> computed at
   *           <tt>stateBasePoint</tt> if available.  If this argument is null,
   *           then it will be computed internally if needed.  The value is
   *           <tt>Teuchos::null</tt>.
   * \param  DfDp
   *           [in] The matrix <tt>d(f)/d(p(p_index))</tt> computed at
   *           <tt>stateBasePoint</tt> if available.  If this argument is null,
   *           then it will be computed internally if needed.  The value is
   *           <tt>Teuchos::null</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>!is_null(W_tilde)</tt>
   * </ul>
   *
   * This function must be called after <tt>intializeStructure()</tt> and
   * before <tt>evalModel()</tt> is called.  After this function is called,
   * <tt>*this</tt> model is fully initialized and ready to compute the
   * requested outputs.
   */
  void initializeState(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &stateBasePoint,
    const Teuchos::RefCountPtr<const Thyra::LinearOpWithSolveBase<Scalar> > &W_tilde,
    const Scalar &coeff_x_dot,
    const Scalar &coeff_x,
    const Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> > &DfDx_dot = Teuchos::null,
    const Teuchos::RefCountPtr<const Thyra::MultiVectorBase<Scalar> > &DfDp = Teuchos::null
    );

  // 2007/05/22: rabartl: ToDo: Add an InterpolationBufferBase
  // stateInterpBuffer object to the initailizeState(...) function that can be
  // used to get x and x_dot at different points in time t.  Then, modify the
  // logic to recompute all of the needed matrices if t != t_base (as passed
  // in through stateBasePoint).  The values of x(t) and xdot(t) can then be
  // gotten from the stateInterpBuffer object!
  
  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> > get_f_space() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief . */
  Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<Scalar> > create_W() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgs() const;
  /** \brief . */
  void evalModel(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
    ) const;

  //@}

private:

  // /////////////////////////
  // Private data members

  Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > stateModel_;
  int p_index_;
  int np_;
  Teuchos::RefCountPtr<const Thyra::DefaultMultiVectorProductVectorSpace<Scalar> > s_bar_space_;
  Teuchos::RefCountPtr<const Thyra::DefaultMultiVectorProductVectorSpace<Scalar> > f_sens_space_;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues_;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> stateBasePoint_;

  mutable Teuchos::RefCountPtr<const Thyra::LinearOpWithSolveBase<Scalar> > W_tilde_;
  mutable Scalar coeff_x_dot_;
  mutable Scalar coeff_x_;
  mutable Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> > DfDx_dot_;
  mutable Teuchos::RefCountPtr<const Thyra::MultiVectorBase<Scalar> > DfDp_;


  // /////////////////////////
  // Private member functions

  void wrapNominalValuesAndBounds();

  void computeDerivativeMatrices(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &point
    ) const;
  
};


// /////////////////////////////////
// Implementations


// Constructors/Intializers/Accessors


template<class Scalar>
ForwardSensitivityModelEvaluator<Scalar>::ForwardSensitivityModelEvaluator()
  : p_index_(0), np_(-1)
{}


template<class Scalar>
void ForwardSensitivityModelEvaluator<Scalar>::initializeStructure(
  const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &stateModel,
  const int p_index
  )
{

  typedef Thyra::ModelEvaluatorBase MEB;

  //
  // Validate input
  //

  TEST_FOR_EXCEPT( is_null(stateModel) );
  TEST_FOR_EXCEPTION(
    !( 0 <= p_index && p_index < stateModel->Np() ), std::logic_error,
    "Error, p_index does not fall in the range [0,"<<(stateModel->Np()-1)<<"]!" );
  // ToDo: Validate support for DfDp!

  //
  // Set the input objects
  //

  stateModel_ = stateModel;
  p_index_ = p_index;
  np_ = stateModel_->get_p_space(p_index)->dim();

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

}


template<class Scalar>
Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> >
ForwardSensitivityModelEvaluator<Scalar>::getStateModel() const
{
  return stateModel_;
}


template<class Scalar>
int ForwardSensitivityModelEvaluator<Scalar>::get_p_index() const
{
  return p_index_;
}


template<class Scalar>
void ForwardSensitivityModelEvaluator<Scalar>::initializeState(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &stateBasePoint,
  const Teuchos::RefCountPtr<const Thyra::LinearOpWithSolveBase<Scalar> > &W_tilde,
  const Scalar &coeff_x_dot,
  const Scalar &coeff_x,
  const Teuchos::RefCountPtr<const Thyra::LinearOpBase<Scalar> > &DfDx_dot,
  const Teuchos::RefCountPtr<const Thyra::MultiVectorBase<Scalar> > &DfDp
  )
{

  typedef Thyra::ModelEvaluatorBase MEB;

#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(
    is_null(stateModel_), std::logic_error,
    "Error, you must call intializeStructure(...) before you call initializeState(...)"
    );
  TEST_FOR_EXCEPT( is_null(stateBasePoint.get_x()) );
  TEST_FOR_EXCEPT( is_null(stateBasePoint.get_x_dot()) );
  TEST_FOR_EXCEPT( is_null(stateBasePoint.get_p(p_index_)) );
  // What about the other parameter values?  We really can't say anything
  // about these and we can't check them.  They can be null just fine.
  if (!is_null(W_tilde)) {
    THYRA_ASSERT_VEC_SPACES("",*W_tilde->domain(),*stateModel_->get_x_space());
    THYRA_ASSERT_VEC_SPACES("",*W_tilde->range(),*stateModel_->get_f_space());
  }
  if (!is_null(DfDx_dot)) {
    THYRA_ASSERT_VEC_SPACES("",*DfDx_dot->domain(),*stateModel_->get_x_space());
    THYRA_ASSERT_VEC_SPACES("",*DfDx_dot->range(),*stateModel_->get_f_space());
  }
  if (!is_null(DfDp)) {
    THYRA_ASSERT_VEC_SPACES("",*DfDp->domain(),*stateModel_->get_p_space(p_index_));
    THYRA_ASSERT_VEC_SPACES("",*DfDp->range(),*stateModel_->get_f_space());
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


// Public functions overridden from ModelEvaulator


template<class Scalar>
Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> >
ForwardSensitivityModelEvaluator<Scalar>::get_x_space() const
{
  return s_bar_space_;
}


template<class Scalar>
Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> >
ForwardSensitivityModelEvaluator<Scalar>::get_f_space() const
{
  return f_sens_space_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
ForwardSensitivityModelEvaluator<Scalar>::getNominalValues() const
{
  return nominalValues_;
}


template<class Scalar>
Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<Scalar> >
ForwardSensitivityModelEvaluator<Scalar>::create_W() const
{
  return Thyra::multiVectorLinearOpWithSolve<Scalar>();
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
ForwardSensitivityModelEvaluator<Scalar>::createInArgs() const
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


template<class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
ForwardSensitivityModelEvaluator<Scalar>::createOutArgs() const
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
void ForwardSensitivityModelEvaluator<Scalar>::evalModel(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  ) const
{

  using Teuchos::Array;
  using Teuchos::RefCountPtr;
  using Teuchos::rcp_dynamic_cast;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::VerboseObjectTempState<Thyra::ModelEvaluatorBase> VOTSME;

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_GEN_BEGIN(
    "ForwardSensitivityModelEvaluator", inArgs, outArgs, Teuchos::null );

  //
  // Update the derivative matrices if they are not already updated for the
  // given time!.
  //

  computeDerivativeMatrices(inArgs);

  //
  // InArgs
  //

  RefCountPtr<const Thyra::DefaultMultiVectorProductVector<Scalar> >
    s_bar = rcp_dynamic_cast<const Thyra::DefaultMultiVectorProductVector<Scalar> >(
      inArgs.get_x().assert_not_null(), true
      );
  RefCountPtr<const Thyra::DefaultMultiVectorProductVector<Scalar> >
    s_bar_dot = rcp_dynamic_cast<const Thyra::DefaultMultiVectorProductVector<Scalar> >(
      inArgs.get_x_dot().assert_not_null(), true
      );
  const Scalar
    alpha = inArgs.get_alpha();
  const Scalar
    beta = inArgs.get_beta();

  RefCountPtr<const Thyra::MultiVectorBase<Scalar> >
    S = s_bar->getMultiVector();
  RefCountPtr<const Thyra::MultiVectorBase<Scalar> >
    S_dot = s_bar_dot->getMultiVector();
  
  //
  // OutArgs
  //

  RefCountPtr<Thyra::DefaultMultiVectorProductVector<Scalar> >
    f_sens = rcp_dynamic_cast<Thyra::DefaultMultiVectorProductVector<Scalar> >(
      outArgs.get_f(), true
      );
  RefCountPtr<Thyra::DefaultMultiVectorLinearOpWithSolve<Scalar> >
    W_sens = rcp_dynamic_cast<Thyra::DefaultMultiVectorLinearOpWithSolve<Scalar> >(
      outArgs.get_W(), true
      );

  RefCountPtr<Thyra::MultiVectorBase<Scalar> >
    F_sens = f_sens->getNonconstMultiVector().assert_not_null();

  //
  // Compute the requested functions
  //

  if(!is_null(F_sens)) {
    // S_diff =  -(coeff_x_dot/coeff_x)*S + S_dot
    Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >
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
    Vp_V( &*F_sens, *DfDp_ );
  }

  if(!is_null(W_sens)) {
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
void ForwardSensitivityModelEvaluator<Scalar>::wrapNominalValuesAndBounds()
{

  using Teuchos::RefCountPtr;
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
void ForwardSensitivityModelEvaluator<Scalar>::computeDerivativeMatrices(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &point
  ) const
{

  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::VerboseObjectTempState<MEB> VOTSME;

  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();

  const Scalar
    t_base = stateBasePoint_.get_t(),
    t = point.get_t();

  TEST_FOR_EXCEPTION(
    t != t_base, std::logic_error,
    "Error, t="<<t<<" is not equal to t_base="<<t_base<<" and we dont' handle this yet!"
    );

  if (is_null(W_tilde_)) {
    TEST_FOR_EXCEPT("ToDo: compute W_tilde from scratch!");
  }
  
  if ( is_null(DfDx_dot_) || is_null(DfDp_) ) {

    MEB::InArgs<Scalar> inArgs = stateBasePoint_;
    MEB::OutArgs<Scalar> outArgs = stateModel_->createOutArgs();
    
    Teuchos::RefCountPtr<Thyra::LinearOpBase<Scalar> > DfDx_dot_compute;
    if (is_null(DfDx_dot_)) {
      DfDx_dot_compute = stateModel_->create_W_op();
      inArgs.set_alpha(1.0);
      inArgs.set_beta(0.0);
      outArgs.set_W_op(DfDx_dot_compute);
    }

    Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > DfDp_compute;
    if (is_null(DfDp_)) {
      DfDp_compute = Thyra::create_DfDp_mv(
        *stateModel_,p_index_,
        MEB::DERIV_MV_BY_COL
        ).getMultiVector();
      outArgs.set_DfDp(
        p_index_,
        MEB::Derivative<Scalar>(DfDp_compute,MEB::DERIV_MV_BY_COL)
        );
    }
    
    VOTSME stateModel_outputTempState(stateModel_,out,verbLevel);
    stateModel_->evalModel(inArgs,outArgs);
    if (!is_null(DfDx_dot_compute))
      DfDx_dot_ = DfDx_dot_compute;
    if (!is_null(DfDp_compute))
      DfDp_ = DfDp_compute;
  
  }

}


} // namespace Rythmos


#endif // RYTHMOS_FORWARD_SENSITIVITY_MODEL_EVALUATOR_HPP
