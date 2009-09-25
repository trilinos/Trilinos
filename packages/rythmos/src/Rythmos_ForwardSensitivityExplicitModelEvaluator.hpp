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

#ifndef RYTHMOS_FORWARD_SENSITIVITY_EXPLICIT_MODEL_EVALUATOR_HPP
#define RYTHMOS_FORWARD_SENSITIVITY_EXPLICIT_MODEL_EVALUATOR_HPP


#include "Rythmos_IntegratorBase.hpp"
#include "Rythmos_ForwardSensitivityModelEvaluatorBase.hpp"
#include "Thyra_ModelEvaluator.hpp" // Interface
#include "Thyra_StateFuncModelEvaluatorBase.hpp" // Implementation
#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"
#include "Thyra_DefaultMultiVectorLinearOpWithSolve.hpp"


namespace Rythmos {


/** \brief Explicit forward sensitivity transient <tt>ModelEvaluator</tt>
 * subclass.
 *
 * This class provides a very general implemenation of a linear forward
 * sensitivity model evaluator for an explicit ODE.
 *
 * \section Rythmos_ForwardSensitivityExplicitModelEvaluator_intro_sec Introduction
 *
 * The form of the parameterized state equation is:

 \verbatim

   x_dot(t) = f(x(t),{p_l},t), over t = [t0,tf] 

   x(t0) = x_init(p)

 \endverbatim

 * As shown above, the parameters are assumed to be steady state and can enter
 * through the intial condition and/or through the ODE equation itself.
 *
 * The explicit forward sensitivity equations written in multi-vector form
 * are:

 \verbatim

   S_dot = d(f)/d(x)*S + d(f)/d(p) = 0, over t = [t0,tf]

   S(t0) = d(x_init)/d(p)

 \endverbatim

 * where <tt>S</tt> is a multi-vector with <tt>np</tt> columns where each
 * column <tt>S(:,j) = d(x)/d(p_j)</tt> is the sensitivity of <tt>x(t)</tt>
 * with respect to the <tt>p_j</tt> parameter.  The sensitivity parameter
 * subvector <tt>p</tt> here is really just one of the parameter subvectors in
 * the state equation.  This index of the parameter subvector for which the
 * sensitivity equations are written for is given by <tt>p_index</tt>.  Note
 * that above <tt>d(f)/d(x)</tt> and <tt>d(f)/d(p</tt>
 * are both evaluated at the solution to the state equation
 * <tt>(x(t),t)</tt> and are not functions of <tt>S_dot</tt> or
 * <tt>S</tt>.
 *
 * Since the model evaluator interface must be expressed in vector form, the
 * multi-vector form of the forward sensitivity equations is flattened out
 * into:

 \verbatim

   s_bar_dot(t) = f_sens(s_bar(t),{p_l},t), over t = [t0,tf]

   s_bar(t0) = s_bar_init

 \endverbatim

 * where

 \verbatim

   s_bar = [ S(:,0); S(:,0); ...; S(:,np-1) ]

            [ d(f)/d(x)*S(:,0) + d(f)/d(p(0))       ]
            [ d(f)/d(x)*S(:,1) + d(f)/d(p(1))       ]
   f_sens = [ ...                                   ]
            [ d(f)/d(x)*S(:,np-1) + d(f)/d(p(np-1)) ]

   s_bar_init = [ d(x_init)/d(p(0)); d(x_init)/d(p(1)); ...; d(x_init)/d(p(np-1)) ]

 \endverbatim

 * The product vector <tt>s_bar</tt> is represented as a specialized
 * <tt>Thyra::ProductVectorBase</tt> subclass object with <tt>np</tt> "blocks"
 * in terms of a single <tt>Thyra::MultiVectorBase</tt> object (which has
 * <tt>np</tt> columns).
 *
 * ToDo: Finish documention!
 */
template<class Scalar>
class ForwardSensitivityExplicitModelEvaluator
  : virtual public ForwardSensitivityModelEvaluatorBase<Scalar>
{
public:

  /** \name Constructors/Intializers/Accessors */
  //@{

  /** \brief . */
  ForwardSensitivityExplicitModelEvaluator();

  //@}

  /** \name Public functions overridden from ForwardSensitivityModelEvaluatorBase. */
  //@{

  /** \brief . */
  void initializeStructure(
    const RCP<const Thyra::ModelEvaluator<Scalar> > &stateModel,
    const int p_index
    );
  
  /** \brief . */
  void initializeStructureInitCondOnly(
    const RCP<const Thyra::ModelEvaluator<Scalar> >& stateModel,
    const RCP<const Thyra::VectorSpaceBase<Scalar> >& p_space
    );

  /** \brief . */
  RCP<const Thyra::ModelEvaluator<Scalar> >
  getStateModel() const;

  /** \brief . */
  RCP<Thyra::ModelEvaluator<Scalar> >
  getNonconstStateModel() const;
  
  /** \brief . */
  int get_p_index() const;
  
  /** \brief . */
  RCP<const Thyra::DefaultMultiVectorProductVectorSpace<Scalar> >
  get_s_bar_space() const;
  
  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_sens_space() const;

  /** \brief Initialize full state for a single point in time.
   */
  void initializePointState(
      Ptr<StepperBase<Scalar> > stateStepper,
      bool forceUpToDateW
      );

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
  int np_;

  RCP<const Thyra::DefaultMultiVectorProductVectorSpace<Scalar> > s_bar_space_;
  RCP<const Thyra::DefaultMultiVectorProductVectorSpace<Scalar> > f_sens_space_;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues_;

  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar> stateBasePoint_;

  mutable RCP<Thyra::LinearOpBase<Scalar> > DfDx_;
  mutable RCP<Thyra::MultiVectorBase<Scalar> > DfDp_;

  mutable RCP<Thyra::LinearOpBase<Scalar> > DfDx_compute_;
  mutable RCP<Thyra::MultiVectorBase<Scalar> > DfDp_compute_;

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
RCP<ForwardSensitivityExplicitModelEvaluator<Scalar> >
forwardSensitivityExplicitModelEvaluator()
{
  RCP<ForwardSensitivityExplicitModelEvaluator<Scalar> > fseme =
    rcp(new ForwardSensitivityExplicitModelEvaluator<Scalar> );
  return fseme;
}


template<class Scalar>
ForwardSensitivityExplicitModelEvaluator<Scalar>::ForwardSensitivityExplicitModelEvaluator()
  : p_index_(0), np_(-1)
{}



// Public functions overridden from ForwardSensitivityModelEvaluatorBase


template<class Scalar>
void ForwardSensitivityExplicitModelEvaluator<Scalar>::initializeStructure(
  const RCP<const Thyra::ModelEvaluator<Scalar> > &stateModel,
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

  this->wrapNominalValuesAndBounds();

  //
  // Wipe out matrix storage
  //

  stateBasePoint_ = MEB::InArgs<Scalar>();
  DfDx_ = Teuchos::null;
  DfDp_ = Teuchos::null;
  DfDx_compute_ = Teuchos::null;
  DfDp_compute_ = Teuchos::null;

}


template<class Scalar>
void ForwardSensitivityExplicitModelEvaluator<Scalar>::initializeStructureInitCondOnly(
  const RCP<const Thyra::ModelEvaluator<Scalar> >& stateModel,
  const RCP<const Thyra::VectorSpaceBase<Scalar> >& p_space
  )
{
  TEST_FOR_EXCEPT_MSG(true, "ToDo: Implement initializeStructureInitCondOnly()!" );
}


template<class Scalar>
RCP<const Thyra::ModelEvaluator<Scalar> >
ForwardSensitivityExplicitModelEvaluator<Scalar>::getStateModel() const
{
  return stateModel_;
}


template<class Scalar>
RCP<Thyra::ModelEvaluator<Scalar> >
ForwardSensitivityExplicitModelEvaluator<Scalar>::getNonconstStateModel() const
{
  return Teuchos::null;
}


template<class Scalar>
int ForwardSensitivityExplicitModelEvaluator<Scalar>::get_p_index() const
{
  return p_index_;
}


template<class Scalar>
RCP<const Thyra::DefaultMultiVectorProductVectorSpace<Scalar> >
ForwardSensitivityExplicitModelEvaluator<Scalar>::get_s_bar_space() const
{
  return s_bar_space_;
}


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
ForwardSensitivityExplicitModelEvaluator<Scalar>::get_p_sens_space() const
{
  return stateModel_->get_p_space(p_index_);
}


template<class Scalar>
void ForwardSensitivityExplicitModelEvaluator<Scalar>::initializePointState(
    Ptr<StepperBase<Scalar> > stateStepper,
    bool forceUpToDateW
    )
{
  TEUCHOS_ASSERT( Teuchos::nonnull(stateStepper) );
#ifdef RYTHMOS_DEBUG
  TEST_FOR_EXCEPTION(
    is_null(stateModel_), std::logic_error,
    "Error, you must call intializeStructure(...) before you call initializePointState(...)"
    );
#endif // RYTHMOS_DEBUG

  Scalar curr_t = stateStepper->getStepStatus().time;
  RCP<const Thyra::VectorBase<Scalar> > x;
  x = get_x(*stateStepper,curr_t);
#ifdef RYTHMOS_DEBUG
  TEST_FOR_EXCEPT( Teuchos::is_null(x) );
#endif // RYTHMOS_DEBUG
      
  stateBasePoint_ = stateStepper->getInitialCondition(); // set parameters
  stateBasePoint_.set_x( x );
  stateBasePoint_.set_t( curr_t );

  // Set whatever derivatives where passed in.  If an input in null, then the
  // member will be null and the null linear operators will be computed later
  // just in time.

  wrapNominalValuesAndBounds();
  
}


// Public functions overridden from ModelEvaulator


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
ForwardSensitivityExplicitModelEvaluator<Scalar>::get_x_space() const
{
  return s_bar_space_;
}


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
ForwardSensitivityExplicitModelEvaluator<Scalar>::get_f_space() const
{
  return f_sens_space_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
ForwardSensitivityExplicitModelEvaluator<Scalar>::getNominalValues() const
{
  return nominalValues_;
}


template<class Scalar>
RCP<Thyra::LinearOpWithSolveBase<Scalar> >
ForwardSensitivityExplicitModelEvaluator<Scalar>::create_W() const
{
  return Thyra::multiVectorLinearOpWithSolve<Scalar>();
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
ForwardSensitivityExplicitModelEvaluator<Scalar>::createInArgs() const
{
  TEUCHOS_ASSERT( !is_null(stateModel_) );
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgs<Scalar> stateModelInArgs = stateModel_->createInArgs();
  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports( MEB::IN_ARG_x );
  inArgs.setSupports( MEB::IN_ARG_t );
  inArgs.setSupports( MEB::IN_ARG_beta,
    stateModelInArgs.supports(MEB::IN_ARG_beta) );
  return inArgs;
}


// Private functions overridden from ModelEvaulatorDefaultBase


template<class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
ForwardSensitivityExplicitModelEvaluator<Scalar>::createOutArgsImpl() const
{
  TEUCHOS_ASSERT( !is_null(stateModel_) );
  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::OutArgs<Scalar> stateModelOutArgs = stateModel_->createOutArgs();
  MEB::OutArgsSetup<Scalar> outArgs;

  outArgs.setModelEvalDescription(this->description());

  outArgs.setSupports(MEB::OUT_ARG_f);

  return outArgs;

}


template<class Scalar>
void ForwardSensitivityExplicitModelEvaluator<Scalar>::evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  ) const
{

  using Teuchos::rcp_dynamic_cast;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::VerboseObjectTempState<Thyra::ModelEvaluatorBase> VOTSME;

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_GEN_BEGIN(
    "ForwardSensitivityExplicitModelEvaluator", inArgs, outArgs, Teuchos::null );

  //
  // Update the derivative matrices if they are not already updated for the
  // given time!.
  //
  
  {
#ifdef ENABLE_RYTHMOS_TIMERS
    TEUCHOS_FUNC_TIME_MONITOR("Rythmos:ForwardSensitivityExplicitModelEvaluator::evalModel: computeMatrices");
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
  RCP<const Thyra::MultiVectorBase<Scalar> >
    S = s_bar->getMultiVector();
  
  //
  // OutArgs
  //

  RCP<Thyra::DefaultMultiVectorProductVector<Scalar> >
    f_sens = rcp_dynamic_cast<Thyra::DefaultMultiVectorProductVector<Scalar> >(
      outArgs.get_f(), true
      );

  RCP<Thyra::MultiVectorBase<Scalar> >
    F_sens = f_sens->getNonconstMultiVector().assert_not_null();

  //
  // Compute the requested functions
  //

  if(!is_null(F_sens)) {

#ifdef ENABLE_RYTHMOS_TIMERS
    TEUCHOS_FUNC_TIME_MONITOR("Rythmos:ForwardSensitivityExplicitModelEvaluator::evalModel: computeSens");
#endif
    
    // Form the residual:  df/dx * S + df/dp
    // F_sens = df/dx * S
    Thyra::apply(
      *DfDx_, Thyra::NOTRANS,
      *S, &*F_sens,
      ST::one(), ST::zero()
      );
    // F_sens += d(f)/d(p)
    Vp_V( &*F_sens, *DfDp_ );
  }
  
  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_END();

}


// private


template<class Scalar>
void ForwardSensitivityExplicitModelEvaluator<Scalar>::wrapNominalValuesAndBounds()
{
  TEUCHOS_ASSERT( !is_null(stateModel_) );

  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  // nominalValues_.clear(); // ToDo: Implement this!

  nominalValues_.set_t(stateModel_->getNominalValues().get_t());
  
  // 2007/05/22: rabartl: Note: Currently there is not much of a reason to set
  // an initial condition here since the initial condition for the
  // sensitivities is really being set in the ForwardSensitivityStepper
  // object!  In the future, a more general use of this class might benefit
  // from setting the initial condition here.

  // 2009/07/20: tscoffe/ccober:  This is the future.  We're going to use this
  // in a more general way, so we need legitimate nominal values now.
  RCP<VectorBase<Scalar> > s_bar_ic = Thyra::createMember(this->get_x_space());
  Thyra::V_S(&*s_bar_ic,ST::zero());
  nominalValues_.set_x(s_bar_ic);
}


template<class Scalar>
void ForwardSensitivityExplicitModelEvaluator<Scalar>::computeDerivativeMatrices(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &point
  ) const
{
  TEUCHOS_ASSERT( !is_null(stateModel_) );

  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::VerboseObjectTempState<MEB> VOTSME;

  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();

  MEB::InArgs<Scalar> inArgs = stateBasePoint_;
  MEB::OutArgs<Scalar> outArgs = stateModel_->createOutArgs();
  
  if (is_null(DfDx_)) {
    DfDx_ = stateModel_->create_W_op();
  }
  if (inArgs.supports(MEB::IN_ARG_beta)) {
    inArgs.set_beta(1.0);
  }
  outArgs.set_W_op(DfDx_);

  if (is_null(DfDp_)) {
    DfDp_ = Thyra::create_DfDp_mv(
      *stateModel_,p_index_,
      MEB::DERIV_MV_BY_COL
      ).getMultiVector();
  }
  outArgs.set_DfDp(
    p_index_,
    MEB::Derivative<Scalar>(DfDp_,MEB::DERIV_MV_BY_COL)
    );
  
  VOTSME stateModel_outputTempState(stateModel_,out,verbLevel);
  stateModel_->evalModel(inArgs,outArgs);
  

}


} // namespace Rythmos


#endif // RYTHMOS_FORWARD_SENSITIVITY_EXPLICIT_MODEL_EVALUATOR_HPP
