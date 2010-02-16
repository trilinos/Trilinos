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

#ifndef RYTHMOS_ADJOINT_MODEL_EVALUATOR_HPP
#define RYTHMOS_ADJOINT_MODEL_EVALUATOR_HPP


#include "Rythmos_IntegratorBase.hpp"
#include "Thyra_ModelEvaluator.hpp" // Interface
#include "Thyra_StateFuncModelEvaluatorBase.hpp" // Implementation
#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultAdjointLinearOpWithSolve.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Teuchos_implicit_cast.hpp"
#include "Teuchos_Assert.hpp"


namespace Rythmos {


/** \brief Standard concrete adjoint ModelEvaluator for time-constant mass
 * matrix models.
 *
 * \section Rythmos_AdjointModelEvaluator_Overview_sec Overview
 *
 * This concrete ModelEvalautor subclass takes any suitable ModelEvalautor
 * object and creates the adjoint model for use by any appropriate time
 * integration method..
 *
 * When the mass matrix <tt>d(f)/d(x_dot)</tt> is not a function of <tt>t</tt>
 * (which this class assumes), then the adjoint model can be represented as:
 
 \verbatim

     d(f)/d(x_dot)^T * lambda_dot - d(f)/d(x)^T * lambda + d(g)/d(x)^T

 \endverbatim

 * This model is stated in reverse time <tt>t_bar <:
 * [0,t_final-t_initial]</tt> (with <tt>d/d(t) = -d/d(t_bar)</tt>) which
 * results in the new adjoint formuation
 
 \verbatim

   f_bar(x_bar_dot, x_bar, t_bar)
     = d(f)/d(x_dot)^T * lambda_rev_dot + d(f)/d(x)^T * lambda - d(g)/d(x)^T

 \endverbatim

 * Where:<ul>
 *
 * <li> <tt>t_bar <: [0,t_final-t_initial]</tt> is reverse time defined so
 * that <tt>t = t_final - t_bar<tt> and <tt>d/d(t) = -d/d(t_bar)</tt>
 *
 * <li> <tt>x_bar = lambda</tt> is the original adjoint
 *
 * <li> <tt>x_bar_dot = lambda_rev_dot</tt> is the reverse-time adjoint time
 * derivative where <tt>lambda_dot = -lambda_rev_dot</tt>.
 *
 * <li> <tt>d(f)/d(x_dot)</tt> and <tt>d(f)/d(x)</tt> are evaluated at
 * <tt>x_dot(t_final-t_bar)</tt> and <tt>x(t_final-t_bar)</tt>.
 *
 * The forward state values <tt>x</tt> and <tt>x_dot</tt> are given through an
 * <tt>InterpolationBufferBase</tt> object that is provided by the client.
 *
 * </ul>
 *
 * <b>WARNING!</b> When interacting with this interface you must take note
 * that reverse time is being used as defined above!  This is especially
 * important if you are going to use lambda_dot for anything.  You have been
 * warned!
 *
 * \section Rythmos_AdjointModelEvaluator_ImplementationNotes_sec Implementation Notes
 *
 * Here, we describe how the residual of the adjoint system
 * <tt>f_bar(...)</tt> is actually computed from the capabilities of the
 * underlying forward model.
 *
 * First, note that

 \verbatim

   W_bar = alpha_bar * d(f_bar)/d(x_bar_dot) + beta_bar * d(f_bar)/d(x_bar)

         = alpha_bar * d(f)/d(x_dot)^T + beta_bar * d(f)/d(x)^T

 \endverbatim

 * This means that <tt>W_bar</tt> can be computed directly as
 * <tt>W_bar_adj</tt> on the underlying forward state ModelEvaluator object as:

 \verbatim

   W_bar_adj = alpha_bar * d(f)/d(x_dot) + beta_bar * d(f)/d(x)

 \endverbatim

 * by passing in <tt>alpha = alpha_bar</tt> and <tt>beta = beta_bar</tt>.  We
 * then use the subclass <tt>Thyra::DefaultAdjointLinearOpWithSolve</tt> to
 * create <tt>W_bar = adjoint(W_bar_adj)</tt> and that is it.
 *
 * Now, given that the client will request the form of <tt>W_bar =
 * adjoint(W_bar_adj)</tt> described above, we want to use this
 * <tt>W_bar_adj</tt> object in computing the adjoint equation residual
 * <tt>f_bar</tt>.  To see how to do this, note that from the above definition
 * of <tt>W_bar</tt> that we have:

 \verbatim

   d(f)/d(x)^T = 1/beta_bar * W_bar_adj^T
     - alpha_bar/beta_bar * d(f)/d(x_dot)^T

 \endverbatim

 * By using the above equation for <tt>d(f)/d(x)^T</tt>, we can eliminate
 * <tt>d(f)/d(x)</tt> from <tt>f_bar</tt> to yield:
 
 \verbatim

   f_bar = d(f)/d(x_dot)^T * lambda_hat + 1/beta_bar * W_bar_adj^T * lambda
           - d(g)/d(x)^T

      where:

         lambda_hat = lambda_rev_dot - alpha_bar/beta_bar * lambda

 \endverbatim

 * Above, we see we need to compute <tt>d(f)/d(x_dot)</tt> sperately from
 * <tt>W_bar_adj</tt> from the underlying forward ModelEvaluator.  Note that
 * for many forward models, that <tt>d(f)/d(x_dot)</tt> will actually be
 * constant and can be computed up front and reused throughout.
 *
 * \todo Add support to response function derivative source d(g)/d(x)^T.
 *
 * \todo Add support for more than one adjoint through the
 * DefaultMultiVectorProductVector[Space] sublasses.
 *
 * \todo Add functionality to the Thyra::ModelEvaluator::OutArgs class to
 * communicate the dependence of a function on its input arguments.  We need
 * to know the exact dependance of <tt>f(...)<tt> on <tt>x_dot</tt>,
 * <tt>x</tt>, and <tt>t</tt> to know if this class can be used and what
 * shortcuts can be used with it.
 */
template<class Scalar>
class AdjointModelEvaluator
  : virtual public Thyra::StateFuncModelEvaluatorBase<Scalar>
{
public:

  /** \name Constructors/Intializers/Accessors */
  //@{

  /** \brief . */
  AdjointModelEvaluator();

  /** \brief Set the underlying forward model and base point. */
  void setFwdStateModel(
    const RCP<const Thyra::ModelEvaluator<Scalar> > &fwdStateModel,
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &basePoint );

  /** \brief Set the forward time range that this adjoint model will be
   * defined over.
   */
  void setFwdTimeRange( const TimeRange<Scalar> &fwdTimeRange );

  /** \brief Set the interpolation buffer that will return values of the state
   * solution <tt>x</tt> and <tt>x_dot</tt> at various points <tt>t</tt> as
   * the adjoint is integrated backwards in time.
   *
   * NOTE: If the model is linear in <tt>x</tt> and <tt>x_dot</tt>, then the
   * client can avoid setting this interpolation buffer since it will never be
   * called.
   *
   * NOTE: This object need be fully initialized at this point.  It only needs
   * to be fully initialized before this->evalModel(...) is called.  This just
   * sets up the basic link to this object.
   */
  void setFwdStateSolutionBuffer(
    const RCP<const InterpolationBufferBase<Scalar> > &fwdStateSolutionBuffer );
  
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
  RCP<Thyra::LinearOpBase<Scalar> > create_W_op() const;
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
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs_bar,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs_bar
    ) const;

  //@}

private:

  // /////////////////////////
  // Private data members

  RCP<const Thyra::ModelEvaluator<Scalar> > fwdStateModel_;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> basePoint_;
  TimeRange<Scalar> fwdTimeRange_;
  RCP<const InterpolationBufferBase<Scalar> > fwdStateSolutionBuffer_;

  mutable bool isInitialized_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar> prototypeInArgs_bar_;
  mutable Thyra::ModelEvaluatorBase::OutArgs<Scalar> prototypeOutArgs_bar_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar> adjointNominalValues_;
  mutable RCP<Thyra::LinearOpBase<Scalar> > my_W_bar_adj_op_;
  mutable RCP<Thyra::LinearOpBase<Scalar> > my_d_f_d_x_dot_op_;

  // /////////////////////////
  // Private member functions

  // Just-in-time initialization function
  void initialize() const;

};


/** \brief Nonmember constructor.
 *
 * \relates AdjointModelEvaluator
 */
template<class Scalar>
RCP<AdjointModelEvaluator<Scalar> >
adjointModelEvaluator(
  const RCP<const Thyra::ModelEvaluator<Scalar> > &fwdStateModel,
  const TimeRange<Scalar> &fwdTimeRange
  )
{
  RCP<AdjointModelEvaluator<Scalar> >
    adjointModel = Teuchos::rcp(new AdjointModelEvaluator<Scalar>);
  adjointModel->setFwdStateModel(fwdStateModel, fwdStateModel->getNominalValues());
  adjointModel->setFwdTimeRange(fwdTimeRange);
  return adjointModel;
}


// /////////////////////////////////
// Implementations


// Constructors/Intializers/Accessors


template<class Scalar>
AdjointModelEvaluator<Scalar>::AdjointModelEvaluator()
  :isInitialized_(false)
{}


template<class Scalar>
void AdjointModelEvaluator<Scalar>::setFwdStateModel(
  const RCP<const Thyra::ModelEvaluator<Scalar> > &fwdStateModel,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &basePoint
  )
{
  TEST_FOR_EXCEPT(is_null(fwdStateModel));
  fwdStateModel_ = fwdStateModel;
  basePoint_ = basePoint;
  isInitialized_ = false;
}


template<class Scalar>
void AdjointModelEvaluator<Scalar>::setFwdTimeRange(
  const TimeRange<Scalar> &fwdTimeRange )
{
  fwdTimeRange_ = fwdTimeRange;
}


template<class Scalar>
void AdjointModelEvaluator<Scalar>::setFwdStateSolutionBuffer(
  const RCP<const InterpolationBufferBase<Scalar> > &fwdStateSolutionBuffer )
{
  TEST_FOR_EXCEPT(is_null(fwdStateSolutionBuffer));
  fwdStateSolutionBuffer_ = fwdStateSolutionBuffer;
}


// Public functions overridden from ModelEvaulator


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
AdjointModelEvaluator<Scalar>::get_x_space() const
{
  initialize();
  return fwdStateModel_->get_f_space();
}


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
AdjointModelEvaluator<Scalar>::get_f_space() const
{
  initialize();
  return fwdStateModel_->get_x_space();
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
AdjointModelEvaluator<Scalar>::getNominalValues() const
{
  initialize();
  return adjointNominalValues_;
}


template<class Scalar>
RCP<Thyra::LinearOpWithSolveBase<Scalar> >
AdjointModelEvaluator<Scalar>::create_W() const
{
  initialize();
  return Thyra::nonconstAdjointLows<Scalar>(fwdStateModel_->create_W());
}


template<class Scalar>
RCP<Thyra::LinearOpBase<Scalar> >
AdjointModelEvaluator<Scalar>::create_W_op() const
{
  initialize();
  return Thyra::nonconstAdjoint<Scalar>(fwdStateModel_->create_W_op());
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
AdjointModelEvaluator<Scalar>::createInArgs() const
{
  initialize();
  return prototypeInArgs_bar_;
}


// Private functions overridden from ModelEvaulatorDefaultBase


template<class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
AdjointModelEvaluator<Scalar>::createOutArgsImpl() const
{
  initialize();
  return prototypeOutArgs_bar_;
}


template<class Scalar>
void AdjointModelEvaluator<Scalar>::evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs_bar,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs_bar
  ) const
{

  using Teuchos::rcp_dynamic_cast;
  using Teuchos::describe;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Thyra::DefaultScaledAdjointLinearOp<Scalar> DSALO;
  typedef Thyra::DefaultAdjointLinearOpWithSolve<Scalar> DALOWS;
  typedef Teuchos::VerboseObjectTempState<Thyra::ModelEvaluatorBase> VOTSME;

  //
  // A) Header stuff
  //

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_GEN_BEGIN(
    "AdjointModelEvaluator", inArgs_bar, outArgs_bar, Teuchos::null );

  initialize();

  VOTSME fwdStateModel_outputTempState(fwdStateModel_,out,verbLevel);

  //const bool trace = includesVerbLevel(verbLevel, Teuchos::VERB_LOW);
  const bool dumpAll = includesVerbLevel(localVerbLevel, Teuchos::VERB_EXTREME);

  //
  // B) Unpack the input and output arguments to see what we have to compute
  //

  // B.1) InArgs

  const Scalar t_bar = inArgs_bar.get_t();
  const RCP<const Thyra::VectorBase<Scalar> >
    lambda_rev_dot = inArgs_bar.get_x_dot().assert_not_null(), // x_bar_dot
    lambda = inArgs_bar.get_x().assert_not_null(); // x_bar
  const Scalar alpha_bar = inArgs_bar.get_alpha();
  const Scalar beta_bar = inArgs_bar.get_beta();

  if (dumpAll) {
    *out << "\nlambda_rev_dot = " << describe(*lambda_rev_dot, Teuchos::VERB_EXTREME);
    *out << "\nlambda = " << describe(*lambda, Teuchos::VERB_EXTREME);
    *out << "\nalpha_bar = " << alpha_bar << "\n";
    *out << "\nbeta_bar = " << beta_bar << "\n";
  }

  // B.2) OutArgs

  const RCP<Thyra::VectorBase<Scalar> > f_bar = outArgs_bar.get_f();

  RCP<DALOWS> W_bar;
  if (outArgs_bar.supports(MEB::OUT_ARG_W))
    W_bar = rcp_dynamic_cast<DALOWS>(outArgs_bar.get_W(), true);

  RCP<DSALO> W_bar_op;
  if (outArgs_bar.supports(MEB::OUT_ARG_W_op))
    W_bar_op = rcp_dynamic_cast<DSALO>(outArgs_bar.get_W_op(), true);

  if (dumpAll) {
    if (!is_null(W_bar)) {
      *out << "\nW_bar = " << describe(*W_bar, Teuchos::VERB_EXTREME);
    }
    if (!is_null(W_bar_op)) {
      *out << "\nW_bar_op = " << describe(*W_bar_op, Teuchos::VERB_EXTREME);
    }
  }
  
  //
  // C) Evaluate the needed quantities from the underlying forward Model
  //

  MEB::InArgs<Scalar> fwdInArgs = fwdStateModel_->createInArgs();

  // C.1) Set the required input arguments

  fwdInArgs = basePoint_;

  if (!is_null(fwdStateSolutionBuffer_)) {
    const Scalar t = fwdTimeRange_.length() - t_bar;
    RCP<const Thyra::VectorBase<Scalar> > x, x_dot;
    get_x_and_x_dot<Scalar>( *fwdStateSolutionBuffer_, t,
      outArg(x), outArg(x_dot) );
    fwdInArgs.set_x(x);
    fwdInArgs.set_x_dot(x);
  }
  else {
    // If we don't have an IB object to get the state from, we will assume
    // that the problem is linear and, therefore, we can pass in any old value
    // of x, x_dot, and t and get the W_bar_adj object that we need.  For this
    // purpose, we will assume the model's base point will do.

    // 2008/05/14: rabartl: ToDo: Implement real variable dependancy
    // communication support to make sure that this is okay!  If the model is
    // really nonlinear we need to check for this and throw if the user did
    // not set up a fwdStateSolutionBuffer object!
  }


  // C.2) Evaluate W_bar_adj if needed

  RCP<Thyra::LinearOpWithSolveBase<Scalar> > W_bar_adj;
  RCP<Thyra::LinearOpBase<Scalar> > W_bar_adj_op;
  {

    MEB::OutArgs<Scalar> fwdOutArgs = fwdStateModel_->createOutArgs();
    
    // Get or create W_bar_adj or W_bar_adj_op if needed
    if (!is_null(W_bar)) {
      // If we have W_bar, the W_bar_adj was already created in
      // this->create_W()
      W_bar_adj = W_bar->getNonconstOp();
      W_bar_adj_op = W_bar_adj;
    }
    else if (!is_null(W_bar_op)) {
      // If we have W_bar_op, the W_bar_adj_op was already created in
      // this->create_W_op()
      W_bar_adj_op = W_bar_op->getNonconstOp();
    }
    else if (!is_null(f_bar)) {
      TEST_FOR_EXCEPT_MSG(true, "ToDo: Unit test this code!");
      // If the user did not pass in W_bar or W_bar_op, then we need to create
      // our own local LOB form W_bar_adj_op of W_bar_adj in order to evaluate
      // the residual f_bar
      if (is_null(my_W_bar_adj_op_)) {
        my_W_bar_adj_op_ = fwdStateModel_->create_W_op();
      }
      W_bar_adj_op = my_W_bar_adj_op_;
    }
    
    // Set W_bar_adj or W_bar_adj_op on the OutArgs object
    if (!is_null(W_bar_adj)) {
      fwdOutArgs.set_W(W_bar_adj);
    }
    else if (!is_null(W_bar_adj_op)) {
      fwdOutArgs.set_W_op(W_bar_adj_op);
    }
    
    // Set alpha and beta on OutArgs object
    if (!is_null(W_bar_adj) || !is_null(W_bar_adj_op)) {
      fwdInArgs.set_alpha(alpha_bar);
      fwdInArgs.set_beta(beta_bar);
    }
    
    // Evaluate the model
    if (!is_null(W_bar_adj) || !is_null(W_bar_adj_op)) {
      fwdStateModel_->evalModel( fwdInArgs, fwdOutArgs );
    }
    
    // Print the objects if requested
    if (!is_null(W_bar_adj) && dumpAll)
      *out << "\nW_bar_adj = " << describe(*W_bar_adj, Teuchos::VERB_EXTREME);
    if (!is_null(W_bar_adj_op) && dumpAll)
      *out << "\nW_bar_adj_op = " << describe(*W_bar_adj_op, Teuchos::VERB_EXTREME);

  }
  
  // C.3) Evaluate d(f)/d(x_dot) if needed

  RCP<Thyra::LinearOpBase<Scalar> > d_f_d_x_dot_op;
  if (!is_null(f_bar)) {
    if (is_null(my_d_f_d_x_dot_op_)) {
      my_d_f_d_x_dot_op_ = fwdStateModel_->create_W_op();
    }
    d_f_d_x_dot_op = my_d_f_d_x_dot_op_;
    MEB::OutArgs<Scalar> fwdOutArgs = fwdStateModel_->createOutArgs();
    fwdOutArgs.set_W_op(d_f_d_x_dot_op);
    fwdInArgs.set_alpha(ST::one());
    fwdInArgs.set_beta(ST::zero());
    fwdStateModel_->evalModel( fwdInArgs, fwdOutArgs );
    if (dumpAll) {
      *out << "\nd_f_d_x_dot_op = " << describe(*d_f_d_x_dot_op, Teuchos::VERB_EXTREME);
    }
  }

  //
  // D) Evaluate the adjoint equation residual:
  //
  //   f_bar = d(f)/d(x_dot)^T * lambda_hat + 1/beta_bar * W_bar_adj^T * lambda
  //           - d(g)/d(x)^T
  //

  if (!is_null(f_bar)) {

    // D.1) lambda_hat = lambda_rev_dot - alpha_bar/beta_bar * lambda
    const RCP<Thyra::VectorBase<Scalar> >
      lambda_hat = createMember(lambda_rev_dot->space());
    Thyra::V_VpStV<Scalar>( outArg(*lambda_hat),
      *lambda_rev_dot, -alpha_bar/beta_bar, *lambda );
    if (dumpAll)
      *out << "\nlambda_hat = " << describe(*lambda_hat, Teuchos::VERB_EXTREME);

    // D.2) f_bar = d(f)/d(x_dot)^T * lambda_hat
    Thyra::apply<Scalar>( *d_f_d_x_dot_op, Thyra::CONJTRANS, *lambda_hat,
      outArg(*f_bar) );

    // D.3) f_bar += 1/beta_bar * W_bar_adj^T * lambda
    Thyra::apply<Scalar>( *W_bar_adj_op, Thyra::CONJTRANS, *lambda,
      outArg(*f_bar), 1.0/beta_bar, ST::one() );

    // D.4) f_bar += - d(g)/d(x)^T
    // 2008/05/15: rabart: ToDo: Implement once we add support for
    // distributed response functions

    if (dumpAll)
      *out << "\nf_bar = " << describe(*f_bar, Teuchos::VERB_EXTREME);

  }

  if (dumpAll) {
    if (!is_null(W_bar)) {
      *out << "\nW_bar = " << describe(*W_bar, Teuchos::VERB_EXTREME);
    }
    if (!is_null(W_bar_op)) {
      *out << "\nW_bar_op = " << describe(*W_bar_op, Teuchos::VERB_EXTREME);
    }
  }


  //
  // E) Do any remaining post processing
  //

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_END();

}


// private


template<class Scalar>
void AdjointModelEvaluator<Scalar>::initialize() const
{

  typedef Thyra::ModelEvaluatorBase MEB;

  if (isInitialized_)
    return;

  //
  // A) Validate the that forward Model is of the correct form!
  //

  MEB::InArgs<Scalar> fwdStateModelInArgs = fwdStateModel_->createInArgs();
  MEB::OutArgs<Scalar> fwdStateModelOutArgs = fwdStateModel_->createOutArgs();

#ifdef RYTHMOS_DEBUG
  TEUCHOS_ASSERT( fwdStateModelInArgs.supports(MEB::IN_ARG_x_dot) );
  TEUCHOS_ASSERT( fwdStateModelInArgs.supports(MEB::IN_ARG_x) );
  TEUCHOS_ASSERT( fwdStateModelInArgs.supports(MEB::IN_ARG_t) );
  TEUCHOS_ASSERT( fwdStateModelInArgs.supports(MEB::IN_ARG_alpha) );
  TEUCHOS_ASSERT( fwdStateModelInArgs.supports(MEB::IN_ARG_beta) );
  TEUCHOS_ASSERT( fwdStateModelOutArgs.supports(MEB::OUT_ARG_f) );
  TEUCHOS_ASSERT( fwdStateModelOutArgs.supports(MEB::OUT_ARG_W) );
#endif

  //
  // B) Set up the prototypical InArgs and OutArgs
  //

  {
    MEB::InArgsSetup<Scalar> inArgs_bar;
    inArgs_bar.setModelEvalDescription(this->description());
    inArgs_bar.setSupports( MEB::IN_ARG_x_dot );
    inArgs_bar.setSupports( MEB::IN_ARG_x );
    inArgs_bar.setSupports( MEB::IN_ARG_t );
    inArgs_bar.setSupports( MEB::IN_ARG_alpha );
    inArgs_bar.setSupports( MEB::IN_ARG_beta );
    prototypeInArgs_bar_ = inArgs_bar;
  }

  {
    MEB::OutArgsSetup<Scalar> outArgs_bar;
    outArgs_bar.setModelEvalDescription(this->description());
    outArgs_bar.setSupports(MEB::OUT_ARG_f);
    if (fwdStateModelOutArgs.supports(MEB::OUT_ARG_W) ) {
      outArgs_bar.setSupports(MEB::OUT_ARG_W);
      outArgs_bar.set_W_properties(fwdStateModelOutArgs.get_W_properties());
    }
    if (fwdStateModelOutArgs.supports(MEB::OUT_ARG_W_op) ) {
      outArgs_bar.setSupports(MEB::OUT_ARG_W_op);
      outArgs_bar.set_W_properties(fwdStateModelOutArgs.get_W_properties());
    }
    prototypeOutArgs_bar_ = outArgs_bar;
  }

  //
  // D) Set up the nominal values for the adjoint
  //

  // Copy structure
  adjointNominalValues_ = prototypeInArgs_bar_;
  // Just set a zero initial condition for the adjoint
  const RCP<Thyra::VectorBase<Scalar> > zero_lambda_vec =
    createMember(fwdStateModel_->get_f_space());
  V_S( zero_lambda_vec.ptr(), ScalarTraits<Scalar>::zero() );
  adjointNominalValues_.set_x_dot(zero_lambda_vec);
  adjointNominalValues_.set_x(zero_lambda_vec);

  //
  // E) Wipe out other cached objects
  //

  my_W_bar_adj_op_ = Teuchos::null;
  my_d_f_d_x_dot_op_ = Teuchos::null;

  //
  // F) We are initialized!
  //

  isInitialized_ = true;

}


} // namespace Rythmos


#endif // RYTHMOS_ADJOINT_MODEL_EVALUATOR_HPP
