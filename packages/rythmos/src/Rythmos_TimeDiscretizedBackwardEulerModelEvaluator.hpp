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


#ifndef RYTHMOS_TIME_DISCRETIZED_BACKWARD_EULER_MODEL_EVALUATOR_HPP
#define RYTHMOS_TIME_DISCRETIZED_BACKWARD_EULER_MODEL_EVALUATOR_HPP


#include "Rythmos_Types.hpp"
#include "Thyra_StateFuncModelEvaluatorBase.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultBlockedTriangularLinearOpWithSolveFactory.hpp" // Default implementation


namespace Rythmos {


/** \brief .
 *
 * ToDo: Finish Documentation!
 */
template<class Scalar>
class TimeDiscretizedBackwardEulerModelEvaluator
  : virtual public Thyra::StateFuncModelEvaluatorBase<Scalar>
{
public:

  /** \name Constructors/initializers/accessors */
  //@{

  /** \brief . */
  TimeDiscretizedBackwardEulerModelEvaluator();

  //@}

  /** \name Overridden from TimeDiscretizedBackwardEulerModelEvaluatorBase */
  //@{

  /** \brief . */
  void initialize(
    const RCP<const Thyra::ModelEvaluator<Scalar> > &daeModel,
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initCond,
    const Scalar finalTime,
    const int numTimeSteps,
    const RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > &W_bar_factory = Teuchos::null
    );

  //@}

  /** \name Public functions overridden from ModelEvaluator */
  //@{

  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;
  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_space() const;
  /** \brief . */
  RCP<Thyra::LinearOpBase<Scalar> > create_W_op() const;
  /** \breif . */
  RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > get_W_factory() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

  //@}

private:

  /** \name Private functions overridden from ModelEvaluatorDefaultBase */
  //@{

  /** \brief . */
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  /** \brief . */
  void evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs
    ) const;

  //@}

private:

  RCP<const Thyra::ModelEvaluator<Scalar> > daeModel_;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> initCond_;
  Scalar finalTime_;
  int numTimeSteps_;

  Scalar initTime_;
  Scalar delta_t_;

  RCP<const Thyra::ProductVectorSpaceBase<Scalar> > x_bar_space_;
  RCP<const Thyra::ProductVectorSpaceBase<Scalar> > f_bar_space_;
  RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > W_bar_factory_;

};


/** \brief Non-member constructor.
 *
 * \relates TimeDiscretizedBackwardEulerModelEvaluator.
 */
template<class Scalar>
RCP<TimeDiscretizedBackwardEulerModelEvaluator<Scalar> >
timeDiscretizedBackwardEulerModelEvaluator(
  const RCP<const Thyra::ModelEvaluator<Scalar> > &daeModel,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initCond,
  const Scalar finalTime,
  const int numTimeSteps,
  const RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > &W_bar_factory
  )
{
  RCP<TimeDiscretizedBackwardEulerModelEvaluator<Scalar> >
    model(new TimeDiscretizedBackwardEulerModelEvaluator<Scalar>());
  model->initialize(daeModel,initCond,finalTime,numTimeSteps,W_bar_factory);
  return model;
}


// ///////////////////////
// Definition


// Constructors/initializers/accessors


template<class Scalar>
TimeDiscretizedBackwardEulerModelEvaluator<Scalar>::TimeDiscretizedBackwardEulerModelEvaluator()
  :finalTime_(-1.0),
   numTimeSteps_(-1),
   initTime_(0.0),
   delta_t_(-1.0) // Flag for uninitailized!
{}


// Overridden from TimeDiscretizedBackwardEulerModelEvaluatorBase


template<class Scalar>
void TimeDiscretizedBackwardEulerModelEvaluator<Scalar>::initialize(
  const RCP<const Thyra::ModelEvaluator<Scalar> > &daeModel,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initCond,
  const Scalar finalTime,
  const int numTimeSteps,
  const RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > &W_bar_factory
  )
{

  TEST_FOR_EXCEPT(is_null(daeModel));
  TEST_FOR_EXCEPT(is_null(initCond.get_x()));
  TEST_FOR_EXCEPT(is_null(initCond.get_x_dot()));
  TEST_FOR_EXCEPT(finalTime <= initCond.get_t());
  TEST_FOR_EXCEPT(numTimeSteps <= 0);
  // ToDo: Validate that daeModel is of the right form!

  daeModel_ = daeModel;
  initCond_ = initCond;
  finalTime_ = finalTime;
  numTimeSteps_ = numTimeSteps;

  initTime_ = initCond.get_t();
  delta_t_ = (finalTime_ - initTime_) / numTimeSteps_;

  x_bar_space_ = productVectorSpace(daeModel_->get_x_space(),numTimeSteps_);
  f_bar_space_ = productVectorSpace(daeModel_->get_f_space(),numTimeSteps_);

  if (!is_null(W_bar_factory)) {
    W_bar_factory_ = W_bar_factory;
  }
  else {
    W_bar_factory_ =
      Thyra::defaultBlockedTriangularLinearOpWithSolveFactory<Scalar>(
        daeModel_->get_W_factory()
        );
  }
  
}


// Public functions overridden from ModelEvaluator


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
TimeDiscretizedBackwardEulerModelEvaluator<Scalar>::get_x_space() const
{
  return x_bar_space_;
}


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
TimeDiscretizedBackwardEulerModelEvaluator<Scalar>::get_f_space() const
{
  return f_bar_space_;
}


template<class Scalar>
RCP<Thyra::LinearOpBase<Scalar> >
TimeDiscretizedBackwardEulerModelEvaluator<Scalar>::create_W_op() const
{
  // Create the block structure for W_op_bar right away!
  RCP<Thyra::PhysicallyBlockedLinearOpBase<Scalar> >
    W_op_bar = Thyra::defaultBlockedLinearOp<Scalar>();
  W_op_bar->beginBlockFill( f_bar_space_, x_bar_space_ );
  for ( int k = 0; k < numTimeSteps_; ++k ) {
    W_op_bar->setNonconstBlock( k, k, daeModel_->create_W_op() );
    if (k > 0)
      W_op_bar->setNonconstBlock( k, k-1, daeModel_->create_W_op() );
  }
  W_op_bar->endBlockFill();
  return W_op_bar;
}


template<class Scalar>
RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
TimeDiscretizedBackwardEulerModelEvaluator<Scalar>::get_W_factory() const
{
  return W_bar_factory_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
TimeDiscretizedBackwardEulerModelEvaluator<Scalar>::getNominalValues() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  TEST_FOR_EXCEPT(true);
  return MEB::InArgs<Scalar>();
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
TimeDiscretizedBackwardEulerModelEvaluator<Scalar>::createInArgs() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(MEB::IN_ARG_x);
  return inArgs;
}


// Private functions overridden from ModelEvaluatorDefaultBase


template<class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
TimeDiscretizedBackwardEulerModelEvaluator<Scalar>::createOutArgsImpl() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::OutArgs<Scalar> daeOutArgs = daeModel_->createOutArgs();
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.setSupports(MEB::OUT_ARG_f);
  outArgs.setSupports(MEB::OUT_ARG_W_op);
  outArgs.set_W_properties(daeOutArgs.get_W_properties());
  return outArgs;
}


template<class Scalar>
void TimeDiscretizedBackwardEulerModelEvaluator<Scalar>::evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs_bar,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs_bar
  ) const
{


  using Teuchos::rcp_dynamic_cast;
  typedef ScalarTraits<Scalar> ST;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Thyra::VectorBase<Scalar> VB;
  typedef Thyra::ProductVectorBase<Scalar> PVB;
  typedef Thyra::BlockedLinearOpBase<Scalar> BLWB;

/*
  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_GEN_BEGIN(
    "Rythmos::ImplicitRKModelEvaluator",inArgs_bar,outArgs_bar,daeModel_
    );
*/

  TEST_FOR_EXCEPTION( delta_t_ <= 0.0, std::logic_error,
    "Error, you have not initialized this object correctly!" );

  //
  // A) Unwrap the inArgs and outArgs to get at product vectors and block op
  //

  const RCP<const PVB> x_bar = rcp_dynamic_cast<const PVB>(inArgs_bar.get_x(), true);
  const RCP<PVB> f_bar = rcp_dynamic_cast<PVB>(outArgs_bar.get_f(), true);
  RCP<BLWB> W_op_bar = rcp_dynamic_cast<BLWB>(outArgs_bar.get_W_op(), true);

  //
  // B) Assemble f_bar and W_op_bar by looping over stages
  //

  MEB::InArgs<Scalar> daeInArgs = daeModel_->createInArgs();
  MEB::OutArgs<Scalar> daeOutArgs = daeModel_->createOutArgs();
  const RCP<VB> x_dot_i = createMember(daeModel_->get_x_space());
  daeInArgs.setArgs(initCond_);
  
  Scalar t_i = initTime_; // ToDo: Define t_init!

  const Scalar oneOverDeltaT = 1.0/delta_t_;

  for ( int i = 0; i < numTimeSteps_; ++i ) {

    // B.1) Setup the DAE's inArgs for time step eqn f(i) ...
    const RCP<const Thyra::VectorBase<Scalar> >
      x_i = x_bar->getVectorBlock(i),
      x_im1 = ( i==0 ? initCond_.get_x() : x_bar->getVectorBlock(i-1) );
    V_VmV( x_dot_i.ptr(), *x_i, *x_im1 ); // x_dot_i = 1/dt * ( x[i] - x[i-1] )
    Vt_S( x_dot_i.ptr(), oneOverDeltaT ); // ... 
    daeInArgs.set_x_dot( x_dot_i );
    daeInArgs.set_x( x_i );
    daeInArgs.set_t( t_i );
    daeInArgs.set_alpha( oneOverDeltaT );
    daeInArgs.set_beta( 1.0 );

    // B.2) Setup the DAE's outArgs for f(i) and/or W(i,i) ...
    if (!is_null(f_bar))
      daeOutArgs.set_f( f_bar->getNonconstVectorBlock(i) );
    if (!is_null(W_op_bar))
      daeOutArgs.set_W_op(W_op_bar->getNonconstBlock(i,i).assert_not_null());

    // B.3) Compute f_bar(i) and/or W_op_bar(i,i) ...
    daeModel_->evalModel( daeInArgs, daeOutArgs );
    daeOutArgs.set_f(Teuchos::null);
    daeOutArgs.set_W_op(Teuchos::null);
    
    // B.4) Evaluate W_op_bar(i,i-1)
    if ( !is_null(W_op_bar) && i > 0 ) {
      daeInArgs.set_alpha( -oneOverDeltaT );
      daeInArgs.set_beta( 0.0 );
      daeOutArgs.set_W_op(W_op_bar->getNonconstBlock(i,i-1).assert_not_null());
      daeModel_->evalModel( daeInArgs, daeOutArgs );
      daeOutArgs.set_W_op(Teuchos::null);
    }

    //
    t_i += delta_t_;

  }

/*  
  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_END();
*/

}


} // namespace Rythmos


#endif // RYTHMOS_TIME_DISCRETIZED_BACKWARD_EULER_MODEL_EVALUATOR_HPP
