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


#ifndef RYTHMOS_STEPPER_AS_MODEL_EVALUATOR_HPP
#define RYTHMOS_STEPPER_AS_MODEL_EVALUATOR_HPP


#include "Thyra_ResponseOnlyModelEvaluatorBase.hpp"
#include "Thyra_ModelEvaluatorDelegatorBase.hpp"

#include "Rythmos_StepperBase.hpp"
#include "Rythmos_IntegratorBase.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_Assert.hpp"


namespace Rythmos {


/** \brief Concrete <tt>Thyra::ModelEvaluator</tt> subclass that takes a
 * parameterized stepper and turns it into a model evaluator <tt>(p,t) ==>
 * g</tt>.
 *
 * ToDo: Finish Documentation!
 */
template<class Scalar>
class StepperAsModelEvaluator
  : virtual public Thyra::ResponseOnlyModelEvaluatorBase<Scalar>
{
public:

  /** \name Constructors, Initialization, Misc. */
  //@{

  /** \brief . */
  StepperAsModelEvaluator();

  /** \brief . */
  void initialize(
    const RCP<StepperBase<Scalar> > &stepper,
    const RCP<IntegratorBase<Scalar> > &integrator,
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
    );

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  int Np() const;
  /** \brief . */
  int Ng() const;
  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> >
  get_p_space(int l) const;
  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> >
  get_g_space(int j) const;
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
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs
    ) const;

  //@}

private:

  // //////////////////////
  // Private types

  typedef Array<RCP<const Thyra::VectorSpaceBase<Scalar> > > SpaceArray_t;
  

  // //////////////////////
  // Private data members

  RCP<StepperBase<Scalar> > stepper_;
  RCP<IntegratorBase<Scalar> > integrator_;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> initialCondition_;

  int Np_;
  int Ng_;

  SpaceArray_t g_space_;
  SpaceArray_t p_space_;

  mutable  Thyra::ModelEvaluatorBase::InArgs<Scalar> currentInitialCondition_;

};


/** \brief Nonmember constructor. */
template<class Scalar>
RCP<StepperAsModelEvaluator<Scalar> >
stepperAsModelEvaluator(
  const RCP<StepperBase<Scalar> > &stepper,
  const RCP<IntegratorBase<Scalar> > &integrator,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
  )
{
  using Teuchos::rcp;
  RCP<StepperAsModelEvaluator<Scalar> >
    stepperAsModelEvaluator = rcp(new StepperAsModelEvaluator<Scalar>());
  stepperAsModelEvaluator->initialize(stepper,integrator,initialCondition);
  return stepperAsModelEvaluator;
}


// ////////////////////////
// Definitions


// Constructors, Initialization, Misc.


template<class Scalar>
StepperAsModelEvaluator<Scalar>::StepperAsModelEvaluator()
  : Np_(0), Ng_(0)
{}


template<class Scalar>
void StepperAsModelEvaluator<Scalar>::initialize(
  const RCP<StepperBase<Scalar> > &stepper,
  const RCP<IntegratorBase<Scalar> > &integrator,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
  )
{

#ifdef RYTHMOS_DEBUG
  TEST_FOR_EXCEPT(is_null(stepper));
  TEST_FOR_EXCEPT(is_null(stepper->getModel()));
  TEST_FOR_EXCEPT(is_null(integrator));
#endif
  stepper_ = stepper;
  integrator_ = integrator;
  initialCondition_ = initialCondition;
  currentInitialCondition_ = initialCondition;

  const RCP<const Thyra::ModelEvaluator<Scalar> >
    stepperModel = stepper_->getModel();

  Np_ = stepperModel->Np();
  p_space_.clear();
  for ( int l = 0; l < Np_; ++l ) {
    p_space_.push_back(stepperModel->get_p_space(l));
  }

  Ng_ = 1;
  g_space_.clear();
  g_space_.push_back(stepper_->getModel()->get_x_space());

}


// Public functions overridden from ModelEvaulator


template<class Scalar>
int StepperAsModelEvaluator<Scalar>::Np() const
{
  return Np_;
}


template<class Scalar>
int StepperAsModelEvaluator<Scalar>::Ng() const
{
  return Ng_;
}


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
StepperAsModelEvaluator<Scalar>::get_p_space(int l) const
{
#ifdef RYTHMOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( l, 0, Np_ );
#endif
  return p_space_[l];
}


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
StepperAsModelEvaluator<Scalar>::get_g_space(int j) const
{
#ifdef RYTHMOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( j, 0, Ng_ );
#endif
  return g_space_[j];
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
StepperAsModelEvaluator<Scalar>::createInArgs() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(Np_);
  inArgs.setSupports(MEB::IN_ARG_t);
  return inArgs;
}


// Private functions overridden from ModelEvaulatorDefaultBase


template<class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
StepperAsModelEvaluator<Scalar>::createOutArgsImpl() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(Np_,Ng_);
  return outArgs;
}


template<class Scalar>
void StepperAsModelEvaluator<Scalar>::evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs
  ) const
{

  using Teuchos::as;
  using Teuchos::describe;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Teuchos::VerboseObjectTempState<InterpolationBufferBase<Scalar> > VOTSSB;

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_GEN_BEGIN(
    "Rythmos::StepperAsModelEvaluator", inArgs, outArgs, Teuchos::null
    );
  VOTSSB integrator_outputTempState(integrator_,out,incrVerbLevel(verbLevel,-1));
  //VOTSSB stepper_outputTempState(stepper_,out,incrVerbLevel(verbLevel,-1));

  // InArgs

  const Scalar finalTime = inArgs.get_t();
  for ( int l = 0; l < Np_; ++l ) {
    currentInitialCondition_.set_p(l,inArgs.get_p(l));
  }

  // OutArgs

  RCP<Thyra::VectorBase<Scalar> >
    g_out = outArgs.get_g(0);

  TEST_FOR_EXCEPT_MSG(
    is_null(g_out), "You must ask for g(0) when you call this function!"
    );

#ifdef RYTHMOS_DEBUG

  THYRA_ASSERT_VEC_SPACES(
    "StepperAsModelEvaluator<Scalar>::evalModel(...)",
    *g_out->space(), *stepper_->get_x_space() );

#endif
  
  // Set up the integrator

  stepper_->setInitialCondition(currentInitialCondition_);
  integrator_->setStepper(stepper_,finalTime);

  // Compute the desired response

  if (!is_null(g_out)) {

    // Get x and xdot at the end time
    Array<Scalar> time_vec = Teuchos::tuple<Scalar>(finalTime);
    Array<RCP<const Thyra::VectorBase<Scalar> > > x_vec;
    integrator_->getFwdPoints( time_vec, &x_vec, 0, 0 );
    
    Thyra::V_V( &*g_out, *x_vec[0] );

  }

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_END();

}


} // namespace Rythmos


#endif // RYTHMOS_STEPPER_AS_MODEL_EVALUATOR_HPP
