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

#include "Rythmos_StepperBase.hpp"
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
    const Teuchos::RefCountPtr<StepperBase<Scalar> > &stepper,
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
    );

  // ToDo: Replace this with a real time integrator!
  STANDARD_MEMBER_COMPOSITION_MEMBERS(int,numTimeSteps);

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  int Np() const;
  /** \brief . */
  int Ng() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> >
  get_p_space(int l) const;
  /** \brief . */
  Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> >
  get_g_space(int j) const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgs() const;
  /** \brief . */
  void evalModel(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs
    ) const;

  //@}

private:

  // //////////////////////
  // Private types

  typedef Teuchos::Array<Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> > > SpaceArray_t;
  

  // //////////////////////
  // Private data members

  Teuchos::RefCountPtr<StepperBase<Scalar> > stepper_;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> initialCondition_;

  int Np_;
  int Ng_;

  SpaceArray_t g_space_;
  SpaceArray_t p_space_;

  mutable  Thyra::ModelEvaluatorBase::InArgs<Scalar> currentInitialCondition_;

};


/** \brief Nonmember constructor. */
template<class Scalar>
Teuchos::RefCountPtr<StepperAsModelEvaluator<Scalar> >
stepperAsModelEvaluator(
  const Teuchos::RefCountPtr<StepperBase<Scalar> > &stepper,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
  )
{
  using Teuchos::RefCountPtr; using Teuchos::rcp;
  RefCountPtr<StepperAsModelEvaluator<Scalar> >
    stepperAsModelEvaluator = rcp(new StepperAsModelEvaluator<Scalar>());
  stepperAsModelEvaluator->initialize(stepper,initialCondition);
  return stepperAsModelEvaluator;
}


// ////////////////////////
// Definitions


// Constructors, Initialization, Misc.


template<class Scalar>
StepperAsModelEvaluator<Scalar>::StepperAsModelEvaluator()
  :numTimeSteps_(-1),
   Np_(0),
   Ng_(0)
{}


template<class Scalar>
void StepperAsModelEvaluator<Scalar>::initialize(
  const Teuchos::RefCountPtr<StepperBase<Scalar> > &stepper,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(is_null(stepper));
  TEST_FOR_EXCEPT(is_null(stepper->getModel()));
#endif
  stepper_ = stepper;
  initialCondition_ = initialCondition;
  currentInitialCondition_ = initialCondition;

  const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> >
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
Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> >
StepperAsModelEvaluator<Scalar>::get_p_space(int l) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_INTEGRAL_IN_RANGE( int, l, 0, Np_ );
#endif
  return p_space_[l];
}


template<class Scalar>
Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> >
StepperAsModelEvaluator<Scalar>::get_g_space(int j) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_INTEGRAL_IN_RANGE( int, j, 0, Ng_ );
#endif
  return g_space_[j];
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
StepperAsModelEvaluator<Scalar>::createInArgs() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.set_Np(Np_);
  inArgs.setSupports(MEB::IN_ARG_t);
  return inArgs;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
StepperAsModelEvaluator<Scalar>::createOutArgs() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.set_Np_Ng(Np_,Ng_);
  return outArgs;
}


template<class Scalar>
void StepperAsModelEvaluator<Scalar>::evalModel(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs
  ) const
{

  using Teuchos::as;
  using Teuchos::RefCountPtr;
  using Teuchos::describe;

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_GEN_BEGIN(
    "Rythmos::StepperAsModelEvaluator", inArgs, outArgs, Teuchos::null
    );

  // InArgs

  const Scalar finalTime = inArgs.get_t();
  for ( int l = 0; l < Np_; ++l ) {
    currentInitialCondition_.set_p(l,inArgs.get_p(l));
  }

  // OutArgs

  RefCountPtr<Thyra::VectorBase<Scalar> >
    g_out = outArgs.get_g(0);

  TEST_FOR_EXCEPT(
    is_null(g_out) && "You must ask for g(0) when you call this function!"
    );

#ifdef TEUCHOS_DEBUG

  THYRA_ASSERT_VEC_SPACES(
    "StepperAsModelEvaluator<Scalar>::evalModel(...)",
    *g_out->space(), *stepper_->get_x_space() );

#endif
  
  // Do the integration
  
  Scalar t0 = currentInitialCondition_.get_t();
  Scalar dt = (finalTime-t0)/numTimeSteps();
  Scalar time = t0;

  stepper_->setInitialCondition(currentInitialCondition_);

  for ( int i = 1 ; i <= numTimeSteps() ; ++i ) {

    if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW) )
      *out << "\ntime step = " << i << ", time = " << time << ":\n";
    
    double dt_taken = stepper_->takeStep(dt,Rythmos::FIXED_STEP);

    TEST_FOR_EXCEPTION(
      dt_taken != dt, std::logic_error,
      "Error, stepper took step of dt = " << dt_taken 
      << " when asked to take step of dt = " << dt << "\n"
      );
    
    time += dt_taken;

    StepStatus<Scalar>
      stepStatus = stepper_->getStepStatus();
    
    RefCountPtr<const Thyra::VectorBase<Scalar> >
      solution = stepStatus.solution,
      solutionDot = stepStatus.solutionDot;
    
    if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_EXTREME) ) {
      *out << "\nsolution = \n" << describe(*solution,verbLevel);
      *out << "\nsolutionDot = \n" << describe(*solutionDot,verbLevel);
    }
    
  }

  // Retrieve the final solution

  if (!is_null(g_out)) {
    Thyra::assign(
      &*g_out,
      *stepper_->getStepStatus().solution
      );
  }

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_END();

}


} // namespace Rythmos


#endif // RYTHMOS_STEPPER_AS_MODEL_EVALUATOR_HPP
