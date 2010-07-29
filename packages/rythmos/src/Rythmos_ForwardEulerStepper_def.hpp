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

#ifndef Rythmos_FORWARDEULER_STEPPER_DEF_H
#define Rythmos_FORWARDEULER_STEPPER_DEF_H

#include "Rythmos_ForwardEulerStepper_decl.hpp"
#include "Rythmos_StepperHelpers.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"


namespace Rythmos {


// --------------------------------------------------------------------
// ForwardEulerStepperMomento definitions:
// --------------------------------------------------------------------

template<class Scalar>
ForwardEulerStepperMomento<Scalar>::ForwardEulerStepperMomento() 
{}

template<class Scalar>
ForwardEulerStepperMomento<Scalar>::~ForwardEulerStepperMomento() 
{}

template<class Scalar>
void ForwardEulerStepperMomento<Scalar>::serialize(
        const StateSerializerStrategy<Scalar>& stateSerializer,
        std::ostream& oStream
        ) const
{
  using Teuchos::is_null;
  TEUCHOS_ASSERT( !is_null(model_) );
  RCP<VectorBase<Scalar> > sol_vec = solution_vector_;
  if (is_null(sol_vec)) {
    sol_vec = Thyra::createMember(model_->get_x_space());
  }
  RCP<VectorBase<Scalar> > res_vec = residual_vector_;
  if (is_null(res_vec)) {
    res_vec = Thyra::createMember(model_->get_f_space());
  }
  RCP<VectorBase<Scalar> > sol_vec_old = solution_vector_old_;
  if (is_null(sol_vec_old)) {
    sol_vec_old = Thyra::createMember(model_->get_x_space());
  }
  stateSerializer.serializeVectorBase(*sol_vec,oStream);
  stateSerializer.serializeVectorBase(*res_vec,oStream);
  stateSerializer.serializeVectorBase(*sol_vec_old,oStream);
  stateSerializer.serializeScalar(t_,oStream);
  stateSerializer.serializeScalar(t_old_,oStream);
  stateSerializer.serializeScalar(dt_,oStream);
  stateSerializer.serializeInt(numSteps_,oStream);
  stateSerializer.serializeBool(isInitialized_,oStream);
  stateSerializer.serializeBool(haveInitialCondition_,oStream);
  RCP<ParameterList> pl = parameterList_;
  if (Teuchos::is_null(pl)) {
    pl = Teuchos::parameterList();
  }
  stateSerializer.serializeParameterList(*pl,oStream);
}

template<class Scalar>
void ForwardEulerStepperMomento<Scalar>::deSerialize(
        const StateSerializerStrategy<Scalar>& stateSerializer,
        std::istream& iStream
        )
{
  using Teuchos::outArg;
  using Teuchos::is_null;
  TEUCHOS_ASSERT( !is_null(model_) );
  if (is_null(solution_vector_)) {
    solution_vector_ = Thyra::createMember(*model_->get_x_space());
  }
  if (is_null(residual_vector_)) {
    residual_vector_ = Thyra::createMember(*model_->get_f_space());
  }
  if (is_null(solution_vector_old_)) {
    solution_vector_old_ = Thyra::createMember(*model_->get_x_space());
  }
  stateSerializer.deSerializeVectorBase(outArg(*solution_vector_),iStream);
  stateSerializer.deSerializeVectorBase(outArg(*residual_vector_),iStream);
  stateSerializer.deSerializeVectorBase(outArg(*solution_vector_old_),iStream);
  stateSerializer.deSerializeScalar(outArg(t_),iStream);
  stateSerializer.deSerializeScalar(outArg(t_old_),iStream);
  stateSerializer.deSerializeScalar(outArg(dt_),iStream);
  stateSerializer.deSerializeInt(outArg(numSteps_),iStream);
  stateSerializer.deSerializeBool(outArg(isInitialized_),iStream);
  stateSerializer.deSerializeBool(outArg(haveInitialCondition_),iStream);
  if (is_null(parameterList_)) {
    parameterList_ = Teuchos::parameterList();
  }
  stateSerializer.deSerializeParameterList(outArg(*parameterList_),iStream);
}

template<class Scalar>
RCP<MomentoBase<Scalar> > ForwardEulerStepperMomento<Scalar>::clone() const
{
  RCP<ForwardEulerStepperMomento<Scalar> > m = rcp(new ForwardEulerStepperMomento<Scalar>());
  m->set_solution_vector(solution_vector_);
  m->set_residual_vector(residual_vector_);
  m->set_solution_vector_old(solution_vector_old_);
  m->set_t(t_);
  m->set_t_old(t_old_);
  m->set_dt(dt_);
  m->set_numSteps(numSteps_);
  m->set_isInitialized(isInitialized_);
  m->set_haveInitialCondition(haveInitialCondition_);
  m->set_parameterList(parameterList_);
  if (!Teuchos::is_null(this->getMyParamList())) {
    m->setParameterList(Teuchos::parameterList(*(this->getMyParamList())));
  }
  m->set_model(model_);
  m->set_basePoint(basePoint_);
  // How do I copy the VerboseObject data?  
  // 07/10/09 tscoffe:  Its not set up in Teuchos to do this yet
  return m;
}

template<class Scalar>
void ForwardEulerStepperMomento<Scalar>::set_solution_vector(const RCP<const VectorBase<Scalar> >& solution_vector )
{ 
  solution_vector_ = Teuchos::null;
  if (!Teuchos::is_null(solution_vector)) {
    solution_vector_ = solution_vector->clone_v(); 
  }
}

template<class Scalar>
RCP<VectorBase<Scalar> > ForwardEulerStepperMomento<Scalar>::get_solution_vector() const
{ 
  return solution_vector_; 
}

template<class Scalar>
void ForwardEulerStepperMomento<Scalar>::set_residual_vector(const RCP<const VectorBase<Scalar> >& residual_vector)
{ 
  residual_vector_ = Teuchos::null;
  if (!Teuchos::is_null(residual_vector)) {
    residual_vector_ = residual_vector->clone_v(); 
  }
}

template<class Scalar>
RCP<VectorBase<Scalar> > ForwardEulerStepperMomento<Scalar>::get_residual_vector() const
{ 
  return residual_vector_; 
}

template<class Scalar>
void ForwardEulerStepperMomento<Scalar>::set_solution_vector_old(const RCP<const VectorBase<Scalar> >& solution_vector_old )
{ 
  solution_vector_old_ = Teuchos::null;
  if (!Teuchos::is_null(solution_vector_old)) {
    solution_vector_old_ = solution_vector_old->clone_v(); 
  }
}

template<class Scalar>
RCP<VectorBase<Scalar> > ForwardEulerStepperMomento<Scalar>::get_solution_vector_old() const
{ 
  return solution_vector_old_; 
}

template<class Scalar>
void ForwardEulerStepperMomento<Scalar>::set_t(const Scalar & t)
{ 
  t_ = t; 
}

template<class Scalar>
Scalar ForwardEulerStepperMomento<Scalar>::get_t() const
{ 
  return t_; 
}

template<class Scalar>
void ForwardEulerStepperMomento<Scalar>::set_t_old(const Scalar & t_old)
{ 
  t_old_ = t_old; 
}

template<class Scalar>
Scalar ForwardEulerStepperMomento<Scalar>::get_t_old() const
{ 
  return t_old_; 
}

template<class Scalar>
void ForwardEulerStepperMomento<Scalar>::set_dt(const Scalar & dt)
{ 
  dt_ = dt; 
}

template<class Scalar>
Scalar ForwardEulerStepperMomento<Scalar>::get_dt() const
{ 
  return dt_; 
}

template<class Scalar>
void ForwardEulerStepperMomento<Scalar>::set_numSteps(const int & numSteps)
{ 
  numSteps_ = numSteps; 
}

template<class Scalar>
int ForwardEulerStepperMomento<Scalar>::get_numSteps() const
{ 
  return numSteps_; 
}

template<class Scalar>
void ForwardEulerStepperMomento<Scalar>::set_isInitialized(const bool & isInitialized)
{ 
  isInitialized_ = isInitialized; 
}

template<class Scalar>
bool ForwardEulerStepperMomento<Scalar>::get_isInitialized() const
{ 
  return isInitialized_; 
}

template<class Scalar>
void ForwardEulerStepperMomento<Scalar>::set_haveInitialCondition(const bool & haveInitialCondition)
{ 
  haveInitialCondition_ = haveInitialCondition; 
}

template<class Scalar>
bool ForwardEulerStepperMomento<Scalar>::get_haveInitialCondition() const
{ 
  return haveInitialCondition_; 
}

template<class Scalar>
void ForwardEulerStepperMomento<Scalar>::set_parameterList(const RCP<const ParameterList>& pl)
{ 
  parameterList_ = Teuchos::null;
  if (!Teuchos::is_null(pl)) {
    parameterList_ = Teuchos::parameterList(*pl); 
  }
}

template<class Scalar>
RCP<ParameterList> ForwardEulerStepperMomento<Scalar>::get_parameterList() const
{ 
  return parameterList_; 
}

template<class Scalar>
void ForwardEulerStepperMomento<Scalar>::setParameterList(const RCP<ParameterList>& paramList)
{ 
  this->setMyParamList(paramList); 
}

template<class Scalar>
RCP<const ParameterList> ForwardEulerStepperMomento<Scalar>::getValidParameters() const
{ 
  return Teuchos::null; 
}

template<class Scalar>
void ForwardEulerStepperMomento<Scalar>::set_model(const RCP<const Thyra::ModelEvaluator<Scalar> >& model)
{ 
  model_ = model; 
}

template<class Scalar>
RCP<const Thyra::ModelEvaluator<Scalar> > ForwardEulerStepperMomento<Scalar>::get_model() const
{ 
  return model_; 
}

template<class Scalar>
void ForwardEulerStepperMomento<Scalar>::set_basePoint(const RCP<const Thyra::ModelEvaluatorBase::InArgs<Scalar> >& basePoint)
{ 
  basePoint_ = basePoint; 
}

template<class Scalar>
RCP<const Thyra::ModelEvaluatorBase::InArgs<Scalar> > ForwardEulerStepperMomento<Scalar>::get_basePoint() const
{ 
  return basePoint_; 
}

// --------------------------------------------------------------------
// ForwardEulerStepper definitions:
// --------------------------------------------------------------------

// Nonmember constructor
template<class Scalar>
RCP<ForwardEulerStepper<Scalar> > forwardEulerStepper()
{
  RCP<ForwardEulerStepper<Scalar> > stepper = rcp(new ForwardEulerStepper<Scalar>());
  return stepper;
}

// Nonmember constructor
template<class Scalar>
RCP<ForwardEulerStepper<Scalar> > forwardEulerStepper(
    const RCP<Thyra::ModelEvaluator<Scalar> >& model
    )
{
  RCP<ForwardEulerStepper<Scalar> > stepper = forwardEulerStepper<Scalar>();
  {
    RCP<StepperBase<Scalar> > stepperBase = 
      Teuchos::rcp_dynamic_cast<StepperBase<Scalar> >(stepper,true);
    setStepperModel(Teuchos::inOutArg(*stepperBase),model);
  }
  return stepper;
}

template<class Scalar>
ForwardEulerStepper<Scalar>::ForwardEulerStepper()
{
  this->defaultInitializAll_();
  dt_ = Teuchos::ScalarTraits<Scalar>::zero();
  numSteps_ = 0;
}

template<class Scalar>
void ForwardEulerStepper<Scalar>::defaultInitializAll_()
{
  model_ = Teuchos::null;
  solution_vector_ = Teuchos::null;
  residual_vector_ = Teuchos::null;
  t_ = ST::nan();
  dt_ = ST::nan();
  t_old_ = ST::nan();
  solution_vector_old_ = Teuchos::null;
  //basePoint_;
  numSteps_ = -1;
  haveInitialCondition_ = false;
  parameterList_ = Teuchos::null;
  isInitialized_ = false;
}

template<class Scalar>
void ForwardEulerStepper<Scalar>::initialize_()
{
  if (!isInitialized_) {
    TEST_FOR_EXCEPTION( is_null(model_), std::logic_error,
       "Error!  Please set a model on the stepper.\n" 
       );
    residual_vector_ = Thyra::createMember(model_->get_f_space());
    isInitialized_ = true;
  }
}

template<class Scalar>
ForwardEulerStepper<Scalar>::~ForwardEulerStepper()
{
}

template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> > ForwardEulerStepper<Scalar>::get_x_space() const
{
  TEST_FOR_EXCEPTION(!haveInitialCondition_,std::logic_error,"Error, attempting to call get_x_space before setting an initial condition!\n");
  return(solution_vector_->space());
}

template<class Scalar>
Scalar ForwardEulerStepper<Scalar>::takeStep(Scalar dt, StepSizeType flag)
{
  TEST_FOR_EXCEPTION( !haveInitialCondition_, std::logic_error,
     "Error!  Attempting to call takeStep before setting an initial condition!\n" 
     );
  this->initialize_();
  if (flag == STEP_TYPE_VARIABLE) { 
    // print something out about this method not supporting automatic variable step-size
    return(-ST::one());
  }
  //Thyra::eval_f<Scalar>(*model_,*solution_vector_,t_+dt,&*residual_vector_);
  eval_model_explicit<Scalar>(*model_,basePoint_,*solution_vector_,t_+dt,Teuchos::outArg(*residual_vector_));

  // solution_vector_old_ = solution_vector_
  Thyra::V_V(Teuchos::outArg(*solution_vector_old_),*solution_vector_);
  // solution_vector = solution_vector + dt*residual_vector
  Thyra::Vp_StV(solution_vector_.ptr(),dt,*residual_vector_); 
  t_old_ = t_;
  t_ += dt;
  dt_ = dt;
  numSteps_++;

  return(dt);
}

template<class Scalar>
const StepStatus<Scalar> ForwardEulerStepper<Scalar>::getStepStatus() const
{
  StepStatus<Scalar> stepStatus;

  if (!haveInitialCondition_) {
    stepStatus.stepStatus = STEP_STATUS_UNINITIALIZED;
  } 
  else if (numSteps_ == 0) {
    stepStatus.stepStatus = STEP_STATUS_UNKNOWN;
    stepStatus.order = this->getOrder();
    stepStatus.time = t_;
    stepStatus.solution = solution_vector_;
  } 
  else {
    stepStatus.stepStatus = STEP_STATUS_CONVERGED;
    stepStatus.stepSize = dt_; 
    stepStatus.order = this->getOrder();
    stepStatus.time = t_;
    stepStatus.stepLETValue = Scalar(-ST::one()); 
    stepStatus.solution = solution_vector_;
    stepStatus.residual = residual_vector_;
  }

  return(stepStatus);
}

template<class Scalar>
std::string ForwardEulerStepper<Scalar>::description() const
{
  std::string name = "Rythmos::ForwardEulerStepper";
  return(name);
}

template<class Scalar>
void ForwardEulerStepper<Scalar>::describe(
      Teuchos::FancyOStream            &out,
      const Teuchos::EVerbosityLevel   verbLevel
      ) const
{
  if ( (static_cast<int>(verbLevel) == static_cast<int>(Teuchos::VERB_DEFAULT) ) ||
       (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW)     )
     ) {
    out << description() << "::describe" << std::endl;
    out << "model = " << model_->description() << std::endl;
  } else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW)) {
  } else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_MEDIUM)) {
  } else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_HIGH)) {
    out << "model = " << std::endl;
    model_->describe(out,verbLevel);
    out << "solution_vector = " << std::endl;
    solution_vector_->describe(out,verbLevel);
    out << "residual_vector = " << std::endl;
    residual_vector_->describe(out,verbLevel);
  }
}

template<class Scalar>
void ForwardEulerStepper<Scalar>::addPoints(
    const Array<Scalar>& time_vec
    ,const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x_vec
    ,const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec
    )
{
  TEST_FOR_EXCEPTION(true,std::logic_error,"Error, addPoints is not implemented for ForwardEulerStepper.\n");
}

template<class Scalar>
void ForwardEulerStepper<Scalar>::getPoints(
    const Array<Scalar>& time_vec
    ,Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >* x_vec
    ,Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec
    ,Array<ScalarMag>* accuracy_vec) const
{
  TEUCHOS_ASSERT( haveInitialCondition_ );
  using Teuchos::constOptInArg;
  using Teuchos::null;
  defaultGetPoints<Scalar>(
      t_old_,
      constOptInArg(*solution_vector_old_),
      Ptr<const VectorBase<Scalar> >(null),
      t_,
      constOptInArg(*solution_vector_),
      Ptr<const VectorBase<Scalar> >(null),
      time_vec,
      ptr(x_vec),
      ptr(xdot_vec),
      ptr(accuracy_vec),
      Ptr<InterpolatorBase<Scalar> >(null)
      );
}

template<class Scalar>
TimeRange<Scalar> ForwardEulerStepper<Scalar>::getTimeRange() const
{
  if (!haveInitialCondition_) {
    return(invalidTimeRange<Scalar>());
  } else {
    return(TimeRange<Scalar>(t_old_,t_));
  }
}

template<class Scalar>
void ForwardEulerStepper<Scalar>::getNodes(Array<Scalar>* time_vec) const
{
  TEUCHOS_ASSERT( time_vec != NULL );
  time_vec->clear();
  if (!haveInitialCondition_) {
    return; 
  } else {
    time_vec->push_back(t_old_);
  }
  if (numSteps_ > 0) {
    time_vec->push_back(t_);
  }
}

template<class Scalar>
void ForwardEulerStepper<Scalar>::removeNodes(Array<Scalar>& time_vec) 
{
  TEST_FOR_EXCEPTION(true,std::logic_error,"Error, removeNodes is not implemented for ForwardEulerStepper.\n");
}

template<class Scalar>
int ForwardEulerStepper<Scalar>::getOrder() const
{
  return(1);
}

template <class Scalar>
void ForwardEulerStepper<Scalar>::setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList)
{
  TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParametersAndSetDefaults(*this->getValidParameters());
  parameterList_ = paramList;
  Teuchos::readVerboseObjectSublist(&*parameterList_,this);
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList> ForwardEulerStepper<Scalar>::getNonconstParameterList()
{
  return(parameterList_);
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList> ForwardEulerStepper<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_param_list = parameterList_;
  parameterList_ = Teuchos::null;
  return(temp_param_list);
}


template<class Scalar>
RCP<const Teuchos::ParameterList>
ForwardEulerStepper<Scalar>::getValidParameters() const
{
  using Teuchos::ParameterList;
  static RCP<const ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<ParameterList> pl = Teuchos::parameterList();
    Teuchos::setupVerboseObjectSublist(&*pl);
    validPL = pl;
  }
  return validPL;
}


template<class Scalar>
void ForwardEulerStepper<Scalar>::setModel(const RCP<const Thyra::ModelEvaluator<Scalar> >& model)
{
  TEST_FOR_EXCEPT( is_null(model) );
  assertValidModel( *this, *model );
  model_ = model;
}


template<class Scalar>
void ForwardEulerStepper<Scalar>::setNonconstModel(const RCP<Thyra::ModelEvaluator<Scalar> >& model)
{
  this->setModel(model); // TODO 09/09/09 tscoffe:  use ConstNonconstObjectContainer!
}


template<class Scalar>
Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
ForwardEulerStepper<Scalar>::getModel() const
{
  return model_;
}

template<class Scalar>
Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >
ForwardEulerStepper<Scalar>::getNonconstModel() 
{
  return Teuchos::null;
}

template<class Scalar>
void ForwardEulerStepper<Scalar>::setInitialCondition(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
    )
{
  basePoint_ = initialCondition;
  t_ = initialCondition.get_t();
  t_old_ = t_;
  solution_vector_ = initialCondition.get_x()->clone_v();
  solution_vector_old_ = solution_vector_->clone_v();
  haveInitialCondition_ = true;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar> 
ForwardEulerStepper<Scalar>::getInitialCondition() const
{
  return basePoint_;
}


template<class Scalar>
bool ForwardEulerStepper<Scalar>::supportsCloning() const
{
  return true;
}

template<class Scalar>
RCP<StepperBase<Scalar> >
ForwardEulerStepper<Scalar>::cloneStepperAlgorithm() const
{

  // Just use the interface to clone the algorithm in a basically
  // uninitialized state

  RCP<StepperBase<Scalar> >
    stepper = Teuchos::rcp(new ForwardEulerStepper<Scalar>());

  if (!is_null(model_)) {
    setStepperModel(Teuchos::inOutArg(*stepper),model_); // Shallow copy is okay!
  }

  if (!is_null(parameterList_)) {
    stepper->setParameterList(Teuchos::parameterList(*parameterList_));
  }

  return stepper;

}

template<class Scalar>
RCP<const MomentoBase<Scalar> >
ForwardEulerStepper<Scalar>::getMomento() const
{
  RCP<ForwardEulerStepperMomento<Scalar> > momento = Teuchos::rcp(new ForwardEulerStepperMomento<Scalar>());
  momento->set_solution_vector(solution_vector_);
  momento->set_solution_vector_old(solution_vector_old_);
  momento->set_residual_vector(residual_vector_);
  momento->set_t(t_);
  momento->set_t_old(t_old_);
  momento->set_dt(dt_);
  momento->set_numSteps(numSteps_);
  momento->set_isInitialized(isInitialized_);
  momento->set_haveInitialCondition(haveInitialCondition_);
  momento->set_parameterList(parameterList_);
  momento->set_model(model_);
  RCP<Thyra::ModelEvaluatorBase::InArgs<Scalar> > bp = Teuchos::rcp(new Thyra::ModelEvaluatorBase::InArgs<Scalar>(basePoint_));
  momento->set_basePoint(bp);
  return momento;
}

template<class Scalar>
void ForwardEulerStepper<Scalar>::setMomento( const Ptr<const MomentoBase<Scalar> >& momentoPtr ) 
{ 
  Ptr<const ForwardEulerStepperMomento<Scalar> > feMomentoPtr = 
    Teuchos::ptr_dynamic_cast<const ForwardEulerStepperMomento<Scalar> >(momentoPtr,true);
  TEST_FOR_EXCEPTION( Teuchos::is_null(feMomentoPtr->get_model()), std::logic_error,
      "Error!  Rythmos::ForwardEulerStepper::setMomento:  The momento must have a valid model through set_model(...) prior to calling ForwardEulerStepper::setMomento(...)."
      );
  TEST_FOR_EXCEPTION( Teuchos::is_null(feMomentoPtr->get_basePoint()), std::logic_error,
      "Error!  Rythmos::ForwardEulerStepper::setMomento:  The momento must have a valid base point through set_basePoint(...) prior to calling ForwardEulerStepper::setMomento(...)."
      );
  model_ = feMomentoPtr->get_model();
  basePoint_ = *(feMomentoPtr->get_basePoint());
  const ForwardEulerStepperMomento<Scalar>& feMomento = *feMomentoPtr;
  solution_vector_ = feMomento.get_solution_vector();
  solution_vector_old_ = feMomento.get_solution_vector_old();
  residual_vector_ = feMomento.get_residual_vector();
  t_ = feMomento.get_t();
  t_old_ = feMomento.get_t_old();
  dt_ = feMomento.get_dt();
  numSteps_ = feMomento.get_numSteps();
  isInitialized_ = feMomento.get_isInitialized();
  haveInitialCondition_ = feMomento.get_haveInitialCondition();
  parameterList_ = feMomento.get_parameterList();
  this->checkConsistentState_();
}

template<class Scalar>
void ForwardEulerStepper<Scalar>::checkConsistentState_()
{
  if (isInitialized_) {
    TEUCHOS_ASSERT( !Teuchos::is_null(model_) );
    TEUCHOS_ASSERT( !Teuchos::is_null(residual_vector_) );
  }
  if (haveInitialCondition_) {
    TEUCHOS_ASSERT( !ST::isnaninf(t_) );
    TEUCHOS_ASSERT( !ST::isnaninf(t_old_) );
    TEUCHOS_ASSERT( !Teuchos::is_null(solution_vector_) );
    TEUCHOS_ASSERT( !Teuchos::is_null(solution_vector_old_) );
    TEUCHOS_ASSERT( t_ >= basePoint_.get_t() );
    TEUCHOS_ASSERT( t_old_ >= basePoint_.get_t() );
  }
  if (numSteps_ > 0) {
    TEUCHOS_ASSERT(isInitialized_);
    TEUCHOS_ASSERT(haveInitialCondition_);
  }
}

// 
// Explicit Instantiation macro
//
// Must be expanded from within the Rythmos namespace!
//

#define RYTHMOS_FORWARD_EULER_STEPPER_INSTANT(SCALAR) \
  \
  template class ForwardEulerStepperMomento< SCALAR >; \
  template class ForwardEulerStepper< SCALAR >; \
  \
  template RCP<ForwardEulerStepper< SCALAR > > forwardEulerStepper(); \
  template RCP<ForwardEulerStepper< SCALAR > > forwardEulerStepper( \
      const RCP<Thyra::ModelEvaluator< SCALAR > >& model \
      ); 


} // namespace Rythmos

#endif //Rythmos_FORWARDEULER_STEPPER_DEF_H
