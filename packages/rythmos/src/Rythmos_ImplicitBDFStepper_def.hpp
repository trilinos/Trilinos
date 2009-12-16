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

#ifndef Rythmos_IMPLICITBDF_STEPPER_DEF_H
#define Rythmos_IMPLICITBDF_STEPPER_DEF_H

#include "Rythmos_ImplicitBDFStepper_decl.hpp"
#include "Rythmos_StepperHelpers.hpp"
#include "Rythmos_ImplicitBDFStepperStepControl.hpp"

namespace Rythmos {

// ////////////////////////////
// Defintions

// Nonmember constructor
template<class Scalar>
RCP<ImplicitBDFStepper<Scalar> > implicitBDFStepper() {
  RCP<ImplicitBDFStepper<Scalar> > stepper = rcp(new ImplicitBDFStepper<Scalar>() );
  return stepper;
}

template<class Scalar>
RCP<ImplicitBDFStepper<Scalar> > implicitBDFStepper(
  const RCP<Thyra::ModelEvaluator<Scalar> >& model,
  const RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
  const RCP<Teuchos::ParameterList>& parameterList
  )
{
  RCP<ImplicitBDFStepper<Scalar> > stepper = Teuchos::rcp(new ImplicitBDFStepper<Scalar>(model,solver,parameterList));
  return stepper;
}

// Constructors, intializers, Misc.


template<class Scalar>
ImplicitBDFStepper<Scalar>::ImplicitBDFStepper()
{
  this->defaultInitializeAll_();
  haveInitialCondition_ = false;
  isInitialized_=false;
}


template<class Scalar>
ImplicitBDFStepper<Scalar>::ImplicitBDFStepper(
  const RCP<Thyra::ModelEvaluator<Scalar> >& model
  ,const RCP<Thyra::NonlinearSolverBase<Scalar> >& solver
  ,const RCP<Teuchos::ParameterList>& parameterList
  )
{
  this->defaultInitializeAll_();
  this->setParameterList(parameterList);
  // Now we instantiate the model and the solver
  setModel(model);
  setSolver(solver);
  haveInitialCondition_ = false;
  isInitialized_=false;
}


template<class Scalar>
ImplicitBDFStepper<Scalar>::ImplicitBDFStepper(
  const RCP<Thyra::ModelEvaluator<Scalar> >& model
  ,const RCP<Thyra::NonlinearSolverBase<Scalar> >& solver
  )
{
  this->defaultInitializeAll_();
  // Now we instantiate the model and the solver
  setModel(model);
  setSolver(solver);
  haveInitialCondition_ = false;
  isInitialized_=false;
}

template<class Scalar>
const Thyra::VectorBase<Scalar>& 
  ImplicitBDFStepper<Scalar>::getxHistory(int index) const
{
  TEST_FOR_EXCEPTION(!isInitialized_,std::logic_error,
      "Error, attempting to call getxHistory before initialization!\n");
  TEST_FOR_EXCEPT( !(( 0 <= index ) && ( index <= maxOrder_ )) );
  TEST_FOR_EXCEPT( !( index <= usedOrder_+1 ) );
  return(*(xHistory_[index]));
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::setStepControlStrategy(const RCP<StepControlStrategyBase<Scalar> >& stepControl)
{
  TEST_FOR_EXCEPTION(stepControl == Teuchos::null,std::logic_error,"Error, stepControl == Teuchos::null!\n");
  stepControl_ = stepControl;    
}

template<class Scalar>
RCP<StepControlStrategyBase<Scalar> > ImplicitBDFStepper<Scalar>::getNonconstStepControlStrategy() 
{
  return(stepControl_);
}

template<class Scalar>
RCP<const StepControlStrategyBase<Scalar> > ImplicitBDFStepper<Scalar>::getStepControlStrategy() const
{
  return(stepControl_);
}


// Overridden from SolverAcceptingStepperBase


template<class Scalar>
void ImplicitBDFStepper<Scalar>::setSolver(const RCP<Thyra::NonlinearSolverBase<Scalar> > &solver)
{
  TEST_FOR_EXCEPT(solver == Teuchos::null)
    solver_ = solver;
}


template<class Scalar>
RCP<Thyra::NonlinearSolverBase<Scalar> >
ImplicitBDFStepper<Scalar>::getNonconstSolver()
{
  return (solver_);
}


template<class Scalar>
RCP<const Thyra::NonlinearSolverBase<Scalar> >
ImplicitBDFStepper<Scalar>::getSolver() const
{
  return (solver_);
}


// Overridden from StepperBase


template<class Scalar>
bool ImplicitBDFStepper<Scalar>::isImplicit() const
{
  return true;
}

template<class Scalar>
bool ImplicitBDFStepper<Scalar>::supportsCloning() const
{
  return true;
}


template<class Scalar>
RCP<StepperBase<Scalar> >
ImplicitBDFStepper<Scalar>::cloneStepperAlgorithm() const
{

  // Just use the interface to clone the algorithm in an basically
  // uninitialized state

  RCP<ImplicitBDFStepper<Scalar> >
    stepper = Teuchos::rcp(new ImplicitBDFStepper<Scalar>());

  if (!is_null(model_))
    stepper->setModel(model_); // Shallow copy is okay!

  if (!is_null(solver_))
    stepper->setSolver(solver_->cloneNonlinearSolver().assert_not_null());
  
  if (!is_null(parameterList_))
    stepper->setParameterList(Teuchos::parameterList(*parameterList_));

  if (!is_null(stepControl_)) {
    if (stepControl_->supportsCloning())
      stepper->setStepControlStrategy(
        stepControl_->cloneStepControlStrategyAlgorithm().assert_not_null());
  }
  
  // At this point, we should have a valid algorithm state.  What might be
  // missing is the initial condition if it was not given in *model_ but was
  // set explicitly.  However, the specification for this function does not
  // guarantee that the full state will be copied in any case!

  return stepper;

}


template<class Scalar>
void ImplicitBDFStepper<Scalar>::setModel(
  const RCP<const Thyra::ModelEvaluator<Scalar> >& model
  )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  TEST_FOR_EXCEPT( is_null(model) );
  assertValidModel( *this, *model );
  model_ = model;
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::setNonconstModel(
  const RCP<Thyra::ModelEvaluator<Scalar> >& model
  )
{
  this->setModel(model); // TODO 09/09/09 tscoffe:  use ConstNonconstObjectContainer!
}


template<class Scalar>
RCP<const Thyra::ModelEvaluator<Scalar> >
ImplicitBDFStepper<Scalar>::getModel() const
{
  return model_;
}


template<class Scalar>
RCP<Thyra::ModelEvaluator<Scalar> >
ImplicitBDFStepper<Scalar>::getNonconstModel()
{
  return Teuchos::null;
}


template<class Scalar>
void ImplicitBDFStepper<Scalar>::setInitialCondition(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
  )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Thyra::ModelEvaluatorBase MEB;
  TEST_FOR_EXCEPT(is_null(initialCondition.get_x()));
  TEST_FOR_EXCEPT(is_null(initialCondition.get_x_dot()));
  basePoint_ = initialCondition;
  xn0_ = initialCondition.get_x()->clone_v();
  xpn0_ = initialCondition.get_x_dot()->clone_v(); 
  time_ = initialCondition.get_t();
  // Generate vectors for use in calculations
  x_dot_base_ = createMember(xpn0_->space());
  V_S(&*x_dot_base_,ST::zero());
  ee_ = createMember(xn0_->space());
  V_S(&*ee_,ST::zero());
  // x history
  xHistory_.clear();
  xHistory_.push_back(xn0_->clone_v());
  xHistory_.push_back(xpn0_->clone_v());

  haveInitialCondition_ = true;
  if (isInitialized_) {
    initialize_();
  }
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar> 
ImplicitBDFStepper<Scalar>::getInitialCondition() const
{
  return basePoint_;
}


template<class Scalar>
Scalar ImplicitBDFStepper<Scalar>::takeStep(Scalar dt, StepSizeType stepType)
{
  
#ifdef ENABLE_RYTHMOS_TIMERS
  TEUCHOS_FUNC_TIME_MONITOR("Rythmos::ImplicitBDFStepper::takeStep");
#endif
  
  using Teuchos::as;
  using Teuchos::incrVerbLevel;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename Thyra::ModelEvaluatorBase::InArgs<Scalar>::ScalarMag ScalarMag;
  typedef Thyra::NonlinearSolverBase<Scalar> NSB;
  typedef Teuchos::VerboseObjectTempState<NSB> VOTSNSB;

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"takeStep");
  VOTSNSB solver_outputTempState(solver_,out,incrVerbLevel(verbLevel,-1));

  if ( !is_null(out) && as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW) ) {
    *out
      << "\nEntering " << this->Teuchos::Describable::description()
      << "::takeStep("<<dt<<","<<toString(stepType)<<") ...\n"; 
  }

  if (!isInitialized_) {
    initialize_(); 
  }

  stepControl_->setOStream(out);
  stepControl_->setVerbLevel(verbLevel);
  stepControl_->setRequestedStepSize(*this,dt,stepType);

  AttemptedStepStatusFlag status;
  while (1) {
    // Set up problem coefficients (and handle first step)
    Scalar hh_old = hh_;
    int desiredOrder;
    stepControl_->nextStepSize(*this,&hh_,&stepType,&desiredOrder);
    TEST_FOR_EXCEPT(!((1 <= desiredOrder) && (desiredOrder <= maxOrder_)));
    TEST_FOR_EXCEPT(!(desiredOrder <= usedOrder_+1));
    currentOrder_ = desiredOrder;
    if (numberOfSteps_ == 0) {
      psi_[0] = hh_;
      if (nef_ == 0) {
        Vt_S(&*xHistory_[1],hh_);
      } else {
        Vt_S(&*xHistory_[1],hh_/hh_old);
      }
    }
    this->updateCoeffs_();
    // compute predictor
    obtainPredictor_();
    // solve nonlinear problem (as follows)
    
    //
    // Setup the nonlinear equations:
    //
    //   f_bar( x_dot_coeff * x_bar + x_dot_base, x_coeff * x_bar + x_base, t_base ) = 0
    //   x_dot_coeff = -alpha_s/dt
    //   x_dot_base = x_prime_pred + (alpha_s/dt) * x_pred
    //   x_coeff = 1
    //   x_base = 0
    //   t_base = tn+dt
    //
    Scalar coeff_x_dot = Scalar(-ST::one())*alpha_s_/hh_;
    V_StVpStV( &*x_dot_base_, ST::one(), *xpn0_, alpha_s_/hh_, *xn0_ );
    if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_EXTREME) ) {
      *out << "model_ = " << std::endl;
      model_->describe(*out,verbLevel);
      *out << "basePoint_ = " << std::endl;
      basePoint_.describe(*out,verbLevel);
      *out << "coeff_x_dot = " << coeff_x_dot << std::endl;
      *out << "x_dot_base_ = " << std::endl;
      x_dot_base_->describe(*out,verbLevel);
      *out << "time_+hh_ = " << time_+hh_ << std::endl;
      *out << "xn0_ = " << std::endl;
      xn0_->describe(*out,verbLevel);
    }
    neModel_.initializeSingleResidualModel(
      model_, basePoint_,
      coeff_x_dot, x_dot_base_,
      ST::one(), Teuchos::null,
      time_+hh_,
      xn0_
      );
    //
    // Solve the implicit nonlinear system to a tolerance of ???
    // 
    if(solver_->getModel().get()!=&neModel_) {
      solver_->setModel( Teuchos::rcp(&neModel_,false) );
    }
    /* // Rythmos::TimeStepNonlinearSolver uses a built in solveCriteria, so you can't pass one in.
    // I believe this is the correct solveCriteria for IDA though.
    Thyra::SolveMeasureType nonlinear_solve_measure_type(Thyra::SOLVE_MEASURE_NORM_RESIDUAL,Thyra::SOLVE_MEASURE_ONE); 
    ScalarMag tolerance = relErrTol_/ScalarMag(10.0); // This should be changed to match the condition in IDA
    Thyra::SolveCriteria<Scalar> nonlinearSolveCriteria(nonlinear_solve_measure_type, tolerance);
    Thyra::SolveStatus<Scalar> nonlinearSolveStatus = solver_->solve( &*xn0_, &nonlinearSolveCriteria, &*delta_ ); 
    */
    if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_EXTREME) ) {
      *out << "xn0 = " << std::endl;
      xn0_->describe(*out,verbLevel);
      *out << "ee = " << std::endl;
      ee_->describe(*out,verbLevel);
    }
    Thyra::SolveStatus<Scalar> nonlinearSolveStatus = solver_->solve( &*xn0_, NULL, &*ee_ ); 
    if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_EXTREME) ) {
      *out << "xn0 = " << std::endl;
      xn0_->describe(*out,verbLevel);
      *out << "ee = " << std::endl;
      ee_->describe(*out,verbLevel);
    }
    // In the above solve, on input *xn0_ is the initial guess that comes from
    // the predictor.  On output, *xn0_ is the solved for solution and *ee_ is
    // the difference computed from the intial guess in *xn0_ to the final
    // solved value of *xn0_.  This is needed for basic numerical stability.
    if (nonlinearSolveStatus.solveStatus == Thyra::SOLVE_STATUS_CONVERGED)  {
      newtonConvergenceStatus_ = 0;
    }
    else {
      newtonConvergenceStatus_ = -1;
    }

    // check error and evaluate LTE
    stepControl_->setCorrection(*this,xn0_,ee_,newtonConvergenceStatus_);
    bool stepPass = stepControl_->acceptStep(*this,&LETvalue_);
    
    if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
      *out << "xn0_ = " << std::endl;
      xn0_->describe(*out,verbLevel);
      *out << "xpn0_ = " << std::endl;
      xpn0_->describe(*out,verbLevel);
      *out << "ee_ = " << std::endl;
      ee_->describe(*out,verbLevel);
      for (int i=0; i<std::max(2,maxOrder_); ++i) {
        *out << "xHistory_[" << i << "] = " << std::endl;
        xHistory_[i]->describe(*out,verbLevel);
      }
    }

    // Check LTE here:
    if (!stepPass) { // stepPass = false
      stepLETStatus_ = STEP_LET_STATUS_FAILED;
      status = stepControl_->rejectStep(*this);
      nef_++;
      if (status == CONTINUE_ANYWAY) {
        break;
      } else {
        restoreHistory_();
      }
    } else { // stepPass = true
      stepLETStatus_ = STEP_LET_STATUS_PASSED;
      break;
    }
  }

  // 08/22/07 the history array must be updated before stepControl_->completeStep.
  completeStep_();   
  stepControl_->completeStep(*this);

  if ( !is_null(out) && as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW) ) {
    *out
      << "\nLeaving " << this->Teuchos::Describable::description()
      << "::takeStep("<<dt<<","<<toString(stepType)<<") ...\n"; 
  }

  return(usedStep_);

}


template<class Scalar>
const StepStatus<Scalar> ImplicitBDFStepper<Scalar>::getStepStatus() const
{

  // 2007/08/24: rabartl: We agreed that getStepStatus() would be free
  // so I have commented out removed all code that is not free

  typedef Teuchos::ScalarTraits<Scalar> ST;
  StepStatus<Scalar> stepStatus;
  if (!isInitialized_) {
    stepStatus.message = "This stepper is uninitialized.";
    stepStatus.stepStatus = STEP_STATUS_UNINITIALIZED;
    stepStatus.stepSize = Scalar(-ST::one());
    stepStatus.order = -1;
    stepStatus.time = Scalar(-ST::one());
    return(stepStatus);
  }

  if (numberOfSteps_ > 0) {
    stepStatus.stepStatus = STEP_STATUS_CONVERGED; 
  } else {
    stepStatus.stepStatus = STEP_STATUS_UNKNOWN;
  }
  stepStatus.stepLETStatus = stepLETStatus_;
  stepStatus.stepSize = usedStep_; 
  stepStatus.order = usedOrder_;
  stepStatus.time = time_;
  stepStatus.stepLETValue = LETvalue_; 
  stepStatus.solution = xHistory_[0];
  stepStatus.solutionDot = Teuchos::null; // This is *not* free!
  stepStatus.residual = Teuchos::null; // This is *not* free!

  return(stepStatus);

}


// Overridden from InterpolationBufferBase


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
ImplicitBDFStepper<Scalar>::get_x_space() const
{
  //TEST_FOR_EXCEPTION(!isInitialized_,std::logic_error,"Error, attempting to call get_x_space before initialization!\n");
  return ( !is_null(model_) ? model_->get_x_space() : Teuchos::null );
}


template<class Scalar>
void ImplicitBDFStepper<Scalar>::addPoints(
  const Array<Scalar>& time_vec,
  const Array<RCP<const Thyra::VectorBase<Scalar> > >& x_vec,
  const Array<RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec
  )
{
  TEST_FOR_EXCEPTION(true,std::logic_error,
    "Error, addPoints is not implemented for ImplicitBDFStepper.\n");
}


template<class Scalar>
TimeRange<Scalar> ImplicitBDFStepper<Scalar>::getTimeRange() const
{
  if ( !isInitialized_ && haveInitialCondition_ )
    return timeRange<Scalar>(time_,time_);
  if ( !isInitialized_ && !haveInitialCondition_ )
    return invalidTimeRange<Scalar>();
  return timeRange<Scalar>(time_-usedStep_,time_);
}


template<class Scalar>
void ImplicitBDFStepper<Scalar>::getPoints(
  const Array<Scalar>& time_vec
  ,Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec
  ,Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec
  ,Array<ScalarMag>* accuracy_vec) const
{
  using Teuchos::as;
  using Teuchos::constOptInArg;
  using Teuchos::null;
  using Teuchos::ptr;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;

  TEUCHOS_ASSERT(haveInitialCondition_);
  // Only do this if we're being called pre-initialization to get the IC.
  if ( (numberOfSteps_ == -1) && 
       (time_vec.length() == 1) &&
       (compareTimeValues<Scalar>(time_vec[0],time_)==0) ) {
    defaultGetPoints<Scalar>(
        time_, constOptInArg(*xn0_), constOptInArg(*xpn0_),
        time_, constOptInArg(*xn0_), constOptInArg(*xpn0_),
        time_vec, ptr(x_vec), ptr(xdot_vec), ptr(accuracy_vec),
        null
        );
    return;
  }
  TEUCHOS_ASSERT(isInitialized_);
#ifdef ENABLE_RYTHMOS_TIMERS
  TEUCHOS_FUNC_TIME_MONITOR("Rythmos::ImplicitBDFStepper::getPoints");
#endif
  if (x_vec)
    x_vec->clear();
  if (xdot_vec)
    xdot_vec->clear();
  for (Teuchos::Ordinal i=0 ; i<time_vec.size() ; ++i) {
    RCP<Thyra::VectorBase<Scalar> >
      x_temp = createMember(xn0_->space());
    RCP<Thyra::VectorBase<Scalar> >
      xdot_temp = createMember(xn0_->space());
    ScalarMag accuracy = -ST::zero();
    interpolateSolution_(
      time_vec[i], &*x_temp, &*xdot_temp,
      accuracy_vec ? &accuracy : 0
      );
    if (x_vec)
      x_vec->push_back(x_temp);
    if (xdot_vec)
      xdot_vec->push_back(xdot_temp);
    if (accuracy_vec)
      accuracy_vec->push_back(accuracy);
  }
  if ( as<int>(this->getVerbLevel()) >= as<int>(Teuchos::VERB_HIGH) ) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"getPoints");
    *out << "Passing out the interpolated values:" << std::endl;
    for (Teuchos::Ordinal i=0; i<time_vec.size() ; ++i) {
      *out << "time_[" << i << "] = " << time_vec[i] << std::endl;
      if (x_vec) {
        *out << "x_vec[" << i << "] = " << std::endl;
        (*x_vec)[i]->describe(*out,this->getVerbLevel());
      }
      if (xdot_vec) {
        *out << "xdot_vec[" << i << "] = ";
        if ( (*xdot_vec)[i] == Teuchos::null) {
          *out << "Teuchos::null" << std::endl;
        }
        else {
          *out << std::endl << Teuchos::describe(*(*xdot_vec)[i],this->getVerbLevel());
        }
      }
      if (accuracy_vec)
        *out << "accuracy[" << i << "] = " << (*accuracy_vec)[i] << std::endl;
    }
  }
}


template<class Scalar>
void ImplicitBDFStepper<Scalar>::getNodes(Array<Scalar>* time_vec) const
{
  TEUCHOS_ASSERT( time_vec != NULL );
  time_vec->clear();
  if (!haveInitialCondition_) {
    return;
  }
  if (numberOfSteps_ > 0) {
    time_vec->push_back(time_-usedStep_);
  }
  time_vec->push_back(time_);
}


template<class Scalar>
void ImplicitBDFStepper<Scalar>::removeNodes(Array<Scalar>& time_vec) 
{
  TEST_FOR_EXCEPTION(true,std::logic_error,
    "Error, removeNodes is not implemented for ImplicitBDFStepper.\n");
}


template<class Scalar>
int ImplicitBDFStepper<Scalar>::getOrder() const
{
  if (!isInitialized_) {
    return(-1);
  }
  return(usedOrder_);
}


// Overridden from Teuchos::ParameterListAcceptor


template<class Scalar>
void ImplicitBDFStepper<Scalar>::setParameterList(
  RCP<Teuchos::ParameterList> const& paramList
  )
{
  TEST_FOR_EXCEPT(paramList == Teuchos::null);
  paramList->validateParameters(*this->getValidParameters(),0);
  parameterList_ = paramList;
  Teuchos::readVerboseObjectSublist(&*parameterList_,this);
}


template<class Scalar>
RCP<Teuchos::ParameterList> ImplicitBDFStepper<Scalar>::getNonconstParameterList()
{
  return(parameterList_);
}


template<class Scalar>
RCP<Teuchos::ParameterList>
ImplicitBDFStepper<Scalar>::unsetParameterList()
{
  RCP<Teuchos::ParameterList> temp_param_list = parameterList_;
  parameterList_ = Teuchos::null;
  return(temp_param_list);
}


template<class Scalar>
RCP<const Teuchos::ParameterList>
ImplicitBDFStepper<Scalar>::getValidParameters() const
{

  static RCP<Teuchos::ParameterList> validPL;

  if (is_null(validPL)) {

    RCP<Teuchos::ParameterList>
      pl = Teuchos::parameterList();

    pl->sublist(RythmosStepControlSettings_name);

    Teuchos::setupVerboseObjectSublist(&*pl);

    validPL = pl;

  }

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"getValidParameters");
  if (Teuchos::as<int>(verbLevel) == Teuchos::VERB_HIGH) {
    *out << "Setting up valid parameterlist." << std::endl;
    validPL->print(*out);
  }

  return (validPL);
  
}


// Overridden from Teuchos::Describable


template<class Scalar>
std::string ImplicitBDFStepper<Scalar>::description() const
{
  std::ostringstream out;
  out << this->Teuchos::Describable::description();
  const TimeRange<Scalar> timeRange = this->getTimeRange();
  if (timeRange.isValid())
    out << " (timeRange="<<timeRange<<")";
  else
    out << " (This stepper is not initialized yet)";
  out << std::endl;
  return out.str();
}


template<class Scalar>
void ImplicitBDFStepper<Scalar>::describe(
  Teuchos::FancyOStream &out,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{

  using Teuchos::as;

  if (!isInitialized_) {
    out << this->description();
    return;
  }
  
  if ( (as<int>(verbLevel) == as<int>(Teuchos::VERB_DEFAULT) ) ||
    (as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW)     )
    )
  {
    out << this->description() << std::endl;
    out << "model_ = " << Teuchos::describe(*model_,verbLevel);
    out << "solver_ = " << Teuchos::describe(*solver_,verbLevel);
    out << "neModel_ = " << Teuchos::describe(neModel_,verbLevel);
  }
  if (as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW)) {
    out << "time_ = " << time_ << std::endl;
    out << "hh_ = " << hh_ << std::endl;
    out << "currentOrder_ = " << currentOrder_ << std::endl;
  }
  if (as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH)) {
    out << "xn0_ = " << Teuchos::describe(*xn0_,verbLevel);
    out << "xpn0_ = " << Teuchos::describe(*xpn0_,verbLevel);
    out << "x_dot_base_ = " << Teuchos::describe(*x_dot_base_,verbLevel);
    for (int i=0 ; i < std::max(2,maxOrder_) ; ++i) {
      out << "xHistory_[" << i << "] = "
          << Teuchos::describe(*xHistory_[i],verbLevel);
    }
    out << "ee_ = " << Teuchos::describe(*ee_,verbLevel);
  }
}


// private


// 2007/08/24: rabartl: Belos: We really should initialize all of this data in
// a member initialization list but since there are like three constructors
// this would mean that we would have to duplicate code (which is error prone)
// or use a macro (which is not easy to debug).  We really should remove all
// but the default constructor which then would set this data once in the
// initialization list.

template<class Scalar>
void ImplicitBDFStepper<Scalar>::defaultInitializeAll_()
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  const Scalar nan = ST::nan(), one = ST::one(), zero = ST::zero();
  // Initialize some data members to their rightful default values
  haveInitialCondition_ = false;
  isInitialized_ = false;
  currentOrder_ = 1;
  usedOrder_ = 1;
  usedStep_ = zero;
  // Initialize the rest of the private data members to invalid values to
  // avoid uninitialed memory
  time_ = nan;
  hh_ = nan;
  maxOrder_ = -1;
  LETvalue_ = -one;
  stepLETStatus_ = STEP_LET_STATUS_UNKNOWN;
  alpha_s_ = -one;
  numberOfSteps_ = -1;
  nef_ = -1;
  nscsco_ = -1;
  newtonConvergenceStatus_ = -1;
}


template<class Scalar>
void ImplicitBDFStepper<Scalar>::obtainPredictor_()
{

  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  if (!isInitialized_) {
    return;
  }

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"obtainPredictor_");
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "currentOrder_ = " << currentOrder_ << std::endl;
  }
  
  // prepare history array for prediction
  for (int i=nscsco_;i<=currentOrder_;++i) {
    Vt_S(&*xHistory_[i],beta_[i]);
  }
  
  // evaluate predictor
  V_V(&*xn0_,*xHistory_[0]);
  V_S(&*xpn0_,ST::zero());
  for (int i=1;i<=currentOrder_;++i) {
    Vp_V(&*xn0_,*xHistory_[i]);
    Vp_StV(&*xpn0_,gamma_[i],*xHistory_[i]);
  }
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "xn0_ = " << std::endl;
    xn0_->describe(*out,verbLevel);
    *out << "xpn0_ = " << std::endl;
    xpn0_->describe(*out,verbLevel);
  }
}


template<class Scalar>
void ImplicitBDFStepper<Scalar>::interpolateSolution_(
  const Scalar& timepoint,
  Thyra::VectorBase<Scalar>* x_ptr,
  Thyra::VectorBase<Scalar>* xdot_ptr,
  ScalarMag* accuracy_ptr
  ) const
{

  typedef std::numeric_limits<Scalar> NL;
  typedef Teuchos::ScalarTraits<Scalar> ST;

#ifdef RYTHMOS_DEBUG
  TEST_FOR_EXCEPTION(
    !isInitialized_,std::logic_error,
    "Error, attempting to call interpolateSolution before initialization!\n");
  const TimeRange<Scalar> currTimeRange = this->getTimeRange();
  TEST_FOR_EXCEPTION(
    !currTimeRange.isInRange(timepoint), std::logic_error,
    "Error, timepoint = " << timepoint << " is not in the time range "
    << currTimeRange << "!" );
#endif

  const Scalar tn = time_;
  const int kused = usedOrder_;

  // order of interpolation
  int kord = kused;
  if ( (kused == 0) || (timepoint == tn) )  {
    kord = 1;
  }

  // Initialize interploation
  Thyra::V_V(x_ptr,*xHistory_[0]);
  Thyra::V_S(xdot_ptr,ST::zero());
  
  // Add history array contributions
  const Scalar delt = timepoint - tn;
  Scalar c = ST::one(); // coefficient for interpolation of x
  Scalar d = ST::zero(); // coefficient for interpolation of xdot
  Scalar gam = delt/psi_[0]; // coefficient for interpolation
  for (int j=1 ; j <= kord ; ++j) {
    d = d*gam + c/psi_[j-1];
    c = c*gam;
    gam = (delt + psi_[j-1])/psi_[j];
    Thyra::Vp_StV(x_ptr,c,*xHistory_[j]);
    Thyra::Vp_StV(xdot_ptr,d,*xHistory_[j]);
  }

  // Set approximate accuracy
  if (accuracy_ptr) {
    *accuracy_ptr = pow(usedStep_,kord);
  }
  
}


template<class Scalar>
void ImplicitBDFStepper<Scalar>::updateHistory_()
{

  using Teuchos::as;

  // Save Newton correction for potential order increase on next step.
  if (usedOrder_ < maxOrder_)  {
    assign( &*xHistory_[usedOrder_+1], *ee_ );
  }
  // Update history arrays
  Vp_V( &*xHistory_[usedOrder_], *ee_ );
  for (int j=usedOrder_-1;j>=0;j--) {
    Vp_V( &*xHistory_[j], *xHistory_[j+1] );
  }
  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"updateHistory_");
  if (as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    for (int i=0;i<std::max(2,maxOrder_);++i) {
      *out << "xHistory_[" << i << "] = " << std::endl;
      xHistory_[i]->describe(*out,verbLevel);
    }
  }

}


template<class Scalar>
void ImplicitBDFStepper<Scalar>::restoreHistory_()
{

  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  // undo preparation of history array for prediction
  for (int i=nscsco_;i<=currentOrder_;++i) {
    Vt_S( &*xHistory_[i], ST::one()/beta_[i] );
  }
  for (int i=1;i<=currentOrder_;++i) {
    psi_[i-1] = psi_[i] - hh_;
  }
  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"restoreHistory_");
  if (as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    for (int i=0;i<maxOrder_;++i) {
      *out << "psi_[" << i << "] = " << psi_[i] << std::endl;
    }
    for (int i=0;i<maxOrder_;++i) {
      *out << "xHistory_[" << i << "] = " << std::endl;
      xHistory_[i]->describe(*out,verbLevel);
    }
  }

} 


template<class Scalar>
void ImplicitBDFStepper<Scalar>::updateCoeffs_()
{

  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  // If the number of steps taken with constant order and constant stepsize is
  // more than the current order + 1 then we don't bother to update the
  // coefficients because we've reached a constant step-size formula.  When
  // this is is not true, then we update the coefficients for the variable
  // step-sizes. 
  if ((hh_ != usedStep_) || (currentOrder_ != usedOrder_)) {
    nscsco_ = 0;
  }
  nscsco_ = std::min(nscsco_+1,usedOrder_+2);
  if (currentOrder_+1 >= nscsco_) {
    beta_[0] = ST::one();
    alpha_[0] = ST::one();
    Scalar temp1 = hh_;
    gamma_[0] = ST::zero();
    for (int i=1;i<=currentOrder_;++i) {
      Scalar temp2 = psi_[i-1];
      psi_[i-1] = temp1;
      beta_[i] = beta_[i-1]*psi_[i-1]/temp2;
      temp1 = temp2 + hh_;
      alpha_[i] = hh_/temp1;
      gamma_[i] = gamma_[i-1]+alpha_[i-1]/hh_;
    }
    psi_[currentOrder_] = temp1;
  }
  alpha_s_ = ST::zero();
  for (int i=0;i<currentOrder_;++i) {
    alpha_s_ = alpha_s_ - Scalar(ST::one()/(i+ST::one()));
  }
  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"updateCoeffs_");
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    for (int i=0;i<=maxOrder_;++i) {
      *out << "alpha_[" << i << "] = " << alpha_[i] << std::endl;
      *out << "beta_[" << i << "] = " << beta_[i] << std::endl;
      *out << "gamma_[" << i << "] = " << gamma_[i] << std::endl;
      *out << "psi_[" << i << "] = " << psi_[i] << std::endl;
      *out << "alpha_s_ = " << alpha_s_ << std::endl;
    }
  }
}


template<class Scalar>
void ImplicitBDFStepper<Scalar>::initialize_()
{

  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  using Thyra::createMember;

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  const bool doTrace = (as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH));
  Teuchos::OSTab ostab(out,1,"initialize_");

  if (doTrace) {
    *out
      << "\nEntering " << this->Teuchos::Describable::description()
      << "::initialize_()...\n";
  }

  TEST_FOR_EXCEPT(model_ == Teuchos::null);
  TEST_FOR_EXCEPT(solver_ == Teuchos::null);
  TEUCHOS_ASSERT(haveInitialCondition_);

  // Initialize Parameter List if none provided.
  if (parameterList_ == Teuchos::null) {
    RCP<Teuchos::ParameterList> emptyParameterList = Teuchos::rcp(new Teuchos::ParameterList);
    this->setParameterList(emptyParameterList);
  }

  // Initialize StepControl
  if (stepControl_ == Teuchos::null) {
    RCP<ImplicitBDFStepperStepControl<Scalar> > implicitBDFStepperStepControl =
      Teuchos::rcp(new ImplicitBDFStepperStepControl<Scalar>());
    RCP<Teuchos::ParameterList> stepControlPL = 
      Teuchos::sublist(parameterList_, RythmosStepControlSettings_name);
    implicitBDFStepperStepControl->setParameterList(stepControlPL);
    this->setStepControlStrategy(implicitBDFStepperStepControl);
  }
  stepControl_->initialize(*this);

  maxOrder_ = stepControl_->getMaxOrder(); // maximum order
  TEST_FOR_EXCEPTION(
      !((1 <= maxOrder_) && (maxOrder_ <= 5)), std::logic_error,
      "Error, maxOrder returned from stepControl_->getMaxOrder() = " << maxOrder_ << " is outside range of [1,5]!\n"
      );

  Scalar zero = ST::zero();

  currentOrder_ = 1; // Current order of integration
  usedOrder_ = 1;  // order used in current step (used after currentOrder_ is updated)
  usedStep_ = zero;
  nscsco_  =  0;
  LETvalue_ = zero;

  alpha_.clear();  // $\alpha_j(n)=h_n/\psi_j(n)$ coefficient used in local error test
  // note:   $h_n$ = current step size, n = current time step
  gamma_.clear();  // calculate time derivative of history array for predictor 
  beta_.clear();   // coefficients used to evaluate predictor from history array
  psi_.clear();    // $\psi_j(n) = t_n-t_{n-j}$ intermediary variable used to 
  // compute $\beta_j(n):$
  for (int i=0 ; i<=maxOrder_ ; ++i) {
    alpha_.push_back(zero);
    beta_.push_back(zero);
    gamma_.push_back(zero);
    psi_.push_back(zero);
  }
  alpha_s_=Scalar(-ST::one());  // $\alpha_s$ fixed-leading coefficient of this BDF method
  hh_=zero;
  numberOfSteps_=0;   // number of total time integration steps taken
  nef_ = 0;

  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "alpha_s_ = " << alpha_s_ << std::endl;
    for (int i=0 ; i<=maxOrder_ ; ++i) {
      *out << "alpha_[" << i << "] = " << alpha_[i] << std::endl;
      *out << "beta_[" << i << "] = " << beta_[i] << std::endl;
      *out << "gamma_[" << i << "] = " << gamma_[i] << std::endl;
      *out << "psi_[" << i << "] = " << psi_[i] << std::endl;
    }
    *out << "numberOfSteps_ = " << numberOfSteps_ << std::endl;
  }

  // setInitialCondition initialized xHistory with xn0, xpn0.  Now we add the rest of the vectors.
  // Store maxOrder_+1 vectors
  for (int i=2 ; i<=maxOrder_ ; ++i) {
    xHistory_.push_back(createMember(xn0_->space())); 
    V_S(&*xHistory_[i],zero);
  }

  isInitialized_ = true;

  if (doTrace) {
    *out
      << "\nLeaving " << this->Teuchos::Describable::description()
      << "::initialize_()...\n";
  }

}


template<class Scalar>
void ImplicitBDFStepper<Scalar>::completeStep_()
{

  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;

#ifdef RYTHMOS_DEBUG
  TEST_FOR_EXCEPT(ST::isnaninf(hh_));
#endif  

  numberOfSteps_ ++;
  nef_ = 0;
  usedStep_ = hh_;
  usedOrder_ = currentOrder_;
  time_ += hh_;
  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"completeStep_");

  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "numberOfSteps_ = " << numberOfSteps_ << std::endl;
    *out << "time_ = " << time_ << std::endl;
  }
  
  // 03/22/04 tscoffe:  Note that updating the history has nothing to do with
  // the step-size and everything to do with the newton correction vectors.
  updateHistory_();

}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::setStepControlData(const StepperBase<Scalar> & stepper)
{
  if (!isInitialized_) {
    initialize_();
  }
  stepControl_->setStepControlData(stepper);
}

// 
// Explicit Instantiation macro
//
// Must be expanded from within the Rythmos namespace!
//

#define RYTHMOS_IMPLICITBDF_STEPPER_INSTANT(SCALAR) \
  \
  template class ImplicitBDFStepper< SCALAR >; \
  \
  template RCP< ImplicitBDFStepper< SCALAR > > \
  implicitBDFStepper();  \
  \
  template RCP< ImplicitBDFStepper< SCALAR > > \
  implicitBDFStepper( \
    const RCP<Thyra::ModelEvaluator< SCALAR > >& model, \
    const RCP<Thyra::NonlinearSolverBase< SCALAR > >& solver, \
    const RCP<Teuchos::ParameterList>& parameterList \
    ); \


} // namespace Rythmos


#endif //Rythmos_IMPLICITBDF_STEPPER_DEF_H
