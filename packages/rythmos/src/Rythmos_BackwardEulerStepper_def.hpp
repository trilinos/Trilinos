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

#ifndef Rythmos_BACKWARD_EULER_STEPPER_DEF_H
#define Rythmos_BACKWARD_EULER_STEPPER_DEF_H

#include "Rythmos_BackwardEulerStepper_decl.hpp"
#include "Rythmos_DataStore.hpp"
#include "Rythmos_LinearInterpolator.hpp"
#include "Rythmos_InterpolatorBaseHelpers.hpp"
#include "Rythmos_StepperHelpers.hpp"

#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_AssertOp.hpp"
#include "Thyra_TestingTools.hpp"

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_as.hpp"

namespace Rythmos {


template<class Scalar>
RCP<BackwardEulerStepper<Scalar> >
backwardEulerStepper(
  const RCP<Thyra::ModelEvaluator<Scalar> >& model,
  const RCP<Thyra::NonlinearSolverBase<Scalar> >& solver
  )
{
  return Teuchos::rcp(new BackwardEulerStepper<Scalar>(model, solver));
}

template<class Scalar>
RCP<BackwardEulerStepper<Scalar> >
backwardEulerStepper()
{
  return Teuchos::rcp(new BackwardEulerStepper<Scalar>());
}


// ////////////////////////////
// Defintions


// Constructors, intializers, Misc.


template<class Scalar>
BackwardEulerStepper<Scalar>::BackwardEulerStepper()
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  this->defaultInitializeAll_();
  numSteps_ = 0;
  dt_ = ST::zero();
}


template<class Scalar>
BackwardEulerStepper<Scalar>::BackwardEulerStepper(
  const RCP<Thyra::ModelEvaluator<Scalar> >& model,
  const RCP<Thyra::NonlinearSolverBase<Scalar> >& solver
  )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  this->defaultInitializeAll_();
  dt_ = ST::zero();
  numSteps_ = 0;
  setModel(model);
  setSolver(solver);
}

template<class Scalar>
void BackwardEulerStepper<Scalar>::defaultInitializeAll_()
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  isInitialized_ = false;
  haveInitialCondition_ = false;
  model_ = Teuchos::null;
  solver_ = Teuchos::null;
  scaled_x_old_ = Teuchos::null;
  x_dot_old_ = Teuchos::null;
  // basePoint_;
  x_ = Teuchos::null;
  x_dot_ = Teuchos::null;
  t_ = ST::nan();
  t_old_ = ST::nan();
  dt_ = ST::nan();
  numSteps_ = -1;
  neModel_ = Teuchos::null;
  parameterList_ = Teuchos::null;
  interpolator_ = Teuchos::null;
}

// Overridden from InterpolatorAcceptingObjectBase

template<class Scalar>
void BackwardEulerStepper<Scalar>::setInterpolator(
  const RCP<InterpolatorBase<Scalar> >& interpolator
  )
{
#ifdef RYTHMOS_DEBUG
  TEST_FOR_EXCEPT(is_null(interpolator));
#endif
  interpolator_ = interpolator;
  isInitialized_ = false;
}

template<class Scalar>
RCP<InterpolatorBase<Scalar> >
BackwardEulerStepper<Scalar>::getNonconstInterpolator()
{
  return interpolator_;
}

template<class Scalar>
RCP<const InterpolatorBase<Scalar> >
BackwardEulerStepper<Scalar>::getInterpolator() const
{
  return interpolator_;
}

template<class Scalar>
RCP<InterpolatorBase<Scalar> >
BackwardEulerStepper<Scalar>::unSetInterpolator()
{
  RCP<InterpolatorBase<Scalar> > temp_interpolator = interpolator_;
  interpolator_ = Teuchos::null;
  return(temp_interpolator);
  isInitialized_ = false;
}


// Overridden from SolverAcceptingStepperBase


template<class Scalar>
void BackwardEulerStepper<Scalar>::setSolver(
  const RCP<Thyra::NonlinearSolverBase<Scalar> > &solver
  )
{
  using Teuchos::as;

  TEST_FOR_EXCEPTION(solver == Teuchos::null, std::logic_error,
      "Error!  Thyra::NonlinearSolverBase RCP passed in through BackwardEulerStepper::setSolver is null!"
      );

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"BES::setSolver");
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "solver = " << solver->description() << std::endl;
  }

  solver_ = solver;

  isInitialized_ = false;

}


template<class Scalar>
RCP<Thyra::NonlinearSolverBase<Scalar> >
BackwardEulerStepper<Scalar>::getNonconstSolver()
{
  return solver_;
}


template<class Scalar>
RCP<const Thyra::NonlinearSolverBase<Scalar> >
BackwardEulerStepper<Scalar>::getSolver() const
{
  return solver_;
}


// Overridden from StepperBase
 

template<class Scalar>
bool BackwardEulerStepper<Scalar>::supportsCloning() const
{
  return true;
}


template<class Scalar>
RCP<StepperBase<Scalar> >
BackwardEulerStepper<Scalar>::cloneStepperAlgorithm() const
{
  RCP<BackwardEulerStepper<Scalar> >
    stepper = Teuchos::rcp(new BackwardEulerStepper<Scalar>);
  stepper->isInitialized_ = isInitialized_;
  stepper->model_ = model_; // Model is stateless so shallow copy is okay!
  if (!is_null(solver_))
    stepper->solver_ = solver_->cloneNonlinearSolver().assert_not_null();
  if (!is_null(x_))
    stepper->x_ = x_->clone_v().assert_not_null();
  if (!is_null(scaled_x_old_))
    stepper->scaled_x_old_ = scaled_x_old_->clone_v().assert_not_null();
  stepper->t_ = t_;
  stepper->t_old_ = t_old_;
  stepper->dt_ = dt_;
  stepper->numSteps_ = numSteps_;
  if (!is_null(neModel_))
    stepper->neModel_
    = Teuchos::rcp(new Rythmos::SingleResidualModelEvaluator<Scalar>);
  if (!is_null(parameterList_))
    stepper->parameterList_ = parameterList(*parameterList_);
  if (!is_null(interpolator_))
    stepper->interpolator_
      = interpolator_->cloneInterpolator().assert_not_null(); // ToDo: Implement cloneInterpolator()
  return stepper;
}

template<class Scalar>
bool BackwardEulerStepper<Scalar>::isImplicit() const
{
  return true;
}

template<class Scalar>
void BackwardEulerStepper<Scalar>::setModel(
  const RCP<const Thyra::ModelEvaluator<Scalar> >& model
  )
{

  using Teuchos::as;

  TEST_FOR_EXCEPT( is_null(model) );
  assertValidModel( *this, *model );

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"BES::setModel");
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "model = " << model->description() << std::endl;
  }
  model_ = model;

  // Wipe out x.  This will either be set thorugh setInitialCondition(...) or
  // it will be taken from the model's nominal vlaues!
//  x_ = Teuchos::null;
//  scaled_x_old_ = Teuchos::null;
//  x_dot_ = Teuchos::null;
//  x_dot_old_ = Teuchos::null;

//  isInitialized_ = false;
//  haveInitialCondition_ = setDefaultInitialConditionFromNominalValues<Scalar>(
//    *model_, Teuchos::ptr(this) );
  
}

template<class Scalar>
void BackwardEulerStepper<Scalar>::setNonconstModel(
  const RCP<Thyra::ModelEvaluator<Scalar> >& model
  )
{
  this->setModel(model); // TODO:  09/09/09 tscoffe:  Use ConstNonconstObjectContainer here!
}

template<class Scalar>
RCP<const Thyra::ModelEvaluator<Scalar> >
BackwardEulerStepper<Scalar>::getModel() const
{
  return model_;
}


template<class Scalar>
RCP<Thyra::ModelEvaluator<Scalar> >
BackwardEulerStepper<Scalar>::getNonconstModel() 
{
  return Teuchos::null;
}


template<class Scalar>
void BackwardEulerStepper<Scalar>::setInitialCondition(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
  )
{

  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Thyra::ModelEvaluatorBase MEB;

  basePoint_ = initialCondition;

  // x

  RCP<const Thyra::VectorBase<Scalar> >
    x_init = initialCondition.get_x();

#ifdef RYTHMOS_DEBUG
  TEST_FOR_EXCEPTION(
    is_null(x_init), std::logic_error,
    "Error, if the client passes in an intial condition to setInitialCondition(...),\n"
    "then x can not be null!" );
#endif

  x_ = x_init->clone_v();

  // x_dot

  RCP<const Thyra::VectorBase<Scalar> >
    x_dot_init = initialCondition.get_x_dot();

  if (!is_null(x_dot_init)) {
    x_dot_ = x_dot_init->clone_v();
  }
  else {
    x_dot_ = createMember(x_->space());
    assign(x_dot_.ptr(),ST::zero());
  }
  
  // t
  
  t_ = initialCondition.get_t();

  t_old_ = t_;

  // x_old 

  scaled_x_old_ = x_->clone_v();

  // x_dot_old
  
  x_dot_old_ = x_dot_->clone_v();


  haveInitialCondition_ = true;

}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
BackwardEulerStepper<Scalar>::getInitialCondition() const
{
  return basePoint_;
}


template<class Scalar>
Scalar BackwardEulerStepper<Scalar>::takeStep(Scalar dt, StepSizeType stepSizeType)
{

  using Teuchos::as;
  using Teuchos::incrVerbLevel;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Thyra::NonlinearSolverBase<Scalar> NSB;
  typedef Teuchos::VerboseObjectTempState<NSB> VOTSNSB;

  initialize();

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"BES::takeStep");
  VOTSNSB solver_outputTempState(solver_,out,incrVerbLevel(verbLevel,-1));

  if ( !is_null(out) && as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW) ) {
    *out
      << "\nEntering " << Teuchos::TypeNameTraits<BackwardEulerStepper<Scalar> >::name()
      << "::takeStep("<<dt<<","<<toString(stepSizeType)<<") ...\n"; 
  }

  dt_ = dt;

  if ((stepSizeType == STEP_TYPE_VARIABLE) || (dt == ST::zero())) {
    if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW) )
      *out << "\nThe arguments to takeStep are not valid for BackwardEulerStepper at this time." << std::endl;
    // print something out about this method not supporting automatic variable step-size
    return(Scalar(-ST::one()));
  }
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "\ndt = " << dt << std::endl;
  }


  //
  // Setup the nonlinear equations:
  //
  //   f( (1/dt)* x + (-1/dt)*x_old), x, t ) = 0
  //

  V_StV( scaled_x_old_.ptr(), Scalar(-ST::one()/dt), *x_ );
  t_old_ = t_;
  if(!neModel_.get()) {
    neModel_ = Teuchos::rcp(new Rythmos::SingleResidualModelEvaluator<Scalar>());
  }
  neModel_->initializeSingleResidualModel(
    model_, basePoint_,
    Scalar(ST::one()/dt), scaled_x_old_,
    ST::one(), Teuchos::null,
    t_old_+dt,
    Teuchos::null
    );
  if( solver_->getModel().get() != neModel_.get() ) {
    solver_->setModel(neModel_);
  }
  // 2007/05/18: rabartl: ToDo: Above, set the stream and the verbosity level
  // on solver_ so that we an see what it is doing!

  //
  // Solve the implicit nonlinear system to a tolerance of ???
  //
  
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW) ) {
    *out << "\nSolving the implicit backward-Euler timestep equation ...\n";
  }

  Thyra::SolveStatus<Scalar>
    neSolveStatus = solver_->solve(&*x_);

  // In the above solve, on input *x_ is the old value of x for the previous
  // time step which is used as the initial guess for the solver.  On output,
  // *x_ is the converged timestep solution.
 
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW) ) {
    *out << "\nOutput status of nonlinear solve:\n" << neSolveStatus;
  }

  // 2007/05/18: rabartl: ToDo: Above, get the solve status from the above
  // solve and at least print warning message if the solve fails!  Actually,
  // you should most likely thrown an exception if the solve fails or return
  // false if appropriate

  //
  // Update the step
  //

  assign( x_dot_old_.ptr(), *x_dot_ );

  // x_dot = (1/dt)*x - (1/dt)*x_old 
  V_StV( x_dot_.ptr(), Scalar(ST::one()/dt), *x_ );
  Vp_StV( x_dot_.ptr(), Scalar(ST::one()), *scaled_x_old_ );

  t_ += dt;

  numSteps_++;

  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "\nt_old_ = " << t_old_ << std::endl;
    *out << "\nt_ = " << t_ << std::endl;
  }

#ifdef RYTHMOS_DEBUG
  // 04/14/09 tscoffe: This code should be moved to StepperValidator

  if ( includesVerbLevel(verbLevel,Teuchos::VERB_LOW) )
    *out << "\nChecking to make sure that solution and the interpolated solution are the same! ...\n";

  {

    typedef ScalarTraits<Scalar> ST;
    typedef ScalarTraits<ScalarMag> SMT;
    
    Teuchos::OSTab tab(out);

    const StepStatus<Scalar> stepStatus = this->getStepStatus();

    RCP<const Thyra::VectorBase<Scalar> >
      x = stepStatus.solution,
      xdot = stepStatus.solutionDot;

    Array<Scalar> time_vec = Teuchos::tuple(stepStatus.time);
    Array<RCP<const Thyra::VectorBase<Scalar> > > x_vec, xdot_vec;
    this->getPoints(time_vec,&x_vec,&xdot_vec,0);

    RCP<const Thyra::VectorBase<Scalar> >
      x_interp = x_vec[0],
      xdot_interp = xdot_vec[0];

    TEST_FOR_EXCEPT(
      !Thyra::testRelNormDiffErr(
        "x", *x, "x_interp", *x_interp,
        "2*epsilon", ScalarMag(100.0*SMT::eps()),
        "2*epsilon", ScalarMag(100.0*SMT::eps()),
        includesVerbLevel(verbLevel,Teuchos::VERB_HIGH) ? out.get() : 0
        )
      );

    TEST_FOR_EXCEPT(
      !Thyra::testRelNormDiffErr(
        "xdot", *xdot, "xdot_interp", *xdot_interp,
        "2*epsilon", ScalarMag(100.0*SMT::eps()),
        "2*epsilon", ScalarMag(100.0*SMT::eps()),
        includesVerbLevel(verbLevel,Teuchos::VERB_HIGH) ? out.get() : 0
        )
      );

  }

  // 2007/07/25: rabartl: ToDo: Move the above test into a helper function so
  // that it can be used from lots of different places!

#endif // RYTHMOS_DEBUG

  if ( !is_null(out) && as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW) ) {
    *out
      << "\nLeaving " << Teuchos::TypeNameTraits<BackwardEulerStepper<Scalar> >::name()
      << "::takeStep(...) ...\n"; 
  }

  return(dt);

}


template<class Scalar>
const StepStatus<Scalar> BackwardEulerStepper<Scalar>::getStepStatus() const
{

  typedef Teuchos::ScalarTraits<Scalar> ST;

  StepStatus<Scalar> stepStatus; // Defaults to unknown status

  if (!isInitialized_) {
    stepStatus.stepStatus = STEP_STATUS_UNINITIALIZED;
  }
  else if (numSteps_ > 0) {
    stepStatus.stepStatus = STEP_STATUS_CONVERGED; 
  }
  // else unknown

  stepStatus.stepSize = dt_;
  stepStatus.order = 1;
  stepStatus.time = t_;
  stepStatus.solution = x_;
  stepStatus.solutionDot = x_dot_;

  return(stepStatus);

}


// Overridden from InterpolationBufferBase


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
BackwardEulerStepper<Scalar>::get_x_space() const
{
  return ( !is_null(model_) ? model_->get_x_space() : Teuchos::null );
}


template<class Scalar>
void BackwardEulerStepper<Scalar>::addPoints(
  const Array<Scalar>& time_vec,
  const Array<RCP<const Thyra::VectorBase<Scalar> > >& x_vec,
  const Array<RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec
  )
{

  typedef Teuchos::ScalarTraits<Scalar> ST;
  using Teuchos::as;

  initialize();

#ifdef RYTHMOS_DEBUG
  TEST_FOR_EXCEPTION(
    time_vec.size() == 0, std::logic_error,
    "Error, addPoints called with an empty time_vec array!\n");
#endif // RYTHMOS_DEBUG

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"BES::setPoints");

  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "time_vec = " << std::endl;
    for (int i=0 ; i<Teuchos::as<int>(time_vec.size()) ; ++i) {
      *out << "time_vec[" << i << "] = " << time_vec[i] << std::endl;
    }
  }
  else if (time_vec.size() == 1) {
    int n = 0;
    t_ = time_vec[n];
    t_old_ = t_;
    Thyra::V_V(x_.ptr(),*x_vec[n]);
    Thyra::V_V(scaled_x_old_.ptr(),*x_);
  }
  else {
    int n = time_vec.size()-1;
    int nm1 = time_vec.size()-2;
    t_ = time_vec[n];
    t_old_ = time_vec[nm1];
    Thyra::V_V(x_.ptr(),*x_vec[n]);
    Scalar dt = t_ - t_old_;
    Thyra::V_StV(scaled_x_old_.ptr(),Scalar(-ST::one()/dt),*x_vec[nm1]);
  }
}


template<class Scalar>
TimeRange<Scalar> BackwardEulerStepper<Scalar>::getTimeRange() const
{
  if ( !isInitialized_ && haveInitialCondition_ )
    return timeRange<Scalar>(t_,t_);
  if ( !isInitialized_ && !haveInitialCondition_ )
    return invalidTimeRange<Scalar>();
  return timeRange<Scalar>(t_old_,t_);
}


template<class Scalar>
void BackwardEulerStepper<Scalar>::getPoints(
  const Array<Scalar>& time_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
  Array<ScalarMag>* accuracy_vec
  ) const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  using Teuchos::constOptInArg;
  using Teuchos::ptr;
#ifdef RYTHMOS_DEBUG
  TEUCHOS_ASSERT(haveInitialCondition_);
#endif // RYTHMOS_DEBUG
  RCP<Thyra::VectorBase<Scalar> > x_temp = x_;
  if (compareTimeValues(t_old_,t_)!=0) {
    Scalar dt = t_ - t_old_;
    x_temp = scaled_x_old_->clone_v();
    Thyra::Vt_S(x_temp.ptr(),Scalar(-ST::one()*dt));  // undo the scaling
  }
  defaultGetPoints<Scalar>(
      t_old_, constOptInArg(*x_temp), constOptInArg(*x_dot_old_),
      t_, constOptInArg(*x_), constOptInArg(*x_dot_),
      time_vec, ptr(x_vec), ptr(xdot_vec), ptr(accuracy_vec),
      ptr(interpolator_.get())
      );

  /*
  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef Teuchos::ScalarTraits<ScalarMag> SMT;
  typename DataStore<Scalar>::DataStoreVector_t ds_nodes;
  typename DataStore<Scalar>::DataStoreVector_t ds_out;

#ifdef RYTHMOS_DEBUG
  TEST_FOR_EXCEPT(!haveInitialCondition_);
  TEST_FOR_EXCEPT( 0 == x_vec );
#endif

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"BES::getPoints");
  if ( !is_null(out) && as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW) ) {
    *out
      << "\nEntering " << Teuchos::TypeNameTraits<BackwardEulerStepper<Scalar> >::name()
      << "::getPoints(...) ...\n"; 
  }
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    for (int i=0 ; i<Teuchos::as<int>(time_vec.size()) ; ++i) {
      *out << "time_vec[" << i << "] = " << time_vec[i] << std::endl;
    }
    *out << "I can interpolate in the interval [" << t_old_ << "," << t_ << "]." << std::endl;
  }
  if (t_old_ != t_) {
    if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
      *out << "Passing two points to interpolator:  " << t_old_ << " and " << t_ << std::endl;
    }
    DataStore<Scalar> ds_temp;
    Scalar dt = t_ - t_old_;
#ifdef RYTHMOS_DEBUG
    TEST_FOR_EXCEPT(
      !Thyra::testRelErr(
        "dt", dt, "dt_", dt_,
        "1e+4*epsilon", ScalarMag(1e+4*SMT::eps()),
        "1e+2*epsilon", ScalarMag(1e+2*SMT::eps()),
        as<int>(verbLevel) >= as<int>(Teuchos::VERB_MEDIUM) ? out.get() : 0
        )
      );
#endif
    RCP<Thyra::VectorBase<Scalar> >
      x_temp = scaled_x_old_->clone_v();
    Thyra::Vt_S(&*x_temp,Scalar(-ST::one()*dt));  // undo the scaling
    ds_temp.time = t_old_;
    ds_temp.x = x_temp;
    ds_temp.xdot = x_dot_old_;
    ds_temp.accuracy = ScalarMag(dt);
    ds_nodes.push_back(ds_temp);
    ds_temp.time = t_;
    ds_temp.x = x_;
    ds_temp.xdot = x_dot_;
    ds_temp.accuracy = ScalarMag(dt);
    ds_nodes.push_back(ds_temp);
  }
  else {
    if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
      *out << "Passing one point to interpolator:  " << t_ << std::endl;
    }
    DataStore<Scalar> ds_temp;
    ds_temp.time = t_;
    ds_temp.x = x_;
    ds_temp.xdot = x_dot_;
    ds_temp.accuracy = ScalarMag(ST::zero());
    ds_nodes.push_back(ds_temp);
  }
  interpolate<Scalar>(*interpolator_,rcp(&ds_nodes,false),time_vec,&ds_out);
  Array<Scalar> time_out;
  dataStoreVectorToVector(ds_out,&time_out,x_vec,xdot_vec,accuracy_vec);
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "Passing out the interpolated values:" << std::endl;
    for (int i=0; i<Teuchos::as<int>(time_out.size()) ; ++i) {
      *out << "time[" << i << "] = " << time_out[i] << std::endl;
      *out << "x_vec[" << i << "] = " << std::endl;
      (*x_vec)[i]->describe(*out,Teuchos::VERB_EXTREME);
      if (xdot_vec) {
        if ( (*xdot_vec)[i] == Teuchos::null) {
          *out << "xdot_vec[" << i << "] = Teuchos::null" << std::endl;
        }
        else {
          *out << "xdot_vec[" << i << "] = " << std::endl;
          (*xdot_vec)[i]->describe(*out,Teuchos::VERB_EXTREME);
        }
      }
      if(accuracy_vec) {
        *out << "accuracy[" << i << "] = " << (*accuracy_vec)[i] << std::endl;
      }
    }
  }
  if ( !is_null(out) && as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW) ) {
    *out
      << "Leaving " << Teuchos::TypeNameTraits<BackwardEulerStepper<Scalar> >::name()
      << "::getPoints(...) ...\n"; 
  }
  */

}


template<class Scalar>
void BackwardEulerStepper<Scalar>::getNodes(Array<Scalar>* time_vec) const
{
  using Teuchos::as;

  TEUCHOS_ASSERT( time_vec != NULL );

  time_vec->clear();

  if (!haveInitialCondition_) {
    return;
  }

  time_vec->push_back(t_old_);

  if (numSteps_ > 0) {
    time_vec->push_back(t_);
  }

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"BES::getNodes");
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << this->description() << std::endl;
    for (int i=0 ; i<Teuchos::as<int>(time_vec->size()) ; ++i) {
      *out << "time_vec[" << i << "] = " << (*time_vec)[i] << std::endl;
    }
  }
}


template<class Scalar>
void BackwardEulerStepper<Scalar>::removeNodes(Array<Scalar>& time_vec) 
{
  initialize();
  using Teuchos::as;
  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"BES::removeNodes");
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "time_vec = " << std::endl;
    for (int i=0 ; i<Teuchos::as<int>(time_vec.size()) ; ++i) {
      *out << "time_vec[" << i << "] = " << time_vec[i] << std::endl;
    }
  }
  TEST_FOR_EXCEPTION(true,std::logic_error,"Error, removeNodes is not implemented for BackwardEulerStepper at this time.\n");
  // TODO:
  // if any time in time_vec matches t_ or t_old_, then do the following:
  // remove t_old_:  set t_old_ = t_ and set scaled_x_old_ = x_
  // remove t_:  set t_ = t_old_ and set x_ = -dt*scaled_x_old_
}


template<class Scalar>
int BackwardEulerStepper<Scalar>::getOrder() const
{
  return(1);
}


// Overridden from Teuchos::ParameterListAcceptor


template <class Scalar>
void BackwardEulerStepper<Scalar>::setParameterList(
  RCP<Teuchos::ParameterList> const& paramList
  )
{
  TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParametersAndSetDefaults(*this->getValidParameters());
  parameterList_ = paramList;
  Teuchos::readVerboseObjectSublist(&*parameterList_,this);
}


template <class Scalar>
RCP<Teuchos::ParameterList>
BackwardEulerStepper<Scalar>::getNonconstParameterList()
{
  return(parameterList_);
}


template <class Scalar>
RCP<Teuchos::ParameterList>
BackwardEulerStepper<Scalar>::unsetParameterList()
{
  RCP<Teuchos::ParameterList>
    temp_param_list = parameterList_;
  parameterList_ = Teuchos::null;
  return(temp_param_list);
}


template<class Scalar>
RCP<const Teuchos::ParameterList>
BackwardEulerStepper<Scalar>::getValidParameters() const
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


// Overridden from Teuchos::Describable


template<class Scalar>
void BackwardEulerStepper<Scalar>::describe(
  Teuchos::FancyOStream &out,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  using Teuchos::as;
  Teuchos::OSTab tab(out);
  if (!isInitialized_) {
    out << this->description() << " : This stepper is not initialized yet" << std::endl;
    return;
  }
  if (
    as<int>(verbLevel) == as<int>(Teuchos::VERB_DEFAULT)
    ||
    as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW)
    )
  {
    out << this->description() << "::describe:" << std::endl;
    out << "model = " << model_->description() << std::endl;
    out << "solver = " << solver_->description() << std::endl;
    if (neModel_ == Teuchos::null) {
      out << "neModel = Teuchos::null" << std::endl;
    } else {
      out << "neModel = " << neModel_->description() << std::endl;
    }
  }
  else if (as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW)) {
    out << "t_ = " << t_ << std::endl;
    out << "t_old_ = " << t_old_ << std::endl;
  }
  else if (as<int>(verbLevel) >= as<int>(Teuchos::VERB_MEDIUM)) {
  }
  else if (as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH)) {
    out << "model_ = " << std::endl;
    model_->describe(out,verbLevel);
    out << "solver_ = " << std::endl;
    solver_->describe(out,verbLevel);
    if (neModel_ == Teuchos::null) {
      out << "neModel = Teuchos::null" << std::endl;
    } else {
      out << "neModel = " << std::endl;
      neModel_->describe(out,verbLevel);
    }
    out << "x_ = " << std::endl;
    x_->describe(out,verbLevel);
    out << "scaled_x_old_ = " << std::endl;
    scaled_x_old_->describe(out,verbLevel);
  }
}


// private


template <class Scalar>
void BackwardEulerStepper<Scalar>::initialize()
{

  if (isInitialized_)
    return;

  TEST_FOR_EXCEPT(is_null(model_));
  TEST_FOR_EXCEPT(is_null(solver_));
  TEST_FOR_EXCEPT(!haveInitialCondition_);

#ifdef RYTHMOS_DEBUG
  THYRA_ASSERT_VEC_SPACES(
    "Rythmos::BackwardEulerStepper::initialize(...)",
    *x_->space(), *model_->get_x_space() );
#endif // RYTHMOS_DEBUG

  if ( is_null(interpolator_) ) {
    // If an interpolator has not been explicitly set, then just create
    // a default linear interpolator.
    interpolator_ = linearInterpolator<Scalar>();
    // 2007/05/18: rabartl: ToDo: Replace this with a Hermete interplator
    // when it is implementated!
  }

  isInitialized_ = true;

}

template<class Scalar>
RCP<const MomentoBase<Scalar> >
BackwardEulerStepper<Scalar>::getMomento() const
{
  RCP<BackwardEulerStepperMomento<Scalar> > momento = Teuchos::rcp(new BackwardEulerStepperMomento<Scalar>());
  momento->set_scaled_x_old(scaled_x_old_);
  momento->set_x_dot_old(x_dot_old_);
  momento->set_x(x_);
  momento->set_x_dot(x_dot_);
  momento->set_t(t_);
  momento->set_t_old(t_old_);
  momento->set_dt(dt_);
  momento->set_numSteps(numSteps_);
  momento->set_isInitialized(isInitialized_);
  momento->set_haveInitialCondition(haveInitialCondition_);
  momento->set_parameterList(parameterList_);
  momento->set_basePoint(basePoint_);
  momento->set_neModel(neModel_);
  momento->set_interpolator(interpolator_);
  return momento;
}

template<class Scalar>
void BackwardEulerStepper<Scalar>::setMomento(
    const Ptr<const MomentoBase<Scalar> >& momentoPtr,
    const RCP<Thyra::ModelEvaluator<Scalar> >& model,
    const RCP<Thyra::NonlinearSolverBase<Scalar> >& solver
    ) 
{ 
  Ptr<const BackwardEulerStepperMomento<Scalar> > feMomentoPtr = 
    Teuchos::ptr_dynamic_cast<const BackwardEulerStepperMomento<Scalar> >(momentoPtr,true);
  const BackwardEulerStepperMomento<Scalar>& feMomento = *feMomentoPtr;
  model_ = model;
  solver_ = solver;
  scaled_x_old_ = feMomento.get_scaled_x_old();
  x_dot_old_ = feMomento.get_x_dot_old();
  x_ = feMomento.get_x();
  x_dot_ = feMomento.get_x_dot();
  t_ = feMomento.get_t();
  t_old_ = feMomento.get_t_old();
  dt_ = feMomento.get_dt();
  numSteps_ = feMomento.get_numSteps();
  isInitialized_ = feMomento.get_isInitialized();
  haveInitialCondition_ = feMomento.get_haveInitialCondition();
  parameterList_ = feMomento.get_parameterList();
  basePoint_ = feMomento.get_basePoint();
  neModel_ = feMomento.get_neModel();
  interpolator_ = feMomento.get_interpolator();
  this->checkConsistentState_();
}

template<class Scalar>
void BackwardEulerStepper<Scalar>::checkConsistentState_()
{
  if (isInitialized_) {
    TEUCHOS_ASSERT( !Teuchos::is_null(model_) );
    TEUCHOS_ASSERT( !Teuchos::is_null(solver_) );
    TEUCHOS_ASSERT( haveInitialCondition_ );
    TEUCHOS_ASSERT( !Teuchos::is_null(interpolator_) );
  }
  if (haveInitialCondition_) {
    // basePoint_ should be defined
    typedef Teuchos::ScalarTraits<Scalar> ST;
    TEUCHOS_ASSERT( !ST::isnaninf(t_) );
    TEUCHOS_ASSERT( !ST::isnaninf(t_old_) );
    TEUCHOS_ASSERT( !Teuchos::is_null(scaled_x_old_) );
    TEUCHOS_ASSERT( !Teuchos::is_null(x_dot_old_) );
    TEUCHOS_ASSERT( !Teuchos::is_null(x_) );
    TEUCHOS_ASSERT( !Teuchos::is_null(x_dot_) );
    TEUCHOS_ASSERT( t_ >= basePoint_.get_t() );
    TEUCHOS_ASSERT( t_old_ >= basePoint_.get_t() );
  }
  if (numSteps_ > 0) {
    TEUCHOS_ASSERT(isInitialized_);
  }
}


// 
// Explicit Instantiation macro
//
// Must be expanded from within the Rythmos namespace!
//

#define RYTHMOS_BACKWARD_EULER_STEPPER_INSTANT(SCALAR) \
  \
  template class BackwardEulerStepper< SCALAR >; \
  \
  template RCP< BackwardEulerStepper< SCALAR > > \
  backwardEulerStepper( \
    const RCP<Thyra::ModelEvaluator< SCALAR > >& model, \
    const RCP<Thyra::NonlinearSolverBase< SCALAR > >& solver \
      ); \
  template RCP< BackwardEulerStepper< SCALAR > > \
  backwardEulerStepper(); 
   



} // namespace Rythmos

#endif //Rythmos_BACKWARD_EULER_STEPPER_DEF_H
