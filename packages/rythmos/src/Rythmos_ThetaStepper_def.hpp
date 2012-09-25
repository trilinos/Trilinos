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

#ifndef Rythmos_THETA_STEPPER_DEF_H
#define Rythmos_THETA_STEPPER_DEF_H

#include "Rythmos_ConfigDefs.h"
#ifdef HAVE_RYTHMOS_EXPERIMENTAL

#include "Rythmos_ThetaStepper_decl.hpp"

namespace Rythmos {


/** \brief Nonmember constructor.
 *
 * \relates ThetaStepper
 */
template<class Scalar>
RCP<ThetaStepper<Scalar> >
thetaStepper(
  const RCP<Thyra::ModelEvaluator<Scalar> >& model,
  const RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
  RCP<Teuchos::ParameterList>& parameterList
  )
{
  Teuchos::RCP<ThetaStepper<Scalar> > stepper = 
    Teuchos::rcp(new ThetaStepper<Scalar>());
  stepper->setParameterList(parameterList);
  stepper->setModel(model);
  stepper->setSolver(solver);

  return stepper;
}

// ////////////////////////////
// Defintions


// Constructors, intializers, Misc.


template<class Scalar>
ThetaStepper<Scalar>::ThetaStepper()
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  this->defaultInitializeAll_();
  Scalar zero = ST::zero();
  t_ = -ST::one();
  t_old_ = zero;
  dt_ = zero;
  dt_old_ = zero;
  numSteps_ = 0;
}

template<class Scalar>
void ThetaStepper<Scalar>::defaultInitializeAll_()
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  isInitialized_ = false;
  haveInitialCondition_ = false;
  model_ = Teuchos::null;
  solver_ = Teuchos::null;
  //basePoint_;
  x_ = Teuchos::null;
  x_old_ = Teuchos::null;
  x_pre_ = Teuchos::null;
  x_dot_ = Teuchos::null;
  x_dot_old_ = Teuchos::null;
  x_dot_really_old_ = Teuchos::null;
  x_dot_base_ = Teuchos::null;
  t_ = ST::nan();
  t_old_ = ST::nan();
  dt_ = ST::nan();
  dt_old_ = ST::nan();
  numSteps_ = -1;
  thetaStepperType_ = INVALID_THETA_STEPPER_TYPE;
  theta_ = ST::nan();
  neModel_ = Teuchos::null;
  parameterList_ = Teuchos::null;
  interpolator_ = Teuchos::null;
  predictor_corrector_begin_after_step_ = -1;
  default_predictor_order_ = -1;
}

template<class Scalar>
bool ThetaStepper<Scalar>::isImplicit() const
{
  return true;
}


template<class Scalar>
void ThetaStepper<Scalar>::setInterpolator(
  const RCP<InterpolatorBase<Scalar> >& interpolator
  )
{
#ifdef HAVE_RYTHMOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(is_null(interpolator));
#endif
  interpolator_ = interpolator;
  isInitialized_ = false;
}

template<class Scalar>
RCP<InterpolatorBase<Scalar> >
  ThetaStepper<Scalar>::getNonconstInterpolator()
{
  return interpolator_;
}

template<class Scalar>
RCP<const InterpolatorBase<Scalar> >
  ThetaStepper<Scalar>::getInterpolator() const
{
  return interpolator_;
}


template<class Scalar>
RCP<InterpolatorBase<Scalar> >
ThetaStepper<Scalar>::unSetInterpolator()
{
  RCP<InterpolatorBase<Scalar> > temp_interpolator = interpolator_;
  interpolator_ = Teuchos::null;
  return(temp_interpolator);
  isInitialized_ = false;
}


// Overridden from SolverAcceptingStepperBase


template<class Scalar>
void ThetaStepper<Scalar>::setSolver(
  const RCP<Thyra::NonlinearSolverBase<Scalar> > &solver
  )
{
  using Teuchos::as;

  TEUCHOS_TEST_FOR_EXCEPTION(solver == Teuchos::null, std::logic_error,
      "Error!  Thyra::NonlinearSolverBase RCP passed in through ThetaStepper::setSolver is null!"
      );

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"TS::setSolver");
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "solver = " << solver->description() << std::endl;
  }

  solver_ = solver;

  isInitialized_ = false;

}


template<class Scalar>
RCP<Thyra::NonlinearSolverBase<Scalar> >
ThetaStepper<Scalar>::getNonconstSolver()
{
  return solver_;
}


template<class Scalar>
RCP<const Thyra::NonlinearSolverBase<Scalar> >
ThetaStepper<Scalar>::getSolver() const
{
  return solver_;
}


// Overridden from StepperBase
 

template<class Scalar>
bool ThetaStepper<Scalar>::supportsCloning() const
{
  return true;
}

template<class Scalar>
RCP<StepperBase<Scalar> >
ThetaStepper<Scalar>::cloneStepperAlgorithm() const
{
  RCP<ThetaStepper<Scalar> >
    stepper = Teuchos::rcp(new ThetaStepper<Scalar>);
  stepper->isInitialized_ = isInitialized_;
  stepper->model_ = model_; // Model is stateless so shallow copy is okay!

  if (!is_null(solver_))
    stepper->solver_ = solver_->cloneNonlinearSolver().assert_not_null();

  stepper->basePoint_ = basePoint_;

  if (!is_null(x_))
    stepper->x_ = x_->clone_v().assert_not_null();
  if (!is_null(x_old_))
    stepper->x_old_ = x_old_->clone_v().assert_not_null();
  if (!is_null(x_pre_))
    stepper->x_pre_ = x_pre_->clone_v().assert_not_null();

  if (!is_null(x_dot_))
    stepper->x_dot_ = x_dot_->clone_v().assert_not_null();
  if (!is_null(x_dot_old_))
    stepper->x_dot_old_ = x_dot_old_->clone_v().assert_not_null();
  if (!is_null(x_dot_really_old_))
    stepper->x_dot_really_old_ = x_dot_really_old_->clone_v().assert_not_null();
  if (!is_null(x_dot_base_))
    stepper->x_dot_base_ = x_dot_base_->clone_v().assert_not_null();

  stepper->t_ = t_;
  stepper->t_old_ = t_old_;

  stepper->dt_ = dt_;
  stepper->dt_old_ = dt_old_;

  stepper->numSteps_ = numSteps_;

  stepper->thetaStepperType_ = thetaStepperType_;
  stepper->theta_ = theta_;
  stepper->predictor_corrector_begin_after_step_ = predictor_corrector_begin_after_step_;
  stepper->default_predictor_order_ = default_predictor_order_;

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
void ThetaStepper<Scalar>::setModel(
  const RCP<const Thyra::ModelEvaluator<Scalar> >& model
  )
{

  using Teuchos::as;

  TEUCHOS_TEST_FOR_EXCEPT( is_null(model) );
  assertValidModel( *this, *model );

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"TS::setModel");
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "model = " << model->description() << std::endl;
  }
  model_ = model;

  // Wipe out x.  This will either be set thorugh setInitialCondition(...) or
  // it will be taken from the model's nominal vlaues!
  x_ = Teuchos::null;
  x_old_ = Teuchos::null;
  x_pre_ = Teuchos::null;

  x_dot_ = Teuchos::null;
  x_dot_old_ = Teuchos::null;
  x_dot_really_old_ = Teuchos::null;
  x_dot_base_ = Teuchos::null;

  isInitialized_ = false;
  haveInitialCondition_ = setDefaultInitialConditionFromNominalValues<Scalar>(
    *model_, Teuchos::ptr(this) );
  
}

template<class Scalar>
void ThetaStepper<Scalar>::setNonconstModel(
  const RCP<Thyra::ModelEvaluator<Scalar> >& model
  )
{
  this->setModel(model); // TODO 09/09/09 tscoffe:  use ConstNonconstObjectContainer!
}


template<class Scalar>
RCP<const Thyra::ModelEvaluator<Scalar> >
ThetaStepper<Scalar>::getModel() const
{
  return model_;
}


template<class Scalar>
RCP<Thyra::ModelEvaluator<Scalar> >
ThetaStepper<Scalar>::getNonconstModel() 
{
  TEUCHOS_TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}


template<class Scalar>
void ThetaStepper<Scalar>::setInitialCondition(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
  )
{

  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Thyra::ModelEvaluatorBase MEB;

  basePoint_ = initialCondition;

  // x

  RCP<const Thyra::VectorBase<Scalar> >
    x_init = initialCondition.get_x();

#ifdef HAVE_RYTHMOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
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

  dt_old_ = 0.0;

  // x pre
  
  x_pre_ = x_->clone_v();

  // x old

  x_old_ = x_->clone_v();

  // x dot base

  x_dot_base_ = x_->clone_v();

  // x dot old

  x_dot_old_ = x_dot_->clone_v();

  // x dot really old

  x_dot_really_old_ = x_dot_->clone_v();

  haveInitialCondition_ = true;

}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar> 
ThetaStepper<Scalar>::getInitialCondition() const
{
  return basePoint_;
}


template<class Scalar>
Scalar ThetaStepper<Scalar>::takeStep(Scalar dt, StepSizeType stepSizeType)
{

  using Teuchos::as;
  using Teuchos::incrVerbLevel;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Thyra::NonlinearSolverBase<Scalar> NSB;
  typedef Teuchos::VerboseObjectTempState<NSB> VOTSNSB;

  initialize_();

  // DEBUG
  //this->setOverridingVerbLevel(Teuchos::VERB_EXTREME);
  //this->setVerbLevel(Teuchos::VERB_EXTREME);

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"TS::takeStep");
  VOTSNSB solver_outputTempState(solver_,out,incrVerbLevel(verbLevel,-1));

  if ( !is_null(out) && as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW) ) {
    *out
      << "\nEntering " << Teuchos::TypeNameTraits<ThetaStepper<Scalar> >::name()
      << "::takeStep("<<dt<<","<<toString(stepSizeType)<<") ...\n"; 
  }

  dt_ = dt;
  V_StV( x_old_.ptr(), Scalar(ST::one()), *x_ );
  V_StV( x_dot_really_old_.ptr(), Scalar(ST::one()), *x_dot_old_ );
  V_StV( x_dot_old_.ptr(),        Scalar(ST::one()), *x_dot_ );

  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_EXTREME) ) {
    *out << "\nSetting dt_ and old data ..." << std::endl;
    *out << "\ndt_ = " << dt_;
    *out << "\nx_old_ = " << *x_old_;
    *out << "\nx_dot_old_ = " << *x_dot_old_;
    *out << "\nx_dot_really_old_ = " << *x_dot_really_old_;
  }

  if ((stepSizeType == STEP_TYPE_VARIABLE) || (dt == ST::zero())) {
    if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW) )
      *out << "\nThe arguments to takeStep are not valid for ThetaStepper at this time." << std::endl;
    return(Scalar(-ST::one()));
  }
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "\ndt = " << dt << std::endl;
    *out << "\nstepSizeType = " << stepSizeType << std::endl;
  }

  // compute predictor
  obtainPredictor_();

  //
  // Setup the nonlinear equations:
  //
  //   substitute:
  // 
  //   x_dot = ( 1/(theta*dt) )*x + ( -1/(theta*dt) )*x_old + ( -(1-theta)/theta )*x_dot_old
  //   
  //   f( x_dot, x, t ) = 0
  //

  const double theta = (numSteps_+1>=predictor_corrector_begin_after_step_) ? theta_ : 1.0;

  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_EXTREME) ) {
    *out << "\nSetting x_dot_base_ ..." << std::endl;
    *out << "\ntheta = " << theta;
    *out << "\nx_ = " << *x_;
    *out << "\nx_dot_old_ = " << *x_dot_old_;
  }

  const Scalar coeff_x_dot = Scalar(ST::one()/(theta*dt)); 
  const Scalar coeff_x = ST::one();

  const Scalar x_coeff = Scalar(-coeff_x_dot);
  const Scalar x_dot_old_coeff = Scalar( -(ST::one()-theta)/theta);

  V_StV( x_dot_base_.ptr(), x_coeff, *x_old_ );
  Vp_StV( x_dot_base_.ptr(), x_dot_old_coeff, *x_dot_old_);

  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_EXTREME) ) {
    *out << "\nx_dot_base_ = " << *x_dot_base_;
  }

  if(!neModel_.get()) {
    neModel_ = Teuchos::rcp(new Rythmos::SingleResidualModelEvaluator<Scalar>());
  }

  neModel_->initializeSingleResidualModel(
    model_, 
    basePoint_,
    coeff_x_dot,
    x_dot_base_,
    coeff_x,
    Teuchos::null, // x_base
    t_+dt, // t_base
    Teuchos::null // x_bar_init
    );
  if( solver_->getModel().get() != neModel_.get() ) {
    solver_->setModel(neModel_);
  }
  
  solver_->setVerbLevel(this->getVerbLevel());

  //
  // Solve the implicit nonlinear system to a tolerance of ???
  //
  
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW) ) {
    *out << "\nSolving the implicit theta-stepper timestep equation"
	 << " with theta = " << theta << "\n";
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
  // Update the step data
  //

  // x_dot = ( 1/(theta*dt) )*x + ( -1/(theta*dt) )*x_old + ( -(1-theta)/theta )*x_dot_old

  V_StV(  x_dot_.ptr(), Scalar( ST::one()/(theta*dt)), *x_ );
  Vp_StV( x_dot_.ptr(), Scalar(-ST::one()/(theta*dt)), *x_old_ );
  Vp_StV( x_dot_.ptr(), Scalar( -(ST::one()-theta)/theta), *x_dot_old_ );

  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_EXTREME) ) {
    *out << "\nUpdating x_dot_ ...\n";
    *out << "\nx_dot_ = " << *x_dot_ << std::endl;
  }

  t_old_ = t_;
  dt_old_ = dt;
  t_ += dt;

  numSteps_++;

  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "\nt_old_ = " << t_old_ << std::endl;
    *out << "\nt_ = " << t_ << std::endl;
  }

#ifdef HAVE_RYTHMOS_DEBUG

  if ( includesVerbLevel(verbLevel,Teuchos::VERB_LOW) )
    *out << "\nChecking to make sure that solution and the interpolated solution are the same! ...\n";

  {

    typedef ScalarTraits<Scalar> ST;
    typedef typename ST::magnitudeType ScalarMag;
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

    TEUCHOS_TEST_FOR_EXCEPT(
      !Thyra::testRelNormDiffErr(
        "x", *x, "x_interp", *x_interp,
        "2*epsilon", ScalarMag(100.0*SMT::eps()),
        "2*epsilon", ScalarMag(100.0*SMT::eps()),
        includesVerbLevel(verbLevel,Teuchos::VERB_HIGH) ? out.get() : 0
        )
      );

    TEUCHOS_TEST_FOR_EXCEPT(
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

#endif // HAVE_RYTHMOS_DEBUG

  if ( !is_null(out) && as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW) ) {
    *out
      << "\nLeaving " << Teuchos::TypeNameTraits<ThetaStepper<Scalar> >::name()
      << "::takeStep(...) ...\n"; 
  }

  return(dt);

}


template<class Scalar>
const StepStatus<Scalar> ThetaStepper<Scalar>::getStepStatus() const
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
ThetaStepper<Scalar>::get_x_space() const
{
  return ( !is_null(model_) ? model_->get_x_space() : Teuchos::null );
}


template<class Scalar>
void ThetaStepper<Scalar>::addPoints(
  const Array<Scalar>& time_vec,
  const Array<RCP<const Thyra::VectorBase<Scalar> > >& x_vec,
  const Array<RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec
  )
{

  typedef Teuchos::ScalarTraits<Scalar> ST;
  using Teuchos::as;

  initialize_();

#ifdef HAVE_RYTHMOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    time_vec.size() == 0, std::logic_error,
    "Error, addPoints called with an empty time_vec array!\n");
#endif // HAVE_RYTHMOS_DEBUG

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"TS::setPoints");

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
    Thyra::V_V(x_dot_base_.ptr(),*x_);
  }
  else {
    int n = time_vec.size()-1;
    int nm1 = time_vec.size()-2;
    t_ = time_vec[n];
    t_old_ = time_vec[nm1];
    Thyra::V_V(x_.ptr(),*x_vec[n]);
    Scalar dt = t_ - t_old_;
    Thyra::V_StV(x_dot_base_.ptr(),Scalar(-ST::one()/dt),*x_vec[nm1]);
  }
}


template<class Scalar>
TimeRange<Scalar> ThetaStepper<Scalar>::getTimeRange() const
{
  if ( !isInitialized_ && haveInitialCondition_ )
    return timeRange<Scalar>(t_,t_);
  if ( !isInitialized_ && !haveInitialCondition_ )
    return invalidTimeRange<Scalar>();
  return timeRange<Scalar>(t_old_,t_);
}


template<class Scalar>
void ThetaStepper<Scalar>::getPoints(
  const Array<Scalar>& time_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
  Array<ScalarMag>* accuracy_vec
  ) const
{
  using Teuchos::constOptInArg;
  using Teuchos::ptr;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  TEUCHOS_ASSERT(haveInitialCondition_);

  RCP<Thyra::VectorBase<Scalar> > x_temp = x_;
  if (compareTimeValues(t_old_,t_)!=0) {
    Scalar dt = t_ - t_old_;
    x_temp = x_dot_base_->clone_v();
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

#ifdef HAVE_RYTHMOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(!haveInitialCondition_);
  TEUCHOS_TEST_FOR_EXCEPT( 0 == x_vec );
#endif

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"TS::getPoints");
  if ( !is_null(out) && as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW) ) {
    *out
      << "\nEntering " << Teuchos::TypeNameTraits<ThetaStepper<Scalar> >::name()
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
#ifdef HAVE_RYTHMOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPT(
      !Thyra::testRelErr(
        "dt", dt, "dt_", dt_,
        "1e+4*epsilon", ScalarMag(1e+4*SMT::eps()),
        "1e+2*epsilon", ScalarMag(1e+2*SMT::eps()),
        as<int>(verbLevel) >= as<int>(Teuchos::VERB_MEDIUM) ? out.get() : 0
        )
      );
#endif
    ds_temp.time = t_old_;
    ds_temp.x = x_old_;
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
      if (x_vec) {
        if ( (*x_vec)[i] == Teuchos::null) {
          *out << "x_vec[" << i << "] = Teuchos::null" << std::endl;
	}
	else {
	  *out << "time[" << i << "] = " << time_out[i] << std::endl;
	  *out << "x_vec[" << i << "] = " << std::endl;
	  (*x_vec)[i]->describe(*out,Teuchos::VERB_EXTREME);
	}
      }
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
      << "Leaving " << Teuchos::TypeNameTraits<ThetaStepper<Scalar> >::name()
      << "::getPoints(...) ...\n"; 
  }
  */

}


template<class Scalar>
void ThetaStepper<Scalar>::getNodes(Array<Scalar>* time_vec) const
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
  Teuchos::OSTab ostab(out,1,"TS::getNodes");
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << this->description() << std::endl;
    for (int i=0 ; i<Teuchos::as<int>(time_vec->size()) ; ++i) {
      *out << "time_vec[" << i << "] = " << (*time_vec)[i] << std::endl;
    }
  }
}


template<class Scalar>
void ThetaStepper<Scalar>::removeNodes(Array<Scalar>& time_vec) 
{
  initialize_();
  using Teuchos::as;
  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"TS::removeNodes");
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "time_vec = " << std::endl;
    for (int i=0 ; i<Teuchos::as<int>(time_vec.size()) ; ++i) {
      *out << "time_vec[" << i << "] = " << time_vec[i] << std::endl;
    }
  }
  TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Error, removeNodes is not implemented for ThetaStepper at this time.\n");
  // TODO:
  // if any time in time_vec matches t_ or t_old_, then do the following:
  // remove t_old_:  set t_old_ = t_ and set x_dot_base_ = x_
  // remove t_:  set t_ = t_old_ and set x_ = -dt*x_dot_base_
}


template<class Scalar>
int ThetaStepper<Scalar>::getOrder() const
{
  return (thetaStepperType_==ImplicitEuler) ? 1 : 2;
}


// Overridden from Teuchos::ParameterListAcceptor


template <class Scalar>
void ThetaStepper<Scalar>::setParameterList(
  RCP<Teuchos::ParameterList> const& paramList
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParametersAndSetDefaults(*this->getValidParameters());
  parameterList_ = paramList;
  Teuchos::readVerboseObjectSublist(&*parameterList_,this);

  RCP<ParameterList> pl_theta = Teuchos::sublist(parameterList_, RythmosStepControlSettings_name);

  std::string thetaStepperTypeString = 
    Teuchos::getParameter<std::string>(*pl_theta, ThetaStepperType_name);
  
  if (thetaStepperTypeString == "Implicit Euler")
    thetaStepperType_ = ImplicitEuler;
  else if (thetaStepperTypeString == "Trapezoid")
    thetaStepperType_ = Trapezoid;
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
		       "Value of " << ThetaStepperType_name << " = " << thetaStepperTypeString 
		       << " is invalid for Rythmos::ThetaStepper");

  default_predictor_order_ = 
    Teuchos::getParameter<int>(*pl_theta, PredictorOrder_name);
  
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  if ( Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_HIGH) ) {
    *out << ThetaStepperType_name << " = " << thetaStepperTypeString << std::endl;
  }
}


template <class Scalar>
RCP<Teuchos::ParameterList>
ThetaStepper<Scalar>::getNonconstParameterList()
{
  return(parameterList_);
}


template <class Scalar>
RCP<Teuchos::ParameterList>
ThetaStepper<Scalar>::unsetParameterList()
{
  RCP<Teuchos::ParameterList>
    temp_param_list = parameterList_;
  parameterList_ = Teuchos::null;
  return(temp_param_list);
}


template<class Scalar>
RCP<const Teuchos::ParameterList>
ThetaStepper<Scalar>::getValidParameters() const
{
  using Teuchos::ParameterList;

  static RCP<const ParameterList> validPL;

  if (is_null(validPL)) {

    RCP<ParameterList> pl_top_level = Teuchos::parameterList();

    RCP<ParameterList> pl = Teuchos::sublist(pl_top_level, RythmosStepControlSettings_name);

    pl->set<std::string> ( ThetaStepperType_name, ThetaStepperType_default,
        "Name of Stepper Type in Theta Stepper"
        );

    pl->set<int> ( PredictorOrder_name, PredictorOrder_default,
        "Order of Predictor in Theta Stepper, can be 0, 1, 2"
        );

    Teuchos::setupVerboseObjectSublist(&*pl_top_level);
    validPL = pl_top_level;
  }
  return validPL;
}


// Overridden from Teuchos::Describable


template<class Scalar>
void ThetaStepper<Scalar>::describe(
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
    out << "x_dot_base_ = " << std::endl;
    x_dot_base_->describe(out,verbLevel);
  }
}


// private


template <class Scalar>
void ThetaStepper<Scalar>::initialize_()
{

  typedef Teuchos::ScalarTraits<Scalar> ST;

  if (isInitialized_)
    return;

  TEUCHOS_TEST_FOR_EXCEPT(is_null(model_));
  TEUCHOS_TEST_FOR_EXCEPT(is_null(solver_));
  TEUCHOS_TEST_FOR_EXCEPT(!haveInitialCondition_);

#ifdef HAVE_RYTHMOS_DEBUG
  THYRA_ASSERT_VEC_SPACES(
    "Rythmos::ThetaStepper::setInitialCondition(...)",
    *x_->space(), *model_->get_x_space() );
#endif // HAVE_RYTHMOS_DEBUG

  if ( is_null(interpolator_) ) {
    // If an interpolator has not been explicitly set, then just create
    // a default linear interpolator.
    interpolator_ = Teuchos::rcp(new LinearInterpolator<Scalar> );
    // 2007/05/18: rabartl: ToDo: Replace this with a Hermete interplator
    // when it is implementated!
  }

  if (thetaStepperType_ == ImplicitEuler)
  {
    theta_ = 1.0;
    predictor_corrector_begin_after_step_ = 2;
  }
  else
  {
    theta_ = 0.5;
    predictor_corrector_begin_after_step_ = 3;
  }

  isInitialized_ = true;

}

template<class Scalar>
void ThetaStepper<Scalar>::obtainPredictor_()
{

  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  if (!isInitialized_) {
    return;
  }

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "Obtaining predictor..." << std::endl;
  }

  const int preferred_predictor_order = std::min(default_predictor_order_, thetaStepperType_ + 1);
  const int max_predictor_order_at_this_timestep = std::max(0, numSteps_);

  const int predictor_order = std::min(preferred_predictor_order, max_predictor_order_at_this_timestep);

  switch (predictor_order) 
  {
    case 0:
      V_StV(x_pre_.ptr(), Scalar(ST::one()), *x_old_);
      break;
    case 1:
    {
      V_StV(x_pre_.ptr(), Scalar(ST::one()), *x_old_);

      TEUCHOS_TEST_FOR_EXCEPT (dt_ <= 0.0);

      Vp_StV(x_pre_.ptr(), dt_, *x_dot_old_);
    }
    break;
    case 2:
    {
      V_StV(x_pre_.ptr(), Scalar(ST::one()), *x_old_);

      TEUCHOS_TEST_FOR_EXCEPT (dt_ <= 0.0);
      TEUCHOS_TEST_FOR_EXCEPT (dt_old_ <= 0.0);

      const Scalar coeff_x_dot_old = (0.5 * dt_) * (2.0 + dt_/dt_old_);
      const Scalar coeff_x_dot_really_old = - (0.5 * dt_) * (dt_/dt_old_);
      
      if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
	*out << "x_dot_old_ = " << *x_dot_old_ << std::endl;
	*out << "x_dot_really_old_ = " << *x_dot_really_old_ << std::endl;
      }
      
      Vp_StV( x_pre_.ptr(), coeff_x_dot_old, *x_dot_old_);
      Vp_StV( x_pre_.ptr(), coeff_x_dot_really_old, *x_dot_really_old_);
    }
    break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
			 "Invalid predictor order " << predictor_order << ". Valid values are 0, 1, and 2.");
  }
  
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "x_pre_ = " << *x_pre_ << std::endl;
  }

  // copy to current solution
  V_StV(x_.ptr(), Scalar(ST::one()), *x_pre_);
}

// 
// Explicit Instantiation macro
//
// Must be expanded from within the Rythmos namespace!
//

#define RYTHMOS_THETA_STEPPER_INSTANT(SCALAR) \
  \
  template class ThetaStepper< SCALAR >; \
  \
  template RCP< ThetaStepper< SCALAR > > \
  thetaStepper( \
    const RCP<Thyra::ModelEvaluator< SCALAR > >& model, \
    const RCP<Thyra::NonlinearSolverBase< SCALAR > >& solver, \
    RCP<Teuchos::ParameterList>& parameterList \
      );
   

} // namespace Rythmos

#endif // HAVE_RYTHMOS_EXPERIMENTAL

#endif //Rythmos_THETA_STEPPER_DEF_H
