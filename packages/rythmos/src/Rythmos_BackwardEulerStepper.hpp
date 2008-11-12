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

#ifndef Rythmos_BACKWARD_EULER_STEPPER_H
#define Rythmos_BACKWARD_EULER_STEPPER_H

#include "Rythmos_StepperBase.hpp"
#include "Rythmos_DataStore.hpp"
#include "Rythmos_LinearInterpolator.hpp"
#include "Rythmos_InterpolatorBaseHelpers.hpp"
#include "Rythmos_SingleResidualModelEvaluator.hpp"
#include "Rythmos_SolverAcceptingStepperBase.hpp"
#include "Rythmos_StepperHelpers.hpp"

#include "Thyra_VectorBase.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_AssertOp.hpp"
#include "Thyra_NonlinearSolverBase.hpp"
#include "Thyra_TestingTools.hpp"

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_as.hpp"


namespace Rythmos {


/** \brief Simple concrete stepper subclass implementing an implicit backward
 * Euler method.
 *
 * This class exists primarily as a simple example of an implicit time stepper
 * and as a vehicle for experimentation.  The <tt>ImplicitBDFStepper</tt> also
 * implements backward Euler and is a more powerful stepper class.  This class
 * does not implement a local truncation error test and therefore also does
 * not handle the automatic step size selection.  Therefore, if you need these
 * features, you should really use the <tt>ImplicitBDFStepper</tt> class.
 */
template<class Scalar>
class BackwardEulerStepper : virtual public SolverAcceptingStepperBase<Scalar>
{
public:
  
  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
  
  /** \name Constructors, intializers, Misc. */
  //@{

  /** \brief . */
  BackwardEulerStepper();
  
  /** \brief . */
  BackwardEulerStepper(
    const RCP<const Thyra::ModelEvaluator<Scalar> > &model,
    const RCP<Thyra::NonlinearSolverBase<Scalar> > &solver
    );
  
  /** \brief . */
  void setInterpolator(RCP<InterpolatorBase<Scalar> > interpolator);
  
  /** \brief . */
  RCP<InterpolatorBase<Scalar> > unsetInterpolator();

  //@}

  /** \name Overridden from SolverAcceptingStepperBase */
  //@{

  /** \brief . */
  void setSolver(
    const RCP<Thyra::NonlinearSolverBase<Scalar> > &solver
    );

  /** \brief . */
  RCP<Thyra::NonlinearSolverBase<Scalar> >
  getNonconstSolver();

  /** \brief . */
  RCP<const Thyra::NonlinearSolverBase<Scalar> >
  getSolver() const;

  //@}

  /** \name Overridden from StepperBase */
  //@{
 
  /** \brief Returns true. */
  bool supportsCloning() const;

  /** \brief Creates copies of all internal data (including the parameter
   * list) except the model which is assumed to stateless.
   *
   * If a shallow copy of the model is not appropirate for some reasone, then
   * the client can simply reset the model using
   * <tt>returnVal->setModel()</tt>.
   */
  RCP<StepperBase<Scalar> > cloneStepperAlgorithm() const;

  /** \brief . */
  void setModel(const RCP<const Thyra::ModelEvaluator<Scalar> > &model);
  
  /** \brief . */
  RCP<const Thyra::ModelEvaluator<Scalar> >
  getModel() const;

  /** \brief . */
  void setInitialCondition(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
    );

  /** \brief . */
  Scalar takeStep(Scalar dt, StepSizeType flag);
  
  /** \brief . */
  const StepStatus<Scalar> getStepStatus() const;
  
  //@}

  /** \name Overridden from InterpolationBufferBase */
  //@{

  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> >
  get_x_space() const;

  /** \brief . */
  void addPoints(
    const Array<Scalar>& time_vec,
    const Array<RCP<const Thyra::VectorBase<Scalar> > >& x_vec,
    const Array<RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec
    );
  
  /** \brief . */
  TimeRange<Scalar> getTimeRange() const;
  
  /** \brief . */
  void getPoints(
    const Array<Scalar>& time_vec,
    Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
    Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
    Array<ScalarMag>* accuracy_vec
    ) const;
  
  /** \brief . */
  void getNodes(Array<Scalar>* time_vec) const;
  
  /** \brief . */
  void removeNodes(Array<Scalar>& time_vec);

  /** \brief . */
  int getOrder() const;

  //@}
  
  /** \name Overridden from Teuchos::ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(RCP<Teuchos::ParameterList> const& paramList);
  
  /** \brief . */
  RCP<Teuchos::ParameterList> getNonconstParameterList();
  
  /** \brief . */
  RCP<Teuchos::ParameterList> unsetParameterList();
  
  /** \brief. */
  RCP<const Teuchos::ParameterList> getValidParameters() const;
 
  //@}

  /** \name Overridden from Teuchos::Describable */
  //@{
  
  /** \brief . */
  void describe(
    Teuchos::FancyOStream  &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;

  //@}

private:

  // ///////////////////////
  // Private date members

  bool isInitialized_;
  bool haveInitialCondition_;
  RCP<const Thyra::ModelEvaluator<Scalar> > model_;
  RCP<Thyra::NonlinearSolverBase<Scalar> > solver_;
  RCP<Thyra::VectorBase<Scalar> > scaled_x_old_;
  RCP<Thyra::VectorBase<Scalar> > x_dot_old_;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> basePoint_;
  RCP<Thyra::VectorBase<Scalar> > x_;
  RCP<Thyra::VectorBase<Scalar> > x_dot_;
  Scalar t_;
  Scalar t_old_;

  Scalar dt_;
  int numSteps_;

  RCP<Rythmos::SingleResidualModelEvaluator<Scalar> >  neModel_;

  RCP<Teuchos::ParameterList> parameterList_;

  RCP<InterpolatorBase<Scalar> > interpolator_;


  // //////////////////////////
  // Private member functions

  void initialize();

};


/** \brief Nonmember constructor.
 *
 * \relates BackwardEulerStepper
 */
template<class Scalar>
RCP<BackwardEulerStepper<Scalar> >
backwardEulerStepper(
  const RCP<const Thyra::ModelEvaluator<Scalar> > &model,
  const RCP<Thyra::NonlinearSolverBase<Scalar> > &solver
  )
{
  return Teuchos::rcp(new BackwardEulerStepper<Scalar>(model, solver));
}


// ////////////////////////////
// Defintions


// Constructors, intializers, Misc.


template<class Scalar>
BackwardEulerStepper<Scalar>::BackwardEulerStepper()
  :isInitialized_(false),
   haveInitialCondition_(false),
   t_(-1.0),
   t_old_(0.0),
   dt_(0.0),
   numSteps_(0)
{}


template<class Scalar>
BackwardEulerStepper<Scalar>::BackwardEulerStepper(
  const RCP<const Thyra::ModelEvaluator<Scalar> > &model,
  const RCP<Thyra::NonlinearSolverBase<Scalar> > &solver
  )
  :isInitialized_(false),
   haveInitialCondition_(false),
   t_(-1.0),
   t_old_(0.0),
   dt_(0.0),
   numSteps_(0)
{
  setModel(model);
  setSolver(solver);
}


template<class Scalar>
void BackwardEulerStepper<Scalar>::setInterpolator(
  RCP<InterpolatorBase<Scalar> > interpolator
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(is_null(interpolator));
#endif
  interpolator_ = interpolator;
  isInitialized_ = false;
}


template<class Scalar>
RCP<InterpolatorBase<Scalar> >
BackwardEulerStepper<Scalar>::unsetInterpolator()
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

  TEST_FOR_EXCEPT(solver == Teuchos::null)

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
void BackwardEulerStepper<Scalar>::setModel(
  const RCP<const Thyra::ModelEvaluator<Scalar> > &model
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
  x_ = Teuchos::null;
  scaled_x_old_ = Teuchos::null;
  x_dot_ = Teuchos::null;
  x_dot_old_ = Teuchos::null;

  isInitialized_ = false;
  haveInitialCondition_ = setDefaultInitialConditionFromNominalValues<Scalar>(
    *model_, Teuchos::ptr(this) );
  
}


template<class Scalar>
RCP<const Thyra::ModelEvaluator<Scalar> >
BackwardEulerStepper<Scalar>::getModel() const
{
  return model_;
}


template<class Scalar>
void BackwardEulerStepper<Scalar>::setInitialCondition(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
  )
{

  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Thyra::ModelEvaluatorBase MEB;

  TEST_FOR_EXCEPT( is_null(model_) );

  basePoint_ = initialCondition;

  // x

  RCP<const Thyra::VectorBase<Scalar> >
    x_init = initialCondition.get_x();

#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(
    is_null(x_init), std::logic_error,
    "Error, if the client passes in an intial condition to setInitialCondition(...),\n"
    "then x can not be null!" );
  THYRA_ASSERT_VEC_SPACES(
    "Rythmos::BackwardEulerStepper::setInitialCondition(...)",
    *x_init->space(), *model_->get_x_space() );
#endif

  x_ = x_init->clone_v();

  // x_dot

  x_dot_ = createMember(model_->get_x_space());

  RCP<const Thyra::VectorBase<Scalar> >
    x_dot_init = initialCondition.get_x_dot();

  if (!is_null(x_dot_init))
    assign(&*x_dot_,*x_dot_init);
  else
    assign(&*x_dot_,ST::zero());
  
  // t
  
  t_ = initialCondition.get_t();

  t_old_ = t_;

  haveInitialCondition_ = true;

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

  V_StV( &*scaled_x_old_, Scalar(-ST::one()/dt), *x_ );
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

  assign( &*x_dot_old_, *x_dot_ );

  // x_dot = (1/dt)*x - (1/dt)*x_old 
  V_StV( &*x_dot_, Scalar(ST::one()/dt), *x_ );
  Vp_StV( &*x_dot_, Scalar(ST::one()), *scaled_x_old_ );

  t_ += dt;

  numSteps_++;

  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "\nt_old_ = " << t_old_ << std::endl;
    *out << "\nt_ = " << t_ << std::endl;
  }

#ifdef TEUCHOS_DEBUG

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

#endif

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

#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(
    time_vec.size() == 0, std::logic_error,
    "Error, addPoints called with an empty time_vec array!\n");
#endif // TEUCHOS_DEBUG

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
    Thyra::V_V(&*x_,*x_vec[n]);
    Thyra::V_V(&*scaled_x_old_,*x_);
  }
  else {
    int n = time_vec.size()-1;
    int nm1 = time_vec.size()-2;
    t_ = time_vec[n];
    t_old_ = time_vec[nm1];
    Thyra::V_V(&*x_,*x_vec[n]);
    Scalar dt = t_ - t_old_;
    Thyra::V_StV(&*scaled_x_old_,Scalar(-ST::one()/dt),*x_vec[nm1]);
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

  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef Teuchos::ScalarTraits<ScalarMag> SMT;
  typename DataStore<Scalar>::DataStoreVector_t ds_nodes;
  typename DataStore<Scalar>::DataStoreVector_t ds_out;

#ifdef TEUCHOS_DEBUG
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
#ifdef TEUCHOS_DEBUG
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

}


template<class Scalar>
void BackwardEulerStepper<Scalar>::getNodes(Array<Scalar>* time_vec) const
{
  using Teuchos::as;


#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(!haveInitialCondition_);
#endif

  time_vec->clear();
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

  if ( is_null(interpolator_) ) {
    // If an interpolator has not been explicitly set, then just create
    // a default linear interpolator.
    interpolator_ = Teuchos::rcp(new LinearInterpolator<Scalar> );
    // 2007/05/18: rabartl: ToDo: Replace this with a Hermete interplator
    // when it is implementated!
  }

  if (is_null(scaled_x_old_))
    scaled_x_old_ = createMember(model_->get_x_space());

  if (is_null(x_dot_))
    x_dot_ = createMember(model_->get_x_space());

  if (is_null(x_dot_old_))
    x_dot_old_ = createMember(model_->get_x_space());

  // Note: above, we don't need to actually initialize x_dot or x_dot_old
  // since these will be initialized after each step

  isInitialized_ = true;

}


} // namespace Rythmos

#endif //Rythmos_BACKWARD_EULER_STEPPER_H
