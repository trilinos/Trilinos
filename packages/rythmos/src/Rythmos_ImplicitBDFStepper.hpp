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

#ifndef Rythmos_IMPLICITBDF_STEPPER_H
#define Rythmos_IMPLICITBDF_STEPPER_H

#include "Rythmos_StepperBase.hpp"
#include "Rythmos_SingleResidualModelEvaluator.hpp"
#include "Rythmos_SolverAcceptingStepperBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_NonlinearSolverBase.hpp"
#include "Thyra_SolveSupportTypes.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_as.hpp"


namespace Rythmos {


/** \brief . */
enum BDFactionFlag { ACTION_UNSET, ACTION_LOWER, ACTION_MAINTAIN, ACTION_RAISE };


/** \brief . */
enum BDFstatusFlag { PREDICT_AGAIN, CONTINUE_ANYWAY, REP_ERR_FAIL, REP_CONV_FAIL };


/** \brief . */
template<class Scalar>
class ImplicitBDFStepper : virtual public SolverAcceptingStepperBase<Scalar>
{
public:

  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
  
  /** \brief Constructors, intializers, Misc. */
  //@{

  /** \brief . */
  ImplicitBDFStepper();

  /** \brief . */
  ImplicitBDFStepper(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >  &model
    ,const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >  &solver
    );

  /** \brief . */
  ImplicitBDFStepper(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >  &model
    ,const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >  &solver
    ,Teuchos::RCP<Teuchos::ParameterList> &parameterList
    );

  /** \brief . */
  bool ErrWtVecSet(
    Thyra::VectorBase<Scalar> *w, 
    const Thyra::VectorBase<Scalar> &y
    );
  // 2007/06/11: rabartl : ToDo: Change above name to errWtVecSet

  /** \brief . */
  Scalar WRMSNorm(
    const Thyra::VectorBase<Scalar> &w
    , const Thyra::VectorBase<Scalar> &y
    ) const;
  // 2007/06/11: rabartl : ToDo: Change above name to wrmNorm

  /** \brief . */
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > get_solution() const;

  /** \brief . */
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > get_residual() const;

  //@}

  /** \name Overridden from SolverAcceptingStepperBase */
  //@{

  /** \brief . */
  void setSolver(
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > &solver
    );

  /** \brief . */
  Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >
  getNonconstSolver();

  /** \brief . */
  Teuchos::RCP<const Thyra::NonlinearSolverBase<Scalar> >
  getSolver() const;

  //@}

  /** \brief Overridden from StepperBase */
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
  Teuchos::RCP<StepperBase<Scalar> > cloneStepper() const;

  /** \brief . */
  void setModel(const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &model);

  /** \brief . */
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
  getModel() const;

  /** \brief . */
  void setInitialCondition(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
    );

  /** \brief . */
  Scalar takeStep(Scalar dt, StepSizeType flag);

  /** \brief . */
  const StepStatus<Scalar> getStepStatus() const;

  /** \brief . */
  void setStepControlData(const StepperBase<Scalar> & stepper);

  //@}

  /** \name Overridden from InterpolationBufferBase */
  //@{

  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
  get_x_space() const;

  /** \brief . */
  bool setPoints(
    const Array<Scalar>& time_vec
    ,const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x_vec
    ,const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec
    ,const Array<ScalarMag> & accuracy_vec 
    );

  /** \brief . */
  bool setRange(
    const TimeRange<Scalar>& range,
    const InterpolationBufferBase<Scalar> & IB
    );

  /** \brief . */
  TimeRange<Scalar> getTimeRange() const;
    
  /** \brief . */
  bool getPoints(
    const Array<Scalar>& time_vec
    ,Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >* x_vec
    ,Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec
    ,Array<ScalarMag>* accuracy_vec
    ) const;

  /** \brief . */
  bool getNodes(Array<Scalar>* time_vec) const;

  /** \brief . */
  bool removeNodes(Array<Scalar>& time_vec);

  /** \brief . */
  int getOrder() const;

  //@}
  
  /** \name Overridden from Teuchos::ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);

  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> getParameterList();

  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();

  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  //@}

  /** \name Overridden from Teuchos::Describable */
  //@{

  /** \brief . */
  void describe(
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;

  //@}

private:

  //
  // Private data members
  //

  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > model_;
  Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver_;
  Rythmos::SingleResidualModelEvaluator<Scalar>   neModel_;

  Teuchos::RCP<Thyra::VectorBase<Scalar> > xn0_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > xpn0_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > x_dot_base_;
  Array<Teuchos::RCP<Thyra::VectorBase<Scalar> > > xHistory_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > ee_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > delta_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > residual_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > errWtVec_;

  Scalar time_;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> basePoint_;

  //typedef typename Thyra::ModelEvaluatorBase::InArgs<Scalar>::ScalarMag ScalarMag;
  ScalarMag relErrTol_; // relative error tolerance
  ScalarMag absErrTol_; // absolute error tolerance
  Scalar hh_;        // Current step-size
  int currentOrder_; // Current order of integration
  int oldOrder_;     // previous order of integration
  int maxOrder_;     // maximum order = std::min(5,user option maxord) - see below.
  int usedOrder_;    // order used in current step (used after currentOrder is updated)
  Scalar alpha_s_;    // $\alpha_s$ fixed-leading coefficient of this BDF method
  std::vector<Scalar> alpha_;    // $\alpha_j(n)=h_n/\psi_j(n)$ coefficient used in local error test
  // note:   $h_n$ = current step size, n = current time step
  Scalar alpha_0_;     // $-\sum_{j=1}^k \alpha_j(n)$ coefficient used in local error test
  Scalar cj_ ;        // $-\alpha_s/h_n$ coefficient used in local error test
  Scalar ck_ ;        // local error coefficient gamma[0] = 0; 
  Scalar ck_enorm_;   // ck * enorm
  EStepLETStatus stepLETStatus_; // Local Error Test Status
  std::vector<Scalar> gamma_;    // $\gamma_j(n)=\sum_{l=1}^{j-1}1/\psi_l(n)$ coefficient used to
  // calculate time derivative of history array for predictor 
  std::vector<Scalar> beta_;     // coefficients used to evaluate predictor from history array
  std::vector<Scalar> psi_;      // $\psi_j(n) = t_n-t_{n-j}$ intermediary variable used to 
  // compute $\beta_j(n)$
  std::vector<Scalar> sigma_;    // $\sigma_j(n) = \frac{h_n^j(j-1)!}{\psi_1(n)*\cdots *\psi_j(n)}$
  int numberOfSteps_;// number of total time integration steps taken
  int  nef_;
  Scalar usedStep_;
  int nscsco_;
  Scalar Ek_;
  Scalar Ekm1_;
  Scalar Ekm2_;
  Scalar Ekp1_;
  Scalar Est_;
  Scalar Tk_;
  Scalar Tkm1_;
  Scalar Tkm2_;
  Scalar Tkp1_;
  int newOrder_;
  bool initialPhase_;
  Scalar stopTime_;
  bool constantStepSize_;
  bool haveInitialCondition_;
  bool isInitialized_;

  // Magic Numbers:
  Scalar h0_safety_;
  Scalar h0_max_factor_;
  Scalar h_phase0_incr_;
  Scalar h_max_inv_;
  Scalar Tkm1_Tk_safety_;
  Scalar Tkp1_Tk_safety_;
  Scalar r_factor_;
  Scalar r_safety_;
  Scalar r_fudge_;
  Scalar r_min_;
  Scalar r_max_;
  Scalar r_hincr_test_;
  Scalar r_hincr_;
  int    max_LET_fail_;
  Scalar minTimeStep_;
  Scalar maxTimeStep_;

  int newtonConvergenceStatus_;

  Teuchos::RCP<Teuchos::ParameterList> parameterList_;

  //
  // Private member functions
  //

  void setDefaultMagicNumbers_(Teuchos::ParameterList &magicNumberList);
  void getInitialCondition_();
  void obtainPredictor_();
  bool interpolateSolution_(
    const Scalar& timepoint,
    Thyra::VectorBase<Scalar>* x_ptr_,
    Thyra::VectorBase<Scalar>* xdot_ptr_,
    ScalarMag* accuracy_ptr_
    ) const;
  void updateHistory_();
  void restoreHistory_();
  void updateCoeffs_();
  void initialize_();
  Scalar checkReduceOrder_();
  BDFstatusFlag rejectStep_();
  void completeStep_();

};


// ////////////////////////////
// Defintions


// Constructors, intializers, Misc.


template<class Scalar>
ImplicitBDFStepper<Scalar>::ImplicitBDFStepper()
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  out->setMaxLenLinePrefix(30);
  haveInitialCondition_ = false;
  isInitialized_=false;
}


template<class Scalar>
ImplicitBDFStepper<Scalar>::ImplicitBDFStepper(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &model
  ,const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > &solver
  ,Teuchos::RCP<Teuchos::ParameterList> &parameterList
  )
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  out->setMaxLenLinePrefix(30);
  this->setParameterList(parameterList);
  // Now we instantiate the model and the solver
  setModel(model);
  setSolver(solver);
  haveInitialCondition_ = false;
  isInitialized_=false;
}


template<class Scalar>
ImplicitBDFStepper<Scalar>::ImplicitBDFStepper(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &model
  ,const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > &solver
  )
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  out->setMaxLenLinePrefix(30);
  // Now we instantiate the model and the solver
  setModel(model);
  setSolver(solver);
  haveInitialCondition_ = false;
  isInitialized_=false;
}


template<class Scalar>
bool ImplicitBDFStepper<Scalar>::ErrWtVecSet(Thyra::VectorBase<Scalar> *w_in, const Thyra::VectorBase<Scalar> &y)
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  Thyra::VectorBase<Scalar> &w = *w_in;
  Thyra::abs(&w,y);
  Vt_S(&w,relErrTol_);
  Vp_S(&w,absErrTol_);
  reciprocal(&w,w);
  Vt_StV(&w,ST::one(),w); // We square w because of how weighted norm_2 is computed.
  // divide by N to get RMS norm
  int N = y.space()->dim();
  Vt_S(&w,Scalar(1.0/N));
  // Now you can compute WRMS norm as:
  // Scalar WRMSnorm = norm_2(w,y); // WRMS norm of y with respect to weights w.
  return true; // This should be updated to reflect success of reciprocal
}


template<class Scalar>
Scalar ImplicitBDFStepper<Scalar>::WRMSNorm(
  const Thyra::VectorBase<Scalar> &w, const Thyra::VectorBase<Scalar> &y
  ) const
{
  return(norm_2(w,y));
}


// Overridden from SolverAcceptingStepperBase


template<class Scalar>
void ImplicitBDFStepper<Scalar>::setSolver(const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > &solver)
{
  TEST_FOR_EXCEPT(solver == Teuchos::null)
    solver_ = solver;
}


template<class Scalar>
Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >
ImplicitBDFStepper<Scalar>::getNonconstSolver()
{
  return (solver_);
}


template<class Scalar>
Teuchos::RCP<const Thyra::NonlinearSolverBase<Scalar> >
ImplicitBDFStepper<Scalar>::getSolver() const
{
  return (solver_);
}


// Overridden from StepperBase


template<class Scalar>
bool ImplicitBDFStepper<Scalar>::supportsCloning() const
{
  return true;
}


template<class Scalar>
Teuchos::RCP<StepperBase<Scalar> >
ImplicitBDFStepper<Scalar>::cloneStepper() const
{

  // Just use the interface to clone the algorithm in an basically
  // uninitialized state

  Teuchos::RCP<ImplicitBDFStepper<Scalar> >
    stepper = Teuchos::rcp(new ImplicitBDFStepper<Scalar>());

  if (!is_null(model_))
    stepper->setModel(model_); // Shallow copy is okay!

  if (!is_null(solver_))
    stepper->setSolver(solver_->cloneNonlinearSolver().assert_not_null());
  
  if (!is_null(parameterList_))
    stepper->setParameterList(Teuchos::parameterList(*parameterList_));

  // At this point, we should have a valid algorithm state.  What might be
  // missing is the initial condition if it was not given in *model_ but was
  // set explicitly.  However, the specification for this function does not
  // guarantee that the full state will be copied in any case!

  return stepper;

}


template<class Scalar>
void ImplicitBDFStepper<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &model
  )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  TEST_FOR_EXCEPT(model == Teuchos::null)
    model_= model;
}


template<class Scalar>
Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
ImplicitBDFStepper<Scalar>::getModel() const
{
  return model_;
}


template<class Scalar>
void ImplicitBDFStepper<Scalar>::setInitialCondition(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
  )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  TEST_FOR_EXCEPT(initialCondition.get_x()==Teuchos::null);
  TEST_FOR_EXCEPT(initialCondition.get_x_dot()==Teuchos::null);
  basePoint_ = initialCondition;
  xn0_ = initialCondition.get_x()->clone_v();
  xpn0_ = initialCondition.get_x_dot()->clone_v(); 
  typedef Thyra::ModelEvaluatorBase MEB;
  time_ = initialCondition.get_t();
  haveInitialCondition_ = true;
  if (isInitialized_) {
    initialize_();
  }
}


template<class Scalar>
Scalar ImplicitBDFStepper<Scalar>::takeStep(Scalar dt, StepSizeType stepType)
{

  using Teuchos::as;
  using Teuchos::incrVerbLevel;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename Thyra::ModelEvaluatorBase::InArgs<Scalar>::ScalarMag ScalarMag;
  typedef Thyra::NonlinearSolverBase<Scalar> NSB;
  typedef Teuchos::VerboseObjectTempState<NSB> VOTSNSB;

  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
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

  if (stepType == FIXED_STEP) {
    constantStepSize_ = true;
    if (dt != ST::zero()) {
      hh_ = dt;
      if (numberOfSteps_ == 0) {
        psi_[0] = hh_;
        cj_ = 1/psi_[0];
      }
    }
    if (hh_ == ST::zero()) {
      return(Scalar(-ST::one()));
    }
  }
  else {
    constantStepSize_ = false;
  }

  if ((stepType == VARIABLE_STEP) && (dt != ST::zero())) {
    maxTimeStep_ = dt;
    h_max_inv_ = Scalar(ST::one()/maxTimeStep_);
  }
  if (hh_ > maxTimeStep_) {
    hh_ = maxTimeStep_;
  }
  
  // Set Error weight std::vector once for each step.
  ErrWtVecSet(&*errWtVec_,*xHistory_[0]);

  BDFstatusFlag status;
  while (1) {
    // Set up problem coefficients (and handle first step)
    updateCoeffs_();
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
    Thyra::SolveStatus<Scalar> nonlinearSolveStatus = solver_->solve( &*xn0_, NULL, &*ee_ ); 
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
    Scalar enorm = checkReduceOrder_();
    ck_enorm_ = ck_*enorm;
    
    if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
      *out << "xn0_ = " << std::endl;
      xn0_->describe(*out,verbLevel);
      *out << "ee_ = " << std::endl;
      ee_->describe(*out,verbLevel);
      *out << "delta_ = " << std::endl;
      delta_->describe(*out,verbLevel);
      for (int i=0; i<std::max(2,maxOrder_); ++i) {
        *out << "xHistory_[" << i << "] = " << std::endl;
        xHistory_[i]->describe(*out,verbLevel);
      }
      *out << "ck_ = " << ck_ << std::endl;
      *out << "enorm = " << enorm << std::endl;
      *out << "Local Truncation Error Check: (ck*enorm) < 1:  (" << ck_enorm_ << ") <?= 1" << std::endl;
    }

    // Check LTE here:
    if ((ck_enorm_) > ST::one()) {
      stepLETStatus_ = STEP_LET_STATUS_FAILED;
      status = rejectStep_(); 
    }
    else {
      stepLETStatus_ = STEP_LET_STATUS_PASSED;
      break;
    }
    if (status == CONTINUE_ANYWAY) {
      break;
    }
    if (status == REP_ERR_FAIL) {
      return(Scalar(-ST::one()));
    }
  }

  completeStep_();  

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
  typedef Teuchos::ScalarTraits<Scalar> ST;
  StepStatus<Scalar> stepStatus;
  if (!isInitialized_) {
    stepStatus.message = "This stepper is uninitialized.";
    stepStatus.stepStatus = STEP_STATUS_UNINITIALIZED;
    stepStatus.stepSize = Scalar(-ST::one());
    stepStatus.order = -1;
    stepStatus.time = Scalar(-ST::one());
    if (model_ != Teuchos::null) {
      stepStatus.solution = Thyra::createMember(model_->get_x_space());
    }
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
  stepStatus.stepLETValue = ck_enorm_; 
  stepStatus.solution = xHistory_[0];
  stepStatus.solutionDot = xHistory_[1];
  stepStatus.residual = residual_;

  return(stepStatus);

}


// Overridden from InterpolationBufferBase


template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
ImplicitBDFStepper<Scalar>::get_x_space() const
{
  return ( !is_null(model_) ? model_->get_x_space() : Teuchos::null );
}


template<class Scalar>
bool ImplicitBDFStepper<Scalar>::setPoints(
  const Array<Scalar>& time_vec,
  const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x_vec,
  const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec,
  const Array<ScalarMag> & accuracy_vec 
  )
{
  return(false);
}


template<class Scalar>
bool ImplicitBDFStepper<Scalar>::setRange(
  const TimeRange<Scalar>& range,
  const InterpolationBufferBase<Scalar>& IB)
{
  return(false);
}


template<class Scalar>
TimeRange<Scalar> ImplicitBDFStepper<Scalar>::getTimeRange() const
{
  return invalidTimeRange<Scalar>();
}


template<class Scalar>
bool ImplicitBDFStepper<Scalar>::getPoints(
  const Array<Scalar>& time_vec
  ,Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >* x_vec
  ,Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec
  ,Array<ScalarMag>* accuracy_vec) const
{
  using Teuchos::as;
  bool status;
  for (unsigned int i=0 ; i<time_vec.size() ; ++i) {
    Teuchos::RCP<Thyra::VectorBase<Scalar> > x_temp = xn0_->clone_v();
    Teuchos::RCP<Thyra::VectorBase<Scalar> > xdot_temp = xn0_->clone_v();
    ScalarMag accuracy;
    status = interpolateSolution_(time_vec[i],&*x_temp,&*xdot_temp,&accuracy);
    if (!status) {
      return(status);
    }
    x_vec->push_back(x_temp);
    xdot_vec->push_back(xdot_temp);
    accuracy_vec->push_back(accuracy);
  }
  if ( as<int>(this->getVerbLevel()) >= as<int>(Teuchos::VERB_HIGH) ) {
    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"getPoints");
    *out << "Passing out the interpolated values:" << std::endl;
    for (unsigned int i=0; i<time_vec.size() ; ++i) {
      *out << "time_[" << i << "] = " << time_vec[i] << std::endl;
      *out << "x_vec[" << i << "] = " << std::endl;
      (*x_vec)[i]->describe(*out,this->getVerbLevel());
      if ( (*xdot_vec)[i] == Teuchos::null) {
        *out << "xdot_vec[" << i << "] = Teuchos::null" << std::endl;
      } else {
        *out << "xdot_vec[" << i << "] = " << std::endl;
        (*xdot_vec)[i]->describe(*out,this->getVerbLevel());
      }
      *out << "accuracy[" << i << "] = " << (*accuracy_vec)[i] << std::endl;
    }
  }
  return(status);
}


template<class Scalar>
bool ImplicitBDFStepper<Scalar>::getNodes(Array<Scalar>* time_vec) const
{
  if (!isInitialized_) {
    return(false);
  }
  if (numberOfSteps_ > 0) {
    time_vec->push_back(time_-usedStep_);
  }
  time_vec->push_back(time_);
  return(true);
}


template<class Scalar>
bool ImplicitBDFStepper<Scalar>::removeNodes(Array<Scalar>& time_vec) 
{
  return(false);
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
  Teuchos::RCP<Teuchos::ParameterList> const& paramList
  )
{

  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  TEST_FOR_EXCEPT(paramList == Teuchos::null);
  parameterList_ = paramList;

  Teuchos::readVerboseObjectSublist(&*parameterList_,this);

  maxOrder_ = parameterList_->get("maxOrder",int(5)); // maximum order
  maxOrder_ = std::max(1,std::min(maxOrder_,5)); // 1 <= maxOrder <= 5

  relErrTol_ = parameterList_->get( "relErrTol", Scalar(1.0e-4) );
  absErrTol_ = parameterList_->get( "absErrTol", Scalar(1.0e-6) );
  constantStepSize_ = parameterList_->get( "constantStepSize", false );
  stopTime_ = parameterList_->get( "stopTime", Scalar(1.0) );

  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"setParameterList");
  out->precision(15);

  setDefaultMagicNumbers_(parameterList_->sublist("magicNumbers"));

  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "maxOrder_ = " << maxOrder_ << std::endl;
    *out << "relErrTol  = " << relErrTol_  << std::endl;
    *out << "absErrTol  = " << absErrTol_  << std::endl;
    *out << "constantStepSize_  = " << constantStepSize_  << std::endl;
    *out << "stopTime_  = " << stopTime_  << std::endl;
  }

}


template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList> ImplicitBDFStepper<Scalar>::getParameterList()
{
  return(parameterList_);
}


template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
ImplicitBDFStepper<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_param_list = parameterList_;
  parameterList_ = Teuchos::null;
  return(temp_param_list);
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
ImplicitBDFStepper<Scalar>::getValidParameters() const
{

  static Teuchos::RCP<Teuchos::ParameterList> validPL;

  if (is_null(validPL)) {

    Teuchos::RCP<Teuchos::ParameterList>
      pl = Teuchos::parameterList();

    pl->set<int>   ( "maxOrder",         5              );
    pl->set<Scalar>( "relErrTol",        Scalar(1.0e-4) );
    pl->set<Scalar>( "absErrTol",        Scalar(1.0e-6) );
    pl->set<bool>  ( "constantStepSize", false          );
    pl->set<Scalar>( "stopTime",         Scalar(10.0)   );

    Teuchos::ParameterList
      &magicNumberList = pl->sublist("magicNumbers");
    magicNumberList.set<Scalar>( "h0_safety",      Scalar(2.0)     );
    magicNumberList.set<Scalar>( "h0_max_factor",  Scalar(0.001)   );
    magicNumberList.set<Scalar>( "h_phase0_incr",  Scalar(2.0)     );
    magicNumberList.set<Scalar>( "h_max_inv",      Scalar(0.0)     );
    magicNumberList.set<Scalar>( "Tkm1_Tk_safety", Scalar(2.0)     );
    magicNumberList.set<Scalar>( "Tkp1_Tk_safety", Scalar(0.5)     );
    magicNumberList.set<Scalar>( "r_factor",       Scalar(0.9)     );
    magicNumberList.set<Scalar>( "r_safety",       Scalar(2.0)     );
    magicNumberList.set<Scalar>( "r_fudge",        Scalar(0.0001)  );
    magicNumberList.set<Scalar>( "r_min",          Scalar(0.125)   );
    magicNumberList.set<Scalar>( "r_max",          Scalar(0.9)     );
    magicNumberList.set<Scalar>( "r_hincr_test",   Scalar(2.0)     );
    magicNumberList.set<Scalar>( "r_hincr",        Scalar(2.0)     );
    magicNumberList.set<int>   ( "max_LET_fail",   int(15)         );
    magicNumberList.set<Scalar>( "minTimeStep",    Scalar(0.0)     );
    magicNumberList.set<Scalar>( "maxTimeStep",    Scalar(10.0)    ); 

    Teuchos::setupVerboseObjectSublist(&*pl);

    validPL = pl;

  }

  return (validPL);
  
}


// Overridden from Teuchos::Describable


template<class Scalar>
void ImplicitBDFStepper<Scalar>::describe(
  Teuchos::FancyOStream &out,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{

  using Teuchos::as;

  if (!isInitialized_) {
    out << this->description() << "This stepper is not initialized yet" << std::endl;
    return;
  }

  if ( (as<int>(verbLevel) == as<int>(Teuchos::VERB_DEFAULT) ) ||
    (as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW)     )
    ) {
    out << this->description() << "::describe" << std::endl;
    out << "model_ = " << model_->description() << std::endl;
    out << "solver_ = " << solver_->description() << std::endl;
    out << "neModel_ = " << neModel_.description() << std::endl;
  }
  else if (as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW)) {
    out << "time_ = " << time_ << std::endl;
    out << "hh_ = " << hh_ << std::endl;
    out << "currentOrder_ = " << currentOrder_ << std::endl;
  }
  else if (as<int>(verbLevel) >= as<int>(Teuchos::VERB_MEDIUM)) {
  }
  else if (as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH)) {
    out << "model_ = " << std::endl;
    model_->describe(out,verbLevel);
    out << "solver_ = " << std::endl;
    solver_->describe(out,verbLevel);
    out << "neModel_ = " << std::endl;
    neModel_.describe(out,verbLevel);
    out << "xn0_ = " << std::endl;
    xn0_->describe(out,verbLevel);
    out << "xpn0_ = " << std::endl;
    xpn0_->describe(out,verbLevel);
    out << "x_dot_base_ = " << std::endl;
    x_dot_base_->describe(out,verbLevel);
    out << "xHistory_ = " << std::endl;
    for (int i=0 ; i < std::max(2,maxOrder_) ; ++i) {
      out << "xHistory_[" << i << "] = " << std::endl;
      xHistory_[i]->describe(out,verbLevel);
    }
    out << "ee_ = " << std::endl;
    ee_->describe(out,verbLevel);
    out << "delta_ = " << std::endl;
    delta_->describe(out,verbLevel);
    out << "residual_ = " << std::endl;
    residual_->describe(out,verbLevel);
    out << "errWtVec_ = " << std::endl;
    errWtVec_->describe(out,verbLevel);
  }
}


// private


template<class Scalar>
void ImplicitBDFStepper<Scalar>::setDefaultMagicNumbers_(
  Teuchos::ParameterList &magicNumberList)
{

  using Teuchos::as;

  // Magic Number Defaults:
  h0_safety_      = magicNumberList.get( "h0_safety",      Scalar(2.0)     );
  h0_max_factor_  = magicNumberList.get( "h0_max_factor",  Scalar(0.001)   );
  h_phase0_incr_  = magicNumberList.get( "h_phase0_incr",  Scalar(2.0)     );
  h_max_inv_      = magicNumberList.get( "h_max_inv",      Scalar(0.0)     );
  Tkm1_Tk_safety_ = magicNumberList.get( "Tkm1_Tk_safety", Scalar(2.0)     );
  Tkp1_Tk_safety_ = magicNumberList.get( "Tkp1_Tk_safety", Scalar(0.5)     );
  r_factor_       = magicNumberList.get( "r_factor",       Scalar(0.9)     );
  r_safety_       = magicNumberList.get( "r_safety",       Scalar(2.0)     );
  r_fudge_        = magicNumberList.get( "r_fudge",        Scalar(0.0001)  );
  r_min_          = magicNumberList.get( "r_min",          Scalar(0.125)   );
  r_max_          = magicNumberList.get( "r_max",          Scalar(0.9)     );
  r_hincr_test_   = magicNumberList.get( "r_hincr_test",   Scalar(2.0)     );
  r_hincr_        = magicNumberList.get( "r_hincr",        Scalar(2.0)     );
  max_LET_fail_   = magicNumberList.get( "max_LET_fail",   int(15)         );
  minTimeStep_    = magicNumberList.get( "minTimeStep",    Scalar(0.0)     );
  maxTimeStep_    = magicNumberList.get( "maxTimeStep",    Scalar(10.0)    ); 

  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"setDefaultMagicNumbers_");
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "h0_safety_ = " << h0_safety_ << std::endl;
    *out << "h0_max_factor_ = " << h0_max_factor_ << std::endl;
    *out << "h_phase0_incr_ = " << h_phase0_incr_ << std::endl;
    *out << "h_max_inv_ = " << h_max_inv_ << std::endl;
    *out << "Tkm1_Tk_safety_ = " << Tkm1_Tk_safety_  << std::endl;
    *out << "Tkp1_Tk_safety_ = " << Tkp1_Tk_safety_ << std::endl;
    *out << "r_factor_ = " << r_factor_ << std::endl;
    *out << "r_safety_ = " << r_safety_ << std::endl;
    *out << "r_fudge_ = " << r_fudge_ << std::endl;
    *out << "r_min_ = " << r_min_ << std::endl;
    *out << "r_max_ = " << r_max_ << std::endl;
    *out << "r_hincr_test_ = " << r_hincr_test_ << std::endl;
    *out << "r_hincr_ = " << r_hincr_ << std::endl;
    *out << "max_LET_fail_ = " << max_LET_fail_ << std::endl;
    *out << "minTimeStep_ = " << minTimeStep_ << std::endl;
    *out << "maxTimeStep_ = " << maxTimeStep_ << std::endl;
  }

}


template<class Scalar>
void ImplicitBDFStepper<Scalar>::getInitialCondition_()
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  if (!haveInitialCondition_) {
    TEST_FOR_EXCEPT(model_->getNominalValues().get_x()==Teuchos::null);
    TEST_FOR_EXCEPT(model_->getNominalValues().get_x_dot()==Teuchos::null);
    basePoint_ = model_->getNominalValues();
    xn0_ = basePoint_.get_x()->clone_v();
    xpn0_ = basePoint_.get_x_dot()->clone_v(); 
    time_ = basePoint_.get_t();
    haveInitialCondition_ = true;
  }
}


template<class Scalar>
void ImplicitBDFStepper<Scalar>::obtainPredictor_()
{

  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  if (!isInitialized_) {
    return;
  }

  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
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
bool ImplicitBDFStepper<Scalar>::interpolateSolution_(
  const Scalar& timepoint
  ,Thyra::VectorBase<Scalar>* x_ptr_
  ,Thyra::VectorBase<Scalar>* xdot_ptr_
  ,ScalarMag* accuracy_ptr_
  ) const
{
  if (!isInitialized_) {
    return(false);
  }
  typedef Teuchos::ScalarTraits<Scalar> ST;
  Scalar tfuzz;   // fuzz factor to check for valid output time
  Scalar tp;      // approximately t{n-1}
  Scalar delt;    // distance between timepoint and time
  Scalar c = ST::one(); // coefficient for interpolation of x
  Scalar d = ST::zero(); // coefficient for interpolation of xdot
  Scalar gam;     // coefficient for interpolation
  int kord;       // order of interpolation
  Scalar tn = time_;
  Scalar hused = usedStep_;
  int kused = usedOrder_;
  Scalar uround = ST::zero();  // unit round-off (set to zero for now)

  tfuzz = 100 * uround * (tn + hh_);
  tp = tn - hused - tfuzz;
  if ( (timepoint - tp)*hh_ < ST::zero() )  {
    return(false);
  }
  if ( timepoint - (time_-tfuzz) > ST::zero() ) {
    return(false);
  }

  Thyra::V_V(x_ptr_,*xHistory_[0]);
  Thyra::V_S(xdot_ptr_,ST::zero());
  kord = kused;
  if ( (kused == 0) || (timepoint == tn) )  {
    kord = 1;
  }

  delt = timepoint - tn;
  gam = delt/psi_[0];
  for (int j=1 ; j <= kord ; ++j) {
    d = d*gam + c/psi_[j-1];
    c = c*gam;
    gam = (delt + psi_[j-1])/psi_[j];
    Thyra::Vp_StV(x_ptr_,c,*xHistory_[j]);
    Thyra::Vp_StV(xdot_ptr_,d,*xHistory_[j]);
  }
  *accuracy_ptr_ = pow(usedStep_,kord);
  return(true);
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
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
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
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
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
    sigma_[0] = ST::one();
    gamma_[0] = ST::zero();
    for (int i=1;i<=currentOrder_;++i) {
      Scalar temp2 = psi_[i-1];
      psi_[i-1] = temp1;
      beta_[i] = beta_[i-1]*psi_[i-1]/temp2;
      temp1 = temp2 + hh_;
      alpha_[i] = hh_/temp1;
      sigma_[i] = Scalar(i+1)*sigma_[i-1]*alpha_[i];
      gamma_[i] = gamma_[i-1]+alpha_[i-1]/hh_;
    }
    psi_[currentOrder_] = temp1;
    alpha_s_ = ST::zero();
    alpha_0_ = ST::zero();
    for (int i=0;i<currentOrder_;++i) {
      alpha_s_ = alpha_s_ - Scalar(ST::one()/(i+ST::one()));
      alpha_0_ = alpha_0_ - alpha_[i];
    }
    cj_ = -alpha_s_/hh_;
    ck_ = std::abs(alpha_[currentOrder_]+alpha_s_-alpha_0_);
    ck_ = std::max(ck_,alpha_[currentOrder_]);
  }
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"updateCoeffs_");
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    for (int i=0;i<=maxOrder_;++i) {
      *out << "alpha_[" << i << "] = " << alpha_[i] << std::endl;
      *out << "beta_[" << i << "] = " << beta_[i] << std::endl;
      *out << "sigma_[" << i << "] = " << sigma_[i] << std::endl;
      *out << "gamma_[" << i << "] = " << gamma_[i] << std::endl;
      *out << "psi_[" << i << "] = " << psi_[i] << std::endl;
      *out << "alpha_s_ = " << alpha_s_ << std::endl;
      *out << "alpha_0_ = " << alpha_0_ << std::endl;
      *out << "cj_ = " << cj_ << std::endl;
      *out << "ck_ = " << ck_ << std::endl;
    }
  }

}


template<class Scalar>
void ImplicitBDFStepper<Scalar>::initialize_()
{

  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  using Thyra::createMember;

  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
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

  if (parameterList_ == Teuchos::null) {
    Teuchos::RCP<Teuchos::ParameterList> emptyParameterList;
    this->setParameterList(emptyParameterList);
  }

  currentOrder_=1; // Current order of integration
  oldOrder_=1; // previous order of integration
  usedOrder_=1;  // order used in current step (used after currentOrder_ is updated)
  alpha_s_=Scalar(-ST::one());  // $\alpha_s$ fixed-leading coefficient of this BDF method
  alpha_.clear();  // $\alpha_j(n)=h_n/\psi_j(n)$ coefficient used in local error test
  // note:   $h_n$ = current step size, n = current time step
  gamma_.clear();  // calculate time derivative of history array for predictor 
  beta_.clear();   // coefficients used to evaluate predictor from history array
  psi_.clear();    // $\psi_j(n) = t_n-t_{n-j}$ intermediary variable used to 
  // compute $\beta_j(n;$
  sigma_.clear();  // $\sigma_j(n) = \frac{h_n^j(j-1)!}{\psi_1(n)*\cdots *\psi_j(n)}$
  Scalar zero = ST::zero();
  for (int i=0 ; i<=maxOrder_ ; ++i) {
    alpha_.push_back(zero);
    beta_.push_back(zero);
    gamma_.push_back(zero);
    psi_.push_back(zero);
    sigma_.push_back(zero);
  }
  alpha_0_=zero;   // $-\sum_{j=1}^k \alpha_j(n)$ coefficient used in local error test
  cj_=zero;      // $-\alpha_s/h_n$ coefficient used in local error test
  ck_=zero;      // local error coefficient gamma_[0] = 0; // $\gamma_j(n)=\sum_{l=1}^{j-1}1/\psi_l(n)$ coefficient used to
  hh_=zero;
  numberOfSteps_=0;   // number of total time integration steps taken
  nef_=0;
  usedStep_=zero;
  nscsco_=0;
  Ek_=zero;
  Ekm1_=zero;
  Ekm2_=zero;
  Ekp1_=zero;
  Est_=zero;
  Tk_=zero;
  Tkm1_=zero;
  Tkm2_=zero;
  Tkp1_=zero;
  newOrder_=currentOrder_;
  initialPhase_=true;

  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "currentOrder_ = " << currentOrder_ << std::endl;
    *out << "oldOrder_ = " << oldOrder_ << std::endl;
    *out << "usedOrder_ = " << usedOrder_ << std::endl;
    *out << "alpha_s_ = " << alpha_s_ << std::endl;
    for (int i=0 ; i<=maxOrder_ ; ++i) {
      *out << "alpha_[" << i << "] = " << alpha_[i] << std::endl;
      *out << "beta_[" << i << "] = " << beta_[i] << std::endl;
      *out << "gamma_[" << i << "] = " << gamma_[i] << std::endl;
      *out << "psi_[" << i << "] = " << psi_[i] << std::endl;
      *out << "sigma_[" << i << "] = " << sigma_[i] << std::endl;
    }
    *out << "alpha_0_ = " << alpha_0_ << std::endl;
    *out << "cj_ = " << cj_ << std::endl;
    *out << "ck_ = " << ck_ << std::endl;
    *out << "numberOfSteps_ = " << numberOfSteps_ << std::endl;
    *out << "nef_ = " << nef_ << std::endl;
    *out << "usedStep_ = " << usedStep_ << std::endl;
    *out << "nscsco_ = " << nscsco_ << std::endl;
    *out << "Ek_ = " << Ek_ << std::endl;
    *out << "Ekm1_ = " << Ekm1_ << std::endl;
    *out << "Ekm2_ = " << Ekm2_ << std::endl;
    *out << "Ekp1_ = " << Ekp1_ << std::endl;
    *out << "Est_ = " << Est_ << std::endl;
    *out << "Tk_ = " << Tk_ << std::endl;
    *out << "Tkm1_ = " << Tkm1_ << std::endl;
    *out << "Tkm2_ = " << Tkm2_ << std::endl;
    *out << "Tkp1_ = " << Tkp1_ << std::endl;
    *out << "newOrder_ = " << newOrder_ << std::endl;
    *out << "initialPhase_ = " << initialPhase_ << std::endl;
  }

  this->getInitialCondition_();
  // Generate vectors for use in calculations
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > 
    x_space = model_->get_x_space();
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > 
    f_space = model_->get_f_space();
  x_dot_base_ = createMember(x_space);
  ee_ = createMember(x_space);
  delta_ = createMember(x_space);
  V_S(&*delta_,zero); // tscoffe 06/14/07:  Initialize delta_ 
  residual_ = createMember(f_space);
  errWtVec_ = createMember(x_space); 
  ErrWtVecSet(&*errWtVec_,*xn0_);
  xHistory_.clear();
  xHistory_.push_back(xn0_->clone_v());
  xHistory_.push_back(xpn0_->clone_v());
  // Store maxOrder_+1 vectors
  for (int i=2 ; i<=maxOrder_ ; ++i) {
    xHistory_.push_back(createMember(x_space)); 
    V_S(&*xHistory_[i],zero);
  }

  // Choose initial step-size
  Scalar time_to_stop = stopTime_ - time_;
  Scalar currentTimeStep;
  if (constantStepSize_) {
    currentTimeStep = hh_;
    //currentTimeStep = 0.1 * time_to_stop;
    //currentTimeStep = std::min(hh_, currentTimeStep);
  } else {
    // compute an initial step-size based on rate of change in the solution initially
    Scalar ypnorm = WRMSNorm(*errWtVec_,*xHistory_[1]);
    if (ypnorm > zero) { // time-dependent DAE
      currentTimeStep = std::min(h0_max_factor_*std::abs(time_to_stop),std::sqrt(2.0)/(h0_safety_*ypnorm));
    } else { // non-time-dependent DAE
      currentTimeStep = h0_max_factor_*std::abs(time_to_stop);
    }
    // choose std::min of user specified value and our value:
    if (hh_ > zero) {
      currentTimeStep = std::min(hh_, currentTimeStep);
    }
    // check for maximum step-size:
    Scalar rh = std::abs(currentTimeStep)*h_max_inv_; 
    if (rh>1.0) currentTimeStep = currentTimeStep/rh;
  }
  hh_ = currentTimeStep;
  if (doTrace) {
    *out << "\nhh_ = " << hh_ << std::endl;
  }
  
  // x history
  assign(&*xHistory_[0],*xn0_);
  V_S(&*xHistory_[1],zero);

  // Coefficient initialization 
  psi_[0] = hh_;
  cj_ = 1/psi_[0];

  isInitialized_ = true;

  if (doTrace) {
    *out
      << "\nLeaving " << this->Teuchos::Describable::description()
      << "::initialize_()...\n";
  }

}


template<class Scalar>
Scalar ImplicitBDFStepper<Scalar>::checkReduceOrder_()
{

  using Teuchos::as;

  // This routine puts its output in newOrder_
  
  // This routine changes the following variables:
  //    Ek, Tk, Est, newOrder, dsDae.delta_x, dsDae.delta_q,
  //    Ekm1, Tkm1, Ekm2, Tkm2 

  // This routine reads but does not change the following variables:
  //    currentOrder_, sigma_, dsDae.newtonCorrectionPtr, dsDae.qNewtonCorrectionPtr,
  //    dsDae.errWtVecPtr, dsDae.qErrWtVecPtr, dsDae.xHistory, dsDae.qHistory

  // 2007/06/11: rabartl: This class has some much private data that it might
  // be good to partition it into two or more classes.  Statements like the
  // above will be hard to maintain and suggest a need for some partitioning.

  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"checkReduceOrder_");
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "sigma_[" << currentOrder_ << "] = " << sigma_[currentOrder_] << std::endl;
    *out << "ee_ = " << std::endl;
    ee_->describe(*out,verbLevel);
    *out << "errWtVec_ = " << std::endl;
    errWtVec_->describe(*out,verbLevel);
  }

  Scalar enorm = WRMSNorm(*errWtVec_,*ee_);
  Ek_ = sigma_[currentOrder_]*enorm;
  Tk_ = Scalar(currentOrder_+1)*Ek_;
  Est_ = Ek_;
  newOrder_ = currentOrder_;
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "currentOrder_ = " << currentOrder_ << std::endl;
    *out << "Ek_ = " << Ek_ << std::endl;
    *out << "Tk_ = " << Tk_ << std::endl;
    *out << "enorm = " << enorm << std::endl;
  }
  if (currentOrder_>1) {
    V_VpV(&*delta_,*xHistory_[currentOrder_],*ee_);
    Ekm1_ = sigma_[currentOrder_-1]*WRMSNorm(*errWtVec_,*delta_);
    Tkm1_ = currentOrder_*Ekm1_;
    if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
      *out << "Ekm1_ = " << Ekm1_ << std::endl;
      *out << "Tkm1_ = " << Tkm1_ << std::endl;
    }
    if (currentOrder_>2) {
      Vp_V(&*delta_,*xHistory_[currentOrder_-1]);
      Ekm2_ = sigma_[currentOrder_-2]*WRMSNorm(*errWtVec_,*delta_);
      Tkm2_ = (currentOrder_-1)*Ekm2_;
      if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
        *out << "Ekm2_ = " << Ekm2_ << std::endl;
        *out << "Tkm2_ = " << Tkm2_ << std::endl;
      }
      if (std::max(Tkm1_,Tkm2_)<=Tk_) {
        newOrder_--;
        Est_ = Ekm1_;
      }
    }
    else if (Tkm1_ <= Tkm1_Tk_safety_ * Tk_) {
      newOrder_--;
      Est_ = Ekm1_;
    }
  }
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "Est_ = " << Est_ << std::endl;
    *out << "newOrder_= " << newOrder_ << std::endl;
  }
  return(enorm);
}


template<class Scalar>
BDFstatusFlag ImplicitBDFStepper<Scalar>::rejectStep_()
{

  using Teuchos::as;

  // This routine puts its output in newTimeStep and newOrder

  // This routine changes the following variables:
  //    initialPhase, nef, psi, newTimeStep,
  //    newOrder, currentOrder_, currentTimeStep, dsDae.xHistory,
  //    dsDae.qHistory, nextTimePt, 
  //    currentTimeStepSum, nextTimePt

  // This routine reads but does not change the following variables:
  //    r_factor, r_safety, Est_, r_fudge_, r_min_, r_max_,
  //    minTimeStep_, maxTimeStep_, time, stopTime_ 

  // Only update the time step if we are NOT running constant stepsize.
  bool adjustStep = (!constantStepSize_);

  Scalar newTimeStep = hh_;
  Scalar rr = 1.0; // step size ratio = new step / old step
  nef_++;
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"rejectStep_");
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "adjustStep = " << adjustStep << std::endl;
    *out << "nef_ = " << nef_ << std::endl;
  }
  if (nef_ >= max_LET_fail_)  {
    std::cerr << "Rythmos_Stepper_ImplicitBDF::rejectStep_:  " 
         << "  Maximum number of local error test failures.  " << std::endl;
    return(REP_ERR_FAIL);
  }
  initialPhase_ = false;
  if (adjustStep) {
    if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
      *out << "initialPhase_ = " << initialPhase_ << std::endl;
    }
    restoreHistory_();
    // restore psi_
//    for (int i=1;i<=currentOrder_;++i)
//      psi_[i-1] = psi_[i] - hh_;

    if ((newtonConvergenceStatus_ < 0)) {
      /// 11/11/05 erkeite:  If the Newton solver fails, don't 
      // rely on the error estimate - it may be full of Nan's.
      rr = r_min_;
      newTimeStep = rr * hh_;

      if (nef_ > 2) {
        newOrder_ = 1;//consistent with block below.
      }
    } else {
      // 03/11/04 tscoffe:  Here is the block for choosing order & 
      // step-size when the local error test FAILS (but Newton 
      // succeeded). 
      if (nef_ == 1) { // first local error test failure
        rr = r_factor_*pow(r_safety_*Est_+r_fudge_,-1.0/(newOrder_+1.0));
        rr = std::max(r_min_,std::min(r_max_,rr));
        newTimeStep = rr * hh_;
      } else if (nef_ == 2) { // second failure
        rr = r_min_;
        newTimeStep = rr * hh_;
      } else if (nef_ > 2) { // third and later failures
        newOrder_ = 1;
        rr = r_min_;
        newTimeStep = rr * hh_;
      }
    }
    if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
      *out << "rr = " << rr << std::endl;
      *out << "newOrder_ = " << newOrder_ << std::endl;
    }
    currentOrder_ = newOrder_;
    if (numberOfSteps_ == 0) { // still first step
      psi_[0] = newTimeStep;
      Vt_S(&*xHistory_[1],rr);
      if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
        *out << "numberOfSteps_ == 0:" << std::endl;
        *out << "psi_[0] = " << psi_[0] << std::endl;
        *out << "xHistory_[1] = " << std::endl;
        xHistory_[1]->describe(*out,verbLevel);
      }
    }
  } else if (!adjustStep) {
    std::cerr << "Rythmos_Stepper_ImplicitBDF::rejectStep_:  "
         << "Warning: Local error test failed with constant step-size." << std::endl;
  }

  BDFstatusFlag return_status = PREDICT_AGAIN;

  // If the step needs to be adjusted:
  if (adjustStep) {
    newTimeStep = std::max(newTimeStep, minTimeStep_);
    newTimeStep = std::min(newTimeStep, maxTimeStep_);

    Scalar nextTimePt = time_ + newTimeStep;

    if (nextTimePt > stopTime_) {
      nextTimePt  = stopTime_;
      newTimeStep = stopTime_ - time_;
    }

    hh_ = newTimeStep;
  } else { // if time step is constant for this step:
    Scalar nextTimePt = time_ + hh_;

    if (nextTimePt > stopTime_) {
      nextTimePt      = stopTime_;
      hh_ = stopTime_ - time_;
    }
    return_status = CONTINUE_ANYWAY;
  }
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "hh_ = " << hh_ << std::endl;
  }
  return(return_status);
}


template<class Scalar>
void ImplicitBDFStepper<Scalar>::completeStep_()
{

  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  numberOfSteps_ ++;
  nef_ = 0;
  time_ += hh_;
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"completeStep_");

  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "numberOfSteps_ = " << numberOfSteps_ << std::endl;
    *out << "nef_ = " << nef_ << std::endl;
    *out << "time_ = " << time_ << std::endl;
  }
  
  // Only update the time step if we are NOT running constant stepsize.
  bool adjustStep = (!constantStepSize_);

  Scalar newTimeStep = hh_;
  Scalar rr = ST::one(); // step size ratio = new step / old step
  // 03/11/04 tscoffe:  Here is the block for choosing order & step-size when
  // the local error test PASSES (and Newton succeeded). 
  int orderDiff = currentOrder_ - usedOrder_;
  usedOrder_ = currentOrder_;
  usedStep_ = hh_;
  if ((newOrder_ == currentOrder_-1) || (currentOrder_ == maxOrder_)) {
    // If we reduced our order or reached std::max order then move to the next phase
    // of integration where we don't automatically double the step-size and
    // increase the order.
    initialPhase_ = false;
  }
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "initialPhase_ = " << initialPhase_ << std::endl;
  }
  if (initialPhase_) {
    currentOrder_++;
    newTimeStep = h_phase0_incr_ * hh_;
    if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
      *out << "currentOrder_ = " << currentOrder_ << std::endl;
      *out << "newTimeStep = " << newTimeStep << std::endl;
    }
  } else { // not in the initial phase of integration
    BDFactionFlag action = ACTION_UNSET;
    if (newOrder_ == currentOrder_-1) {
      action = ACTION_LOWER;
    } else if (newOrder_ == maxOrder_) {
      action = ACTION_MAINTAIN;
    } else if ((currentOrder_+1>=nscsco_) || (orderDiff == 1)) {
      // If we just raised the order last time then we won't raise it again
      // until we've taken currentOrder_+1 steps at order currentOrder_.
      action = ACTION_MAINTAIN;
    } else { // consider changing the order 
      V_StVpStV(&*delta_,ST::one(),*ee_,Scalar(-ST::one()),*xHistory_[currentOrder_+1]);
      Tkp1_ = WRMSNorm(*errWtVec_,*delta_);
      Ekp1_ = Tkp1_/(currentOrder_+2);
      if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
        *out << "delta_ = " << std::endl;
        delta_->describe(*out,verbLevel);
        *out << "Tkp1_ = ||delta_||_WRMS = " << Tkp1_ << std::endl;
        *out << "Ekp1_ = " << Ekp1_ << std::endl;
      }
      if (currentOrder_ == 1) {
        if (Tkp1_ >= Tkp1_Tk_safety_ * Tk_) {
          action = ACTION_MAINTAIN;
        } else {
          action = ACTION_RAISE;
        }
      } else {
        if (Tkm1_ <= std::min(Tk_,Tkp1_)) {
          action = ACTION_LOWER;
        } else if (Tkp1_ >= Tk_) {
          action = ACTION_MAINTAIN;
        } else {
          action = ACTION_RAISE;
        }
      }
    }
    if (action == ACTION_RAISE) {
      currentOrder_++;
      Est_ = Ekp1_;
    } else if (action == ACTION_LOWER) {
      currentOrder_--;
      Est_ = Ekm1_;
    }
    newTimeStep = hh_;
    rr = pow(r_safety_*Est_+r_fudge_,-1.0/(currentOrder_+1.0));
    if (rr >= r_hincr_test_) {
      rr = r_hincr_;
      newTimeStep = rr*hh_;
    } else if (rr <= 1) {
      rr = std::max(r_min_,std::min(r_max_,rr));
      newTimeStep = rr*hh_;
    }
    if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
      *out << "Est_ = " << Est_ << std::endl;
      *out << "currentOrder_ = " << currentOrder_ << std::endl;
      *out << "rr  = " << rr << std::endl;
      *out << "newTimeStep = " << newTimeStep << std::endl;
    }
  }
  // 03/22/04 tscoffe:  Note that updating the history has nothing to do with
  // the step-size and everything to do with the newton correction vectors.
  updateHistory_();

  // 12/01/05 tscoffe:  This is necessary to avoid currentTimeStep == 0 right
  // before a breakpoint.  So I'm checking to see if time is identically
  // equal to stopTime_, in which case we are right before a breakpoint and we
  // should not adjust currentStepSize because that would result in
  // currentStepSize == 0.
  if (time_ < stopTime_) {
    // If the step needs to be adjusted:
    if (adjustStep) {
      newTimeStep = std::max(newTimeStep, minTimeStep_);
      newTimeStep = std::min(newTimeStep, maxTimeStep_);

      Scalar nextTimePt = time_ + newTimeStep;

      if (nextTimePt > stopTime_) {
        nextTimePt  = stopTime_;
        newTimeStep = stopTime_ - time_;
      }

      hh_ = newTimeStep;
    } else { // if time step is constant for this step:
      Scalar nextTimePt = time_ + hh_;

      if (nextTimePt > stopTime_) {
        nextTimePt      = stopTime_;
        hh_ = stopTime_ - time_;
      }
    }
  }
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "hh_ = " << hh_ << std::endl;
  }
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::setStepControlData(const StepperBase<Scalar> & stepper)
{
  const ImplicitBDFStepper<Scalar>& bdfstepper = Teuchos::dyn_cast<const ImplicitBDFStepper<Scalar> >(stepper);
  currentOrder_ = bdfstepper.getOrder();
}


} // namespace Rythmos


#endif //Rythmos_IMPLICITBDF_STEPPER_H
