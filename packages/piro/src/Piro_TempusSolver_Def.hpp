// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#include "Piro_TempusSolver.hpp"

#include "Piro_ObserverToTempusIntegrationObserverAdapter.hpp"
#include "Piro_ValidPiroParameters.hpp"
#include "Piro_MatrixFreeDecorator.hpp"

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Assert.hpp"

#include "Thyra_DefaultAddedLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_DefaultModelEvaluatorWithSolveFactory.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

#include "Piro_InvertMassMatrixDecorator.hpp"

#ifdef HAVE_PIRO_IFPACK2
#include "Thyra_Ifpack2PreconditionerFactory.hpp"
#include "Tpetra_CrsMatrix.hpp"
#endif

#ifdef HAVE_PIRO_MUELU
#include <Thyra_MueLuPreconditionerFactory.hpp>
#include "Stratimikos_MueLuHelpers.hpp"
#endif

#ifdef HAVE_PIRO_NOX
#  include "Thyra_NonlinearSolver_NOX.hpp"
#endif

//#define DEBUG_OUTPUT

#include <string>
#include <stdexcept>
#include <iostream>

template <typename Scalar>
Piro::TempusSolver<Scalar>::TempusSolver(
    const Teuchos::RCP<Teuchos::ParameterList> &appParams,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &in_model,
    const Teuchos::RCP<Piro::ObserverBase<Scalar> > &piroObserver):
  TransientSolver<Scalar>(in_model), 
  out_(Teuchos::VerboseObjectBase::getDefaultOStream()),
  isInitialized_(false),
  piroObserver_(piroObserver),
  supports_x_dotdot_(false),
  initial_state_reset_{false}
{
  std::string sens_method_string = appParams->get("Sensitivity Method","None");
  this->setSensitivityMethod(sens_method_string); 
  sens_method_ = this->getSensitivityMethod(); 
  std::string jacobianSource = appParams->get("Jacobian Operator", "Have Jacobian");
  if (jacobianSource == "Matrix-Free") {
    Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > model;
    if (appParams->isParameter("Matrix-Free Perturbation")) {
      model = Teuchos::rcp(new Piro::MatrixFreeDecorator<Scalar>(in_model,
                           appParams->get<double>("Matrix-Free Perturbation")));
    }
    else model = Teuchos::rcp(new Piro::MatrixFreeDecorator<Scalar>(in_model));
    initialize(appParams, model);
  }
  else {
    initialize(appParams, in_model);
  }
}

template <typename Scalar>
void Piro::TempusSolver<Scalar>::initialize(
    const Teuchos::RCP<Teuchos::ParameterList> &appParams,
    const Teuchos::RCP< Thyra::ModelEvaluator<Scalar> > &in_model) 
{
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  model_ = in_model;  
  num_p_ = in_model->Np();
  num_g_ = in_model->Ng();

  //
  *out_ << "\nA) Get the base parameter list ...\n";
  //

  if (appParams->isSublist("Tempus")) {

    RCP<Teuchos::ParameterList> tempusPL = sublist(appParams, "Tempus", true);
    abort_on_failure_ = tempusPL->get<bool>("Abort on Failure", true); 

    RCP<Teuchos::ParameterList> integratorPL = sublist(tempusPL, "Tempus Integrator", true);
    //IKT, 10/31/16, FIXME: currently there is no Verbosity Sublist in Tempus, but
    //Curt will add this at some point.  When this option is added, set Verbosity
    //based on that sublist, rather than hard-coding it here.
    solnVerbLevel_ = Teuchos::VERB_DEFAULT;

    RCP<Teuchos::ParameterList> timeStepControlPL = Teuchos::null; 
    RCP<Teuchos::ParameterList> albTimeStepControlPL = Teuchos::null; 
    if (tempusPL->isSublist("Albany Time Step Control Options")) {
      *out_ << "\n    Using 'Albany Time Step Control Options'.\n";
      abort_on_fail_at_min_dt_ = true; 
      albTimeStepControlPL = sublist(tempusPL, "Albany Time Step Control Options"); 
      if (integratorPL->isSublist("Time Step Control")) {
        TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, 
            "\n Error!  You are attempting to specify 'Albany Time Step Control Options' and 'Time Step Control Strategy' \n "
             << "parameter lists.  Please pick one of these to use, and re-run.\n"); 
      } 
      if (abort_on_failure_ == true) {
        TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, 
            "\n Error!  'Abort on Failure = true' is invalid option when using 'Albany Time Step Control Options'\n"
            << "Please re-run with 'Abort on Failure = false' or use Tempus Time Step Control.\n"); 
      }
      abort_on_failure_ = false;  
      t_initial_ = albTimeStepControlPL->get<Scalar>("Initial Time", 0.0);
      t_final_ = albTimeStepControlPL->get<Scalar>("Final Time");
      Scalar dt_initial; 
      dt_initial = albTimeStepControlPL->get<Scalar>("Initial Time Step");
      Scalar dt_min = albTimeStepControlPL->get<Scalar>("Minimum Time Step", dt_initial);
      Scalar dt_max = albTimeStepControlPL->get<Scalar>("Maximum Time Step", dt_initial);
      Scalar reduc_factor = albTimeStepControlPL->get<Scalar>("Reduction Factor", 1.0);
      Scalar ampl_factor = albTimeStepControlPL->get<Scalar>("Amplification Factor", 1.0);
      timeStepControlPL = sublist(integratorPL, "Time Step Control", false);
      timeStepControlPL->set<Scalar>("Initial Time", t_initial_); 
      timeStepControlPL->set<Scalar>("Final Time", t_final_); 
      timeStepControlPL->set<Scalar>("Initial Time Step", dt_initial); 
      timeStepControlPL->set<Scalar>("Minimum Time Step", dt_min); 
      timeStepControlPL->set<Scalar>("Maximum Time Step", dt_max); 
      timeStepControlPL->set<int>("Initial Time Index", 0); 
      timeStepControlPL->set<int>("Final Time Index", 1.0e6);
      timeStepControlPL->set<std::string>("Integrator Step Type", "Variable"); 
      RCP<Teuchos::ParameterList> timeStepControlStrategyPL = sublist(timeStepControlPL, "Time Step Control Strategy"); 
      timeStepControlStrategyPL->set<std::string>("Time Step Control Strategy List", "basic_vs"); 
      RCP<Teuchos::ParameterList> basic_vs_PL = sublist(timeStepControlStrategyPL, "basic_vs"); 
      basic_vs_PL->set<std::string>("Name", "Basic VS"); 
      basic_vs_PL->set<Scalar>("Reduction Factor", reduc_factor); 
      basic_vs_PL->set<Scalar>("Amplification Factor", ampl_factor); 
      //The following options are such that dt is reduced if solve fails,
      //and amplified if solve is successful, to be consistent with what is done in Albany
      //for quasistatics, and for Schwarz. 
      basic_vs_PL->set<Scalar>("Minimum Value Monitoring Function", 1.0e20); 
      basic_vs_PL->set<Scalar>("Maximum Value Monitoring Function", 1.0e20); 
    }
    else { 
      *out_ << "\n    Using Tempus 'Time Step Control'.\n";
      RCP<Teuchos::ParameterList> timeStepControlPL = sublist(integratorPL, "Time Step Control", true);
      t_initial_ = timeStepControlPL->get<Scalar>("Initial Time", 0.0);
      t_final_ = timeStepControlPL->get<Scalar>("Final Time", t_initial_);
    }
    //*out_ << "tempusPL = " << *tempusPL << "\n";
    RCP<Teuchos::ParameterList> stepperPL = sublist(tempusPL, "Tempus Stepper", true);
    //*out_ << "stepperPL = " << *stepperPL << "\n";
    const std::string stepperType = stepperPL->get<std::string>("Stepper Type", "Backward Euler");
    //*out_ << "Stepper Type = " << stepperType << "\n";

    //
    // *out_ << "\nB) Create the Stratimikos linear solver factory ...\n";
    //
    // This is the linear solve strategy that will be used to solve for the
    // linear system with the W.
    //
    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;

#ifdef HAVE_PIRO_IFPACK2
    typedef Thyra::PreconditionerFactoryBase<double> Base;
    typedef Thyra::Ifpack2PreconditionerFactory<Tpetra::CrsMatrix<double> > Impl;
    linearSolverBuilder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), "Ifpack2");
#endif
#ifdef HAVE_PIRO_MUELU
    Stratimikos::enableMueLu(linearSolverBuilder);
#endif

    linearSolverBuilder.setParameterList(sublist(tempusPL, "Stratimikos", true));
    tempusPL->validateParameters(*getValidTempusParameters(),0);
    RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory = createLinearSolveStrategy(linearSolverBuilder);

    //
    *out_ << "\nC) Create and initalize the forward model ...\n";

    //
    // C.1) Create the underlying Thyra::ModelEvaluator
    // already constructed as "model". Decorate if needed.
    // IKT, 12/8/16: it may be necessary to expand the list of conditions
    // below, as more explicit schemes get added to Tempus
    // Explicit time-integrators for 1st order ODEs 
    if (
      stepperType == "Forward Euler" ||
      stepperType == "RK Forward Euler" ||
      stepperType == "RK1" || 
      stepperType == "RK Explicit 4 Stage" ||
      stepperType == "RK Explicit 3/8 Rule" ||
      stepperType == "RK Explicit 4 Stage 3rd order by Runge" ||
      stepperType == "RK Explicit 5 Stage 3rd order by Kinnmark and Gray"||
      stepperType == "RK Explicit 3 Stage 3rd order" ||
      stepperType == "RK Explicit 3 Stage 3rd order TVD" ||
      stepperType == "RK Explicit 3 Stage 3rd order by Heun" ||
      stepperType == "RK Explicit 2 Stage 2nd order by Runge" ||
      stepperType == "RK Explicit Midpoint" || 
      stepperType == "RK Explicit Trapezoidal" ||
      stepperType == "Heuns Method" ||
      stepperType == "Bogacki-Shampine 3(2) Pair" ||
      stepperType == "SSPERK22" || 
      stepperType == "SSPRK2" ||
      stepperType == "SSPERK33" || 
      stepperType == "SSPRK3" ||
      stepperType == "SSPERK54" ||
      stepperType == "General ERK" ) {

      Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > origModel = model_;
      model_ = Teuchos::rcp(new Piro::InvertMassMatrixDecorator<Scalar>(
      sublist(tempusPL,"Stratimikos", true), origModel, tempusPL->get("Constant Mass Matrix", false), tempusPL->get("Lump Mass Matrix", false),false));
    }

    //Explicit time-integrators for 2nd order ODEs
    //IKT, FIXME: fill this in as more explicit integrators for 2nd order ODEs are added to Tempus.
    else if (stepperType == "Newmark Explicit a-Form") {
      Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > origModel = model_;
      model_ = Teuchos::rcp(new Piro::InvertMassMatrixDecorator<Scalar>(
        sublist(tempusPL,"Stratimikos", true), origModel, tempusPL->get("Constant Mass Matrix", false), tempusPL->get("Lump Mass Matrix", false),true));
    }
    // C.2) Create the Thyra-wrapped ModelEvaluator

    thyraModel_ = rcp(new Thyra::DefaultModelEvaluatorWithSolveFactory<Scalar>(model_, lowsFactory));

    const RCP<const Thyra::VectorSpaceBase<double> > x_space = thyraModel_->get_x_space();

    //
    *out_ << "\nD) Create the stepper and integrator for the forward problem ...\n";

    //Create Tempus integrator with observer using tempusPL, model_ and sensitivity method
    piroTempusIntegrator_ = Teuchos::rcp(new Piro::TempusIntegrator<Scalar>(tempusPL, model_, sens_method_));
    this->setPiroTempusIntegrator(piroTempusIntegrator_);  

    //Get stepper from integrator
    fwdStateStepper_ = piroTempusIntegrator_->getStepper();

    //Set observer
    supports_x_dotdot_ = model_->createInArgs().supports(Thyra::ModelEvaluatorBase::IN_ARG_x_dot_dot);
    setObserver();  

  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
        appParams->isSublist("Tempus"),
        Teuchos::Exceptions::InvalidParameter,
        "\n Error! Piro::TempusSolver: must have Tempus sublist ");

  }

  isInitialized_ = true;
}

template <typename Scalar>
Piro::TempusSolver<Scalar>::TempusSolver(
    const Teuchos::RCP<Piro::TempusIntegrator<Scalar> > &stateIntegrator,
    const Teuchos::RCP<Tempus::Stepper<Scalar> > &stateStepper,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > &timeStepSolver,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &underlyingModel,
    Scalar finalTime,
    const std::string sens_method_string, 
    Teuchos::EVerbosityLevel verbosityLevel) :
  TransientSolver<Scalar>(underlyingModel), 
  piroTempusIntegrator_(stateIntegrator),
  fwdStateStepper_(stateStepper),
  fwdTimeStepSolver_(timeStepSolver),
  model_(underlyingModel),
  t_initial_(0.0),
  t_final_(finalTime),
  num_p_(model_->Np()),
  num_g_(model_->Ng()),
  out_(Teuchos::VerboseObjectBase::getDefaultOStream()),
  solnVerbLevel_(verbosityLevel),
  isInitialized_(true)
{
  if (fwdStateStepper_->getModel() != underlyingModel) {
    fwdStateStepper_->setModel(underlyingModel);
  }
  this->setSensitivityMethod(sens_method_string); 
  sens_method_ = this->getSensitivityMethod(); 
  this->setPiroTempusIntegrator(piroTempusIntegrator_);  
}

template <typename Scalar>
Piro::TempusSolver<Scalar>::TempusSolver(
    const Teuchos::RCP<Piro::TempusIntegrator<Scalar> > &stateIntegrator,
    const Teuchos::RCP<Tempus::Stepper<Scalar> > &stateStepper,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > &timeStepSolver,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &underlyingModel,
    Scalar initialTime,
    Scalar finalTime,
    const std::string sens_method_string, 
    Teuchos::EVerbosityLevel verbosityLevel) :
  TransientSolver<Scalar>(underlyingModel), 
  piroTempusIntegrator_(stateIntegrator),
  fwdStateStepper_(stateStepper),
  fwdTimeStepSolver_(timeStepSolver),
  model_(underlyingModel),
  t_initial_(initialTime),
  t_final_(finalTime),
  num_p_(model_->Np()),
  num_g_(model_->Ng()),
  out_(Teuchos::VerboseObjectBase::getDefaultOStream()),
  solnVerbLevel_(verbosityLevel),
  isInitialized_(true)
{
  if (fwdStateStepper_->getModel() != underlyingModel) {
    fwdStateStepper_->setModel(underlyingModel);
  }
  this->setSensitivityMethod(sens_method_string); 
  sens_method_ = this->getSensitivityMethod(); 
  this->setPiroTempusIntegrator(piroTempusIntegrator_);  
}

template <typename Scalar>
void Piro::TempusSolver<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const 
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Set initial time and initial condition 
  Thyra::ModelEvaluatorBase::InArgs<Scalar> state_ic = model_->getNominalValues();
  Teuchos::RCP<const Thyra::VectorBase<Scalar>> xinit, xdotinit, xdotdotinit; 
  if(t_initial_ > 0.0 && state_ic.supports(Thyra::ModelEvaluatorBase::IN_ARG_t)) {
    state_ic.set_t(t_initial_);
    //If initial state has not been reset, get the initial state from ME in args
    if (!initial_state_reset_) { 
      if (state_ic.supports(Thyra::ModelEvaluatorBase::IN_ARG_x)) {
        xinit = state_ic.get_x();
      }
      if (state_ic.supports(Thyra::ModelEvaluatorBase::IN_ARG_x_dot)) {
        xdotinit = state_ic.get_x_dot();
      }
      if (state_ic.supports(Thyra::ModelEvaluatorBase::IN_ARG_x_dot_dot)) {
        xdotdotinit = state_ic.get_x_dot_dot();
      }
      piroTempusIntegrator_->initializeSolutionHistory(t_initial_, xinit, xdotinit, xdotdotinit); 
      //Reset observer.  This is necessary for correct observation of solution
      //since initializeSolutionHistory modifies the solutionHistory object.
      setObserver();
#ifdef DEBUG_OUTPUT
      if (xinit != Teuchos::null) { 
        *out_ << "\n*** Piro::TempusSolver::evalModelImpl xinit at time = " << t_initial_ << " ***\n";
        Teuchos::Range1D range;
        RTOpPack::ConstSubVectorView<Scalar> xinitv;
        xinit->acquireDetachedView(range, &xinitv);
        auto xinita = xinitv.values();
        for (auto i = 0; i < xinita.size(); ++i) *out_ << xinita[i] << " ";
        *out_ << "\n*** Piro::TempusSolver::evalModelImpl xinit at time = " << t_initial_ << " ***\n";
      }
      if (xdotinit != Teuchos::null) { 
        *out_ << "\n*** Piro::TempusSolver::evalModelImpl xdotinit at time = " << t_initial_ << " ***\n";
        Teuchos::Range1D range;
        RTOpPack::ConstSubVectorView<Scalar> xdotinitv;
        xdotinit->acquireDetachedView(range, &xdotinitv);
        auto xdotinita = xdotinitv.values();
        for (auto i = 0; i < xdotinita.size(); ++i) *out_ << xdotinita[i] << " ";
        *out_ << "\n*** Piro::TempusSolver::evalModelImpl xdotinit at time = " << t_initial_ << " ***\n";
      }
      if (xdotdotinit != Teuchos::null) { 
        *out_ << "\n*** Piro::TempusSolver::evalModelImpl xdotdotinit at time = " << t_initial_ << " ***\n";
        Teuchos::Range1D range;
        RTOpPack::ConstSubVectorView<Scalar> xdotdotinitv;
        xdotdotinit->acquireDetachedView(range, &xdotdotinitv);
        auto xdotdotinita = xdotdotinitv.values();
        for (auto i = 0; i < xdotdotinita.size(); ++i) *out_ << xdotdotinita[i] << " ";
        *out_ << "\n*** Piro::TempusSolver::evalModelImpl xdotdotinit at time = " << t_initial_ << " ***\n";
      }
#endif
    }
  }
  
  // Set parameters as part of initial conditions
  for (int l = 0; l < num_p_; ++l) {
    auto p_in = inArgs.get_p(l); 
    if (Teuchos::nonnull(p_in)) {
      state_ic.set_p(l, p_in);
    }
  }

  //*out_ << "\nstate_ic:\n" << Teuchos::describe(state_ic, solnVerbLevel_);

  //
  *out_ << "\nE) Solve the forward problem ...\n";
  //
  RCP<const Thyra::VectorBase<Scalar> > finalSolution;
  RCP<const Tempus::SolutionState<Scalar> > solutionState;
  RCP<const Tempus::SolutionHistory<Scalar> > solutionHistory;
    
  *out_ << "T final requested: " << t_final_ << " \n";

  piroTempusIntegrator_->advanceTime(t_final_);
  double time = piroTempusIntegrator_->getTime();
  *out_ << "T final actual: " << time << "\n";
 
  Scalar diff = 0.0; 
  if (abs(t_final_) == 0) diff = abs(time-t_final_);
  else diff = abs(time-t_final_)/abs(t_final_);  
  if (diff > 1.0e-10) {
    if (abort_on_failure_ == true) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true,
        Teuchos::Exceptions::InvalidParameter,
        "\n Error! Piro::TempusSolver: time-integrator did not make it to final time " <<
        "specified in Input File.  Final time in input file is " << t_final_ <<
        ", whereas actual final time is " << time << ".  If you'd like to " <<
        "suppress this exception, run with 'Abort on Failure' set to 'false' in " << 
        "Tempus sublist.\n" );
    }
    else {
       *out_ << "\n WARNING: Piro::TempusSolver did not make it to final time, but "
            << "solver will not abort since you have specified 'Abort on Failure' = 'false'.\n"; 
    }
  }

  solutionHistory = piroTempusIntegrator_->getSolutionHistory();
  auto numStates = solutionHistory->getNumStates();
  solutionState = (*solutionHistory)[numStates-1];
  //Get final solution from solutionHistory.
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;
  Teuchos::RCP<const Thyra::VectorBase<Scalar>> x = solutionState->getX(); 
  Teuchos::RCP<const DMVPV> X = Teuchos::rcp_dynamic_cast<const DMVPV>(x);
  finalSolution = (sens_method_ == NONE) ? x : X->getMultiVector()->col(0);

  if (Teuchos::VERB_MEDIUM <= solnVerbLevel_) {
    *out_ << "Final Solution\n" << *finalSolution << "\n";
  }


  // As post-processing step, calculate responses at final solution
  Thyra::ModelEvaluatorBase::InArgs<Scalar> modelInArgs = model_->createInArgs();
  
  modelInArgs.set_x(finalSolution);
  for (int l=0; l < num_p_; ++l) { 
    auto p_in = inArgs.get_p(l); 
    modelInArgs.set_p(l, p_in);
  }
  //Set time to be final time at which the solve occurs (< t_final_ in the case we don't make it to t_final_).
  //IKT: get final time from solutionHistory workingSpace, which is different than how it is done in Piro::RythmosSolver class.
  //IKT, 11/1/16, FIXME? workingState pointer is null right now, so the following
  //code is commented out for now.  Use t_final_ and soln_dt in set_t instead for now.
  /*RCP<Tempus::SolutionState<Scalar> > workingState = solutionHistory->getWorkingState();
  const Scalar time = workingState->getTime();
  const Scalar dt   = workingState->getTimeStep();
  const Scalar t = time + dt;
  modelInArgs.set_t(t);*/
  const Scalar soln_dt = solutionState->getTimeStep();
  modelInArgs.set_t(t_final_ - soln_dt);

  //Calculate responses and sensitivities 
  this->evalConvergedModelResponsesAndSensitivities(modelInArgs, outArgs);

}


template <typename Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
Piro::TempusSolver<Scalar>::getValidTempusParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> validPL =
    Teuchos::rcp(new Teuchos::ParameterList("ValidTempusSolverParams"));
  validPL->sublist("Tempus", false, "");
  validPL->sublist("Albany Time Step Control Options", false, "");
  validPL->sublist("Albany Time Step Control Options", false, "").set<Scalar>("Initial Time", 0.0, "");
  validPL->sublist("Albany Time Step Control Options", false, "").set<Scalar>("Initial Time Step", 0.1, "");
  validPL->sublist("Albany Time Step Control Options", false, "").set<Scalar>("Minimum Time Step", 0.1, "");
  validPL->sublist("Albany Time Step Control Options", false, "").set<Scalar>("Maximum Time Step", 0.1, "");
  validPL->sublist("Albany Time Step Control Options", false, "").set<Scalar>("Final Time", 1.0, "");
  validPL->sublist("Albany Time Step Control Options", false, "").set<Scalar>("Reduction Factor", 1.0, "");
  validPL->sublist("Albany Time Step Control Options", false, "").set<Scalar>("Amplification Factor", 1.0, "");
  validPL->sublist("Stratimikos", false, "");
  validPL->sublist("NonLinear Solver", false, "");
  //validPL->set<std::string>("Verbosity Level", "", "");
  validPL->set<bool>("Lump Mass Matrix", false, "Boolean to tell code whether to lump mass matrix");
  validPL->set<bool>("Constant Mass Matrix", false, "Boolean to tell code if mass matrix is constant in time");
  validPL->set<bool>("Abort on Failure", true, "");
  validPL->set<std::string>("Integrator Name", "Tempus Integrator", "");
  validPL->sublist("Tempus Integrator", false, "");
  validPL->sublist("Tempus Stepper", false, "");
  validPL->sublist("Time Step Control", false, "");
  return validPL;
}

template <typename Scalar>
void Piro::TempusSolver<Scalar>::
addStepperFactory(const std::string & stepperName,const Teuchos::RCP<Piro::TempusStepperFactory<Scalar> > & factory)
{
  stepperFactories_[stepperName] = factory;
}

template <typename Scalar>
void Piro::TempusSolver<Scalar>::
addStepControlFactory(const std::string & stepControlName,
                      const Teuchos::RCP<Piro::TempusStepControlFactory<Scalar>> & step_control_strategy)
{
  stepControlFactories_[stepControlName] = step_control_strategy;
}

template <typename Scalar>
void Piro::TempusSolver<Scalar>::
setStartTime(const Scalar start_time)
{
  Teuchos::RCP<const Tempus::TimeStepControl<Scalar> > tsc_const = piroTempusIntegrator_->getTimeStepControl();
  Teuchos::RCP<Tempus::TimeStepControl<Scalar> > tsc = Teuchos::rcp_const_cast<Tempus::TimeStepControl<Scalar> >(tsc_const); 
  t_initial_ = start_time;  
  tsc->setInitTime(start_time);
} 

template <typename Scalar>
Scalar Piro::TempusSolver<Scalar>::
getStartTime() const
{
  Teuchos::RCP<const Tempus::TimeStepControl<Scalar> > tsc = piroTempusIntegrator_->getTimeStepControl();
  Scalar start_time = tsc->getInitTime(); 
  return start_time; 
} 

template <typename Scalar>
void Piro::TempusSolver<Scalar>::
setFinalTime(const Scalar final_time)
{
  Teuchos::RCP<const Tempus::TimeStepControl<Scalar> > tsc_const = piroTempusIntegrator_->getTimeStepControl();
  Teuchos::RCP<Tempus::TimeStepControl<Scalar> > tsc = Teuchos::rcp_const_cast<Tempus::TimeStepControl<Scalar> >(tsc_const); 
  t_final_ = final_time; 
  tsc->setFinalTime(final_time); 
} 

template <typename Scalar>
Scalar Piro::TempusSolver<Scalar>::
getFinalTime() const
{
  Teuchos::RCP<const Tempus::TimeStepControl<Scalar> > tsc = piroTempusIntegrator_->getTimeStepControl();
  Scalar final_time = tsc->getFinalTime(); 
  return final_time; 
} 

template <typename Scalar>
void Piro::TempusSolver<Scalar>::
setInitTimeStep(const Scalar init_time_step)
{
  Teuchos::RCP<const Tempus::TimeStepControl<Scalar> > tsc_const = piroTempusIntegrator_->getTimeStepControl();
  Teuchos::RCP<Tempus::TimeStepControl<Scalar> > tsc = Teuchos::rcp_const_cast<Tempus::TimeStepControl<Scalar> >(tsc_const); 
  tsc->setInitTimeStep(init_time_step); 
} 


template <typename Scalar>
Scalar Piro::TempusSolver<Scalar>::
getInitTimeStep() const
{
  Teuchos::RCP<const Tempus::TimeStepControl<Scalar> > tsc = piroTempusIntegrator_->getTimeStepControl();
  auto init_time_step = tsc->getInitTimeStep(); 
  return init_time_step; 
} 
template <typename Scalar>
void Piro::TempusSolver<Scalar>::
setObserver() const
{
  Teuchos::RCP<Tempus::IntegratorObserverBasic<Scalar> > observer = Teuchos::null;
  if (Teuchos::nonnull(piroObserver_)) {
    //Get solutionHistory from integrator
    const Teuchos::RCP<const Tempus::SolutionHistory<Scalar> > solutionHistory = piroTempusIntegrator_->getSolutionHistory();
    const Teuchos::RCP<const Tempus::TimeStepControl<Scalar> > timeStepControl = piroTempusIntegrator_->getTimeStepControl();
    //Create Tempus::IntegratorObserverBasic object
    observer = Teuchos::rcp(new ObserverToTempusIntegrationObserverAdapter<Scalar>(solutionHistory,
                                timeStepControl, piroObserver_, supports_x_dotdot_, abort_on_fail_at_min_dt_, sens_method_));
  }
  if (Teuchos::nonnull(observer)) {
    //Set observer in integrator
    piroTempusIntegrator_->clearObservers();
    piroTempusIntegrator_->setObserver(observer);
    fwdStateStepper_->initialize();
    //Reinitialize everything in integrator class, since we have changed the observer.
    piroTempusIntegrator_->initialize();
  }
}

template <typename Scalar>
void Piro::TempusSolver<Scalar>::
setInitialState(Scalar t0,
      Teuchos::RCP<Thyra::VectorBase<Scalar> > x0,
      Teuchos::RCP<Thyra::VectorBase<Scalar> > xdot0,
      Teuchos::RCP<Thyra::VectorBase<Scalar> > xdotdot0) 
{
   piroTempusIntegrator_->initializeSolutionHistory(t0, x0, xdot0, xdotdot0); 
   //Reset observer.  This is necessary for correct observation of solution
   //since initializeSolutionHistory modifies the solutionHistory object.
   setObserver();
   initial_state_reset_ = true;  
}

template <typename Scalar>
void Piro::TempusSolver<Scalar>::
setInitialGuess(Teuchos::RCP< const Thyra::VectorBase<Scalar> > initial_guess) 
{
   fwdStateStepper_->setInitialGuess(initial_guess); 
   fwdStateStepper_->initialize();
}

template <typename Scalar>
Teuchos::RCP<Tempus::SolutionHistory<Scalar> > Piro::TempusSolver<Scalar>::
getSolutionHistory() const
{
  Teuchos::RCP<const Tempus::SolutionHistory<Scalar> > soln_history_const = piroTempusIntegrator_->getSolutionHistory();
  Teuchos::RCP<Tempus::SolutionHistory<Scalar> > soln_history = Teuchos::rcp_const_cast<Tempus::SolutionHistory<Scalar> >(soln_history_const); 
  return soln_history;
}
  

template <typename Scalar>
Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > Piro::TempusSolver<Scalar>::
getSolver() const
{
  return fwdStateStepper_->getSolver(); 
}


template <typename Scalar>
Tempus::Status Piro::TempusSolver<Scalar>::
getTempusIntegratorStatus() const
{
  return piroTempusIntegrator_->getStatus(); 
}


template <typename Scalar>
Teuchos::RCP<Piro::TempusSolver<Scalar> >
Piro::tempusSolver(
    const Teuchos::RCP<Teuchos::ParameterList> &appParams,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &in_model,
    const Teuchos::RCP<Piro::ObserverBase<Scalar> > &piroObserver)
{
  Teuchos::RCP<Teuchos::FancyOStream> out_(Teuchos::VerboseObjectBase::getDefaultOStream());
  return Teuchos::rcp(new TempusSolver<Scalar>(appParams, in_model, piroObserver));
}


