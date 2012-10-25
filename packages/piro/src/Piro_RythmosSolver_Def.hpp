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

#include <cmath>

#include "Piro_RythmosSolver.hpp"
#include "Piro_ValidPiroParameters.hpp"

#include "Rythmos_BackwardEulerStepper.hpp"
#include "Rythmos_ExplicitRKStepper.hpp"
#include "Rythmos_ImplicitBDFStepper.hpp"
#include "Rythmos_SimpleIntegrationControlStrategy.hpp"
#include "Rythmos_RampingIntegrationControlStrategy.hpp"
#include "Rythmos_ForwardSensitivityStepper.hpp"
#include "Rythmos_ImplicitBDFStepperRampingStepControl.hpp"
#include "Rythmos_StepperAsModelEvaluator.hpp"
#include "Rythmos_CompositeIntegrationObserver.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

#include "Rythmos_IntegratorBuilder.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Assert.hpp"
#include "Piro_RythmosNOX_RowSumUpdater.hpp"

#ifdef Piro_ENABLE_NOX
#  include "Thyra_NonlinearSolver_NOX.hpp"
#endif

#include "Thyra_ScaledModelEvaluator.hpp"

template <typename Scalar>
Piro::RythmosSolver<Scalar>::RythmosSolver(Teuchos::RCP<Teuchos::ParameterList> in_appParams,
                          Teuchos::RCP< Thyra::ModelEvaluatorDefaultBase<Scalar> > in_model,
                          Teuchos::RCP<Rythmos::IntegrationObserverBase<Scalar> > in_observer) :
  appParams(in_appParams),
  model(in_model),
  observer(in_observer)
{
  // For dumping default parameters from Rythmos
  {
    //Rythmos::IntegratorBuilder<double> b;
    //std::cout << *(b.getValidParameters()) << std::endl;
    //Teuchos::writeParameterListToXmlFile(*b.getValidParameters(), "sample.xml");
  }

  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  out = Teuchos::VerboseObjectBase::getDefaultOStream();

  num_p = model->createInArgs().Np();
  num_g = model->createOutArgs().Ng();

//   TEUCHOS_TEST_FOR_EXCEPTION(num_p > 1, Teuchos::Exceptions::InvalidParameter,
//                      std::endl << "Error in Piro::RythmosSolver " <<
//                      "Not Implemented for Np>1 : " << num_p << std::endl);
//   TEUCHOS_TEST_FOR_EXCEPTION(num_g > 1, Teuchos::Exceptions::InvalidParameter,
//                      std::endl << "Error in Piro::RythmosSolver " <<
//                      "Not Implemented for Ng>1 : " << num_g << std::endl);

  //
  *out << "\nA) Get the base parameter list ...\n";
  //

  RCP<Teuchos::ParameterList> rythmosPL = sublist(appParams, "Rythmos", true);
  rythmosPL->validateParameters(*getValidRythmosParameters(),0);

  {
    const std::string verbosity = rythmosPL->get("Verbosity Level", "VERB_DEFAULT");
    solnVerbLevel = Teuchos::VERB_DEFAULT;
    if      (verbosity == "VERB_NONE")    solnVerbLevel = Teuchos::VERB_NONE;
    else if (verbosity == "VERB_LOW")     solnVerbLevel = Teuchos::VERB_LOW;
    else if (verbosity == "VERB_MEDIUM")  solnVerbLevel = Teuchos::VERB_MEDIUM;
    else if (verbosity == "VERB_HIGH")    solnVerbLevel = Teuchos::VERB_HIGH;
    else if (verbosity == "VERB_EXTREME") solnVerbLevel = Teuchos::VERB_EXTREME;
  }

  t_final = rythmosPL->get("Final Time", 0.1);
  
  const std::string stepperType = rythmosPL->get("Stepper Type", "Backward Euler");
  
  //
  *out << "\nC) Create and initalize the forward model ...\n";
  //
      
  *out << "\nD) Create the stepper and integrator for the forward problem ...\n";
  //
  
  if (rythmosPL->get<std::string>("Nonlinear Solver Type") == "Rythmos") {
    Teuchos::RCP<Rythmos::TimeStepNonlinearSolver<double> > rythmosTimeStepSolver = 
      Rythmos::timeStepNonlinearSolver<double>();
    if (rythmosPL->getEntryPtr("NonLinear Solver")) {
      RCP<Teuchos::ParameterList> nonlinePL =
	sublist(rythmosPL, "NonLinear Solver", true);
      rythmosTimeStepSolver->setParameterList(nonlinePL);
    }
    fwdTimeStepSolver = rythmosTimeStepSolver;
  }
  else if (rythmosPL->get<std::string>("Nonlinear Solver Type") == "NOX") {
#ifdef Piro_ENABLE_NOX
    Teuchos::RCP<Thyra::NOXNonlinearSolver> nox_solver =  Teuchos::rcp(new Thyra::NOXNonlinearSolver);
    Teuchos::RCP<Teuchos::ParameterList> nox_params = Teuchos::rcp(new Teuchos::ParameterList);
    *nox_params = appParams->sublist("NOX");
    nox_solver->setParameterList(nox_params);
    fwdTimeStepSolver = nox_solver;
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,"Requested NOX solver for a Rythmos Transient solve, Trilinos was not built with NOX enabled.  Please rebuild Trilinos or use the native Rythmos nonlinear solver.");
#endif
    
  }
  
  if (stepperType == "Backward Euler") {
    fwdStateStepper = Rythmos::backwardEulerStepper<Scalar> (model, fwdTimeStepSolver);
    fwdStateStepper->setParameterList(sublist(rythmosPL, "Rythmos Stepper", true));
  }
  else if (stepperType == "Explicit RK") {
    fwdStateStepper = Rythmos::explicitRKStepper<Scalar>(model);
    fwdStateStepper->setParameterList(sublist(rythmosPL, "Rythmos Stepper", true));
  }
  else if (stepperType == "BDF") {
    Teuchos::RCP<Teuchos::ParameterList> BDFparams = 
      Teuchos::sublist(rythmosPL, "Rythmos Stepper", true);
    Teuchos::RCP<Teuchos::ParameterList> BDFStepControlPL =
      Teuchos::sublist(BDFparams,"Step Control Settings");
    
    fwdStateStepper = Teuchos::rcp( new Rythmos::ImplicitBDFStepper<Scalar>(model,fwdTimeStepSolver,BDFparams) );
    fwdStateStepper->setInitialCondition(model->getNominalValues());
    
  }
  else 
    TEUCHOS_TEST_FOR_EXCEPTION( true, Teuchos::Exceptions::InvalidParameter,
				std::endl << "Error! Piro::Epetra::RythmosSolver: Invalid Steper Type: "
				<< stepperType << std::endl);
  
  // Step control strategy
  {
    // If the stepper can accept a step control strategy, then attempt to build one.
    RCP<Rythmos::StepControlStrategyAcceptingStepperBase<Scalar> > scsa_stepper = 
      Teuchos::rcp_dynamic_cast<Rythmos::StepControlStrategyAcceptingStepperBase<Scalar> >(fwdStateStepper);

    if ( nonnull(scsa_stepper) ) {
      std::string step_control_strategy = rythmosPL->get("Step Control Strategy Type", "None");

      if (step_control_strategy == "None") {
	// don't do anything, stepper will build default
      }
      else if (step_control_strategy == "ImplicitBDFRamping") {
	
	const RCP<Rythmos::ImplicitBDFStepperRampingStepControl<Scalar> > rscs = 
	  rcp(new Rythmos::ImplicitBDFStepperRampingStepControl<Scalar>);
	
	const RCP<ParameterList> p = parameterList(rythmosPL->sublist("Rythmos Step Control Strategy"));
	rscs->setParameterList(p);

	scsa_stepper->setStepControlStrategy(rscs);
      }
      else {
	TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,"Error! Piro::Epetra::RythmosSolver: Invalid step control strategy type: " << step_control_strategy << std::endl);
      }
      
    }

  }
  
  {
    RCP<Teuchos::ParameterList>
      integrationControlPL = sublist(rythmosPL, "Rythmos Integration Control", true);
    
    RCP<Rythmos::DefaultIntegrator<Scalar> > defaultIntegrator;
    
    if (rythmosPL->get("Rythmos Integration Control Strategy", "Simple") == "Simple") {
      defaultIntegrator = Rythmos::controlledDefaultIntegrator<Scalar>(Rythmos::simpleIntegrationControlStrategy<Scalar>(integrationControlPL));
    }
    else if(rythmosPL->get<std::string>("Rythmos Integration Control Strategy") == "Ramping") {
      defaultIntegrator = Rythmos::controlledDefaultIntegrator<Scalar>(Rythmos::rampingIntegrationControlStrategy<Scalar>(integrationControlPL));  
    }

    fwdStateIntegrator = defaultIntegrator;
  }

  fwdStateIntegrator->setParameterList(sublist(rythmosPL, "Rythmos Integrator", true));
  
  if (observer != Teuchos::null) 
    fwdStateIntegrator->setIntegrationObserver(observer);
  
}

template <typename Scalar>
Piro::RythmosSolver<Scalar>::~RythmosSolver()
{
}

template <typename Scalar>
Teuchos::RCP<const Rythmos::IntegratorBase<Scalar> > 
Piro::RythmosSolver<Scalar>::getRythmosIntegrator() const
{
  return fwdStateIntegrator;
}

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::RythmosSolver<Scalar>::get_p_space(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l >= num_p || l < 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error in Piro::RythmosSolver::get_p_map():  " <<
                     "Invalid parameter index l = " <<
                     l << std::endl);
  return model->get_p_space(l);
}

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::RythmosSolver<Scalar>::get_g_space(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(j > num_g || j < 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error in Piro::RythmosSolver::get_g_map():  " <<
                     "Invalid response index j = " <<
                     j << std::endl);

  if      (j < num_g) return model->get_g_space(j);
  else return model->get_x_space(); // j == num_g
}

template<typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::RythmosSolver<Scalar>::getNominalValues() const
{
  Thyra::ModelEvaluatorBase::InArgs<Scalar> result = this->createInArgs();
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> modelNominalValues = model->getNominalValues();
  for (int l = 0; l < num_p; ++l) {
    result.set_p(l, modelNominalValues.get_p(l));
  }
  return result;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar> Piro::RythmosSolver<Scalar>::createInArgs() const
{
  Thyra::ModelEvaluatorBase::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(num_p);
  return inArgs;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar> Piro::RythmosSolver<Scalar>::createOutArgsImpl() const
{
  Thyra::ModelEvaluatorBase::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());

  // Ng is 1 bigger then model-Ng so that the solution vector can be an outarg
  outArgs.set_Np_Ng(num_p, num_g+1);

  Thyra::ModelEvaluatorBase::OutArgs<Scalar> model_outargs = model->createOutArgs();
  for (int i=0; i<num_g; i++)
    for (int j=0; j<num_p; j++)
      if (!model_outargs.supports( Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, i, j).none())
        outArgs.setSupports( Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, i, j,
             Thyra::ModelEvaluatorBase::DerivativeSupport( Thyra::ModelEvaluatorBase::DERIV_MV_BY_COL));

  return outArgs;
}

template <typename Scalar>
void Piro::RythmosSolver<Scalar>::evalModelImpl(
       const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
       const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs ) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Parse InArgs

  RCP<const Thyra::VectorBase<Scalar> > p_in;
  if (num_p > 0) p_in = inArgs.get_p(0);

  // Parse OutArgs: always 1 extra
  RCP< Thyra::VectorBase<Scalar> > g_out; 
  if (num_g > 0) g_out = outArgs.get_g(0); 
  RCP< Thyra::VectorBase<Scalar> > gx_out = outArgs.get_g(num_g); 

  // Parse out-args for sensitivity calculation
  RCP< Thyra::MultiVectorBase<Scalar> > dgdp_out;
  if (num_p>0 && num_g>0)
    dgdp_out = outArgs.get_DgDp(0,0).getMultiVector();

  RCP<const Thyra::VectorBase<Scalar> > finalSolution;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> state_ic = model->getNominalValues();

  // Set paramters p_in as part of initial conditions
  if (num_p > 0) state_ic.set_p(0,p_in);

  //*out << "\nstate_ic:\n" << Teuchos::describe(state_ic,Teuchos::VERB_NONE);
  *out << "\nstate_ic:\n" << Teuchos::describe(state_ic,solnVerbLevel);

  if (dgdp_out == Teuchos::null) {
      //
      *out << "\nE) Solve the forward problem ...\n";
      //
  
      fwdStateStepper->setInitialCondition(state_ic);
      fwdStateIntegrator->setStepper(fwdStateStepper, t_final, true);
  
      Teuchos::Array<RCP<const Thyra::VectorBase<Scalar> > > x_final_array;
      fwdStateIntegrator->getFwdPoints(
        Teuchos::tuple<Scalar>(t_final), &x_final_array, NULL, NULL
        );
      finalSolution = x_final_array[0];
  
      if (Teuchos::VERB_MEDIUM <= solnVerbLevel)
         std::cout << "Final Solution\n" << *finalSolution << std::endl;

     // As post-processing step, calc responses at final solution
     Thyra::ModelEvaluatorBase::InArgs<Scalar>  model_inargs = model->createInArgs();
     Thyra::ModelEvaluatorBase::OutArgs<Scalar> model_outargs = model->createOutArgs();
     model_inargs.set_x(finalSolution);
     if (num_p > 0)  model_inargs.set_p(0, p_in);
     if (g_out != Teuchos::null) {
       Thyra::put_scalar(0.0,g_out.ptr());
       model_outargs.set_g(0, g_out);
     }

     model->evalModel(model_inargs, model_outargs);

   }
   else {//Doing sensitivities
      //
      *out << "\nE) Solve the forward problem with Sensitivities...\n";
      //

      RCP<Rythmos::ForwardSensitivityStepper<Scalar> > stateAndSensStepper =
        Rythmos::forwardSensitivityStepper<Scalar>();
      stateAndSensStepper->initializeSyncedSteppers(
          model, 0, model->getNominalValues(),
          fwdStateStepper, fwdTimeStepSolver);

      //
      // Set the initial condition for the state and forward sensitivities
      //

      RCP< Thyra::VectorBase<Scalar> > s_bar_init
        = createMember(stateAndSensStepper->getFwdSensModel()->get_x_space());
      assign( s_bar_init.ptr(), 0.0 );
      RCP< Thyra::VectorBase<Scalar> > s_bar_dot_init
        = createMember(stateAndSensStepper->getFwdSensModel()->get_x_space());
      assign( s_bar_dot_init.ptr(), 0.0 );
      // Above, I believe that these are the correct initial conditions for
      // s_bar and s_bar_dot given how the EpetraExt::DiagonalTransientModel
      // is currently implemented!

      RCP<const Rythmos::StateAndForwardSensitivityModelEvaluator<Scalar> >
        stateAndSensModel = stateAndSensStepper->getStateAndFwdSensModel();

      Thyra::ModelEvaluatorBase::InArgs<Scalar>
        state_and_sens_ic = stateAndSensStepper->getModel()->createInArgs();

      // Copy time, parameters etc.
      state_and_sens_ic.setArgs(state_ic);
      // Set initial condition for x_bar = [ x; s_bar ]
      state_and_sens_ic.set_x(
        stateAndSensModel->create_x_bar_vec(state_ic.get_x(),s_bar_init)
        );
      // Set initial condition for x_bar_dot = [ x_dot; s_bar_dot ]
      state_and_sens_ic.set_x_dot(
        stateAndSensModel->create_x_bar_vec(state_ic.get_x_dot(),s_bar_dot_init)
        );

//      *out << "\nstate_and_sens_ic:\n" << Teuchos::describe(state_and_sens_ic,Teuchos::VERB_DEFAULT);

      stateAndSensStepper->setInitialCondition(state_and_sens_ic);

      //
      // Use a StepperAsModelEvaluator to integrate the state+sens
      //

      RCP<Rythmos::StepperAsModelEvaluator<Scalar> >
        stateAndSensIntegratorAsModel = Rythmos::stepperAsModelEvaluator(
          Teuchos::rcp_implicit_cast<Rythmos::StepperBase<Scalar> >(stateAndSensStepper),
          Teuchos::rcp_implicit_cast<Rythmos::IntegratorBase<Scalar> >(fwdStateIntegrator),
          state_and_sens_ic
          );

      *out << "\nUse the StepperAsModelEvaluator to integrate state + sens x_bar(p,t_final) ... \n";

      RCP< Thyra::VectorBase<Scalar> > x_bar_final;

      Teuchos::OSTab tab(out);

      x_bar_final = createMember(stateAndSensIntegratorAsModel->get_g_space(0));

      eval_g(
        *stateAndSensIntegratorAsModel,
        0, *state_ic.get_p(0),
        t_final,
        0, &*x_bar_final
        );

      *out
        << "\nx_bar_final = x_bar(p,t_final) evaluated using "
        << "stateAndSensIntegratorAsModel:\n"
        << Teuchos::describe(*x_bar_final,solnVerbLevel);

     // As post-processing step, calc responses and gradient at final solution
     finalSolution = Thyra::productVectorBase<Scalar>(x_bar_final)->getVectorBlock(0);

      *out << "\nF) Check the solution to the forward problem ...\n";

     // Extract sensitivity vectors into Epetra_MultiVector
     Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > dxdp =
         Teuchos::rcp_dynamic_cast<const Thyra::DefaultMultiVectorProductVector<Scalar> >
           (Thyra::productVectorBase<Scalar>(x_bar_final)->getVectorBlock(1))
             ->getMultiVector();;

     Thyra::assign(dgdp_out.ptr(), 0.0);

/*********************  NEED TO CONVERT TO THYRA *******************

     Teuchos::RCP<Epetra_MultiVector> dgdx
          = Teuchos::rcp(new Epetra_MultiVector(finalSolution->Map(),
                                                   dgdp_out->GlobalLength()));
     Teuchos::Array<int> p_indexes =
       outArgs.get_DgDp(0,0).getDerivativeMultiVector().getParamIndexes();
 
     EpetraExt::ModelEvaluator::DerivativeMultiVector dmv_dgdp(dgdp_out,
                                                               DERIV_MV_BY_COL,
                                                               p_indexes);
 
     EpetraExt::ModelEvaluator::InArgs model_inargs = model->createInArgs();
     EpetraExt::ModelEvaluator::OutArgs model_outargs = model->createOutArgs();
     model_inargs.set_x(finalSolution);
     model_inargs.set_p(0, p_in);

     if (g_out != Teuchos::null) {
       g_out->PutScalar(0.0);
       model_outargs.set_g(0, g_out);
     }
     model_outargs.set_DgDp(0,0,dmv_dgdp);
     model_outargs.set_DgDx(0,dgdx);

     model->evalModel(model_inargs, model_outargs);

 
     // (3) Calculate dg/dp = dg/dx*dx/dp + dg/dp
     // This may be the transpose of what we want since we specified
     // we want dg/dp by column in createOutArgs().
     // In this case just interchange the order of dgdx and dxdp
     // We should really probably check what the underlying ME does

     if (Teuchos::VERB_MEDIUM <= solnVerbLevel) cout << " dgdx \n" << *dgdx << endl;
     if (Teuchos::VERB_MEDIUM <= solnVerbLevel) cout << " dxdp \n" << *dxdp << endl;

     dgdp_out->Multiply('T', 'N', 1.0, *dgdx, *dxdp, 1.0);
*********************/

   }

   // return the final solution as an additional g-vector, if requested
   if (gx_out != Teuchos::null)  Thyra::copy(*finalSolution, gx_out.ptr());
}

template <typename Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
Piro::RythmosSolver<Scalar>::getValidRythmosParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> validPL =
     Teuchos::rcp(new Teuchos::ParameterList("ValidRythmosParams"));

  validPL->set<std::string>("Nonlinear Solver Type", "");

  Teuchos::setStringToIntegralParameter<int>(
    "Nonlinear Solver Type",
    "Rythmos",
    "Determines which nonlinear solver to use.",
    Teuchos::tuple<std::string>("Rythmos","NOX"),
    validPL.get()
    );

  validPL->sublist("NonLinear Solver", false, "");
  validPL->sublist("Rythmos Builder", false, "");
  

  validPL->set<double>("Final Time", 1.0, "");
  validPL->sublist("Rythmos Stepper", false, "");
  validPL->sublist("Rythmos Integrator", false, "");
  validPL->set<std::string>("Rythmos Integration Control Strategy", "Simple", "");
  validPL->set<std::string>("Step Control Strategy Type", "None", "");
  validPL->sublist("Rythmos Step Control Strategy", false, "");
  validPL->sublist("Rythmos Integration Control", false, "");
  validPL->sublist("Stratimikos", false, "");
  validPL->set<std::string>("Verbosity Level", "", "");
  validPL->set<std::string>("Stepper Type", "", "");

  validPL->set<double>("Alpha", 1.0, "");
  validPL->set<double>("Beta", 1.0, "");
  validPL->set<double>("Max State Error", 1.0, "");
  validPL->set<std::string>("Name", "", "");
  validPL->set<bool>("Invert Mass Matrix", false, "");

  return validPL;
}
