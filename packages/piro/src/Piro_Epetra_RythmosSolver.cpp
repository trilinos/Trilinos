// @HEADER
// ************************************************************************
// 
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// 
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#include <cmath>

#include "Piro_Epetra_RythmosSolver.hpp"
#include "Piro_Epetra_InvertMassMatrixDecorator.hpp"
#include "Piro_ValidPiroParameters.hpp"
#include "Piro_Epetra_MatrixFreeDecorator.hpp"

#include "EpetraExt_ModelEvaluator.h"

#include "Rythmos_BackwardEulerStepper.hpp"
#include "Rythmos_ExplicitRKStepper.hpp"
#include "Rythmos_SimpleIntegrationControlStrategy.hpp"
#include "Rythmos_ForwardSensitivityStepper.hpp"
#include "Rythmos_StepperAsModelEvaluator.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"


Piro::Epetra::RythmosSolver::RythmosSolver(Teuchos::RCP<Teuchos::ParameterList> appParams_,
                          Teuchos::RCP<EpetraExt::ModelEvaluator> model_,
                          Teuchos::RCP<Rythmos::IntegrationObserverBase<Scalar> > observer ) :
  appParams(appParams_),
  model(model_)
{
  //appParams->validateParameters(*Piro::getValidPiroParameters(),0);

  using Teuchos::RCP;
  using Teuchos::rcp;

      out = Teuchos::VerboseObjectBase::getDefaultOStream();

  // Allow for Matrix-Free implementation
  string jacobianSource = appParams->get("Jacobian Operator", "Have Jacobian");
  if (jacobianSource == "Matrix-Free") {
    if (appParams->isParameter("Matrix-Free Perturbation")) {
      model = Teuchos::rcp(new Piro::Epetra::MatrixFreeDecorator(model,
                           appParams->get<double>("Matrix-Free Perturbation")));
    }
    else model = Teuchos::rcp(new Piro::Epetra::MatrixFreeDecorator(model));
  }


  num_p = model->createInArgs().Np();
  num_g = model->createOutArgs().Ng();

  TEST_FOR_EXCEPTION(num_p > 1, Teuchos::Exceptions::InvalidParameter,
                     std::endl << "Error in Piro::Epetra::RythmosSolver " <<
                     "Not Implemented for Np>1 : " << num_p << std::endl);
  TEST_FOR_EXCEPTION(num_g > 1, Teuchos::Exceptions::InvalidParameter,
                     std::endl << "Error in Piro::Epetra::RythmosSolver " <<
                     "Not Implemented for Ng>1 : " << num_g << std::endl);

      //
      *out << "\nA) Get the base parameter list ...\n";
      //

      RCP<Teuchos::ParameterList> rythmosPL = sublist(appParams, "Rythmos", true);
      rythmosPL->validateParameters(*getValidRythmosParameters(),0);

      {
        const string verbosity = rythmosPL->get("Verbosity Level", "VERB_DEFAULT");
        solnVerbLevel = Teuchos::VERB_DEFAULT;
        if      (verbosity == "VERB_NONE")    solnVerbLevel = Teuchos::VERB_NONE;
        else if (verbosity == "VERB_LOW")     solnVerbLevel = Teuchos::VERB_LOW;
        else if (verbosity == "VERB_MEDIUM")  solnVerbLevel = Teuchos::VERB_MEDIUM;
        else if (verbosity == "VERB_HIGH")    solnVerbLevel = Teuchos::VERB_HIGH;
        else if (verbosity == "VERB_EXTREME") solnVerbLevel = Teuchos::VERB_EXTREME;
      }

      const int numTimeSteps = rythmosPL->get("Num Time Steps", 10);
      const Scalar t_init = 0.0;
      t_final = rythmosPL->get("Final Time", 0.1);
      
      const Rythmos::TimeRange<Scalar> fwdTimeRange(t_init, t_final);
      const Scalar delta_t = t_final / (double) numTimeSteps;
      *out << "\ndelta_t = " << delta_t;

      const string stepperType = rythmosPL->get("Stepper Type", "Backward Euler");

      //
      *out << "\nB) Create the Stratimikos linear solver factory ...\n";
      //
      // This is the linear solve strategy that will be used to solve for the
      // linear system with the W.
      //
      
      Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
      linearSolverBuilder.setParameterList(sublist(rythmosPL, "Stratimikos", true));
      RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >
	W_factory = createLinearSolveStrategy(linearSolverBuilder);
      
      //
      *out << "\nC) Create and initalize the forward model ...\n";
      //
      
      // C.1) Create the underlying EpetraExt::ModelEvaluator
      
      // already constructed as "model". Decorate if needed.
      
      if (stepperType == "Explicit RK") {
        if (rythmosPL->get("Invert Mass Matrix", false)) {
          Teuchos::RCP<EpetraExt::ModelEvaluator> origModel = model;
          bool lump=rythmosPL->get("Lump Mass Matrix", false);
          model = Teuchos::rcp(new Piro::Epetra::InvertMassMatrixDecorator(
                   sublist(rythmosPL,"Stratimikos", true), origModel));
        }
      }
      
      // C.2) Create the Thyra-wrapped ModelEvaluator
      
      fwdStateModel = epetraModelEvaluator(model, W_factory);
      
      const RCP<const Thyra::VectorSpaceBase<Scalar> >
	x_space = fwdStateModel->get_x_space();
      
      //
      *out << "\nD) Create the stepper and integrator for the forward problem ...\n";
      //

      fwdTimeStepSolver = Rythmos::timeStepNonlinearSolver<double>();

      if (rythmosPL->getEntryPtr("NonLinear Solver")) {
        RCP<Teuchos::ParameterList> nonlinePL =
           sublist(rythmosPL, "NonLinear Solver", true);
        fwdTimeStepSolver->setParameterList(nonlinePL);
      }

      if (stepperType == "Backward Euler")
        fwdStateStepper = Rythmos::backwardEulerStepper<double>
           (fwdStateModel, fwdTimeStepSolver);
      else if (stepperType == "Explicit RK")
        fwdStateStepper = Rythmos::explicitRKStepper<double>(fwdStateModel);
      else 
        TEST_FOR_EXCEPTION( true, Teuchos::Exceptions::InvalidParameter,
                     std::endl << "Error! Piro::Epetra::RythmosSolver: Invalid Steper Type: "
                     << stepperType << std::endl);

      fwdStateStepper->setParameterList(sublist(rythmosPL, "Rythmos Stepper", true));

      {
        RCP<Teuchos::ParameterList>
          integrationControlPL = sublist(rythmosPL, "Rythmos Integration Control", true);
        // set default to false
        bool var = integrationControlPL->get( "Take Variable Steps", false );
        integrationControlPL->set( "Fixed dt", Teuchos::as<double>(delta_t) );

       // RCP<Rythmos::IntegratorBase<Scalar> >
        RCP<Rythmos::DefaultIntegrator<Scalar> >
          defaultIntegrator = Rythmos::controlledDefaultIntegrator<Scalar>(
            Rythmos::simpleIntegrationControlStrategy<Scalar>(integrationControlPL)
            );
        fwdStateIntegrator = defaultIntegrator;
      }
      fwdStateIntegrator->setParameterList(sublist(rythmosPL, "Rythmos Integrator", true));

   if (observer != Teuchos::null) 
     fwdStateIntegrator->setIntegrationObserver(observer);

}

Piro::Epetra::RythmosSolver::~RythmosSolver()
{
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::RythmosSolver::get_x_map() const
{
  Teuchos::RCP<const Epetra_Map> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::RythmosSolver::get_f_map() const
{
  Teuchos::RCP<const Epetra_Map> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::RythmosSolver::get_p_map(int l) const
{
  TEST_FOR_EXCEPTION(l >= num_p || l < 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error in Piro::Epetra::RythmosSolver::get_p_map():  " <<
                     "Invalid parameter index l = " <<
                     l << std::endl);
  return model->get_p_map(l);
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::RythmosSolver::get_g_map(int j) const
{
  TEST_FOR_EXCEPTION(j > num_g || j < 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error in Piro::Epetra::RythmosSolver::get_g_map():  " <<
                     "Invalid response index j = " <<
                     j << std::endl);

  if      (j < num_g) return model->get_g_map(j);
  else if (j == num_g) return model->get_x_map();
}

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::RythmosSolver::get_x_init() const
{
  Teuchos::RCP<const Epetra_Vector> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::RythmosSolver::get_p_init(int l) const
{
  TEST_FOR_EXCEPTION(l >= num_p || l < 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error in Piro::Epetra::RythmosSolver::get_p_init():  " <<
                     "Invalid parameter index l = " <<
                     l << std::endl);
  return model->get_p_init(l);
}

EpetraExt::ModelEvaluator::InArgs Piro::Epetra::RythmosSolver::createInArgs() const
{
  //return underlyingME->createInArgs();
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(num_p);
//  inArgs.setSupports(IN_ARG_x,true);
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs Piro::Epetra::RythmosSolver::createOutArgs() const
{
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());

  // Ng is 1 bigger then model-Ng so that the solution vector can be an outarg
  outArgs.set_Np_Ng(num_p, num_g+1);

  EpetraExt::ModelEvaluator::OutArgs model_outargs = model->createOutArgs();
  for (int i=0; i<num_g; i++)
    for (int j=0; j<num_p; j++)
      if (!model_outargs.supports(OUT_ARG_DgDp, i, j).none())
        outArgs.setSupports(OUT_ARG_DgDp, i, j,
                            DerivativeSupport(DERIV_MV_BY_COL));

  return outArgs;
}

void Piro::Epetra::RythmosSolver::evalModel( const InArgs& inArgs,
                                     const OutArgs& outArgs ) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Parse InArgs

  RCP<const Epetra_Vector> p_in;
  if (num_p > 0) p_in = inArgs.get_p(0);

  // Parse OutArgs: always 1 extra
  RCP<Epetra_Vector> g_out; 
  if (num_g > 0) g_out = outArgs.get_g(0); 
  RCP<Epetra_Vector> gx_out = outArgs.get_g(num_g); 

  // Parse out-args for sensitivity calculation
  RCP<Epetra_MultiVector> dgdp_out;
  if (num_p>0 && num_g>0)
    dgdp_out = outArgs.get_DgDp(0,0).getMultiVector();

  RCP<const Epetra_Vector> finalSolution;

  Thyra::ModelEvaluatorBase::InArgs<Scalar>
     state_ic = fwdStateModel->getNominalValues();

  // Set paramters p_in as part of initial conditions
  if (num_p > 0)
    state_ic.set_p(0,Thyra::create_Vector(p_in, fwdStateModel->get_p_space(0)));

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
      const RCP<const Thyra::VectorBase<Scalar> > x_final = x_final_array[0];
  
      finalSolution = Thyra::get_Epetra_Vector(*model->get_x_map(), x_final);

      if (Teuchos::VERB_MEDIUM <= solnVerbLevel)
         cout << "Final Solution\n" << *finalSolution << std::endl;

     // As post-processing step, calc responses at final solution
     EpetraExt::ModelEvaluator::InArgs model_inargs = model->createInArgs();
     EpetraExt::ModelEvaluator::OutArgs model_outargs = model->createOutArgs();
     model_inargs.set_x(finalSolution);
     if (num_p > 0)  model_inargs.set_p(0, p_in);
     if (g_out != Teuchos::null) {
       g_out->PutScalar(0.0);
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
          fwdStateModel, 0, fwdStateModel->getNominalValues(),
          fwdStateStepper, fwdTimeStepSolver);

      //
      // Set the initial condition for the state and forward sensitivities
      //

      RCP<Thyra::VectorBase<Scalar> > s_bar_init
        = createMember(stateAndSensStepper->getFwdSensModel()->get_x_space());
      assign( &*s_bar_init, 0.0 );
      RCP<Thyra::VectorBase<Scalar> > s_bar_dot_init
        = createMember(stateAndSensStepper->getFwdSensModel()->get_x_space());
      assign( &*s_bar_dot_init, 0.0 );
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

      RCP<Thyra::VectorBase<Scalar> > x_bar_final;

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
     const RCP<const Thyra::VectorBase<Scalar> > x_final =
           Thyra::productVectorBase<Scalar>(x_bar_final)->getVectorBlock(0);
     finalSolution = Thyra::get_Epetra_Vector(*model->get_x_map(), x_final);

      *out << "\nF) Check the solution to the forward problem ...\n";

     // Extract sensitivity vectors into Epetra_MultiVector
     Teuchos::RCP<const Epetra_MultiVector> dxdp =
       Thyra::get_Epetra_MultiVector(*model->get_x_map(),
         Teuchos::rcp_dynamic_cast<const Thyra::DefaultMultiVectorProductVector<Scalar> >
           (Thyra::productVectorBase<Scalar>(x_bar_final)->getVectorBlock(1))
             ->getMultiVector() );;

     dgdp_out->PutScalar(0.0);

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

   }

   // return the final solution as an additional g-vector, if requested
   if (gx_out != Teuchos::null)  *gx_out = *finalSolution;
}

Teuchos::RCP<const Teuchos::ParameterList>
Piro::Epetra::RythmosSolver::getValidRythmosParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> validPL =
     Teuchos::rcp(new Teuchos::ParameterList("ValidRythmosParams"));;
  validPL->sublist("NonLinear Solver", false, "");
  validPL->set<int>("Num Time Steps", 0, "");
  validPL->set<double>("Final Time", 1.0, "");
  validPL->sublist("Rythmos Stepper", false, "");
  validPL->sublist("Rythmos Integrator", false, "");
  validPL->sublist("Rythmos Integration Control", false, "");
  validPL->sublist("Stratimikos", false, "");
  validPL->set<std::string>("Verbosity Level", "", "");
  validPL->set<std::string>("Stepper Type", "", "");

  validPL->set<double>("Alpha", 1.0, "");
  validPL->set<double>("Beta", 1.0, "");
  validPL->set<double>("Max State Error", 1.0, "");
  validPL->set<std::string>("Name", "", "");
  validPL->set<bool>("Invert Mass Matrix", false, "");
  validPL->set<bool>("Lump Mass Matrix", false, "");
  validPL->set<std::string>("Stepper Method", "", "");
  return validPL;
}

