//@HEADER

// ***********************************************************************
//
//                     Rythmos Package
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

#include "EpetraExt_DiagonalTransientModel.hpp"
#include "Rythmos_BackwardEulerStepper.hpp"
#include "Rythmos_ImplicitBDFStepper.hpp"
#include "Rythmos_ImplicitRKStepper.hpp"
#include "Rythmos_ForwardSensitivityStepper.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"
#include "Rythmos_DefaultIntegrator.hpp"
#include "Rythmos_SimpleIntegrationControlStrategy.hpp"
#include "Rythmos_StepperAsModelEvaluator.hpp"
#include "Rythmos_RKButcherTableauBuilder.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_DirectionalFiniteDiffCalculator.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_VerbosityLevelCommandLineProcessorHelpers.hpp"
#include "Teuchos_as.hpp"

#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif // HAVE_MPI


namespace {


const std::string TimeStepNonlinearSolver_name = "TimeStepNonlinearSolver";

const std::string Stratimikos_name = "Stratimikos";

const std::string DiagonalTransientModel_name = "DiagonalTransientModel";

const std::string RythmosStepper_name = "Rythmos Stepper";

const std::string RythmosIntegrator_name = "Rythmos Integrator";

const std::string FdCalc_name = "FD Calc";


Teuchos::RCP<const Teuchos::ParameterList>
getValidParameters()
{
  using Teuchos::RCP; using Teuchos::ParameterList;
  static RCP<const ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->sublist(TimeStepNonlinearSolver_name);
    pl->sublist(Stratimikos_name);
    pl->sublist(DiagonalTransientModel_name);
    pl->sublist(RythmosStepper_name);
    pl->sublist(RythmosIntegrator_name);
    pl->sublist(FdCalc_name);
    validPL = pl;
  }
  return validPL;
}


} // namespace


int main(int argc, char *argv[])
{

  using std::endl;
  typedef double Scalar;
  typedef double ScalarMag;
  using Teuchos::describe;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_implicit_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::as;
  using Teuchos::ParameterList;
  using Teuchos::CommandLineProcessor;
  typedef Teuchos::ParameterList::PrintOptions PLPrintOptions;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Thyra::DefaultMultiVectorProductVectorSpace<Scalar> DMVPVS;
  using Thyra::productVectorBase;
  
  bool result, success = true;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  RCP<Epetra_Comm> epetra_comm;
#ifdef HAVE_MPI
  epetra_comm = rcp( new Epetra_MpiComm(MPI_COMM_WORLD) );
#else
  epetra_comm = rcp( new Epetra_SerialComm );
#endif // HAVE_MPI

  RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {

    //
    // Read commandline options
    //

    CommandLineProcessor clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);

    std::string paramsFileName = "";
    clp.setOption( "params-file", &paramsFileName,
      "File name for XML parameters" );

    std::string extraParamsString = "";
    clp.setOption( "extra-params", &extraParamsString,
      "Extra XML parameters" );

    std::string extraParamsFile = "";
    clp.setOption( "extra-params-file", &extraParamsFile, "File containing extra parameters in XML format.");

    double maxStateError = 1e-6;
    clp.setOption( "max-state-error", &maxStateError,
      "The maximum allowed error in the integrated state in relation to the exact state solution" );

    double finalTime = 1e-3;
    clp.setOption( "final-time", &finalTime,
      "Final integration time (initial time is 0.0)" );

    int numTimeSteps = 10;
    clp.setOption( "num-time-steps", &numTimeSteps,
      "Number of (fixed) time steps.  If <= 0.0, then variable time steps are taken" );

    bool useBDF = false;
    clp.setOption( "use-BDF", "use-BE", &useBDF,
      "Use BDF or Backward Euler (BE)" );

    bool useIRK = false;
    clp.setOption( "use-IRK", "use-other", &useIRK,
      "Use IRK or something" );

    bool doFwdSensSolve = false;
    clp.setOption( "fwd-sens-solve", "state-solve", &doFwdSensSolve,
      "Do the forward sensitivity solve or just the state solve" );

    bool doFwdSensErrorControl = false;
    clp.setOption( "fwd-sens-err-cntrl", "no-fwd-sens-err-cntrl", &doFwdSensErrorControl,
      "Do error control on the forward sensitivity solve or not" );

    double maxRestateError = 0.0;
    clp.setOption( "max-restate-error", &maxRestateError,
      "The maximum allowed error between the state integrated by itself verses integrated along with DxDp" );

    double maxSensError = 1e-4;
    clp.setOption( "max-sens-error", &maxSensError,
      "The maximum allowed error in the integrated sensitivity in relation to"
      " the finite-difference sensitivity" );

    Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_DEFAULT;
    setVerbosityLevelOption( "verb-level", &verbLevel,
      "Top-level verbosity level.  By default, this gets deincremented as you go deeper into numerical objects.",
      &clp );

    bool testExactSensitivity = false;
    clp.setOption( "test-exact-sens", "no-test-exact-sens", &testExactSensitivity,
      "Test the exact sensitivity with finite differences or not." );

    bool dumpFinalSolutions = false;
    clp.setOption(
      "dump-final-solutions", "no-dump-final-solutions", &dumpFinalSolutions,
      "Determine if the final solutions are dumpped or not." );
    
    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;
    
    if ( Teuchos::VERB_DEFAULT == verbLevel )
      verbLevel = Teuchos::VERB_LOW;

    const Teuchos::EVerbosityLevel
      solnVerbLevel = ( dumpFinalSolutions ? Teuchos::VERB_EXTREME : verbLevel );

    //
    // Get the base parameter list that all other parameter lists will be read
    // from.
    //
    
    RCP<ParameterList>
      paramList = Teuchos::parameterList();
    if (paramsFileName.length())
      updateParametersFromXmlFile( paramsFileName, paramList.ptr() );
    if(extraParamsFile.length())
      Teuchos::updateParametersFromXmlFile( "./"+extraParamsFile, paramList.ptr() );
    if (extraParamsString.length())
      updateParametersFromXmlString( extraParamsString, paramList.ptr() );

    if (testExactSensitivity) {
      paramList->sublist(DiagonalTransientModel_name).set("Exact Solution as Response",true);
    }

    paramList->validateParameters(*getValidParameters(),0); // Only validate top level lists!

    //
    // Create the Stratimikos linear solver factory.
    //
    // This is the linear solve strategy that will be used to solve for the
    // linear system with the W.
    //

    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
    linearSolverBuilder.setParameterList(sublist(paramList,Stratimikos_name));
    RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >
      W_factory = createLinearSolveStrategy(linearSolverBuilder);
    
    //
    // Create the underlying EpetraExt::ModelEvaluator
    //

    RCP<EpetraExt::DiagonalTransientModel>
      epetraStateModel = EpetraExt::diagonalTransientModel(
        epetra_comm,
        sublist(paramList,DiagonalTransientModel_name)
        );

    *out <<"\nepetraStateModel valid options:\n";
    epetraStateModel->getValidParameters()->print(
      *out, PLPrintOptions().indent(2).showTypes(true).showDoc(true)
      );
    
    //
    // Create the Thyra-wrapped ModelEvaluator
    //
    
    RCP<Thyra::ModelEvaluator<double> >
      stateModel = epetraModelEvaluator(epetraStateModel,W_factory);

    *out << "\nParameter names = " << *stateModel->get_p_names(0) << "\n";
    
    //
    // Create the Rythmos stateStepper
    //

    RCP<Rythmos::TimeStepNonlinearSolver<double> >
      nonlinearSolver = Rythmos::timeStepNonlinearSolver<double>();
    RCP<ParameterList>
      nonlinearSolverPL = sublist(paramList,TimeStepNonlinearSolver_name);
    nonlinearSolverPL->get("Default Tol",1e-3*maxStateError); // Set default if not set
    nonlinearSolver->setParameterList(nonlinearSolverPL);

    RCP<Rythmos::StepperBase<Scalar> > stateStepper;

    if (useBDF) {
      stateStepper = rcp(
        new Rythmos::ImplicitBDFStepper<double>(
          stateModel, nonlinearSolver
          )
        );
    }
    else if (useIRK) {
      // We need a separate LOWSFB object for the IRK stepper
      RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >
        irk_W_factory = createLinearSolveStrategy(linearSolverBuilder);
      RCP<Rythmos::RKButcherTableauBase<double> > irkbt = Rythmos::createRKBT<double>("Backward Euler");
      stateStepper = Rythmos::implicitRKStepper<double>(
        stateModel, nonlinearSolver, irk_W_factory, irkbt
        );
    }
    else {
      stateStepper = rcp(
        new Rythmos::BackwardEulerStepper<double>(
          stateModel, nonlinearSolver
          )
        );
    }

    *out <<"\nstateStepper:\n" << describe(*stateStepper,verbLevel);
    *out <<"\nstateStepper valid options:\n";
    stateStepper->getValidParameters()->print(
      *out, PLPrintOptions().indent(2).showTypes(true).showDoc(true)
      );

    stateStepper->setParameterList(sublist(paramList,RythmosStepper_name));

    //
    // Setup finite difference objects that will be used for tests
    //

    Thyra::DirectionalFiniteDiffCalculator<Scalar> fdCalc;
    fdCalc.setParameterList(sublist(paramList,FdCalc_name));
    fdCalc.setOStream(out);
    fdCalc.setVerbLevel(verbLevel);

    //
    // Use a StepperAsModelEvaluator to integrate the state
    //

    const MEB::InArgs<Scalar>
      state_ic = stateModel->getNominalValues();
    *out << "\nstate_ic:\n" << describe(state_ic,verbLevel);

    RCP<Rythmos::IntegratorBase<Scalar> > integrator;
    {
      RCP<ParameterList>
        integratorPL = sublist(paramList,RythmosIntegrator_name);
      integratorPL->set( "Take Variable Steps", as<bool>(numTimeSteps < 0) );
      integratorPL->set( "Fixed dt", as<double>((finalTime - state_ic.get_t())/numTimeSteps) );
      RCP<Rythmos::IntegratorBase<Scalar> >
        defaultIntegrator = Rythmos::controlledDefaultIntegrator<Scalar>(
          Rythmos::simpleIntegrationControlStrategy<Scalar>(integratorPL)
          );
      integrator = defaultIntegrator;
    }

    RCP<Rythmos::StepperAsModelEvaluator<Scalar> >
      stateIntegratorAsModel = Rythmos::stepperAsModelEvaluator(
        stateStepper, integrator, state_ic
        );
    stateIntegratorAsModel->setVerbLevel(verbLevel);
    
    *out << "\nUse the StepperAsModelEvaluator to integrate state x(p,finalTime) ... \n";
    
    RCP<Thyra::VectorBase<Scalar> > x_final;

    {
      
      Teuchos::OSTab tab(out);

      x_final = createMember(stateIntegratorAsModel->get_g_space(0));
      
      eval_g(
        *stateIntegratorAsModel,
        0, *state_ic.get_p(0),
        finalTime,
        0, &*x_final
        );

      *out
        << "\nx_final = x(p,finalTime) evaluated using stateIntegratorAsModel:\n"
        << describe(*x_final,solnVerbLevel);

    }

    //
    // Test the integrated state against the exact analytical state solution
    //

    RCP<const Thyra::VectorBase<Scalar> >
      exact_x_final = create_Vector(
        epetraStateModel->getExactSolution(finalTime),
        stateModel->get_x_space()
        );
    
    result = Thyra::testRelNormDiffErr(
      "exact_x_final", *exact_x_final, "x_final", *x_final,
      "maxStateError", maxStateError, "warningTol", 1.0, // Don't warn
      &*out, solnVerbLevel
      );
    if (!result) success = false;

    //
    // Solve and test the forward sensitivity computation
    //
      
    if (doFwdSensSolve) {

      //
      // Create the forward sensitivity stepper
      //
      
      RCP<Rythmos::ForwardSensitivityStepper<Scalar> > stateAndSensStepper = 
        Rythmos::forwardSensitivityStepper<Scalar>();
      if (doFwdSensErrorControl) {
        stateAndSensStepper->initializeDecoupledSteppers(
          stateModel, 0, stateModel->getNominalValues(),
          stateStepper, nonlinearSolver,
          integrator->cloneIntegrator(), finalTime
          );
      }
      else {
        stateAndSensStepper->initializeSyncedSteppers(
          stateModel, 0, stateModel->getNominalValues(),
          stateStepper, nonlinearSolver
          );
        // The above call will result in stateStepper and nonlinearSolver being
        // cloned.  This helps to ensure consistency between the state and
        // sensitivity computations!
      }

      //
      // Set the initial condition for the state and forward sensitivities
      //

      RCP<Thyra::VectorBase<Scalar> > s_bar_init
        = createMember(stateAndSensStepper->getFwdSensModel()->get_x_space());
      assign( s_bar_init.ptr(), 0.0 );
      RCP<Thyra::VectorBase<Scalar> > s_bar_dot_init
        = createMember(stateAndSensStepper->getFwdSensModel()->get_x_space());
      assign( s_bar_dot_init.ptr(), 0.0 );
      // Above, I believe that these are the correct initial conditions for
      // s_bar and s_bar_dot given how the EpetraExt::DiagonalTransientModel
      // is currently implemented!

      RCP<const Rythmos::StateAndForwardSensitivityModelEvaluator<Scalar> >
        stateAndSensModel = stateAndSensStepper->getStateAndFwdSensModel();

      MEB::InArgs<Scalar>
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

      *out << "\nstate_and_sens_ic:\n" << describe(state_and_sens_ic,verbLevel);
 
      stateAndSensStepper->setInitialCondition(state_and_sens_ic);

      //
      // Use a StepperAsModelEvaluator to integrate the state+sens
      //
    
      RCP<Rythmos::StepperAsModelEvaluator<Scalar> >
        stateAndSensIntegratorAsModel = Rythmos::stepperAsModelEvaluator(
          rcp_implicit_cast<Rythmos::StepperBase<Scalar> >(stateAndSensStepper),
          integrator, state_and_sens_ic
          );
      stateAndSensIntegratorAsModel->setVerbLevel(verbLevel);
    
      *out << "\nUse the StepperAsModelEvaluator to integrate state + sens x_bar(p,finalTime) ... \n";
    
      RCP<Thyra::VectorBase<Scalar> > x_bar_final;

      {
      
        Teuchos::OSTab tab(out);

        x_bar_final = createMember(stateAndSensIntegratorAsModel->get_g_space(0));
      
        eval_g(
          *stateAndSensIntegratorAsModel,
          0, *state_ic.get_p(0),
          finalTime,
          0, &*x_bar_final
          );

        *out
          << "\nx_bar_final = x_bar(p,finalTime) evaluated using stateAndSensIntegratorAsModel:\n"
          << describe(*x_bar_final,solnVerbLevel);

      }

      //
      // Test that the state computed above is same as computed initially!
      //

      *out << "\nChecking that x(p,finalTime) computed as part of x_bar above is the same ...\n";

      {

        Teuchos::OSTab tab(out);

        RCP<const Thyra::VectorBase<Scalar> >
          x_in_x_bar_final = productVectorBase<Scalar>(x_bar_final)->getVectorBlock(0);

        result = Thyra::testRelNormDiffErr<Scalar>(
          "x_final", *x_final,
          "x_in_x_bar_final", *x_in_x_bar_final,
          "maxRestateError", maxRestateError,
          "warningTol", 1.0, // Don't warn
          &*out, solnVerbLevel
          );
        if (!result) success = false;

      }

      //
      // Compute DxDp using finite differences
      //

      *out << "\nApproximating DxDp(p,t) using directional finite differences of integrator for x(p,t) ...\n";
    
      RCP<Thyra::MultiVectorBase<Scalar> > DxDp_fd_final;

      {
      
        Teuchos::OSTab tab(out);
      
      
        MEB::InArgs<Scalar>
          fdBasePoint = stateIntegratorAsModel->createInArgs();
      
        fdBasePoint.set_t(finalTime);
        fdBasePoint.set_p(0,stateModel->getNominalValues().get_p(0));
      
        DxDp_fd_final = createMembers(
          stateIntegratorAsModel->get_g_space(0),
          stateIntegratorAsModel->get_p_space(0)->dim()
          );
      
        typedef Thyra::DirectionalFiniteDiffCalculatorTypes::SelectedDerivatives
          SelectedDerivatives;
      
        MEB::OutArgs<Scalar> fdOutArgs =
          fdCalc.createOutArgs(
            *stateIntegratorAsModel,
            SelectedDerivatives().supports(MEB::OUT_ARG_DgDp,0,0)
            );
        fdOutArgs.set_DgDp(0,0,DxDp_fd_final);
      
        // Silence the model evaluators that are called.  The fdCal object
        // will show all of the inputs and outputs for each call.
        stateStepper->setVerbLevel(Teuchos::VERB_NONE);
        stateIntegratorAsModel->setVerbLevel(Teuchos::VERB_NONE);
      
        fdCalc.calcDerivatives(
          *stateIntegratorAsModel, fdBasePoint,
          stateIntegratorAsModel->createOutArgs(), // Don't bother with function value
          fdOutArgs
          );
        
        *out
          << "\nFinite difference DxDp_fd_final = DxDp(p,finalTime): "
          << describe(*DxDp_fd_final,solnVerbLevel);

      }

      //
      // Test that the integrated sens and the F.D. sens are similar
      //

      *out << "\nChecking that integrated DxDp(p,finalTime) and finite-diff DxDp(p,finalTime) are similar ...\n";

      {

        Teuchos::OSTab tab(out);

        RCP<const Thyra::VectorBase<Scalar> >
          DxDp_vec_final = Thyra::productVectorBase<Scalar>(x_bar_final)->getVectorBlock(1);
      
        RCP<const Thyra::VectorBase<Scalar> >
          DxDp_fd_vec_final = Thyra::multiVectorProductVector(
            rcp_dynamic_cast<const Thyra::DefaultMultiVectorProductVectorSpace<Scalar> >(
              DxDp_vec_final->range()
              ),
            DxDp_fd_final
            );

        result = Thyra::testRelNormDiffErr(
          "DxDp_vec_final", *DxDp_vec_final,
          "DxDp_fd_vec_final", *DxDp_fd_vec_final,
          "maxSensError", maxSensError,
          "warningTol", 1.0, // Don't warn
          &*out, solnVerbLevel
          );
        if (!result) success = false;

      }

    }
    
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,success);

  if(success)
    *out << "\nEnd Result: TEST PASSED" << endl;
  else
    *out << "\nEnd Result: TEST FAILED" << endl;
  
  return ( success ? 0 : 1 );

} // end main() [Doxygen looks for this!]

