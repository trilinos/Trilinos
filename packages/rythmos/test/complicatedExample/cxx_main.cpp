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

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Version.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif // HAVE_MPI

#include "ExampleApplication.hpp"

// Includes for Rythmos:
#include "Rythmos_ConfigDefs.h"
#include "Rythmos_Version.h"
//#include "ExampleApplicationRythmosInterface.hpp"
#include "Rythmos_ForwardEulerStepper.hpp"
#include "Rythmos_BackwardEulerStepper.hpp"
#include "Rythmos_ExplicitRKStepper.hpp"
#include "Rythmos_ImplicitBDFStepper.hpp"
// DefaultIntegrator includes: 
#include "Rythmos_InterpolationBuffer.hpp"
#include "Rythmos_LinearInterpolator.hpp"
#include "Rythmos_HermiteInterpolator.hpp"
#include "Rythmos_DefaultIntegrator.hpp"
#include "Rythmos_SimpleIntegrationControlStrategy.hpp"
// Includes for Thyra:
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_LinearNonlinearSolver.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"
#include "Thyra_TestingTools.hpp"

#ifdef HAVE_RYTHMOS_NOX
#  include "Thyra_NonlinearSolver_NOX.hpp"
#endif

// Includes for Stratimikos:
#ifdef HAVE_RYTHMOS_STRATIMIKOS
#  include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#endif

#include <string>

// Includes for Teuchos:
#include "Teuchos_as.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_VerbosityLevelCommandLineProcessorHelpers.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

enum EMethod { METHOD_FE, METHOD_BE, METHOD_ERK, METHOD_BDF };
enum EStepMethod { STEP_TYPE_FIXED, STEP_TYPE_VARIABLE };

int main(int argc, char *argv[])
{

  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Array;

  bool verbose = false; // verbosity level.
  bool result, success = true; // determine if the run was successfull

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

#ifdef HAVE_MPI
  MPI_Comm mpiComm = MPI_COMM_WORLD;
#endif // HAVE_MPI

  try { // catch exceptions

    bool useNOX = false; // Use NOX?
    double lambda_min = -0.9;   // min ODE coefficient
    double lambda_max = -0.01;  // max ODE coefficient
    double coeff_s = 0.0;  // Default is no forcing term
    std::string lambda_fit = "linear"; // Lambda model
    int numElements = 1; // number of elements in vector
    double x0 = 10.0; // ODE initial condition
    double finalTime = 1.0; // ODE final time
    int N = 10;  // number of steps to take
    const int num_methods = 4;
    const EMethod method_values[] = { METHOD_FE, METHOD_BE, METHOD_ERK, METHOD_BDF };
    const char * method_names[] = { "FE", "BE", "ERK", "BDF" };
    EMethod method_val = METHOD_ERK;
    const int num_step_methods = 2;
    const EStepMethod step_method_values[] = { STEP_TYPE_FIXED, STEP_TYPE_VARIABLE };
    const char * step_method_names[] = { "fixed", "variable" };
    EStepMethod step_method_val = STEP_TYPE_FIXED;
    double maxError = 1e-6;
    double maxRestepError = 1.0e4*Teuchos::ScalarTraits<double>::prec();
    bool version = false;  // display version information 
    double reltol = 1.0e-2;
    double abstol = 1.0e-4;
    int maxOrder = 5;
    bool useIntegrator = false;
    int buffersize = 100;
    Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_LOW;

    // Parse the command-line options:
    Teuchos::CommandLineProcessor  clp(false); // Don't throw exceptions
    clp.addOutputSetupOptions(true);
#ifdef HAVE_RYTHMOS_STRATIMIKOS
    Stratimikos::DefaultLinearSolverBuilder lowsfCreator;
    lowsfCreator.setupCLP(&clp);
#endif
    clp.setOption( "x0", &x0, "Constant ODE initial condition." );
    clp.setOption( "T", &finalTime, "Final time for simulation." );
    clp.setOption( "lambda_min", &lambda_min, "Lower bound for ODE coefficient");
    clp.setOption( "lambda_max", &lambda_max, "Upper bound for ODE coefficient");
    clp.setOption( "lambda_fit", &lambda_fit, "Lambda model:  random, linear");
    clp.setOption( "force_coeff", &coeff_s, "Forcing term coefficient");
    clp.setOption( "numelements", &numElements, "Problem size");
    clp.setOption( "method", &method_val, num_methods, method_values, method_names, "Integration method" );
    clp.setOption( "stepmethod", &step_method_val, num_step_methods, step_method_values, step_method_names, "Stepping method" );
    clp.setOption( "numsteps", &N, "Number of integration steps to take" );
    clp.setOption( "maxerror", &maxError, "Maximum error" );
    clp.setOption( "max-restep-error", &maxRestepError, "Maximum error between solutions of two steppers" );
    clp.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not" );
    clp.setOption( "version", "run", &version, "Version of this code" );
    clp.setOption( "reltol", &reltol, "Relative Error Tolerance" );
    clp.setOption( "abstol", &abstol, "Absolute Error Tolerance" );
    clp.setOption( "maxorder", &maxOrder, "Maximum Implicit BDF order" );
    clp.setOption( "useintegrator", "normal", &useIntegrator, "Use DefaultIntegrator as integrator" );
    clp.setOption( "buffersize", &buffersize, "Number of solutions to store in InterpolationBuffer" );
    clp.setOption( "useNOX", "noNOX", &useNOX, "Use NOX as the nonlinear solver" );
    setVerbosityLevelOption( "verb-level", &verbLevel, "Overall verbosity level.", &clp );

    Teuchos::CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

#ifdef HAVE_RYTHMOS_STRATIMIKOS
    lowsfCreator.readParameters(out.get());
    *out << "\nThe parameter list after being read in:\n";
    lowsfCreator.getParameterList()->print(*out,2,true,false);
#endif

    if (version) // Display version information and exit.
    {
      *out << Rythmos::version() << std::endl;
      *out << "complicatedExample Version 0.1 - 06/23/05" << std::endl;
      return(0);
    }

    if (lambda_min > lambda_max)
    {
      *out << "lamba_min must be less than lambda_max" << std::endl;
      return(1);
    }

    if (finalTime <= 0.0)
    {
      *out << "Final simulation time must be > 0.0." << std::endl;
      return(1);
    }

    *out << std::setprecision(15);
    
    // Set up the parameter list for the application:
    Teuchos::ParameterList params;
    bool implicitFlag = ((method_val==METHOD_BE) | (method_val==METHOD_BDF));
    //*out << "implicitFlag = " << implicitFlag << std::endl;
    params.set( "Implicit", implicitFlag );
    params.set( "Lambda_min", lambda_min );
    params.set( "Lambda_max", lambda_max );
    params.set( "Lambda_fit", lambda_fit );
    params.set( "NumElements", numElements );
    params.set( "x0", x0 );
    params.set( "Coeff_s", coeff_s );
#ifdef HAVE_MPI
    RCP<Epetra_Comm> epetra_comm_ptr_ = rcp( new Epetra_MpiComm(mpiComm) );
#else
    RCP<Epetra_Comm> epetra_comm_ptr_ = rcp( new Epetra_SerialComm  );
#endif // HAVE_MPI

    // Create the factory for the LinearOpWithSolveBase object
    RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
      W_factory;
    if((method_val == METHOD_BE) | (method_val == METHOD_BDF)) {
#ifdef HAVE_RYTHMOS_STRATIMIKOS
      W_factory = lowsfCreator.createLinearSolveStrategy("");
      *out
        << "\nCreated a LinearOpWithSolveFactory described as:\n"
        << Teuchos::describe(*W_factory,Teuchos::VERB_MEDIUM);
#endif
    }

    // create interface to problem
    RCP<ExampleApplication>
      epetraModel = rcp(new ExampleApplication(epetra_comm_ptr_, params));
    RCP<Thyra::ModelEvaluator<double> >
      model = rcp(new Thyra::EpetraModelEvaluator(epetraModel,W_factory));
    Thyra::ModelEvaluatorBase::InArgs<double> model_ic = model->getNominalValues();
    RCP<ExampleApplication>
      epetraModelSlave = rcp(new ExampleApplication(epetra_comm_ptr_, params));
    RCP<Thyra::ModelEvaluator<double> >
      modelSlave = rcp(new Thyra::EpetraModelEvaluator(epetraModelSlave,W_factory));
    Thyra::ModelEvaluatorBase::InArgs<double> modelSlave_ic = modelSlave->getNominalValues();

    // Create Stepper object depending on command-line input
    std::string method;
    RCP<Rythmos::StepperBase<double> > stepper_ptr;
    RCP<Rythmos::StepperBase<double> > stepperSlave_ptr;
    if ( method_val == METHOD_ERK ) {
      stepper_ptr = Rythmos::explicitRKStepper<double>(model);
      stepper_ptr->setInitialCondition(model_ic);
      RCP<Teuchos::ParameterList> ERKparams = Teuchos::parameterList();
      ERKparams->sublist("VerboseObject").set(
        "Verbosity Level",
        Teuchos::getVerbosityLevelParameterValueName(verbLevel)
        );
      stepper_ptr->setParameterList(ERKparams);
      method = "Explicit Runge-Kutta of order 4";
      step_method_val = STEP_TYPE_FIXED;
    }
    else if (method_val == METHOD_FE) {
      stepper_ptr = Rythmos::forwardEulerStepper<double>(model);
      stepper_ptr->setInitialCondition(model_ic);
      RCP<Teuchos::ParameterList> FEparams = Teuchos::parameterList();
      FEparams->sublist("VerboseObject").set(
        "Verbosity Level",
        Teuchos::getVerbosityLevelParameterValueName(verbLevel)
        );
      stepper_ptr->setParameterList(FEparams);
      method = "Forward Euler";
      step_method_val = STEP_TYPE_FIXED;
    }
    else if ((method_val == METHOD_BE) | (method_val == METHOD_BDF)) {
      RCP<Thyra::NonlinearSolverBase<double> >
        nonlinearSolver;
      RCP<Thyra::NonlinearSolverBase<double> >
        nonlinearSolverSlave;
      if (useNOX) {
#ifdef HAVE_RYTHMOS_NOX
        RCP<Thyra::NOXNonlinearSolver>
          _nonlinearSolver = rcp(new Thyra::NOXNonlinearSolver);
        RCP<Thyra::NOXNonlinearSolver>
          _nonlinearSolverSlave = rcp(new Thyra::NOXNonlinearSolver);
        nonlinearSolver = _nonlinearSolver;
        nonlinearSolverSlave = _nonlinearSolverSlave;
#else // HAVE_RYTHMOS_NOX
        TEST_FOR_EXCEPT_MSG(true, "Error: NOX is not enabled!");
#endif // HAVE_RYTHMOS_NOX
      } 
      else {
        RCP<Rythmos::TimeStepNonlinearSolver<double> >
          _nonlinearSolver = rcp(new Rythmos::TimeStepNonlinearSolver<double>());
        RCP<Rythmos::TimeStepNonlinearSolver<double> >
          _nonlinearSolverSlave = rcp(new Rythmos::TimeStepNonlinearSolver<double>());
        RCP<Teuchos::ParameterList>
          nonlinearSolverPL = Teuchos::parameterList();
        nonlinearSolverPL->set("Default Tol",double(1e-3*maxError));
        _nonlinearSolver->setParameterList(nonlinearSolverPL);
        _nonlinearSolverSlave->setParameterList(nonlinearSolverPL);
        nonlinearSolver = _nonlinearSolver;
        nonlinearSolverSlave = _nonlinearSolverSlave;
      } 
      if (method_val == METHOD_BE) {
        stepper_ptr = rcp(
          new Rythmos::BackwardEulerStepper<double>(model,nonlinearSolver));
        RCP<Teuchos::ParameterList>
          BEparams = Teuchos::parameterList();
        BEparams->sublist("VerboseObject").set(
          "Verbosity Level",
          Teuchos::getVerbosityLevelParameterValueName(verbLevel)
          );
        stepper_ptr->setParameterList(BEparams);
        stepper_ptr->setInitialCondition(model_ic);
        method = "Backward Euler";
        step_method_val = STEP_TYPE_FIXED;
      } 
      else {
        RCP<Teuchos::ParameterList>
          BDFparams = Teuchos::parameterList();
        RCP<Teuchos::ParameterList> BDFStepControlPL =
          Teuchos::sublist(BDFparams,"Step Control Settings");

        BDFparams->sublist("VerboseObject").set(
          "Verbosity Level",
          Teuchos::getVerbosityLevelParameterValueName(verbLevel)
          );

        BDFStepControlPL->sublist("VerboseObject").set(
          "Verbosity Level",
          Teuchos::getVerbosityLevelParameterValueName(verbLevel)
          );
        BDFStepControlPL->set( "stopTime", finalTime );
        BDFStepControlPL->set( "maxOrder", maxOrder );
        BDFStepControlPL->set( "relErrTol", reltol );
        BDFStepControlPL->set( "absErrTol", abstol );
        stepper_ptr = rcp(
          new Rythmos::ImplicitBDFStepper<double>(model,nonlinearSolver,BDFparams)
          );
        stepper_ptr->setInitialCondition(model_ic);
        stepperSlave_ptr = rcp(
          new Rythmos::ImplicitBDFStepper<double>(modelSlave,nonlinearSolverSlave,BDFparams)
          );
        stepperSlave_ptr->setInitialCondition(modelSlave_ic);
        method = "Implicit BDF";
        // step_method_val setting is left alone in this case
      }
    }
    else {
      TEST_FOR_EXCEPT(true);
    }
    Rythmos::StepperBase<double> &stepper = *stepper_ptr;
    if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
      stepper.describe(*out,verbLevel);
    }

    int numSteps = 0;
    double t0 = 0.0;
    double dt = (finalTime-t0)/N;
    double time = t0;

    RCP<const Thyra::VectorBase<double> > x_computed_thyra_ptr;
    if (step_method_val == STEP_TYPE_FIXED)
    {
      if (useIntegrator)
      {
        // Create integrator using stepper and linear interpolation buffer:
        RCP<Rythmos::InterpolatorBase<double> > 
          linearInterpolator = Rythmos::linearInterpolator<double>();
        RCP<Rythmos::InterpolationBuffer<double> > 
          IB = Rythmos::interpolationBuffer<double>(linearInterpolator,buffersize);
        RCP<Teuchos::ParameterList> 
          IBParams = Teuchos::parameterList();
        IBParams->sublist("VerboseObject").set(
          "Verbosity Level",
          Teuchos::getVerbosityLevelParameterValueName(verbLevel)
          );
        IB->setParameterList(IBParams);
        // Set up fixed-step-size integration:
        RCP<Teuchos::ParameterList> 
          integratorParams = Teuchos::parameterList();
        integratorParams->set( "Take Variable Steps", false );
        integratorParams->set( "Fixed dt", dt );
        integratorParams->sublist("VerboseObject").set(
          "Verbosity Level",
          Teuchos::getVerbosityLevelParameterValueName(verbLevel)
          );
        RCP<Rythmos::IntegratorBase<double> >
          defaultIntegrator = Rythmos::controlledDefaultIntegrator<double>(
            Rythmos::simpleIntegrationControlStrategy<double>(integratorParams)
            );
        RCP<Rythmos::DefaultIntegrator<double> > dIntegrator = 
          Teuchos::rcp_dynamic_cast<Rythmos::DefaultIntegrator<double> >(defaultIntegrator,true);
        Rythmos::DefaultIntegrator<double> &integrator = *dIntegrator;
        integrator.setStepper(stepper_ptr,finalTime);
        integrator.setTrailingInterpolationBuffer(IB);
        // Ask for desired time value:
        Array<double> time_vals;
        for (int i=0 ; i<=N ; ++i)
          time_vals.push_back(i*dt);
        Array<RCP<const Thyra::VectorBase<double> > > x_vec;
        Array<RCP<const Thyra::VectorBase<double> > > xdot_vec;
        Array<double> accuracy_vec;
        integrator.getPoints(time_vals,&x_vec,&xdot_vec,&accuracy_vec);
        // Get solution out of stepper:
        x_computed_thyra_ptr = x_vec.back();
      }
      else
      {
        // Integrate forward with fixed step sizes:
        for (int i=1 ; i<=N ; ++i)
        {
          double dt_taken = stepper.takeStep(dt,Rythmos::STEP_TYPE_FIXED);
          time += dt_taken;
          numSteps++;
          if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) )
            stepper.describe(*out,verbLevel);
          if (dt_taken != dt)
          {
            cerr << "Error, stepper took step of dt = " << dt_taken 
              << " when asked to take step of dt = " << dt << std::endl;
            break;
          }
        }
        // Get solution out of stepper:
        Rythmos::StepStatus<double> stepStatus = stepper.getStepStatus();
        x_computed_thyra_ptr = stepStatus.solution;
      }
    }
    else // (step_method_val == STEP_TYPE_VARIABLE)
    {
      if (useIntegrator)
      {
        // Set up fixed-step-size integration:
        // Create integrator using stepper and interpolation buffer:
        //RCP<Rythmos::InterpolatorBase<double> > 
        //  linearInterpolator = Rythmos::linearInterpolator<double>();
        //RCP<Rythmos::InterpolationBuffer<double> > 
        //  IB = Rythmos::interpolationBuffer<double>(linearInterpolator,buffersize);
        RCP<Rythmos::InterpolatorBase<double> > 
          hermiteInterpolator = Rythmos::hermiteInterpolator<double>();
        RCP<Rythmos::InterpolationBuffer<double> > 
          IB = Rythmos::interpolationBuffer<double>(hermiteInterpolator,buffersize);
        RCP<Teuchos::ParameterList> 
          IBParams = Teuchos::parameterList();
        IBParams->sublist("VerboseObject").set(
          "Verbosity Level",
          Teuchos::getVerbosityLevelParameterValueName(verbLevel)
          );
        IB->setParameterList(IBParams);
        RCP<Teuchos::ParameterList> 
          integratorParams = Teuchos::parameterList();
        integratorParams->set( "Take Variable Steps", true );
        integratorParams->sublist("VerboseObject").set(
          "Verbosity Level",
          Teuchos::getVerbosityLevelParameterValueName(verbLevel)
          );
        RCP<Rythmos::IntegratorBase<double> >
          defaultIntegrator = Rythmos::controlledDefaultIntegrator<double>(
            Rythmos::simpleIntegrationControlStrategy<double>(integratorParams)
            );
        RCP<Rythmos::DefaultIntegrator<double> > dIntegrator = 
          Teuchos::rcp_dynamic_cast<Rythmos::DefaultIntegrator<double> >(defaultIntegrator,true);
        Rythmos::DefaultIntegrator<double> &integrator = *dIntegrator;
        integrator.setStepper(stepper_ptr,finalTime);
        integrator.setTrailingInterpolationBuffer(IB);
        // Ask for desired time value:
        Array<double> time_vals;
        for (int i=0 ; i<=N ; ++i)
          time_vals.push_back(i*dt);
        Array<RCP<const Thyra::VectorBase<double> > > x_vec;
        Array<RCP<const Thyra::VectorBase<double> > > xdot_vec;
        Array<double> accuracy_vec;
        integrator.getPoints(time_vals,&x_vec,&xdot_vec,&accuracy_vec);
        // Get solution out of stepper:
        x_computed_thyra_ptr = x_vec.back();
      }
      else
      {
        // HIGH output data structures:
        // Create a place to store the computed solutions
        Rythmos::StepStatus<double> stepStatus = stepper.getStepStatus();
        x_computed_thyra_ptr = stepStatus.solution;
        // Convert Thyra::VectorBase to Epetra_Vector
        // Create a place to store the exact numerical solution
        RCP<Epetra_Vector> x_numerical_exact_ptr = rcp(new Epetra_Vector(*epetraModel->get_x_map()));
        Epetra_Vector& x_numerical_exact = *x_numerical_exact_ptr;
        // Create a place to store the relative difference:
        RCP<Epetra_Vector> x_rel_diff_ptr = rcp(new Epetra_Vector(*epetraModel->get_x_map()));
        Epetra_Vector& x_rel_diff = *x_rel_diff_ptr;
        // get lambda from the problem:
        RCP<const Epetra_Vector> lambda_ptr = epetraModel->get_coeff();
        const Epetra_Vector &lambda = *lambda_ptr;

        while (time < finalTime)
        {
          double dt_taken = stepper.takeStep(0.0,Rythmos::STEP_TYPE_VARIABLE);
          if (method_val == METHOD_BDF) {
            stepperSlave_ptr->setStepControlData(stepper);
            double slave_dt_taken = stepperSlave_ptr->takeStep(dt_taken,Rythmos::STEP_TYPE_FIXED);
            // Check that returned dt matches exactly
            TEST_FOR_EXCEPT(dt_taken != slave_dt_taken);
            Rythmos::StepStatus<double> stepStatusMaster = stepper.getStepStatus();
            Rythmos::StepStatus<double> stepStatusSlave = stepperSlave_ptr->getStepStatus();
            // Check that the stepStatus matches exactly:
            TEST_FOR_EXCEPT(stepStatusMaster.stepLETStatus != stepStatusSlave.stepLETStatus);
            TEST_FOR_EXCEPT(stepStatusMaster.stepSize != stepStatusSlave.stepSize);
            TEST_FOR_EXCEPT(stepStatusMaster.time != stepStatusSlave.time);
            if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
              *out << "Master order = " << stepStatusMaster.order << endl;
              *out << " Slave order = " << stepStatusSlave.order << endl;
              *out << "Master LET Value = " << stepStatusMaster.stepLETValue << endl;
              *out << " Slave LET Value = " << stepStatusSlave.stepLETValue << endl;
            }
            TEST_FOR_EXCEPTION(stepStatusMaster.order != stepStatusSlave.order, std::logic_error,
                "Error, stepStatusMaster.order = " << stepStatusMaster.order << 
                " != stepStatusSlave.order = " << stepStatusSlave.order << "!");
            // We will allow a difference of some multiplier of machine epsilon:
            const double
              normLETDiff = std::abs(stepStatusMaster.stepLETValue - stepStatusSlave.stepLETValue);
            TEST_FOR_EXCEPTION(
              normLETDiff > maxRestepError, std::logic_error,
              "Error, normLETDiff = " << normLETDiff << " > maxRestepError = " << maxRestepError << "!" );
            if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
              *out << "normLETDiff = " << normLETDiff << std::endl;
            }
            // Create a non-const Thyra VectorBase to use as a temp vector
            RCP<Thyra::VectorBase<double> > vec_temp = stepStatusSlave.solution->clone_v();
            // Check that the solution matches exactly
            Thyra::V_StVpStV<double>(&*vec_temp,1.0,*stepStatusMaster.solution,-1.0,*stepStatusSlave.solution);
            double normSolutionDiff = Thyra::norm_inf<double>(*vec_temp);
            if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
              *out << "normSolutionDiff = " << normSolutionDiff << endl;
            }
            const double eps = 1.0e6*Teuchos::ScalarTraits<double>::prec();
            TEST_FOR_EXCEPTION(normSolutionDiff > eps, std::logic_error,
                "Error, normSolutionDiff = " << normSolutionDiff << " > eps = " << eps << "!");
            // Check that solution dot matches exactly
            {

              RCP<const Thyra::VectorBase<double> >
                master_x_dot = get_xdot(stepper,stepStatusMaster.time),
                slave_x_dot = get_xdot(*stepperSlave_ptr,stepStatusSlave.time);
              Thyra::V_VmV<double>(&*vec_temp,*master_x_dot,*slave_x_dot);
              double normSolutionDotDiff = Thyra::norm_inf<double>(*vec_temp);
              if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
                *out << "normSolutionDotDiff = " << normSolutionDotDiff << endl;
              }
              TEST_FOR_EXCEPTION(normSolutionDotDiff > eps, std::logic_error,
                  "Error, normSolutionDotDiff = " << normSolutionDotDiff << " > eps = " << eps << "!");
            }
            // Do not check that the residual matches because the residual isn't stored in ImplicitBDFStepper.
          }
          numSteps++;
          if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) )
            stepper.describe(*out,verbLevel);
          if (dt_taken < 0)
          {
            *out << "Error, stepper failed for some reason with step taken = " << dt_taken << endl;
            break;
          }
          if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) )
          {
            // Get solution out of stepper:
            stepStatus = stepper.getStepStatus();
            x_computed_thyra_ptr = stepStatus.solution;
            // Convert Thyra::VectorBase to Epetra_Vector
            RCP<const Epetra_Vector>
              x_computed_ptr = Thyra::get_Epetra_Vector(*(epetraModel->get_x_map()),x_computed_thyra_ptr);
            if ((method_val == METHOD_BDF) && (maxOrder == 1))
            {
              int myN = x_numerical_exact.MyLength();
              if (numSteps == 1) // First step
              {
                for (int i=0 ; i<myN ; ++i)
                  x_numerical_exact[i] = x0;
              }
              for (int i=0 ; i<myN ; ++i)
                x_numerical_exact[i] = ( x_numerical_exact[i]
                    +dt_taken*epetraModel->evalR(time+dt_taken,lambda[i],coeff_s))
                    /(1-lambda[i]*dt_taken);
              for (int i=0 ; i<myN ; ++i)
                x_rel_diff[i] = (x_numerical_exact[i]-(*x_computed_ptr)[i])/x_numerical_exact[i];
              if (myN == 1)
                *out << "Computed x(" << time+dt_taken << ") = " << (*x_computed_ptr)[0] 
                    << "  Numerical Exact = " << x_numerical_exact[0] 
                    << "  Rel Diff = " << x_rel_diff[0] << std::endl;
              else
              {
                for (int i=0 ; i<myN ; ++i)
                  *out << "Computed x_" << i << "(" << time+dt_taken << ") = " << (*x_computed_ptr)[i] 
                      << "  Numerical Exact = " << x_numerical_exact[i] 
                      << "  Rel Diff = " << x_rel_diff[i] <<  std::endl;
              }
            }
            else
            {
              // compute exact answer
              RCP<const Epetra_Vector>
                x_star_ptr = epetraModel->getExactSolution(time);
              const Epetra_Vector& x_star = *x_star_ptr;
              int myN = x_computed_ptr->MyLength();
              for (int i=0 ; i<myN ; ++i)
                x_rel_diff[i] = (x_star[i]-(*x_computed_ptr)[i])/x_star[i];
              if (myN == 1)
                *out << "Computed x(" << time+dt_taken << ") = " << (*x_computed_ptr)[0] 
                    << "  Exact = " << x_star[0] 
                    << "  Rel Diff = " << x_rel_diff[0] << std::endl;
              else
              {
                for (int i=0 ; i<myN ; ++i)
                  *out << "Computed x_" << i << "(" << time+dt_taken << ") = " << (*x_computed_ptr)[i] 
                      << "  Exact = " << x_star[i] 
                      << "  Rel Diff = " << x_rel_diff[i] << std::endl;
              }
            }
          }
          time += dt_taken;
          *out << "Took stepsize of: " << dt_taken << " time = " << time << endl;
        }
        // Get solution out of stepper:
        stepStatus = stepper.getStepStatus();
        x_computed_thyra_ptr = stepStatus.solution;
      }
    }
    *out << "Integrated to time = " << time << endl;

    // Convert solution from Thyra::VectorBase to Epetra_Vector
    RCP<const Epetra_Vector>
      x_computed_ptr = Thyra::get_Epetra_Vector(*(epetraModel->get_x_map()),x_computed_thyra_ptr);
    const Epetra_Vector &x_computed = *x_computed_ptr;

    // compute exact answer
    RCP<const Epetra_Vector>
      x_star_ptr = epetraModel->getExactSolution(finalTime);
    const Epetra_Vector& x_star = *x_star_ptr;
    
    // get lambda from the problem:
    RCP<const Epetra_Vector> lambda_ptr = epetraModel->get_coeff();
    const Epetra_Vector &lambda = *lambda_ptr;

    // compute numerical exact answer (for FE and BE)
    RCP<const Epetra_Vector> x_numerical_exact_ptr; 
    if (method_val == METHOD_FE) 
    {
      RCP<Epetra_Vector> x_exact_ptr = rcp(new Epetra_Vector(x_star.Map()));
      Epetra_Vector& x_exact = *x_exact_ptr;
      int myN = x_exact.MyLength();
      for ( int i=0 ; i<myN ; ++i)
      {
        x_exact[i] = x0;
        for (int j=1 ; j<=N ; ++j)
        {
          x_exact[i] = (1+lambda[i]*dt)*x_exact[i]+dt*epetraModel->evalR(0+j*dt,lambda[i],coeff_s);
        }
        //x_exact[i] = x0*pow(1+lambda[i]*dt,N);
      }
      x_numerical_exact_ptr = x_exact_ptr;
    } 
    else if (method_val == METHOD_BE) 
    {
      RCP<Epetra_Vector> x_exact_ptr = rcp(new Epetra_Vector(x_star.Map()));
      Epetra_Vector& x_exact = *x_exact_ptr;
      int myN = x_exact.MyLength();
      for ( int i=0 ; i<myN ; ++i)
      {
        x_exact[i] = x0;
        for (int j=1 ; j<=N ; ++j)
        {
          x_exact[i] = (x_exact[i]+dt*epetraModel->evalR(0+j*dt,lambda[i],coeff_s))/(1-lambda[i]*dt);
        }
        //x_exact[i] = x0*pow(1/(1-lambda[i]*dt),N);
      }
      x_numerical_exact_ptr = x_exact_ptr;
    }
    else if (method_val == METHOD_BDF)
    {
      // exact bdf solution here?
    }

    // 06/03/05 tscoffe to get an Epetra_Map associated with an Epetra_Vector:
    // x.Map()
    // to get an Epetra_Comm associated with an Epetra_Vector:
    // x.Comm()
    
    int MyPID = x_computed.Comm().MyPID();
    if (MyPID == 0)
    {
      *out << "Integrating \\dot{x}=\\lambda x from t = " << t0 
                << " to t = " << finalTime << std::endl;
      *out << "using " << method << std::endl;
      *out << "with initial x_0 = " << x0
                << ", \\Delta t = " << dt  << "." << std::endl;
      *out << "Took " << numSteps << " steps." << std::endl;
    }
    int MyLength = x_computed.MyLength();
    if (verbose)
    {
      for (int i=0 ; i<MyLength ; ++i)
      {
        *out << std::setprecision(15);
        *out << "lambda[" << MyPID*MyLength+i << "] = " << lambda[i] << std::endl;
      }
      // Print out computed and exact solutions:
      for (int i=0 ; i<MyLength ; ++i)
      {
        *out << std::setprecision(15);
        *out << "Computed: x[" << MyPID*MyLength+i << "] = ";
        *out << std::setw(20) << x_computed[i] << "\t";
        *out << "Exact: x[" << MyPID*MyLength+i << "] = ";
        *out << std::setw(20) << x_star[i] << std::endl;
      }
    }
    
    // Check numerics against exact numerical method for FE and BE case:
    double numerical_error = 0;
    double numerical_error_mag = 0;
    if (x_numerical_exact_ptr.get())
    {
      const Epetra_Vector& x_numerical_exact = *x_numerical_exact_ptr;
      for ( int i=0 ; i<MyLength ; ++i)
      {
        if (verbose) 
        {
          *out << std::setprecision(15);
          *out << "Computed: x[" << MyPID*MyLength+i << "] = ";
          *out << std::setw(20) << x_computed[i] << "\t";
          *out << "Numerical Exact: x[" << MyPID*MyLength+i << "] = ";
          *out << std::setw(20) << x_numerical_exact[i] << std::endl;
        }
        const double thisError = x_numerical_exact[i]-x_computed[i];
        numerical_error += thisError*thisError;
        numerical_error_mag += x_numerical_exact[i]*x_numerical_exact[i];
      }
      numerical_error = sqrt(numerical_error)/sqrt(numerical_error_mag);
      result = Thyra::testMaxErr(
        "Exact numerical error",numerical_error
        ,"maxError",1.0e-10
        ,"maxWarning",1.0e-9
        ,&*out,""
        );
      if(!result) success = false;
    }

    // Check numerics against exact DE solution:
    double error = 0;
    double errorMag = 0;
    for (int i=0 ; i<MyLength ; ++i)
    {
      const double thisError = x_computed[i]-x_star[i];
      error += thisError*thisError;
      errorMag += x_star[i]*x_star[i];
    }
    error = sqrt(error)/sqrt(errorMag);
    result = Thyra::testMaxErr(
      "Exact DE solution error",error
      ,"maxError",maxError
      ,"maxWarning",10.0*maxError
      ,&*out,""
      );
    if(!result) success = false;

#ifdef HAVE_RYTHMOS_STRATIMIKOS
    // Write the final parameters to file
    if(W_factory.get())
      lowsfCreator.writeParamsFile(*W_factory);
#endif

   } // end try
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,success)

  return success ? 0 : 1;

} // end main() [Doxygen looks for this!]
