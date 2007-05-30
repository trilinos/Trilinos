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
#include "Rythmos_ForwardSensitivityStepper.hpp"
#include "Thyra_DefaultRealLinearSolverBuilder.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"
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
#  include "mpi.h"
#else
#  include "Epetra_SerialComm.h"
#endif // HAVE_MPI

namespace {

const std::string TimeStepNonlinearSolver_name = "TimeStepNonlinearSolver";

const std::string Stratimikos_name = "Stratimikos";

const std::string DiagonalTransientModel_name = "DiagonalTransientModel";

const std::string RythmosStepper_name = "Rythmos Stepper";

Teuchos::RefCountPtr<const Teuchos::ParameterList>
getValidParameters()
{
  using Teuchos::RefCountPtr; using Teuchos::ParameterList;
  static RefCountPtr<const ParameterList> validPL;
  if (is_null(validPL)) {
    RefCountPtr<ParameterList> pl = Teuchos::parameterList();
    pl->sublist(TimeStepNonlinearSolver_name);
    pl->sublist(Stratimikos_name);
    pl->sublist(DiagonalTransientModel_name);
    pl->sublist(RythmosStepper_name);
    validPL = pl;
  }
  return validPL;
}


} // namespace

int main(int argc, char *argv[])
{

  typedef double Scalar;
  typedef double ScalarMag;
  using Teuchos::describe;
  using Teuchos::RefCountPtr;
  using Teuchos::rcp;
  using Teuchos::as;
  using Teuchos::ParameterList;
  using Teuchos::CommandLineProcessor;
  typedef Teuchos::ParameterList::PrintOptions PLPrintOptions;
  typedef Thyra::ModelEvaluatorBase MEB;
  
  bool result, success = true;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  RefCountPtr<Epetra_Comm> epetra_comm;
#ifdef HAVE_MPI
  epetra_comm = rcp( new Epetra_MpiComm(MPI_COMM_WORLD) );
#else
  epetra_comm = rcp( new Epetra_SerialComm );
#endif // HAVE_MPI

  RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {

    //
    // Read commandline options
    //

    CommandLineProcessor clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);

    std::string linearSolverParamsFile = "";
    clp.setOption( "linear-solver-params-file", &linearSolverParamsFile,
      "File name for XML linear solver parameters for Stratimikos" );

    std::string linearSolverExtraParams = "";
    clp.setOption( "linear-solver-extra-params", &linearSolverExtraParams,
      "Extra XML parameter list string for linear solver parameters for Stratimikos" );

    double maxError = 1e-6;
    clp.setOption( "max-error", &maxError,
      "The maximum allowed error ???" );

    double finalTime = 1.0;
    clp.setOption( "final-time", &finalTime,
      "Final integration time (initial time is 0.0)" );

    int numTimeSteps = 10;
    clp.setOption( "num-time-steps", &numTimeSteps,
      "Number of (fixed) time steps." );

    bool doFwdSensSolve = false;
    clp.setOption( "fwd-sens-solve", "state-solve", &doFwdSensSolve,
      "Do the forward sensitivity solve or just the state solve" );

    Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_DEFAULT;
    setVerbosityLevelOption( "verb-level", &verbLevel,
      "Overall verbosity level.", &clp );

    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;
    
    if ( Teuchos::VERB_DEFAULT == verbLevel )
      verbLevel = Teuchos::VERB_LOW;

    //
    // Get the base parameter list that all other parameter lists will be read
    // from.
    //
    
    RefCountPtr<ParameterList>
      paramList = Teuchos::parameterList();
    if (linearSolverParamsFile.length())
      updateParametersFromXmlFile( linearSolverParamsFile, &*paramList );
    if (linearSolverExtraParams.length())
      updateParametersFromXmlString( linearSolverExtraParams, &*paramList );
    paramList->validateParameters(*getValidParameters(),0); // Only validate top level lists!

    //
    // Create the Stratimikos linear solver factory.
    //
    // This is the linear solve strategy that will be used to solve for the
    // linear system with the W.
    //

    Thyra::DefaultRealLinearSolverBuilder linearSolverBuilder;
    linearSolverBuilder.setParameterList(sublist(paramList,Stratimikos_name));
    RefCountPtr<Thyra::LinearOpWithSolveFactoryBase<Scalar> >
      W_factory = createLinearSolveStrategy(linearSolverBuilder);
    
    //
    // Create the underlying EpetraExt::ModelEvaluator
    //

    RefCountPtr<EpetraExt::DiagonalTransientModel>
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
    
    RefCountPtr<Thyra::ModelEvaluator<double> >
      stateModel = epetraModelEvaluator(epetraStateModel,W_factory);
    
    //
    // Create the Rythmos stateStepper
    //

    RefCountPtr<Rythmos::TimeStepNonlinearSolver<double> >
      nonlinearSolver = Teuchos::rcp(new Rythmos::TimeStepNonlinearSolver<double>());
    RefCountPtr<ParameterList>
      nonlinearSolverPL = sublist(paramList,TimeStepNonlinearSolver_name);
    nonlinearSolverPL->get("Default Tol",1e-3*maxError); // Set default if not set
    nonlinearSolver->setParameterList(nonlinearSolverPL);

    RefCountPtr<Rythmos::StepperBase<Scalar> > stateStepper;
    
    if (doFwdSensSolve) {
      stateStepper = rcp(new Rythmos::BackwardEulerStepper<double>());
      // The model and nonlinear solver will be set in the forward stepper!
    }
    else {
      stateStepper = rcp(
        new Rythmos::BackwardEulerStepper<double>(
          stateModel, nonlinearSolver
          )
        );
    }
   
    // ToDo: Above, add a factory interface to create the time stateStepper
    // strategy given a parameter sublist (kind of like with Stratimikos).

    *out <<"\nstateStepper:\n" << describe(*stateStepper,verbLevel);
    *out <<"\nstateStepper valid options:\n";
    stateStepper->getValidParameters()->print(
      *out, PLPrintOptions().indent(2).showTypes(true).showDoc(true)
      );

    stateStepper->setParameterList(sublist(paramList,RythmosStepper_name));

    // The default problem is just the state stepper
    RefCountPtr<Rythmos::StepperBase<Scalar> > stepper = stateStepper;

    //
    // Create the forward sensitivity stepper if needed
    //

    RefCountPtr<Rythmos::ForwardSensitivityStepper<Scalar> > stateAndSensStepper;

    if (doFwdSensSolve) {

      stateAndSensStepper = Rythmos::forwardSensitivityStepper<Scalar>(
        stateModel, 0, stateModel->getNominalValues(),
        stateStepper, nonlinearSolver
        );
      // The above call will result in stateStepper and nonlinearSolver being
      // cloned.  This helps to ensure consistency between the state and
      // sensitivity computations!
      stepper = stateAndSensStepper;

      // Set the initial condition for the forward sensitivities
      RefCountPtr<Thyra::VectorBase<Scalar> > s_bar_init
        = createMember(stateAndSensStepper->getFwdSensModel()->get_x_space());

      assign( &*s_bar_init, 0.0 );
      // WARNING!!!!!!! You have to actually set the right initial condition
      // but you can do that later.

      stateAndSensStepper->setFwdSensInitialCondition(s_bar_init);

    }

    *out <<"\nstepper:\n" << describe(*stepper,verbLevel);
    *out <<"\nstepper->get_s_space():\n" << describe(*stepper->get_x_space(),verbLevel);
    
    //
    // Time step through the problem taking constant time steps, testing the
    // solution as you go ...
    //

    double t0 = 0.0;
    double dt = (finalTime-t0)/numTimeSteps;
    double time = t0;

    stepper->setVerbLevel(verbLevel);

    for ( int i = 1 ; i <= numTimeSteps ; ++i ) {

     if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_MEDIUM) )
       *out << "\ntime step = " << i << ", time = " << time << ":\n";

     Teuchos::OSTab tab(out);

     if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_MEDIUM) )
       *out << "\nTaking step ...\n";

     double dt_taken = stepper->takeStep(dt,Rythmos::FIXED_STEP);

     TEST_FOR_EXCEPTION(
       dt_taken != dt, std::logic_error,
       "Error, stepper took step of dt = " << dt_taken 
       << " when asked to take step of dt = " << dt << "\n"
       );

     time += dt_taken;

     Rythmos::StepStatus<Scalar>
       stepStatus = stepper->getStepStatus();

     RefCountPtr<const Thyra::VectorBase<Scalar> >
       solution = stepStatus.solution,
       solutionDot = stepStatus.solutionDot;

     *out << "\nsolution = \n" << describe(*solution,verbLevel);
     *out << "\nsolutionDot = \n" << describe(*solutionDot,verbLevel);

     RefCountPtr<const Thyra::VectorBase<Scalar> > exact_x, solved_x;
     ScalarMag rel_err;

     exact_x = create_Vector(
       epetraStateModel->getExactSolution(time),
       stateModel->get_x_space()
       );
     
     if (doFwdSensSolve ) {
       solved_x = productVectorBase(solution)->getVectorBlock(0).assert_not_null();
     }
     else {
       solved_x = get_x(*stepper,time);
     }

     rel_err = relErr(*exact_x,*solved_x);
     
     result = rel_err <= maxError;
     
     if (!result) success = false;
     
     if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_MEDIUM) ) {
        *out << "\n||exact_x||inf = " << norm_inf(*exact_x);
        *out << ", ||solved_x||inf = " << norm_inf(*solved_x);
        *out
          << ", relErr(exact_x,solved_x) = " << rel_err
          << " <= " << maxError << " : " << Thyra::passfail(result) << "\n";
        if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
          *out << "\nexact_x =\n" << describe(*exact_x,verbLevel);
          *out << "\nsolved_x =\n" << describe(*solved_x,verbLevel);
        }
     }

    }
    
    if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW) )
      *out << "\nFinal time = " << time << "\n";
    
    Teuchos::OSTab tab(out);

    RefCountPtr<const Thyra::VectorBase<Scalar> >
      solution = get_x(*stepper,time);

    RefCountPtr<const Thyra::VectorBase<Scalar> > exact_x, solved_x;
    ScalarMag rel_err;

    exact_x = create_Vector(
      epetraStateModel->getExactSolution(time),
      stateModel->get_x_space()
      );
    
    if (doFwdSensSolve ) {
      solved_x = productVectorBase(solution)->getVectorBlock(0);
    }
    else {
      solved_x = solution;
    }
    
    rel_err = relErr(*exact_x,*solved_x);
      
    result = rel_err <= maxError;
    
    if (!result) success = false;

    if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW) ) {
    
      *out << "\n||exact_x||inf = " << norm_inf(*exact_x);
      *out << ", ||solved_x||inf = " << norm_inf(*solved_x);
      *out
        << ", relErr(exact_x,solved_x) = " << rel_err
        << " <= " << maxError << " : " << Thyra::passfail(result) << "\n";
      if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
        *out << "\nexact_x =\n" << describe(*exact_x,verbLevel);
        *out << "\nsolved_x =\n" << describe(*solved_x,verbLevel);
      }

    }
    
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,success)

  if(success)
    *out << "\nEnd Result: TEST PASSED" << endl;
  else
    *out << "\nEnd Result: TEST FAILED" << endl;
  
  return ( success ? 0 : 1 );

} // end main() [Doxygen looks for this!]
