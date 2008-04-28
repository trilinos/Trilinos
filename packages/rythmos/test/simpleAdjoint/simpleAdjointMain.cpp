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
#include "Rythmos_AdjointModelEvaluator.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"
#include "Rythmos_DefaultIntegrator.hpp"
#include "Rythmos_SimpleIntegrationControlStrategy.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
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


Teuchos::RCP<const Teuchos::ParameterList>
getValidParameters()
{
  using Teuchos::RCP; using Teuchos::ParameterList;
  static RCP<const ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->sublist(TimeStepNonlinearSolver_name).disableRecursiveValidation();
    pl->sublist(Stratimikos_name).disableRecursiveValidation();
    pl->sublist(DiagonalTransientModel_name).disableRecursiveValidation();
    pl->sublist(RythmosStepper_name).disableRecursiveValidation();
    pl->sublist(RythmosIntegrator_name).disableRecursiveValidation();
    validPL = pl;
  }
  return validPL;
}


} // namespace


//
// Problem:
//
//  Terminal response:
//
//     h(x(t_final)) = 1/2 * x^T * x
//
//  Subject to ODE:
//
//     x_dot = alpha * x, for t in [0,t_final]
//      x(0) = p
//
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

    double finalTime = 1e-3;
    clp.setOption( "final-time", &finalTime,
      "Final integration time (initial time is 0.0)" );

    int numTimeSteps = 10;
    clp.setOption( "num-time-steps", &numTimeSteps,
      "Number of (fixed) time steps.  If <= 0.0, then variable time steps are taken" );

    Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_DEFAULT;
    setVerbosityLevelOption( "verb-level", &verbLevel,
      "Top-level verbosity level.  By default, this gets deincremented as you go deeper into numerical objects.",
      &clp );

    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    //
    // A) Get the base parameter list that all other parameter lists will be
    // read from.
    //
    
    RCP<ParameterList>
      paramList = Teuchos::parameterList();
    if (paramsFileName.length())
      updateParametersFromXmlFile( paramsFileName, &*paramList );

    paramList->validateParameters(*getValidParameters());

    //
    // B) Create the Stratimikos linear solver factory.
    //
    // This is the linear solve strategy that will be used to solve for the
    // linear system with the W.
    //

    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
    linearSolverBuilder.setParameterList(sublist(paramList,Stratimikos_name));
    RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >
      W_factory = createLinearSolveStrategy(linearSolverBuilder);

    //
    // C) Create and initalize the forward model
    //
    
    // C.1) Create the underlying EpetraExt::ModelEvaluator

    RCP<EpetraExt::DiagonalTransientModel>
      epetraStateModel = EpetraExt::diagonalTransientModel(
        epetra_comm,
        sublist(paramList,DiagonalTransientModel_name)
        );

    *out <<"\nepetraStateModel valid options:\n";
    epetraStateModel->getValidParameters()->print(
      *out, PLPrintOptions().indent(2).showTypes(true).showDoc(true)
      );
    
    // C.2) Create the Thyra-wrapped ModelEvaluator
    
    RCP<Thyra::ModelEvaluator<double> >
      fwdStateModel = epetraModelEvaluator(epetraStateModel,W_factory);

    //
    // D) Create the Adjoint ME wrapper object
    //

    RCP<Thyra::ModelEvaluator<double> >
      adjointModel = Rythmos::adjointModelEvaluator<double>(fwdStateModel);


    //
    // E) Create a stepper for the adjoint
    //

    TEST_FOR_EXCEPT(true);

    //
    // F) Set up the initial condition for the adjoint at the final time
    //

    TEST_FOR_EXCEPT(true);

    //
    // G) Integrate the adjoint backwards in time (using inverse time)
    //

    TEST_FOR_EXCEPT(true);

    //
    // H) Test the final adjoint againt exact discrete solution.
    //

    TEST_FOR_EXCEPT(true);

    
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,success);

  if(success)
    *out << "\nEnd Result: TEST PASSED" << endl;
  else
    *out << "\nEnd Result: TEST FAILED" << endl;
  
  return ( success ? 0 : 1 );

} // end main() [Doxygen looks for this!]

