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
#include "Rythmos_TimeDiscretizedBackwardEulerModelEvaluator.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_DampenedNewtonNonlinearSolver.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_VerbosityLevelCommandLineProcessorHelpers.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_as.hpp"


#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif // HAVE_MPI


namespace {


const std::string DiagonalTransientModel_name = "DiagonalTransientModel";

const std::string DAELinearSolver_name = "DAE Linear Solver";

const std::string OverallLinearSolver_name = "Overall Linear Solver";

const std::string NonlinearSolver_name = "Nonlinear Solver";


Teuchos::RCP<const Teuchos::ParameterList>
getValidParameters()
{
  using Teuchos::RCP; using Teuchos::ParameterList;
  static RCP<const ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->sublist(DiagonalTransientModel_name).disableRecursiveValidation();
    pl->sublist(DAELinearSolver_name).disableRecursiveValidation();
    pl->sublist(OverallLinearSolver_name).disableRecursiveValidation();
    pl->sublist(NonlinearSolver_name).disableRecursiveValidation();
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
  using Thyra::createMember;
  using Thyra::createMembers;
  
  bool success = true;

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
    // A) Read commandline options
    //

    CommandLineProcessor clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);

    std::string paramsFileName = "";
    clp.setOption( "params-file", &paramsFileName,
      "File name for XML parameters" );

    std::string extraParamsString = "";
    clp.setOption( "extra-params", &extraParamsString,
      "Extra XML parameter string" );

    Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_DEFAULT;
    setVerbosityLevelOption( "verb-level", &verbLevel,
      "Top-level verbosity level.  By default, this gets deincremented as you go deeper into numerical objects.",
      &clp );

    double finalTime = 1.0;
    clp.setOption( "final-time", &finalTime, "Final time (the inital time)" );

    int numTimeSteps = 2;
    clp.setOption( "num-time-steps", &numTimeSteps, "Number of time steps" );

    bool dumpFinalSolutions = false;
    clp.setOption(
      "dump-final-solutions", "no-dump-final-solutions", &dumpFinalSolutions,
      "Determine if the final solutions are dumpped or not." );

    double maxStateError = 1e-6;
    clp.setOption( "max-state-error", &maxStateError,
      "The maximum allowed error in the integrated state in relation to the exact state solution" );
    
    // ToDo: Read in more parameters

    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;
    
    if ( Teuchos::VERB_DEFAULT == verbLevel )
      verbLevel = Teuchos::VERB_LOW;

    const Teuchos::EVerbosityLevel
      solnVerbLevel = ( dumpFinalSolutions ? Teuchos::VERB_EXTREME : verbLevel );

    //
    // B) Get the base parameter list that all other parameter lists will be
    // read from.
    //
    
    RCP<ParameterList> paramList = Teuchos::parameterList();
    if (paramsFileName.length())
      updateParametersFromXmlFile( paramsFileName, &*paramList );
    if (extraParamsString.length())
      updateParametersFromXmlString( extraParamsString, &*paramList );

    paramList->validateParameters(*getValidParameters());

    //
    // C) Create the Stratimikos linear solver factories.
    //

    // Get the linear solve strategy that will be used to solve for the linear
    // system with the dae's W matrix.
    Stratimikos::DefaultLinearSolverBuilder daeLinearSolverBuilder;
    daeLinearSolverBuilder.setParameterList(sublist(paramList,DAELinearSolver_name));
    RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >
      daeLOWSF = createLinearSolveStrategy(daeLinearSolverBuilder);

    // Get the linear solve strategy that can be used to override the overall
    // linear system solve
    Stratimikos::DefaultLinearSolverBuilder overallLinearSolverBuilder;
    overallLinearSolverBuilder.setParameterList(sublist(paramList,OverallLinearSolver_name));
    RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >
      overallLOWSF = createLinearSolveStrategy(overallLinearSolverBuilder);
    
    //
    // D) Create the underlying EpetraExt::ModelEvaluator
    //

    RCP<EpetraExt::DiagonalTransientModel> epetraDaeModel =
      EpetraExt::diagonalTransientModel(
        epetra_comm,
        sublist(paramList,DiagonalTransientModel_name)
        );

    *out <<"\nepetraDaeModel valid options:\n";
    epetraDaeModel->getValidParameters()->print(
      *out, PLPrintOptions().indent(2).showTypes(true).showDoc(true)
      );
    
    //
    // E) Create the Thyra-wrapped ModelEvaluator
    //
    
    RCP<Thyra::ModelEvaluator<double> > daeModel =
      epetraModelEvaluator(epetraDaeModel,daeLOWSF);

    //
    // F) Create the TimeDiscretizedBackwardEulerModelEvaluator
    //
    
    MEB::InArgs<Scalar> initCond = daeModel->createInArgs();
    initCond.setArgs(daeModel->getNominalValues());

    RCP<Thyra::ModelEvaluator<Scalar> >
      discretizedModel = Rythmos::timeDiscretizedBackwardEulerModelEvaluator<Scalar>(
        daeModel, initCond, finalTime, numTimeSteps, overallLOWSF );
    
    *out << "\ndiscretizedModel = " << describe(*discretizedModel,verbLevel);

    //
    // F) Setup a nonlinear solver and solve the system
    //

    // F.1) Setup a nonlinear solver

    Thyra::DampenedNewtonNonlinearSolver<Scalar> nonlinearSolver;
    nonlinearSolver.setOStream(out);
    nonlinearSolver.setVerbLevel(verbLevel);
    //nonlinearSolver.setParameterList(sublist(paramList,NonlinearSolver_name));
    //2007/11/27: rabartl: ToDo: Implement parameter list handling for
    //DampenedNonlinearSolve so that I can uncomment the above line.
    nonlinearSolver.setModel(discretizedModel);

    // F.2) Solve the system

    RCP<Thyra::VectorBase<Scalar> > 
      x_bar = createMember(discretizedModel->get_x_space());
    V_S( x_bar.ptr(), 0.0 );

    Thyra::SolveStatus<Scalar> solveStatus =
      Thyra::solve( nonlinearSolver, &*x_bar );

    *out << "\nsolveStatus:\n" << solveStatus;

    *out << "\nx_bar = " << describe(*x_bar,solnVerbLevel);
    
    //
    // G) Verify that the solution is correct???
    //

    // Check against the end time exact solution.

    RCP<const Thyra::VectorBase<Scalar> >
      exact_x_final = Thyra::create_Vector(
        epetraDaeModel->getExactSolution(finalTime),
        daeModel->get_x_space()
        );

    RCP<const Thyra::VectorBase<Scalar> > solved_x_final
      = rcp_dynamic_cast<Thyra::ProductVectorBase<Scalar> >(x_bar,true)->getVectorBlock(numTimeSteps-1);
    
    const bool result = Thyra::testRelNormDiffErr(
      "exact_x_final", *exact_x_final, "solved_x_final", *solved_x_final,
      "maxStateError", maxStateError, "warningTol", 1.0, // Don't warn
      &*out, solnVerbLevel
      );
    if (!result) success = false;

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,success);

  if(success)
    *out << "\nEnd Result: TEST PASSED" << endl;
  else
    *out << "\nEnd Result: TEST FAILED" << endl;
  
  return ( success ? 0 : 1 );

} // end main() [Doxygen looks for this!]

