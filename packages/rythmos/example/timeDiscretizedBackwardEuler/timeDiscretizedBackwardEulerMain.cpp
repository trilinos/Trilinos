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
#include "Thyra_DefaultRealLinearSolverBuilder.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
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
    Thyra::DefaultRealLinearSolverBuilder daeLinearSolverBuilder;
    daeLinearSolverBuilder.setParameterList(sublist(paramList,DAELinearSolver_name));
    RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >
      daeLOWSF = createLinearSolveStrategy(daeLinearSolverBuilder);

    // Get the linear solve strategy that can be used to override the overall
    // linear system solve
    Thyra::DefaultRealLinearSolverBuilder overallLinearSolverBuilder;
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
    // F) Create the RHS and solve the overall linear system (which solves the
    // forward problem)
    //

    // We know that this is a linear model so we need to just done one linear
    // solve and we are done!

    // F.1) Setup the input and output objects for (linear) f_bar(x_bar) = 0

    RCP<Thyra::VectorBase<Scalar> > 
      x_bar = createMember(discretizedModel->get_x_space());
    V_S( &*x_bar, 0.0 );
 
    RCP<Thyra::VectorBase<Scalar> > 
      f_bar = createMember(discretizedModel->get_f_space());
 
    RCP<Thyra::LinearOpWithSolveBase<Scalar> > 
      W_bar = discretizedModel->create_W();

    MEB::InArgs<Scalar> inArgs_bar = discretizedModel->createInArgs();
    inArgs_bar.set_x(x_bar);

    MEB::OutArgs<Scalar> outArgs_bar = discretizedModel->createOutArgs();
    outArgs_bar.set_f(f_bar);
    outArgs_bar.set_W(W_bar);
    
    discretizedModel->evalModel(inArgs_bar,outArgs_bar);

    // Solve linear system!

    Thyra::SolveStatus<Scalar> solveStatus =
      Thyra::solve( *W_bar, Thyra::NOTRANS, *f_bar, &*x_bar );

    *out << "\nsolveStatus:\n" << solveStatus;

    *out << "\nx_bar = " << describe(*x_bar,solnVerbLevel);
    
    //
    // G) Verify that the solution is correct???
    //

    TEST_FOR_EXCEPT("ToDo: Implement!");
    
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,success);

  if(success)
    *out << "\nEnd Result: TEST PASSED" << endl;
  else
    *out << "\nEnd Result: TEST FAILED" << endl;
  
  return ( success ? 0 : 1 );

} // end main() [Doxygen looks for this!]

