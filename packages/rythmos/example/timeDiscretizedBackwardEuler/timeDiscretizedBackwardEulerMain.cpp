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


#include "Rythmos_BackwardEulerStepper.hpp"
#include "Rythmos_LinearTimeInvariantModelEvaluator.hpp"
#include "Rythmos_TimeDiscretizedBackwardEulerModelEvaluator.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
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
  
  bool result, success = true;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {

    //
    // A) Read commandline options
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
    if (linearSolverParamsFile.length())
      updateParametersFromXmlFile( linearSolverParamsFile, &*paramList );
    if (linearSolverExtraParams.length())
      updateParametersFromXmlString( linearSolverExtraParams, &*paramList );

    // ToDo: Validate the parameter list

    Thyra::DefaultRealLinearSolverBuilder linearSolverBuilder;
    linearSolverBuilder.setParameterList(paramList);
    RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >
      W_factory = createLinearSolveStrategy(linearSolverBuilder);

    //
    *out << "\nC) Get the data for the problem:\n\n";
    //
    // Get the defining matrices M, K, delta_t, initial conditions etc ...
    //

    // 2007/11/14: rabartl: Initially we will just fake M, K etc ...

    // Create space for LTI matrices
    const int n = 2;
    const RCP<const Thyra::VectorSpaceBase<Scalar> >
      space = Thyra::defaultSpmdVectorSpace<Scalar>(n);
    
    // Intially just create random square multi-vectors for M and K
    RCP<Thyra::MultiVectorBase<Scalar> >
      mvM = createMembers(space,n,"M"),
      mvK = createMembers(space,n,"K");
    assign( &*mvM, 1.0 );
    assign( &*mvK, 2.0 );
    RCP<const Thyra::LinearOpBase<Scalar> >
      M = mvM,
      K = mvK;
    
    *out << "\nM = " << describe(*M,verbLevel);
    *out << "\nK = " << describe(*K,verbLevel);

    //
    *out << "\nD) Create the LinearTimeInvariantModelEvaluator DAE ...\n\n";
    //

    RCP<Thyra::ModelEvaluator<Scalar> >
      ltiModel = Rythmos::linearTimeInvariantModelEvaluator(M,K,W_factory);

    *out << "\nltiModel = " << describe(*ltiModel,verbLevel);

    //
    // E) Create the TimeDiscretizedBackwardEulerModelEvaluator
    //
    
    // E.1) Create the initial condition

    RCP<Thyra::VectorBase<Scalar> >
      x_init = createMember(ltiModel->get_x_space()),
      x_dot_init = createMember(ltiModel->get_x_space());

    assign( &*x_init, 1.0 );
    assign( &*x_dot_init, 0.0 );

    MEB::InArgs<Scalar> initCond = ltiModel->createInArgs();
    initCond.set_x(x_init);
    initCond.set_x_dot(x_dot_init);
    initCond.set_t(0.0);

    RCP<Thyra::ModelEvaluator<Scalar> >
      discretizedModel = Rythmos::timeDiscretizedBackwardEulerModelEvaluator<Scalar>(
        ltiModel, initCond, finalTime, numTimeSteps );
    
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

    Thyra::solve( *W_bar, Thyra::NOTRANS, *f_bar, &*x_bar );

    *out << "\nx_bar = " << describe(*x_bar,verbLevel);
    
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

