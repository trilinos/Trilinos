// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_DefaultInverseLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_TestingTools.hpp"
#include "EpetraExt_CrsMatrixIn.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_TimeMonitor.hpp"
#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif


/**
 * Example program that shows one way to create a physics-based
 * preconditioner using the building blocks of Stratimikos and Thyra
 * implicit linear operators.
 *
 * See the string in setDocString(...) for more detailed documentation ...
 */


namespace {


// Read a Epetra_CrsMatrix from a Matrix market file
Teuchos::RCP<Epetra_CrsMatrix>
readEpetraCrsMatrixFromMatrixMarket(
  const std::string fileName, const Epetra_Comm &comm
  )
{
  Epetra_CrsMatrix *A = 0;
  TEUCHOS_TEST_FOR_EXCEPTION(
    0!=EpetraExt::MatrixMarketFileToCrsMatrix( fileName.c_str(), comm, A ),
    std::runtime_error,
    "Error, could not read matrix file \""<<fileName<<"\"!"
    );
  return Teuchos::rcp(A);
}


// Read an Epetra_CrsMatrix in as a wrapped Thyra::EpetraLinearOp object
Teuchos::RCP<const Thyra::LinearOpBase<double> >
readEpetraCrsMatrixFromMatrixMarketAsLinearOp(
  const std::string fileName, const Epetra_Comm &comm,
  const std::string &label
  )
{
  return Thyra::epetraLinearOp(
    readEpetraCrsMatrixFromMatrixMarket(fileName,comm),label
    );
}


} // namespace


int main(int argc, char* argv[])
{

  using Teuchos::describe;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcp_const_cast;
  using Teuchos::RCP;
  using Teuchos::CommandLineProcessor;
  using Teuchos::ParameterList;
  using Teuchos::sublist;
  typedef ParameterList::PrintOptions PLPrintOptions;
  using Thyra::inverse;
  using Thyra::initializePreconditionedOp;
  using Thyra::initializeOp;
  using Thyra::unspecifiedPrec;
  using Thyra::solve;
  typedef RCP<const Thyra::LinearOpBase<double> > LinearOpPtr;
  typedef RCP<Thyra::VectorBase<double> > VectorPtr;

  bool success = true;
  bool verbose = true;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {

    //
    // Read in options from the command line
    //

    const int numVerbLevels = 6;
    Teuchos::EVerbosityLevel
      verbLevelValues[] =
      {
        Teuchos::VERB_DEFAULT, Teuchos::VERB_NONE,
        Teuchos::VERB_LOW, Teuchos::VERB_MEDIUM,
        Teuchos::VERB_HIGH, Teuchos::VERB_EXTREME
      };
    const char*
      verbLevelNames[] =
      { "default", "none", "low", "medium", "high", "extreme" };

    Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_MEDIUM;
    std::string paramsFile = "";
    std::string extraParamsFile = "";
    std::string baseDir = ".";
    bool useP1Prec = true;
    bool invertP1 = false;
    bool showParams = false;
    double solveTol = 1e-8;

    CommandLineProcessor clp(false); // Don't throw exceptions

    clp.setOption( "verb-level", &verbLevel,
      numVerbLevels, verbLevelValues, verbLevelNames,
      "Verbosity level used for all objects."
      );
    clp.setOption( "params-file", &paramsFile,
      "File containing parameters in XML format.",
      true // Required
      );
    clp.setOption( "extra-params-file", &extraParamsFile,
      "File containing extra parameters in XML format."
      );
    clp.setOption( "base-dir", &baseDir,
      "Base directory for all of the files."
      );
    clp.setOption( "use-P1-prec", "use-algebraic-prec", &useP1Prec,
      "Use the physics-based preconditioner for P2 based on P1 or just an algebraic preconditioner"
      );
    clp.setOption( "invert-P1", "prec-P1-only", &invertP1,
      "In the physics-based preconditioner, use P1's preconditioner or fully invert P1."
      );
    clp.setOption( "solve-tol", &solveTol,
      "Tolerance for the solution to determine success or failure!"
      );
    clp.setOption( "show-params", "no-show-params", &showParams,
      "Show the parameter list or not."
      );
    clp.setDocString(
      "This example program shows an attempt to create a physics-based\n"
      "preconditioner using the building blocks of Stratimikos and Thyra\n"
      "implicit linear operators.  The idea is to use the linear operator\n"
      "for first-order Lagrange finite elements as the preconditioner for the\n"
      "operator using second-order Lagrange finite elements.  The details of the\n"
      "PDE being represented, the mesh being used, or the boundary conditions are\n"
      "beyond the scope of this example.\n"
      "\n"
      "The matrices used in this example are:\n"
      "\n"
      "  P2: The discretized PDE matrix using second-order Lagrange finite elements.\n"
      "  P1: The discretized PDE matrix using first-order Lagrange finite elements.\n"
      "  M22: The mass matrix for the second-order Lagrange finite-element basis functions.\n"
      "  M11: The mass matrix for the first-order Lagrange finite-element basis functions.\n"
      "  M21: A rectangular matrix that uses mixed first- and second-order basis functions\n"
      "       to map from P2 space to P1 space.\n"
      "  M12: A rectangular matrix that uses mixed first- and second-order basis functions\n"
      "       to map from P1 space to P2 space.\n"
      "\n"
      "The above matrices are read from Matrix Market *.mtx files in the directory given by\n"
      "the --base-dir command-line option.\n"
      "\n"
      "The preconditioner operator created in this example program is:\n"
      "\n"
      "   precP2Op = (inv(M22) * M12) * prec(P1) * (inv(M11) * M21)\n"
      "\n"
      "where prec(P1) is either the algebraic preconditioner for P1 (--prec-P1-only)\n"
      "or is a full solve for P1 (--invert-P1).\n"
      "\n"
      "We use Stratimikos to specify the linear solvers and/or algebraic\n"
      "preconditioners and we use the Thyra implicit operators to build the\n"
      "implicitly multiplied linear operators associated with the preconditioner.\n"
      "\n"
      "Warning!  This physics-based preconditioner is singular and can not\n"
      "be used to solve the linear system given a random right-hand side.  However,\n"
      "singular or not, this example shows how one can use Thyra/Stratimikos to quickly\n"
      "try out these types of preconditioning ideas (even if they do not work).\n"
      );

    // Note: Use --help on the command line to see the above documentation

    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;


    //
    *out << "\nA) Reading in Epetra_CrsMatrix objects for P1, P2, M11, M12, M21, and M22 ...\n";
    //

#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif

    LinearOpPtr P1=readEpetraCrsMatrixFromMatrixMarketAsLinearOp(
      baseDir+"/P1.mtx",comm,"P1");
    *out << "\nP1 = " << describe(*P1,verbLevel) << "\n";
    LinearOpPtr P2= readEpetraCrsMatrixFromMatrixMarketAsLinearOp(
      baseDir+"/P2.mtx",comm,"P2");
    *out << "\nP2 = " << describe(*P2,verbLevel) << "\n";
    LinearOpPtr M11=readEpetraCrsMatrixFromMatrixMarketAsLinearOp(
      baseDir+"/M11.mtx",comm,"M11");
    *out << "\nM11 = " << describe(*M11,verbLevel) << "\n";
    LinearOpPtr M22=readEpetraCrsMatrixFromMatrixMarketAsLinearOp(
      baseDir+"/M22.mtx",comm,"M22");
    *out << "\nM22 = " << describe(*M22,verbLevel) << "\n";
    LinearOpPtr M12=readEpetraCrsMatrixFromMatrixMarketAsLinearOp(
      baseDir+"/M12.mtx",comm,"M12");
    *out << "\nM12 = " << describe(*M12,verbLevel) << "\n";
    LinearOpPtr M21=readEpetraCrsMatrixFromMatrixMarketAsLinearOp(
      baseDir+"/M21.mtx",comm,"M21");
    *out << "\nM21 = " << describe(*M21,verbLevel) << "\n";

    // ToDo: Replace the above functions with a general Thyra strategy object
    // to do the reading

    //
    *out << "\nB) Get the preconditioner and/or linear solver strategies to invert M11, M22, P1, and P2 ...\n";
    //

    //
    // Get separate parameter sublists for each square operator separately
    // that specify the type of linear solver and/or preconditioner to use.
    //

    RCP<ParameterList> paramList =
      Teuchos::getParametersFromXmlFile( baseDir+"/"+paramsFile );
    if (extraParamsFile.length()) {
      Teuchos::updateParametersFromXmlFile( baseDir+"/"+extraParamsFile, paramList.ptr() );
    }
    if (showParams) {
      *out << "\nRead in parameter list:\n\n";
      paramList->print(*out,PLPrintOptions().indent(2).showTypes(true));
    }

    Stratimikos::DefaultLinearSolverBuilder M11_linsolve_strategy_builder;
    M11_linsolve_strategy_builder.setParameterList(
      sublist(paramList,"M11 Solver",true) );

    Stratimikos::DefaultLinearSolverBuilder M22_linsolve_strategy_builder;
    M22_linsolve_strategy_builder.setParameterList(
      sublist(paramList,"M22 Solver",true) );

    Stratimikos::DefaultLinearSolverBuilder P1_linsolve_strategy_builder;
    P1_linsolve_strategy_builder.setParameterList(
      sublist(paramList,"P1 Solver",true) );

    Stratimikos::DefaultLinearSolverBuilder P2_linsolve_strategy_builder;
    P2_linsolve_strategy_builder.setParameterList(
      sublist(paramList,"P2 Solver",true) );

    //
    // Create the linear solver and/or preconditioner strategies
    // (i.e. factories)
    //

    // For M11 and M22, we want full linear solver factories with embedded
    // algebraic preconditioner factories.

    RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > M11_linsolve_strategy
      = createLinearSolveStrategy(M11_linsolve_strategy_builder);

    RCP<Thyra::LinearOpWithSolveFactoryBase<double> > M22_linsolve_strategy
      = createLinearSolveStrategy(M22_linsolve_strategy_builder);

    // For P1, we only want its preconditioner factory.

    RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > P1_linsolve_strategy;
    RCP<const Thyra::PreconditionerFactoryBase<double> > P1_prec_strategy;
    if(invertP1)
      P1_linsolve_strategy
        = createLinearSolveStrategy(P1_linsolve_strategy_builder);
    else
      P1_prec_strategy
        = createPreconditioningStrategy(P1_linsolve_strategy_builder);

    // For P2, we only want a linear solver factory.  We will supply the
    // preconditioner ourselves (that is the whole point of this example).

    RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > P2_linsolve_strategy
      = createLinearSolveStrategy(P2_linsolve_strategy_builder);

    //
    *out << "\nC) Create the physics-based preconditioner! ...\n";
    //

    *out << "\nCreating inv(M11) ...\n";
    LinearOpPtr invM11 = inverse(*M11_linsolve_strategy, M11);
    *out << "\ninvM11 = " << describe(*invM11,verbLevel) << "\n";

    *out << "\nCreating inv(M22) ...\n";
    LinearOpPtr invM22 = inverse(*M22_linsolve_strategy, M22);
    *out << "\ninvM22 = " << describe(*invM22,verbLevel) << "\n";

    *out << "\nCreating prec(P1) ...\n";
    LinearOpPtr invP1;
    if(invertP1) {
      *out << "\nCreating prec(P1) as a full solver ...\n";
      invP1 = inverse(*P1_linsolve_strategy, P1);
    }
    else {
      *out << "\nCreating prec(P1) as just an algebraic preconditioner ...\n";
      RCP<Thyra::PreconditionerBase<double> >
        precP1 = prec(*P1_prec_strategy,P1);
      *out << "\nprecP1 = " << describe(*precP1,verbLevel) << "\n";
      invP1 = precP1->getUnspecifiedPrecOp();
    }
    rcp_const_cast<Thyra::LinearOpBase<double> >(
      invP1)->setObjectLabel("invP1"); // Cast to set label ...
    *out << "\ninvP1 = " << describe(*invP1,verbLevel) << "\n";

    LinearOpPtr P2ToP1 = multiply( invM11, M21 );
    *out << "\nP2ToP1 = " << describe(*P2ToP1,verbLevel) << "\n";

    LinearOpPtr P1ToP2 = multiply( invM22, M12 );
    *out << "\nP1ToP2 = " << describe(*P1ToP2,verbLevel) << "\n";

    LinearOpPtr precP2Op = multiply( P1ToP2, invP1, P2ToP1 );
    *out << "\nprecP2Op = " << describe(*precP2Op,verbLevel) << "\n";

    //
    *out << "\nD) Setup the solver for P2 ...\n";
    //

    RCP<Thyra::LinearOpWithSolveBase<double> >
      P2_lows = P2_linsolve_strategy->createOp();
    if(useP1Prec) {
      *out << "\nCreating the solver P2 using the specialized precP2Op\n";
      initializePreconditionedOp<double>( *P2_linsolve_strategy, P2,
        unspecifiedPrec(precP2Op), P2_lows.ptr());
    }
    else {
      *out << "\nCreating the solver P2 using algebraic preconditioner\n";
      initializeOp(*P2_linsolve_strategy, P2, P2_lows.ptr());
    }
    *out << "\nP2_lows = " << describe(*P2_lows, verbLevel) << "\n";

    //
    *out << "\nE) Solve P2 for a random RHS ...\n";
    //

    VectorPtr x = createMember(P2->domain());
    VectorPtr b = createMember(P2->range());
    Thyra::randomize(-1.0, +1.0, b.ptr());
    Thyra::assign(x.ptr(), 0.0); // Must give an initial guess!

    Thyra::SolveStatus<double>
      solveStatus = solve<double>( *P2_lows, Thyra::NOTRANS, *b, x.ptr() );

    *out << "\nSolve status:\n" << solveStatus;

    *out << "\nSolution ||x|| = " << Thyra::norm(*x) << "\n";

    if(showParams) {
      *out << "\nParameter list after use:\n\n";
      paramList->print(*out, PLPrintOptions().indent(2).showTypes(true));
    }

    //
    *out << "\nF) Checking the error in the solution of r=b-P2*x ...\n";
    //

    VectorPtr P2x = Thyra::createMember(b->space());
    Thyra::apply( *P2, Thyra::NOTRANS, *x, P2x.ptr() );
    VectorPtr r = Thyra::createMember(b->space());
    Thyra::V_VmV<double>(r.ptr(), *b, *P2x);

    double
      P2x_nrm = Thyra::norm(*P2x),
      r_nrm = Thyra::norm(*r),
      b_nrm = Thyra::norm(*b),
      r_nrm_over_b_nrm = r_nrm / b_nrm;

    bool result = r_nrm_over_b_nrm <= solveTol;
    if(!result) success = false;

    *out
      << "\n||P2*x|| = " << P2x_nrm << "\n";

    *out
      << "\n||P2*x-b||/||b|| = " << r_nrm << "/" << b_nrm
      << " = " << r_nrm_over_b_nrm << " <= " << solveTol
      << " : " << Thyra::passfail(result) << "\n";

    Teuchos::TimeMonitor::summarize(*out<<"\n");
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose,std::cerr,success)

  if (verbose) {
    if(success)  *out << "\nCongratulations! All of the tests checked out!\n";
    else         *out << "\nOh no! At least one of the tests failed!\n";
  }

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );

}
