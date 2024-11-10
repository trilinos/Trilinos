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
#include "Thyra_LinearOpTester.hpp"
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
 * Example program that reads in a matrix from file then tries to solve it
 * using one linear solver as a preconditioner for a another solver. .
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
  using Teuchos::getParametersFromXmlFile;
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

    CommandLineProcessor clp(false); // Don't throw exceptions

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
    clp.setOption( "verb-level", &verbLevel,
      numVerbLevels, verbLevelValues, verbLevelNames,
      "Verbosity level used for all objects."
      );

    std::string matrixFile = ".";
    clp.setOption( "matrix-file", &matrixFile,
      "Matrix file."
      );

    std::string paramListFile = "";
    clp.setOption( "param-list-file", &paramListFile,
      "Parameter list for preconditioner and solver blocks."
      );

    bool showParams = false;
    clp.setOption( "show-params", "no-show-params", &showParams,
      "Show the parameter list or not."
      );

    bool testPrecIsLinearOp = true;
    clp.setOption( "test-prec-is-linear-op", "test-prec-is-linear-op", &testPrecIsLinearOp,
      "Test if the preconditioner is a linear operator or not."
      );

    double solveTol = 1e-8;
    clp.setOption( "solve-tol", &solveTol,
      "Tolerance for the solution to determine success or failure!"
      );

    clp.setDocString(
      "This example program shows how to use one linear solver (e.g. AztecOO)\n"
      "as a preconditioner for another iterative solver (e.g. Belos).\n"
      );

    // Note: Use --help on the command line to see the above documentation

    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;


    //
    *out << "\nA) Reading in the matrix ...\n";
    //

#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif

    const LinearOpPtr A = readEpetraCrsMatrixFromMatrixMarketAsLinearOp(
      matrixFile, comm, "A");
    *out << "\nA = " << describe(*A,verbLevel) << "\n";

    const RCP<ParameterList> paramList = getParametersFromXmlFile(paramListFile);
    if (showParams) {
      *out << "\nRead in parameter list:\n\n";
      paramList->print(*out, PLPrintOptions().indent(2).showTypes(true));
    }

    //
    *out << "\nB) Get the preconditioner as a forward solver\n";
    //

    const RCP<ParameterList> precParamList = sublist(paramList, "Preconditioner Solver");
    Stratimikos::DefaultLinearSolverBuilder precSolverBuilder;
    precSolverBuilder.setParameterList(precParamList);
    const RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > precSolverStrategy
      = createLinearSolveStrategy(precSolverBuilder);
    //precSolverStrategy->setVerbLevel(verbLevel);

    const LinearOpPtr A_inv_prec = inverse<double>(*precSolverStrategy, A,
      Thyra::SUPPORT_SOLVE_FORWARD_ONLY,
      Teuchos::null, // Use internal solve criteria
      Thyra::IGNORE_SOLVE_FAILURE // Ignore solve failures since this is just a prec
      );
    *out << "\nA_inv_prec = " << describe(*A_inv_prec, verbLevel) << "\n";

    if (testPrecIsLinearOp) {
      *out << "\nTest that the preconditioner A_inv_prec is indeed a linear operator.\n";
      Thyra::LinearOpTester<double> linearOpTester;
      linearOpTester.check_adjoint(false);
      const bool linearOpCheck = linearOpTester.check(*A_inv_prec, out.ptr());
      if (!linearOpCheck) { success = false; }
    }

    //
    *out << "\nC) Create the forward solver using the created preconditioner ...\n";
    //

    const RCP<ParameterList> fwdSolverParamList = sublist(paramList, "Forward Solver");
    Stratimikos::DefaultLinearSolverBuilder fwdSolverSolverBuilder;
    fwdSolverSolverBuilder.setParameterList(fwdSolverParamList);
    const RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > fwdSolverSolverStrategy
      = createLinearSolveStrategy(fwdSolverSolverBuilder);

    const RCP<Thyra::LinearOpWithSolveBase<double> >
      A_lows = fwdSolverSolverStrategy->createOp();

    initializePreconditionedOp<double>( *fwdSolverSolverStrategy, A,
        unspecifiedPrec(A_inv_prec), A_lows.ptr());
    //A_lows->setVerbLevel(verbLevel);
    *out << "\nA_lows = " << describe(*A_lows, verbLevel) << "\n";

    //
    *out << "\nD) Solve the linear system for a random RHS ...\n";
    //

    VectorPtr x = createMember(A->domain());
    VectorPtr b = createMember(A->range());
    Thyra::randomize(-1.0, +1.0, b.ptr());
    Thyra::assign(x.ptr(), 0.0); // Must give an initial guess!

    Thyra::SolveStatus<double>
      solveStatus = solve<double>( *A_lows, Thyra::NOTRANS, *b, x.ptr() );

    *out << "\nSolve status:\n" << solveStatus;

    *out << "\nSolution ||x|| = " << Thyra::norm(*x) << "\n";

    if(showParams) {
      *out << "\nParameter list after use:\n\n";
      paramList->print(*out, PLPrintOptions().indent(2).showTypes(true));
    }

    //
    *out << "\nF) Checking the error in the solution of r=b-A*x ...\n";
    //

    VectorPtr Ax = Thyra::createMember(b->space());
    Thyra::apply( *A, Thyra::NOTRANS, *x, Ax.ptr() );
    VectorPtr r = Thyra::createMember(b->space());
    Thyra::V_VmV<double>(r.ptr(), *b, *Ax);

    double
      Ax_nrm = Thyra::norm(*Ax),
      r_nrm = Thyra::norm(*r),
      b_nrm = Thyra::norm(*b),
      r_nrm_over_b_nrm = r_nrm / b_nrm;

    bool resid_tol_check = ( r_nrm_over_b_nrm <= solveTol );
    if(!resid_tol_check) success = false;

    *out
      << "\n||A*x|| = " << Ax_nrm << "\n";

    *out
      << "\n||A*x-b||/||b|| = " << r_nrm << "/" << b_nrm
      << " = " << r_nrm_over_b_nrm << " <= " << solveTol
      << " : " << Thyra::passfail(resid_tol_check) << "\n";

    Teuchos::TimeMonitor::summarize(*out<<"\n");
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success)

  if (verbose) {
    if(success)  *out << "\nCongratulations! All of the tests checked out!\n";
    else         *out << "\nOh no! At least one of the tests failed!\n";
  }

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
